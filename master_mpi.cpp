#include "master_mpi.h"

// SNAP redefines max and min!
#undef max
#undef min

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <limits>
#include <fstream>
#include <utility>

using namespace std;

MasterNode::MasterNode(SimConfig config)
  : MPINode(config)
  , m_epflSolver(m_network, m_config)
{
  srand(time(NULL));
  initializeGroundTruths();
}

void MasterNode::run()
{
  SimConfig initialConfig = m_config;
  if (m_config.objType != -1) {
    runWithCurrentConfig();
    return;
  }
  for (int obj = EC2; obj <= EPFL_ML; ++obj) {
    m_config = initialConfig;
    m_config.setObjType(static_cast<AlgorithmType>(obj));
    m_config.logfile = TStr::Fmt("%s%d_%s.log", m_config.logfile.CStr(),
        m_config.nodes, algorithmTypeToString(static_cast<AlgorithmType>(obj)));

    runWithCurrentConfig();
  }
}

void MasterNode::runWithCurrentConfig()
{
  cout << "Running with objective: " << m_config.objType << endl;

  time_t startTime = time(NULL);
  SimConfig initialConfig = m_config;

  // Description: results[configStep][groundTruth][observerPercentage]
  vector<vector<result_t>> results[m_config.steps];
  for (int step = 0; step < m_config.steps; ++step, ++m_config) {
    double startTime = time(NULL);
    initializeTests();
    cout << "Nodes were initialized in " <<
        difftime(time(NULL), startTime) << "s." << endl;
    for (size_t idx = 0; idx < m_realizations.size(); ++idx) {
      const GraphHypothesis& realization = m_realizations[idx];
      vector<result_t> incrementalResults = simulate(idx);

      cout << m_config << "-" << realization.getSource() << "\t";

      for (result_t& r : incrementalResults) {
        cout << r.first << "\t";

        // Add distance (in hops) to the source node.
        r.second.push_back(
            m_idToShortestPathsFromSource[idx].GetDat(r.first));

        // Keep a boolean, whether the source was identified or not.
        r.first = realization.getSource() == r.first ? 1 : 0;

        for (size_t i = 0; i < r.second.size(); ++i)
          cout << r.second[i] << "\t";
        cout << endl;
      }

      results[step].push_back(incrementalResults);
    }
    cout << "Time: " << difftime(time(NULL), startTime) << "s\n\n";
  }

  processResults(results, initialConfig);
  cout << "Time: " << difftime(time(NULL), startTime) << "s" << endl;
}

void MasterNode::initializeGroundTruths()
{
  // Read from file, if ground truth is given.
  if (!m_config.groundTruthFile.Empty()) {
    cout << "Reading ground truth from " <<
        m_config.groundTruthFile.CStr() << endl;

    m_realizations.push_back(GraphHypothesis::readHypothesisFromFile(
          m_config.groundTruthFile.CStr()));
  } else {
    // Generate artificially.
    for (int truth = 0; truth < m_config.groundTruths; ++truth) {
      int sourceId = rand() % m_network->GetNodes();
      switch (m_config.infType) {
        case BETA:
          m_realizations.push_back(GraphHypothesis::generateHypothesis(
                m_network, sourceId, m_config.cluster));
          break;
        case GAUSSIAN:
          m_realizations.push_back(
              GraphHypothesis::generateHypothesisUsingGaussianModel(
                m_network, sourceId, m_config.cluster));
          break;
      }
    }
  }

  for (const GraphHypothesis& realization : m_realizations) {
    TIntH idToShortestPathsFromSource;
    TSnap::GetShortPath(m_network, realization.getSource(), idToShortestPathsFromSource);
    m_idToShortestPathsFromSource.push_back(idToShortestPathsFromSource);
  }
}

void MasterNode::computeTestPriors(double currentWeight)
{
  double testPriors[m_config.nodes];
  double nullVals[m_config.nodes];
  memset(nullVals, 0, m_config.nodes * sizeof(*nullVals));

  // Get from workers the node count for each test node.
  MPI::COMM_WORLD.Reduce(nullVals, testPriors, m_config.nodes,
      MPI::DOUBLE, MPI::SUM, MPI_MASTER);

  for (int testNode = 0; testNode < m_config.nodes; ++testNode) {
    testPriors[testNode] /= currentWeight;
    m_testsPrior[testNode] = testPriors[testNode];
  }

  // Broadcast back to workers computed test priors.
  MPI::COMM_WORLD.Bcast(testPriors, m_config.nodes, MPI::DOUBLE, MPI_MASTER);
}

double MasterNode::computeCurrentWeight(double *weightSum) {
  // Compute current mass of all the clusters.
  double currentMass = 0.0;
  double crtSum[2] = {0.0, 0.0};
  MPI::COMM_WORLD.Reduce(&currentMass, &crtSum[0], 1,
      MPI::DOUBLE, MPI::SUM, MPI_MASTER);
  MPI::COMM_WORLD.Reduce(&currentMass, &crtSum[1], 1,
      MPI::DOUBLE, MPI::SUM, MPI_MASTER);

  currentMass = crtSum[0];
  if (weightSum)
    *weightSum = currentMass;
  if (m_config.objType == EC2) // EC2
    currentMass = crtSum[0] * crtSum[0] - crtSum[1];

  // cout << "Weight Sum: " << *weightSum << " Mass: " << currentMass << endl;
  return currentMass;
}

void MasterNode::initializeTests()
{
  if (m_config.objType == RANDOM || m_config.objType == VOI)
    initializeTestVector();
  else if (m_config.objType != EPFL_ML)
    initializeTestHeap();
}

void MasterNode::initializeTestHeap()
{
  double nullVals[m_config.nodes];
  double sums[m_config.objSums][m_config.nodes];

  MPI::Status status;
  memset(sums, 0, m_config.objSums * m_config.nodes * sizeof(**sums));
  memset(nullVals, 0, m_config.nodes * sizeof(*nullVals));

  // Compute current mass of all the clusters.
  double weightSum = 0;
  double currentMass = computeCurrentWeight(&weightSum);
  // Reinitialize test priors
  computeTestPriors(weightSum);

  for (int s = 0; s < m_config.objSums; ++s)
    MPI::COMM_WORLD.Reduce(nullVals, sums[s], m_config.nodes,
        MPI::DOUBLE, MPI::SUM, MPI_MASTER);

  m_tests.clear();
  for (int test = 0; test < m_config.nodes; ++test) {
    double score = 0.0;
    double testPositivePb = m_testsPrior[test];

    if (m_config.objType == EC2) {
      double positiveMass = sums[POSITIVE_SUM][test] * sums[POSITIVE_SUM][test] -
          sums[POSITIVE_DIAG_SUM][test];
      double negativeMass = sums[NEGATIVE_SUM][test] * sums[NEGATIVE_SUM][test] -
          sums[NEGATIVE_DIAG_SUM][test];
      score = currentMass -
          (testPositivePb * positiveMass + (1 - testPositivePb) * negativeMass);
    } else if (m_config.objType == GBS) {   // GBS
      score = testPositivePb * sums[POSITIVE_SUM][test] +
        (1 - testPositivePb) * sums[NEGATIVE_SUM][test];
      score = currentMass - score;
    } else if (m_config.objType == VOI) {   // VOI
      score = sums[0][test];
    } else {
      cout << "Code path should not get here!" << endl;
    }

    GraphTest graphTest(test);
    graphTest.setScore(score);

    m_tests.push_back(graphTest);
  }

  std::make_heap<vector<GraphTest>::iterator, TestCompareFunction>(
      m_tests.begin(), m_tests.end(), TestCompareFunction());
}

void MasterNode::initializeTestVector()
{
  m_tests.clear();
  for (int testNode = 0; testNode < m_config.nodes; ++testNode)
    m_tests.push_back(GraphTest(testNode));
}

vector<result_t> MasterNode::simulate(int realizationIdx)
{
  if (m_config.objType == EPFL_ML)
    return simulateEPFLPolicy(realizationIdx);

  return simulateAdaptivePolicy(realizationIdx);
}

vector<result_t> MasterNode::simulateEPFLPolicy(int realizationIdx)
{
  vector<result_t> results;

  SimConfig config = m_config;
  for (config.testThreshold = 0.01; config.testThreshold < 1.00;
       config.testThreshold += 0.01) {
    m_epflSolver = EPFLSolver(m_network, config);

    vector<pair<double, int>> clusterSortedScores;
    result_t result =
        m_epflSolver.solve(m_realizations[realizationIdx], clusterSortedScores);
    result.second.push_back(computeNDCG(clusterSortedScores, realizationIdx));
    results.push_back(result);
  }

  return results;
}

vector<result_t> MasterNode::simulateAdaptivePolicy(int realizationIdx)
{
  vector<result_t> results;

  vector<GraphTest> tests(m_tests);
  const GraphHypothesis& realization = m_realizations[realizationIdx];

  // Request scores for each of the tests.
  double nextPcnt = 0.01;
  for (int count = 0; count < m_config.nodes; ++count) {
    GraphTest nextTest = selectNextTest(tests);
    cout << count << ". " << nextTest.getNodeId() << " " << nextTest.getScore() << endl;

    // Inform the workers about the selected test.
    int nodeId = nextTest.getNodeId();
    bool outcome = realization.getTestOutcome(nextTest, m_previousTests);
    int infection = realization.getInfectionTime(nextTest.getNodeId());

    MPI::COMM_WORLD.Bcast(&nodeId, 1, MPI::INT, MPI_MASTER);
    MPI::COMM_WORLD.Bcast(&outcome, 1, MPI::BOOL, MPI_MASTER);
    MPI::COMM_WORLD.Bcast(&infection, 1, MPI::INT, MPI_MASTER);

    m_previousTests.push_back(make_pair(infection, nodeId));

    if (count == static_cast<int>(nextPcnt * m_config.nodes)) {
      results.push_back(identifyCluster(realizationIdx));
      nextPcnt += 0.01;
    }

    if (m_config.objType == RANDOM)
      continue;

    for (GraphTest& test : tests)
      test.setScore((1-m_config.eps) * (1-m_config.eps) * test.getScore());
  }

  return results;
}

GraphTest MasterNode::selectRandomTest(vector<GraphTest>& tests)
{
  size_t randomIdx = rand() % tests.size();
  GraphTest test = tests[randomIdx];

  tests.erase(tests.begin() + randomIdx);
  return test;
}

GraphTest MasterNode::selectVOITest(vector<GraphTest>& tests,
    double currentWeight)
{
  double maxTestScore = -numeric_limits<double>::max();
  size_t maxTestIndex = 0;

  computeTestPriors(currentWeight);
  for (size_t i = 0; i < tests.size(); ++i) {
    GraphTest& test = tests[i];
    recomputeTestScore(test, currentWeight, currentWeight);
    if (test.getScore() > maxTestScore) {
      maxTestScore = test.getScore();
      maxTestIndex = i;
    }
  }

  int invalidNode = -1;
  MPI::COMM_WORLD.Bcast(&invalidNode, 1, MPI::INT, MPI_MASTER);

  GraphTest test = tests[maxTestIndex];
  tests.erase(tests.begin() + maxTestIndex);
  return test;
}

GraphTest MasterNode::selectNextTest(vector<GraphTest>& tests)
{
  if (m_config.objType == RANDOM)
    return selectRandomTest(tests);

  double weightSum = 0;
  double currentMass = computeCurrentWeight(&weightSum);
  if (m_config.objType == VOI)
    return selectVOITest(tests, currentMass);

  TestCompareFunction tstCmpFcn;
#if DBG
  int count = 0;
#endif
  do {
    // Remove the top test from the heap.
    GraphTest crtTop = tests.front();
    pop_heap<vector<GraphTest>::iterator, TestCompareFunction>(
        tests.begin(), tests.end(), TestCompareFunction());
    tests.pop_back();

#if DBG
    std::cout << "Top: " << crtTop.getScore() << " Next: " <<
      tests.front().getScore() << std::endl;
#endif

    // Exit early if it's the last element in the heap.
    if (!tests.size()) {
      int invalidNode = -1;
      MPI::COMM_WORLD.Bcast(&invalidNode, 1, MPI::INT, MPI_MASTER);

      return crtTop;
    }

    // Recompute its score and keep it if it stays on top.
    recomputeTestScore(crtTop, weightSum, currentMass);

#if DBG
    std::cout << "Rescored Top: " << crtTop.getScore() << std::endl;
#endif

    if (tstCmpFcn(tests.front(), crtTop)) {
#if DBG
      std::cout << "Pushed " << count << " elems back to heap... " << std::endl;
#endif
      int invalidNode = -1;
      MPI::COMM_WORLD.Bcast(&invalidNode, 1, MPI::INT, MPI_MASTER);

      return crtTop;
    }

    // Otherwise push it back to the heap.
    tests.push_back(crtTop);
#if DBG
    count++;
#endif
    push_heap<vector<GraphTest>::iterator, TestCompareFunction>(
        tests.begin(), tests.end(), TestCompareFunction());
  } while (true);

  cout << "SHOULD NOT BE REACHED" << endl;
  return tests.front();
}

void MasterNode::recomputeTestScoreEC2(
    GraphTest& test, double currentMass, double* nullVals)
{
  double sums[m_config.objSums];
  MPI::COMM_WORLD.Reduce(nullVals, sums, m_config.objSums,
      MPI::DOUBLE, MPI::SUM, MPI_MASTER);

  pair<double, double> mass(
      sums[POSITIVE_SUM] * sums[POSITIVE_SUM] - sums[POSITIVE_DIAG_SUM],
      sums[NEGATIVE_SUM] * sums[NEGATIVE_SUM] - sums[NEGATIVE_DIAG_SUM]);
  double expectedMass = m_testsPrior[test.getNodeId()] * mass.first +
      (1 - m_testsPrior[test.getNodeId()]) * mass.second;

  test.setScore(currentMass - expectedMass);
}

void MasterNode::recomputeTestScoreGBS(
    GraphTest& test, double currentMass, double* nullVals)
{
  double sums[m_config.objSums];
  MPI::COMM_WORLD.Reduce(nullVals, sums, m_config.objSums,
      MPI::DOUBLE, MPI::SUM, MPI_MASTER);

  double expectedMass = m_testsPrior[test.getNodeId()] * sums[POSITIVE_SUM] +
    (1 - m_testsPrior[test.getNodeId()]) * sums[NEGATIVE_SUM];

  test.setScore(currentMass - expectedMass);
}

void MasterNode::recomputeTestScoreVOI(GraphTest& test)
{
  double nullVal = 0;
  double maxClusterMassWithTest;
  MPI::COMM_WORLD.Reduce(&nullVal, &maxClusterMassWithTest, 1,
      MPI::DOUBLE, MPI::MAX, MPI_MASTER);
  test.setScore(maxClusterMassWithTest);
}


void MasterNode::recomputeTestScore(
    GraphTest& test, double weightSum, double currentMass)
{
  int currentTestNode = test.getNodeId();
  MPI::COMM_WORLD.Bcast(&currentTestNode, 1, MPI::INT, MPI_MASTER);

  double nullVals[m_config.objSums];
  memset(nullVals, 0, m_config.objSums * sizeof(*nullVals));

  // Recompute test prior.
  double junk = 0.0;
  double positiveTestPrior;
  MPI::COMM_WORLD.Reduce(&junk, &positiveTestPrior, 1,
        MPI::DOUBLE, MPI::SUM, MPI_MASTER);
  m_testsPrior[test.getNodeId()] = positiveTestPrior / weightSum;
  if (m_config.objType == EC2)
    recomputeTestScoreEC2(test, currentMass, nullVals);
  else if (m_config.objType == GBS)
    recomputeTestScoreGBS(test, currentMass, nullVals);
  else if (m_config.objType == VOI)
    recomputeTestScoreVOI(test);
  else
    cout << "This part of the code should never be reached.." << endl;
}

double MasterNode::computeNDCG(
    const vector<pair<double, int>>& clusterSortedScores,
    int realizationIdx) const
{
  const TIntH& nodeRelevance = m_idToShortestPathsFromSource[realizationIdx];

  int maxDistance = 0;
  for (int i = 0; i < m_config.nodes; ++i)
    maxDistance = max(nodeRelevance.GetDat(i).Val, maxDistance);

  double dcg = maxDistance - nodeRelevance.GetDat(clusterSortedScores[0].second);
  for (int i = 1; i < m_config.ndcgN; ++i)
    dcg += (double) (maxDistance -
        nodeRelevance.GetDat(clusterSortedScores[i].second)) / log2(i+1);

  vector<int> idealRelevanceOrdering;
  for (int i = 0; i < m_config.nodes; ++i)
    idealRelevanceOrdering.push_back(maxDistance - nodeRelevance.GetDat(i).Val);
  sort(idealRelevanceOrdering.begin(), idealRelevanceOrdering.end(),
       greater<int>());

  double idcg = idealRelevanceOrdering[0];
  for (int i = 1; i < m_config.ndcgN; ++i)
    idcg += (double) idealRelevanceOrdering[i] / log2(i+1);

  return dcg / idcg;
}

result_t MasterNode::identifyCluster(int realizationIdx)
{
  int realSource = m_realizations[realizationIdx].getSource();

  // Gather information from all processes.
  int clusterNodes[m_config.nodes];
  double clusterWeight[m_config.nodes];
  vector<pair<double, int>> allClusterWeights;

  double realSourceMass = numeric_limits<double>::max();
  int nodes;
  double mass = 0.0;

  MPI::Status status;
  for (int worker = 1; worker < m_config.mpi.nodes; ++worker) {
    MPI::COMM_WORLD.Recv(&nodes, 1, MPI::INT, worker, 0, status);
    MPI::COMM_WORLD.Recv(&clusterNodes, nodes, MPI::INT, worker, 0, status);
    MPI::COMM_WORLD.Recv(&clusterWeight, nodes, MPI::DOUBLE, worker, 0, status);

    for (int i = 0; i < nodes; ++i) {
      mass += clusterWeight[i];
      allClusterWeights.push_back(make_pair(clusterWeight[i], clusterNodes[i]));
    }
  }

  // Sort to obtain relevance of each node.
  sort(allClusterWeights.begin(), allClusterWeights.end(),
       greater<pair<double, int>>());
  for (size_t i = 0; i < allClusterWeights.size(); ++i) {
    allClusterWeights[i].first = allClusterWeights[i].first * 100 / mass;
    if (allClusterWeights[i].second == realSource)
      realSourceMass = allClusterWeights[i].first;
  }

  // Obtain rank of true solution.
  int rank = 0;
  for (; rank < m_config.nodes; ++rank)
    if (allClusterWeights[rank].second == realSource)
      break;

  pair<int, vector<double>> result;
  // Identified solution (source node id).
  result.first = allClusterWeights[0].second;
  // Confidence in the solution identified (probability)
  result.second.push_back(allClusterWeights[0].first);
  // Difference between the probability of solution and the correct one.
  result.second.push_back(allClusterWeights[0].first - realSourceMass);
  // Rank of the correct solution.
  result.second.push_back(rank);
  // NDCG metric.
  result.second.push_back(computeNDCG(allClusterWeights, realizationIdx));

  return result;
}

void MasterNode::processResults(vector<vector<result_t>> *results,
    SimConfig& initialConfig)
{
  ofstream dumpStream;
  dumpStream.open(initialConfig.logfile.CStr());
  dumpStream << fixed << std::setprecision(3);

  for (int step = 0; step < initialConfig.steps; ++step, ++initialConfig) {
    int vars = results[step][0][0].second.size();

    // For each observer percentage (1%, 2%, ..., 100%), compute averages.
    double pcnt = 0.01;
    for (size_t observers = 0; observers < results[step][0].size();
         ++observers, pcnt += 0.01) {
      double identificationCount = 0.0;
      double averages[vars];

      memset(averages, 0, vars * sizeof(*averages));
      for (int truth = 0; truth < initialConfig.groundTruths; ++truth) {
        identificationCount += results[step][truth][observers].first;
        for (int var = 0; var < vars; ++var) {
          averages[var] += results[step][truth][observers].second[var];
        }
      }

      // Normalize computations.
      for (int var = 0; var < vars; ++var)
        averages[var] /= initialConfig.groundTruths;
      identificationCount /= initialConfig.groundTruths;

      // Compute standard error.
      double stderrs[vars];
      double squaredTruths =
          initialConfig.groundTruths * initialConfig.groundTruths;

      memset(stderrs, 0, vars * sizeof(*stderrs));
      for (int truth = 0; truth < initialConfig.groundTruths; ++truth) {
        for (int var = 0; var < vars; ++var) {
          double tmp = results[step][truth][observers].second[var] - averages[var];
          stderrs[var] += tmp * tmp / squaredTruths;
        }
      }

      for (int var = 0; var < vars; ++var)
        stderrs[var] = sqrt(stderrs[var]);

      // Dump to output streams.
      cout << pcnt << "\t";
      dumpStream << pcnt << "\t";
      for (int var = 0; var < vars; ++var) {
        dumpStream << averages[var] << "\t" << stderrs[var] << "\t";
        cout << averages[var] << "\t" << stderrs[var] << "\t";
      }
      dumpStream << identificationCount << endl;
      cout << identificationCount << endl;
    }
  }

  dumpStream.close();
}
