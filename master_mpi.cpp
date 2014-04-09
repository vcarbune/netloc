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
  , m_epflSolver(m_network, m_config, config.sampling)
{
  srand(time(NULL));
  initializeGroundTruths();
  initializeNodeInfectionTimeMap();
}

void MasterNode::runWithCurrentConfig()
{
  cout << "Running with objective: " <<
      algorithmTypeToString(m_config.objType) << endl;

  time_t startTime = time(NULL);
  SimConfig initialConfig = m_config;

  // Description: results[configStep][groundTruth][observerPercentage]
  vector<vector<result_t>> results[m_config.steps];
  for (int step = 0; step < m_config.steps; ++step, ++m_config) {
    double startTime = time(NULL);
    reset();
    for (size_t idx = 0; idx < m_realizations.size(); ++idx) {
      double startTime = time(NULL);
      const GraphHypothesis& realization = m_realizations[idx];
      vector<result_t> incrementalResults = simulate(idx);
      cout << "Progress: " << 100 * idx / m_realizations.size() << "\% ("
          << difftime(time(NULL), startTime) << "s)" << endl;

      double pcnt = m_config.sampling;
      for (result_t& r : incrementalResults) {
#if DBG
        cout << pcnt << "-" << realization.getSource() << "\t";
        cout << r.first << "\t";
#endif

        // Add distance (in hops) to the source node.
        r.second.push_back(
            m_idToShortestPathsFromSource[idx].GetDat(r.first));

        // Keep a boolean, whether the source was identified or not.
        r.first = realization.getSource() == r.first ? 1 : 0;

#if DBG
        for (size_t i = 0; i < r.second.size(); ++i)
          cout << r.second[i] << "\t";
        cout << endl;
#endif

        pcnt += m_config.sampling;
      }

      results[step].push_back(incrementalResults);
    }
    cout << "Total Time: " << difftime(time(NULL), startTime) << "s\n\n";
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
          m_network, m_config.groundTruthFile.CStr()));
  } else {
    // Generate artificially.
    for (int truth = 0; truth < m_config.groundTruths; ++truth) {
      int sourceId = m_nid[rand() % m_network->GetNodes()];
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
    if (!m_network->IsNode(realization.getSource()))
      cout << "WTF?" << realization.getSource() << endl;
    TSnap::GetShortPath(m_network, realization.getSource(), idToShortestPathsFromSource);
    m_idToShortestPathsFromSource.push_back(idToShortestPathsFromSource);
  }
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

void MasterNode::reset()
{
  // cout << "MasterNode::reset()" << endl;
  if (!m_config.cluster.keep)
    initializeNodeInfectionTimeMap();
  initializeTests();
}

void MasterNode::initializeNodeInfectionTimeMap()
{
  m_histogramInfo.clear();

  int junk = 0;
  MPI::Status status;
  for (int node = 0; node < m_network->GetNodes(); ++node) {
    vector<double> uniqueInfectionTimes;
    if (!m_config.ignoreTime) {
      int infectionTimeVectorSize[m_config.mpi.nodes];
      MPI::COMM_WORLD.Gather(&junk, 1, MPI::INT,
          infectionTimeVectorSize, 1, MPI::INT, MPI_MASTER);

      for (int worker = 1; worker <  m_config.mpi.nodes; ++worker) {
        int items = infectionTimeVectorSize[worker];
        double infectionTimes[items];
        MPI::COMM_WORLD.Recv(
            infectionTimes, items, MPI::DOUBLE, worker, 0, status);

        vector<double> result(items + uniqueInfectionTimes.size());
        vector<double>::iterator it = set_union(
            uniqueInfectionTimes.begin(), uniqueInfectionTimes.end(),
            infectionTimes, infectionTimes + items, result.begin());
        result.resize(it - result.begin());

        uniqueInfectionTimes.swap(result);
      }

      // The infection time vector has the following structure:
      // Nbins, minT1, minT2, ..., minNbins[, infty]
      uniqueInfectionTimes = createHistograms(uniqueInfectionTimes);
    } else {
      uniqueInfectionTimes.push_back(INFECTED_TRUE);
      uniqueInfectionTimes.push_back(INFECTED_FALSE);
    }

#if DBG
    cout << node << "(" << m_nid[node] << "): ";
    for (size_t i = 0; i < uniqueInfectionTimes.size(); ++i)
      cout << uniqueInfectionTimes[i] << " ";
    cout << endl;
#endif

    // Is it really useful to store this in the master node?
    m_histogramInfo.push_back(uniqueInfectionTimes);

    int size = uniqueInfectionTimes.size();
    MPI::COMM_WORLD.Bcast(&size, 1, MPI::INT, MPI_MASTER);
    MPI::COMM_WORLD.Bcast(&uniqueInfectionTimes.front(), size, MPI::DOUBLE, MPI_MASTER);

    m_nodeInfectionTimes.push_back(
        createDiscreteTimeValuesFromHistogramBounds(uniqueInfectionTimes));

#if DBG
    cout << node << "(" << m_nid[node] << "): ";
    for (size_t i = 0; i < m_nodeInfectionTimes[node].size(); ++i)
      cout << m_nodeInfectionTimes[node][i] << " ";
    cout << endl;
#endif
  }
}

void MasterNode::initializeTests()
{
  m_tests.clear();
  if (m_config.objType == EPFL_ML || m_config.objType == EPFL_EC2)
    return;
  else if (m_config.objType == RANDOM || m_config.objType == VOI ||
           m_config.objType == EC2_HIGH)
    initializeTestVector();
  else
    initializeTestHeap();
}

void MasterNode::initializeTestHeap()
{
  // cout << "MasterNode::initializeTestHeap()" << endl;
  double weightSum = 0;
  double currentGraphWeight = computeCurrentWeight(&weightSum);
  for (int node = 0; node < m_config.nodes; ++node) {
    GraphTest test(node);
    recomputeTestScore(test, weightSum, currentGraphWeight);
    m_tests.push_back(test);
  }

  // Inform worker nodes to stop processing requests for recomputations.
  int invalidNode = -1;
  MPI::COMM_WORLD.Bcast(&invalidNode, 1, MPI::INT, MPI_MASTER);

  make_heap<vector<GraphTest>::iterator, TestCompareFunction>(
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
  if (m_config.objType == EPFL_ML || m_config.objType == EPFL_EC2)
    return simulateEPFLPolicy(realizationIdx);

  return simulateAdaptivePolicy(realizationIdx);
}

vector<result_t> MasterNode::simulateEPFLPolicy(int realizationIdx)
{
  vector<result_t> results;
  for (double pcnt = m_config.sampling; pcnt < (1.00 + m_config.sampling);
       pcnt += m_config.sampling) {
    m_epflSolver = EPFLSolver(m_network, m_config, pcnt);
    if (m_config.objType == EPFL_EC2) {
      m_epflSolver.setObserverList(m_ec2observers[realizationIdx], pcnt);
    }

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
  m_previousTests.clear();

  const GraphHypothesis& realization = m_realizations[realizationIdx];

  // Request scores for each of the tests.
  double nextPcnt = m_config.sampling;
  for (int count = 0; count < m_config.nodes; ++count) {
    GraphTest nextTest = selectNextTest(tests);
#if DBG
    cout << count << ". " << m_nid[nextTest.getNodeId()] << " " << nextTest.getScore() << endl;
#endif

    // Inform the workers about the selected test.
    int nodeId = nextTest.getNodeId();
    double score = nextTest.getScore();

    // Use the observer.
    double infection = realization.getInfectionTime(nextTest.getNodeId());

    MPI::COMM_WORLD.Bcast(&nodeId, 1, MPI::INT, MPI_MASTER);
    MPI::COMM_WORLD.Bcast(&infection, 1, MPI::DOUBLE, MPI_MASTER);
    MPI::COMM_WORLD.Bcast(&score, 1, MPI::DOUBLE, MPI_MASTER);

    m_previousTests.push_back(make_pair(infection, nodeId));

    if (fabs(count - nextPcnt * m_config.nodes) < 0.5) {
      results.push_back(identifyCluster(realizationIdx));
      nextPcnt += m_config.sampling;
    }

    if (m_config.objType == RANDOM)
      continue;
  }
  results.push_back(identifyCluster(realizationIdx));

  // Save observers for each ground truth so that they can be used for EPFL.
  if (m_config.objType == EC2) {
    vector<int> observers;
    for (size_t i = 0; i < m_previousTests.size(); ++i)
      observers.push_back(m_previousTests[i].second);
    m_ec2observers.push_back(observers);
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

GraphTest MasterNode::selectHighestDegreeTest(vector<GraphTest>& tests)
{
  size_t idx = 0;
  int highestDegree = m_network->GetNI(tests[idx].getNodeId()).GetOutDeg();

  for (size_t i = 1; i < tests.size(); ++i) {
    int nodeDegree = m_network->GetNI(tests[i].getNodeId()).GetOutDeg();
    if (nodeDegree > highestDegree) {
      idx = i;
      highestDegree = nodeDegree;
    }
  }

  GraphTest test = tests[idx];
  test.setScore(highestDegree);
  tests.erase(tests.begin() + idx);
  return test;
}

GraphTest MasterNode::selectVoITest(vector<GraphTest>& tests)
{
  double weightSum = 0;
  computeCurrentWeight(&weightSum);

  size_t maxIdx = 0;
  for (size_t idx = 0; idx < tests.size(); ++idx) {
    int currentTestNode = tests[idx].getNodeId();
    MPI::COMM_WORLD.Bcast(&currentTestNode, 1, MPI::INT, MPI_MASTER);

    // For each infection time, a different prior, a different mass.
    int totalInfectionTimes = m_nodeInfectionTimes[currentTestNode].size();

    double nullVals[totalInfectionTimes];
    memset(nullVals, 0, totalInfectionTimes * sizeof(*nullVals));

    double priors[totalInfectionTimes];
    MPI::COMM_WORLD.Reduce(
        nullVals, priors, totalInfectionTimes, MPI::DOUBLE, MPI::SUM, MPI_MASTER);
    for (int i = 0; i < totalInfectionTimes; ++i)
      priors[i] = priors[i] / weightSum;

    // Broadcast priors back to workers.
    MPI::COMM_WORLD.Bcast(priors, totalInfectionTimes, MPI::DOUBLE, MPI_MASTER);

    double score;
    double nullVal = 0;
    MPI::COMM_WORLD.Reduce(&nullVal, &score, 1, MPI::DOUBLE, MPI::MAX,
        MPI_MASTER);

    tests[idx].setScore(score);
    if (tests[idx].getScore() > tests[maxIdx].getScore())
      maxIdx = idx;
  }

  int currentTestNode = -1;
  MPI::COMM_WORLD.Bcast(&currentTestNode, 1, MPI::INT, MPI_MASTER);

  GraphTest test = tests[maxIdx];
  tests.erase(tests.begin() + maxIdx);
  return test;
}

GraphTest MasterNode::selectNextTest(vector<GraphTest>& tests)
{
  // cout << "MasterNode::selectNextTest()" << endl;

  if (m_config.objType == RANDOM)
    return selectRandomTest(tests);
  if (m_config.objType == EC2_HIGH)
    return selectHighestDegreeTest(tests);
  if (m_config.objType == VOI)
    return selectVoITest(tests);

  double weightSum = 0;
  double currentMass = computeCurrentWeight(&weightSum);

  TestCompareFunction tstCmpFcn;
#if DBG
  int count = 0;
  cout << "Current Mass: " << currentMass << endl;
#endif
  do {
    // Remove the top test from the heap.
    GraphTest crtTop = tests.front();
    pop_heap<vector<GraphTest>::iterator, TestCompareFunction>(
        tests.begin(), tests.end(), TestCompareFunction());
    tests.pop_back();

    // Exit early if it's the last element in the heap.
    if (!tests.size()) {
      int invalidNode = -1;
      MPI::COMM_WORLD.Bcast(&invalidNode, 1, MPI::INT, MPI_MASTER);

      return crtTop;
    }

    // Recompute its score and keep it if it stays on top.
    recomputeTestScore(crtTop, weightSum, currentMass);
    if (tstCmpFcn(tests.front(), crtTop)) {
#if DBG
      cout << "Heap cost: " << count << endl;
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

void MasterNode::recomputeTestScore(
    GraphTest& test, double weightSum, double currentMass)
{
  // Broadcast current node for which scores are recomputed.
  int currentTestNode = test.getNodeId();
  MPI::COMM_WORLD.Bcast(&currentTestNode, 1, MPI::INT, MPI_MASTER);

  // For each infection time, a different prior, a different mass.
  int totalInfectionTimes = m_nodeInfectionTimes[currentTestNode].size();

  double nullVals[totalInfectionTimes];
  memset(nullVals, 0, totalInfectionTimes * sizeof(*nullVals));

  // Recompute test priors, for each infection time!
  double priors[totalInfectionTimes];
  MPI::COMM_WORLD.Reduce(
      nullVals, priors, totalInfectionTimes, MPI::DOUBLE, MPI::SUM, MPI_MASTER);

  // Recompute masses.
  double mass[totalInfectionTimes];
  double sqMass[totalInfectionTimes];
  MPI::COMM_WORLD.Reduce(
      nullVals, mass, totalInfectionTimes, MPI::DOUBLE, MPI::SUM, MPI_MASTER);
  MPI::COMM_WORLD.Reduce(
      nullVals, sqMass, totalInfectionTimes, MPI::DOUBLE, MPI::SUM, MPI_MASTER);

  double expectedMass = 0.0;
  for (int i = 0; i < totalInfectionTimes; ++i) {
    priors[i] = priors[i] / weightSum;
    if (m_config.objType == EC2)
      expectedMass += priors[i] * (mass[i] * mass[i] - sqMass[i]);
    else if (m_config.objType == GBS)
      expectedMass += priors[i] * mass[i];
    else
      cout << "Warning: not yet supported :(" << endl;
  }
  test.setScore(currentMass - expectedMass);
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
    if (!m_network->IsNode(allClusterWeights[i].second))
      cerr << "Whoa, not in the network!?" << allClusterWeights[i].second;
  }

  // Obtain rank of true solution.
  int rank = 0;
  for (; rank < m_config.nodes; ++rank)
    if (allClusterWeights[rank].second == realSource)
      break;

  pair<int, vector<double>> result;
  // Identified solution (source node id).
  result.first = allClusterWeights[0].second;
  // How much mass was in the solution cluster?
  result.second.push_back(realSourceMass);
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
  dumpStream << fixed << setprecision(3);

  for (int step = 0; step < initialConfig.steps; ++step, ++initialConfig) {
    int vars = results[step][0][0].second.size();

    // For each observer percentage (1%, 2%, ..., 100%), compute averages.
    double pcnt = m_config.sampling;
    for (size_t observers = 0; observers < results[step][0].size();
         ++observers, pcnt += m_config.sampling) {
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
