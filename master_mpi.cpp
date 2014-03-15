#include "master_mpi.h"

// SNAP redefines max and min!
#undef max
#undef min

#include <cmath>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <limits>
#include <fstream>
#include <utility>

MasterNode::MasterNode(SimConfig config)
  : MPINode(config)
{
  initializeGroundTruths();
}

void MasterNode::run()
{
  time_t startTime = time(NULL);
  SimConfig initialConfig = m_config;

  vector<result_t> results[m_config.steps];
  for (int step = 0; step < m_config.steps; ++step, ++m_config) {
    initializeTestHeap();
    double startTime = time(NULL);
    for (const GraphHypothesis& realization : m_realizations) {
      result_t result = simulate(realization);
      // Add distance (in hops) to the source node.
      result.second.push_back(TSnap::GetShortPath(m_network,
            realization.getSource(), result.first));
      // Keep a boolean, whether the source was identified or not.
      result.first = realization.getSource() == result.first ? 1 : 0;
      results[step].push_back(result);

      cout << m_config << "-" << realization.getSource() << "\t";
      for (size_t i = 0; i < result.second.size(); ++i)
        cout << result.second[i] << "\t";
      cout << endl;
    }
    cout << "Time: " << difftime(time(NULL), startTime) << "s\n\n";
  }

  processResults(results, initialConfig);
  cout << "Time: " << difftime(time(NULL), startTime) << "s" << endl;
}

void MasterNode::initializeGroundTruths()
{
  // Read from file, if ground truth is given.
  if (!m_config.groundTruth.Empty()) {
    cout << "Reading ground truth from " << m_config.groundTruth.CStr() << endl;
    m_realizations.push_back(
        GraphHypothesis::readHypothesisFromFile(m_config.groundTruth.CStr()));
    return;
  }

  // Generate artificially otherwise.
  for (int truth = 0; truth < m_config.groundTruths; ++truth)
    m_realizations.push_back(GraphHypothesis::generateHypothesis(
          m_network, rand() % m_network->GetNodes(),
          m_config.cascadeBound, m_config.beta));
}

double MasterNode::computeCurrentMass() {
  // Compute current mass of all the clusters.
  double currentWeight = 0.0;
  double crtSum[2];
  crtSum[0] = 0.0; crtSum[1] = 0.0;
  MPI::COMM_WORLD.Reduce(&currentWeight, &crtSum[0], 1,
      MPI::DOUBLE, MPI::SUM, MPI_MASTER);
  MPI::COMM_WORLD.Reduce(&currentWeight, &crtSum[1], 1,
      MPI::DOUBLE, MPI::SUM, MPI_MASTER);

  currentWeight = crtSum[0];
  if (m_config.objType == 0) // EC2
    currentWeight = crtSum[0] * crtSum[0] - crtSum[1];

  return currentWeight;
}

void MasterNode::initializeTestHeap()
{
  double nullVals[m_config.nodes];
  double sums[m_config.objSums][m_config.nodes];

  MPI::Status status;
  memset(sums, 0, m_config.objSums * m_config.nodes * sizeof(**sums));
  memset(nullVals, 0, m_config.nodes * sizeof(*nullVals));

  // Get from workers the node count for each test node.
  MPI::COMM_WORLD.Reduce(&nullVals, sums[0], m_config.nodes,
      MPI::DOUBLE, MPI::SUM, MPI_MASTER);

  int totalHypothesis = m_config.nodes * m_config.clusterSize;
  for (int testNode = 0; testNode < m_config.nodes; ++testNode) {
    m_testsPrior[testNode] = sums[0][testNode] / totalHypothesis;
    sums[0][testNode] = m_testsPrior[testNode];
  }

  // Broadcast back to workers computed test priors.
  MPI::COMM_WORLD.Bcast(sums[0], m_config.nodes, MPI::DOUBLE, MPI_MASTER);

  // Get from workers positive & negative mass of tests still in the loop.
  if (m_config.objType != 2) {
    for (int s = 0; s < m_config.objSums; ++s)
      MPI::COMM_WORLD.Reduce(&nullVals, sums[s], m_config.nodes,
          MPI::DOUBLE, MPI::SUM, MPI_MASTER);
  } else {
    MPI::COMM_WORLD.Reduce(&nullVals, sums[0], m_config.nodes,
        MPI::DOUBLE, MPI::MAX, MPI_MASTER);
  }

  // Compute current mass of all the clusters.
  double currentWeight = computeCurrentMass();

  // Store on the fly the test node probabilities.
  m_tests.clear();

  for (int test = 0; test < m_config.nodes; ++test) {
    double score = 0.0;
    double testPositivePb = m_testsPrior[test];

    if (m_config.objType == 0) {          // EC2
      double positiveMass = sums[POSITIVE_SUM][test] * sums[POSITIVE_SUM][test] -
          sums[POSITIVE_DIAG_SUM][test];
      double negativeMass = sums[NEGATIVE_SUM][test] * sums[NEGATIVE_SUM][test] -
          sums[NEGATIVE_DIAG_SUM][test];
      score = currentWeight -
          (testPositivePb * positiveMass + (1 - testPositivePb) * negativeMass);
    } else if (m_config.objType == 1) {   // GBS
      score = testPositivePb * sums[POSITIVE_SUM][test] +
        (1 - testPositivePb) * sums[NEGATIVE_SUM][test];
      score = currentWeight - score;
    } else if (m_config.objType == 2) {   // VOI
      score = sums[0][test];
    } else if (m_config.objType == 3) {   // RANDOM
      cout << "Code path should not get here!" << endl;
    }

    GraphTest graphTest(test);
    graphTest.setScore(score);

    m_tests.push_back(graphTest);
  }

  std::make_heap<vector<GraphTest>::iterator, TestCompareFunction>(
      m_tests.begin(), m_tests.end(), TestCompareFunction());
}

pair<int, vector<double>> MasterNode::simulate(
    const GraphHypothesis& realization)
{
  vector<GraphTest> tests(m_tests);

  // Request scores for each of the tests.
  int totalTests = m_config.testThreshold * m_config.nodes;
  for (int count = 0; count < totalTests; ++count) {
    GraphTest nextTest = selectNextTest(tests);
    cout << count << ". " << nextTest.getNodeId() << " " << nextTest.getScore() << endl;

    // Inform the workers about the selected test.
    int nodeId = nextTest.getNodeId();
    bool outcome = realization.getTestOutcome(nextTest);
    int infection = realization.getInfectionTime(nextTest.getNodeId());

    MPI::COMM_WORLD.Bcast(&nodeId, 1, MPI::INT, MPI_MASTER);
    MPI::COMM_WORLD.Bcast(&outcome, 1, MPI::BOOL, MPI_MASTER);
    MPI::COMM_WORLD.Bcast(&infection, 1, MPI::INT, MPI_MASTER);

    for (GraphTest& test : tests)
      test.setScore((1-EPS) * (1-EPS) * test.getScore());
  }

  // Identify the cluster where the mass is concentrated.
  result_t solution = identifyCluster(realization.getSource());
  return solution;
}

GraphTest MasterNode::selectNextTest(vector<GraphTest>& tests)
{
  double currentWeight = computeCurrentMass();
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
    recomputeTestScore(crtTop, currentWeight);

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
    GraphTest& test, double currentWeight, double* nullVals)
{
  double sums[m_config.objSums];
  MPI::COMM_WORLD.Reduce(nullVals, sums, m_config.objSums,
      MPI::DOUBLE, MPI::SUM, MPI_MASTER);

  pair<double, double> mass(
      sums[POSITIVE_SUM] * sums[POSITIVE_SUM] - sums[POSITIVE_DIAG_SUM],
      sums[NEGATIVE_SUM] * sums[NEGATIVE_SUM] - sums[NEGATIVE_DIAG_SUM]);
  double expectedMass = m_testsPrior[test.getNodeId()] * mass.first +
      (1 - m_testsPrior[test.getNodeId()]) * mass.second;

  test.setScore(currentWeight - expectedMass);
}

void MasterNode::recomputeTestScoreGBS(
    GraphTest& test, double currentWeight, double* nullVals)
{
  double sums[m_config.objSums];
  MPI::COMM_WORLD.Reduce(nullVals, sums, m_config.objSums,
      MPI::DOUBLE, MPI::SUM, MPI_MASTER);

  double expectedMass = m_testsPrior[test.getNodeId()] * sums[POSITIVE_SUM] +
    (1 - m_testsPrior[test.getNodeId()]) * sums[NEGATIVE_SUM];

  test.setScore(currentWeight - expectedMass);
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
    GraphTest& test, double currentWeight)
{
  int currentTestNode = test.getNodeId();
  MPI::COMM_WORLD.Bcast(&currentTestNode, 1, MPI::INT, MPI_MASTER);

  double nullVals[m_config.objSums];
  memset(nullVals, 0, m_config.objSums * sizeof(*nullVals));

  if (m_config.objType == 0)
    recomputeTestScoreEC2(test, currentWeight, nullVals);
  else if (m_config.objType == 1)
    recomputeTestScoreGBS(test, currentWeight, nullVals);
  else if (m_config.objType == 2)
    recomputeTestScoreVOI(test);
  else
    cout << "Configuration NOT yet supported!" << endl;
}

result_t MasterNode::identifyCluster(int realSource)
{
  // Gather information from all processes.
  int clusterNodes[m_config.nodes];
  double clusterWeight[m_config.nodes];
  double allClusterWeights[m_config.nodes];

  int nodes;
  double mass = 0.0;

  MPI::Status status;
  for (int worker = 1; worker < m_config.mpi.nodes; ++worker) {
    MPI::COMM_WORLD.Recv(&nodes, 1, MPI::INT, worker, 0, status);
    MPI::COMM_WORLD.Recv(&clusterNodes, nodes, MPI::INT, worker, 0, status);
    MPI::COMM_WORLD.Recv(&clusterWeight, nodes, MPI::DOUBLE, worker, 0, status);
    for (int i = 0; i < nodes; ++i) {
      mass += clusterWeight[i];
      allClusterWeights[clusterNodes[i]] = clusterWeight[i];
    }
  }

  int maxIndex = 0;
  for (int i = 0; i < m_config.nodes; ++i) {
    allClusterWeights[i] = allClusterWeights[i] * 100 / mass;
    if (allClusterWeights[i] > allClusterWeights[maxIndex])
      maxIndex = i;
  }

  int rank = 0;
  for (int i = 0; i < m_config.nodes; ++i) {
    if (allClusterWeights[realSource] < allClusterWeights[i])
      rank++;
  }

  pair<int, vector<double>> result;
  result.first = maxIndex; // identified solution
  result.second.push_back(allClusterWeights[maxIndex]); // solution confidence
  result.second.push_back(result.second[0] - allClusterWeights[realSource]);
  result.second.push_back(rank);

  return result;
}

void MasterNode::processResults(vector<result_t> *results,
    SimConfig& initialConfig)
{
  ofstream dumpStream;
  dumpStream.open(initialConfig.logfile.CStr());
  dumpStream << fixed << std::setprecision(3);

  for (int step = 0; step < initialConfig.steps; ++step, ++initialConfig) {
    int vars = results[step][0].second.size();
    double identificationCount = 0.0;
    double averages[vars];

    memset(averages, 0, vars * sizeof(*averages));
    for (int truth = 0; truth < initialConfig.groundTruths; ++truth) {
      identificationCount += results[step][truth].first;
      for (int var = 0; var < vars; ++var) {
        averages[var] += results[step][truth].second[var];
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
        double tmp = results[step][truth].second[var] - averages[var];
        stderrs[var] += tmp * tmp / squaredTruths;
      }
    }

    for (int var = 0; var < vars; ++var)
      stderrs[var] = sqrt(stderrs[var]);

    // Dump to output streams.
    cout << initialConfig.clusterSize << "\t";
    dumpStream << initialConfig.clusterSize << "\t";
    for (int var = 0; var < vars; ++var) {
      dumpStream << averages[var] << "\t" << stderrs[var] << "\t";
      cout << averages[var] << "\t" << stderrs[var] << "\t";
    }
    dumpStream << identificationCount << endl;
    cout << identificationCount << endl;
  }

  dumpStream.close();
}
