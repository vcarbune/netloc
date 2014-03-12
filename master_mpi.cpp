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

double computeCurrentMass(const SimConfig& config) {
  // Compute current mass of all the clusters.
  double currentWeight = 0.0;
  double crtSum[2];
  crtSum[0] = 0.0; crtSum[1] = 0.0;
  MPI::COMM_WORLD.Reduce(&currentWeight, &crtSum[0], 1,
      MPI::DOUBLE, MPI::SUM, MPI_MASTER);
  MPI::COMM_WORLD.Reduce(&currentWeight, &crtSum[1], 1,
      MPI::DOUBLE, MPI::SUM, MPI_MASTER);

  currentWeight = crtSum[0];
  if (config.objType == 0) { // EC2
    currentWeight = crtSum[0] * crtSum[0] - crtSum[1];
  }

  return currentWeight;
}

void recomputeTestScore(
    GraphTest& test, const SimConfig& config, double currentWeight)
{
  int currentTestNode = test.getNodeId();
  MPI::COMM_WORLD.Bcast(&currentTestNode, 1, MPI::INT, MPI_MASTER);

  double sums[config.objSums];
  double junk[config.objSums];
  memset(junk, 0, config.objSums * sizeof(*junk));

  MPI::COMM_WORLD.Reduce(&junk, sums, config.objSums,
      MPI::DOUBLE, MPI::SUM, MPI_MASTER);

  double totalHypothesis = config.nodes * config.clusterSize;
  double testPositivePb = (double) sums[CONS_HYPO_SUM] / totalHypothesis;
  double score = 0.0;
  if (config.objType == 0) { // EC2
    double positiveMass = sums[POSITIVE_SUM] * sums[POSITIVE_SUM] -
        sums[POSITIVE_DIAG_SUM];
    double negativeMass = sums[NEGATIVE_SUM] * sums[NEGATIVE_SUM] -
        sums[NEGATIVE_DIAG_SUM];
    score = currentWeight -
        (testPositivePb * positiveMass + (1 - testPositivePb) * negativeMass);
  } else if (config.objType == 1) { // GBS
    score = testPositivePb * sums[POSITIVE_SUM] +
      (1 - testPositivePb) * sums[NEGATIVE_SUM];
    score = currentWeight - score;
  } else if (config.objType == 2) { // VOI
    cout << "Configuration NOT supported!" << endl;
  } else if (config.objType == 3) { // RANDOM
    cout << "Code path should not get here!" << endl;
  }

  test.setScore(score);
}

GraphTest selectNextTest(
    vector<GraphTest>& tests, const SimConfig& config) {
  double currentWeight = computeCurrentMass(config);
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
    recomputeTestScore(crtTop, config, currentWeight);

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

vector<GraphTest> buildTestHeap(const SimConfig& config) {
  double junk[config.nodes];
  double sums[config.objSums][config.nodes];

  MPI::Status status;
  memset(sums, 0, config.objSums * config.nodes * sizeof(**sums));
  memset(junk, 0, config.nodes * sizeof(*junk));

#if DBG
  time_t sumTime = time(NULL);
#endif

  // Get from workers positive & negative mass of tests still in the loop.
  for (int s = 0; s < config.objSums; ++s)
    MPI::COMM_WORLD.Reduce(&junk, sums[s], config.nodes,
        MPI::DOUBLE, MPI::SUM, MPI_MASTER);

  // Compute current mass of all the clusters.
  double currentWeight = computeCurrentMass(config);

#if DBG
  cout << "Current Weight: " << currentWeight << endl;
  cout << "Partial scores: " << difftime(time(NULL), sumTime) << "s" << endl;
#endif

  vector<GraphTest> tests;
  int totalHypothesis = config.nodes * config.clusterSize;
  for (int test = 0; test < config.nodes; ++test) {
    double score = 0.0;
    double testPositivePb = (double) sums[CONS_HYPO_SUM][test] / totalHypothesis;

    if (config.objType == 0) { // EC2
      double positiveMass = sums[POSITIVE_SUM][test] * sums[POSITIVE_SUM][test] -
          sums[POSITIVE_DIAG_SUM][test];
      double negativeMass = sums[NEGATIVE_SUM][test] * sums[NEGATIVE_SUM][test] -
          sums[NEGATIVE_DIAG_SUM][test];
      score = currentWeight -
          (testPositivePb * positiveMass + (1 - testPositivePb) * negativeMass);
    } else if (config.objType == 1) { // GBS
      score = testPositivePb * sums[POSITIVE_SUM][test] +
        (1 - testPositivePb) * sums[NEGATIVE_SUM][test];
      score = currentWeight - score;
    } else if (config.objType == 2) { // VOI
      cout << "Configuration NOT supported!" << endl;
    } else if (config.objType == 3) { // RANDOM
      cout << "Code path should not get here!" << endl;
    }

    GraphTest graphTest(test);
    graphTest.setScore(score);

    tests.push_back(graphTest);
  }

#if DBG
  cout << "Built test heap ready" << endl;
#endif
  std::make_heap<vector<GraphTest>::iterator, TestCompareFunction>(
      tests.begin(), tests.end(), TestCompareFunction());
  return tests;
}

pair<int, vector<double>> identifyClusterUsingParallelComputations(
    int realSource, const SimConfig& config)
{
  // Gather information from all processes.
  int clusterNodes[config.nodes];
  double clusterWeight[config.nodes];
  double allClusterWeights[config.nodes];

  int nodes;
  double mass = 0.0;

  MPI::Status status;
  for (int worker = 1; worker < config.mpi.nodes; ++worker) {
    MPI::COMM_WORLD.Recv(&nodes, 1, MPI::INT, worker, 0, status);
    MPI::COMM_WORLD.Recv(&clusterNodes, nodes, MPI::INT, worker, 0, status);
    MPI::COMM_WORLD.Recv(&clusterWeight, nodes, MPI::DOUBLE, worker, 0, status);
    for (int i = 0; i < nodes; ++i) {
      mass += clusterWeight[i];
      allClusterWeights[clusterNodes[i]] = clusterWeight[i];
    }
  }

  int maxIndex = 0;
  for (int i = 0; i < config.nodes; ++i) {
    allClusterWeights[i] = allClusterWeights[i] * 100 / mass;
    if (allClusterWeights[i] > allClusterWeights[maxIndex])
      maxIndex = i;
  }

  int rank = 0;
  for (int i = 0; i < config.nodes; ++i) {
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

pair<int, vector<double>> simulate(
    vector<GraphTest>& tests,
    const GraphHypothesis& realization,
    const SimConfig& config)
{
  // Request scores for each of the tests.
  int totalTests = config.testThreshold * config.nodes;
  for (int count = 0; count < totalTests; ++count) {
    GraphTest nextTest = selectNextTest(tests, config);
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
  pair<int, vector<double>> solution =
      identifyClusterUsingParallelComputations(realization.getSource(), config);

  return solution;
}

void processResults(vector<pair<int, vector<double>>> *results,
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
    double tmp;

    memset(stderrs, 0, vars * sizeof(*stderrs));
    for (int truth = 0; truth < initialConfig.groundTruths; ++truth) {
      for (int var = 0; var < vars; ++var) {
        tmp = results[step][truth].second[var] - averages[var];
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

void initializeGroundTruths(PUNGraph network,
    vector<GraphHypothesis>& realizations,
    const SimConfig& config)
{
  if (!config.groundTruth.Empty()) {
    cout << config.mpi.rank << ": Reading ground truth from " <<
      config.groundTruth.CStr() << endl;
    realizations.push_back(GraphHypothesis::readHypothesisFromFile(
          config.groundTruth.CStr()));
    return;
  }

  for (int truth = 0; truth < config.groundTruths; ++truth) {
    realizations.push_back(GraphHypothesis::generateHypothesis(
          network, rand() % network->GetNodes(),
          config.cascadeBound, config.beta));
  }
}

void startMaster(PUNGraph network, SimConfig config)
{
  time_t startTime = time(NULL);
  cout << fixed << std::setprecision(3);

  vector<GraphHypothesis> realizations;
  initializeGroundTruths(network, realizations, config);

  SimConfig initialConfig = config;
  vector<pair<int, vector<double>>> results[config.steps];
  for (int step = 0; step < config.steps; ++step, ++config) {
    double startTime = time(NULL);

    vector<GraphTest> tests = buildTestHeap(config);
    for (int truth = 0; truth < config.groundTruths; ++truth) {
      vector<GraphTest> tempTests(tests);
      pair<int, vector<double>> result =
          simulate(tempTests, realizations[truth], config);

      // Add distance (in hops) to the source node.
      result.second.push_back(TSnap::GetShortPath(network,
            realizations[truth].getSource(), result.first));

      // Keep a boolean, whether the source was identified or not.
      result.first = realizations[truth].getSource() == result.first ? 1 : 0;

      // Store all the data to process the results later.
      results[step].push_back(result);

      cout << config << "-" << truth << "\t";
      for (size_t i = 0; i < result.second.size(); ++i)
        cout << result.second[i] << "\t";
      cout << endl;
    }
    cout << "Time: " << difftime(time(NULL), startTime) << "s\n\n";
  }

  processResults(results, initialConfig);
  cout << "Time: " << difftime(time(NULL), startTime) << "s" << endl;
}
