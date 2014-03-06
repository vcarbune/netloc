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

pair <int, double> selectNextTestUsingParallelLazyGreedy(
    bool *testWasUsed, const SimConfig& config)
{
  // Gather test nodes and scores from mpi nodes.
  int testNodeId[config.mpi.nodes];
  double testScore[config.mpi.nodes];

  int junk = -1;
  double score = 0.0;

  MPI::COMM_WORLD.Gather(&junk, 1, MPI::INT, &testNodeId, 1, MPI::INT, MPI_MASTER);
  MPI::COMM_WORLD.Gather(&score, 1, MPI::DOUBLE, &testScore, 1, MPI::DOUBLE, MPI_MASTER);

  double maxTestScore = numeric_limits<double>::min();
  int maxTestNode = -1;

  for (int i = 1; i < config.mpi.nodes; ++i) {
#if DBG
    cout << testNodeId[i] << " - " << testScore[i] << endl;
#endif
    if (testNodeId[i] != -1 && testScore[i] > maxTestScore) {
      maxTestNode = testNodeId[i];
      maxTestScore = testScore[i];
    }
  }

  return pair<int, double>(maxTestNode, maxTestScore);
}

pair<int, double> selectNextTestUsingParallelComputations(
    bool *testWasUsed, const SimConfig& config) {
  // Final masses, summed from what was received from each node.
  double junk[config.nodes];
  double sums[config.objSums][config.nodes];
  double crtSum[2];

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
  double currentWeight = 0.0;
  crtSum[0] = 0.0; crtSum[1] = 0.0;
  MPI::COMM_WORLD.Reduce(&currentWeight, &crtSum[0], 1,
      MPI::DOUBLE, MPI::SUM, MPI_MASTER);
  MPI::COMM_WORLD.Reduce(&currentWeight, &crtSum[1], 1,
      MPI::DOUBLE, MPI::SUM, MPI_MASTER);

  currentWeight = crtSum[0];
  if (config.objType == 0) { // EC2
    currentWeight = crtSum[0] * crtSum[0] - crtSum[1];
  }

#if DBG
  cout << "Current Weight: " << currentWeight << endl;
  cout << "Partial scores: " << difftime(time(NULL), sumTime) << "s" << endl;
#endif

  // Generalize to non-EC2 objective functions.
  double maxTestScore = numeric_limits<double>::min();
  int maxTestNode = -1;

  if (config.objType == 3) { // RANDOM
  }

  int totalHypothesis = config.nodes * config.clusterSize;
  for (int test = 0; test < config.nodes; ++test) {
    if (testWasUsed[test])
      continue;

    double score = 0.0;
    double testPositivePb =
        (double) sums[CONS_HYPO_SUM][test] / totalHypothesis;

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
    }

    if (score > maxTestScore) {
      maxTestScore = score;
      maxTestNode = test;
    }
  }

  return pair<int, double>(maxTestNode, maxTestScore);
}

pair<int, vector<double>> identifyClusterUsingLazyGreedyAndVoting(
    int realSource, const SimConfig& config)
{
  // Gather information from all processes.
  double totalMass;
  double zeroValues[config.nodes];
  double clusterMass[config.nodes];

  double junk = 0.0;
  MPI::COMM_WORLD.Reduce(&junk, &totalMass, 1, MPI::DOUBLE, MPI::SUM, MPI_MASTER);

  memset(zeroValues, 0, config.nodes * sizeof(*zeroValues));
  MPI::COMM_WORLD.Reduce(zeroValues, clusterMass, config.nodes,
      MPI::DOUBLE, MPI::SUM, MPI_MASTER);

  int maxIndex = 0;
  for (int i = 0; i < config.nodes; ++i)
    if (clusterMass[maxIndex] < clusterMass[i])
      maxIndex = i;
  double realSourceMass = clusterMass[realSource];

  pair<int, vector<double>> result;
  result.first = maxIndex; // identified solution

  // Mass in solution cluster and difference of mass.
  result.second.push_back(100 * realSourceMass);
  result.second.push_back(100 * (clusterMass[maxIndex] - realSourceMass));

  result.second[0] /= totalMass;
  result.second[1] /= totalMass;

  // Rank of the solution.
  result.second.push_back(0.0);
  for (int i = 0; i < config.nodes; ++i)
    if (realSourceMass < clusterMass[i])
      result.second[2] += 1;

  return result;
}

pair<int, vector<double>> identifyClusterUsingParallelComputations(
    int realSource, const SimConfig& config)
{
  // Gather information from all processes.
  int sourceNodes[config.mpi.nodes];
  double maxMass[config.mpi.nodes];

  double e = 0.0;
  double mass = 0.0;

  MPI::COMM_WORLD.Reduce(&e, &mass, 1, MPI::DOUBLE, MPI::SUM, MPI_MASTER);
  MPI::COMM_WORLD.Gather(&e, 1, MPI::DOUBLE, &maxMass, 1, MPI::DOUBLE, MPI_MASTER);
  MPI::COMM_WORLD.Gather(&e, 1, MPI::INT, &sourceNodes, 1, MPI::INT, MPI_MASTER);

  int maxIndex = 1;
  for (int i = 1; i < config.mpi.nodes; ++i)
    if (maxMass[maxIndex] < maxMass[i])
      maxIndex = i;

  MPI::COMM_WORLD.Bcast(&realSource, 1, MPI::INT, MPI_MASTER);

  double realSourceMass;
  MPI::Status status;
  MPI::COMM_WORLD.Recv(&realSourceMass, config.nodes, MPI::DOUBLE, MPI::ANY_SOURCE, 1, status);

  pair<int, vector<double>> result;
  result.first = sourceNodes[maxIndex]; // identified solution
  result.second.push_back(100 * realSourceMass); // solution confidence
  result.second.push_back(100 * (maxMass[maxIndex] - realSourceMass));

  result.second[0] /= mass;
  result.second[1] /= mass;

  return result;
}

pair<int, vector<double>> simulate(
    const GraphHypothesis& realization,
    const SimConfig& config)
{
  bool testWasUsed[config.nodes];
  for (int i = 0; i < config.nodes; ++i)
    testWasUsed[i] = false;

  // Request scores for each of the tests.
  int totalTests = config.testThreshold * config.nodes;
  for (int count = 0; count < totalTests; ++count) {
    pair<int, double> nextTest;
    if (config.objType == 3) {
      do {
        int maxTestNode = rand() % config.nodes;
        if (!testWasUsed[maxTestNode]) {
          nextTest = pair<int, double>(maxTestNode, 0.314159);
          break;
        }
      } while (true);
    } else if (config.lazy) {
      nextTest = selectNextTestUsingParallelLazyGreedy(testWasUsed, config);
    } else {
      nextTest = selectNextTestUsingParallelComputations(testWasUsed, config);
    }

#if DBG
    cout << count << ". " << nextTest.first << " " << nextTest.second << endl;
#endif
    testWasUsed[nextTest.first] = true;

    // Inform the workers about the selected test.
    GraphTest test(nextTest.first);
    bool outcome = realization.getTestOutcome(test);
    int infection = realization.getInfectionTime(nextTest.first);

#if DBG
    time_t bcastTime = time(NULL);
#endif
    MPI::COMM_WORLD.Bcast(&nextTest.first, 1, MPI::INT, MPI_MASTER);
    MPI::COMM_WORLD.Bcast(&outcome, 1, MPI::BOOL, MPI_MASTER);
    MPI::COMM_WORLD.Bcast(&infection, 1, MPI::INT, MPI_MASTER);
#if DBG
    cout << "Broadcasting took " << difftime(time(NULL), bcastTime) <<
        "s" << endl;
#endif
  }

  // Identify the cluster where the mass is concentrated.
  pair<int, vector<double>> solution = config.lazy ?
    identifyClusterUsingLazyGreedyAndVoting(realization.getSource(), config) :
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

void startMaster(PUNGraph network, SimConfig config)
{
  time_t startTime = time(NULL);
  cout << fixed << std::setprecision(3);

  vector<GraphHypothesis> realizations;
  for (int truth = 0; truth < config.groundTruths; ++truth) {
    realizations.push_back(GraphHypothesis::generateHypothesis(
          network, rand() % network->GetNodes(),
          config.cascadeBound, config.beta));
  }

  SimConfig initialConfig = config;
  vector<pair<int, vector<double>>> results[config.steps];
  for (int step = 0; step < config.steps; ++step, ++config) {
    double startTime = time(NULL);
    for (int truth = 0; truth < config.groundTruths; ++truth) {
      pair<int, vector<double>> result =
          simulate(realizations[truth], config);

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
