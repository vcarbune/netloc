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

#define POSITIVE_SUM      0
#define NEGATIVE_SUM      1
#define CONS_HYPO_SUM     2
#define POSITIVE_DIAG_SUM 3
#define NEGATIVE_DIAG_SUM 5

#define EC2_SUMS 5

pair<int, double> selectNextTest(bool *testWasUsed, const SimConfig& config) {
  // Final masses, summed from what was received from each node.
  double junk[config.nodes];
  double sums[EC2_SUMS][config.nodes];
  double crtSum[2];

  MPI::Status status;
  memset(sums, 0, EC2_SUMS * config.nodes * sizeof(**sums));
  memset(junk, 0, config.nodes * sizeof(*junk));

#if DBG
  time_t sumTime = time(NULL);
#endif

  // Get from workers positive & negative mass of tests still in the loop.
  for (int s = 0; s < EC2_SUMS; ++s) {
    MPI::COMM_WORLD.Reduce(&junk, sums[s], config.nodes,
        MPI::DOUBLE, MPI::SUM, MPI_MASTER);
  }

  // Compute current mass of all the clusters.
  double currentWeight = 0.0;
  MPI::COMM_WORLD.Reduce(&currentWeight, &crtSum[0], 1,
      MPI::DOUBLE, MPI::SUM, MPI_MASTER);
  MPI::COMM_WORLD.Reduce(&currentWeight, &crtSum[1], 1,
      MPI::DOUBLE, MPI::SUM, MPI_MASTER);
  currentWeight = crtSum[0] * crtSum[0] - crtSum[1];

#if DBG
  cout << "Partial scores: " << difftime(time(NULL), sumTime) << "s" << endl;
#endif

  // Generalize to non-EC2 objective functions.
  double maxTestScore = numeric_limits<double>::min();
  int maxTestNode = -1;

  int totalHypothesis = config.nodes * config.clusterSize;
  for (int test = 0; test < config.nodes; ++test) {
    if (testWasUsed[test])
      continue;

    double positiveMass = sums[POSITIVE_SUM][test] * sums[POSITIVE_SUM][test] -
        sums[POSITIVE_DIAG_SUM][test];
    double negativeMass = sums[NEGATIVE_SUM][test] * sums[NEGATIVE_SUM][test] -
        sums[NEGATIVE_DIAG_SUM][test];
    double testPositivePb =
        (double) sums[CONS_HYPO_SUM][test] / totalHypothesis;
    double score = currentWeight -
        (testPositivePb * positiveMass + (1 - testPositivePb) * negativeMass);

    if (score > maxTestScore) {
      maxTestScore = score;
      maxTestNode = test;
    }
  }

  return pair<int, double>(maxTestNode, maxTestScore);
}

pair<int, pair<double, double>> identifyCluster(int realSource, const SimConfig& config)
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

  pair<int, pair<double, double>> result;
  result.first = sourceNodes[maxIndex];               // identified solution
  result.second.first = 100 * realSourceMass / mass;  // solution confidence
  result.second.second = 100 * maxMass[maxIndex] / mass - result.second.first;

  return result;
}

pair<int, pair<double, double>> simulate(
    const GraphHypothesis& realization,
    const SimConfig& config)
{
  bool testWasUsed[config.nodes];
  for (int i = 0; i < config.nodes; ++i)
    testWasUsed[i] = false;

  // Request scores for each of the tests.
  int totalTests = config.testThreshold * config.nodes;
  for (int count = 0; count < totalTests; ++count) {
    pair<int, double> nextTest = selectNextTest(testWasUsed, config);

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
  pair<int, pair<double, double>> solution =
      identifyCluster(realization.getSource(), config);

  return solution;
}

void processResults(vector<pair<int, pair<double, double>>> *results,
    SimConfig& initialConfig)
{
  ofstream dumpStream;
  dumpStream.open(initialConfig.logfile.CStr());
  dumpStream << fixed << std::setprecision(3);

  for (int step = 0; step < initialConfig.steps; ++step, ++initialConfig) {
    double averageMassBuffer = 0.0;
    double averageDiffBuffer = 0.0;
    double identificationCount = 0.0;

    // Compute averages and identification count.
    for (int truth = 0; truth < initialConfig.groundTruths; ++truth) {
      averageMassBuffer += results[step][truth].second.first;
      averageDiffBuffer += results[step][truth].second.second;
      identificationCount += results[step][truth].first;
    }
    averageMassBuffer /= initialConfig.groundTruths;
    averageDiffBuffer /= initialConfig.groundTruths;
    identificationCount /= initialConfig.groundTruths;


    // Compute standard error.
    double stderrMass = 0.0;
    double stderrDiff = 0.0;
    double squaredTruths =
        initialConfig.groundTruths * initialConfig.groundTruths;

    double tmp;
    for (int truth = 0; truth < initialConfig.groundTruths; ++truth) {
      tmp = results[step][truth].second.first - averageMassBuffer;
      stderrMass += tmp * tmp / squaredTruths;
      tmp = results[step][truth].second.second - averageDiffBuffer;
      stderrDiff += tmp * tmp / squaredTruths;
    }

    stderrMass = sqrt(stderrMass);
    stderrDiff = sqrt(stderrDiff);

    // Dump to output streams.
    dumpStream << initialConfig.clusterSize << "\t" <<
      averageMassBuffer << "\t" << stderrMass << "\t" <<
      averageDiffBuffer << "\t" << stderrDiff << "\t" <<
      identificationCount << endl;

    // Dump to cout too, to be included in the mail.
    cout << initialConfig.clusterSize << "\t" <<
      averageMassBuffer << "\t" << stderrMass << "\t" <<
      averageDiffBuffer << "\t" << stderrDiff << "\t" <<
      identificationCount << endl;
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
  vector<pair<int, pair<double, double>>> results[config.steps];
  for (int step = 0; step < config.steps; ++step, ++config) {
    double startTime = time(NULL);
    for (int truth = 0; truth < config.groundTruths; ++truth) {
      pair<int, pair<double, double>> result =
          simulate(realizations[truth], config);

      // Keep a boolean, whether the source was identified or not.
      result.first = realizations[truth].getSource() == result.first ? 1 : 0;

      // Store all the data to process the results later.
      results[step].push_back(result);

      cout << config << "-" << truth << "\t" <<
        result.second.first << "\t" << result.second.second << "\t" << endl;
    }
    cout << "Time: " << difftime(time(NULL), startTime) << "s\n\n";
  }

  processResults(results, initialConfig);
  cout << "Time: " << difftime(time(NULL), startTime) << "s" << endl;
}