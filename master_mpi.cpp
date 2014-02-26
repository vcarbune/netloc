#include "master_mpi.h"

#include <iostream>
#include <limits>

#undef max
#undef min

#define POSITIVE_SUM      0
#define POSITIVE_DIAG_SUM 1
#define NEGATIVE_SUM      2
#define NEGATIVE_DIAG_SUM 3
#define CONS_HYPO_SUM     4

#define EC2_SUMS 5

pair<int, double> selectNextTest(bool *testWasUsed, const SimConfig& config) {
  // Final masses, summed from what was received from each node.
  double buffer[config.nodes];
  double sums[config.nodes][EC2_SUMS];
  double crtSum[2];

  MPI::Status status;
  for (int test = 0; test < config.nodes; ++test)
    for (int s = 0; s < EC2_SUMS; ++s)
      sums[test][s] = 0.0;

  // Get from workers positive & negative mass of tests still in the loop.
  for (int worker = 1; worker < config.mpi.nodes; worker++) {
    for (int s = 0; s < EC2_SUMS; ++s) {
      MPI::COMM_WORLD.Recv(buffer, config.nodes, MPI::DOUBLE, worker, 0, status);
      for (int test = 0; test < config.nodes; ++test)
        sums[test][s] += buffer[test];
    }
  }

  // Compute current mass of all the clusters.
  double currentWeight = 0.0;
  MPI::COMM_WORLD.Reduce(&currentWeight, &crtSum[0], 1,
      MPI::DOUBLE, MPI::SUM, MPI_MASTER);
  MPI::COMM_WORLD.Reduce(&currentWeight, &crtSum[1], 1,
      MPI::DOUBLE, MPI::SUM, MPI_MASTER);
  currentWeight = crtSum[0] * crtSum[0] - crtSum[1];

  // Generalize to non-EC2 objective functions.
  double maxTestScore = numeric_limits<double>::min();
  int maxTestNode = -1;

  int totalHypothesis = config.nodes * config.clusterSize;
  for (int test = 0; test < config.nodes; ++test) {
    if (testWasUsed[test])
      continue;

    double positiveMass = sums[test][POSITIVE_SUM] * sums[test][POSITIVE_SUM] -
        sums[test][POSITIVE_DIAG_SUM];
    double negativeMass = sums[test][NEGATIVE_SUM] * sums[test][NEGATIVE_SUM] -
        sums[test][NEGATIVE_DIAG_SUM];
    double testPositivePb =
        (double) sums[test][CONS_HYPO_SUM] / totalHypothesis;
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
  result.first = sourceNodes[maxIndex];           // identified solution
  result.second.first = maxMass[maxIndex] / mass; // confidence
  result.second.second = result.second.first - realSourceMass / mass;
  return result;
}

pair<int, pair<double, double>> simulate(
    const GraphHypothesis& realization, const SimConfig& config)
{
  bool testWasUsed[config.nodes];
  for (int i = 0; i < config.nodes; ++i)
    testWasUsed[i] = false;

  // Request scores for each of the tests.
  int totalTests = config.testThreshold * config.nodes;
  for (int count = 0; count < totalTests; ++count) {
    pair<int, double> nextTest = selectNextTest(testWasUsed, config);

    cout << count << ". " << nextTest.first << " " << nextTest.second << endl;
    testWasUsed[nextTest.first] = true;

    // Inform the workers about the selected test.
    GraphTest test(nextTest.first);
    bool outcome = realization.getTestOutcome(test);
    int infection = realization.getInfectionTime(nextTest.first);

    MPI::COMM_WORLD.Bcast(&nextTest.first, 1, MPI::INT, MPI_MASTER);
    MPI::COMM_WORLD.Bcast(&outcome, 1, MPI::BOOL, MPI_MASTER);
    MPI::COMM_WORLD.Bcast(&infection, 1, MPI::INT, MPI_MASTER);
  }

  // Identify the cluster where the mass is concentrated.
  pair<int, pair<double, double>> solution =
      identifyCluster(realization.getSource(), config);

  cout << "Expected: " << realization.getSource() << endl;
  cout << "Found: " << solution.first << "(" <<
      solution.second.first << ", " << solution.second.second <<
      ")" << endl;

  return solution;
}

void startMaster(PUNGraph network, SimConfig config)
{
  time_t startTime = time(NULL);
  for (int truth = 0; truth < config.groundTruths; ++truth) {
    GraphHypothesis realization = GraphHypothesis::generateHypothesis(
        network, rand() % network->GetNodes(), config.cascadeBound, config.beta);
    pair<int, pair<double, double>> result = simulate(realization, config);
  }
  cout << "Time: " << difftime(time(NULL), startTime) << "s" << endl;
}
