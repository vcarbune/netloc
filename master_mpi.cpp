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

void startMaster(PUNGraph network, SimConfig config)
{
  time_t startTime = time(NULL);

  vector<GraphHypothesis> realizations;
  for (int truth = 0; truth < config.groundTruths; ++truth) {
    int sourceNode = rand() % network->GetNodes();
    realizations.push_back(GraphHypothesis::generateHypothesis(
          network, sourceNode, config.cascadeBound, config.beta));
  }

  bool testWasUsed[config.nodes];
  for (int i = 0; i < config.nodes; ++i)
    testWasUsed[i] = false;

  // Final masses, summed from what was received from each node.
  double buffer[config.nodes];
  double sums[config.nodes][EC2_SUMS];
  double crtSum[2];

  MPI::Status status;

  // Request scores for each of the tests.
  int totalTests = config.testThreshold * config.nodes;
  for (int count = 0; count < totalTests; ++count) {
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

    // Compute the aggregated score of the tests and select the lowest one.
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

      // Select the test with maximum score (that is, maximum between
      // current mass and the expected mass after the test is run).
      if (score > maxTestScore) {
        maxTestScore = score;
        maxTestNode = test;
      }
    }

    cout << count << ". " << maxTestNode << " - " << maxTestScore << endl;
    testWasUsed[maxTestNode] = true;

    // Broadcast the selected test to all the workers.
    GraphTest test(maxTestNode);
    bool outcome = realizations[0].getTestOutcome(test);
    int infection = realizations[0].getInfectionTime(maxTestNode);

    MPI::COMM_WORLD.Bcast(&maxTestNode, 1, MPI::INT, MPI_MASTER);
    MPI::COMM_WORLD.Bcast(&outcome, 1, MPI::BOOL, MPI_MASTER);
    MPI::COMM_WORLD.Bcast(&infection, 1, MPI::INT, MPI_MASTER);
  }

  // Gather information from all processes.
  double maxClusterMass[config.mpi.nodes];
  double totalMass[config.mpi.nodes];
  int sourceNodes[config.mpi.nodes];
  MPI::COMM_WORLD.Gather(&totalMass, 1, MPI::DOUBLE,
      &totalMass, 1, MPI::DOUBLE, MPI_MASTER);
  MPI::COMM_WORLD.Gather(&maxClusterMass, 1, MPI::DOUBLE,
      &maxClusterMass, 1, MPI::DOUBLE, MPI_MASTER);
  MPI::COMM_WORLD.Gather(&sourceNodes, 1, MPI::INT,
      &sourceNodes, 1, MPI::INT, MPI_MASTER);

  currentMass = 0.0;
  double maxPb = numeric_limits<double>::min();
  int source = -1;
  for (int i = 1; i < config.mpi.nodes; ++i) {
    currentMass += totalMass[i];
    if (maxPb < maxClusterMass[i]) {
      maxPb = maxClusterMass[i];
      source = sourceNodes[i];
    }
  }
  maxPb /= currentMass;

  cout << "True source: " << realizations[0].getSource() << endl;
  cout << "Found source: " << source << "(" << maxPb << ")" << endl;
  cout << "Time: " << difftime(time(NULL), startTime) << "s" << endl;
}
