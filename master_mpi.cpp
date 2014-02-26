#include "master_mpi.h"

#include <iostream>
#include <limits>

#undef max
#undef min

void startMaster(PUNGraph network, SimConfig config)
{
  // Gather message from workers.
  MPI::Status status;
  for (int worker = 1; worker < config.mpi.nodes; worker++) {
    int msg;
    MPI::COMM_WORLD.Recv(&msg, 1, MPI::INT, worker, 0, status);
    cout << "Received from " << worker << ": " << msg << endl;
  }

  time_t startTime = time(NULL);

  GraphHypothesis realization = GraphHypothesis::generateHypothesis(network,
          rand() % network->GetNodes(), config.cascadeBound, config.beta);

  bool testNodeWasRun[config.nodes];
  for (int i = 0; i < config.nodes; ++i)
    testNodeWasRun[i] = false;

  // Final masses, summed from what was received from each node.
  double positiveMassDiagonal[config.nodes],
         tempPositiveMassDiagonal[config.nodes];
  double positiveMass[config.nodes],
         tempPositiveMass[config.nodes];
  double negativeMassDiagonal[config.nodes],
         tempNegativeMassDiagonal[config.nodes];
  double negativeMass[config.nodes],
         tempNegativeMass[config.nodes];
  int testConsistentHypothesis[config.nodes],
      tempTestConsistentHypothesis[config.nodes];
  double currentMassDiagonal = 0.0;
  double tempCurrentMassDiagonal = 0.0;
  double currentMass = 0.0;
  double tempCurrentMass = 0.0;
  double finalScores[config.nodes];

  // Request scores for each of the tests.
  int totalTests = config.testThreshold * config.nodes;
  // MPI::COMM_WORLD.Bcast(&totalTests, 1, MPI::INT, MPI_MASTER);
  for (int count = 0; count < totalTests; ++count) {
    currentMass = 0.0;
    currentMassDiagonal = 0.0;

    // Reset total scores to zero
    for (int test = 0; test < config.nodes; ++test) {
      positiveMassDiagonal[test] = 0.0;
      positiveMass[test] = 0.0;
      negativeMassDiagonal[test] = 0.0;
      negativeMass[test] = 0.0;
      testConsistentHypothesis[test] = 0;
    }

    // Get from workers positive & negative mass of tests still in the loop.
    for (int worker = 1; worker < config.mpi.nodes; worker++) {
      MPI::COMM_WORLD.Recv(&tempPositiveMassDiagonal, config.nodes, MPI::DOUBLE, worker, 0, status);
      MPI::COMM_WORLD.Recv(&tempPositiveMass, config.nodes, MPI::DOUBLE, worker, 0, status);
      MPI::COMM_WORLD.Recv(&tempNegativeMassDiagonal, config.nodes, MPI::DOUBLE, worker, 0, status);
      MPI::COMM_WORLD.Recv(&tempNegativeMass, config.nodes, MPI::DOUBLE, worker, 0, status);
      MPI::COMM_WORLD.Recv(&tempCurrentMassDiagonal, 1, MPI::DOUBLE, worker, 0, status);
      MPI::COMM_WORLD.Recv(&tempCurrentMass, 1, MPI::DOUBLE, worker, 0, status);
      MPI::COMM_WORLD.Recv(&tempTestConsistentHypothesis, config.nodes, MPI::INT, worker, 0, status);

      currentMassDiagonal += tempCurrentMassDiagonal;
      currentMass += tempCurrentMass;

      for (int test = 0; test < config.nodes; ++test) {
        if (testNodeWasRun[test])
          continue;

        positiveMassDiagonal[test] += tempPositiveMassDiagonal[test];
        positiveMass[test] += tempPositiveMass[test];
        negativeMassDiagonal[test] += tempNegativeMassDiagonal[test];
        negativeMass[test] += tempNegativeMass[test];
        testConsistentHypothesis[test] += tempTestConsistentHypothesis[test];
      }
   }

    double currentWeight = currentMass * currentMass - currentMassDiagonal;

    // Compute the aggregated score of the tests and select the lowest one.
    double maxTestScore = numeric_limits<double>::min();
    int maxTestNode = -1;

    int totalHypothesis = config.nodes * config.clusterSize;
    for (int test = 0; test < config.nodes; ++test) {
      if (testNodeWasRun[test])
        continue;

      positiveMass[test] = positiveMass[test] * positiveMass[test] -
          positiveMassDiagonal[test];
      negativeMass[test] = negativeMass[test] * negativeMass[test] -
          negativeMassDiagonal[test];

      double testPositivePb =
          (double) testConsistentHypothesis[test] / totalHypothesis;

      finalScores[test] = testPositivePb * positiveMass[test] +
          (1 - testPositivePb) * negativeMass[test];
      finalScores[test] = currentWeight - finalScores[test];

      // Select the test with maximum score (that is, maximum between
      // current mass and the expected mass after the test is run).
      if (finalScores[test] > maxTestScore || maxTestNode < 0) {
        maxTestScore = finalScores[test];
        maxTestNode = test;
      }
    }

    cout << count << ". " << maxTestNode << " - " << maxTestScore << endl;
    testNodeWasRun[maxTestNode] = true;

    // Broadcast the selected test to all the workers.
    GraphTest test(maxTestNode);
    bool outcome = realization.getTestOutcome(test);
    int infection = realization.getInfectionTime(maxTestNode);

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

  cout << "True source: " << realization.getSource() << endl;
  cout << "Found source: " << source << "(" << maxPb << ")" << endl;
  cout << "Time: " << difftime(time(NULL), startTime) << "s" << endl;
}
