#include "worker_mpi.h"

#include <iostream>

void startWorker(PUNGraph network, SimConfig config)
{
  // Send message to master.
  MPI::Status status;
  MPI::COMM_WORLD.Send(&config.mpi.rank, 1, MPI::INT, MPI_MASTER, 0);

  // Determine the node clusters for which this node is responsible.
  int clustersPerNode = network->GetNodes() / (config.mpi.nodes - 1);
  int startNode = (config.mpi.rank-1) * clustersPerNode;
  int endNode = config.mpi.rank * clustersPerNode;
  if (config.mpi.rank == config.mpi.nodes - 1)
    endNode = network->GetNodes();

  // Generate clusters.
  vector<GraphHypothesisCluster> clusters;
  for (int source = startNode; source < endNode; source++)
    clusters.push_back(GraphHypothesisCluster::generateHypothesisCluster(
          network, source, 1, config.beta, config.cascadeBound,
          config.clusterSize));

  cout << config.mpi.rank << ": " <<  startNode << " -> " << endNode << endl;

  bool testNodeWasRun[config.nodes];
  for (int i = 0; i < network->GetNodes(); ++i)
    testNodeWasRun[i] = false;

  double positiveMassDiagonal[config.nodes];
  double positiveMass[config.nodes];
  double negativeMassDiagonal[config.nodes];
  double negativeMass[config.nodes];
  double currentMassDiagonal;
  double currentMass;
  int testConsistentHypothesis[config.nodes];

  int totalTests;
  MPI::COMM_WORLD.Bcast(&totalTests, 1, MPI::INT, MPI_MASTER);
  for (int count = 0; count < totalTests; ++count) {
    // Compute the positive & negative mass for the remaining testNodes.

    // TODO(vcarbune): threadify!
    for (int testNode = 0; testNode < config.nodes; ++testNode) {
      positiveMassDiagonal[testNode] = 0.0;
      positiveMass[testNode] = 0.0;
      negativeMassDiagonal[testNode] = 0.0;
      negativeMass[testNode] = 0.0;
      currentMassDiagonal = 0.0;
      currentMass = 0.0;
      testConsistentHypothesis[testNode] = 0;

      if (testNodeWasRun[testNode])
        continue;

      for (const GraphHypothesisCluster& cluster : clusters) {
        GraphTest test(testNode);
        pair<double, double> mass = cluster.computeMassWithTest(test);

        positiveMass[testNode] += mass.first;
        negativeMass[testNode] += mass.second;

        positiveMassDiagonal[testNode] += mass.first * mass.first;
        negativeMassDiagonal[testNode] += mass.second * mass.second;

        double crtWeight = cluster.getWeight();
        currentMass += crtWeight;
        currentMassDiagonal += crtWeight * crtWeight;
        testConsistentHypothesis[testNode] += cluster.getNodeCount(testNode);
      }
    }

    // Send to the master node the computed masses.
    MPI::COMM_WORLD.Send(&positiveMassDiagonal, config.nodes, MPI::DOUBLE, MPI_MASTER, 0);
    MPI::COMM_WORLD.Send(&positiveMass, config.nodes, MPI::DOUBLE, MPI_MASTER, 0);
    MPI::COMM_WORLD.Send(&negativeMassDiagonal, config.nodes, MPI::DOUBLE, MPI_MASTER, 0);
    MPI::COMM_WORLD.Send(&negativeMass, config.nodes, MPI::DOUBLE, MPI_MASTER, 0);
    MPI::COMM_WORLD.Send(&currentMassDiagonal, 1, MPI::DOUBLE, MPI_MASTER, 0);
    MPI::COMM_WORLD.Send(&currentMass, 1, MPI::DOUBLE, MPI_MASTER, 0);
    MPI::COMM_WORLD.Send(&testConsistentHypothesis, config.nodes, MPI::INT, MPI_MASTER, 0);

    // Receive from the master node the testNode that was selected to run.
    int selectedNode;
    bool outcome;
    int infectionTime;

    MPI::COMM_WORLD.Bcast(&selectedNode, 1, MPI::INT, MPI_MASTER);
    if (selectedNode > config.nodes)
      return;

    MPI::COMM_WORLD.Bcast(&outcome, 1, MPI::BOOL, MPI_MASTER);
    MPI::COMM_WORLD.Bcast(&infectionTime, 1, MPI::INT, MPI_MASTER);

    testNodeWasRun[selectedNode] = true;

    // Update the inner hypothesis weights, by running the testNode.
    GraphTest test(selectedNode);
    test.setOutcome(outcome);
    test.setInfectionTime(infectionTime);

    for (GraphHypothesisCluster& cluster : clusters)
      cluster.updateMassWithTest(test);
  }
}
