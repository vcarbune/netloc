#include "worker_mpi.h"

#include <cstring>
#include <iostream>

#undef max
#undef min

void simulate(vector<GraphHypothesisCluster>& clusters,
              const SimConfig& config)
{
  bool testWasUsed[config.nodes];
  for (int i = 0; i < config.nodes; ++i)
    testWasUsed[i] = false;

  double currentMassDiagonal;
  double currentMass;

  double massBuffer[config.objSums][config.nodes];

  int totalTests = config.testThreshold * config.nodes;
  for (int count = 0; count < totalTests; ++count) {
    memset(massBuffer, 0, config.objSums * config.nodes * sizeof(**massBuffer));

    // Compute the positive & negative mass for the remaining testNodes.
    for (int testNode = 0; testNode < config.nodes; ++testNode) {
      if (testWasUsed[testNode])
        continue;

      for (const GraphHypothesisCluster& cluster : clusters) {
        GraphTest test(testNode);
        pair<double, double> mass = cluster.computeMassWithTest(test);

        massBuffer[POSITIVE_SUM][testNode] += mass.first;
        massBuffer[NEGATIVE_SUM][testNode] += mass.second;
        massBuffer[POSITIVE_DIAG_SUM][testNode] += mass.first * mass.first;
        massBuffer[NEGATIVE_DIAG_SUM][testNode] += mass.second * mass.second;
        massBuffer[CONS_HYPO_SUM][testNode] += cluster.getNodeCount(testNode);
      }
    }

    currentMassDiagonal = 0.0;
    currentMass = 0.0;

    for (const GraphHypothesisCluster& cluster : clusters) {
      double crtWeight = cluster.getWeight();
      currentMass += crtWeight;
      currentMassDiagonal += crtWeight * crtWeight;
    }

    // Send to the master node the computed masses.
    for (int s = 0; s < config.objSums; ++s)
      MPI::COMM_WORLD.Reduce(massBuffer[s], NULL, config.nodes,
          MPI::DOUBLE, MPI::SUM, MPI_MASTER);

    MPI::COMM_WORLD.Reduce(&currentMass, NULL, 1, MPI::DOUBLE,
        MPI::SUM, MPI_MASTER);
    MPI::COMM_WORLD.Reduce(&currentMassDiagonal, NULL, 1, MPI::DOUBLE,
        MPI::SUM, MPI_MASTER);

    // Receive from the master node the testNode that was selected to run.
    int selectedNode;
    bool outcome;
    int infectionTime;

    MPI::COMM_WORLD.Bcast(&selectedNode, 1, MPI::INT, MPI_MASTER);
    MPI::COMM_WORLD.Bcast(&outcome, 1, MPI::BOOL, MPI_MASTER);
    MPI::COMM_WORLD.Bcast(&infectionTime, 1, MPI::INT, MPI_MASTER);

    testWasUsed[selectedNode] = true;

    // Update the inner hypothesis weights, by running the testNode.
    GraphTest test(selectedNode);
    test.setOutcome(outcome);
    test.setInfectionTime(infectionTime);

    for (GraphHypothesisCluster& cluster : clusters)
      cluster.updateMassWithTest(test);
  }

  // Identify the cluster with the highest mass and send it back
  // to the master process to identify the highest one.
  double maxMass = numeric_limits<double>::min();
  double totalMass = 0.0;
  int sourceNode = -1;

  for (const GraphHypothesisCluster& cluster : clusters) {
    totalMass += cluster.getMass();
    if (cluster.getMass() > maxMass) {
      maxMass = cluster.getMass();
      sourceNode = cluster.getSource();
    }
  }
  MPI::COMM_WORLD.Reduce(&totalMass, NULL, 1, MPI::DOUBLE, MPI::SUM, MPI_MASTER);
  MPI::COMM_WORLD.Gather(&maxMass, 1, MPI::DOUBLE, NULL, 1, MPI::DOUBLE, MPI_MASTER);
  MPI::COMM_WORLD.Gather(&sourceNode, 1, MPI::INT, NULL, 1, MPI::INT, MPI_MASTER);
}

void startWorker(PUNGraph network, SimConfig config)
{
  // Determine the node clusters for which this node is responsible.
  int clustersPerNode = network->GetNodes() / (config.mpi.nodes - 1);
  int startNode = (config.mpi.rank-1) * clustersPerNode;
  int endNode = config.mpi.rank * clustersPerNode;
  if (config.mpi.rank == config.mpi.nodes - 1)
    endNode = network->GetNodes();

  for (int step = 0; step < config.steps; ++step, ++config) {
    // Generate clusters.
    vector<GraphHypothesisCluster> clusters;
    for (int source = startNode; source < endNode; source++)
      clusters.push_back(GraphHypothesisCluster::generateHypothesisCluster(
            network, source, 1, config.beta, config.cascadeBound,
            config.clusterSize));

    for (int truth = 0; truth < config.groundTruths; ++truth) {
      for (GraphHypothesisCluster& cluster : clusters)
        cluster.resetWeight(1);

      simulate(clusters, config);

      int realSource;
      MPI::COMM_WORLD.Bcast(&realSource, 1, MPI::INT, MPI_MASTER);
      if (realSource >= startNode && realSource < endNode) {
        double realSourceMass = clusters[realSource - startNode].getWeight();
        if (clusters[realSource - startNode].getSource() != realSource)
          cout << "!!ERROR!!" << endl;
        MPI::COMM_WORLD.Send(&realSourceMass, 1, MPI::DOUBLE, MPI_MASTER, 1);
      }
    }
  }
}
