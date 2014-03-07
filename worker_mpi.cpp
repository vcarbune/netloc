#include "worker_mpi.h"

#include <cstring>
#include <iostream>

#undef max
#undef min

void computeCurrentWeight(const vector<GraphHypothesisCluster>& clusters,
    const SimConfig& config)
{
  double currentMassDiagonal = 0.0;
  double currentMass = 0.0;

  for (const GraphHypothesisCluster& cluster : clusters) {
    double crtWeight = cluster.getWeight();
    currentMass += crtWeight;
    currentMassDiagonal += crtWeight * crtWeight;
  }

  MPI::COMM_WORLD.Reduce(&currentMass, NULL, 1, MPI::DOUBLE,
      MPI::SUM, MPI_MASTER);
  MPI::COMM_WORLD.Reduce(&currentMassDiagonal, NULL, 1, MPI::DOUBLE,
      MPI::SUM, MPI_MASTER);
}

void reducePartialTestScores(
    const vector<GraphHypothesisCluster>& clusters,
    const SimConfig &config) {
  // Distribute partial sum for the mass.
  computeCurrentWeight(clusters, config);

  // Recompute requests coming from master.
  int currentTestNode = 0;
  MPI::COMM_WORLD.Bcast(&currentTestNode, 1, MPI::INT, MPI_MASTER);

  double massBuffer[config.objSums];
  while (currentTestNode != -1) {
    memset(massBuffer, 0, config.objSums * sizeof(*massBuffer));

    // Compute partial mass scores.
    GraphTest test(currentTestNode);
    for (const GraphHypothesisCluster& cluster : clusters) {
      pair<double, double> mass = cluster.computeMassWithTest(test);
      massBuffer[POSITIVE_SUM] += mass.first;
      massBuffer[NEGATIVE_SUM] += mass.second;
      if (config.objType == 0) { // Only for EC2
        massBuffer[POSITIVE_DIAG_SUM] += mass.first * mass.first;
        massBuffer[NEGATIVE_DIAG_SUM] += mass.second * mass.second;
      }
      massBuffer[CONS_HYPO_SUM] += cluster.getNodeCount(currentTestNode);
    }

    // Return the information to the master node.
    MPI::COMM_WORLD.Reduce(massBuffer, NULL, config.objSums,
        MPI::DOUBLE, MPI::SUM, MPI_MASTER);
    MPI::COMM_WORLD.Bcast(&currentTestNode, 1, MPI::INT, MPI_MASTER);
  }
}

void buildTestHeap(
    const vector<GraphHypothesisCluster>& clusters,
    const SimConfig& config)
{
  double massBuffer[config.objSums][config.nodes];
  memset(massBuffer, 0, config.objSums * config.nodes * sizeof(**massBuffer));

  // Compute the positive & negative mass for the remaining testNodes.
  for (int testNode = 0; testNode < config.nodes; ++testNode) {
    for (const GraphHypothesisCluster& cluster : clusters) {
      GraphTest test(testNode);
      pair<double, double> mass = cluster.computeMassWithTest(test);

      massBuffer[POSITIVE_SUM][testNode] += mass.first;
      massBuffer[NEGATIVE_SUM][testNode] += mass.second;
      if (config.objType == 0) { // Only for EC2
        massBuffer[POSITIVE_DIAG_SUM][testNode] += mass.first * mass.first;
        massBuffer[NEGATIVE_DIAG_SUM][testNode] += mass.second * mass.second;
      }
      massBuffer[CONS_HYPO_SUM][testNode] += cluster.getNodeCount(testNode);
    }
  }

  // Send to the master node the computed masses.
  for (int s = 0; s < config.objSums; ++s)
    MPI::COMM_WORLD.Reduce(massBuffer[s], NULL, config.nodes,
        MPI::DOUBLE, MPI::SUM, MPI_MASTER);

  computeCurrentWeight(clusters, config);
}

void simulate(
    vector<GraphHypothesisCluster>& clusters,
    const SimConfig& config)
{
  int totalTests = config.testThreshold * config.nodes;
  for (int count = 0; count < totalTests; ++count) {
    if (config.objType == 3) {
      // NOTHING; test is selected at random!
    } else {
      reducePartialTestScores(clusters, config);
    }

    // Receive from the master node the testNode that was selected to run.
    int selectedNode;
    bool outcome;
    int infectionTime;

    MPI::COMM_WORLD.Bcast(&selectedNode, 1, MPI::INT, MPI_MASTER);
    MPI::COMM_WORLD.Bcast(&outcome, 1, MPI::BOOL, MPI_MASTER);
    MPI::COMM_WORLD.Bcast(&infectionTime, 1, MPI::INT, MPI_MASTER);

    // Update the inner hypothesis weights, by running the testNode.
    GraphTest test(selectedNode);
    test.setOutcome(outcome);
    test.setInfectionTime(infectionTime);

    for (GraphHypothesisCluster& cluster : clusters)
      cluster.updateMassWithTest(test);
  }

  // Send the cluster masses to the central node.
  int clusterNodes[clusters.size()];
  double clusterWeight[clusters.size()];
  for (size_t node = 0; node < clusters.size(); ++node) {
    clusterNodes[node] = clusters[node].getSource();
    clusterWeight[node] = clusters[node].getWeight();
  }

  int nodes = clusters.size();
  MPI::COMM_WORLD.Send(&nodes, 1, MPI::INT, MPI_MASTER, 0);
  MPI::COMM_WORLD.Send(&clusterNodes, nodes, MPI::INT, MPI_MASTER, 0);
  MPI::COMM_WORLD.Send(&clusterWeight, nodes, MPI::DOUBLE, MPI_MASTER, 0);
}

void startWorker(PUNGraph network, SimConfig config)
{
  // Determine the node clusters for which this node is responsible.
  int clustersPerNode = network->GetNodes() / (config.mpi.nodes - 1);
  int startNode = (config.mpi.rank-1) * clustersPerNode;
  int endNode = network->GetNodes();
  if (config.mpi.rank < config.mpi.nodes - 1)
    endNode = config.mpi.rank * clustersPerNode;

  int mpiClusterStartNode = startNode;
  int mpiClusterEndNode = endNode;

  for (int step = 0; step < config.steps; ++step, ++config) {
    // Generate clusters.
    vector<GraphHypothesisCluster> clusters;
    for (int source = mpiClusterStartNode; source < mpiClusterEndNode; source++)
      clusters.push_back(GraphHypothesisCluster::generateHypothesisCluster(
            network, source, 1, config.beta, config.cascadeBound,
            config.clusterSize));

    buildTestHeap(clusters, config);
    for (int truth = 0; truth < config.groundTruths; ++truth) {
      for (GraphHypothesisCluster& cluster : clusters)
        cluster.resetWeight(1);
      simulate(clusters, config);
    }
  }
}
