#include "worker_mpi.h"

#include <cstring>
#include <iostream>

#undef max
#undef min

WorkerNode::WorkerNode(SimConfig config)
  : MPINode(config)
{
}

void WorkerNode::run()
{
  // No reason for a distributed version EPFL max-likelihood (distributed
  // matrix computations are a different beast, so keeping, for now, algorithm
  // implementation on single machine, running on central MPI node only)
  if (m_config.objType == EPFL_ML) {
    cout << "Note: EPFL runs only on central MPI node" << endl;
    cout << "Worker " << m_config.mpi.rank << " exiting..." << endl;
    return;
  }

  // Determine the node clusters for which this node is responsible.
  int clustersPerNode = m_network->GetNodes() / (m_config.mpi.nodes - 1);
  int startNode = (m_config.mpi.rank-1) * clustersPerNode;
  int endNode = m_network->GetNodes();
  if (m_config.mpi.rank < m_config.mpi.nodes - 1)
    endNode = m_config.mpi.rank * clustersPerNode;

  for (int step = 0; step < m_config.steps; ++step, ++m_config) {
    initializeClusters(startNode, endNode);
    initializeTestHeap();
    for (int truth = 0; truth < m_config.groundTruths; ++truth) {
      simulate();
      for (GraphHypothesisCluster& cluster : m_clusters)
        cluster.resetWeight(1);
    }
  }
}

void WorkerNode::initializeClusters(int startNode, int endNode)
{
  m_clusters.clear();
  for (int src = startNode; src < endNode; src++)
    m_clusters.push_back(GraphHypothesisCluster::generateHypothesisCluster(
          m_network, src, 1, m_config.beta, m_config.cascadeBound,
          m_config.clusterSize));
}

void WorkerNode::initializeTestHeap()
{
  if (m_config.objType == RANDOM)
    return;

  double massBuffer[m_config.objSums][m_config.nodes];
  memset(massBuffer, 0, m_config.objSums * m_config.nodes * sizeof(**massBuffer));

  // Compute partial test priors first and distribute them to master node.
  for (int testNode = 0; testNode < m_config.nodes; ++testNode)
    for (const GraphHypothesisCluster& cluster : m_clusters)
      massBuffer[0][testNode] += cluster.getNodeCount(testNode);
  MPI::COMM_WORLD.Reduce(massBuffer[0], NULL, m_config.nodes,
      MPI::DOUBLE, MPI::SUM, MPI_MASTER);

  // Get test priors computed at the central node.
  MPI::COMM_WORLD.Bcast(massBuffer[0], m_config.nodes, MPI::DOUBLE, MPI_MASTER);
  for (int testNode = 0; testNode < m_config.nodes; ++testNode)
    m_testsPrior[testNode] = massBuffer[0][testNode];

  // Compute the positive & negative mass for the remaining testNodes.
  for (int testNode = 0; testNode < m_config.nodes; ++testNode) {
    for (const GraphHypothesisCluster& cluster : m_clusters) {
      GraphTest test(testNode);
      pair<double, double> mass = cluster.computeMassWithTest(test);

      if (m_config.objType == VOI) {
        double expectedMass = m_testsPrior[testNode] * mass.first +
            (1 - m_testsPrior[testNode]) * mass.second;
        massBuffer[0][testNode] = max(massBuffer[0][testNode], expectedMass);
      } else {
        massBuffer[POSITIVE_SUM][testNode] += mass.first;
        massBuffer[NEGATIVE_SUM][testNode] += mass.second;
        if (m_config.objType == EC2) { // Only for EC2
          massBuffer[POSITIVE_DIAG_SUM][testNode] += mass.first * mass.first;
          massBuffer[NEGATIVE_DIAG_SUM][testNode] += mass.second * mass.second;
        }
      }
    }
  }

  // Send to the master node the computed masses.
  MPI::Op op = MPI::SUM;
  if (m_config.objType == VOI)
    op = MPI::MAX;
  for (int s = 0; s < m_config.objSums; ++s)
    MPI::COMM_WORLD.Reduce(massBuffer[s], NULL, m_config.nodes,
        MPI::DOUBLE, op, MPI_MASTER);
  computeCurrentMass();
}

void WorkerNode::simulate()
{
  int totalTests = m_config.testThreshold * m_config.nodes;
  for (int count = 0; count < totalTests; ++count) {
    recomputePartialTestScores();

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

    for (GraphHypothesisCluster& cluster : m_clusters)
      cluster.updateMassWithTest(test);
  }

  // Send the cluster masses to the central node.
  int clusterNodes[m_clusters.size()];
  double clusterWeight[m_clusters.size()];
  for (size_t node = 0; node < m_clusters.size(); ++node) {
    clusterNodes[node] = m_clusters[node].getSource();
    clusterWeight[node] = m_clusters[node].getWeight();
  }

  int nodes = m_clusters.size();
  MPI::COMM_WORLD.Send(&nodes, 1, MPI::INT, MPI_MASTER, 0);
  MPI::COMM_WORLD.Send(&clusterNodes, nodes, MPI::INT, MPI_MASTER, 0);
  MPI::COMM_WORLD.Send(&clusterWeight, nodes, MPI::DOUBLE, MPI_MASTER, 0);
}

void WorkerNode::computeCurrentMass()
{
  double currentMassDiagonal = 0.0;
  double currentMass = 0.0;

  for (const GraphHypothesisCluster& cluster : m_clusters) {
    double crtWeight = cluster.getWeight();
    currentMass += crtWeight;
    currentMassDiagonal += crtWeight * crtWeight;
  }

  MPI::COMM_WORLD.Reduce(&currentMass, NULL, 1, MPI::DOUBLE,
      MPI::SUM, MPI_MASTER);
  MPI::COMM_WORLD.Reduce(&currentMassDiagonal, NULL, 1, MPI::DOUBLE,
      MPI::SUM, MPI_MASTER);
}

void WorkerNode::recomputePartialTestScores()
{
  // No requests will come for random policy.
  if (m_config.objType == RANDOM)
    return;

  // Distribute partial sum for the mass.
  computeCurrentMass();

  // Recompute requests coming from master.
  int currentTestNode = 0;
  MPI::COMM_WORLD.Bcast(&currentTestNode, 1, MPI::INT, MPI_MASTER);

  double massBuffer[m_config.objSums];
  while (currentTestNode != -1) {
    memset(massBuffer, 0, m_config.objSums * sizeof(*massBuffer));

    // Compute partial mass scores.
    GraphTest test(currentTestNode);
    for (const GraphHypothesisCluster& cluster : m_clusters) {
      pair<double, double> mass = cluster.computeMassWithTest(test);
      if (m_config.objType == VOI) {
        double expectedMass = m_testsPrior[currentTestNode] * mass.first +
            (1 - m_testsPrior[currentTestNode]) * mass.second;
        massBuffer[0] = max(massBuffer[0], expectedMass);
      } else {
        massBuffer[POSITIVE_SUM] += mass.first;
        massBuffer[NEGATIVE_SUM] += mass.second;
        if (m_config.objType == EC2) {
          massBuffer[POSITIVE_DIAG_SUM] += mass.first * mass.first;
          massBuffer[NEGATIVE_DIAG_SUM] += mass.second * mass.second;
        }
      }
    }

    // Return the information to the master node.
    MPI::Op op = MPI::SUM;
    if (m_config.objType == VOI)
      op = MPI::MAX;
    MPI::COMM_WORLD.Reduce(massBuffer, NULL, m_config.objSums,
        MPI::DOUBLE, op, MPI_MASTER);

    // Get the next node to be evaluated.
    MPI::COMM_WORLD.Bcast(&currentTestNode, 1, MPI::INT, MPI_MASTER);
  }
}
