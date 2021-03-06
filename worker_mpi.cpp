#include "worker_mpi.h"

#include <algorithm>
#include <cstring>
#include <iostream>

#undef max
#undef min

using namespace std;

WorkerNode::WorkerNode(SimConfig config)
  : MPINode(config)
{
  // Note: This assumes the network has been read in the superclass constructor.

  // Determine the node clusters for which this node is responsible.
  int clustersPerNode = m_network->GetNodes() / (m_config.mpi.nodes - 1);
  int startNode = (m_config.mpi.rank-1) * clustersPerNode;
  int endNode = m_network->GetNodes();
  if (m_config.mpi.rank < m_config.mpi.nodes - 1)
    endNode = m_config.mpi.rank * clustersPerNode;
  m_nodeRange.first = startNode;
  m_nodeRange.second = endNode;

  // Initialize clusters!
  initializeClusters();
  initializeNodeInfectionTimeMap();
}

void WorkerNode::runWithCurrentConfig()
{
  // No reason for a distributed version EPFL max-likelihood (distributed
  // matrix computations are a different beast, so keeping, for now, algorithm
  // implementation on single machine, running on central MPI node only)
  if (m_config.objType == EPFL_ML || m_config.objType == EPFL_EC2) {
    cout << "Note: EPFLSolver runs only on central MPI node" << endl;
    return;
  }

  cout << m_config.mpi.rank << ": " <<
      GraphHypothesis::GetConsumedBytes() << " bytes" << endl;
  cout << m_config.mpi.rank << ": " <<
      GraphHypothesisCluster::GetConsumedBytes() << " bytes" << endl;

  for (int step = 0; step < m_config.steps; ++step, ++m_config) {
    reset();
    for (int truth = 0; truth < m_config.groundTruths; ++truth) {
      simulate();
      for (GraphHypothesisCluster& cluster : m_clusters)
        cluster.resetWeight(1);
    }
  }
}

void WorkerNode::initializeClusters()
{
  m_clusters.clear();

  // Actual nodes in the network.
  for (int src = m_nodeRange.first; src < m_nodeRange.second; src++) {
    m_clusters.push_back(GraphHypothesisCluster::generateHypothesisCluster(
          m_network, m_nid[src], 1, m_config));
  }

  double avgSize = 0.0;
  for (const GraphHypothesisCluster& c : m_clusters)
    avgSize += c.getAverageHypothesisSize();

  cout << "Worker " << m_config.mpi.rank << ": " <<
      avgSize / m_clusters.size() << endl;
}

void WorkerNode::initializeNodeInfectionTimeMap()
{
  int size;
  m_histogramInfo.clear();
  m_nodeInfectionTimes.clear();

  // Collect the infection times of nodes.
  for (int node = 0; node < m_network->GetNodes(); ++node) {
    vector<double> collectedInfectionTimes;
    if (!m_config.ignoreTime) {
      for (const GraphHypothesisCluster& cluster : m_clusters)
        cluster.collectInfectionTimes(m_nid[node], collectedInfectionTimes);
      sort(collectedInfectionTimes.begin(), collectedInfectionTimes.end());

      // Send to the central node.
      size = collectedInfectionTimes.size();
      MPI::COMM_WORLD.Gather(&size, 1, MPI::INT, NULL, 0, MPI::INT, MPI_MASTER);

      // Send the actual infection times.
      MPI::COMM_WORLD.Send(&collectedInfectionTimes.front(),
          collectedInfectionTimes.size(), MPI::DOUBLE, MPI_MASTER, 0);
    }

    // Gather the centralized vector.
    MPI::COMM_WORLD.Bcast(&size, 1, MPI::INT, MPI_MASTER);

    double infectionTimes[size];
    MPI::COMM_WORLD.Bcast(infectionTimes, size, MPI::DOUBLE, MPI_MASTER);

    m_histogramInfo.push_back(
        vector<double>(infectionTimes, infectionTimes + size));
    m_nodeInfectionTimes.push_back(
        createDiscreteTimeValuesFromHistogramBounds(m_histogramInfo[node]));
  }
}

void WorkerNode::initializeTestHeap()
{
  if (m_config.objType == RANDOM || m_config.objType == VOI ||
      m_config.objType == EC2_HIGH)
    return;

  recomputePartialTestScores();
}

void WorkerNode::reset()
{
  // cout << "WorkerNode::reset()" << endl;
  m_previousTests.clear();
  if (!m_config.cluster.keep) {
    initializeClusters();
    initializeNodeInfectionTimeMap();
  }
  initializeTestHeap();
}

void WorkerNode::simulate()
{
  m_previousTests.clear();

  double nextPcnt = m_config.sampling;
  for (int count = 0; count < m_config.maxObservers; ++count) {
    recomputePartialTestScores();

    // Receive from the master node the testNode that was selected to run.
    int selectedNode;
    double score;
    double infectionTime;

    MPI::COMM_WORLD.Bcast(&selectedNode, 1, MPI::INT, MPI_MASTER);
    MPI::COMM_WORLD.Bcast(&infectionTime, 1, MPI::DOUBLE, MPI_MASTER);
    MPI::COMM_WORLD.Bcast(&score, 1, MPI::DOUBLE, MPI_MASTER);

    // Update the inner hypothesis weights, by running the testNode.
    GraphTest test(m_nid[selectedNode]);
    test.setInfectionTime(infectionTime);

    for (GraphHypothesisCluster& cluster : m_clusters) {
      cluster.updateMassWithTest(m_config.eps, test, m_previousTests);
    }

    m_previousTests.push_back(make_pair(test.getInfectionTime(), test.getNodeId()));

    if (fabs(count - nextPcnt * m_config.nodes) < 0.5) {
      sendClusterData();
      nextPcnt += m_config.sampling;
    }
  }
  sendClusterData();
}

void WorkerNode::sendClusterData()
{
  // Send the cluster masses to the central node.
  int clusterNodes[m_clusters.size()];
  double clusterWeight[m_clusters.size()];
  for (size_t idx = 0; idx < m_clusters.size(); ++idx) {
    clusterNodes[idx] = m_clusters[idx].getSource();
    clusterWeight[idx] = m_clusters[idx].getWeight();
  }

  int nodes = m_clusters.size();
  MPI::COMM_WORLD.Send(&nodes, 1, MPI::INT, MPI_MASTER, 0);
  MPI::COMM_WORLD.Send(&clusterNodes, nodes, MPI::INT, MPI_MASTER, 0);
  MPI::COMM_WORLD.Send(&clusterWeight, nodes, MPI::DOUBLE, MPI_MASTER, 0);
}

void WorkerNode::computeCurrentWeight()
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

void WorkerNode::recomputeVoIScores()
{
  computeCurrentWeight();

  int currentTestNode = 0;
  MPI::COMM_WORLD.Bcast(&currentTestNode, 1, MPI::INT, MPI_MASTER);
  while (currentTestNode != -1) {
    int totalInfectionTimes = m_nodeInfectionTimes[currentTestNode].size();
    double priors[totalInfectionTimes];
    memset(priors, 0, totalInfectionTimes * sizeof(*priors));

    // Cache cluster masses with infection time.
    vector<vector<double>> clusterMassTestInfection;

    // For each infection time, compute mass and prior.
    GraphTest test(m_nid[currentTestNode]);
    for (const GraphHypothesisCluster& cluster : m_clusters) {
      vector<double> clusterMasses;
      for (size_t i = 0; i < m_nodeInfectionTimes[currentTestNode].size(); ++i) {
        test.setInfectionTime(m_nodeInfectionTimes[currentTestNode][i]);
        pair<double, double> vals = cluster.computeMassWithTest(
            m_config.eps, test, m_config.ignoreTime, m_previousTests,
            m_histogramInfo[currentTestNode]);
        priors[i] += vals.first;
        clusterMasses.push_back(vals.second);
      }
      clusterMassTestInfection.push_back(clusterMasses);
    }

    // Centralize priors to master node.
    MPI::COMM_WORLD.Reduce(
      priors, NULL, totalInfectionTimes, MPI::DOUBLE, MPI::SUM, MPI_MASTER);
    MPI::COMM_WORLD.Bcast(priors, totalInfectionTimes, MPI::DOUBLE, MPI_MASTER);

    // Compute maximum score for this test.
    double maxScore = -numeric_limits<double>::max();
    for (size_t c = 0; c < m_clusters.size(); ++c) {
      double expectedWeight = 0;
      for (size_t i = 0; i < m_nodeInfectionTimes[currentTestNode].size(); ++i)
        expectedWeight += priors[i] * clusterMassTestInfection[c][i];
      maxScore = max(maxScore, m_clusters[c].getWeight() - expectedWeight);
    }

    MPI::COMM_WORLD.Reduce(&maxScore, NULL, 1, MPI::DOUBLE, MPI::MAX, MPI_MASTER);
    MPI::COMM_WORLD.Bcast(&currentTestNode, 1, MPI::INT, MPI_MASTER);
  }
}

void WorkerNode::recomputePartialTestScores()
{
  // No requests will come for random policy.
  if (m_config.objType == RANDOM || m_config.objType == EC2_HIGH)
    return;

  if (m_config.objType == VOI) {
    recomputeVoIScores();
    return;
  }

  computeCurrentWeight();

  // Recompute requests coming from master.
  int currentTestNode = 0;
  MPI::COMM_WORLD.Bcast(&currentTestNode, 1, MPI::INT, MPI_MASTER);
  while (currentTestNode != -1) {
    int totalInfectionTimes = m_nodeInfectionTimes[currentTestNode].size();

    double mass[totalInfectionTimes];
    memset(mass, 0, totalInfectionTimes * sizeof(*mass));

    double sqMass[totalInfectionTimes];
    memset(sqMass, 0, totalInfectionTimes * sizeof(*sqMass));

    double priors[totalInfectionTimes];
    memset(priors, 0, totalInfectionTimes * sizeof(*priors));

    // For each infection time, compute mass and prior.
    GraphTest test(m_nid[currentTestNode]);
    for (const GraphHypothesisCluster& cluster : m_clusters) {
      for (size_t i = 0; i < m_nodeInfectionTimes[currentTestNode].size(); ++i) {
        test.setInfectionTime(m_nodeInfectionTimes[currentTestNode][i]);
        pair<double, double> vals = cluster.computeMassWithTest(
            m_config.eps, test, m_config.ignoreTime, m_previousTests,
            m_histogramInfo[currentTestNode]);
        priors[i] += vals.first;
        mass[i] += vals.second;
        sqMass[i] += vals.second * vals.second;
      }
    }

    // Priors, mass, squared mass..
    MPI::COMM_WORLD.Reduce(
        priors, NULL, totalInfectionTimes, MPI::DOUBLE, MPI::SUM, MPI_MASTER);
    MPI::COMM_WORLD.Reduce(
        mass, NULL, totalInfectionTimes, MPI::DOUBLE, MPI::SUM, MPI_MASTER);
    MPI::COMM_WORLD.Reduce(
        sqMass, NULL, totalInfectionTimes, MPI::DOUBLE, MPI::SUM, MPI_MASTER);

    // Get the next node to be evaluated.
    MPI::COMM_WORLD.Bcast(&currentTestNode, 1, MPI::INT, MPI_MASTER);
  }
}
