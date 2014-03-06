#include "worker_mpi.h"

#include <cstring>
#include <iostream>

#undef max
#undef min

void reducePartialTestScores(
    bool *testWasUsed,
    const vector<GraphHypothesisCluster>& clusters,
    const SimConfig &config)
{
  double massBuffer[config.objSums][config.nodes];
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
      if (config.objType == 0) { // Only for EC2
        massBuffer[POSITIVE_DIAG_SUM][testNode] += mass.first * mass.first;
        massBuffer[NEGATIVE_DIAG_SUM][testNode] += mass.second * mass.second;
      }
      massBuffer[CONS_HYPO_SUM][testNode] += cluster.getNodeCount(testNode);
    }
  }

  double currentMassDiagonal = 0.0;
  double currentMass = 0.0;

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
}

GraphTest selectBestLocalTest(vector<GraphTest> &tests,
    const vector<GraphHypothesisCluster>& clusters,
    const SimConfig& config)
{
  // Local LazyGreedy approach.
  GraphTest bestTest = GraphTest(-1);

  if (tests.size())
    bestTest = lazyRescoreTests(tests, clusters, config.objType);

  // Inform the master about the score this test has (how much much it reduces)
  int testNodeId = bestTest.getNodeId();

  // Normalize scores relative to the weight of the local graph.
  double currentWeight = computeGraphWeight(clusters, config.objType);
  double score = 100 * bestTest.getScore() / currentWeight;

  MPI::COMM_WORLD.Gather(&testNodeId, 1, MPI::INT, NULL, 1, MPI::INT, MPI_MASTER);
  MPI::COMM_WORLD.Gather(&score, 1, MPI::DOUBLE, NULL, 1, MPI::DOUBLE, MPI_MASTER);

  return bestTest;
}

void simulate(
    vector<GraphTest> &tests,
    vector<GraphHypothesisCluster>& clusters,
    const SimConfig& config)
{
  GraphTest bestTest(-1); // only lazy greedy

  bool testWasUsed[config.nodes];
  for (int i = 0; i < config.nodes; ++i)
    testWasUsed[i] = false;

  for (GraphTest& test : tests)
    rescoreTest(test, clusters, config.objType);

  make_heap<vector<GraphTest>::iterator, TestCompareFunction>(
      tests.begin(), tests.end(), TestCompareFunction());

  int totalTests = config.testThreshold * config.nodes;
  for (int count = 0; count < totalTests; ++count) {
    if (config.objType == 3) {
      // NOTHING; test is selected at random!
    } else if (config.lazy) {
      bestTest = selectBestLocalTest(tests, clusters, config);
    } else {
      reducePartialTestScores(testWasUsed, clusters, config);
    }

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

    if (config.lazy && !(test == bestTest)) {
      tests.push_back(bestTest);
      push_heap<vector<GraphTest>::iterator, TestCompareFunction>(
          tests.begin(), tests.end(), TestCompareFunction());
    }

    for (GraphTest& test : tests)
      test.setScore((1-EPS) * (1-EPS) * test.getScore());
  }

  // Identify the cluster with the highest mass and send it back
  // to the master process to identify the highest one.
  double maxMass = -numeric_limits<double>::max();
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

  // If lazy greedy is enabled, send back mass for all clusters directly.
  if (config.lazy) {
    double clusterMass[config.nodes];
    for (int i = 0; i < config.nodes; ++i)
      clusterMass[i] = clusters[i].getMass();

    MPI::COMM_WORLD.Reduce(clusterMass, NULL, config.nodes,
        MPI::DOUBLE, MPI::SUM, MPI_MASTER);
  } else {
    MPI::COMM_WORLD.Gather(&maxMass, 1, MPI::DOUBLE, NULL, 1, MPI::DOUBLE, MPI_MASTER);
    MPI::COMM_WORLD.Gather(&sourceNode, 1, MPI::INT, NULL, 1, MPI::INT, MPI_MASTER);
  }
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

  int mpiTestStartNode = startNode;
  int mpiTestEndNode = endNode;

  if (!config.lazy) {
    // Tests don't matter if there's no lazy greedy approach.
    mpiTestStartNode = mpiTestEndNode = 0;
  } else {
    // If we have lazy greedy enabled, then all the clusters are needed.
    mpiClusterStartNode = 0;
    mpiClusterEndNode = network->GetNodes();
  }

  for (int step = 0; step < config.steps; ++step, ++config) {
    // Generate clusters.
    vector<GraphHypothesisCluster> clusters;
    for (int source = mpiClusterStartNode; source < mpiClusterEndNode; source++)
      clusters.push_back(GraphHypothesisCluster::generateHypothesisCluster(
            network, source, 1, config.beta, config.cascadeBound,
            config.clusterSize));

    for (int truth = 0; truth < config.groundTruths; ++truth) {
      for (GraphHypothesisCluster& cluster : clusters)
        cluster.resetWeight(1);

      // Configure the tests for which this MPI slave node is responsible.
      vector<GraphTest> tests;
      for (int node = mpiTestStartNode; node < mpiTestEndNode; ++node)
        tests.push_back(GraphTest(node));

      simulate(tests, clusters, config);

      if (config.lazy)
        continue;

      int realSource;
      MPI::COMM_WORLD.Bcast(&realSource, 1, MPI::INT, MPI_MASTER);
      if (realSource >= mpiClusterStartNode && realSource < mpiClusterEndNode) {
        int realSourceIdx = realSource - mpiClusterStartNode;
        double realSourceMass = clusters[realSourceIdx].getWeight();
        if (clusters[realSourceIdx].getSource() != realSource)
          cout << "!!ERROR!!" << endl;
        MPI::COMM_WORLD.Send(&realSourceMass, 1, MPI::DOUBLE, MPI_MASTER, 1);
      }
    }
  }
}
