/**
 * Main implementation of semi-supervised learning for source localization.
 * Author: Victor Carbune (vcarbune@ethz.ch)
 */

#include <algorithm>
#include <limits>
#include <iostream>
#include <queue>

#include <cassert>
#include <ctime>
#include <cmath>

#include "ec2.h"
#include "hypothesis.h"
#include "test.h"
#include "utils.h"

#include "snap/snap-core/Snap.h"

#define TRIALS 5

// Snap defines it's own macros of max(), min() and this doesn't allow the
// proper use of numeric_limits<int>::min()/max(), therefore undefine them.
#undef max
#undef min

using namespace std;

void runSimulation(const PUNGraph& network,
                   const vector<GraphHypothesisCluster>& clusters,
                   const vector<GraphTest>& tests,
                   vector<vector<double>> *runStats)
{

  double crtScore;
  vector<double> scoreList;

  /*
  scoreList[0] = 0;
  scoreList[1] = numeric_limits<double>::max();
  scoreList[2] = numeric_limits<double>::min();
  */

  for (int trials = 0; trials < TRIALS; ++trials) {
    // Initialize temporary variables.
    vector<GraphHypothesisCluster> tempClusters = clusters;
    vector<GraphTest> tempTests = tests;

    // Select a different realization at each run.
    int index = rand() % tempClusters.size();
    GraphHypothesis realization = tempClusters[index].getRandomHypothesis();

    // Run the simulation with the current configuration.
    crtScore = runEC2<GraphTest, GraphHypothesisCluster, GraphHypothesis>(
        tempTests, tempClusters, realization);

    scoreList.push_back(crtScore);
    /*
    scoreList[0] += crtScore;
    scoreList[1] = min(scoreList[1], crtScore);
    scoreList[2] = max(scoreList[2], crtScore);
    */
  }
  /*
  scoreList[0] /= TRIALS;
  */
  runStats->push_back(scoreList);
}

inline PUNGraph generateNetwork(int nodes, int edges) {
  // To easily swap the generation model later.
  return TSnap::GenRndGnm<PUNGraph>(nodes, edges);
}

inline void generateClusters(vector<GraphHypothesisCluster> *clusters,
                             const PUNGraph& network,
                             int clusterSize = 10,
                             double beta = 0.1,
                             double size = 0.4)
{
  // Generate all possible hypothesis clusters that we want to search through.
  clusters->clear();
  for (int source = 0; source < network->GetNodes(); ++source)
    clusters->push_back(
        GraphHypothesisCluster(network, source, clusterSize, beta, size));
}

inline void generateTests(vector<GraphTest> *tests,
                          const PUNGraph& network)
{
  // Generate all tests that are enough to differentiate between hypothesis.
  tests->clear();
  for (int node = 0; node < network->GetNodes(); ++node) {
    assert(network->IsNode(node));
    tests->push_back(GraphTest(node));
  }
}

void generateSimulationStats(ParameterVariationType simulationParameter,
                             vector<vector<double>> *runStats,
                             vector<double> *runParams)
{
  PUNGraph network;
  vector<GraphHypothesisCluster> clusters;
  vector<GraphTest> tests;

  runStats->clear();
  runParams->clear();

  switch (simulationParameter) {
    case NodeVar:
      for (int nodes = 100; nodes < 1001; nodes += 50) {
        network = generateNetwork(nodes, nodes * log(nodes));
        generateClusters(&clusters, network);
        generateTests(&tests, network);

        runSimulation(network, clusters, tests, runStats);
        runParams->push_back(static_cast<double>(nodes));
        cout << "Running with " << nodes << " nodes... " << endl;
      }
      break;
    case EdgeVar:
      for (int edges = 100; edges < 1001; edges += 50) {
        network = generateNetwork((int) pow(edges, 1.1), edges);
        generateClusters(&clusters, network);
        generateTests(&tests, network);

        runSimulation(network, clusters, tests, runStats);
        runParams->push_back(static_cast<double>(edges));
        cout << "Running with " << edges << " edges... " << endl;
      }
      break;
    case BetaVar:
      for (double beta = 0.1; beta <= 0.5; beta += 0.05) {
        network = generateNetwork(500, 1500);
        generateClusters(&clusters, network, 10, beta, 0.4);
        generateTests(&tests, network);

        runSimulation(network, clusters, tests, runStats);
        runParams->push_back(beta);
        cout << "Running with " << beta << " beta... " << endl;
      }
      break;
    case HypothesisVar:
      for (int clusterSize = 5; clusterSize <= 15; ++clusterSize) {
        network = generateNetwork(500, 1500);
        generateClusters(&clusters, network, clusterSize);
        generateTests(&tests, network);

        runSimulation(network, clusters, tests, runStats);
        runParams->push_back(static_cast<double>(clusterSize));
        cout << "Running with " << clusterSize << " cluster size... " << endl;
      }
      break;
    case CascBoundVar:
      for (double bound = 0.4; bound <= 0.6; bound += 0.05) {
        network = generateNetwork(500, 1500);
        generateClusters(&clusters, network, 10, 0.1, bound);
        generateTests(&tests, network);

        runSimulation(network, clusters, tests, runStats);
        runParams->push_back(bound);
        cout << "Running with " << bound << " cascade size bound... " << endl;
      }
      break;
  }

  cout << endl;
}

void dumpSimulationStats(const vector<double>& runParams,
                         const vector<vector<double>>& runStats,
                         ostream& outstream)
{
  assert(runParams.size() == runStats.size());

  for (size_t trial = 0; trial < runStats.size(); ++trial) {
    outstream << runParams[trial] << "\t";
    for (size_t i = 0; i < runStats[trial].size(); ++i)
      outstream << runStats[trial][i] << "\t";
    outstream << endl;
  }
}

int main(int argc, char *argv[])
{
  vector<double> runParams;
  vector<vector<double>> runStats;

  srand(time(NULL));

  generateSimulationStats(HypothesisVar, &runStats, &runParams);
  dumpSimulationStats(runParams, runStats, cout);

  return 0;
}
