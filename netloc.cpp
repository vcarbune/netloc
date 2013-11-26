/**
 * Main implementation of semi-supervised learning for source localization.
 * Author: Victor Carbune (vcarbune@ethz.ch)
 */

#include <algorithm>
#include <limits>
#include <iostream>
#include <fstream>
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
  }

  runStats->push_back(scoreList);
}

inline void generateNetwork(PUNGraph *network, const SimConfig& config) {
  // To easily swap the generation model later.
  *network = TSnap::GenRndGnm<PUNGraph>(config.nodes, config.edges);
}

inline void generateClusters(vector<GraphHypothesisCluster> *clusters,
                             const PUNGraph& network,
                             const SimConfig& config)
{
  // Generate all possible hypothesis clusters that we want to search through.
  clusters->clear();
  for (int source = 0; source < network->GetNodes(); ++source)
    clusters->push_back(
        GraphHypothesisCluster(network, source,
              config.clusterSize, config.beta, config.cascadeBound));
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

void generateSimulationStats(SimulationType simulationParameter,
                             vector<vector<double>> *runStats,
                             vector<double> *runParams,
                             ostream& fout)
{
  PUNGraph network;
  vector<GraphHypothesisCluster> clusters;
  vector<GraphTest> tests;

  runStats->clear();
  runParams->clear();

  SimConfig config(simulationParameter);
  for (int step = 0; step < 5; step++, ++config) {
    cout << "Running one simulation step ... " << config.getSimParamValue();

    time_t start = time(NULL);

    generateNetwork(&network, config);
    generateClusters(&clusters, network, config);
    generateTests(&tests, network);

    runSimulation(network, clusters, tests, runStats);
    runParams->push_back(static_cast<double>(config.getSimParamValue()));

    // If provided, log each simulation step into an output stream.
    if (fout) {
      fout << config.getSimParamValue() << "\t";
      for (size_t i = 0; i < (*runStats)[step].size(); ++i)
        fout << (*runStats)[step][i] << "\t";
      fout << endl;
    }

    cout << " (" << difftime(time(NULL), start) << " sec)" << endl;
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

  ofstream dumpStream;
  dumpStream.open("dump.log");

  generateSimulationStats(NodeVar, &runStats, &runParams, dumpStream);
  dumpSimulationStats(runParams, runStats, cout);

  dumpStream.close();
  return 0;
}
