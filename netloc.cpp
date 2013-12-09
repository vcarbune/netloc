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
#include <thread>

#include "ec2.h"
#include "hypothesis.h"
#include "test.h"
#include "utils.h"

#include "snap/snap-core/Snap.h"

#define TRIALS 20

// Snap defines its own macros of max(), min() and this doesn't allow the
// proper use of numeric_limits<int>::min()/max(), therefore undefine them.
#undef max
#undef min

using namespace std;

inline void generateNetwork(PUNGraph *network, const SimConfig& config) {
  // To easily swap the generation model later.
  *network = TSnap::GenRndGnm<PUNGraph>(config.nodes, config.edges);
}

inline void generateClusters(vector<GraphHypothesisCluster> *clusters,
                             const PUNGraph& network,
                             const SimConfig& config)
{
  // TODO(vcarbune): threadify

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
  // TODO(vcarbune): threadify

  // Generate all tests that are enough to differentiate between hypothesis.
  tests->clear();
  for (int node = 0; node < network->GetNodes(); ++node) {
    assert(network->IsNode(node));
    tests->push_back(GraphTest(node));
  }
}

void runSimulation(const SimConfig config,
                   ostream& fout,
                   vector<vector<double>> *runStats)
{
  double crtScore;
  vector<double> scoreList;

  PUNGraph network;
  vector<GraphHypothesisCluster> clusters;
  vector<GraphTest> tests;

  generateNetwork(&network, config);
  generateClusters(&clusters, network, config);
  generateTests(&tests, network);

  time_t start = time(NULL);
  int fails = 0;
  scoreList.push_back(config.getSimParamValue());

  for (int trials = 0; trials < TRIALS; ++trials) {
    // Initialize temporary variables.
    vector<GraphHypothesisCluster> tempClusters(clusters);
    vector<GraphTest> tempTests(tests);

    // Select a different realization at each run.
    int index = rand() % tempClusters.size();
    GraphHypothesis realization = tempClusters[index].generateHypothesis(true);

    // Run the simulation with the current configuration.
    crtScore = runEC2<GraphTest, GraphHypothesisCluster, GraphHypothesis>(
        tempTests, tempClusters, realization);

    if (crtScore == network->GetNodes())
      fails++;

    scoreList.push_back(crtScore);
  }

  scoreList.push_back((double) fails / TRIALS);
  runStats->push_back(scoreList);

  cout << "Runtime: " << difftime(time(NULL), start) << endl;

  // If provided, log each simulation step into an output stream.
  if (fout) {
    for (size_t i = 0; i < scoreList.size(); ++i)
      fout << scoreList[i] << "\t";
    fout << endl;
  }
}

void generateSimulationStats(vector<vector<double>> *runStats,
                             SimConfig& config,
                             ostream& fout)
{
  runStats->clear();
  vector<thread> threads;

  for (int step = 0; step < config.steps; step++, ++config) {
    cout << endl << "Current configuration: " << config.getSimParamValue() << endl;
    runSimulation(config, fout, runStats);
  }

  for (size_t i = 0; i < threads.size(); ++i)
    threads[i].join();

  cout << endl;
}

void dumpSimulationStats(const vector<vector<double>>& runStats,
                         ostream& outstream)
{
  for (size_t trial = 0; trial < runStats.size(); ++trial) {
    for (size_t i = 0; i < runStats[trial].size(); ++i)
      outstream << runStats[trial][i] << "\t";
    outstream << endl;
  }
}

int main(int argc, char *argv[])
{
  vector<vector<double>> runStats;
  SimConfig config = getSimConfigFromEnv(argc, argv);

  srand(time(NULL));

  ofstream dumpStream;
  dumpStream.open(config.logfile.CStr());

  generateSimulationStats(&runStats, config, dumpStream);
  dumpSimulationStats(runStats, cout);

  dumpStream.close();
  return 0;
}
