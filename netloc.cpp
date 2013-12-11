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

// Needs to be included after Snap..
#include <future>

#define TRIALS 5
#define THREADS 8

// Snap defines its own macros of max(), min() and this doesn't allow the
// proper use of numeric_limits<int>::min()/max(), therefore undefine them.
#undef max
#undef min

typedef unique_ptr<GraphHypothesisCluster> PGraphHypothesisCluster;

using namespace std;

inline void generateNetwork(PUNGraph *network, const SimConfig& config) {
  // To easily swap the generation model later.
  *network = TSnap::GenRndGnm<PUNGraph>(config.nodes, config.edges);
}

inline void generateClusters(vector<GraphHypothesisCluster> *clusters,
                             const PUNGraph& network,
                             const SimConfig& config)
{
  clusters->clear();

  // TODO(vcarbune): figure out why the code below breaks..
  /*
  vector<future<GraphHypothesisCluster>> futures;
  for (int source = 0; source < network->GetNodes(); ++source) {
    futures.push_back(async(launch::async,
          [](const PUNGraph& net, const SimConfig& cfg, int src) {
              return GraphHypothesisCluster(net, src, cfg.clusterSize,
                cfg.beta, cfg.cascadeBound);
          }, ref(network), ref(config), source));
  }

  for (future<GraphHypothesisCluster>& ftr : futures) {
    ftr.wait();
    clusters->push_back(move(ftr.get()));
  }
  */

  // Generate all possible hypothesis clusters that we want to search through.
  for (int source = 0; source < network->GetNodes(); source++)
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
  time_t start = time(NULL);

  double crtScore;
  vector<double> scoreList;

  PUNGraph network;
  vector<GraphHypothesisCluster> clusters;
  vector<GraphTest> tests;

  generateNetwork(&network, config);
  generateClusters(&clusters, network, config);
  generateTests(&tests, network);

  cout << "Initialization: " << difftime(time(NULL), start) <<
    " seconds " << endl;

  int fails = 0;
  scoreList.push_back(config.getSimParamValue());

  for (int trials = 0; trials < TRIALS; ++trials) {
    // Initialize temporary variables.
    vector<GraphHypothesisCluster> remainingClusters;
    vector<GraphHypothesisCluster> tempClusters(clusters);
    vector<GraphTest> tempTests(tests);

    // Select a different realization at each run.
    int index = rand() % tempClusters.size();
    GraphHypothesis realization = tempClusters[index].generateHypothesis(true);

    // Run the simulation with the current configuration.
    crtScore = runEC2<GraphTest, GraphHypothesisCluster, GraphHypothesis>(
        tempTests, tempClusters, realization, remainingClusters);

    if (!remainingClusters.size()) {
      scoreList.push_back(crtScore);
      fails++;
      continue;
    }

    bool found = false;
    cout << "Correct: " << tempClusters[index].getSource() << "\t";
    cout << "Candidates: ";
    for (const GraphHypothesisCluster& cluster : remainingClusters) {
      if (cluster.getSource() == tempClusters[index].getSource())
        found = true;

      cout << cluster.getSource() << "(" <<
        TSnap::GetShortPath(network, tempClusters[index].getSource(), cluster.getSource()) <<
        " hops)\t";
    }
    cout << "Tests: " << crtScore << endl;

    if (!found)
      fails++;

    scoreList.push_back(crtScore);
  }

  scoreList.push_back((double) fails / TRIALS);
  runStats->push_back(scoreList);

  cout << "Runtime: " << difftime(time(NULL), start) << " seconds " << endl;

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
