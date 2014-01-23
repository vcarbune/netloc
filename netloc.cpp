/**
 * Main implementation of semi-supervised learning for source localization.
 * Author: Victor Carbune (vcarbune@ethz.ch)
 */

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <memory>
#include <limits>
#include <fstream>

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

#define TRIALS 7

// Snap defines its own macros of max(), min() and this doesn't allow the
// proper use of numeric_limits<int>::min()/max(), therefore undefine them.
#undef max
#undef min

using namespace std;

inline void generateNetwork(PUNGraph *network, const SimConfig& config) {
  // Keep the same underlying network if the
  // simulation parameter is not the number of nodes.
  if (!network->Empty()) {
    cout << "Keeping the previously initialized network." << endl;
    return;
  }

  switch (config.networkType) {
    case 0:
      *network = TSnap::ConvertGraph<PUNGraph, PNGraph>(TSnap::GenForestFire(config.nodes, 0.35, 0.32));
      break;
    case 1:
      *network = TSnap::GenPrefAttach(config.nodes, 3);
      break;
    case 2:
      *network = TSnap::GenRndGnm<PUNGraph>(
          config.nodes, config.nodes * log(config.nodes));
      break;
  }

  cout << "Generated network is connected: " << TSnap::GetMxWccSz(*network) << endl;
}

inline void generateClusters(vector<GraphHypothesisCluster> *clusters,
                             const PUNGraph& network,
                             const SimConfig& config)
{
  // Generate all possible hypothesis clusters that we want to search through.
  clusters->clear();
  /*
  int maxOutDegree = numeric_limits<int>::min();
  for (int source = 0; source < network->GetNodes(); source++)
    maxOutDegree = max(maxOutDegree, network->GetNI(source).GetOutDeg());
  */

  for (int source = 0; source < network->GetNodes(); source++) {
    clusters->push_back(GraphHypothesisCluster(
        network, source,
        config.clusterSize, config.beta, config.cascadeBound, 1));
  }
}

inline void generateTests(vector<GraphTest> *tests,
                          const PUNGraph& network)
{
  // Generate all tests that are enough to differentiate between hypothesis.
  tests->clear();
  for (int node = 0; node < network->GetNodes(); ++node)
    tests->push_back(GraphTest(node));
}

void runSimulation(PUNGraph network,
                   const vector<GraphHypothesis>& realizations,
                   const SimConfig config,
                   ostream& fout,
                   vector<vector<double>> *runStats)
{
  time_t start = time(NULL);

  double crtScore;
  vector<double> scoreList;

  if (realizations.size() != TRIALS) {
    cout << "Different number of trial expected" << endl;
    return;
  }

  vector<GraphHypothesisCluster> clusters;
  vector<GraphTest> tests;

  generateClusters(&clusters, network, config);
  generateTests(&tests, network);

  cout << "Initialization: " << difftime(time(NULL), start) <<
    " seconds " << endl;

  int fails = 0;
  scoreList.push_back(config.getSimParamValue());

  for (int trial = 0; trial < TRIALS; ++trial) {
    // Initialize temporary variables.
    vector<int> topClusters;
    vector<GraphHypothesisCluster> tempClusters(clusters);
    vector<GraphTest> tempTests(tests);

    // Select a different realization at each run.
    GraphHypothesis realization = realizations[trial];

    // Run the simulation with the current configuration.
    time_t startTime = time(NULL);
    crtScore = runEC2<GraphTest, GraphHypothesisCluster, GraphHypothesis>(
        tempTests, tempClusters, realization, config.lazy);

    bool found = false;
    cout << "Correct: " << realization.getSource() << "\t";
    cout << "Tests: " << crtScore << "%\t";
    cout << "Time: " << difftime(time(NULL), startTime) << "s\t";
    cout << "Candidates: ";
    sort(tempClusters.begin(), tempClusters.end());
    for (int i = 0; i < config.topN; ++i) {
      int foundSource = tempClusters[i].getSource();
      if (foundSource == realization.getSource())
        found = true;

      cout << tempClusters[i].getSource() << "(" <<
          TSnap::GetShortPath(network, realization.getSource(), foundSource) << " hops, " <<
          tempClusters[i].getWeight() << ")\t";
    }

    if (!found)
      fails++;

    cout << endl;
    scoreList.push_back(crtScore);
  }

  scoreList.push_back(1 - (double) fails / TRIALS);
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

  PUNGraph network;
  generateNetwork(&network, config);

  vector<GraphHypothesis> realizations;
  for (int trial = 0; trial < TRIALS; trial++)
    realizations.push_back(GraphHypothesis::generateHypothesis(network,
          rand() % network->GetNodes(), config.cascadeBound, config.beta));

  for (int step = 0; step < config.steps; step++, ++config) {
    cout << endl << "Current configuration: " << config.getSimParamValue() << endl;
    runSimulation(network, realizations, config, fout, runStats);
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

  cout << fixed << std::setprecision(2);
  dumpStream << fixed << std::setprecision(2);

  generateSimulationStats(&runStats, config, dumpStream);
  dumpSimulationStats(runStats, cout);

  dumpStream.close();
  return 0;
}
