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

#define GROUND_TRUTHS 5

// Snap defines its own macros of max(), min() and this doesn't allow the
// proper use of numeric_limits<int>::min()/max(), therefore undefine them.
#undef max
#undef min

using namespace std;

void runSimulation(PUNGraph network,
                   const vector<GraphHypothesis>& realizations,
                   vector<int>& identificationCount,
                   const SimConfig config,
                   ostream& fout)
{

  double crtScore;
  vector<double> scoreList;
  scoreList.push_back(config.clusterSize);

  time_t startTime = time(NULL);
  vector<GraphHypothesisCluster> clusters;
  vector<GraphTest> tests;

  generateClusters(&clusters, network, config);
  generateTests(&tests, network);

  cout << "Initialization took " << difftime(time(NULL), startTime)
    << " seconds " << endl;

  int fails = 0;
  for (int truth = 0; truth < GROUND_TRUTHS; ++truth) {
    // Initialize temporary variables.
    vector<int> topClusters;
    vector<GraphHypothesisCluster> tempClusters(clusters);
    vector<GraphTest> tempTests(tests);

    // Select a different realization at each run.
    GraphHypothesis realization = realizations[truth];

    // Run the simulation with the current configuration.
    time_t startTime = time(NULL);
    crtScore = runEC2<GraphTest, GraphHypothesisCluster, GraphHypothesis>(
        tempTests, tempClusters, realization, config.lazy,
        config.massThreshold, config.testThreshold);

    bool found = false;
    cout << "Correct: " << realization.getSource() << "\t";
    cout << "Tests: " << crtScore << "%\t";
    cout << "Time: " << difftime(time(NULL), startTime) << "s\t";
    cout << "Candidates: ";

    double mass = 0.0;
    for (const GraphHypothesisCluster& cluster : tempClusters)
      mass += cluster.getWeight();
    for (GraphHypothesisCluster& cluster : tempClusters)
      cluster.normalizeWeight(mass / 100);

    sort(tempClusters.begin(), tempClusters.end());
    for (int i = 0; i < config.topN; ++i) {
      int foundSource = tempClusters[i].getSource();
      if (foundSource == realization.getSource()) {
        identificationCount[truth]++;
        found = true;
      }

      cout << foundSource << "(" <<
          TSnap::GetShortPath(network, realization.getSource(), foundSource) << " hops, " <<
          tempClusters[i].getWeight() << ")\t";
    }

    if (!found)
      fails++;

    cout << endl;

    // Output either the tests exhausted, or mass of the solution at the end.
    switch (config.outputType) {
      case 1:
        for (const GraphHypothesisCluster& hc : tempClusters)
          if (hc.getSource() == realization.getSource()) {
            scoreList.push_back(hc.getWeight());
            break;
          }
        break;
      case 2:
        for (const GraphHypothesisCluster& hc : tempClusters)
          if (hc.getSource() == realization.getSource()) {
            scoreList.push_back(tempClusters[0].getWeight() - hc.getWeight());
            break;
          }
        break;
      default:
        scoreList.push_back(crtScore);
    }
  }

  // Log a summary of success rate.
  scoreList.push_back(1 - (double) fails / GROUND_TRUTHS);

  // If provided, log each simulation step into an output stream.
  if (fout) {
    for (size_t i = 0; i < scoreList.size(); ++i)
      fout << scoreList[i] << "\t";
    fout << endl;
  }
}

void startSimulations(SimConfig& config, ostream& fout)
{
  // Keep network constants over multiple simulation runs.
  PUNGraph network;
  generateNetwork(&network, config);

  // Initialize ground truths (hypothesis for which we are searching).
  vector<GraphHypothesis> realizations;
  for (int truth = 0; truth < GROUND_TRUTHS; truth++)
    realizations.push_back(GraphHypothesis::generateHypothesis(network,
          rand() % network->GetNodes(), config.cascadeBound, config.beta));

  // Keep an identification count of each ground truth.
  vector<int> identificationCount(GROUND_TRUTHS);
  for (int step = 0; step < config.steps; step++, ++config) {
    cout << "Current configuration: " << config << endl;
    runSimulation(network, realizations, identificationCount, config, fout);
  }


  for (int truth = 0; truth < GROUND_TRUTHS; ++truth) {
    fout << (double) identificationCount[truth] / config.steps << " ";
  }

  fout << endl;
  cout << endl;
}

int main(int argc, char *argv[])
{
  SimConfig config = SimConfig::getSimConfigFromEnv(argc, argv);

  srand(time(NULL));

  ofstream dumpStream;
  dumpStream.open(config.logfile.CStr());

  cout << fixed << std::setprecision(2);
  dumpStream << fixed << std::setprecision(2);

  startSimulations(config, dumpStream);

  dumpStream.close();
  return 0;
}
