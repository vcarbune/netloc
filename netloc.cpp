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

// Snap defines its own macros of max(), min() and this doesn't allow the
// proper use of numeric_limits<int>::min()/max(), therefore undefine them.
#undef max
#undef min

#include "Eigen/Dense"

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
  for (int truth = 0; truth < config.groundTruths; ++truth) {
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
        config.massThreshold, config.testThreshold, config.objType);

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
  scoreList.push_back(1 - (double) fails / config.groundTruths);

  // If provided, log each simulation step into an output stream.
  if (fout) {
    for (size_t i = 0; i < scoreList.size(); ++i)
      fout << scoreList[i] << "\t";
    fout << endl;
  }
}

// EPFL Solution
void solveUsingML(PUNGraph network,
    const vector<GraphHypothesis>& realizations,
    vector<int>& identificationCount,
    const SimConfig config,
    ostream& fout)
{
  // Select 34% of observers using high-degree nodes.
  int observers = config.testThreshold * network->GetNodes();
  vector<pair<int, int>> nodeDegrees;
  for (int node = 0; node < network->GetNodes(); ++node)
    nodeDegrees.push_back(
        pair<int, int>(network->GetNI(node).GetOutDeg(), node));
  sort(nodeDegrees.begin(), nodeDegrees.end(), std::greater<pair<int, int>>());

  for (size_t truth = 0; truth < realizations.size(); ++truth) {
    GraphHypothesis realization = realizations[truth];
    // Get infection times for all the observers and keep ascending.
    vector<pair<double, int>> observerNodes;
    for (int o = 0; o < observers; ++o) {
      observerNodes.push_back(pair<double, int>(
            realization.getInfectionTime(nodeDegrees[o].second),
            nodeDegrees[o].second));
    }

    // Everything is related to the reference observer.
    sort(observerNodes.begin(), observerNodes.end());
    int referenceObserver = observerNodes[0].second;
    TIntH idToShortestPathsFromReference;
    TSnap::GetShortPath(network, referenceObserver, idToShortestPathsFromReference);

    // BFS tree from reference observer. Compute LCA w.r.t. to ref. observer.
    PNGraph bfsTree = TSnap::GetBfsTree(network, referenceObserver, true, true);
    int lca[observers][observers];
    for (int k = 1; k < observers; ++k) {
      for (int i = k + 1; i < observers; ++i) {
        int Ni = observerNodes[i].second;
        int Nk = observerNodes[k].second;
        int heightNi = idToShortestPathsFromReference.GetDat(Ni);
        int heightNk = idToShortestPathsFromReference.GetDat(Nk);
        while (heightNi > heightNk && Ni != referenceObserver) {
          // Note: fingers crossed that the same values are in GetInNId array.
          Ni = bfsTree->GetNI(Ni).GetInNId(0);
          heightNi = idToShortestPathsFromReference.GetDat(Ni);
        }
        while (heightNk > heightNi) {
          // Note: fingers crossed that the same values are in GetInNId array.
          Nk = bfsTree->GetNI(Nk).GetInNId(0);
          heightNk = idToShortestPathsFromReference.GetDat(Nk);
        }
        while (Ni != Nk) {
          Nk = bfsTree->GetNI(Nk).GetInNId(0);
          Ni = bfsTree->GetNI(Ni).GetInNId(0);
        }
        if (heightNk != heightNi || Ni != Nk)
          cout << "Something awfully wrong here" << endl;
        lca[k][i] = lca[i][k] = Ni;
      }
    }

    // Delay Covariance
    Eigen::MatrixXf lambda(observers - 1, observers - 1);
    for (int k = 0; k < observers - 1; ++k) {
      for (int i = 0; i < observers - 1; ++i) {
        float value = config.epflSigma * config.epflSigma;
        if (k == i)
          value *= idToShortestPathsFromReference.GetDat(observerNodes[k+1].second);
        else
          value *= idToShortestPathsFromReference.GetDat(lca[k+1][i+1]);
        lambda(k, i) = value;
      }
    }
    Eigen::MatrixXf invLambda(observers - 1, observers - 1);
    invLambda = lambda.inverse();

    // Observed Delay
    Eigen::VectorXf d(observers - 1);
    for (int o = 0; o < observers - 1; ++o)
      d[o] = observerNodes[o+1].first - observerNodes[o].first;

    vector<pair<double, int>> scores;
    for (int s = 0; s < network->GetNodes(); ++s) {
      TIntH idToShortestPathsFromSource;
      TSnap::GetShortPath(network, s, idToShortestPathsFromSource);

      // Deterministic Delay
      Eigen::VectorXf miu_s(observers - 1);
      for (int o = 0; o < observers - 1; ++o)
        miu_s[o] = config.epflMiu *
          (idToShortestPathsFromSource.GetDat(observerNodes[o+1].second) -
           idToShortestPathsFromSource.GetDat(referenceObserver));

      // Estimator value.
      double estimator = miu_s.transpose() * invLambda * (0.5 * miu_s - d);
      scores.push_back(pair<double, int>(estimator, s));
    }
    sort(scores.begin(), scores.end(), std::greater<pair<double, int>>());
    cout << "Expected: " << realization.getSource() <<
      " Found: " << scores[0].second << "(" << scores[0].first << ")" << endl;
    if (realization.getSource() == scores[0].second)
      identificationCount[truth]++;
  }
}

void startSimulations(SimConfig& config, ostream& fout)
{
  // Keep network constants over multiple simulation runs.
  PUNGraph network;
  generateNetwork(&network, config);

  // Initialize ground truths (hypothesis for which we are searching).
  vector<GraphHypothesis> realizations;
  for (int truth = 0; truth < config.groundTruths; truth++) {
    if (config.epflSolver){
      realizations.push_back(
          GraphHypothesis::generateHypothesisUsingGaussianModel(
            network, rand() % network->GetNodes(), config.cascadeBound,
            config.epflMiu, config.epflSigma));
    } else {
      realizations.push_back(GraphHypothesis::generateHypothesis(network,
            rand() % network->GetNodes(), config.cascadeBound, config.beta));
    }
  }

  // Keep an identification count of each ground truth.
  vector<int> identificationCount(config.groundTruths);
  for (int step = 0; step < config.steps; step++, ++config) {
    if (config.epflSolver) {
      solveUsingML(network, realizations, identificationCount, config, fout);
    } else {
      cout << "Current configuration: " << config << endl;
      runSimulation(network, realizations, identificationCount, config, fout);
    }
  }

  double totalCount = 0.0;
  for (int truth = 0; truth < config.groundTruths; ++truth) {
    fout << (double) identificationCount[truth] / config.steps << " ";
    totalCount += identificationCount[truth] / config.steps;
  }
  cout << totalCount / config.groundTruths << endl;

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
