/**
 * Main implementation of semi-supervised learning for source localization.
 * Author: Victor Carbune (vcarbune@ethz.ch)
 */

#include <algorithm>
#include <iostream>
#include <queue>

#include <cassert>
#include <ctime>

#include "ec2.h"
#include "hypothesis.h"
#include "test.h"

#include "snap/snap-core/Snap.h"

using namespace std;

void plotGraphs(const PUNGraph& network, const vector<GraphTest>& tests,
                const GraphHypothesis& realization)
{
  TIntStrH nodeTestLabels;
  TIntStrH nodeInfectionLabels;

  for (size_t i = 0; i < tests.size(); ++i) {
    int id = tests[i].getNodeId();
    nodeTestLabels.AddDat(id, TStr::Fmt("%d", i));
  }

  for (int i = 0; i < network->GetNodes(); ++i) {
    int infTime = realization.getInfectionTime(i);
    if (infTime != -1)
      nodeInfectionLabels.AddDat(i, TStr::Fmt("%d", infTime));
    else
      nodeInfectionLabels.AddDat(i, TStr::Fmt(""));

    if (!nodeTestLabels.IsKey(i))
      nodeTestLabels.AddDat(i, TStr::Fmt(""));
  }

  // Visualize the order in which the nodes have been tested.
  TStr FNameDemo = TStr::Fmt("graph_with_test_order.%s", "png");
  remove(FNameDemo.CStr());
  TSnap::DrawGViz(network, TGVizLayout(2), FNameDemo, "Tests", nodeTestLabels);
  printf("Drawing graph '%s'\n", FNameDemo.CStr());

  // Visualization true infection order.
  FNameDemo = TStr::Fmt("graph_with_infection_order.%s", "png");
  remove(FNameDemo.CStr());
  TSnap::DrawGViz(network, TGVizLayout(2), FNameDemo, "Infection", nodeInfectionLabels);
  printf("Drawing graph '%s'\n", FNameDemo.CStr());
}

int runSimulation(const PUNGraph& network, int H, bool plot)
{
  // Get all possible hypothesis clusters that we want to search through.
  vector<GraphHypothesisCluster> clusters;
  for (int source = 0; source < network->GetNodes(); ++source)
    clusters.push_back(GraphHypothesisCluster(network, source, H));

  // Get true hypothesis.
  int index = rand() % clusters.size();
  GraphHypothesis realization = clusters[index].getRandomHypothesis();

  // Create tests.
  vector<GraphTest> tests;
  for (int node = 0; node < network->GetNodes(); ++node) {
    assert(network->IsNode(node));
    tests.push_back(GraphTest(node));
  }

  vector<GraphTest> testsInOrder =
    runEC2<GraphTest, GraphHypothesisCluster, GraphHypothesis>(
        tests, clusters, realization);

  // Visual plot of how the algorithm works.
  if (plot)
    plotGraphs(network, testsInOrder, realization);

  /*
  int nonEmptyClusters = 0;
  for (GraphHypothesisCluster cluster : clusters)
    if (cluster.countHypothesisAvailable())
      nonEmptyClusters++;

  cout << "Clusters left: " << nonEmptyClusters << endl;
  cout << "True hypothesis cluster: " << index << endl;
  cout << "Tests remaining: " << tests.size() << endl;
  */

  return testsInOrder.size();
}

int main(int argc, char *argv[])
{
  // Get underlying network.
  PUNGraph network = TSnap::GenRndGnm<PUNGraph>(300, 1000);
  double score;
  int runs = 10;

  srand(time(NULL));
  for (int hSize = 1; hSize < 15; ++hSize) {
    score = 0;
    for (int run = 0; run < runs; ++run) {
      score += runSimulation(network, 1, false);
    }
    cout << hSize << "\t" << score / runs << endl;
  }

  return 0;
}
