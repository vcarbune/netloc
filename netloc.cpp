/**
 * Network and Cascades Generator.
 */

#include <iostream>

#include "hypothesis.h"
#include "test.h"
#include "snap/snap-core/Snap.h"

using namespace std;

int main(int argc, char *argv[])
{
  // TODO(vcarbune): Set running envinronment variables.

  // Get underlying network.
  PUNGraph network = TSnap::GenRndGnm<PUNGraph>(100, 100);

  // Get all possible hypothesis clusters that we want to search through.
  vector<HypothesisCluster> clusters;
  for (int source = 0; source < network->GetNodes(); ++source)
    clusters.push_back(HypothesisCluster(network, source, 10));

  // Get true hypothesis.
  Hypothesis realization =
    clusters[rand() % clusters.size()].getRandomHypothesis();

  cout << "Graph " << network->GetNodes() << " " << network->GetEdges() << endl;
  cout << "Hypothesis: " << clusters.size() << " clusters";

  for (size_t i = 0; i < clusters.size(); ++i)
    clusters[i].printState();

  // Rule out hypothesis incompatible with t.
  Test t(0, 1);
  for (size_t i = 0; i < clusters.size(); ++i)
    clusters[i].ruleOutHypothesis(t);

  return 0;
}
