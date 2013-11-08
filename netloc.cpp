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

int main(int argc, char *argv[])
{
  // Get underlying network.
  PUNGraph network = TSnap::GenRndGnm<PUNGraph>(100, 4800);

  // Get all possible hypothesis clusters that we want to search through.
  vector<GraphHypothesisCluster> clusters;
  for (int source = 0; source < network->GetNodes(); ++source)
    clusters.push_back(GraphHypothesisCluster(network, source, 10));

  // Get true hypothesis.
  srand(time(NULL));
  int index = rand() % clusters.size();
  GraphHypothesis realization = clusters[index].getRandomHypothesis();

  // Create tests.
  vector<GraphTest> tests;
  for (int node = 0; node < network->GetNodes(); ++node) {
    assert(network->IsNode(node));
    tests.push_back(GraphTest(node));
  }

  runEC2<GraphTest, GraphHypothesisCluster, GraphHypothesis>(
      tests, clusters, realization);

  cout << "True hypothesis cluster: " << index << endl;
  cout << "Tests remaining: " << tests.size() << endl;

  return 0;
}
