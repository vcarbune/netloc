/**
 * Network and Cascades Generator.
 */

#include "hypothesis.h"
#include "test.h"
#include "snap/snap-core/Snap.h"

using namespace std;

int main(int argc, char *argv[])
{
  // Get underlying network.
  PUNGraph network = TSnap::GenRndGnm<PUNGraph>(10, 10);

  // Get all possible hypothesis clusters.
  vector<HypothesisCluster> clusters;
  for (int source = 0; source < network->GetNodes(); ++source)
    clusters.push_back(HypothesisCluster(network, source));

  // Get true hypothesis.
  Hypothesis realization =
    clusters[rand() % clusters.size()].getRandomHypothesis();

  // Rule out hypothesis incompatible with t.
  Test t(0, 1);
  for (size_t i = 0; i < clusters.size(); ++i)
    clusters[i].ruleOutHypothesis(t);

  return 0;
}
