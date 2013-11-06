/**
 * Network and Cascades Generator.
 */

#include <algorithm>
#include <iostream>
#include <queue>

#include <cassert>
#include <ctime>

#include "hypothesis.h"
#include "test.h"
#include "snap/snap-core/Snap.h"

using namespace std;

void rescoreTests(vector<Test>& tests,
                  const vector<HypothesisCluster>& clusters)
{
  for (size_t test = 0; test < tests.size(); ++test) {
    Test& t = tests[test];

    int countConsistentHypothesis = 0;
    int countTotalHypothesis = 0;
    for (size_t i = 0; i < clusters.size(); ++i) {
      countConsistentHypothesis +=
        clusters[i].countHypothesisConsistentWithTest(t);
      countTotalHypothesis += clusters[i].countTotalHypothesis();
    }

    // Score if test outcome is True
    double positiveScore =
      ((double) countConsistentHypothesis / countTotalHypothesis) *
      ((double) 1 / countConsistentHypothesis);
    // Score if test outcome is False
    double negativeScore =
      ((double) (1 - countConsistentHypothesis) / countTotalHypothesis) *
      ((double) 1 / (countTotalHypothesis - countConsistentHypothesis));

    t.first = positiveScore + negativeScore -
      ((double) 1 / countTotalHypothesis);

    t.first *= -1;
  }
}

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
  srand(time(NULL));
  int index = rand() % clusters.size();
  Hypothesis realization = clusters[index].getRandomHypothesis();

  cout << "Graph " << network->GetNodes() << " " << network->GetEdges() << endl;
  cout << "Hypothesis: " << clusters.size() << " clusters" << endl;

  for (size_t i = 0; i < clusters.size(); ++i) {
    if (clusters[i].countTotalHypothesis())
      clusters[i].printState();
  }

  // Rule out hypothesis incompatible with t.
  vector<Test> tests;
  for (int node = 0; node < network->GetNodes(); ++node) {
    assert(network->IsNode(node));
    tests.push_back(Test(0, node));
  }

  int clustersLeft = clusters.size();
  while (!tests.empty() && clustersLeft != 1) {
    rescoreTests(tests, clusters);
    priority_queue<Test> pq;
    for (vector<Test>::iterator it = tests.begin(); it != tests.end(); ++it) {
      pq.push(*it);
    }

    Test t = pq.top();
    bool outcome = realization.second.first.IsKey(t.second);
    for (size_t i = 0; i < clusters.size(); ++i)
      clusters[i].ruleOutHypothesis(t, outcome);

    tests.erase(find(tests.begin(), tests.end(), t));

    clustersLeft = 0;
    for (size_t i = 0; i < clusters.size(); ++i) {
      if (clusters[i].countTotalHypothesis()) {
        clustersLeft++;
        clusters[i].printState();
      }
    }

    cout << endl;
  }

  cout << "True hypothesis cluster: " << index << endl;

  return 0;
}
