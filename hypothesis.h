/**
 * Data structure used to represent a hypothesis cluster.
 */

#ifndef HYPOTHESIS_H

#include <map>
#include <vector>

#include "test.h"
#include "snap/snap-core/Snap.h"

using namespace std;

typedef map<int, int> Hypothesis;

class HypothesisCluster {
  public:
    HypothesisCluster(PUNGraph, int);

    void ruleOutHypothesis(Test);
    Hypothesis getRandomHypothesis();
  private:
    void generateCascades();

    PUNGraph m_network;
    int m_sourceId;
    vector<Hypothesis> m_hypothesis;
};

#endif // HYPOTHESIS_H
