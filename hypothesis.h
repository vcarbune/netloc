/**
 * Data structure used to represent a hypothesis cluster.
 */

#ifndef HYPOTHESIS_H

#include <map>
#include <vector>

#include "test.h"
#include "snap/snap-core/Snap.h"

using namespace std;

/**
 * One hypothesis is a pair of the form
 *    <weight, <infectionTimeHash, actualCascade>>
 */
typedef pair<double, pair<TIntH, PNGraph>> Hypothesis;

class HypothesisCluster {
  public:
    HypothesisCluster(PUNGraph, int, int);

    void ruleOutHypothesis(Test);
    Hypothesis getRandomHypothesis();
    int getSource() { return m_sourceId; }

    void printState();

  private:
    void generateHypothesisCluster(int);
    Hypothesis generateHypothesis();

    PUNGraph m_network;
    vector<Hypothesis> m_hypothesis;

    int m_sourceId;
};

#endif // HYPOTHESIS_H
