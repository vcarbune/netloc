/**
 * Data structure used to represent a hypothesis cluster.
 * Author: Victor Carbune (vcarbune@ethz.ch)
 */

#ifndef HYPOTHESIS_H_
#define HYPOTHESIS_H_

#include <map>
#include <vector>

#include "test.h"
#include "snap/snap-core/Snap.h"

using namespace std;

/**
 * One hypothesis is a pair of the form
 *    <weight, <infectionTimeHash, actualCascade>>
 */
class GraphHypothesis {
  public:
    GraphHypothesis(PNGraph, TIntH, double);

    int getInfectionTime(int) const;
    bool getTestOutcome(const GraphTest&) const;

    void setWeight(double);

  private:
    double m_weight;
    TIntH m_infectionTimeHash;
    PNGraph m_cascade;
};

class GraphHypothesisCluster {
  public:
    GraphHypothesisCluster(PUNGraph, int, int);

    int countHypothesisConsistentWithTest (const GraphTest&) const;
    int countHypothesisAvailable() const { return m_hypothesis.size(); }

    void removeHypothesisInconsistentWithTest(const GraphTest&);
    GraphHypothesis getRandomHypothesis();

    int getSource() { return m_sourceId; }
    virtual void printState();

  private:
    void generateHypothesisCluster(int);
    GraphHypothesis generateHypothesis();

    PUNGraph m_network;
    vector<GraphHypothesis> m_hypothesis;

    int m_sourceId;
};

#endif // HYPOTHESIS_H_
