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
    int getSize() const { return m_cascade->GetNodes(); }

  private:
    double m_weight;
    TIntH m_infectionTimeHash;
    PNGraph m_cascade;
};

class GraphHypothesisCluster {
  public:
    GraphHypothesisCluster(PUNGraph, int, int);
    GraphHypothesisCluster(PUNGraph, int, int, double, double);

    int countHypothesisConsistentWithTest (const GraphTest&) const;
    int countHypothesisAvailable() const { return m_hypothesis.size(); }

    void removeHypothesisInconsistentWithTest(const GraphTest&);
    GraphHypothesis getRandomHypothesis() const;
    GraphHypothesis generateHypothesis() const;

    int getSource() { return m_sourceId; }
    virtual void printState();

  private:
    void generateHypothesisCluster(int);

    PUNGraph m_network;
    vector<GraphHypothesis> m_hypothesis;

    int m_sourceId;
    double m_beta;
    double m_size;
};

#endif // HYPOTHESIS_H_
