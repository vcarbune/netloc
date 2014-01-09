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
    GraphHypothesis(TIntH, double);

    bool isConsistentWithTest(const GraphTest& test) const;

    int getInfectionTime(int) const;
    bool getTestOutcome(const GraphTest&) const;

    int getSize() const { return m_infectionTimeHash.Len(); }

  private:
    TIntH m_infectionTimeHash;
};

class GraphHypothesisCluster {
  public:
    GraphHypothesisCluster(PUNGraph, int, int);
    GraphHypothesisCluster(PUNGraph, int, int, double, double);

    int countHypothesisConsistentWithTest (const GraphTest&) const;
    int countHypothesisAvailable() const { return m_hypothesis.size(); }

    void removeHypothesisInconsistentWithTest(const GraphTest&);
    GraphHypothesis getRandomHypothesis() const;
    GraphHypothesis generateHypothesis(bool = false) const;

    int getSource() const { return m_sourceId; }
    virtual void printState();

    void setWeight(double);

    // Not used yet.
    void setHopsFromSource(int);
    int getHopsFromSource() const { return m_hops; }

    bool operator< (const GraphHypothesisCluster& o) const {
      return countHypothesisAvailable() < o.countHypothesisAvailable();
    }
  private:
    void generateHypothesisCluster(int);

    PUNGraph m_network;
    vector<GraphHypothesis> m_hypothesis;
    int m_sourceId;
    double m_beta;
    double m_size;
    int m_hops;
    double m_weight;
};

#endif // HYPOTHESIS_H_
