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

class GraphHypothesis {
  public:
    GraphHypothesis(TIntH);

    bool isConsistentWithTest(const GraphTest& test) const;

    int getInfectionTime(int) const;
    bool getTestOutcome(const GraphTest&) const;

    int getSize() const { return m_infectionTimeHash.Len(); }

    double weight;

  private:
    TIntH m_infectionTimeHash;
};

class GraphHypothesisCluster {
  public:
    GraphHypothesisCluster(PUNGraph, int, int, double, double, double);

    void updateMassWithTest(const GraphTest&);
    double computeMassWithTest(const GraphTest&) const;
    double getMass() const { return m_weight; }

    // int countConsistentHypothesis() const { return m_crtConsistentHypothesis; }
    void countConsistentHypothesis(const GraphTest&, int*, int*) const;

    void removeHypothesisInconsistentWithTest(const GraphTest&);
    GraphHypothesis getRandomHypothesis() const;
    GraphHypothesis generateHypothesis(bool = false) const;

    int getSource() const { return m_sourceId; }
    double getWeight() { return m_weight; }

    bool operator< (const GraphHypothesisCluster& o) const {
      return m_weight > o.m_weight;
    }
  private:
    void generateHypothesisCluster(int);
    // void updateConsistentHypothesisCount();

    PUNGraph m_network;
    vector<GraphHypothesis> m_hypothesis;
    int m_sourceId;
    double m_beta;
    double m_size;
    int m_hops;
    double m_weight;

    int m_crtConsistentHypothesis;
};

#endif // HYPOTHESIS_H_
