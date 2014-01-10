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
    GraphHypothesis(TIntH, double);

    bool isConsistentWithTest(const GraphTest& test) const;

    bool isMarkedAsInconsistent() const { return m_isMarkedAsInconsistent; }
    void markAsInconsistent() { m_isMarkedAsInconsistent = true; }

    int getInfectionTime(int) const;
    bool getTestOutcome(const GraphTest&) const;

    int getSize() const { return m_infectionTimeHash.Len(); }

  private:
    TIntH m_infectionTimeHash;
    bool m_isMarkedAsInconsistent;
};

class GraphHypothesisCluster {
  public:
    GraphHypothesisCluster(PUNGraph, int, int, double, double, double);

    void markInconsistentHypothesis(const GraphTest&);
    int countConsistentHypothesis() const { return m_crtConsistentHypothesis; }
    void countConsistentHypothesis(const GraphTest&, int*, int*) const;

    void removeHypothesisInconsistentWithTest(const GraphTest&);
    GraphHypothesis getRandomHypothesis() const;
    GraphHypothesis generateHypothesis(bool = false) const;

    int getSource() const { return m_sourceId; }

    void updateWeight(const vector<GraphTest>&);
    double getWeight() { return m_weight; }
    void setWeight(double);

    // Not used yet.
    void setHopsFromSource(int);
    int getHopsFromSource() const { return m_hops; }

    bool operator< (const GraphHypothesisCluster& o) const {
      return m_weight > o.m_weight;
      // return countHypothesisAvailable() < o.countHypothesisAvailable();
    }
  private:
    void generateHypothesisCluster(int);
    void updateConsistentHypothesisCount();

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
