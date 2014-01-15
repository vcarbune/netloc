/**
 * Data structure used to represent a hypothesis cluster.
 * Author: Victor Carbune (vcarbune@ethz.ch)
 */

#ifndef HYPOTHESIS_H_
#define HYPOTHESIS_H_

#include <map>
#include <vector>
#include <unordered_map>

#include "test.h"
#include "snap/snap-core/Snap.h"

using namespace std;

class GraphHypothesis {
  public:
    explicit GraphHypothesis(unordered_map<int, int>&);

    bool isConsistentWithTest(const GraphTest& test) const;

    int getInfectionTime(int) const;
    bool getTestOutcome(const GraphTest&) const;
    int getSize() const { return m_infectionHash.size(); }

    double weight;

  private:
    unordered_map<int, int> m_infectionHash;
};

class GraphHypothesisCluster {
  public:
    GraphHypothesisCluster(PUNGraph, int, int, double, double, double);

    void updateMassWithTest(const GraphTest&);
    pair<double, double> computeMassWithTest(const GraphTest&) const;
    double getMass() const { return m_weight; }

    int getTotalHypothesis() const { return m_hypothesis.size(); }
    int getNodeCount(int nodeId) const { return m_nodeCount[nodeId]; }

    GraphHypothesis getRandomHypothesis() const;
    GraphHypothesis generateHypothesis(double, double, vector<int>* = NULL) const;

    int getSource() const { return m_sourceId; }
    double getWeight() { return m_weight; }
    void normalizeWeight(double mass) { m_weight /= mass; }

    bool operator< (const GraphHypothesisCluster& o) const {
      return m_weight > o.m_weight;
    }
  private:
    void generateHypothesisCluster(double, double, int);

    PUNGraph m_network;
    vector<GraphHypothesis> m_hypothesis;
    int m_sourceId;
    double m_weight;

    // Counts the number of times a particular node appears in all the generated
    // hypothesis. This can be used to speed up the computation of the probability
    // of a test outcome to be positive.
    vector<int> m_nodeCount;
};

#endif // HYPOTHESIS_H_
