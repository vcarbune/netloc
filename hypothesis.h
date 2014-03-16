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
#undef min
#undef max

class GraphHypothesis {
  public:
    static GraphHypothesis generateHypothesis(
        PUNGraph, int, double, double, std::vector<int>* = NULL);
    static GraphHypothesis generateHypothesisUsingGaussianModel(
        PUNGraph, int, double, double, double, std::vector<int>* = NULL);
    static GraphHypothesis readHypothesisFromFile(const char*);

    void writeHypothesisToFile(const char*);
    bool isConsistentWithTest(const GraphTest& test) const;

    double getInfectionTime(int) const;
    bool getTestOutcome(const GraphTest&) const;

    unsigned int getSize() const { return m_infectionHash.size(); }
    short unsigned int getSource() const { return m_sourceId; }

    double weight;

  private:
    GraphHypothesis(unsigned int, std::unordered_map<int, double>&);

    std::unordered_map<int, double> m_infectionHash;
    unsigned int m_sourceId;
};

class GraphHypothesisCluster {
  public:
    static GraphHypothesisCluster generateHypothesisCluster(
        PUNGraph, int, double, double, int, int);
    void updateMassWithTest(const GraphTest&);
    std::pair<double, double> computeMassWithTest(const GraphTest&) const;
    double getMass() const { return m_weight; }

    int getTotalHypothesis() const { return m_hypothesis.size(); }
    int getNodeCount(int nodeId) const { return m_nodeCount[nodeId]; }

    int getSource() const { return m_sourceId; }
    double getWeight() const { return m_weight; }
    void normalizeWeight(double mass) { m_weight /= mass; }

    void resetWeight(double);

    bool operator< (const GraphHypothesisCluster& o) const {
      return m_weight > o.m_weight;
    }

  private:
    GraphHypothesisCluster(PUNGraph, int, double);

    PUNGraph m_network;
    int m_sourceId;
    double m_weight;

    std::vector<GraphHypothesis> m_hypothesis;

    // Counts the number of times a particular node appears in all the generated
    // hypothesis. This can be used to speed up the computation of the probability
    // of a test outcome to be positive.
    std::vector<int> m_nodeCount;
};

#endif // HYPOTHESIS_H_
