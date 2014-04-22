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
#include "utils.h"

#include "snap/snap-core/Snap.h"
#undef min
#undef max


class GraphHypothesis : public Counter<GraphHypothesis> {
  public:
    /* Generators */
    static GraphHypothesis generateHypothesis(
        PUNGraph, int, const HypothesisClusterConfig, bool = true);
    static GraphHypothesis generateHypothesisUsingGaussianModel(
        PUNGraph, int, const HypothesisClusterConfig, bool = true);

    /* Input/Output */
    static GraphHypothesis readHypothesisFromFile(PUNGraph, const char*);
    void writeHypothesisToFile(const char*);

    /* Tests */
    bool isConsistentWithTest(
        const GraphTest&, const std::vector<std::pair<double, int>>&) const;
    double getInfectionTime(int) const;
    bool getTestOutcome(const GraphTest&) const;

    /* Properties */
    int getSize() const { return m_infectionHash.size(); }
    int getSource() const { return m_sourceId; }

    double weight;

  private:
    GraphHypothesis(unsigned int, std::unordered_map<int, double>&, double);
    bool isConsistentWithPreviousTests(
        const GraphTest&, const std::vector<std::pair<double, int>>&) const;

    std::unordered_map<int, double> m_infectionHash;
    int m_sourceId;
    double m_maxInfectionTime;
};

class GraphHypothesisCluster : public Counter<GraphHypothesisCluster>{
  public:
    static GraphHypothesisCluster generateHypothesisCluster(
        PUNGraph, int, double, const SimConfig& config);

    /* Tests */
    void updateMassWithTest(const double,
        const GraphTest&, const std::vector<std::pair<double, int>>&);
    std::pair<double, double> computeMassWithTest(
        const double, const GraphTest&, bool,
        const std::vector<std::pair<double, int>>&,
        const std::vector<double>&) const;

    /* Properties */
    int getSource() const { return m_sourceId; }
    double getWeight() const { return m_weight; }
    double getAverageHypothesisSize() const;

    void resetWeight(double);
    void multiplyWeights(double);

    /* Time Info */
    void collectInfectionTimes(int, std::vector<double>&) const;

    bool operator< (const GraphHypothesisCluster& o) const {
      return m_weight > o.m_weight;
    }

  private:
    GraphHypothesisCluster(PUNGraph, int, double);

    PUNGraph m_network;
    int m_sourceId;
    double m_weight;

    std::vector<GraphHypothesis> m_hypothesis;
};

template<class T>
int Counter<T>::count = 0;

#endif
