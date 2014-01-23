/**
 * Various utility functions to do good stuff.
 *
 * Author (Victor Carbune): vcarbune@ethz.ch
 */

#ifndef UTILS_H_
#define UTILS_H_

#include "hypothesis.h"
#include "test.h"
#include "snap/snap-core/Snap.h"

#include <vector>

// TODO(vcarbune): currently unused, refactor or remove.
void plotGraphs(const PUNGraph& network, const std::vector<GraphTest>& tests,
                const GraphHypothesis& realization);

enum SimulationType {
  BetaVar = 0,  // beta: [0.1, 0.5]
  ClusterVar,  // hypothesis per cluster: [1, 15]
  CascBoundVar    // max cascade size: [0.4, 0.8]
};

struct SimConfig {
  SimConfig(SimulationType type)
    : nodes(500)
    , edges(500 * log(500))
    , clusterSize(100)
    , beta(0.05)
    , cascadeBound(0.3)
    , steps(1)
    , lazy(false)
    , networkType(1)
    , outputType(0)
    , testThreshold(0.25)
    , massThreshold(0.50)
    , m_type(type)
  {
    if (type == BetaVar || type == ClusterVar || type == CascBoundVar) {
      // The base network on which the simulation is done is fixed throughout
      // so in the initial configuration we should have it enough.
      nodes = 500;
      edges = 1500;
    }

    if (type == ClusterVar) {
      clusterSize = 5;
    }
  }

  SimConfig& operator+=(int inc) {
    while (inc)
      ++(*this), inc--;

    return *this;
  }

  SimConfig& operator++() {
    switch (m_type) {
      case BetaVar:
        beta += 0.05;
        break;
      case ClusterVar:
        clusterSize *= 2;
        break;
      case CascBoundVar:
        cascadeBound += 0.05;
        break;
    }

    return *this;
  }

  double getSimParamValue() const {
    switch (m_type) {
      case BetaVar:
        return beta;
      case ClusterVar:
        return static_cast<double>(clusterSize);
      case CascBoundVar:
        return static_cast<double>(cascadeBound);
    }

    return -1.0;
  }

  SimulationType getSimType() const { return m_type; }

  int nodes;
  int edges;
  int clusterSize;      // cluster size, number of hypothesis in a cluster
  double beta;          // edge selection probability
  double cascadeBound;  // cascade size, percentage of the network
  int steps;
  int topN;             // keep topN possible solutions
  TStr logfile;
  bool lazy;
  int networkType;
  int outputType;       // probability of source, or number of tests.
  double testThreshold; // stop ec2 after a percentage of tests have run.
  double massThreshold; // stop ec2 if there's a cluster with more mass.

private:
  SimulationType m_type;
};

SimConfig getSimConfigFromEnv(int argc, char *argv[]);

#endif // UTILS_H_
