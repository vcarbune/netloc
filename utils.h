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
  NodeVar,  // nodes: G(V, VlogV)
  EdgeVar,  // edge: G(sqrt(E), E)
  BetaVar,  // beta: [0.1, 0.5]
  HypothesisVar,  // hypothesis per cluster: [1, 15]
  CascBoundVar    // max cascade size: [0.4, 0.8]
};

struct SimConfig {
  SimConfig(SimulationType type)
    : nodes(100)
    , edges(200)
    , beta(0.1)
    , clusterSize(10)
    , cascadeBound(0.4)
    , m_type(type)
  {
    if (type == BetaVar || type == HypothesisVar || type == CascBoundVar) {
      // The base network on which the simulation is done is fixed throughout
      // so in the initial configuration we should have it enough.
      nodes = 500;
      edges = 1500;
    }

    if (type == HypothesisVar) {
      clusterSize = 5;
    }
  }

  SimConfig& operator++() {
    switch (m_type) {
      case NodeVar:
        nodes += 50;
        edges = nodes * log(nodes);
        break;
      case EdgeVar:
        edges += 50;
        nodes = pow(edges, 1.1);
        break;
      case BetaVar:
        beta += 0.05;
        break;
      case HypothesisVar:
        clusterSize += 5;
        break;
      case CascBoundVar:
        cascadeBound += 0.05;
        break;
    }

    return *this;
  }

  double getSimParamValue() const {
    switch (m_type) {
      case NodeVar:
        return static_cast<double>(nodes);
      case EdgeVar:
        return static_cast<double>(edges);
      case BetaVar:
        return beta;
      case HypothesisVar:
        return static_cast<double>(clusterSize);
      case CascBoundVar:
        return static_cast<double>(cascadeBound);
    }

    return -1.0;
  }

  int nodes;
  int edges;
  int beta;
  int clusterSize;  // cluster size, number of hypothesis in a cluster
  int cascadeBound; // cascade size, percentage of the network

private:
  SimulationType m_type;
};

#endif // UTILS_H_
