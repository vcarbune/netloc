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

enum ParameterVariationType {
  NodeVar,  // nodes: G(V, VlogV)
  EdgeVar,  // edge: G(sqrt(E), E)
  BetaVar,  // beta: [0.1, 0.5]
  HypothesisVar, // hypothesis per cluster: [1, 15]
  CascBoundVar // max cascade size: [0.4, 0.8]
};

#endif // UTILS_H_
