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

void plotGraphs(const PUNGraph& network, const std::vector<GraphTest>& tests,
                const GraphHypothesis& realization);

#endif // UTILS_H_
