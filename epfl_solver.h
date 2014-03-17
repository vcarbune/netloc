/**
 * EPFL Solver class contains the implementation of the Maximum-Likelihood
 * solution of the problem, according to the paper
 *
 * "Locating the Source of Diffusion in Large-Scale Networks"
 * Link:  http://infoscience.epfl.ch/record/181841/files/PRL.pdf
 *
 * Within the present thesis, it is used as a comparison method.
 * Author (Victor Carbune): vcarbune@ethz.ch
 */
#include "utils.h"

#include "snap/snap-core/Snap.h"
#undef min
#undef max

#ifndef EPFL_SOLVER_H_
#define EPFL_SOLVER_H_

class EPFLSolver {
  public:
    EPFLSolver(PUNGraph, SimConfig);
    result_t solve(const GraphHypothesis&);

  private:
    std::vector<int> m_observerNodes;

    PUNGraph m_network;
    SimConfig m_config;
};

#endif
