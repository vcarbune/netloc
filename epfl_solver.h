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
#include "hypothesis.h"
#include "utils.h"

#include "snap/snap-core/Snap.h"
#undef min
#undef max

#ifndef EPFL_SOLVER_H_
#define EPFL_SOLVER_H_

// pair: (mean, variance)
typedef std::pair<double, double> gaussian_t;

class EPFLSolver {
  public:
    EPFLSolver(PUNGraph, SimConfig, double);
    result_t solve(const GraphHypothesis&,
        std::vector<std::pair<double, int>>&);

    void setObserverList(const std::vector<int>&, double);
    size_t countObservers() { return m_observerNodes.size(); }

  private:
    gaussian_t computeGaussianMoments(
        const std::vector<std::pair<double, int>>&);

    std::vector<int> m_observerNodes;
    gaussian_t m_gaussianParams;
    PUNGraph m_network;
    SimConfig m_config;
};

#endif
