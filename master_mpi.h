/**
 * Master Job for MPI.
 * Author (Victor Carbune): vcarbune@ethz.ch
 */

#include "mpi.h"
#include "test.h"
#include "hypothesis.h"
#include "utils.h"
#include "epfl_solver.h"

#include "snap/snap-core/Snap.h"

#ifndef MASTER_MPI_H_
#define MASTER_MPI_H_

class MasterNode : public MPINode {
  public:
    MasterNode(SimConfig);
    virtual void run();

  private:
    /* Initialize */
    void initializeGroundTruths();
    double computeCurrentWeight();

    void initializeTests();
    void initializeTestHeap();
    void initializeTestVector();

    /* Simulate */
    result_t simulate(int);
    result_t simulateAdaptivePolicy(int);
    result_t simulateEPFLPolicy(int);

    GraphTest selectRandomTest(std::vector<GraphTest>&);
    GraphTest selectNextTest(std::vector<GraphTest>&);

    void recomputeTestScoreEC2(GraphTest&, double, double*);
    void recomputeTestScoreGBS(GraphTest&, double, double*);
    void recomputeTestScoreVOI(GraphTest&);
    void recomputeTestScore(GraphTest&, double);

    /* Results */
    result_t identifyCluster(int, const TIntH&);
    void processResults(std::vector<result_t>*, SimConfig&);

    /* Variables */
    std::vector<GraphHypothesis> m_realizations;
    std::vector<TIntH> m_idToShortestPathsFromSource;

    std::vector<GraphTest> m_tests;
    std::vector<std::pair<double, int>> m_previousTests;

    EPFLSolver m_epflSolver;
};

#endif
