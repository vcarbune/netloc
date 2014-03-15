/**
 * Master Job for MPI.
 * Author (Victor Carbune): vcarbune@ethz.ch
 */

#include "mpi.h"
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
    double computeCurrentMass();

    void initializeTests();
    void initializeTestHeap();
    void initializeTestVector();

    /* Simulate */
    result_t simulate(const GraphHypothesis&);
    result_t simulateAdaptivePolicy(const GraphHypothesis&);
    result_t simulateEPFLPolicy(const GraphHypothesis&);

    GraphTest selectRandomTest(vector<GraphTest>&);
    GraphTest selectNextTest(vector<GraphTest>&);

    void recomputeTestScoreEC2(GraphTest&, double, double*);
    void recomputeTestScoreGBS(GraphTest&, double, double*);
    void recomputeTestScoreVOI(GraphTest&);
    void recomputeTestScore(GraphTest&, double);

    /* Results */
    result_t identifyCluster(int);
    void processResults(vector<result_t>*, SimConfig&);

    /* Variables */
    vector<GraphHypothesis> m_realizations;
    vector<GraphTest> m_tests;

    EPFLSolver m_epflSolver;
};

#endif
