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

  private:
    virtual void runWithCurrentConfig();

    /* Initialize */
    void initializeGroundTruths();
    double computeCurrentWeight(double*);

    void reset();
    void initializeNodeInfectionTimeMap();
    void initializeTests();
    void initializeTestHeap();
    void initializeTestVector();

    /* Simulate */
    std::vector<result_t> simulate(int);
    std::vector<result_t> simulateAdaptivePolicy(int);
    std::vector<result_t> simulateEPFLPolicy(int);

    GraphTest selectRandomTest(std::vector<GraphTest>&);
    GraphTest selectVOITest(std::vector<GraphTest>&, double);
    GraphTest selectNextTest(std::vector<GraphTest>&);

    void recomputeTestScoreEC2(GraphTest&, double, double*);
    void recomputeTestScoreGBS(GraphTest&, double, double*);
    void recomputeTestScoreVOI(GraphTest&);
    void recomputeTestScore(GraphTest&, double, double);

    /* Results */
    double computeNDCG(const std::vector<std::pair<double, int>>&, int) const;
    result_t identifyCluster(int);
    void processResults(std::vector<std::vector<result_t>>*, SimConfig&);

    /* Variables */
    std::vector<GraphHypothesis> m_realizations;
    std::vector<TIntH> m_idToShortestPathsFromSource;

    std::vector<GraphTest> m_tests;
    std::vector<std::pair<double, int>> m_previousTests;

    EPFLSolver m_epflSolver;

    /* m_ec2observers[TRUTH][PCNT]; */
    std::vector<std::vector<int>> m_ec2observers;
    std::vector<std::vector<double>> m_nodeInfectionTimes;
};

#endif
