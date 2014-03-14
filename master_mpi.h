/**
 * Master Job for MPI.
 * Author (Victor Carbune): vcarbune@ethz.ch
 */

#include "mpi.h"
#include "utils.h"

#include "snap/snap-core/Snap.h"

#ifndef MASTER_MPI_H_
#define MASTER_MPI_H_

typedef pair<int, vector<double>> result_t;

class MasterNode : public MPINode {
  public:
    MasterNode(SimConfig);
    virtual void run();

  private:
    /* Initialize */
    void initializeGroundTruths();
    double computeCurrentMass();
    void rebuildTestHeap();

    /* Simulate */
    result_t simulate(const GraphHypothesis&);
    GraphTest selectNextTest(vector<GraphTest>&);
    void recomputeTestScore(GraphTest&, double);

    /* Results */
    result_t identifyCluster(int);
    void processResults(vector<result_t>*, SimConfig&);

    vector<GraphHypothesis> m_realizations;
    vector<GraphTest> m_tests;
};

#endif
