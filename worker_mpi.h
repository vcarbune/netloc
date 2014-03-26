/**
 * Worker Job for MPI.
 * Author (Victor Carbune): vcarbune@ethz.ch
 */

#include "mpi.h"
#include "hypothesis.h"
#include "test.h"
#include "utils.h"

#include "snap/snap-core/Snap.h"

#ifndef WORKER_MPI_H_
#define WORKER_MPI_H_

class WorkerNode : public MPINode {
  public:
    WorkerNode(SimConfig);
    virtual void run();

  private:
    void runWithCurrentConfig();

    /* Initialize */
    void initializeClusters();
    void computeTestPriors();
    void initializeTestHeap();

    /* Simulate */
    void reset();
    void simulate();

    /* Computational Methods */
    void computeCurrentWeight();
    void recomputePartialTestScores();

    std::vector<GraphHypothesisCluster> m_clusters;
    std::vector<std::pair<double, int>> m_previousTests;
    std::pair<int, int> m_nodeRange;
};

#endif
