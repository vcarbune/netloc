/**
 * Worker Job for MPI.
 * Author (Victor Carbune): vcarbune@ethz.ch
 */

#include "mpi.h"
#include "utils.h"

#include "snap/snap-core/Snap.h"

#ifndef WORKER_MPI_H_
#define WORKER_MPI_H_

class WorkerNode : public MPINode {
  public:
    WorkerNode(SimConfig);
    virtual void run();

  private:
    /* Initialize */
    void initializeClusters(int, int);
    void initializeTestHeap();

    /* Simulate */
    void simulate();

    /* Computational Methods */
    void computeCurrentMass();
    void recomputePartialTestScores();

    std::vector<GraphHypothesisCluster> m_clusters;
};

#endif
