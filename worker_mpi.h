/**
 * Slave Job for MPI.
 * Author (Victor Carbune): vcarbune@ethz.ch
 */

#include "mpi.h"
#include "utils.h"

#include "snap/snap-core/Snap.h"

#ifndef WORKER_MPI_H_
#define WORKER_MPI_H_

void startWorker(PUNGraph, SimConfig);

#endif
