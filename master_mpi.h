/**
 * Master Job for MPI.
 * Author (Victor Carbune): vcarbune@ethz.ch
 */

#include "mpi.h"
#include "utils.h"

#include "snap/snap-core/Snap.h"

#ifndef MASTER_MPI_H_
#define MASTER_MPI_H_

void startMaster(PUNGraph network, SimConfig config);

#endif
