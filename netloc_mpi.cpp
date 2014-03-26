#include <iostream>
#include <iomanip>

#include "master_mpi.h"
#include "worker_mpi.h"

using namespace std;

int main(int argc, char *argv[]) {
  MPI::Init(argc, argv);

  // Get initialization parameters.
  SimConfig config = SimConfig::getSimConfigFromEnv(argc, argv,
      MPI::COMM_WORLD.Get_rank() != MPI_MASTER);
  config.mpi.nodes = MPI::COMM_WORLD.Get_size();
  config.mpi.rank = MPI::COMM_WORLD.Get_rank();
  MPI::COMM_WORLD.Barrier();

  // cout << fixed << std::setprecision(3);

  // The network is read directly from a file previously generated.
  PUNGraph network;
  if (config.mpi.rank == MPI_MASTER) {
    MasterNode master(config);
    master.run();
  } else {
    WorkerNode worker(config);
    worker.run();
  }

  MPI::Finalize();

  return 0;
}
