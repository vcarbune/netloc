#include <iostream>

#include "master_mpi.h"
#include "worker_mpi.h"

using namespace std;

int main(int argc, char *argv[]) {
  MPI::Init(argc, argv);

  srand(time(NULL));

  SimConfig config = SimConfig::getSimConfigFromEnv(argc, argv);
  config.mpi.nodes = MPI::COMM_WORLD.Get_size();
  config.mpi.rank = MPI::COMM_WORLD.Get_rank();
  MPI::COMM_WORLD.Barrier();

  // The network is read directly from a file previously generated.
  if (config.netinFile.Empty()) {
    MPI::Finalize();
    cerr << "Network file MUST be provided!" << endl;
    return 0;
  }

  // Read the network.
  PUNGraph network;
  { TFIn FIn(config.netinFile); network = TUNGraph::Load(FIn); }
  config.nodes = network->GetNodes();
  cout << config.mpi.rank << ": " <<
      "Loaded network (N=" << config.nodes << ") from file..." << endl;

  if (config.mpi.rank == MPI_MASTER)
    startMaster(network, config);
  else
    startWorker(network, config);

  MPI::Finalize();

  return 0;
}
