/**
 * Various utility functions.
 *
 * Author (Victor Carbune): vcarbune@ethz.ch
 */

#ifndef UTILS_H_
#define UTILS_H_

#include "hypothesis.h"
#include "test.h"
#include "snap/snap-core/Snap.h"

#include <vector>

#define MPI_MASTER 0

/* MPI DEFINES */
#define RANDOM_SUMS       0
#define REGULAR_SUMS      2
#define EC2_SUMS          4

#define POSITIVE_SUM      0
#define NEGATIVE_SUM      1
#define POSITIVE_DIAG_SUM 2
#define NEGATIVE_DIAG_SUM 3

struct MPIConfig {
  int rank;
  int nodes;
};

class SimConfig {
  public:
    static SimConfig getSimConfigFromEnv(int argc, char *argv[]);

    // Parameters changing from one step of the simulation to the other.
    SimConfig& operator++();

    int nodes;
    int networkType;
    int clusterSize;      // cluster size, number of hypothesis in a cluster
    double beta;          // edge selection probability
    double cascadeBound;  // cascade size, percentage of the network
    int steps;
    int topN;             // keep topN possible solutions
    TStr logfile;
    bool lazy;
    int outputType;       // probability of source or number of tests.
    double testThreshold; // stop ec2 after a percentage of tests have run.
    double massThreshold; // stop ec2 if there's a cluster with more mass.

    TStr netinFile;
    TStr netoutFile;

    int objType;
    int objSums;

    TStr groundTruth;
    int groundTruths;

    bool epflSolver;      // run EPFL algorithm on the given network.
    double epflMiu;
    double epflSigma;

    MPIConfig mpi;
  private:
    SimConfig() {};
};
ostream& operator<<(ostream& os, const SimConfig& config);

// Base class for MPI Nodes.
class MPINode {
  public:
    MPINode(SimConfig);
    virtual void run() = 0;

  protected:
    void readNetwork();

    PUNGraph m_network;
    SimConfig m_config;

    vector<double> m_testsPrior;
};

// Initialization
void generateNetwork(PUNGraph *network, SimConfig& config);
void generateClusters(vector<GraphHypothesisCluster> *clusters,
                      const PUNGraph network,
                      const SimConfig& config);
void generateTests(vector<GraphTest> *tests,
                   const PUNGraph network);

#endif // UTILS_H_
