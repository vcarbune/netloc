/**
 * Various utility functions.
 *
 * Author (Victor Carbune): vcarbune@ethz.ch
 */

#include <vector>

#include "test.h"

#include "snap/snap-core/Snap.h"
#undef min
#undef max

#ifndef UTILS_H_
#define UTILS_H_

#define MPI_MASTER 0
#define DBG 0

/* MPI DEFINES */
#define RANDOM_SUMS       0
#define REGULAR_SUMS      2
#define EC2_SUMS          4
#define GBS_SUMS          2
#define VOI_SUMS          1

#define POSITIVE_SUM      0
#define NEGATIVE_SUM      1
#define POSITIVE_DIAG_SUM 2
#define NEGATIVE_DIAG_SUM 3

typedef std::pair<int, std::vector<double>> result_t;

struct MPIClusterConfig {
  int nodes;  // total number of nodes initialized and used by MPI.
  int rank;   // rank of the current node.
};

struct HypothesisClusterConfig {
  int size;         // number of hypothesis per node.
  int simulations;  // number of simulations for hypothesis in cluster.
  double beta;      // edge activation probability.
  double bound;     // percentage bound of the total network.
  double miu;       // for gaussian-generated cascades.
  double sigma;
};

enum AlgorithmType {
  EC2 = 0,    // Adaptive     Edge Cutting Equivalence Class
  GBS,        // Adaptive     Generalized Binary Search
  VOI,        // Adaptive     Value of Information
  RANDOM,     // Non-adaptive Random Selection
  EPFL_ML     // Non-adaptive EPFL Approach (using Maximum Likelihood)
};

enum InfectionType {
  BETA = 0,   // SSI          Default model
  GAUSSIAN    // Gaussian Distributed Delay (joint distribution for theta)
};

const char* algorithmTypeToString(AlgorithmType);

class SimConfig {
  public:
    static SimConfig getSimConfigFromEnv(int argc, char *argv[], bool = true);

    void setObjType(AlgorithmType);

    int nodes;
    int steps;
    TStr logfile;
    double testThreshold; // stop ec2 after a percentage of tests have run.

    TStr netinFile;

    AlgorithmType objType;
    int objSums;

    InfectionType infType;

    TStr groundTruthFile;
    int groundTruths;

    // Noisy Model Re-weighing (tolerance to noisy measurements).
    double eps;
    int ndcgN;

    MPIClusterConfig mpi;
    HypothesisClusterConfig cluster;

    // Parameters changing from one step of the simulation to the other.
    SimConfig& operator++();

  private:
    SimConfig() {};
};
std::ostream& operator<<(std::ostream& os, const SimConfig& config);

// Base class for MPI Nodes.
class MPINode {
  public:
    MPINode(SimConfig);
    virtual void run() = 0;

  protected:
    void readNetwork();

    PUNGraph m_network;
    SimConfig m_config;

    std::vector<double> m_testsPrior;
};

#endif // UTILS_H_
