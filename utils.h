/**
 * Various utility functions.
 *
 * Author (Victor Carbune): vcarbune@ethz.ch
 */

#include <vector>

#include "hypothesis.h"
#include "test.h"

#include "snap/snap-core/Snap.h"
#undef min
#undef max

#ifndef UTILS_H_
#define UTILS_H_

#define MPI_MASTER 0
#define DBG 0

/* Noisy Model Re-Weighing*/
#define EPS 0.05


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

struct MPIConfig {
  int rank;
  int nodes;
};

enum AlgorithmType {
  EC2 = 0,    // Adaptive     Edge Cutting Equivalence Class
  GBS,        // Adaptive     Generalized Binary Search
  VOI,        // Adaptive     Value of Information
  RANDOM,     // Non-adaptive Random Selection
  EPFL_ML     // Non-adaptive EPFL Approach (using Maximum Likelihood)
};

class SimConfig {
  public:
    static SimConfig getSimConfigFromEnv(int argc, char *argv[], bool = true);

    // Parameters changing from one step of the simulation to the other.
    SimConfig& operator++();

    int nodes;
    int clusterSize;      // cluster size, number of hypothesis in a cluster
    double beta;          // edge selection probability
    double cascadeBound;  // cascade size, percentage of the network
    int steps;
    TStr logfile;
    double testThreshold; // stop ec2 after a percentage of tests have run.

    TStr netinFile;

    AlgorithmType objType;
    int objSums;

    TStr groundTruthFile;
    int groundTruths;

    double epflMiu;
    double epflSigma;

    MPIConfig mpi;
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
