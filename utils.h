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
#define DBG 1

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

/* MASS RESET THRESHOLD */
#define RESET_SCORE_THRESHOLD         1
#define RESET_SCORE_FACTOR            1

typedef std::pair<int, std::vector<double>> result_t;

struct MPIClusterConfig {
  int nodes;  // total number of nodes initialized and used by MPI.
  int rank;   // rank of the current node.
};

struct HypothesisClusterConfig {
  int size;         // number of hypothesis per node.
  double beta;      // edge activation probability.
  double bound;     // percentage bound of the total network.
  double cbound;    // percentage bound of each hypothesis in cluster.
  double miu;       // for gaussian-generated cascades.
  double sigma;
  bool keep;
};

enum AlgorithmType {
  EC2 = 0,    // Adaptive     Edge Cutting Equivalence Class
  GBS,        // Adaptive     Generalized Binary Search
  RANDOM,     // Non-adaptive Random Selection
  EPFL_ML,    // Non-adaptive EPFL Approach (using Highest Degree observers)
  EPFL_EC2,   // Adaptive     EPFL Approach (using EC2 selected observers)
  VOI         // Adaptive     Value of Information
};

enum InfectionType {
  BETA = 0,   // SIS          Default model
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
    TStr netinFile;

    AlgorithmType objType;
    int objSums;

    InfectionType infType;

    TStr groundTruthFile;
    int groundTruths;

    // Noisy Model Re-weighing (tolerance to noisy measurements).
    double eps;
    double sampling = 0.05;
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
    void run();

  protected:
    void readNetwork();
    virtual void runWithCurrentConfig() = 0;

    PUNGraph m_network;
    SimConfig m_config;

    std::vector<double> m_testsPrior;
    TIntV m_nid;
};

#endif // UTILS_H_
