/**
 * Various utility functions.
 */

#include "utils.h"

using namespace std;

SimConfig& SimConfig::operator++() {
  clusterSize *= 2;
  return *this;
}

SimConfig SimConfig::getSimConfigFromEnv(int argc, char *argv[])
{
  Env = TEnv(argc, argv, TNotify::StdNotify);
  Env.PrepArgs(TStr::Fmt("NetLoc. Build: %s, %s. Time: %s",
        __TIME__, __DATE__, TExeTm::GetCurTm()), 1, true);

  const TInt paramNodes = Env.GetIfArgPrefixInt(
      "-n=", 500, "Network size");
  const TInt paramClusterSize = Env.GetIfArgPrefixInt(
      "-c=", 500, "Cluster size");
  const double paramCascadeSize = Env.GetIfArgPrefixFlt(
      "-s=", 0.4, "Ground truth size");
  const double paramBeta = Env.GetIfArgPrefixFlt(
      "-b=", 0.06, "Activation probability on edges)");
  const TInt paramSteps = Env.GetIfArgPrefixInt(
      "-steps=", 1, "Number of simulation steps.");
  const TStr dumpFile = Env.GetIfArgPrefixStr(
      "-dump=", "dump.log", "File where to dump the output");
  const TStr networkInFile = Env.GetIfArgPrefixStr(
      "-netin=", "", "File where to load the network from");
  const TBool paramLazy = Env.GetIfArgPrefixBool(
      "-lazy=", false, "Lazy evaluation");
  const TInt paramNetworkType = Env.GetIfArgPrefixInt(
      "-type=", 0, "Network Type: "
                   "ForestFire - 0, Barabasi-Albert - 1, Erdos-Renyi - 2");
  const TInt paramOutputType = Env.GetIfArgPrefixInt(
      "-output=", 0, "Output Type: Tests - 0, Probability - 1, Diff - 2");
  const double paramTestThreshold = Env.GetIfArgPrefixFlt(
      "-testthr=", 0.25, "Tests Threshold (%)");
  const TInt paramObjType = Env.GetIfArgPrefixInt(
      "-obj=", 0, "Objective Type: "
                  "EC2 - 0, GBS - 1, VoI - 2");
  const TStr paramGroundTruth = Env.GetIfArgPrefixStr(
      "-truth=", "", "The ground truth to search for");
  const TInt paramGroundTruths = Env.GetIfArgPrefixInt(
      "-truths=", 20, "The total number of ground truths");
  const TBool paramEpfl = Env.GetIfArgPrefixInt(
      "-epfl=", false, "EPFL solution?");

  SimConfig config;

  config.nodes = paramNodes.Val;
  config.clusterSize = paramClusterSize.Val;
  config.steps = paramSteps.Val;
  config.cascadeBound = paramCascadeSize;
  config.beta = paramBeta;
  config.logfile = dumpFile;
  config.lazy = paramLazy;
  config.networkType = paramNetworkType;
  config.outputType = paramOutputType;
  config.testThreshold = paramTestThreshold;
  config.netinFile = networkInFile;
  config.objType = paramObjType;
  switch (paramObjType) {
    case 0:
      config.objSums = EC2_SUMS;
      break;
    case 3:
      config.objSums = RANDOM_SUMS;
      break;
    default:
      config.objSums = REGULAR_SUMS;
  }

  config.groundTruth = paramGroundTruth;
  config.groundTruths = paramGroundTruths;
  if (!paramGroundTruth.Empty())
    config.groundTruths = 1;

  config.epflSolver = paramEpfl;
  config.epflMiu = 8;
  config.epflSigma = 2;

  return config;
}

ostream& operator<<(ostream& os, const SimConfig& config)
{
  os << config.clusterSize;
  return os;
}

MPINode::MPINode(SimConfig config)
  : m_config(config)
{
  readNetwork();
  m_testsPrior.resize(m_config.nodes);
}

void MPINode::readNetwork()
{
  if (m_config.netinFile.Empty()) {
    cerr << "Network file MUST be provided!" << endl;
    return;
  }

  TFIn FIn(m_config.netinFile);
  m_network = TUNGraph::Load(FIn);

  m_config.nodes = m_network->GetNodes();
  cout << m_config.mpi.rank << ": " <<
      "Loaded network (N=" << m_config.nodes << ") from file..." << endl;
}

void generateNetwork(PUNGraph *network, SimConfig& config) {
  if (!config.netinFile.Empty()) {
    { TFIn FIn(config.netinFile); *network = TUNGraph::Load(FIn); }
    config.nodes = (*network)->GetNodes();
    cout << "Loaded network (N=" << config.nodes <<
        ") from file..." << endl;
    return;
  }
  switch (config.networkType) {
    case 0:
      *network = TSnap::ConvertGraph<PUNGraph, PNGraph>(TSnap::GenForestFire(config.nodes, 0.35, 0.32));
      break;
    case 1:
      *network = TSnap::GenPrefAttach(config.nodes, 3);
      break;
    case 2:
      *network = TSnap::GenRndGnm<PUNGraph>(
          config.nodes, config.nodes * log(config.nodes));
      break;
  }

  cout << "Generated network is connected: " << TSnap::GetMxWccSz(*network) << endl;
}

void generateClusters(vector<GraphHypothesisCluster> *clusters,
                      const PUNGraph network,
                      const SimConfig& config)
{
  // Generate all possible hypothesis clusters that we want to search through.
  clusters->clear();
  for (int source = 0; source < network->GetNodes(); source++) {
    clusters->push_back(GraphHypothesisCluster::generateHypothesisCluster(
        network, source, 1,
        config.beta, config.cascadeBound, config.clusterSize));
  }
}

void generateTests(vector<GraphTest> *tests,
                   const PUNGraph network)
{
  // Generate all tests that are enough to differentiate between hypothesis.
  tests->clear();
  for (int node = 0; node < network->GetNodes(); ++node)
    tests->push_back(GraphTest(node));
}
