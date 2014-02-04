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
  Env.PrepArgs(TStr::Fmt("NetLoc. Build: %s, %s. Time: %s", __TIME__, __DATE__, TExeTm::GetCurTm()));

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
  const TInt paramKeepTopN = Env.GetIfArgPrefixInt(
      "-topN=", 1, "Keep topN solutions");
  const TStr dumpFile = Env.GetIfArgPrefixStr(
      "-dump=", "dump.log", "File where to dump the output");
  const TBool paramLazy = Env.GetIfArgPrefixBool(
      "-lazy=", false, "Lazy evaluation");
  const TInt paramNetworkType = Env.GetIfArgPrefixInt(
      "-type=", 0, "Network Type: "
                   "ForestFire - 0, Barabasi-Albert - 1, Erdos-Renyi - 2");
  const TInt paramOutputType = Env.GetIfArgPrefixInt(
      "-output=", 0, "Output Type: Tests - 0, Probability - 1, Diff - 2");
  const double paramTestThreshold = Env.GetIfArgPrefixFlt(
      "-testthr=", 0.25, "Tests Threshold (%)");
  const double paramMassThreshold = Env.GetIfArgPrefixFlt(
      "-massthr=", 1.00, "Mass Threshold (%)");

  SimConfig config;

  config.nodes = paramNodes.Val;
  config.clusterSize = paramClusterSize.Val;
  config.steps = paramSteps.Val;
  config.cascadeBound = paramCascadeSize;
  config.beta = paramBeta;
  config.logfile = dumpFile;
  config.topN = paramKeepTopN;
  config.lazy = paramLazy;
  config.networkType = paramNetworkType;
  config.outputType = paramOutputType;
  config.testThreshold = paramTestThreshold;
  config.massThreshold = paramMassThreshold;

  return config;
}

ostream& operator<<(ostream& os, const SimConfig& config)
{
  os << config.clusterSize;
  return os;
}

void generateNetwork(PUNGraph *network, const SimConfig& config) {
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
    clusters->push_back(GraphHypothesisCluster(
        network, source,
        config.clusterSize, config.beta, config.cascadeBound, 1));
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
