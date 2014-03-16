/**
 * Various utility functions.
 */

#include "utils.h"

using namespace std;

SimConfig& SimConfig::operator++() {
  clusterSize *= 2;
  return *this;
}

SimConfig SimConfig::getSimConfigFromEnv(int argc, char *argv[], bool silent)
{
  Env = TEnv(argc, argv, silent ? NULL : TNotify::StdNotify);
  if (!silent)
    Env.PrepArgs(TStr::Fmt("NetLoc. Build: %s, %s. Time: %s",
          __TIME__, __DATE__, TExeTm::GetCurTm()));

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
  const double paramTestThreshold = Env.GetIfArgPrefixFlt(
      "-testthr=", 0.25, "Tests Threshold (%)");
  const TInt paramObjType = Env.GetIfArgPrefixInt(
      "-obj=", 0, "Objective Type: "
                  "EC2 - 0, GBS - 1, VoI - 2");
  const TStr paramGroundTruth = Env.GetIfArgPrefixStr(
      "-truth=", "", "The ground truth to search for");
  const TInt paramGroundTruths = Env.GetIfArgPrefixInt(
      "-truths=", 20, "The total number of ground truths");

  SimConfig config;

  config.clusterSize = paramClusterSize.Val;
  config.steps = paramSteps.Val;
  config.cascadeBound = paramCascadeSize;
  config.beta = paramBeta;
  config.logfile = dumpFile;
  config.testThreshold = paramTestThreshold;
  config.netinFile = networkInFile;

  config.objType = static_cast<AlgorithmType>(paramObjType.Val);
  switch (config.objType) {
    case EC2:
      config.objSums = EC2_SUMS;
      break;
    case GBS:
      config.objSums = GBS_SUMS;
      break;
    case VOI:
      config.objSums = VOI_SUMS;
      break;
    case RANDOM:
      config.objSums = RANDOM_SUMS;
      break;
    case EPFL_ML:
      config.objSums = REGULAR_SUMS;
  }

  config.groundTruth = paramGroundTruth;
  config.groundTruths = paramGroundTruths;
  if (!paramGroundTruth.Empty())
    config.groundTruths = 1;

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
}
