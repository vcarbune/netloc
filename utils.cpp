/**
 * Various utility functions.
 */

#include <fstream>
#include <iostream>
#include "utils.h"

using namespace std;

const char* algorithmTypeToString(AlgorithmType obj) {
  switch(obj) {
    case EC2:
      return "ec2";
    case GBS:
      return "gbs";
    case VOI:
      return "voi";
    case RANDOM:
      return "rand";
    case EPFL_ML:
      return "epfl";
  }

  return "huh";
}

SimConfig& SimConfig::operator++() {
  cluster.size *= 2;
  return *this;
}

ostream& operator<<(ostream& os, const SimConfig& config)
{
  os << config.cluster.size;
  return os;
}


SimConfig SimConfig::getSimConfigFromEnv(int argc, char *argv[], bool silent)
{
  Env = TEnv(argc, argv, silent ? NULL : TNotify::StdNotify);
  if (!silent)
    Env.PrepArgs(TStr::Fmt("NetLoc. Build: %s, %s. Time: %s",
          __TIME__, __DATE__, TExeTm::GetCurTm()));

  /* EC2 Model Parameters */
  const TInt paramClusterSize = Env.GetIfArgPrefixInt(
      "-c=", 500, "Cluster size");
  const double paramEps = Env.GetIfArgPrefixFlt(
      "-eps=", 0.01, "EPS Noisy Measurement Probability");
  const double paramCascadeSize = Env.GetIfArgPrefixFlt(
      "-bound=", 1.0, "Ground truth size");
  const double paramBeta = Env.GetIfArgPrefixFlt(
      "-beta=", 0.06, "Activation probability on edges)");
  const TInt paramSteps = Env.GetIfArgPrefixInt(
      "-steps=", 1, "Number of simulation steps.");
  const TInt paramObjType = Env.GetIfArgPrefixInt(
      "-obj=", -1, "Objective Type: All(-1), "
                  "EC2 - 0, GBS - 1, VoI - 2, Random - 3, EPFL - 4");
  const TInt paramSimulations = Env.GetIfArgPrefixInt(
      "-sim=", 5, "Simulation steps for cascades");

  /* The infection model of the ground truth */
  const TInt paramInfectionType = Env.GetIfArgPrefixInt(
      "-inf=", 0, "Infection Type: 0 - Default, 1 - Gaussian");

  /* Logging and Input */
  const TStr paramOutputLog = Env.GetIfArgPrefixStr(
      "-dump=", "dump.log", "File where to dump the output");
  const TStr networkInFile = Env.GetIfArgPrefixStr(
      "-netin=", "", "File where to load the network from");
  const TStr paramGroundTruthFile = Env.GetIfArgPrefixStr(
      "-truth=", "", "The ground truth to search for");
  const TInt paramGroundTruths = Env.GetIfArgPrefixInt(
      "-truths=", 20, "The total number of ground truths");

  SimConfig config;

  config.cluster.size = paramClusterSize.Val;
  config.cluster.simulations = paramSimulations;
  config.cluster.beta = paramBeta;
  config.cluster.bound = paramCascadeSize;

  config.eps = paramEps;

  config.steps = paramSteps.Val;
  config.logfile = paramOutputLog;
  config.netinFile = networkInFile;

  /* Fixed configuration for Gaussian generated */
  config.cluster.miu = 8;
  config.cluster.sigma = 2;

  /* Keep clusters from one iteration to another. */
  config.cluster.keep = true;

  config.setObjType(static_cast<AlgorithmType>(paramObjType.Val));
  config.infType = static_cast<InfectionType>(paramInfectionType.Val);
  if (config.infType == GAUSSIAN)
    config.cluster.beta = 1.00 / config.cluster.miu;

  config.groundTruthFile = paramGroundTruthFile;
  config.groundTruths = paramGroundTruths;
  if (!paramGroundTruthFile.Empty())
    config.groundTruths = 1;

  return config;
}

void SimConfig::setObjType(AlgorithmType objType)
{
  this->objType = objType;
  switch (this->objType) {
    case EC2:
      objSums = EC2_SUMS;
      break;
    case GBS:
      objSums = GBS_SUMS;
      break;
    case VOI:
      objSums = VOI_SUMS;
      break;
    case RANDOM:
      objSums = RANDOM_SUMS;
      break;
    case EPFL_ML:
      objSums = REGULAR_SUMS;
  }
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
  m_config.ndcgN = static_cast<int>(0.02 * m_network->GetNodes());
}
