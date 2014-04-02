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
      return "epfl_high";
    case EPFL_EC2:
      return "epfl_ec2";
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
      "-eps=", 0.02, "EPS Noisy Measurement Probability");
  const double paramCascadeSize = Env.GetIfArgPrefixFlt(
      "-bound=", 1.0, "Ground truth size");
  const double paramArtifCascadeSize = Env.GetIfArgPrefixFlt(
      "-cbound=", 0.20, "Artif truths size");
  const double paramBeta = Env.GetIfArgPrefixFlt(
      "-beta=", 0.06, "Activation probability on edges)");
  const TInt paramSteps = Env.GetIfArgPrefixInt(
      "-steps=", 1, "Number of simulation steps.");
  const TInt paramObjType = Env.GetIfArgPrefixInt(
      "-obj=", -1, "Objective Type: All(-1), "
                  "EC2 - 0, GBS - 1, VoI - 2, Random - 3, EPFL - 4");


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
  config.cluster.beta = paramBeta;
  config.cluster.bound = paramCascadeSize;
  config.cluster.cbound = paramArtifCascadeSize;

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
    case EPFL_EC2:
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
  m_network->GetNIdV(m_nid);

  m_config.nodes = m_network->GetNodes();
  m_config.ndcgN = static_cast<int>(0.05 * m_network->GetNodes());
}

void MPINode::run()
{
  SimConfig initialConfig = m_config;
  vector<AlgorithmType> objectives;

  if (m_config.objType == EPFL_EC2) {
    // First run classical EPFL approach, using highest degree nodes.
    objectives.push_back(EPFL_ML);
    // First run EC2 and store observers.
    objectives.push_back(EC2);
    // Then run EPFL approach using the stored EC2 observers.
    objectives.push_back(EPFL_EC2);
  } else if (m_config.objType == -1) {
    // Run with all the objectives and compare everything.
    for (int obj = EC2; obj <= EPFL_EC2; ++obj)
      objectives.push_back(static_cast<AlgorithmType>(obj));
  } else {
    // Run only with the specified objective.
    objectives.push_back(static_cast<AlgorithmType>(m_config.objType));
  }

  for (AlgorithmType objType : objectives) {
    // Reset configuration.
    m_config = initialConfig;
    // Set objective properly.
    m_config.setObjType(static_cast<AlgorithmType>(objType));
    // Set logfile properly.
    m_config.logfile = TStr::Fmt("%s%d_%s.log", m_config.logfile.CStr(),
        m_config.nodes, algorithmTypeToString(static_cast<AlgorithmType>(objType)));

    runWithCurrentConfig();
  }
}
