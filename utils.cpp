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
    case EC2_HIGH:
      return "ec2_high";
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

vector<double> createHistograms(const vector<double>& uniqueInfectionTimes)
{
  vector<double> histogramInfo;

  double minTime = uniqueInfectionTimes[0];
  double maxTime = uniqueInfectionTimes[uniqueInfectionTimes.size() - 1];
  if (maxTime == INFECTED_FALSE) {
    if (uniqueInfectionTimes.size() > 2)
      maxTime = uniqueInfectionTimes[uniqueInfectionTimes.size() - 2];
    else
      maxTime = minTime + 0.1;
  }

  int totalBins = sqrt(uniqueInfectionTimes.size() - 1);
  histogramInfo.push_back(totalBins);

  for (int bin = 0; bin <= totalBins; ++bin) {
    double binMinTime = ((double) bin / totalBins) * (maxTime - minTime);
    histogramInfo.push_back(binMinTime);
  }

  // We always consider the non-infected case too, but in the getHistogramBin()

  return histogramInfo;
}

vector<double> createDiscreteTimeValuesFromHistogramBounds(
    const std::vector<double>& histogramInfo)
{
  // Always use the left-most value of the bin.
  vector<double> times;
  int bins = histogramInfo[0];

  // Maybe we consider just boolean values.
  if (bins == INFECTED_TRUE)
    times.push_back(INFECTED_TRUE);

  for (int bin = 1; bin <= bins; ++bin)
    times.push_back((double) (histogramInfo[bin] + histogramInfo[bin+1]) / 2);
  times.push_back(INFECTED_FALSE);

  return times;
}

int getHistogramBin(double time, const vector<double>& histogramInfo)
{
  int totalBins = histogramInfo[0];
  if (time == INFECTED_FALSE)
    return totalBins + 1;

  for (int bin = 1; bin <= totalBins; ++bin)
    if (time >= histogramInfo[bin] && time < histogramInfo[bin+1])
      return bin;

  // If it's in the last bin, the upper limit.
  if (histogramInfo[totalBins+1] == time)
    return totalBins;

  // Just warn...
  cout << "This shouldn't happen " << time << endl;
  for (size_t i = 0; i < histogramInfo.size(); ++i)
    cout << histogramInfo[i] << " ";
  cout << endl;

  return totalBins + 1;
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
  const bool paramIgnoreTime = Env.GetIfArgPrefixBool(
      "-time=", false,
      "Ignore Time: Boolean Vals (true), Discrete Values(false)");

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
  config.ignoreTime = paramIgnoreTime;
  config.objType = static_cast<AlgorithmType>(paramObjType.Val);
  config.infType = static_cast<InfectionType>(paramInfectionType.Val);
  if (config.infType == GAUSSIAN) {
    config.cluster.beta = 1.00 / config.cluster.miu;
  }

  config.groundTruthFile = paramGroundTruthFile;
  config.groundTruths = paramGroundTruths;
  if (!paramGroundTruthFile.Empty())
    config.groundTruths = 1;

  return config;
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
  cout << "Fraction of nodes in largest WCC: " <<
      TSnap::GetMxWccSz(m_network) << endl;
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
    objectives.push_back(VOI);
    objectives.push_back(GBS);
    objectives.push_back(RANDOM);
    objectives.push_back(EC2);
    objectives.push_back(EC2_HIGH);
    objectives.push_back(EPFL_ML);
    objectives.push_back(EPFL_EC2);
    objectives.push_back(VOI);
  } else {
    // Run only with the specified objective.
    objectives.push_back(static_cast<AlgorithmType>(m_config.objType));
  }

  for (size_t i = 0; i < objectives.size(); ++i) {
    // Reset configuration.
    m_config = initialConfig;
    // Set objective properly.
    m_config.objType = objectives[i];
    // Set logfile properly.
    m_config.logfile = TStr::Fmt("%s%d_%s.log", m_config.logfile.CStr(),
        m_config.nodes, algorithmTypeToString(objectives[i]));

    runWithCurrentConfig();
  }
}
