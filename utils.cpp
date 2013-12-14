/**
 * Various utility functions.
 */

#include "utils.h"

using namespace std;

void plotGraphs(const PUNGraph& network, const vector<GraphTest>& tests,
                const GraphHypothesis& realization)
{
  TIntStrH nodeTestLabels;
  TIntStrH nodeInfectionLabels;

  for (size_t i = 0; i < tests.size(); ++i) {
    int id = tests[i].getNodeId();
    nodeTestLabels.AddDat(id, TStr::Fmt("%d", i));
  }

  for (int i = 0; i < network->GetNodes(); ++i) {
    int infTime = realization.getInfectionTime(i);
    if (infTime != -1)
      nodeInfectionLabels.AddDat(i, TStr::Fmt("%d", infTime));
    else
      nodeInfectionLabels.AddDat(i, TStr::Fmt(""));

    if (!nodeTestLabels.IsKey(i))
      nodeTestLabels.AddDat(i, TStr::Fmt(""));
  }

  // Visualize the order in which the nodes have been tested.
  TStr FNameDemo = TStr::Fmt("graph_with_test_order.%s", "png");
  remove(FNameDemo.CStr());
  TSnap::DrawGViz(network, TGVizLayout(2), FNameDemo, "Tests", nodeTestLabels);
  printf("Drawing graph '%s'\n", FNameDemo.CStr());

  // Visualization true infection order.
  FNameDemo = TStr::Fmt("graph_with_infection_order.%s", "png");
  remove(FNameDemo.CStr());
  TSnap::DrawGViz(network, TGVizLayout(2), FNameDemo, "Infection", nodeInfectionLabels);
  printf("Drawing graph '%s'\n", FNameDemo.CStr());
}

SimConfig getSimConfigFromEnv(int argc, char *argv[])
{
  Env = TEnv(argc, argv, TNotify::StdNotify);
  Env.PrepArgs(TStr::Fmt("NetLoc. Build: %s, %s. Time: %s", __TIME__, __DATE__, TExeTm::GetCurTm()));

  const TInt paramSimulation = Env.GetIfArgPrefixInt(
      "-sim=", NodeVar, "Simulation Type"
          "(NodeVar - 0, BetaVar - 1, HypothesisVar - 2, CascBoundVar - 3)");
  const TInt paramNodes = Env.GetIfArgPrefixInt(
      "-n=", 0, "Network size");
  const TInt paramClusterSize = Env.GetIfArgPrefixInt(
      "-c=", 100, "Cluster size");
  const double paramCascadeSize = Env.GetIfArgPrefixFlt(
      "-s=", 0.3, "True cascade size");
  const double paramBeta = Env.GetIfArgPrefixFlt(
      "-b=", 0.05, "Activation probability on edges)");
  const TInt paramSteps = Env.GetIfArgPrefixInt(
      "-steps=", 1, "Number of simulation steps.");
  const TInt paramStartStep = Env.GetIfArgPrefixInt(
      "-start=", 0, "Simulation step");
  const TInt paramKeepTopN = Env.GetIfArgPrefixInt(
      "-topN=", 1, "Keep topN solutions");
  const TStr dumpFile = Env.GetIfArgPrefixStr(
      "-dump=", "dump.log", "File where to dump the output");

  SimConfig config(static_cast<SimulationType>(paramSimulation.Val));

  if (paramNodes.Val) {
    config.nodes = paramNodes.Val;
    config.edges = config.nodes * log(config.nodes);
  }

  config.clusterSize = paramClusterSize.Val;
  config.steps = paramSteps.Val;
  config.cascadeBound = paramCascadeSize;
  config.beta = paramBeta;
  config.logfile = dumpFile;
  config.topN = paramKeepTopN;

  config += paramStartStep.Val;

  return config;
}
