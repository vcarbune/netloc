/**
 * 
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
  Env.PrepArgs(TStr::Fmt("InfMax. Build: %s, %s. Time: %s", __TIME__, __DATE__, TExeTm::GetCurTm()));

  const TInt paramSimulation = Env.GetIfArgPrefixInt(
      "-sim=", NodeVar, "Simulation Type"
          "(NodeVar - 0, EdgeVar - 1, BetaVar - 2, HypothesisVar - 3, CascBoundVar - 4)");
  const TInt paramNodes = Env.GetIfArgPrefixInt(
      "-n=", 0, "Network Size (number of nodes)");
  const TInt paramClusterSize = Env.GetIfArgPrefixInt(
      "-c=", 10, "Cluster Size");
  const double paramCascadeSize = Env.GetIfArgPrefixFlt(
      "-s=", 0.4, "Cascade Size (percentage of network size)");
  const double paramBeta = Env.GetIfArgPrefixFlt(
      "-b=", 0.1, "Beta (activation probability on edges)");
  const TInt paramSteps = Env.GetIfArgPrefixInt(
      "-steps=", 5, "Number of simulation steps.");
  const TInt paramStartStep = Env.GetIfArgPrefixInt(
      "-start=", 0, "Simulation step (could be used for resuming)");
  const TStr dumpFile = Env.GetIfArgPrefixStr(
      "-dump=", "dump.log", "File where to dump the output");

  SimConfig config(static_cast<SimulationType>(paramSimulation.Val));

  if (paramNodes.Val)
    config.nodes = paramNodes.Val;

  config.clusterSize = paramClusterSize.Val;
  config.steps = paramSteps.Val;
  config.cascadeBound = paramCascadeSize;
  config.beta = paramBeta;
  config.logfile = dumpFile;

  config += paramStartStep.Val;

  return config;
}
