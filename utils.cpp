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
