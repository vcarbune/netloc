/**
 * Just a simple example to get used to SNAP.
 */
#include <iostream>
#include <ctime>

#include "snap/snap-core/Snap.h"

using namespace std;

void populateGraph(PNEGraph g) {
  srand48(time(NULL));;

  int nodes = 10 + (long) (drand48() * 10);
  int edges = 5 + (long) (drand48() * 5);

  for (int i = 0; i < nodes; ++i)
    g->AddNode(i);

  for (int i = 0; i < edges; ++i)
    g->AddEdge((long) (drand48() * nodes), (long) (drand48() * nodes));

}

void printGraph(PNEGraph g) {
  cout << "Nodes: " << g->GetNodes() << endl;
  for (TNEGraph::TNodeI it = g->BegNI(); it < g->EndNI(); it++)
    cout << "Node " << it.GetId() << ": " << it.GetDeg() << ", " <<
      it.GetInDeg() << endl;

  cout << endl;

  cout << "Edges: " << g->GetEdges() << endl;
  for (TNEGraph::TEdgeI it = g->BegEI(); it < g->EndEI(); it++)
    cout << "Edge " << it.GetId() << ": (" << it.GetSrcNId() << ", " <<
      it.GetDstNId() << ")" << endl;
}

void saveGraph(PNEGraph g, const char *filename) {
  TFOut fout(filename);
  g->Save(fout);
  fout.Flush();
}


int main(int argc, char* argv[]) {
  PNEGraph graph = PNEGraph::New();

  populateGraph(graph);
  printGraph(graph);
  saveGraph(graph, "graph.dat");

  return 0;
}
