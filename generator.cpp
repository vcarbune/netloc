/**
 * Dumb simple example to get used to SNAP generators.
 */
#include <iostream>
#include <ctime>

#include "snap/snap-core/Snap.h"

using namespace std;

int main(int argc, char* argv[]) {
  PUNGraph fullGraph = TSnap::GenFull<PUNGraph>(10);
  cout << fullGraph->GetNodes() << " " << fullGraph->GetEdges() << endl;

  PUNGraph starGraph = TSnap::GenStar<PUNGraph>(10);
  cout << starGraph->GetNodes() << " " << starGraph->GetEdges() << endl;

  return 0;
}
