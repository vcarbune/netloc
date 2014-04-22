/**
 * Codered Converted.
 * Author: Victor Carbune (vcarbune@ethz.ch)
 */
#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>

#include "../snap/snap-core/Snap.h"
#undef max
#undef min

using namespace std;

int main(int argc, char* argv[]) {
  Env = TEnv(argc, argv, TNotify::StdNotify);
  Env.PrepArgs(TStr::Fmt("Codered to SNAP Converter. Build: %s, %s. Time: %s",
        __TIME__, __DATE__, TExeTm::GetCurTm()), 1, true);

  const TStr aslinksInputFile = Env.GetIfArgPrefixStr(
      "-aslinksInputFile=",
      "data/codered/20010713.link.v4",
      "File containing edge lists between AS numbers");
  const TStr snapNetworkOutputFile = Env.GetIfArgPrefixStr(
      "-snapNetworkOutputFile=", "data/codered/coderedNetwork.dat",
      "File where to write the output network as a snap binary");

  // Build the network from the most popular websites.
  PUNGraph graph = TSnap::LoadEdgeList<PUNGraph>(aslinksInputFile, 0, 1);
  // Save network.
  { TFOut FOut(snapNetworkOutputFile); graph->Save(FOut); }
  return 0;
}
