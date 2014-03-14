/**
 * Synthetic Data Generator for evaluating performance of source localization.
 * Author: Victor Carbune (vcarbune@ethz.ch)
 */
#include <ctime>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>

#include "../snap/snap-core/Snap.h"
#undef max
#undef min

#include "../hypothesis.h"

using namespace std;

int main(int argc, char* argv[]) {
  Env = TEnv(argc, argv, TNotify::StdNotify);
  Env.PrepArgs(TStr::Fmt("NetLoc Synthetic Data Generator."
       "Build: %s, %s. Time: %s",  __TIME__, __DATE__,
       TExeTm::GetCurTm()), 1, true);

  // Generate Network
  const TInt networkSize = Env.GetIfArgPrefixInt(
      "-networkSize=", 500, "Network Size");
  const TInt networkType = Env.GetIfArgPrefixInt(
      "-networkType=", 0,
      "Network Type (0 - forest, 1 - barabasi, 2 - erdos)");

  PUNGraph network;
  const char* networkTypeStr = "unknown";
  switch (networkType) {
    case 0:
      network = TSnap::ConvertGraph<PUNGraph,
              PNGraph>(TSnap::GenForestFire(networkSize, 0.35, 0.32));
      networkTypeStr = "forest";
      break;
    case 1:
      network = TSnap::GenPrefAttach(networkSize, 3);
      networkTypeStr = "barabasi";
      break;
    case 2:
      network = TSnap::GenRndGnm<PUNGraph>(
          networkSize, networkSize * log(networkSize.Val));
      networkTypeStr = "erdos";
      break;
  }

  cout << "Percentage of nodes in largest connected component: " <<
    TSnap::GetMxWccSz(network) << endl;
  cout << "Saving network to file..." << endl;

  const TStr snapNetworkOutputFile = Env.GetIfArgPrefixStr(
      "-snapNetworkOutputFile=",
      TStr::Fmt("data/synthetic/%s%d.dat", networkTypeStr, networkSize),
      "File where to write the output network as a snap binary");
  { TFOut FOut(snapNetworkOutputFile); network->Save(FOut); }

  // Generate Ground Truth on top of the network.
  const double groundTruthSize = Env.GetIfArgPrefixFlt(
      "-groundTruthSize=", 0.4, "Ground Truth Size");
  const double  edgeActivationProbability = Env.GetIfArgPrefixFlt(
      "-beta=", 0.06, "Activation probability on edges");
  GraphHypothesis realization = GraphHypothesis::generateHypothesis(
      network, rand() % network->GetNodes(), groundTruthSize,
      edgeActivationProbability);

  const TStr groundTruthOutputFile = Env.GetIfArgPrefixStr(
      "-groundTruthOutputFile=",
      TStr::Fmt("data/synthetic/%s%d_gt.txt", networkTypeStr, networkSize),
      "File where to write the ground truth");
  realization.writeHypothesisToFile(groundTruthOutputFile.CStr());

  return 0;
}
