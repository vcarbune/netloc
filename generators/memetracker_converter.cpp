/**
 * Memetracker Converter. Original data is taken from the MemeTracker project.
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

#include <boost/network/uri.hpp>
#include <boost/network/uri/uri_io.hpp>

using namespace boost::network;
using namespace std;

int main(int argc, char* argv[]) {
  Env = TEnv(argc, argv, TNotify::StdNotify);
  Env.PrepArgs(TStr::Fmt("Memetracker Converter. Build: %s, %s. Time: %s",
        __TIME__, __DATE__, TExeTm::GetCurTm()), 1, true);
  const TStr netinfInputFile = Env.GetIfArgPrefixStr(
      "-netinfInputFile=",
      "data/memetracker/InfoNet5000Q1000NEXP.txt",
      "File containing the inferred network, used as underlying structure");
  const TStr clusteredCascadesInputFile = Env.GetIfArgPrefixStr(
      "-clusteredCascadesInputFile=",
      "data/memetracker/clust-qt08080902w3mfq5.txt",
      "File containing clustered instances of memes and their path on the web");
  const TInt cascadeInputLine = Env.GetIfArgPrefixInt(
      "-cascadeInputLine=", 75926,
      "Specific cascade to read from clusteredCascadesInputFile");

  const TStr snapNetworkOutputFile = Env.GetIfArgPrefixStr(
      "-snapNetworkOutputFile=", "data/memetracker/memeNetwork.dat",
      "File where to write the output network as a snap binary");
  const TStr groundTruthOutputFile = Env.GetIfArgPrefixStr(
      "-groundTruthOutputFile=", "data/memetracker/memeGroundTruth.txt",
      "File where to write the ground truth");

  // Build the network from the most popular websites.
  ifstream inputfile(netinfInputFile.CStr());
  string line;
  getline(inputfile, line); // read header line from file
  map<string, unsigned int> urlToNodeIdHash;
  PUNGraph graph = PUNGraph::New();
  unsigned int nodeId = 0;
  while (getline(inputfile, line)) {
    istringstream iss(line);
    string idx, src, dst;
    iss >> idx >> src >> dst;
    if (urlToNodeIdHash.find(src) == urlToNodeIdHash.end()) {
      urlToNodeIdHash[src] = nodeId;
      graph->AddNode(nodeId);
      nodeId++;
    }
    if (urlToNodeIdHash.find(dst) == urlToNodeIdHash.end()) {
      urlToNodeIdHash[dst] = nodeId;
      graph->AddNode(nodeId);
      nodeId++;
    }
    graph->AddEdge(urlToNodeIdHash[src], urlToNodeIdHash[dst]);
  }
  // Save network.
  { TFOut FOut(snapNetworkOutputFile); graph->Save(FOut); }

  // Read one memetracker entry.
  ifstream memetracker(clusteredCascadesInputFile.CStr());
  for (int i = 0; i < cascadeInputLine; ++i)
    getline(memetracker, line);
  getline(memetracker, line);
  int entries, dummyInt;
  istringstream iss(line);
  iss >> dummyInt >> entries;

  cout << "Building cascade for ";
  while (!iss.eof()) {
    string phrase;
    iss >> phrase;
    cout << phrase << " ";
  }
  cout << endl;

  // Dump cascade to some file.
  ofstream dumpStream;
  dumpStream.open(groundTruthOutputFile.CStr());

  string dummy, url;
  map<string, unsigned int> infectionTimeHash;
  unsigned int infectionTime = 0;
  for (int i = 0; i < entries; ++i) {
    // Read through each "infected" URL.
    getline(memetracker, line);
    istringstream iss(line);
    // These fields of the cascade entry ar not important.
    iss >> dummy >> dummy >> dummy >> dummy;
    iss >> url;
    // Parse the URL and identify the host website.
    uri::uri instance(url);
    assert(instance.is_valid());
    // If node not in network or already infected, skip.
    if (urlToNodeIdHash.find(instance.host()) == urlToNodeIdHash.end() ||
        infectionTimeHash.find(instance.host()) != infectionTimeHash.end())
      continue;
    infectionTimeHash[instance.host()] = infectionTime++;
    // Dump as pair of <nodeId, infectionTime>.
    dumpStream << urlToNodeIdHash[instance.host()] << " " <<
      infectionTimeHash[instance.host()] << endl;
  }
  return 0;
}
