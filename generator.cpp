/**
 * Dumb simple example to get used to SNAP generators.
 */
#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>

#include "snap/snap-core/Snap.h"

#undef max
#undef min

#include <iostream>
#include <boost/network/uri.hpp>
#include <boost/network/uri/uri_io.hpp>
using namespace boost::network;

using namespace std;

int main(int argc, char* argv[]) {
  Env = TEnv(argc, argv, TNotify::StdNotify);
  Env.PrepArgs(TStr::Fmt("NetLoc. Build: %s, %s. Time: %s",
        __TIME__, __DATE__, TExeTm::GetCurTm()), 1, true);

  const TStr networkFile = Env.GetIfArgPrefixStr(
      "-network=", "meme.log", "File from where to read.");
  const TStr memeTrackerFile = Env.GetIfArgPrefixStr(
      "-meme=", "tracker.log", "File from where to read cascades.");
  const TStr cascadeoutFile = Env.GetIfArgPrefixStr(
      "-cascadeOut=", "cascade.txt", "Write actual cascade there");
  const TStr netoutFile = Env.GetIfArgPrefixStr(
      "-netOut=", "memeNetwork.dat", "Write actual network there");

  ifstream inputfile(networkFile.CStr());
  string line;

  // Build the network from the most popular websites.
  getline(inputfile, line); // read header line from file
  map<string, unsigned int> urlToNodeIdHash;

  PUNGraph graph = PUNGraph::New();
  unsigned int nodeId = 0;
  while (getline(inputfile, line)) {
    istringstream iss(line);

    string idx, src, dst;
    iss >> idx >> src >> dst;
    cout << idx << " " << src << " " << dst << endl;

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
  { TFOut FOut(netoutFile); graph->Save(FOut); }

  // Build one cascade
  ifstream memetracker(memeTrackerFile.CStr());

  // jump straight to an interestesting entry.
  for (int i = 0; i < 75926; ++i)
    getline(memetracker, line);

  // Read first meme entry.
  getline(memetracker, line);
  int entries, dummy;
  string phrase;

  istringstream iss(line);
  iss >> dummy >> entries;

  cout << "Building cascade for ";
  while (!iss.eof()) {
    iss >> phrase;
    cout << phrase << " ";
  }
  cout << endl;

  // Dump cascade to some file.
  ofstream dumpStream;
  dumpStream.open(cascadeoutFile.CStr());

  map<string, unsigned int> siteInfectionTimeHash;
  unsigned int infectionTime = 0;
  for (int i = 0; i < entries; ++i) {
    getline(memetracker, line);
    istringstream iss(line);

    string date1, date2;
    iss >> date1 >> date2;

    string dummy;
    iss >> dummy >> dummy;

    string url;
    iss >> url;


    uri::uri instance(url);
    assert(instance.is_valid());

    // If node not in network or already infected, skip.
    if (urlToNodeIdHash.find(instance.host()) == urlToNodeIdHash.end())
      continue;
    if (siteInfectionTimeHash.find(instance.host()) != siteInfectionTimeHash.end())
      continue;

    cout << date1 << " " << date2 << " " << url;
    cout << "host: " << instance.host() << endl;

    siteInfectionTimeHash[instance.host()] = infectionTime++;
    dumpStream << urlToNodeIdHash[instance.host()] << " " <<
      siteInfectionTimeHash[instance.host()] << endl;
  }

  return 0;
}
