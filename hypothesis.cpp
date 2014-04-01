/**
 * Data structure implementation to represent a hypothesis cluster.
 * Author: Victor Carbune (vcarbune@ethz.ch)
 */

#include "hypothesis.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <numeric>
#include <random>
#include <sstream>
#include <queue>

#include <iostream>

using namespace std;

GraphHypothesis::GraphHypothesis(unsigned int sourceId,
                                 unordered_map<int, double>& infectionTime,
                                 double maxInfectionTime)
  : m_sourceId(sourceId)
  , m_maxInfectionTime(maxInfectionTime)
{
  m_infectionHash.swap(infectionTime);
}

bool GraphHypothesis::isConsistentWithTest(
    const GraphTest& test, const vector<pair<double,int>>& prevTests) const
{
  double infectionTime = test.getInfectionTime();

  // The test wasn't validated with the realization (so it's value is set to
  // true/false, depending on the computation that needs to be done).
  if (fabs(infectionTime - INFECTED_UNDEFINED) < 1E-1)
    return test.getOutcome() == this->getTestOutcome(test);

  if (fabs(infectionTime - INFECTED_FALSE) < 1E-1)
    return !this->getTestOutcome(test);

  // There are two ways of treating this:
  // - by default down-weight, because the cluster is definitely not the source.
  // - by default do nothing, because the point couldn't be reached in this  instance.
  /*
  if (m_maxInfectionTime < infectionTime)
    return true;
  */
  return this->getTestOutcome(test, prevTests);
}

double GraphHypothesis::getInfectionTime(int nodeId) const
{
  if (m_infectionHash.find(nodeId) != m_infectionHash.end())
    return m_infectionHash.at(nodeId);

  return INFECTED_FALSE;
}

bool GraphHypothesis::getTestOutcome(const GraphTest& test) const
{
  return m_infectionHash.find(test.getNodeId()) != m_infectionHash.end();
}

bool GraphHypothesis::getTestOutcome(
    const GraphTest& test, const vector<pair<double, int>>& prevTests) const
{
  // Make sure the test keeps order with respect to other tests.
  for (size_t i = 0; i < prevTests.size(); ++i) {
    if (getInfectionTime(prevTests[i].second) == INFECTED_FALSE ||
        prevTests[i].first == INFECTED_FALSE)
      continue;

    double realizationDiffTime = test.getInfectionTime() - prevTests[i].first;
    double cascadeDiffTime = getInfectionTime(test.getNodeId()) -
        getInfectionTime(prevTests[i].second);

    // If they have different signs, then the order is not preserved.
    if (realizationDiffTime * cascadeDiffTime < 0)
      return false;
  }

  return true;
}

/**
 * Generates one cascade, on top of the underlying network structure.
 * Extended from snap/examples/cascades.
 */
GraphHypothesis GraphHypothesis::generateHypothesis(
    PUNGraph network, int sourceId, const HypothesisClusterConfig config,
    bool isTrueHypothesis)
{
  unsigned int trueCascadeSize = config.bound * network->GetNodes();
  unsigned int artifCascadeSize = config.cbound * network->GetNodes();
  unsigned int maxCascadeSize = isTrueHypothesis ? trueCascadeSize : artifCascadeSize;
  double maxInfectionTime = 1;

  // Add the source node (fixed for this cluster).
  unordered_map<int, double> infectionTime;
  infectionTime[sourceId] = 0;

  // Make sure worker nodes have different seeds.
  // TInt::Rnd.PutSeed(0);
  for (unsigned int run = 0; run < maxCascadeSize; run++) {
    maxInfectionTime++;
    for (const auto& p : infectionTime) {
      const TUNGraph::TNodeI& crtIt = network->GetNI(p.first);
      for (int neighbour = 0; neighbour < crtIt.GetOutDeg(); ++neighbour) {
        if (TInt::Rnd.GetUniDev() > config.beta) // Flip a coin!
          continue;

        unsigned int neighbourId = crtIt.GetOutNId(neighbour);
        if (infectionTime.find(neighbourId) != infectionTime.end())
          continue;

        infectionTime[neighbourId] = maxInfectionTime;
        if (infectionTime.size() == maxCascadeSize)
          return GraphHypothesis(sourceId, infectionTime, maxInfectionTime);
      }
    }
  }
  return GraphHypothesis(sourceId, infectionTime, maxInfectionTime);
}

/**
 * Generate one cascade using the EPFL Gaussian diffusion model.
 * See paper "Locating the Source of Diffusion in Large-Scale Networks".
 */
GraphHypothesis GraphHypothesis::generateHypothesisUsingGaussianModel(
    PUNGraph network, int sourceId, const HypothesisClusterConfig cluster,
    bool)
{
  // Assign propagation delays to all edges using a normal distribution.
  std::default_random_engine generator;

  // miu/sigma = 4 (synthetic data)
  std::normal_distribution<double> d(cluster.miu, cluster.sigma);

  unordered_map<int, double> infectionTime;
  infectionTime[sourceId] = 0;
  double maxInfectionTime = 0;

  // bfs to mark infection times.
  queue<int> q;
  q.push(sourceId);
  while (!q.empty()) {
    int crtNode = q.front();
    q.pop();

    const TUNGraph::TNodeI& crtIt = network->GetNI(crtNode);
    for (int neighbour = 0; neighbour < crtIt.GetOutDeg(); ++neighbour) {
      unsigned int neighbourId = crtIt.GetOutNId(neighbour);
      if (infectionTime.find(neighbourId) != infectionTime.end())
        continue;

      q.push(neighbourId);
      infectionTime[neighbourId] = infectionTime[crtNode] + d(generator);
      maxInfectionTime = max(infectionTime[neighbourId], maxInfectionTime);
    }
  }

  return GraphHypothesis(sourceId, infectionTime, maxInfectionTime);
}

/**
 * Generates a hypothesis on top of a weighted directed graph,
 * starting from a given source id.
 */
/*
GraphHypothesis GraphHypothesis::generateHypothesisUsingWeightedGraph(
  const TNodeEDatNet<int, double>& network, int sourceId,
  const HypothesisClusterConfig cluster)
{
}
*/

GraphHypothesis GraphHypothesis::readHypothesisFromFile(
    PUNGraph network,
    const char* filename)
{
  int node;
  int maxInfectionStep = -1;
  int infectionTime;
  int srcNode = -1;
  int skippedNodes = 0;
  unordered_map<int, double> infectionTimeMap;

  ifstream inputfile(filename);
  for(string line; getline(inputfile, line); ) {
    istringstream iss(line);
    iss >> node >> infectionTime;
    if (!network->IsNode(node)) {
      skippedNodes++;
      continue;
    }
    infectionTimeMap[node] = infectionTime;
    maxInfectionStep = max(infectionTime, maxInfectionStep);
    if (srcNode == -1 && infectionTime == 0.00) {
      srcNode = node;
      cout << "SOURCE: " << srcNode << endl;
    }
  }
  inputfile.close();
  cout << "Warning: " << skippedNodes << " nodes were not in the network." << endl;
  return GraphHypothesis(srcNode, infectionTimeMap, maxInfectionStep);
}

void GraphHypothesis::writeHypothesisToFile(const char* filename)
{
  ofstream outputStream;
  outputStream.open(filename);
  for (auto it = m_infectionHash.begin(); it != m_infectionHash.end(); ++it)
    outputStream << it->first << " " << it->second << endl;
  outputStream.close();
}

GraphHypothesisCluster::GraphHypothesisCluster(PUNGraph network,
                                               int sourceId,
                                               double weight)
  : m_network(network)
  , m_sourceId(sourceId)
  , m_weight(weight)
{
}

/**
 * Internal generator for all the hypothesis in the cluster.
 */
GraphHypothesisCluster GraphHypothesisCluster::generateHypothesisCluster(
    PUNGraph network, int source, double weight, const SimConfig& config)
{
  GraphHypothesisCluster cluster(network, source, weight);
  for (int h = 0; h < config.cluster.size; h++) {
    cluster.m_hypothesis.push_back(GraphHypothesis::generateHypothesis(
          network, source, config.cluster, false));
    cluster.m_hypothesis[h].weight = weight / config.cluster.size;
  }
  return cluster;
}

void GraphHypothesisCluster::updateMassWithTest(const double eps,
    const GraphTest& test, const vector<pair<double,int>>& prevTests)
{
  double weight = 0.0;
  for (GraphHypothesis& h : m_hypothesis) {
    h.weight *= (h.isConsistentWithTest(test, prevTests) ? (1-eps) : eps);
    weight += h.weight;
  }
  m_weight = weight;
}

pair<double, double> GraphHypothesisCluster::computeMassWithTest(
    double& positiveTestPrior,
    const double eps, const GraphTest& test,
    const vector<pair<double, int>>& prevTests) const
{
  double positiveMass = 0.0;
  double negativeMass = 0.0;
  for (const GraphHypothesis& h : m_hypothesis) {
    bool outcome = h.isConsistentWithTest(test, prevTests);
    positiveMass += h.weight * (outcome ? (1-eps) : eps);
    negativeMass += h.weight * (outcome ? eps : (1-eps));
    if (outcome)
      positiveTestPrior += h.weight;
  }
  return pair<double, double>(positiveMass, negativeMass);
}

double GraphHypothesisCluster::getAverageHypothesisSize() const
{
  double size = 0.0;
  for (const GraphHypothesis& h : m_hypothesis)
    size += h.getSize();
  return size / m_hypothesis.size();
}

void GraphHypothesisCluster::resetWeight(double weight)
{
  double hWeight = (double) weight / m_hypothesis.size();
  for (GraphHypothesis& h : m_hypothesis)
    h.weight = hWeight;
  m_weight = weight;
}
