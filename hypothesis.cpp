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
                                 int infectionStep)
  : m_sourceId(sourceId)
  , m_infectionStep(infectionStep)
{
  m_infectionHash.swap(infectionTime);
}

bool GraphHypothesis::isConsistentWithTest(
    const GraphTest& test, const vector<pair<double,int>>& prevTests) const {
  // The test wasn't validated with the realization (so it's value is set to
  // true/false, depending on the computation that needs to be done).
  if (test.getInfectionTime() == INFECTED_UNDEFINED)
    return test.getOutcome() == this->getTestOutcome(test);

  if (test.getInfectionTime() == INFECTED_FALSE)
    return !this->getTestOutcome(test);

  double infectionTime = test.getInfectionTime();
  return this->getTestOutcome(test, prevTests) ||
      m_infectionStep < infectionTime;
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
  if (getInfectionTime(test.getNodeId()) == INFECTED_FALSE)
      return false;

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
  int runTimes = isTrueHypothesis ? trueCascadeSize : config.simulations;
  int infectionStep = 1;

  // Add the source node (fixed for this cluster).
  unordered_map<int, double> infectionTime;
  infectionTime[sourceId] = 0;

  // Make sure worker nodes have different seeds.
  TInt::Rnd.PutSeed(0);
  for (int run = 0; run < runTimes; run++) {
    infectionStep++;
    for (const auto& p : infectionTime) {
      const TUNGraph::TNodeI& crtIt = network->GetNI(p.first);
      for (int neighbour = 0; neighbour < crtIt.GetOutDeg(); ++neighbour) {
        if (TInt::Rnd.GetUniDev() > config.beta) // Flip a coin!
          continue;

        unsigned int neighbourId = crtIt.GetOutNId(neighbour);
        if (infectionTime.find(neighbourId) != infectionTime.end())
          continue;

        infectionTime[neighbourId] = infectionStep;
        if (infectionTime.size() == trueCascadeSize)
          return GraphHypothesis(sourceId, infectionTime, infectionStep);
      }
    }
  }
  return GraphHypothesis(sourceId, infectionTime, infectionStep);
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
    }
  }

  return GraphHypothesis(sourceId, infectionTime, infectionTime.size());
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

GraphHypothesis GraphHypothesis::readHypothesisFromFile(const char* filename)
{
  int node;
  int maxInfectionStep = -1;
  int infectionTime;
  int srcNode = -1;
  unordered_map<int, double> infectionTimeMap;

  ifstream inputfile(filename);
  for(string line; getline(inputfile, line); ) {
    istringstream iss(line);
    iss >> node >> infectionTime;
    infectionTimeMap[node] = infectionTime;
    maxInfectionStep = max(infectionTime, maxInfectionStep);
    if (srcNode == -1)
      srcNode = node;
  }
  inputfile.close();

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
