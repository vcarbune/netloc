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
                                 unordered_map<int, double>& infectionTime)
  : m_sourceId(sourceId)
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
      m_infectionHash.size() < infectionTime;
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

  // Add the source node (fixed for this cluster).
  unordered_map<int, double> infectionTime;
  infectionTime[sourceId] = 0;

  // Make sure worker nodes have different seeds.
  TInt::Rnd.PutSeed(0);
  for (int run = 0; run < runTimes; run++) {
    for (const auto& p : infectionTime) {
      const TUNGraph::TNodeI& crtIt = network->GetNI(p.first);
      for (int neighbour = 0; neighbour < crtIt.GetOutDeg(); ++neighbour) {
        if (TInt::Rnd.GetUniDev() > config.beta) // Flip a coin!
          continue;

        unsigned int neighbourId = crtIt.GetOutNId(neighbour);
        if (infectionTime.find(neighbourId) != infectionTime.end())
          continue;

        infectionTime[neighbourId] = infectionTime.size();
        if (infectionTime.size() == trueCascadeSize)
          return GraphHypothesis(sourceId, infectionTime);
      }
    }
  }
  return GraphHypothesis(sourceId, infectionTime);
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
  std::random_device rd;
  std::mt19937 gen(rd());

  // miu/sigma = 4 (synthetic data)
  std::normal_distribution<> d(cluster.miu, cluster.sigma);

  unordered_map<int, double> infectionTime;
  infectionTime[sourceId] = 0;

  unordered_map<int, int> visited;
  queue<int> q;

  //int infectedNodes = 0;
  q.push(sourceId);

  // bfs to mark infection times.
  while (!q.empty()) {
    int crtNode = q.front();
    q.pop();

    const TUNGraph::TNodeI& crtIt = network->GetNI(crtNode);
    for (int neighbour = 0; neighbour < crtIt.GetOutDeg(); ++neighbour) {
      unsigned int neighbourId = crtIt.GetOutNId(neighbour);
      if (infectionTime.find(neighbourId) != infectionTime.end())
        continue;

      q.push(neighbourId);
      infectionTime[neighbourId] = infectionTime[crtNode] + d(gen);
    }
  }

  return GraphHypothesis(sourceId, infectionTime);
}

GraphHypothesis GraphHypothesis::readHypothesisFromFile(const char* filename)
{
  int node; 
  int infectionTime;
  int srcNode = -1;
  unordered_map<int, double> infectionTimeMap;

  ifstream inputfile(filename);
  for(string line; getline(inputfile, line); ) {
    istringstream iss(line);
    iss >> node >> infectionTime;
    infectionTimeMap[node] = infectionTime;
    if (srcNode == -1)
      srcNode = node;
  }
  inputfile.close();

  return GraphHypothesis(srcNode, infectionTimeMap);
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
    if (config.infType == BETA) {
      cluster.m_hypothesis.push_back(GraphHypothesis::generateHypothesis(
            network, source, config.cluster, false));
    } else if (config.infType == GAUSSIAN) {
      cluster.m_hypothesis.push_back(
          GraphHypothesis::generateHypothesisUsingGaussianModel(
            network, source, config.cluster, false));
    }
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

void GraphHypothesisCluster::resetWeight(double weight)
{
  m_weight = weight;
  double hWeight = (double) weight / m_hypothesis.size();
  for (size_t i = 0; i < m_hypothesis.size(); ++i)
    m_hypothesis[i].weight = hWeight;
}
