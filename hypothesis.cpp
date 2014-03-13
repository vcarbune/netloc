/**
 * Data structure implementation to represent a hypothesis cluster.
 * Author: Victor Carbune (vcarbune@ethz.ch)
 */

#include <cassert>
#include <cmath>
#include <fstream>
#include <numeric>
#include <random>
#include <sstream>

#include <queue>

#include "hypothesis.h"

#define INITIAL_RUNS 4

GraphHypothesis::GraphHypothesis(unsigned int sourceId,
                                 unordered_map<int, double>& infectionTime)
  : m_sourceId(sourceId)
{
  m_infectionHash.swap(infectionTime);
}

bool GraphHypothesis::isConsistentWithTest(const GraphTest& test) const {
  if (test.getInfectionTime() == INFECTED_UNDEFINED)
    return test.getOutcome() == this->getTestOutcome(test);

  if (test.getInfectionTime() == INFECTED_FALSE)
    return !this->getTestOutcome(test);

  double infectionTime = test.getInfectionTime();
  return this->getTestOutcome(test) || m_infectionHash.size() < infectionTime;
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

/**
 * Generates one cascade, on top of the underlying network structure.
 * Extended from snap/examples/cascades.
 */
GraphHypothesis GraphHypothesis::generateHypothesis(PUNGraph network,
    int sourceId, double size, double beta, vector<int> *nodeCount)
{
  bool isTrueHypothesis = nodeCount == NULL;

  unsigned int cascadeSize = size * network->GetNodes();
  int runTimes = isTrueHypothesis ? cascadeSize : INITIAL_RUNS;

  // Add the source node (fixed for this cluster).
  unordered_map<int, double> infectionTime;
  infectionTime[sourceId] = 0;

  for (int run = 0; run < runTimes; run++) {
    for (const auto& p : infectionTime) {
      const TUNGraph::TNodeI& crtIt = network->GetNI(p.first);
      for (int neighbour = 0; neighbour < crtIt.GetOutDeg(); ++neighbour) {
        if (TInt::Rnd.GetUniDev() > beta) // Flip a coin!
            continue;

        unsigned int neighbourId = crtIt.GetOutNId(neighbour);
        if (infectionTime.find(neighbourId) != infectionTime.end())
          continue;

        if (nodeCount)
          (*nodeCount)[neighbourId]++;

        infectionTime[neighbourId] = infectionTime.size();
        if (infectionTime.size() == cascadeSize)
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
    PUNGraph network, int sourceId, double size,
    double miu, double sigma, vector<int> *nodeCount)
{
  bool isTrueHypothesis = nodeCount == NULL;
  unsigned int cascadeSize = size * network->GetNodes();

  // Assign propagation delays to all edges using a normal distribution.
  std::random_device rd;
  std::mt19937 gen(rd());

  // miu/sigma = 4 (synthetic data)
  std::normal_distribution<> d(miu, sigma);

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

GraphHypothesisCluster::GraphHypothesisCluster(PUNGraph network,
                                               int sourceId,
                                               double weight)
  : m_network(network)
  , m_sourceId(sourceId)
  , m_weight(weight)
{
  m_nodeCount.resize(network->GetNodes(), 0);
}

/**
 * Internal generator for all the hypothesis in the cluster.
 */
GraphHypothesisCluster GraphHypothesisCluster::generateHypothesisCluster(
    PUNGraph network, int source, double weight,
    double beta, int size, int clusterSize)
{
  GraphHypothesisCluster cluster(network, source, weight);
  for (int h = 0; h < clusterSize; h++) {
    cluster.m_hypothesis.push_back(
        GraphHypothesis::generateHypothesis(network, source, size, beta,
            &cluster.m_nodeCount));
    cluster.m_hypothesis[h].weight = weight / clusterSize;
  }
  return cluster;
}

void GraphHypothesisCluster::updateMassWithTest(const GraphTest& test)
{
  double weight = 0.0;
  for (GraphHypothesis& h : m_hypothesis) {
    h.weight *= (h.isConsistentWithTest(test) ? (1-EPS) : EPS);
    weight += h.weight;
  }
  m_weight = weight;
}

pair<double, double> GraphHypothesisCluster::computeMassWithTest(
    const GraphTest& test) const
{
  double positiveMass = 0.0;
  double negativeMass = 0.0;
  for (const GraphHypothesis& h : m_hypothesis) {
    bool outcome = h.getTestOutcome(test);
    positiveMass += h.weight * (outcome ? (1-EPS) : EPS);
    negativeMass += h.weight * (outcome ? EPS : (1-EPS));
  }
  return pair<double, double>(positiveMass, negativeMass);
}

void GraphHypothesisCluster::resetWeight(double weight)
{
  m_weight = weight;
  double hWeight = weight / m_hypothesis.size();
  for (size_t i = 0; i < m_hypothesis.size(); ++i)
    m_hypothesis[i].weight = hWeight;
}
