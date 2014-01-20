/**
 * Data structure implementation to represent a hypothesis cluster.
 * Author: Victor Carbune (vcarbune@ethz.ch)
 */

#include <cassert>
#include <iostream>

#include "hypothesis.h"

#define INITIAL_RUNS 5
#define EPS 0.05

GraphHypothesis::GraphHypothesis(unordered_map<int, int>& infectionTime)
{
  m_infectionHash.swap(infectionTime);
}

bool GraphHypothesis::isConsistentWithTest(const GraphTest& test) const {
  return test.getOutcome() == this->getTestOutcome(test);
}

int GraphHypothesis::getInfectionTime(int nodeId) const
{
  if (m_infectionHash.find(nodeId) != m_infectionHash.end())
    return m_infectionHash.at(nodeId);

  return -1;
}

bool GraphHypothesis::getTestOutcome(const GraphTest& test) const
{
  return m_infectionHash.find(test.getNodeId()) != m_infectionHash.end();
}

GraphHypothesisCluster::GraphHypothesisCluster(PUNGraph network,
                                               int sourceId,
                                               int clusterSize,
                                               double beta,
                                               double size,
                                               double weight)
  : m_network(network)
  , m_sourceId(sourceId)
  , m_weight(weight)
{
  for (int i = 0; i < network->GetNodes(); ++i)
    m_nodeCount.push_back(0);
  generateHypothesisCluster(size, beta, clusterSize);
}

/**
 * Internal generator for all the hypothesis in the cluster.
 */
void GraphHypothesisCluster::generateHypothesisCluster(
    double size, double beta, int clusterSize)
{
  // TODO(vcarbune): threadify
  for (int h = 0; h < clusterSize; h++) {
    m_hypothesis.push_back(generateHypothesis(size, beta, &m_nodeCount));
    m_hypothesis[h].weight = m_weight / clusterSize;
  }
}

/**
 * Generates one cascade, on top of the underlying network structure.
 * Extended from snap/examples/cascades.
 *
 * TODO(vcarbune): Do we want to sample the cascade (e.g. removing nodes?)
 */
GraphHypothesis GraphHypothesisCluster::generateHypothesis(
    double size, double beta, vector<int> *nodeCount) const
{
  bool isTrueHypothesis = nodeCount == NULL;

  unsigned int cascadeSize = size * m_network->GetNodes();
  int runTimes = isTrueHypothesis ? cascadeSize : INITIAL_RUNS;

  // Add the source node (fixed for this cluster).
  unordered_map<int, int> infectionTime;
  infectionTime[m_sourceId] = 0;

  for (int run = 0; run < runTimes; run++) {
    for (const auto& p : infectionTime) {
      const TUNGraph::TNodeI& crtIt = m_network->GetNI(p.first);
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
          return GraphHypothesis(infectionTime);
      }
    }
  }

  return GraphHypothesis(infectionTime);
}

void GraphHypothesisCluster::updateMassWithTest(const GraphTest& test)
{
  m_weight = 0;
  for (GraphHypothesis& h : m_hypothesis) {
    h.weight *= (h.isConsistentWithTest(test) ? (1-EPS) : EPS);
    m_weight += h.weight;
  }
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
