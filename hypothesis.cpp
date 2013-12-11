/**
 * Data structure implementation to represent a hypothesis cluster.
 * Author: Victor Carbune (vcarbune@ethz.ch)
 */

#include <cassert>
#include <iostream>

#include "hypothesis.h"

#define INITIAL_RUNS 10

GraphHypothesis::GraphHypothesis(TIntH hash, double weight)
  : m_weight(weight)
  , m_infectionTimeHash(hash)
{
}

bool GraphHypothesis::isConsistentWithTest(const GraphTest& test) const {
  if (test.getInfectionTime() == -1)
    return this->getTestOutcome(test) == test.getOutcome();

  return m_infectionTimeHash.Len() < test.getInfectionTime() ||
    test.getOutcome() == this->getTestOutcome(test);
}

int GraphHypothesis::getInfectionTime(int nodeId) const
{
  if (m_infectionTimeHash.IsKey(nodeId))
    return m_infectionTimeHash.GetDat(nodeId);

  return -1;
}

bool GraphHypothesis::getTestOutcome(const GraphTest& test) const
{
  return m_infectionTimeHash.IsKey(test.getNodeId());
}

void GraphHypothesis::setWeight(double weight)
{
  m_weight = weight;
}

GraphHypothesisCluster::GraphHypothesisCluster(PUNGraph network,
                                               int sourceId,
                                               int maxHypothesis,
                                               double beta,
                                               double size)
  : m_network(network)
  , m_sourceId(sourceId)
  , m_beta(beta)
  , m_size(size)
{
  generateHypothesisCluster(maxHypothesis);
}

void GraphHypothesisCluster::printState()
{
  std::cout << "Cluster with source " << m_sourceId << " has " <<
    m_hypothesis.size() << " hypothesis left " << std::endl;
}

/**
 * Internal generator for all the hypothesis in the cluster.
 */
void GraphHypothesisCluster::generateHypothesisCluster(int maxHypothesis)
{
  // TODO(vcarbune): threadify
  for (int h = 0; h < maxHypothesis; h++)
    m_hypothesis.push_back(generateHypothesis());
}


/**
 * Generates one cascade, on top of the underlying network structure.
 * Extended from snap/examples/cascades.
 *
 * TODO(vcarbune): Do we want to sample the cascade (e.g. removing nodes?)
 */
GraphHypothesis GraphHypothesisCluster::generateHypothesis(
    bool isTrueHypothesis) const
{
  // TODO(vcarbune): discuss about these values.
  int cascadeSize = m_size * m_network->GetNodes();
  int runTimes = isTrueHypothesis ? cascadeSize : INITIAL_RUNS;

  TIntH nodeInfectionTime;

  // Add the source node (fixed for this cluster).
  nodeInfectionTime.AddDat(m_sourceId, nodeInfectionTime.Len());

  for (int run = 0; run < runTimes; run++) {
    for (TIntH::TIter I = nodeInfectionTime.BegI();
         I < nodeInfectionTime.EndI();
         I++) {
      const TUNGraph::TNodeI& crtIt = m_network->GetNI(I->Key());
      for (int neighbour = 0; neighbour < crtIt.GetOutDeg(); ++neighbour) {
        if (TInt::Rnd.GetUniDev() > m_beta) // Flip a coin!
            continue;

        int neighbourId = crtIt.GetOutNId(neighbour);
        if (nodeInfectionTime.IsKey((neighbourId)))
            continue;

        nodeInfectionTime.AddDat(neighbourId, nodeInfectionTime.Len());
        if (nodeInfectionTime.Len() == cascadeSize)
          return GraphHypothesis(nodeInfectionTime, 0);
      }
    }
  }

  // cout << (double) nodeInfectionTime.Len() / m_network->GetNodes() << endl;
  return GraphHypothesis(nodeInfectionTime, 0);
}

int GraphHypothesisCluster::countHypothesisConsistentWithTest (
    const GraphTest& test) const
{
  int total = 0;
  for (const GraphHypothesis& h : m_hypothesis)
    if (h.isConsistentWithTest(test))
      total++;

  return total;
}


/**
 * Eliminates all the cascades that are not possible,
 * considering the outcome of Test t to be true.
 */
void GraphHypothesisCluster::removeHypothesisInconsistentWithTest(const GraphTest& t)
{
  vector<GraphHypothesis> tmp;
  tmp.swap(m_hypothesis);

  for (size_t i = 0; i < tmp.size(); ++i)
    if (tmp[i].isConsistentWithTest(t))
      m_hypothesis.push_back(tmp[i]);
}

GraphHypothesis GraphHypothesisCluster::getRandomHypothesis() const
{
  return m_hypothesis[(long) (drand48() * m_hypothesis.size())];
}
