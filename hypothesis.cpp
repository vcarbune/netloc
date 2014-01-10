/**
 * Data structure implementation to represent a hypothesis cluster.
 * Author: Victor Carbune (vcarbune@ethz.ch)
 */

#include <cassert>
#include <iostream>

#include "hypothesis.h"

#define INITIAL_RUNS 10

GraphHypothesis::GraphHypothesis(TIntH hash, double weight)
  : m_infectionTimeHash(hash)
  , m_isMarkedAsInconsistent(false)
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

GraphHypothesisCluster::GraphHypothesisCluster(PUNGraph network,
                                               int sourceId,
                                               int maxHypothesis,
                                               double beta,
                                               double size)
  : m_network(network)
  , m_sourceId(sourceId)
  , m_beta(beta)
  , m_size(size)
  , m_hops(0)
  , m_weight(0)
{
  generateHypothesisCluster(maxHypothesis);
  updateConsistentHypothesisCount();
}

void GraphHypothesisCluster::setHopsFromSource(int hops)
{
  m_hops = hops;
}

void GraphHypothesisCluster::updateWeight(const vector<GraphTest>& A)
{
  int inconsistent = 0;
  for (const GraphTest& t : A) {
    for (const GraphHypothesis& h : m_hypothesis) {
      if (!h.isConsistentWithTest(t))
        inconsistent++;
    }
  }

  m_weight = pow((double) 0.25/0.75, inconsistent);
}

void GraphHypothesisCluster::setWeight(double weight)
{
  m_weight = weight;
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

void GraphHypothesisCluster::updateConsistentHypothesisCount()
{
  m_crtConsistentHypothesis = 0;
  for (const GraphHypothesis &h : m_hypothesis)
    if (!h.isMarkedAsInconsistent())
      m_crtConsistentHypothesis++;

  // cout << "Consistency percentage: " << (double) m_crtConsistentHypothesis / m_hypothesis.size() << endl;
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

  return GraphHypothesis(nodeInfectionTime, 0);
}

void GraphHypothesisCluster::countConsistentHypothesis (
    const GraphTest& test, int *withTest, int *withPrevTests) const
{
  for (const GraphHypothesis& h : m_hypothesis) {
    if (h.isMarkedAsInconsistent())
      continue;

    (*withPrevTests)++;
    if (h.isConsistentWithTest(test))
      (*withTest)++;
  }
}

void GraphHypothesisCluster::markInconsistentHypothesis(const GraphTest& t)
{
  for (GraphHypothesis& h : m_hypothesis)
    if (!h.isConsistentWithTest(t))
      h.markAsInconsistent();

  updateConsistentHypothesisCount();
}

GraphHypothesis GraphHypothesisCluster::getRandomHypothesis() const
{
  return m_hypothesis[(long) (drand48() * m_hypothesis.size())];
}
