/**
 * Data structure implementation to represent a hypothesis cluster.
 * Author: Victor Carbune (vcarbune@ethz.ch)
 */

#include <cassert>
#include <iostream>

#include "hypothesis.h"

#define INITIAL_RUNS 5
#define EPS 0.90

GraphHypothesis::GraphHypothesis(TIntH hash)
  : m_infectionTimeHash(hash)
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
                                               double size,
                                               double weight)
  : m_network(network)
  , m_sourceId(sourceId)
  , m_beta(beta)
  , m_size(size)
  , m_hops(0)
  , m_weight(weight)
{
  generateHypothesisCluster(maxHypothesis);
  // updateConsistentHypothesisCount();
}

/**
 * Internal generator for all the hypothesis in the cluster.
 */
void GraphHypothesisCluster::generateHypothesisCluster(int maxHypothesis)
{
  // TODO(vcarbune): threadify
  for (int h = 0; h < maxHypothesis; h++) {
    m_hypothesis.push_back(generateHypothesis());
    m_hypothesis[h].weight = m_weight / maxHypothesis;
  }
}

/*
void GraphHypothesisCluster::updateConsistentHypothesisCount()
{
  m_crtConsistentHypothesis = 0;
  for (const GraphHypothesis &h : m_hypothesis)
    if (h.weight)
      m_crtConsistentHypothesis++;

  // cout << "Consistency percentage: " << (double) m_crtConsistentHypothesis / m_hypothesis.size() << endl;
}
*/

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
          return GraphHypothesis(nodeInfectionTime);
      }
    }
  }

  return GraphHypothesis(nodeInfectionTime);
}

void GraphHypothesisCluster::updateMassWithTest(const GraphTest& test)
{
  m_weight = 0;
  for (GraphHypothesis& h : m_hypothesis) {
    h.weight *= (h.isConsistentWithTest(test) ? EPS : (1-EPS));
    m_weight += h.weight;
  }

  // updateConsistentHypothesisCount();
}

double GraphHypothesisCluster::computeMassWithTest(const GraphTest& test) const
{
  double mass = 0.0;
  for (const GraphHypothesis& h : m_hypothesis)
    mass += h.weight * (h.isConsistentWithTest(test) ? EPS : (1-EPS));

  return mass;
}

void GraphHypothesisCluster::countConsistentHypothesis (
    const GraphTest& test, int *withTest, int *withPrevTests) const
{
  (*withPrevTests) += m_hypothesis.size();
  for (const GraphHypothesis& h : m_hypothesis)
    if (h.isConsistentWithTest(test))
      (*withTest)++;
}

GraphHypothesis GraphHypothesisCluster::getRandomHypothesis() const
{
  return m_hypothesis[(long) (drand48() * m_hypothesis.size())];
}
