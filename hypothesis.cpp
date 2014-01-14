/**
 * Data structure implementation to represent a hypothesis cluster.
 * Author: Victor Carbune (vcarbune@ethz.ch)
 */

#include <cassert>
#include <iostream>

#include "hypothesis.h"

#define INITIAL_RUNS 5
#define EPS 0.05

GraphHypothesis::GraphHypothesis(TIntH hash)
  : m_infectionTimeHash(hash)
{
}

bool GraphHypothesis::isConsistentWithTest(const GraphTest& test) const {
  bool outcome = test.getOutcome() == this->getTestOutcome(test);

  if (test.getInfectionTime() == -1)
    return outcome;

  return outcome || m_infectionTimeHash.Len() < test.getInfectionTime();
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

  int cascadeSize = size * m_network->GetNodes();
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
        if (TInt::Rnd.GetUniDev() > beta) // Flip a coin!
            continue;

        unsigned int neighbourId = crtIt.GetOutNId(neighbour);
        if (nodeInfectionTime.IsKey((neighbourId)))
            continue;

        if (nodeCount)
          (*nodeCount)[neighbourId]++;

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
    h.weight *= (h.isConsistentWithTest(test) ? (1-EPS) : EPS);
    m_weight += h.weight;
  }
}

double GraphHypothesisCluster::computeMassWithTest(const GraphTest& test) const
{
  double mass = 0.0;
  for (const GraphHypothesis& h : m_hypothesis)
    mass += h.weight * (h.isConsistentWithTest(test) ? (1-EPS) : EPS);

  return mass;
}

GraphHypothesis GraphHypothesisCluster::getRandomHypothesis() const
{
  return m_hypothesis[(long) (drand48() * m_hypothesis.size())];
}
