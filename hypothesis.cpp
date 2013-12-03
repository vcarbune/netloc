/**
 * Data structure implementation to represent a hypothesis cluster.
 * Author: Victor Carbune (vcarbune@ethz.ch)
 */

#include <cassert>
#include <iostream>

#include "hypothesis.h"

#define INFECTION_RUNS 10

GraphHypothesis::GraphHypothesis(PNGraph cascade, TIntH hash, double weight)
  : m_weight(weight)
  , m_infectionTimeHash(hash)
  , m_cascade(cascade)
{
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

GraphHypothesisCluster::GraphHypothesisCluster(PUNGraph network, int sourceId,
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

GraphHypothesisCluster::GraphHypothesisCluster(PUNGraph network, int sourceId,
                                               int maxHypothesis)
  : GraphHypothesisCluster(network, sourceId, maxHypothesis, 0.05, 0.4)
{
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
  for (int h = 0; h < maxHypothesis; ++h)
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
  PUNGraph weaklyConnectedComponents = TSnap::GetMxWcc(m_network);

  // TODO(vcarbune): discuss about these values.
  int cascadeSize = m_size * m_network->GetNodes();
  int runTimes = isTrueHypothesis ? cascadeSize : INFECTION_RUNS;

  PNGraph casc = TNGraph::New();
  TIntH nodeInfectionTime;

  // Add the source node (fixed for this cluster).
  casc->AddNode(m_sourceId);
  nodeInfectionTime.AddDat(m_sourceId, nodeInfectionTime.Len());

  for (int run = 0; run < runTimes; run++) {
    TIntV cascade;
    casc->GetNIdV(cascade);

    for (int node = 0; node < cascade.Len(); ++node) {
      const TUNGraph::TNodeI crtIt = m_network->GetNI(cascade[node]);

      for (int neighbour = 0; neighbour < crtIt.GetOutDeg(); ++neighbour) {
        int neighbourId = crtIt.GetOutNId(neighbour);

        // If the node is already in the cascade, or random activation fails,
        // move to the next neighbour.
        if (casc->IsNode(neighbourId) || TInt::Rnd.GetUniDev() > m_beta)
          continue;


        casc->AddNode(neighbourId);
        casc->AddEdge(crtIt.GetId(), neighbourId);

        assert(!nodeInfectionTime.IsKey(neighbourId));
        nodeInfectionTime.AddDat(neighbourId, nodeInfectionTime.Len());

        if (casc->GetNodes() == cascadeSize)
          return GraphHypothesis(casc, nodeInfectionTime, 0);
      }
    }
  }

  return GraphHypothesis(casc, nodeInfectionTime, 0);
}

int GraphHypothesisCluster::countHypothesisConsistentWithTest (
    const GraphTest& test) const
{
  int total = 0;

  for (GraphHypothesis h : m_hypothesis)
    if (h.getTestOutcome(test) == test.getOutcome())
      total++;

  return total;
}


/**
 * Eliminates all the cascades that are not possible,
 * considering the outcome of Test t to be true.
 */
void GraphHypothesisCluster::removeHypothesisInconsistentWithTest(const GraphTest& t)
{
  vector<GraphHypothesis>::iterator it;
  for (it = m_hypothesis.begin(); it != m_hypothesis.end();) {
    if (it->getTestOutcome(t) == t.getOutcome())
      ++it;
    else
      it = m_hypothesis.erase(it);
  }
}

GraphHypothesis GraphHypothesisCluster::getRandomHypothesis() const
{
  return m_hypothesis[(long) (drand48() * m_hypothesis.size())];
}
