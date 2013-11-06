/**
 * Data structure implementation to represent a hypothesis cluster.
 * Author: Victor Carbune (vcarbune@ethz.ch)
 */

#include <iostream>

#include "hypothesis.h"

GraphHypothesis::GraphHypothesis(PNGraph cascade, TIntH hash, double weight)
  : m_weight(weight)
  , m_infectionTimeHash(hash)
  , m_cascade(cascade)
{
}

bool GraphHypothesis::getTestOutcome(const GraphTest& test) const
{
  return m_infectionTimeHash.IsKey(test.getNodeId());
}

void GraphHypothesis::setWeight(double weight)
{
  m_weight = weight;
}

GraphHypothesisCluster::GraphHypothesisCluster(PUNGraph network, int sourceId, int maxHypothesis)
  : m_network(network)
  , m_sourceId(sourceId)
{
  generateHypothesisCluster(maxHypothesis);
}

void GraphHypothesisCluster::printState()
{
  std::cout << "Cluster with source " << m_sourceId << " has " <<
    m_hypothesis.size() << " hypothesis left " << std::endl;

  /*
  for (size_t i = 0; i < m_hypothesis.size(); ++i) {
    std::cout << "Cascade " << i << ": " <<
      (m_hypothesis[i].second.second->GetNodes()) << " (nodes) " <<
      (m_hypothesis[i].second.second->GetEdges()) << " (edges) " <<
      std::endl;
  }
  */
}

/**
 * Internal generator for all the hypothesis in the cluster.
 */
void GraphHypothesisCluster::generateHypothesisCluster(int maxHypothesis)
{
  for (int h = 0; h < maxHypothesis; ++h) {
    m_hypothesis.push_back(generateHypothesis());
  }
}


/**
 * Generates one cascade, on top of the underlying network structure.
 * Extended from snap/examples/cascades.
 *
 * TODO(vcarbune): Do we want to sample the cascade (e.g. removing nodes?)
 */
GraphHypothesis GraphHypothesisCluster::generateHypothesis()
{
  PUNGraph weaklyConnectedComponents = TSnap::GetMxWcc(m_network);

  // TODO(vcarbune): discuss about these values.
  double beta = 0.1;
  int cascadeSize = 0.4 * m_network->GetNodes();
  int runTimes = cascadeSize;

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
        if (casc->IsNode(neighbourId) || TInt::Rnd.GetUniDev() > beta)
          continue;

        casc->AddNode(neighbourId);
        nodeInfectionTime.AddDat(neighbourId, nodeInfectionTime.Len());
        casc->AddEdge(crtIt.GetId(), neighbourId);

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

GraphHypothesis GraphHypothesisCluster::getRandomHypothesis()
{
  return m_hypothesis[(long) (drand48() * m_hypothesis.size())];
}
