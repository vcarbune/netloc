/**
 * Data structure implementation to represent a hypothesis cluster.
 */

#include <iostream>

#include "hypothesis.h"


HypothesisCluster::HypothesisCluster(PUNGraph network, int sourceId, int maxHypothesis)
  : m_network(network)
  , m_sourceId(sourceId)
{
  generateHypothesisCluster(maxHypothesis);
}

void HypothesisCluster::printState()
{
  std::cout << "Cluster with source " << m_sourceId << std::endl;
  for (size_t i = 0; i < m_hypothesis.size(); ++i) {
    std::cout << "Cascade " << i << ": " <<
      (m_hypothesis[i].second.second->GetNodes()) << " (nodes) " <<
      (m_hypothesis[i].second.second->GetEdges()) << " (edges) " <<
      std::endl;
  }
}

/**
 * Internal generator for all the hypothesis in the cluster.
 */
void HypothesisCluster::generateHypothesisCluster(int maxHypothesis)
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
Hypothesis HypothesisCluster::generateHypothesis()
{
  PUNGraph weaklyConnectedComponents = TSnap::GetMxWcc(m_network);

  // TODO(vcarbune): discuss about these values.
  double beta = 0.1;
  int cascadeSize = 0.4 * m_network->GetNodes();
  int runTimes = cascadeSize;

  Hypothesis h;
  h.second.second = TNGraph::New();

  PNGraph casc = h.second.second;
  TIntH& nodeInfectionTime = h.second.first;

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
          return h;
      }
    }
  }

  return h;
}

/**
 * Eliminates all the cascades that are not possible,
 * considering the outcome of Test t to be true.
 */
void HypothesisCluster::ruleOutHypothesis(Test t)
{
  vector<Hypothesis>::iterator it;
  for (it = m_hypothesis.begin(); it != m_hypothesis.end();) {
    TIntH& nodeInfectionTime = (*it).second.first;
    if (!nodeInfectionTime.IsKey(t.second))
      ++it;
    else
      it = m_hypothesis.erase(it);
  }
}

Hypothesis HypothesisCluster::getRandomHypothesis()
{
  return m_hypothesis[(long) (drand48() * m_hypothesis.size())];
}
