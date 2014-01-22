/**
 * Data structure implementation to represent a hypothesis cluster.
 * Author: Victor Carbune (vcarbune@ethz.ch)
 */

#include <cassert>
#include <numeric>

#include "hypothesis.h"

#define INITIAL_RUNS 7

GraphHypothesis::GraphHypothesis(unsigned int nodes,
                                 unordered_map<int, int>& infectionTime)
{
  m_infectionHash.swap(infectionTime);
}

bool GraphHypothesis::isConsistentWithTest(const GraphTest& test) const {
  if (test.getInfectionTime() == INFECTED_UNDEFINED)
    return test.getOutcome() == this->getTestOutcome(test);

  if (test.getInfectionTime() == INFECTED_FALSE)
    return !this->getTestOutcome(test);

  unsigned int infectionTime =
    static_cast<unsigned int>(test.getInfectionTime());
  return this->getTestOutcome(test) || m_infectionHash.size() < infectionTime;
}

int GraphHypothesis::getInfectionTime(int nodeId) const
{
  if (m_infectionHash.find(nodeId) != m_infectionHash.end())
    return m_infectionHash.at(nodeId);

  return INFECTED_FALSE;
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
  m_nodeCount.resize(network->GetNodes());

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
          return GraphHypothesis(m_network->GetNodes(), infectionTime);
      }
    }
  }

  return GraphHypothesis(m_network->GetNodes(), infectionTime);
}

void GraphHypothesisCluster::updateMassWithTest(const GraphTest& test)
{
  m_weight = 0;
  for (GraphHypothesis& h : m_hypothesis) {
    h.weight *= (h.isConsistentWithTest(test) ? (1-EPS) : EPS);
    m_weight += h.weight;
  }

  /*
  double oldMass = m_weight;
  double pcnt = (oldMass - m_weight) / oldMass;
  if (fabs(pcnt - EPS) > 1E-3)
    cout << "Diff is " << pcnt << endl;
  */
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
