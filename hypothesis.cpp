/**
 * Data structure implementation to represent a hypothesis cluster.
 */

#include "hypothesis.h"

HypothesisCluster::HypothesisCluster(PUNGraph network, int sourceId)
  : m_network(network)
  , m_sourceId(sourceId)
{
  generateCascades();
}

/**
 * Generates all the cascades among which search will be conducted,
 * on top of the underlying fixed network structure set in the constructor.
 */
void HypothesisCluster::generateCascades()
{
  // TODO(vcarbune): Cascade generation according to various models.
  map<int, int> cascade;
  for (int t = 1; t < 300; ++t) {
    cascade[(long) (drand48() * m_network->GetNodes())] = t;
  }

  // The source is always infected at first time.
  cascade[0] = m_sourceId;
  m_hypothesis.push_back(cascade);
}

/**
 * Eliminates all the cascades that are not possible,
 * considering the outcome of Test t to be true.
 */
void HypothesisCluster::ruleOutHypothesis(Test t)
{
  vector<Hypothesis>::iterator it;
  for (it = m_hypothesis.begin(); it != m_hypothesis.end();) {
    if ((*it).find(t.second) == (*it).end())
      ++it;
    else
      it = m_hypothesis.erase(it);
  }
}

Hypothesis HypothesisCluster::getRandomHypothesis()
{
  return m_hypothesis[(long) (drand48() * m_hypothesis.size())];
}
