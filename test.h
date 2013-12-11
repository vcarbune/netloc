/**
 * Data structure used to abstract away a test instance.
 * Author: Victor Carbune (vcarbune@ethz.ch)
 */

#ifndef GRAPH_TEST_H_
#define GRAPH_TEST_H_

#include <utility>

#include "ec2.h"

using namespace std;

/**
 * One test is a simple container of score and nodeId. The meaning of the two items is:
 * Score  - the expected reduction in weight from the prior to the posterior
 *          distribution, given that the test runs.
 * Node   - the node tested (and further assumed to be infected)
 */

class GraphTest : public Test {
  public:
    explicit GraphTest(int source) : Test()
        , m_nodeId(source)
        , m_infectionTime(-1)
    {
      setOutcome(true);
    }

    int getNodeId() const { return m_nodeId; }
    int getInfectionTime() const { return m_infectionTime; }

    void setInfectionTime(int infectionTime) {
      m_infectionTime = infectionTime;
    }

    bool operator==(const GraphTest& o) {
      return m_nodeId == o.m_nodeId;
    }

  private:
    int m_nodeId;
    int m_infectionTime;
};

#endif // GRAPH_TEST_H_
