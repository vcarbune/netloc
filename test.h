/**
 * Data structure used to abstract away a test instance.
 * Author: Victor Carbune (vcarbune@ethz.ch)
 */

#ifndef GRAPH_TEST_H_
#define GRAPH_TEST_H_

#include <limits>
#include <utility>

#define INFECTED_FALSE        -2
#define INFECTED_UNDEFINED    -1

/**
 * Base Class (Generic)
 */
class Test {
  public:
    Test()
      : m_score(std::numeric_limits<double>::max())
      , m_outcome(true) {}

    void setOutcome(bool outcome) { m_outcome = outcome; }
    bool getOutcome() const { return m_outcome; }

    void setScore(double score) { m_score = score; }
    double getScore() const { return m_score; }

    virtual bool operator==(const Test& o) const { return false; }

  private:
    double m_score;
    bool m_outcome;
};

class TestCompareFunction {
  public:
    bool operator() (const Test& p, const Test& q) const {
      return p.getScore() < q.getScore();
    }
};

/**
 * GraphTest is a particular instance used to solve
 * the source localization problem.
 */

class GraphTest : public Test {
  public:
    explicit GraphTest(int source) : Test()
        , m_nodeId(source)
        , m_infectionTime(INFECTED_UNDEFINED) {}

    int getNodeId() const { return m_nodeId; }
    double getInfectionTime() const { return m_infectionTime; }

    void setInfectionTime(double infectionTime) {
      m_infectionTime = infectionTime;
    }

    bool operator==(const GraphTest& o) {
      return m_nodeId == o.m_nodeId;
    }

  private:
    int m_nodeId;
    double m_infectionTime;
};

#endif // GRAPH_TEST_H_
