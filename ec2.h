/**
 * This is an abstract implementation of EC2.
 * Source: http://arxiv.org/pdf/1010.3091v1.pdf
 *
 * Author: Victor Carbune (vcarbune@ethz.ch)
 */

#ifndef EC2_H_
#define EC2_H_

#include <iostream>
#include <limits>

#include <queue>
#include <vector>

// Abstract classes that can be implemented in order to use runEC2.

class Test {
  public:
    void setOutcome(bool outcome) { m_outcome = outcome; }
    bool getOutcome() const { return m_outcome; }

    void setScore(double score) { m_score = score; }
    double getScore() const { return m_score; }

    virtual bool operator==(const Test& o) const { return false; }
    virtual bool operator<(const Test& o) const {
      return getScore() < o.getScore();
    }

  private:
    double m_score;
    bool m_outcome;
};

class Hypothesis {
  public:
    virtual bool getTestOutcome(const Test&) const = 0;
    virtual void setWeight(double) = 0;
};

class HypothesisCluster {
  public:
    virtual int countHypothesisConsistentWithTest(const Test&) const = 0;
    virtual int countHypothesisAvailable() const = 0;

    virtual void removeHypothesisInconsistentWithTest(const Test&) = 0;
    virtual void printState() const = 0;
};

class TestCompareFunction {
  public:
    bool operator() (const Test& p, const Test& q) const {
      return p.getScore() < q.getScore();
    }
};


// Actual EC2 template implementation.

template <class TTest, class THypothesisCluster>
inline void rescoreTest(TTest& test,
                 const std::vector<THypothesisCluster>& clusters,
                 int countTotalHypothesis)
{
  int countConsistentHypothesis = 0;

  // Count hypothesis matching test.
  for (const THypothesisCluster& cluster : clusters)
    countConsistentHypothesis +=
      cluster.countHypothesisConsistentWithTest(test);

  // Score if test outcome is True
  double positiveScore =
    ((double) countConsistentHypothesis / countTotalHypothesis) *
    ((double) 1 / countConsistentHypothesis);

  // Score if test outcome is False
  double negativeScore =
    ((double) (1 - countConsistentHypothesis)) *
    ((double) 1 / (countTotalHypothesis - countConsistentHypothesis));

  double score = positiveScore + negativeScore -
    ((double) 1 / countTotalHypothesis);

  test.setScore(score);
}

template <class TTest, class THypothesisCluster>
void rescoreTests(std::vector<TTest>& tests,
                  const std::vector<THypothesisCluster>& clusters)
{
  int countTotalHypothesis = 0;
  for (const THypothesisCluster& cluster : clusters)
      countTotalHypothesis += cluster.countHypothesisAvailable();

  double prevScore = tests.front().getScore();
  rescoreTest(tests.front(), clusters, countTotalHypothesis);
  if (prevScore <= tests.front().getScore()) {
    tests.front().setScore(prevScore);
    return;
  }

  for (TTest& t : tests)
    rescoreTest(t, clusters, countTotalHypothesis);

  std::make_heap<typename std::vector<TTest>::iterator, TestCompareFunction>(
      tests.begin(), tests.end(), TestCompareFunction());
}

template <class TTest, class THypothesisCluster, class THypothesis>
TTest runOneEC2Step(std::vector<TTest>& tests,
                    std::vector<THypothesisCluster>& clusters,
                    const THypothesis& realization)
{
  rescoreTests(tests, clusters);

  TTest& test = tests.front();
  test.setOutcome(realization.getTestOutcome(test));

  for (std::size_t i = 0; i < clusters.size(); ++i)
    clusters[i].removeHypothesisInconsistentWithTest(test);

  std::pop_heap<typename std::vector<TTest>::iterator, TestCompareFunction>(
      tests.begin(), tests.end(), TestCompareFunction());
  tests.pop_back();

  return test;
}

template <class TTest, class THypothesisCluster, class THypothesis>
size_t runEC2(std::vector<TTest>& tests,
              std::vector<THypothesisCluster>& clusters,
              const THypothesis& realization)
{
  typename std::vector<TTest>::iterator it;
  std::vector<TTest> testRunOrder;

  rescoreTests(tests, clusters);

  int clustersLeft = clusters.size();
  while (!tests.empty() && clustersLeft != 1) {
    TTest t = runOneEC2Step<TTest, THypothesisCluster, THypothesis>(
        tests, clusters, realization);
    testRunOrder.push_back(t);

    clustersLeft = 0;
    for (std::size_t i = 0; i < clusters.size(); ++i)
      if (clusters[i].countHypothesisAvailable())
        clustersLeft++;
  }

  assert (clustersLeft == 1 || tests.empty());
  return testRunOrder.size();
}


#endif // EC2_H_
