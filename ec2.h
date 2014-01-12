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
#include <cmath>
#include <cstdio>

#include <queue>
#include <vector>
#include <thread>

#define WORK_THREADS 8

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
};

class HypothesisCluster {
  public:
    virtual int countConsistentHypothesis() const = 0;
    virtual void countConsistentHypothesis(const Test&, int*, int*) const = 0;

    virtual void markInconsistentHypothesis(const Test&) = 0;
    virtual void setWeight(double) = 0;
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
    const std::vector<THypothesisCluster>& clusters)
{
  int prevConsistentHypothesis = 0;
  int testConsistentHypothesis = 0;

  double positiveMass = 0.0;
  double positiveDiagonalMass = 0.0;
  double negativeMass = 0.0;
  double negativeDiagonalMass = 0.0;

  for (const THypothesisCluster& cluster : clusters) {
    cluster.countConsistentHypothesis(test,
        &testConsistentHypothesis, &prevConsistentHypothesis);

    test.setOutcome(true);
    double positiveClusterMass = cluster.computeMassWithTest(test);

    test.setOutcome(false);
    double negativeClusterMass = cluster.computeMassWithTest(test);

    positiveMass += positiveClusterMass;
    negativeMass += negativeClusterMass;

    positiveDiagonalMass -= positiveClusterMass * positiveClusterMass;
    negativeDiagonalMass -= negativeClusterMass * negativeClusterMass;
  }

  positiveMass = positiveMass * positiveMass - positiveDiagonalMass;
  negativeMass = negativeMass * negativeMass - negativeDiagonalMass;

  double testPositivePb = (double) testConsistentHypothesis / prevConsistentHypothesis;
  double score = testPositivePb * positiveMass +
    (1 - testPositivePb) * negativeMass;

  test.setScore(score);
}

template <class TTest, class THypothesisCluster>
void rescoreTests(std::vector<TTest>& tests,
                  const std::vector<THypothesisCluster>& clusters)
{
  double prevScore = tests.front().getScore();
  rescoreTest(tests.front(), clusters);

  double relativeChange = prevScore - tests.front().getScore();
  relativeChange /= prevScore;

  if (relativeChange < 10 * 1E-2) {
      tests.front().setScore(prevScore);
      return;
  }

  /*
  std::cout << prevScore << " now is " << tests.front().getScore() <<
    " " << tests.front().getNodeId() << std::endl;
  */

  size_t test = 0;
  for (; test < tests.size(); test += WORK_THREADS) {
    std::vector<std::thread> threads;

    for (int t = 0; t < WORK_THREADS && test + t < tests.size(); ++t)
      threads.push_back(std::thread(rescoreTest<TTest, THypothesisCluster>,
            std::ref(tests[test + t]), std::ref(clusters)));

    for (std::thread& t : threads)
      t.join();
  }

  test -= WORK_THREADS;
  for (; test < tests.size(); test++)
    rescoreTest(tests[test], clusters);

  std::make_heap<typename std::vector<TTest>::iterator, TestCompareFunction>(
      tests.begin(), tests.end(), TestCompareFunction());
}

template<class THypothesisCluster, class TTest>
void threadWrapperClusterCleanup(THypothesisCluster& cluster,
    const TTest& test, bool* isNewlyEmptied)
{
  // TODO(vcarbune): remove this, hypothesis will never dissapear.
  *isNewlyEmptied = false;
  if (!cluster.countConsistentHypothesis())
    return;

  cluster.updateMassWithTest(test);

  // TODO(vcarbune): remove this, hypothesis will never dissapear.
  if (!cluster.countConsistentHypothesis())
    *isNewlyEmptied = true;
}

template <class TTest, class THypothesisCluster, class THypothesis>
TTest runOneEC2Step(std::vector<TTest>& tests,
                    std::vector<THypothesisCluster>& clusters,
                    std::vector<int>& removedClusters,
                    const THypothesis& realization)
{
  rescoreTests(tests, clusters);

  TTest test = tests.front();
  test.setOutcome(realization.getTestOutcome(test));
  test.setInfectionTime(realization.getInfectionTime(test.getNodeId()));

  /* Original: remove inconsistent hypothesis from all clusters. */
  bool emptiedClusters[clusters.size()];
  std::size_t i = 0;
  for (; i < clusters.size(); i += WORK_THREADS) {
    std::vector<std::thread> threads;

    for (int t = 0; t < WORK_THREADS && i + t < clusters.size(); ++t)
      threads.push_back(
          std::thread(threadWrapperClusterCleanup<THypothesisCluster, TTest>,
            std::ref(clusters[i+t]), std::ref(test), &emptiedClusters[i+t]));

    for (std::thread& t : threads)
      t.join();
  }

  for (; i < clusters.size(); i++) {
    threadWrapperClusterCleanup<THypothesisCluster, TTest>(clusters[i], test,
        &emptiedClusters[i]);
  }

  for (i = 0; i < clusters.size(); i++) {
    if (emptiedClusters[i])
      removedClusters.push_back(clusters[i].getSource());
  }

  std::pop_heap<typename std::vector<TTest>::iterator, TestCompareFunction>(
      tests.begin(), tests.end(), TestCompareFunction());
  tests.pop_back();

  return test;
}

template <class TTest, class THypothesisCluster, class THypothesis>
size_t runEC2(std::vector<TTest>& tests,
              std::vector<THypothesisCluster>& clusters,
              const THypothesis& realization,
              int topN,
              std::vector<int>& topClusters)
{
  typename std::vector<TTest>::iterator it;
  std::vector<TTest> testRunOrder;

  while (!tests.empty() && topClusters.size() < clusters.size()) {
    TTest t = runOneEC2Step<TTest, THypothesisCluster, THypothesis>(
        tests, clusters, topClusters, realization);
    testRunOrder.push_back(t);

    /*
    std::cout << testRunOrder.size() << ". " << t.getScore() << " - id:" <<
      t.getNodeId() << std::endl;
    */

    if (t.getScore() == 0)
      break;
  }

  /*
  sort(clusters.begin(), clusters.end());
  for (THypothesisCluster& cluster : clusters) {
    std::cout << cluster.getSource() << " --> " << cluster.getWeight() <<
      std::endl;
  }

  for (int i = 0; i < topN; ++i) {
    topClusters.push_back(clusters[i].getSource());
    std::cout << clusters[i].getSource() << ", " << clusters[i].getWeight() << " " << std::endl;
  }
  */

  return testRunOrder.size();
}

#endif // EC2_H_
