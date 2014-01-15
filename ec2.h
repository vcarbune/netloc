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
#include <utility>

#define WORK_THREADS 8
#define MASS_THRESHOLD 0.07

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
  // TODO(vcarbune): Redefine interface...
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
  int totalHypothesis = 0;
  int testConsistentHypothesis = 0;

  double positiveMass = 0.0;
  double positiveDiagonalMass = 0.0;
  double negativeMass = 0.0;
  double negativeDiagonalMass = 0.0;

  for (const THypothesisCluster& cluster : clusters) {
    testConsistentHypothesis += cluster.getNodeCount(test.getNodeId());
    totalHypothesis += cluster.getTotalHypothesis();

    std::pair<double, double> mass = cluster.computeMassWithTest(test);

    positiveMass += mass.first;
    negativeMass += mass.second;

    positiveDiagonalMass -= mass.first * mass.first;
    negativeDiagonalMass -= mass.second * mass.second;
  }

  positiveMass = positiveMass * positiveMass - positiveDiagonalMass;
  negativeMass = negativeMass * negativeMass - negativeDiagonalMass;

  double testPositivePb = (double) testConsistentHypothesis / totalHypothesis;
  test.setScore(testPositivePb * positiveMass + (1 - testPositivePb) * negativeMass);
}

template <class TTest, class THypothesisCluster>
void rescoreTests(std::vector<TTest>& tests,
                  const std::vector<THypothesisCluster>& clusters)
{
  double prevScore = tests.front().getScore();
  rescoreTest(tests.front(), clusters);

  double relativeChange = prevScore - tests.front().getScore();
  relativeChange /= prevScore;

  if (relativeChange < 0.3) {
      tests.front().setScore(prevScore);
      return;
  }

  std::vector<std::thread> threads;
  for (size_t test = 0; test < tests.size(); test += WORK_THREADS) {
    for (int t = 0; t < WORK_THREADS && test + t < tests.size(); ++t)
      threads.push_back(std::thread(rescoreTest<TTest, THypothesisCluster>,
            std::ref(tests[test + t]), std::ref(clusters)));
    for (std::thread& t : threads)
      t.join();
    threads.clear();
  }

  std::make_heap<typename std::vector<TTest>::iterator, TestCompareFunction>(
      tests.begin(), tests.end(), TestCompareFunction());
}

template <class TTest, class THypothesisCluster, class THypothesis>
TTest runOneEC2Step(std::vector<TTest>& tests,
                    std::vector<THypothesisCluster>& clusters,
                    const THypothesis& realization)
{
  rescoreTests(tests, clusters);

  TTest test = tests.front();
  test.setOutcome(realization.getTestOutcome(test));
  test.setInfectionTime(realization.getInfectionTime(test.getNodeId()));

  std::vector<std::thread> threads;
  for (std::size_t i = 0; i < clusters.size(); i += WORK_THREADS) {
    for (int t = 0; t < WORK_THREADS && i + t < clusters.size(); ++t)
      threads.push_back(std::thread(&THypothesisCluster::updateMassWithTest,
          &clusters[i+t], std::ref(test)));
    for (std::thread& t : threads)
      t.join();
    threads.clear();
  }

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

  double mass = 0.0;
  bool shouldStop = false;

  while (!tests.empty() && !shouldStop) {
    TTest t = runOneEC2Step<TTest, THypothesisCluster, THypothesis>(
        tests, clusters, realization);
    testRunOrder.push_back(t);

    mass = 0.0;
    for (unsigned i = 0; i < clusters.size(); ++i)
      mass += clusters[i].getMass();

    for (THypothesisCluster& cluster : clusters)
      if (cluster.getWeight() / mass > MASS_THRESHOLD)
        shouldStop = true;
  }

  // Normalize weights and probabilities
  for (THypothesisCluster& cluster : clusters)
    cluster.normalizeWeight(mass);

  // The score is the number of tests used.
  return testRunOrder.size();
}

#endif // EC2_H_
