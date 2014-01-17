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
TTest rescoreTests(std::vector<TTest>& tests,
                  const std::vector<THypothesisCluster>& clusters)
{
  int count = 0;
  std::cout << ":: " << tests.front().getNodeId() << std::endl;
  double score = tests.front().getScore();
  do {
    // Remove the top test from the heap.
    TTest crtTop = tests.front();
    std::pop_heap<typename std::vector<TTest>::iterator, TestCompareFunction>(
        tests.begin(), tests.end(), TestCompareFunction());
    tests.pop_back();

    // Exit early if it's the last element in the heap.
    if (!tests.size())
      return crtTop;

    // Recompute its score and keep it if it stays on top.
    rescoreTest(crtTop, clusters);
    if (crtTop.getScore() >= tests.front().getScore()) {
      std::cout << ":: " << crtTop.getNodeId() << std::endl;
      std::cout << (score - crtTop.getScore()) / score << std::endl;
      std::cout << "Pushed " << count << " elems back to heap... " << std::endl;
      return crtTop;
    }

    // Otherwise push it back to the heap.
    tests.push_back(crtTop); count++;
    std::push_heap<typename std::vector<TTest>::iterator, TestCompareFunction>(
        tests.begin(), tests.end(), TestCompareFunction());
  } while (true);

  /*
  TTest secondTop = tests.front();
  // std::cout << prevScore << " now " << crtTop.getScore() << " - " << crtTop.getNodeId() << " " << tests.front().getScore() << " - " << tests.front().getNodeId() << std::endl;
  if (crtTop.getScore() >= tests.front().getScore()) {
    std::cout << "yeaaaah" << std::endl;
    return crtTop;
  } else {
    crtTop.setScore(prevScore);
    tests.push_back(crtTop);
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

  TTest newTop = tests.front();
  std::pop_heap<typename std::vector<TTest>::iterator, TestCompareFunction>(
      tests.begin(), tests.end(), TestCompareFunction());
  tests.pop_back();

  // std::cout << newTop.getScore() << " - " << newTop.getNodeId() << " " << tests.front().getScore() << " - " << tests.front().getNodeId() << std::endl;
  if (newTop == crtTop && secondTop == tests.front())
    std::cout << "same order was kept.. " << std::endl;
  return newTop;
  */
}

template <class TTest, class THypothesisCluster, class THypothesis>
TTest runOneEC2Step(std::vector<TTest>& tests,
                    std::vector<THypothesisCluster>& clusters,
                    const THypothesis& realization)
{
  TTest test = rescoreTests(tests, clusters);
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

  // Initially assign score to all tests and create the heap from the vector.
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

  while (!tests.empty() && !shouldStop) {
    TTest t = runOneEC2Step<TTest, THypothesisCluster, THypothesis>(
        tests, clusters, realization);
    testRunOrder.push_back(t);

    std::cout << testRunOrder.size() << ". " << t.getNodeId() << " --> " << t.getScore() << std::endl;

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
