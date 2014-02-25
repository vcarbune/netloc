/**
 * This is an implementation of EC2.
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

#define DBG 1
#define WORK_THREADS 16

#define EPS 0.05

#undef min
#undef max

// Abstract classes that can be implemented in order to use runEC2.
// TODO(vcarbune): Not really very useful anymore.

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


template <class TTest, class THypothesisCluster>
inline void rescoreTestVoI(TTest& test,
                           const std::vector<THypothesisCluster>& clusters)
{
  double testConsistentHypothesis = 0.0;
  int totalHypothesis = 0;

  for (const THypothesisCluster& cluster : clusters) {
    testConsistentHypothesis += cluster.getNodeCount(test.getNodeId());
    totalHypothesis += cluster.getTotalHypothesis();
  }

  double testPositivePb = (double) testConsistentHypothesis / totalHypothesis;
  double score = 0.0;
  for (const THypothesisCluster& cluster : clusters) {
    std::pair<double, double> mass = cluster.computeMassWithTest(test);
    double cMass = testPositivePb * mass.first +
        (1 - testPositivePb) * mass.second;
    score = std::max(cMass, score);
  }
  test.setScore(score);
}

template <class TTest, class THypothesisCluster>
inline void rescoreTestGBS(TTest& test,
                           const std::vector<THypothesisCluster>& clusters)
{
  double testConsistentHypothesis = 0.0;
  int totalHypothesis = 0;

  for (const THypothesisCluster& cluster : clusters) {
    testConsistentHypothesis += cluster.getNodeCount(test.getNodeId());
    totalHypothesis += cluster.getTotalHypothesis();
  }

  double totalMass = 0.0;
  double testPositivePb = (double) testConsistentHypothesis / totalHypothesis;
  for (const THypothesisCluster& cluster : clusters) {
    std::pair<double, double> mass = cluster.computeMassWithTest(test);
    double cMass = testPositivePb * mass.first +
        (1 - testPositivePb) * mass.second;
    totalMass += cMass;
  }

  test.setScore(-totalMass);
}

template <class TTest, class THypothesisCluster>
inline void rescoreTestEC2(TTest& test,
                           const std::vector<THypothesisCluster>& clusters)
{
  int totalHypothesis = 0;
  int testConsistentHypothesis = 0;

  double positiveMass = 0.0;
  double positiveDiagonalMass = 0.0;
  double negativeMass = 0.0;
  double negativeDiagonalMass = 0.0;
  double currentMass = 0.0;
  double currentDiagonalMass = 0.0;

  for (const THypothesisCluster& cluster : clusters) {
    testConsistentHypothesis += cluster.getNodeCount(test.getNodeId());
    totalHypothesis += cluster.getTotalHypothesis();

    std::pair<double, double> mass = cluster.computeMassWithTest(test);

    positiveMass += mass.first;
    negativeMass += mass.second;

    positiveDiagonalMass += mass.first * mass.first;
    negativeDiagonalMass += mass.second * mass.second;

    double crtMass = cluster.getWeight();
    currentMass += crtMass;
    currentDiagonalMass += crtMass * crtMass;
  }

  currentMass = currentMass * currentMass - currentDiagonalMass;
  positiveMass = positiveMass * positiveMass - positiveDiagonalMass;
  negativeMass = negativeMass * negativeMass - negativeDiagonalMass;

  double testPositivePb = (double) testConsistentHypothesis / totalHypothesis;
  double expectedMass =
      testPositivePb * positiveMass + (1 - testPositivePb) * negativeMass;

  std::cout << totalHypothesis << std::endl;
  test.setScore(currentMass - expectedMass);
}

template <class TTest, class THypothesisCluster>
inline void rescoreTest(TTest& test,
                        const std::vector<THypothesisCluster>& clusters,
                        const int objType)
{
  switch(objType) {
    case 1:
      return rescoreTestGBS(test, clusters);
    case 2:
      return rescoreTestVoI(test, clusters);
    default:
      return rescoreTestEC2(test, clusters);
  }
}


template <class TTest, class THypothesisCluster>
TTest lazyRescoreTests(std::vector<TTest>& tests,
                       const std::vector<THypothesisCluster>& clusters,
                       const int objType)
{
  TestCompareFunction tstCmpFcn;
#if DBG
  int count = 0;
#endif
  do {
    // Remove the top test from the heap.
    TTest crtTop = tests.front();
    std::pop_heap<typename std::vector<TTest>::iterator, TestCompareFunction>(
        tests.begin(), tests.end(), TestCompareFunction());
    tests.pop_back();

#if DBG
    std::cout << "Top: " << crtTop.getScore() << " Next: " <<
      tests.front().getScore() << std::endl;
#endif

    // Exit early if it's the last element in the heap.
    if (!tests.size())
      return crtTop;

    // Recompute its score and keep it if it stays on top.
    rescoreTest(crtTop, clusters, objType);

#if DBG
    std::cout << "Rescored Top: " << crtTop.getScore() << std::endl;
#endif

    if (tstCmpFcn(tests.front(), crtTop)) {
#if DBG
      std::cout << "Pushed " << count << " elems back to heap... " << std::endl;
#endif
      return crtTop;
    }

    // Otherwise push it back to the heap.
    tests.push_back(crtTop);
#if DBG
    count++;
#endif
    std::push_heap<typename std::vector<TTest>::iterator, TestCompareFunction>(
        tests.begin(), tests.end(), TestCompareFunction());
  } while (true);

  std::cout << "SHOULD NOT BE REACHED" << std::endl;
  return tests.front();
}

template <class TTest, class THypothesisCluster>
TTest completeRescoreTests(std::vector<TTest>& tests,
                           const std::vector<THypothesisCluster>& clusters,
                           const int objType)
{
  TestCompareFunction tstCmpFcn;

  // Rescore all the tests on multiple threads.
  std::vector<std::thread> threads;
  for (size_t test = 0; test < tests.size(); test += WORK_THREADS) {
    for (int t = 0; t < WORK_THREADS && test + t < tests.size(); ++t)
      threads.push_back(std::thread(rescoreTest<TTest, THypothesisCluster>,
          std::ref(tests[test + t]), std::ref(clusters), objType));
    for (std::thread& t : threads)
      t.join();
    threads.clear();
  }

  // Just get the top scored element from the vector (no heap required).
  TTest top = tests.front();
  for (const TTest& test : tests)
    if (tstCmpFcn(top, test))
      top = test;

  tests.erase(std::find(tests.begin(), tests.end(), top));
  return top;
}

template <class TTest, class THypothesisCluster, class THypothesis>
TTest runOneEC2Step(std::vector<TTest>& tests,
                    std::vector<THypothesisCluster>& clusters,
                    const THypothesis& realization,
                    bool lazyEval,
                    const int objType)
{
  // Select the best test at this point.
  TTest test = lazyEval ? lazyRescoreTests(tests, clusters, objType) :
      completeRescoreTests(tests, clusters, objType);
  test.setOutcome(realization.getTestOutcome(test));
  test.setInfectionTime(realization.getInfectionTime(test.getNodeId()));

  // Update mass of each cluster considering the test outcome.
  std::vector<std::thread> threads;
  for (std::size_t i = 0; i < clusters.size(); i += WORK_THREADS) {
    for (int t = 0; t < WORK_THREADS && i + t < clusters.size(); ++t)
      threads.push_back(std::thread(&THypothesisCluster::updateMassWithTest,
          &clusters[i+t], std::ref(test)));
    for (std::thread& t : threads)
      t.join();
    threads.clear();
  }

  for (TTest& test : tests)
    test.setScore((1-EPS) * (1-EPS) * test.getScore());

  return test;
}

template <class TTest, class THypothesisCluster, class THypothesis>
double runEC2(std::vector<TTest>& tests,
              std::vector<THypothesisCluster>& clusters,
              const THypothesis& realization,
              bool lazyEval,
              const double massThreshold,
              const double testThreshold,
              const int objType)
{
  typename std::vector<TTest>::iterator it;
  std::vector<TTest> testRunOrder;

  unsigned int totalTests = tests.size();

  // Initially assign score to all tests and create the heap from the vector.
  std::vector<std::thread> threads;
  for (size_t test = 0; test < tests.size(); test += WORK_THREADS) {
    for (int t = 0; t < WORK_THREADS && test + t < tests.size(); ++t)
      threads.push_back(std::thread(rescoreTest<TTest, THypothesisCluster>,
            std::ref(tests[test + t]), std::ref(clusters), objType));
    for (std::thread& t : threads)
      t.join();
    threads.clear();
  }

  std::make_heap<typename std::vector<TTest>::iterator, TestCompareFunction>(
      tests.begin(), tests.end(), TestCompareFunction());

  double mass = 0.0;
  double percentageMass = 0.0;
  double percentageTests = 0.0;

  while (percentageMass < massThreshold && percentageTests < testThreshold) {
    TTest t = runOneEC2Step<TTest, THypothesisCluster, THypothesis>(
        tests, clusters, realization, lazyEval, objType);
    testRunOrder.push_back(t);

#if DBG
    std::cout << testRunOrder.size() << ". " << t.getNodeId() << " --> " << t.getScore() << std::endl;
#endif

    mass = 0.0;
    for (unsigned i = 0; i < clusters.size(); ++i)
      mass += clusters[i].getMass();

    percentageMass = 0.0;
    for (THypothesisCluster& cluster : clusters) {
      percentageMass = std::max(cluster.getWeight() / mass, percentageMass);
    }

    percentageTests = (double) testRunOrder.size() / totalTests;
  }

  return 100 * percentageTests;
}

#endif // EC2_H_
