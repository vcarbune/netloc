/**
 * Maximum-Likelihood Solution (see header file comments for details).
 * Author (Victor Carbune): vcarbune@ethz.ch
 */

#include "epfl_solver.h"

// SNAP redefines max and min!
#undef max
#undef min

#include <algorithm>
#include <iostream>

#include "Eigen/Dense"

using namespace std;

EPFLSolver::EPFLSolver(PUNGraph network, SimConfig config, double pcnt)
  : m_network(network)
  , m_config(config)
{
  // Keep top % of nodes as observers based on highest-degree policy.
  int observers = pcnt * m_network->GetNodes();

  // Sort by degree (using default pair comparator).
  vector<pair<int, int>> degrees;
  for (int node = 0; node < m_network->GetNodes(); ++node)
    degrees.push_back(pair<int, int>(m_network->GetNI(node).GetOutDeg(), node));
  sort(degrees.begin(), degrees.end(), greater<pair<int, int>>());

  // Store observer nodes locally.
  for (int node = 0; node < observers; ++node)
    m_observerNodes.push_back(degrees[node].second);
}

void EPFLSolver::setObserverList(const vector<int>& observers, double pcnt)
{
  int observersCount = pcnt * m_network->GetNodes();
  m_observerNodes.clear();
  for (int node = 0; node < observersCount; ++node) {
    m_observerNodes.push_back(observers[node]);
    // cout << observers[node] << " ";
  }
  // cout << endl;
}

result_t EPFLSolver::solve(const GraphHypothesis& realization,
    vector<pair<double, int>>& clusterSortedScores)
{
  // Get infection times for all the observers and keep ascending.
  vector<pair<double, int>> observers;
  for (int observer : m_observerNodes) {
    // If chosen observers are not infected, we just skip them.
    if (realization.getInfectionTime(observer) == INFECTED_FALSE)
      continue;

    observers.push_back(
        pair<double, int>(realization.getInfectionTime(observer), observer));
  }

  if (observers.size() == 0) {
    result_t empty;

    empty.first = 0;
    empty.second.push_back(0);      // mass
    empty.second.push_back(100);    // diff
    empty.second.push_back(0);      // correct mass

    for (int s = 0; s < m_config.nodes; ++s)
      clusterSortedScores.push_back(pair<double, int>(0.0, s));

    return empty;
  }

  // Everything is related to the reference observer (first infected).
  sort(observers.begin(), observers.end());
  int referenceObserver = observers[0].second;

  // Compute gaussian moments.
  gaussian_t moments;
  if (m_config.infType == BETA) {
    double beta = m_config.cluster.beta;
    moments = make_pair(1.00/beta, sqrt((1.00 - beta)/(beta*beta)));
    cout << moments.first << " -- " << moments.second << endl;
  } else {
    moments = make_pair(m_config.cluster.miu, m_config.cluster.beta);
    // moments = computeGaussianMoments(observers);
  }

  TIntH idToShortestPathsFromReference;
  TSnap::GetShortPath(m_network, referenceObserver, idToShortestPathsFromReference);

  // BFS tree from reference observer. Compute LCA w.r.t. to ref. observer.
  PNGraph bfsTree = TSnap::GetBfsTree(m_network, referenceObserver, true, true);
  int lca[observers.size()][observers.size()];
  for (size_t k = 1; k < observers.size(); ++k) {
    for (size_t i = k + 1; i < observers.size(); ++i) {
      int Ni = observers[i].second;
      int Nk = observers[k].second;
      int heightNi = idToShortestPathsFromReference.GetDat(Ni);
      int heightNk = idToShortestPathsFromReference.GetDat(Nk);

      while (heightNi > heightNk && Ni != referenceObserver) {
        // Note: fingers crossed that the same values are in GetInNId array.
        Ni = bfsTree->GetNI(Ni).GetInNId(0);
        heightNi = idToShortestPathsFromReference.GetDat(Ni);
      }
      while (heightNk > heightNi) {
        // Note: fingers crossed that the same values are in GetInNId array.
        Nk = bfsTree->GetNI(Nk).GetInNId(0);
        heightNk = idToShortestPathsFromReference.GetDat(Nk);
      }
      while (Ni != Nk) {
        Nk = bfsTree->GetNI(Nk).GetInNId(0);
        Ni = bfsTree->GetNI(Ni).GetInNId(0);
      }
      if (heightNk != heightNi || Ni != Nk)
        cout << "Something awfully wrong here" << endl;
      lca[k][i] = lca[i][k] = Ni;
    }
  }

  // Delay Covariance
  Eigen::MatrixXf lambda(observers.size() - 1, observers.size() - 1);
  for (size_t k = 0; k < observers.size() - 1; ++k) {
    for (size_t i = 0; i < observers.size() - 1; ++i) {
      float value = moments.second * moments.second;
      if (k == i)
        value *= idToShortestPathsFromReference.GetDat(observers[k+1].second);
      else
        value *= idToShortestPathsFromReference.GetDat(lca[k+1][i+1]);
      lambda(k, i) = value;
    }
  }

  Eigen::MatrixXf invLambda(observers.size() - 1, observers.size() - 1);
  invLambda = lambda.inverse();

  // Observed Delay
  Eigen::VectorXf d(observers.size() - 1);
  for (size_t o = 0; o < observers.size() - 1; ++o)
    d[o] = observers[o+1].first - observers[0].first;

  for (int s = 0; s < m_network->GetNodes(); ++s) {
    TIntH idToShortestPathsFromSource;
    TSnap::GetShortPath(m_network, s, idToShortestPathsFromSource);

    // Deterministic Delay
    Eigen::VectorXf miu_s(observers.size() - 1);
    for (size_t o = 0; o < observers.size() - 1; ++o) {
      miu_s[o] = moments.first *
        (idToShortestPathsFromSource.GetDat(observers[o+1].second) -
         idToShortestPathsFromSource.GetDat(referenceObserver));
    }

    // Estimator value.
    double estimator = miu_s.transpose() * invLambda * (d - 0.5 * miu_s);
    clusterSortedScores.push_back(pair<double, int>(estimator, s));
  }
  sort(clusterSortedScores.begin(), clusterSortedScores.end(),
       std::greater<pair<double, int>>());

  int realSourceIdx = 0;
  for (size_t i = 0; i < clusterSortedScores.size(); ++i) {
    if (clusterSortedScores[i].second == realization.getSource())
      realSourceIdx = i;
  }

  result_t result;
  result.first = clusterSortedScores[0].second;          // identified solution
  result.second.push_back(clusterSortedScores[0].first); // solution score?
  result.second.push_back(clusterSortedScores[0].first -
      clusterSortedScores[realSourceIdx].first);
  result.second.push_back(realSourceIdx);                // rank

  return result;
}

gaussian_t EPFLSolver::computeGaussianMoments(
    const vector<pair<double, int>>& observers)
{
  gaussian_t moments(0.0, 0.0);
  vector<double> data;

  // Mean
  for (size_t i = 1; i < observers.size(); ++i) {
    data.push_back(observers[i].first - observers[0].first);
    moments.first += data[i-1];
  }
  moments.first /= data.size();

  // Standard Deviation
  for (size_t i = 0; i < data.size(); ++i)
    moments.second += (data[i] - moments.first) * (data[i] - moments.first);
  moments.second = sqrt(moments.second / (data.size() - 1));

  // cout << "Gaussian(" << moments.first << ", " << moments.second << ")\n";
  return moments;
}

