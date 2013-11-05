/**
 * Data structure used to abstract away a test instance.
 */

#include <utility>

using namespace std;

/**
 * One test is a pair of (score, nodeId). The meaning of the two items is:
 * Score  - the expected reduction in weight from the prior to the posterior
 *          distribution, given that the test runs.
 * Node   - the node tested (and further assumed to be infected)
 *
 * Using the implicit STL pair structure for now, as it easily allows us to
 * store a priority queue with all the tests and pick the one that cuts the
 * most edges in EC2.
 *
 * TODO(vcarbune): Make this an explicit class for readability
 * and extension purposes.
 */
typedef pair<double, int> Test;
