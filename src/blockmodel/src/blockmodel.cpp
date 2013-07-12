/* vim:set ts=4 sw=4 sts=4 et: */

#include <algorithm>
#include <cassert>
#include <cmath>
#include <stdexcept>
#include <block/blockmodel.h>
#include <igraph/cpp/vertex_selector.h>

using namespace igraph;

namespace {
    double binary_entropy(double prob) {
        if (prob <= 0 || prob >= 1)
            return 0.0;
        return prob * std::log(prob) + (1 - prob) * std::log(1 - prob);
    }
}


/***************************************************************************/

void PointMutation::perform(Blockmodel& model) const {
    assert(model.getType(vertex) == from);
    model.setType(vertex, to);
}

void PointMutation::reverse() {
    std::swap(from, to);
}

PointMutation PointMutation::reversed() const {
    return PointMutation(vertex, to, from);
}

void PointMutation::undo(Blockmodel& model) const {
    assert(model.getType(vertex) == to);
    model.setType(vertex, from);
}

/***************************************************************************/

void Blockmodel::getEdgeCountsFromAffectedGroupsAfter(
        const PointMutation& mutation,
        igraph::Vector& countsFrom, igraph::Vector& countsTo) const {
    m_edgeCounts.getRow(mutation.from, countsFrom);
    m_edgeCounts.getCol(mutation.to, countsTo);

    if (mutation.from == mutation.to)
        return;

    // Get the neighbors of the affected vertex
    Vector neighbors = m_pGraph->neighbors(mutation.vertex);

    // Adjust countsFrom and countsTo accordingly
    for (Vector::const_iterator it = neighbors.begin();
         it != neighbors.end(); it++) {
        long otherType = m_types[*it];
        countsFrom[otherType]--;
        countsTo[otherType]++;
        if (otherType == mutation.from) {
            countsFrom[otherType]--;
            countsFrom[mutation.to]++;
        } else if (otherType == mutation.to) {
            countsTo[otherType]++;
            countsTo[mutation.from]--;
        }
    }
}

double Blockmodel::getLogLikelihoodIncrease(
        const PointMutation& mutation) {
    double result = -getLogLikelihood();
    mutation.perform(*this);
    result += getLogLikelihood();
    mutation.undo(*this);
    return result;
}

long int Blockmodel::getTotalEdgesBetweenGroups(int type1, int type2) const {
    if (type1 == type2)
        return (m_typeCounts[type1] - 1) * m_typeCounts[type1];
    return m_typeCounts[type1] * m_typeCounts[type2];
}

void Blockmodel::getTotalEdgesFromGroup(int type, Vector& result) const {
    result = m_typeCounts[type] * m_typeCounts;
    result[type] -= m_typeCounts[type];
}

void Blockmodel::getTotalEdgesFromAffectedGroupsAfter(
        const PointMutation& mutation,
        Vector& countsFrom, Vector& countsTo) const {
    countsFrom = m_typeCounts[mutation.from] * m_typeCounts;
    countsFrom[mutation.from] -= m_typeCounts[mutation.from];

    if (mutation.from == mutation.to) {
        countsTo = countsFrom;
        return;
    }

    countsTo = m_typeCounts[mutation.to] * m_typeCounts;
    countsTo[mutation.to] -= m_typeCounts[mutation.to];

    // We take into account the point mutation from here
    countsFrom -= m_typeCounts;
    countsFrom[mutation.from] -= m_typeCounts[mutation.from];
    countsTo += m_typeCounts;
    countsTo[mutation.to] += m_typeCounts[mutation.to];

    countsFrom[mutation.from] += 2;
    countsFrom[mutation.to] += m_typeCounts[mutation.from]-1;
    countsTo[mutation.from] -= m_typeCounts[mutation.to]+1;
}

void Blockmodel::randomize(MersenneTwister& rng) {
    for (Vector::iterator it = m_types.begin(); it != m_types.end(); it++)
        *it = rng.randint(m_numTypes);
    // Invalidate the log-likelihood cache
    invalidateCache();
    // Recount the edges
    recountEdges();
}

void Blockmodel::recountEdges() {
    if (m_pGraph == NULL)
        return;

    long n = m_pGraph->vcount();
    Vector edgelist = m_pGraph->getEdgelist();
    Vector::const_iterator it = edgelist.begin();

    m_typeCounts.fill(0);
    for (long i = 0; i < n; i++) {
        m_typeCounts[m_types[i]]++;
    }

    m_edgeCounts.fill(0);
    while (it != edgelist.end()) {
        long type1 = m_types[*(it++)];
        long type2 = m_types[*(it++)];
        m_edgeCounts(type1,type2) += 1;
        m_edgeCounts(type2,type1) += 1;
    }
}

void Blockmodel::setGraph(igraph::Graph* graph) {
    size_t oldSize = m_types.size();

    m_pGraph = graph;
    if (m_pGraph != NULL) {
        m_types.resize(m_pGraph->vcount());
        for (; oldSize < m_types.size(); oldSize++)
            m_types[oldSize] = 0;
        recountEdges();
    }
}

void Blockmodel::setNumTypes(int numTypes) {
    if (numTypes <= 0)
        throw std::runtime_error("must have at least one type");

    m_numTypes = numTypes;
    m_typeCounts.resize(numTypes);
    if (m_pGraph != NULL)
        m_typeCounts[0] += m_pGraph->vcount() - m_typeCounts.sum();

    m_edgeCounts.resize(numTypes, numTypes);
    recountEdges();
}

void Blockmodel::setType(long index, int newType) {
    // Save the old type
    int oldType = m_types[index];
    if (oldType == newType)
        return;
    // Get the neighbors of the affected vertex
    Vector neighbors = m_pGraph->neighbors(index);
    // Adjust the edge counts and the type counts
    // Here we assume that there are no loop edges
    m_typeCounts[oldType]--; m_typeCounts[newType]++;
    for (Vector::const_iterator it = neighbors.begin();
         it != neighbors.end(); it++) {
        long otherType = m_types[*it];
        m_edgeCounts(oldType, otherType)--;
        m_edgeCounts(otherType, oldType)--;
        m_edgeCounts(newType, otherType)++;
        m_edgeCounts(otherType, newType)++;
    }
    // Set the type of the vertex to the new type
    m_types[index] = newType;
    // Invalidate the log-likelihood cache
    invalidateCache();
}

void Blockmodel::setTypes(const Vector& types) {
    if (types.min() < 0)
        throw std::runtime_error("negative type index found in type vector");
    if (types.max() >= m_numTypes)
        throw std::runtime_error("too large type index found in type vector");

    m_types = types;
    // Invalidate the log-likelihood cache
    invalidateCache();
    // Recount the edges
    recountEdges();
}

/***************************************************************************/

Graph UndirectedBlockmodel::generate(MersenneTwister& rng) const {
    long int n = m_types.size();
    const Matrix probs = getProbabilities();
    Graph graph(n);
    Vector edges;

    for (long int v1 = 0; v1 < n; v1++) {
        int type1 = m_types[v1];
        for (long int v2 = v1+1; v2 < n; v2++) {
            if (rng.random() < probs(type1, m_types[v2])) {
                edges.push_back(v1);
                edges.push_back(v2);
            }
        }
    }

    graph.addEdges(edges);
    return graph;
}

double UndirectedBlockmodel::recalculateLogLikelihood() const {
    double den, result = 0.0;

    for (int i = 0; i < m_numTypes; i++) {
        if (m_typeCounts[i] == 0)
            continue;

        for (int j = i+1; j < m_numTypes; j++) {
            den = getTotalEdgesBetweenGroups(i, j);
            if (den == 0)
                continue;
            result += den * binary_entropy(m_edgeCounts(i, j) / den);
        }

        den = getTotalEdgesBetweenGroups(i, i);
        if (den == 0)
            continue;
        result += (den/2) * binary_entropy(m_edgeCounts(i, i) / den);
    }

    m_logLikelihood = result;
    return result;
}

double UndirectedBlockmodel::getLogLikelihoodIncrease(
        const PointMutation& mutation) {
    double result = 0.0;

    if (mutation.from == mutation.to)
        return result;

    /* Get the old probabilities and counts */

    Vector oldProbsFrom(m_numTypes),  oldProbsTo(m_numTypes);
    Vector oldCountsFrom(m_numTypes), oldCountsTo(m_numTypes);

    getTotalEdgesFromGroup(mutation.from, oldCountsFrom);
    getTotalEdgesFromGroup(mutation.to,   oldCountsTo);
    getProbabilitiesFromGroup(mutation.from, oldProbsFrom);
    getProbabilitiesFromGroup(mutation.to,   oldProbsTo);

    /* Get the new probabilities and counts */

    Vector newCountsFrom(m_numTypes), newCountsTo(m_numTypes);
    Vector newProbsFrom(m_numTypes),  newProbsTo(m_numTypes);

    getTotalEdgesFromAffectedGroupsAfter(mutation, newCountsFrom, newCountsTo);
    getEdgeCountsFromAffectedGroupsAfter(mutation, newProbsFrom,  newProbsTo);
    for (int i = 0; i < m_numTypes; i++) {
        if (newCountsFrom[i] > 0)
            newProbsFrom[i] /= newCountsFrom[i];
        if (newCountsTo[i] > 0)
            newProbsTo[i] /= newCountsTo[i];
    }

    /* Correct the diagonal elements in oldCounts* and newCounts* as they
     * contain twice the number of edges there */
    oldCountsFrom[mutation.from] /= 2;
    oldCountsTo[mutation.to] /= 2;
    newCountsFrom[mutation.from] /= 2;
    newCountsTo[mutation.to] /= 2;

    /* Calculate the increase */
    for (int i = 0; i < m_numTypes; i++) {
        result += newCountsFrom[i] * binary_entropy(newProbsFrom[i]);
        result -= oldCountsFrom[i] * binary_entropy(oldProbsFrom[i]);
        result += newCountsTo  [i] * binary_entropy(newProbsTo  [i]);
        result -= oldCountsTo  [i] * binary_entropy(oldProbsTo  [i]);
    }
    result -= newCountsFrom[mutation.to] *
              binary_entropy(newProbsFrom[mutation.to]);
    result += oldCountsFrom[mutation.to] *
              binary_entropy(oldProbsFrom[mutation.to]);

    return result;
}

double UndirectedBlockmodel::getProbability(int type1, int type2) const {
    if (m_pGraph == NULL)
        return m_probabilities(type1, type2);

    double possibleEdges = getTotalEdgesBetweenGroups(type1, type2);

    if (possibleEdges == 0)
        return 0;

    return m_edgeCounts(type1, type2) / possibleEdges;
}

Matrix UndirectedBlockmodel::getProbabilities() const {
    Matrix result(m_numTypes, m_numTypes);
    getProbabilities(result);
    return result;
}

void UndirectedBlockmodel::getProbabilities(Matrix& result) const {
    if (m_pGraph == NULL) {
        result = m_probabilities;
        return;
    }

    result.resize(m_numTypes, m_numTypes);

    for (int i = 0; i < m_numTypes; i++) {
        double den;

        if (m_typeCounts[i] == 0) {
            /* numerator is strongly zero */
            for (int j = i; j < m_numTypes; j++) {
                result(i, j) = result(j, i) = 0;
            }
            continue;
        }

        for (int j = i; j < m_numTypes; j++) {
            den = getTotalEdgesBetweenGroups(i, j);
            result(i, j) = result(j, i) =
                (den == 0) ? 0 : (m_edgeCounts(i, j) / den);
        }
    }
}

Vector UndirectedBlockmodel::getProbabilitiesFromGroup(int type) const {
    Vector result;
    getProbabilitiesFromGroup(type, result);
    return result;
}

void UndirectedBlockmodel::getProbabilitiesFromGroup(int type,
        igraph::Vector& result) const {
    getTotalEdgesFromGroup(type, result);
    for (int i = 0; i < m_numTypes; i++) {
        if (result[i] > 0)
            result[i] = m_edgeCounts(type, i) / result[i];
    }
}

void UndirectedBlockmodel::setNumTypes(int numTypes) {
    Blockmodel::setNumTypes(numTypes);
    if (m_pGraph != NULL)
        m_probabilities.resize(0, 0);
    else
        m_probabilities.resize(numTypes, numTypes);
}

void UndirectedBlockmodel::setProbabilities(const Matrix& p) {
    if (!p.isSymmetric())
        throw new std::runtime_error("probability matrix must be symmetric");
    m_probabilities = p;
}

/***************************************************************************/

namespace {
    /* Auxiliary function to generate Poisson-distributed random numbers */
    double poisson(MersenneTwister& rng, double lambda) {
        double p0 = std::exp(-lambda), p = 1;
        int k = 0;
        do {
            k++;
            p *= rng.random();
        } while (p > p0);
        return k-1;
    }
}

Graph DegreeCorrectedUndirectedBlockmodel::generate(
        MersenneTwister& rng) const {
    long int n = m_types.size();
    const Matrix rates = getRates();
    const Vector stickinesses = getStickinesses();

    Graph graph(n);
    Vector edges;

    for (long int v1 = 0; v1 < n; v1++) {
        double s1 = stickinesses[v1];
        int type1 = m_types[v1];
        for (long int v2 = v1+1; v2 < n; v2++) {
            double s2 = stickinesses[v2];
            double rate = s1 * s2 * rates(type1, m_types[v2]);
            int numEdges = poisson(rng, rate);
            for (int i = 0; i < numEdges; i++) {
                edges.push_back(v1);
                edges.push_back(v2);
            }
        }
    }

    graph.addEdges(edges);
    return graph;
}

/* Auxiliary functions for getLogLikelihoodIncrease */
namespace {
	inline double b(double x) {
		return x == 0 ? 0 : x*(std::log(x) - 1);
	}
}

double DegreeCorrectedUndirectedBlockmodel::getLogLikelihoodIncrease(
        const PointMutation& mutation) {
	int r = mutation.from, s = mutation.to;

	if (r == s)
		return 0.0;

	double result = 0.0, degree;
	Vector k(m_numTypes);

    // Get the neighbors of the affected vertex
    Vector neighbors = m_pGraph->neighbors(mutation.vertex);

	// Calculate k, the number of edges between the vertex and other
	// vertices of a given type
    for (Vector::const_iterator it = neighbors.begin();
         it != neighbors.end(); it++)
        k[m_types[*it]]++;
	degree = neighbors.size();

	// Calculate the difference
	for (int t = 0; t < m_numTypes; t++) {
		double diff = 0.0;
		if (t == r) {
			if (k[r] > 0) {
				diff += b(m_edgeCounts(r, r) - 2 * k[r]);
				diff -= b(m_edgeCounts(r, r));
			}
			if (k[r] != k[s]) {
				diff += b(m_edgeCounts(r, s) + k[r] - k[s]);
				diff -= b(m_edgeCounts(r, s));
			}
			diff /= 2;
		} else if (t == s) {
			if (k[r] != k[s]) {
				diff += b(m_edgeCounts(r, s) + k[r] - k[s]);
				diff -= b(m_edgeCounts(r, s));
			}
			if (k[s] > 0) {
				diff += b(m_edgeCounts(s, s) + 2 * k[s]);
				diff -= b(m_edgeCounts(s, s));
			}
			diff /= 2;
		} else if (k[t] > 0) {
			diff += b(m_edgeCounts(r, t) - k[t]);
			diff -= b(m_edgeCounts(r, t));
			diff += b(m_edgeCounts(s, t) + k[t]);
			diff -= b(m_edgeCounts(s, t));
		}
		result += diff;
	}

	// Correction for the m_degrees * logTheta term.
	// Affected are all the nodes in groups r and s
	result += b(m_sumOfDegreesByType[r]) - b(m_sumOfDegreesByType[r] - degree);
	result += b(m_sumOfDegreesByType[s]) - b(m_sumOfDegreesByType[s] + degree);
	
	return result;
}

double DegreeCorrectedUndirectedBlockmodel::recalculateLogLikelihood() const {
    Vector logTheta = m_degrees;

    for (size_t i = 0; i < logTheta.size(); i++) {
        logTheta[i] /= m_sumOfDegreesByType[m_types[i]];
        logTheta[i] = std::log(logTheta[i]);
    }

    double result = m_degrees * logTheta;

    for (int i = 0; i < m_numTypes; i++) {
        double ec;

        ec = m_edgeCounts(i, i);
        if (ec > 0) {
            result += 0.5 * ec * (std::log(ec) - 1);
		}

        for (int j = i+1; j < m_numTypes; j++) {
            ec = m_edgeCounts(i, j);
            if (ec > 0)
                result += ec * (std::log(ec) - 1);
        }
    }

    m_logLikelihood = result;
    return result;
}

double DegreeCorrectedUndirectedBlockmodel::getRate(int i, int j) const {
    return m_edgeCounts(i, j);
}

void DegreeCorrectedUndirectedBlockmodel::getRates(Matrix& result) const {
    result = m_edgeCounts;
}

Matrix DegreeCorrectedUndirectedBlockmodel::getRates() const {
    Matrix result;
    getRates(result);
    return result;
}

void DegreeCorrectedUndirectedBlockmodel::getStickinesses(Vector& result) const {
    if (m_pGraph) {
        /* If we have a graph, estimate stickinesses from degrees */
        m_pGraph->degree(&result, V(m_pGraph));
        for (size_t i = 0; i < result.size(); i++)
            result[i] /= m_sumOfDegreesByType[m_types[i]];
    } else {
        /* Use pre-defined stickiness values */
        result = m_stickinesses;
    }
}

Vector DegreeCorrectedUndirectedBlockmodel::getStickinesses() const {
    Vector result;
    getStickinesses(result);
    return result;
}

void DegreeCorrectedUndirectedBlockmodel::recountEdges() {
    if (m_pGraph == NULL)
        return;

    Blockmodel::recountEdges();

    long n = m_pGraph->vcount();

    m_sumOfDegreesByType.fill(0);
    for (long i = 0; i < n; i++) {
        m_sumOfDegreesByType[m_types[i]] += m_degrees[i];
    }
}

void DegreeCorrectedUndirectedBlockmodel::setStickinesses(
        const igraph::Vector& stickinesses) {
    if (m_pGraph != NULL)
        throw std::runtime_error("cannot set stickiness values for a model "
                "with an associated graph");
    m_stickinesses = stickinesses;
}

void DegreeCorrectedUndirectedBlockmodel::setGraph(igraph::Graph* graph) {
    if (graph != NULL)
        graph->degree(&m_degrees, VertexSelector::All());
    else
        m_degrees.clear();

    Blockmodel::setGraph(graph);
}

void DegreeCorrectedUndirectedBlockmodel::setNumTypes(int numTypes) {
    m_sumOfDegreesByType.resize(numTypes);
    Blockmodel::setNumTypes(numTypes);
}

void DegreeCorrectedUndirectedBlockmodel::setRates(const igraph::Matrix& r) {
    if (!r.isSymmetric())
        throw new std::runtime_error("rate matrix must be symmetric");
    if (m_pGraph != NULL)
        throw new std::runtime_error("cannot set rate matrix when a graph is linked "
                                     "to the model");
    m_edgeCounts = r;

    // Invalidate the log-likelihood cache
    invalidateCache();
    // Recount the edges
    recountEdges();
}

void DegreeCorrectedUndirectedBlockmodel::setType(long index, int newType) {
    int oldType = m_types[index];
    if (m_types[index] == newType)
        return;

    double degree = m_degrees[index];
    m_sumOfDegreesByType[oldType] -= degree;
    Blockmodel::setType(index, newType);
    m_sumOfDegreesByType[newType] += degree;
}


