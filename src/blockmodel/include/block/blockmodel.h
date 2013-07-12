/* vim:set ts=4 sw=4 sts=4 et: */

#ifndef BLOCKMODEL_BLOCKMODEL_H
#define BLOCKMODEL_BLOCKMODEL_H

#include <cmath>
#include <stdexcept>
#include <igraph/cpp/graph.h>
#include <igraph/cpp/matrix.h>
#include <igraph/cpp/vector.h>
#include <mtwister/mt.h>

class Blockmodel;

/// Simple struct representing a point mutation of a blockmodel
/**
 * A point mutation is a step that moves vertex i from group k to group l
 */
class PointMutation {
	friend class Blockmodel;

public:
	/// The index of the vertex being moved
	int vertex;
	/// The group in which the vertex is before the move
	int from;
	/// The group in which the vertex is after the move
	int to;

	/// Constructor
	PointMutation(int vertex, int from, int to) :
		vertex(vertex), from(from), to(to) {}

    /// Reverses the mutation
    void reverse();

    /// Returns a reversed copy of the mutation
    PointMutation reversed() const;

private:
    /// Performs the mutation on the given model
    void perform(Blockmodel& model) const;

    /// Undoes the mutation on the given model
    void undo(Blockmodel& model) const;
};

/// Abstract base class for various types of blockmodels
class Blockmodel {
protected:
    /// Pointer to a graph to which this model will be fitted
    igraph::Graph* m_pGraph;

    /// The number of types in the model
    int m_numTypes;

    /// Vector storing the vertex types
    igraph::Vector m_types;

    /// Vector storing the number of vertices of a given type
    igraph::Vector m_typeCounts;

    /// Matrix storing the number of edges between pairs of vertex types
    /**
     * Due to some optimizations, the diagonal of the matrix actually stores
     * twice the number of edges going between vertices of the same type.
     */
    igraph::Matrix m_edgeCounts;
    
    /// Cached value of the log-likelihood
    mutable double m_logLikelihood;

public:
    /// Constructs a new blockmodel not associated with any given graph
    explicit Blockmodel()
        : m_pGraph(0), m_numTypes(0), m_types(),
		  m_typeCounts(), m_edgeCounts(), m_logLikelihood(1) {
    }
	
	/// Copies a blockmodel to another one
	/**
	 * This method is implemented explicitly instead of by the means of a
	 * copy constructor or assignment operator to avoid the temptation of
	 * assigning dereferenced \ref Blockmodel pointers to each other, which
	 * would copy the members of the \em base class (i.e. \ref Blockmodel)
	 * but not necessarily the members of the actual \ref Blockmodel
	 * subclass.
	 */
	virtual void assignFrom(const Blockmodel* other) {
		*this = *other;
	}

	/// Copies a blockmodel to another one (auto_ptr variant)
	/**
	 * \see assignFrom(const Blockmodel*)
	 */
	void assignFrom(const std::auto_ptr<Blockmodel>& other) {
		assignFrom(other.get());
	}

    /// Convenience function to create a blockmodel instance
    /**
     * This function could not have been made a constructor as it calls virtual
     * functions, which is forbidden in constructors.
     */
    template <typename ModelType>
    static ModelType create(igraph::Graph* pGraph, int numTypes) {
        ModelType result;
        result.setGraph(pGraph);
        result.setNumTypes(numTypes);
        return result;
    }

    /// Generates a new graph according to the current parameters of the blockmodel
    igraph::Graph generate() const {
        MersenneTwister rng;
        return generate(rng);
    }

    /// Generates a new graph according to the current parameters of the blockmodel
    virtual igraph::Graph generate(MersenneTwister& rng) const = 0;

    /// Returns the number of edges between the two given groups
    long getEdgeCount(long ri, long ci) const {
        return m_edgeCounts(ri, ci);
    }

	/**
	 * Returns the actual number of edges after a point mutation
	 * between the two affected groups and others.
     *
     * Edges within the same group are counted twice.
	 *
	 * \param       mutation    the point mutation we are planning to perform
	 * \param[out]  countsFrom  the actual number of edges between group
	 *                          \c mutation.from and other groups after the move
	 * \param[out]  countsTo    the actual number of edges between group
	 *                          \c mutation.to and other groups after the move
	 */
	void getEdgeCountsFromAffectedGroupsAfter(const PointMutation& mutation,
			igraph::Vector& countsFrom, igraph::Vector& countsTo) const;

    /// Returns the whole edge count matrix
    igraph::Matrix getEdgeCounts() const {
        return m_edgeCounts;
    }

	/// Returns the probability of the given edge in the model
	virtual double getEdgeProbability(int v1, int v2) = 0;

    /// Returns a pointer to the graph associated to the model (const)
    const igraph::Graph* getGraph() const {
        return m_pGraph;
    }

    /// Returns the log-likelihood of the model
    double getLogLikelihood() const {
        if (m_logLikelihood < 0)
            return m_logLikelihood;
        return recalculateLogLikelihood();
    }

    /// Returns the increase in the log-likelihood of the model after a point mutation
    /**
     * The default is a dumb implementation which actually performs the move and
     * calculates the log-likelihood difference. Subclasses may override this
     * method if a more efficient solution can be provided.
     */
    virtual double getLogLikelihoodIncrease(const PointMutation& mutation);

    /// Returns the number of observations in this model
    long getNumObservations() const {
        return (m_pGraph->vcount() * (m_pGraph->vcount()-1) / 2);
    }

    /// Returns the number of free parameters in this model
    virtual int getNumParameters() const = 0;

    /// Returns the number of types in this model
    int getNumTypes() const {
        return m_numTypes;
    }

    /// Counts how many edges could there be (theoretically) from a group to others
    void getTotalEdgesFromGroup(int type, igraph::Vector& result) const;

    /// Counts how many edges could there be (theoretically) between the two groups
    long int getTotalEdgesBetweenGroups(int type1, int type2) const;

	/**
	 * Returns the theoretically possible number of edges after a point mutation
	 * between the two affected groups and the others.
	 *
     * Edges within the same group are counted twice.
	 *
	 * \param       mutation    the point mutation we are planning to perform
	 * \param[out]  countsFrom  the theoretically possible number of edges between
	 *                          group \c mutation.from and other groups after the move
	 * \param[out]  countsTo    the theoretically possible number of edges between
	 *                          group \c mutation.to and other groups after the move
	 */
	void getTotalEdgesFromAffectedGroupsAfter(const PointMutation& mutation,
			igraph::Vector& countsFrom, igraph::Vector& countsTo) const;

    /// Returns the type of the given vertex
    int getType(long index) const {
        return m_types[index];
    }

    /// Returns how many vertices have the given type
    long getTypeCount(long index) const {
        return m_typeCounts[index];
    }

    /// Returns the whole type count vector
    igraph::Vector getTypeCounts() const {
        return m_typeCounts;
    }

    /// Returns the whole type vector
    igraph::Vector getTypes() const {
        return m_types;
    }
    
    /// Returns the number of vertices in the model
    size_t getVertexCount() const {
        return m_types.size();
    }

    /// Performs the given mutation on the model
    virtual void performMutation(const PointMutation& mutation) {
		mutation.perform(*this);
	}

	/// Randomizes the current configuration of the model
	virtual void randomize(MersenneTwister& rng);

    /// Returns the log-likelihood of the model (with forced recalculation)
    virtual double recalculateLogLikelihood() const = 0;

    /// Sets the graph associated to the model
    /**
     * If the graph is not NULL, the type vector will be resized to the number
     * of vertices in the graph and the edge counts will be re-calculated.
     */
    virtual void setGraph(igraph::Graph* graph);

    /// Sets the number of types
    /**
     * This method should be called only after construction as it will
     * re-create the edge count matrix and the type count vector
     */
    virtual void setNumTypes(int numTypes);

    /// Sets the type of a single vertex
    virtual void setType(long index, int newType);

    /// Sets the types of multiple vertices
    void setTypes(const igraph::Vector& types);

protected:
    /// Recounts the edges and updates m_typeCounts and m_edgeCounts
    virtual void recountEdges();

    /// Invalidates the cached log-likelihood value
    void invalidateCache() {
        m_logLikelihood = 1.0;
    }
};

/// Class representing an undirected blockmodel
class UndirectedBlockmodel : public Blockmodel {
private:
    /// Matrix storing the edge probabilities between pairs of vertex types
    /**
     * This is used only if there is no graph associated to the model (i.e.
     * m_pGraph is NULL).
     */
    igraph::Matrix m_probabilities;

public:
    /// Constructs a new undirected blockmodel not associated with any given graph
    explicit UndirectedBlockmodel() : Blockmodel() {}

	virtual void assignFrom(const Blockmodel* other) {
		*this = dynamic_cast<const UndirectedBlockmodel&>(*other);
	}

    /// Generates a new graph according to the current parameters of the blockmodel
    virtual igraph::Graph generate(MersenneTwister& rng) const;

	/// Returns the probability of the given edge in the model
	virtual double getEdgeProbability(int v1, int v2) {
		return getProbability(m_types[v1], m_types[v2]);
	}

    /// Returns the increase in the log-likelihood of the model after a point mutation
    virtual double getLogLikelihoodIncrease(const PointMutation& mutation);

    /// Returns the number of free parameters in this model
	virtual int getNumParameters() const {
        return (m_numTypes * (m_numTypes+1) / 2.) + m_types.size() + 1;
    }

    /// Returns the estimated probability between the two given types
    double getProbability(int type1, int type2) const;

    /// Returns the estimated probability matrix
    igraph::Matrix getProbabilities() const;

    /// Returns the estimated probability matrix
    void getProbabilities(igraph::Matrix& result) const;

    /// Returns one row of the estimated probability matrix
    igraph::Vector getProbabilitiesFromGroup(int type) const;

    /// Returns one row of the estimated probability matrix
    void getProbabilitiesFromGroup(int type, igraph::Vector& result) const;

    /// Returns the log-likelihood of the model (with forced recalculation)
    virtual double recalculateLogLikelihood() const;

    /// Sets the number of types
    /**
     * This method should be called only after construction as it will
     * re-create the edge count matrix and the type count vector
     */
    virtual void setNumTypes(int numTypes);

    /// Sets the probability of an edge between the two given groups
    /**
     * This method works only if there is no graph associated to the model
     * (i.e. m_pGraph is NULL).
     */
    void setProbability(int type1, int type2, double p) {
        m_probabilities(type1, type2) = m_probabilities(type2, type1) = p;
    }

    /// Sets the probability of edges between vertex groups.
    /**
     * This method works only if there is no graph associated to the model
     * (i.e. m_pGraph is NULL).
     */
    void setProbabilities(const igraph::Matrix& p);
};

/// Class representing an undirected degree-corrected blockmodel
class DegreeCorrectedUndirectedBlockmodel : public Blockmodel {
private:
    /// Degree vector for the graph
	/**
	 * This vector stores the actual degrees of the vertices in the graph if
	 * the model is associated to a graph.
	 */
    igraph::Vector m_degrees;

	/// Drift counter for the iterative re-calculation of log-likelihood
	double m_driftCounter;

	/// Stickiness vector for the vertices
	/**
	 * For models without graphs, this vector must be set to the stickiness values
	 * of each vertex. For models associated to a graph, this vector is unused and
	 * empty.
	 */
	igraph::Vector m_stickinesses;

    /// Sum of degrees for vertices in a given group
    igraph::Vector m_sumOfDegreesByType;

public:
    /**
     * \brief Constructs a new undirected degree-corrected blockmodel not
     *        associated with any given graph
     */
    explicit DegreeCorrectedUndirectedBlockmodel()
        : Blockmodel(), m_degrees(), m_driftCounter(0), m_stickinesses(),
		m_sumOfDegreesByType() {}

	virtual void assignFrom(const Blockmodel* other) {
		*this = dynamic_cast<const DegreeCorrectedUndirectedBlockmodel&>(*other);
	}

    /// Generates a new graph according to the current parameters of the blockmodel
    virtual igraph::Graph generate(MersenneTwister& rng) const;

	/// Returns the probability of the given edge in the model
	virtual double getEdgeProbability(int v1, int v2) {
		int type1 = m_types[v1], type2 = m_types[v2];
		double lambda = m_edgeCounts(type1, type2);

		if (m_pGraph == NULL) {
			/* No graph associated */
			lambda *= m_stickinesses[v1] * m_stickinesses[v2];
		} else {
			/* We have a graph, so we use m_degrees */
			lambda *= m_degrees[v1] / m_sumOfDegreesByType[type1];
			lambda *= m_degrees[v2] / m_sumOfDegreesByType[type2];
		}

		return 1 - std::exp(-lambda);
	}

    /// Returns the increase in the log-likelihood of the model after a point mutation
    virtual double getLogLikelihoodIncrease(const PointMutation& mutation);

    /// Returns the number of free parameters in this model
	virtual int getNumParameters() const {
        return (m_numTypes * (m_numTypes+1) / 2.) + 2 * m_types.size() + 1;
    }

    /// Returns the rate parameter between the two given types in this model
    double getRate(int type1, int type2) const;

    /// Returns the matrix of rate parameters for this model
    igraph::Matrix getRates() const;

    /// Returns the matrix of rate parameters for this model
    void getRates(igraph::Matrix& result) const;

    /// Returns the stickiness values of all the vertices in this model
    igraph::Vector getStickinesses() const;

    /// Returns the stickiness values of all the vertices in this model
    void getStickinesses(igraph::Vector& result) const;

    /// Performs the given mutation on the model
	/**
	 * This method is overridden from the parent because it is more efficient
	 * to calculate the new log-likelihood using \ref getLogLikelihoodIncrease
	 * instead of re-calculating it completely.
	 *
	 * However, after every 8192 steps, we re-calculate the log-likelihood
	 * completely to avoid the accumulation of numerical errors.
	 */
    virtual void performMutation(const PointMutation& mutation) {
		if (m_driftCounter < 8192) {
			double oldLogLikelihood = m_logLikelihood;
			double increase = getLogLikelihoodIncrease(mutation);
			Blockmodel::performMutation(mutation);
			m_logLikelihood = oldLogLikelihood + increase;
			m_driftCounter++;
		} else {
			Blockmodel::performMutation(mutation);
			m_driftCounter = 0;
		}
	}

    /// Returns the log-likelihood of the model (with forced recalculation)
    virtual double recalculateLogLikelihood() const;

    /// Recounts the edges and updates m_typeCounts, m_edgeCounts and m_sumOfDegreesByType
    virtual void recountEdges();

    /// Sets the graph associated to the model
    /**
     * If the graph is not NULL, the type vector will be resized to the number
     * of vertices in the graph and the edge counts will be re-calculated.
     * The cached degree vector will also be re-calculated.
     */
    virtual void setGraph(igraph::Graph* graph);

	/// Sets the stickiness values of vertices in the model
	/**
	 * This method should be called only if the model is not associated to
	 * a graph.
	 */
	void setStickinesses(const igraph::Vector& stickinesses);

    /// Sets the number of types
    /**
     * This method should be called only after construction as it will
     * re-create the edge count matrix and the type count vector
     */
    virtual void setNumTypes(int numTypes);

    /// Sets the Poisson rates of vertex group pairs.
    /**
     * This method works only if there is no graph associated to the model
     * (i.e. m_pGraph is NULL).
     */
    void setRates(const igraph::Matrix& r);

    /// Sets the type of a single vertex
    virtual void setType(long index, int newType);
};
#endif
