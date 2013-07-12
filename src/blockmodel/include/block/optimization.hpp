/* vim:set ts=4 sw=4 sts=4 et: */

#ifndef BLOCKMODEL_OPTIMIZATION_HPP
#define BLOCKMODEL_OPTIMIZATION_HPP

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <block/blockmodel.h>
#include <block/math.hpp>
#include <igraph/cpp/types.h>
#include <mtwister/mt.h>

/// Abstract optimization strategy class
template <typename Model>
class OptimizationStrategy {
protected:
    /// The number of steps taken
    int m_stepCount;

public:
    /// Constructs an optimization strategy not attached to any model
    OptimizationStrategy() : m_stepCount(0) {}

    /// Returns the number of steps taken so far
    int getStepCount() const {
        return m_stepCount;
    }

    /// Runs the optimization procedure
    void optimize(Model* model) {
        while (step(model));
    }

    /// Runs one step of the optimization strategy
    /**
     * Returns true if the state changed, false otherwise.
     */
    virtual bool step(Model* model) = 0;

    /// Increases the step counter by one
    void stepDone() { m_stepCount++; }
};

/// Greedy optimization strategy for an undirected blockmodel
/**
 * This strategy takes a blockmodel and repeatedly moves vertices between
 * groups in a way that maximizes the local contribution of a vertex to
 * the global likelihood.
 */
template <typename Model>
class GreedyStrategy : public OptimizationStrategy<Model> {
public:
    virtual bool step(Model* pModel) {
        long int n = pModel->getGraph()->vcount();
        int k = pModel->getNumTypes();
        igraph::Vector newTypes(pModel->getTypes());

        for (long int i = 0; i < n; i++) {
            double bestLogL = 0;
            double newLogL;

            for (int j = 0; j < k; j++) {
                PointMutation mutation(i, newTypes[i], j);
                newLogL = pModel->getLogLikelihoodIncrease(mutation);
                if (newLogL > bestLogL) {
                    bestLogL = newLogL;
                    newTypes[i] = j;
                }
            }
        }

        double oldLogL = pModel->getLogLikelihood();
        pModel->setTypes(newTypes);
        return (pModel->getLogLikelihood() > oldLogL);
    }
};

/// Optimization strategy that uses a random number generator
template <typename Model>
class RandomizedOptimizationStrategy : public OptimizationStrategy<Model> {
protected:
    /// The random generator used by the strategy
    std::auto_ptr<MersenneTwister> m_pRng;

public:
    /// Constructor
    RandomizedOptimizationStrategy() : m_pRng(new MersenneTwister()) {}

    /// Returns the random generator used by the strategy
    MersenneTwister* getRNG() {
        return m_pRng.get();
    }

    /// Sets the random generator used by the strategy
    void setRNG(MersenneTwister* pRng) {
        m_pRng.reset(pRng);
    }
};

/// Metropolis-Hastings algorithm for a blockmodel
/**
 * In each step, a vertex is selected randomly and a random new group is
 * proposed. There are two possibilities:
 *
 * 1. The new group improves the log-likelihood of the model. In this case,
 *    the new group is accepted unconditionally.
 *
 * 2. The new group does not improve the log-likelihood of the model. In
 *    this case, the new group is accepted with a probability equal to
 *    the likelihood ratio of the new and the old configuration. If the
 *    new group is rejected, the same sample will be returned.
 */
class MetropolisHastingsStrategy : public RandomizedOptimizationStrategy<Blockmodel> {
private:
    /// The moving average that tracks the acceptance ratio
    MovingAverage<bool> m_acceptanceRatio;

    /// Whether the last proposal was accepted or not
    bool m_lastProposalAccepted;

public:
    /// Constructor
    MetropolisHastingsStrategy() : RandomizedOptimizationStrategy<Blockmodel>(),
        m_acceptanceRatio(1000), m_lastProposalAccepted(false) {
    }

    /// Returns the acceptance ratio
    float getAcceptanceRatio() const {
        return m_acceptanceRatio.value();
    }

    /// Advances the Markov chain by one step
    virtual bool step(Blockmodel* pModel) {
        int i = m_pRng->randint(pModel->getGraph()->vcount());
        int newType = m_pRng->randint(pModel->getNumTypes());
        PointMutation mutation(i, pModel->getType(i), newType);
        double logLDiff = pModel->getLogLikelihoodIncrease(mutation);

        m_lastProposalAccepted =
            (logLDiff >= 0) || (m_pRng->random() <= std::exp(logLDiff));
        if (m_lastProposalAccepted)
            pModel->performMutation(mutation);

        m_acceptanceRatio.push_back(m_lastProposalAccepted);

        stepDone();

        return true;
    }

    /// Returns whether the last proposal was accepted or not
    bool wasLastProposalAccepted() const {
        return m_lastProposalAccepted;
    }
};

/// Gibbs sampling for a blockmodel
/**
 * In each step, a vertex is selected randomly and a new group is
 * proposed according to the log-likelihood conditioned on all but the
 * selected vertex.
 */
class GibbsSamplingStrategy : public RandomizedOptimizationStrategy<Blockmodel> {
public:
    /// Constructor
    GibbsSamplingStrategy() : RandomizedOptimizationStrategy<Blockmodel>() {}

    /// Returns the acceptance ratio
    float getAcceptanceRatio() const { return 1.0; }

    /// Advances the Markov chain by one step
    virtual bool step(Blockmodel* pModel) {
        int i = m_pRng->randint(pModel->getGraph()->vcount());
        int oldType = pModel->getType(i);
        long int k = pModel->getNumTypes();
        igraph::Vector logLs(k);

        // TODO: maybe this can be calculated more efficiently?
        for (int j = 0; j < k; j++) {
            PointMutation mutation(i, oldType, j);
            logLs[j] = pModel->getLogLikelihoodIncrease(mutation);
        }

        // Subtract the minimum log-likelihood from the log-likelihoods
        logLs -= logLs.min();

        // Run exp(), get cumulative sum
        igraph::real_t (*igraph_exp)(igraph::real_t) = &std::exp;
        std::transform(logLs.begin(), logLs.end(), logLs.begin(), igraph_exp);
        std::partial_sum(logLs.begin(), logLs.end(), logLs.begin());

        // Select a new type based on the log-likelihood distribution
        logLs.binsearch(m_pRng->random() * logLs.back(), &k);
        pModel->setType(i, k);

        stepDone();

        return true;
    }

    /// Returns whether the last proposal was accepted or not
    bool wasLastProposalAccepted() const { return true; }
};

#endif
