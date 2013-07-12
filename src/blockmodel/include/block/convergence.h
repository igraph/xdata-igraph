/* vim:set ts=4 sw=4 sts=4 et: */

#ifndef BLOCKMODEL_CONVERGENCE_H
#define BLOCKMODEL_CONVERGENCE_H

#include <block/statistics.h>
#include <igraph/cpp/vector.h>

/// Abstract convergence criterion class
/**
 * Convergence criterion classes are used to decide whether a given Markov
 * chain has converged to its stationary distribution or not. There are a
 * number of different approaches for defining convergence; in this library
 * we use approaches that make their decisions based on a block of K samples
 * from the Markov chain. The check() method can be used to check for
 * convergence; the result is true if the chain has converged, false if not.
 *
 * The convergence criterion classes are usually stateful; e.g., they usually
 * compare the current block with the previous one(s). The reset() method
 * can be used to reset the criterion to the initial state.
 *
 * Convergence criterion classes may also return a short, one-line summary of
 * their state if requested via the report() method. By default, the report()
 * method returns an empty string.
 */
class ConvergenceCriterion {
public:
    virtual bool check(const igraph::Vector& samples) = 0;
    virtual std::string report() const { return ""; }
    virtual void reset() {}
};

/// Converge criterion for Markov chains based on entropy estimation
/**
 * This class decides that the Markov chain has converged if the difference
 * between the estimated entropy (i.e. the average log-likelihood) of the
 * last two blocks is less than a given threshold (1.0 by default).
 */
class EntropyConvergenceCriterion : public ConvergenceCriterion {
private:
    /// Estimated entropy of the last block
    double m_lastEntropy;

    /// The entropy difference threshold
    float m_threshold;

public:
    explicit EntropyConvergenceCriterion(float threshold = 1.0)
        : ConvergenceCriterion(), m_threshold(threshold) {
        reset();
    }

    virtual bool check(const igraph::Vector& samples);
    virtual std::string report() const;
    virtual void reset();
};

/// Converge criterion for Markov chains based on a Mann-Whitney U test
/**
 * This class decides that the Markov chain has converged if the last two
 * blocks of log-likelihood values pass the Mann-Whitney U test at a
 * given significance level (say, 95%).
 */
class MannWhitneyConvergenceCriterion : public ConvergenceCriterion {
private:
    /// Log-likelihood values of the last block
    igraph::Vector m_prevSamples;

    /// The p-value threshold
    float m_significance;

    /// A Mann-Whitney test instance
    MannWhitneyTest m_test;

public:
    explicit MannWhitneyConvergenceCriterion(float significance = 0.1)
        : ConvergenceCriterion(), m_prevSamples(),
          m_significance(significance) {
        reset();
    }

    virtual bool check(const igraph::Vector& samples);
    virtual std::string report() const;
    virtual void reset();
};

#endif

