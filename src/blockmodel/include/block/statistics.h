/* vim:set ts=4 sw=4 sts=4 et: */

#ifndef BLOCKMODEL_STATISTICS_H
#define BLOCKMODEL_STATISTICS_H

#include <igraph/cpp/vector.h>

/// Alternative hypothesis types
typedef enum {
    NOT_EQUAL, GREATER_THAN, LESS_THAN
} H1;

/// Returns the area under the normal curve to the left of the given z value
double z_probability(double z);

/// Two-sample Mann-Whitney U test
class MannWhitneyTest {
private:
    /// The alternative hypothesis
    H1 m_alternative;

    /// Size of sample A
    unsigned long m_nA;

    /// Size of sample B
    unsigned long m_nB;

    /// The tie correction value that was applied
    double m_tieCorrection;

    /// The tolerance used in floating point comparisons
    double m_tolerance;

    /// The value of the test statistic
    double m_u;

public:
    /// Constructs a two-sample Mann-Whitney test with no samples
    explicit MannWhitneyTest(H1 alternative=NOT_EQUAL, double tolerance=1e-7) :
        m_alternative(alternative), m_tieCorrection(),
        m_tolerance(tolerance), m_u() {
        igraph::Vector empty;
        setData(empty, empty);
    }

    /// Constructs a two-sample Mann-Whitney test with the given samples
    MannWhitneyTest(const igraph::Vector& xA, const igraph::Vector& xB,
            H1 alternative=NOT_EQUAL, double tolerance=1e-7) :
        m_alternative(alternative), m_tieCorrection(),
        m_tolerance(tolerance), m_u() {
        setData(xA, xB);
    }

    /// Returns the p-value of the test
    double getP() const;

    /// Returns the value of the test statistic
    double getTestStatistic() const {
        return m_u;
    }

    /// Returns the tie correction that was applied to the p-value
    double getTieCorrection() const {
        return m_tieCorrection;
    }

    /// Sets the data used in the two-sample test
    void setData(const igraph::Vector& xA, const igraph::Vector& xB);

    /// Returns the size of the first sample
    unsigned long sizeA() const {
        return m_nA;
    }

    /// Returns the size of the second sample
    unsigned long sizeB() const {
        return m_nB;
    }
};

#endif
