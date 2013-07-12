/* vim:set ts=4 sw=4 sts=4 et: */

#include <algorithm>
#include <cmath>
#include <numeric>
#include <block/statistics.h>
#include <block/util.hpp>

using igraph::Vector;

double z_probability(double z) {
    const double Z_MAX = 6.0;
    double x, y, w;

    if (z == 0.0)
        x = 0.0;
    else {
        y = 0.5 * std::fabs(z);
        if (y >= Z_MAX * 0.5)
            x = 1.0;
        else if (y < 1.0) {
            w = y * y;
            x = ((((((((0.000124818987 * w
                        -0.001075204047) * w +0.005198775019) * w
                      -0.019198292004) * w +0.059054035642) * w
                    -0.151968751364) * w +0.319152932694) * w
                  -0.531923007300) * w +0.797884560593) * y * 2.0;
        } else {
            y = y - 2.0;
            x = (((((((((((((-0.000045255659 * y
                        +0.000152529290) * y -0.000019538132) * y
                      -0.000676904986) * y +0.001390604284) * y
                    -0.000794620820) * y -0.002034254874) * y
                  +0.006549791214) * y -0.010557625006) * y
                +0.011630447319) * y -0.009279453341) * y
              +0.005353579108) * y -0.002141268741) * y
            +0.000535310849) * y +0.999936657524;
		}
    }
		
    if (z > 0.0)
        return ((x + 1.0) * 0.5);
    return (1.0-x) * 0.5;
}

void MannWhitneyTest::setData(const Vector& xA, const Vector& xB) {
    unsigned long n;
    double uA;
    double nAnB;

    m_nA = xA.size(); m_nB = xB.size();
    n = m_nA + m_nB;
    nAnB = ((double)m_nA) * m_nB;

    Vector joined = xA;
    joined.append(xB);

    Vector ranks;
    ranks.resize(joined.size());
    rank_vector(joined.begin(), joined.end(), ranks.begin(), m_tolerance);

    uA = nAnB - std::accumulate(ranks.begin(), ranks.begin() + m_nA, 0.0);
    uA += m_nA / 2.0 * (m_nA + 1);

    if (m_alternative == NOT_EQUAL) {
        m_u = std::min(uA, nAnB - uA);
    } else {
        m_u = uA;
    }

    // Calculate tie correction value
    m_tieCorrection = 0.0;
    if (ranks.size() > 1) {
        ranks.sort();

        Vector::iterator curr = ranks.begin(), end = ranks.end();
        Vector::iterator next = curr + 1;
        while (next != end) {
            if (*curr == *next) {
                int nties = 1;
                while (next != end && *curr == *next) {
                    curr++; next++; nties++;
                }
                m_tieCorrection += std::pow(nties, 3.0) - nties;
            }
            if (next != end) {
                curr++; next++;
            }
        }
    }
    m_tieCorrection = std::sqrt(1 - m_tieCorrection / (std::pow(n, 3.0) - n)); 
}

double MannWhitneyTest::getP() const {
    double nAnB = ((double)m_nA) * m_nB;
    double sd = std::sqrt(m_tieCorrection / 12.0 * nAnB * (m_nA + m_nB + 1));
    double z = (m_u - nAnB / 2.0) / sd;

    if (m_alternative == NOT_EQUAL) {
        z = -std::fabs(z);
        return 2 * z_probability(z);
    } else if (m_alternative == LESS_THAN) {
        return 1.0 - z_probability(z);
    } else {
        return z_probability(z);
    }
}
