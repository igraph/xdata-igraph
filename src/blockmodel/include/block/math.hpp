/* vim:set ts=4 sw=4 sts=4 et: */

#ifndef BLOCKMODEL_MATH_H
#define BLOCKMODEL_MATH_H

#include <cmath>
#include <cstdio>
#include <numeric>
#include <vector>

#ifndef isnan
template <typename T>
inline bool block_isnan(T x) {
    return x != x;
}

#define isnan(x) block_isnan(x)
#endif

/// Calculates the moving average of some time series
template <typename T>
class MovingAverage {
private:
    /// Vector to store the values
    std::vector<T> m_data;
    
    /// Current window size
    float m_windowSize;

    /// Iterator pointing to the next element of the vector
    typename std::vector<T>::iterator m_next;

    /// Iterator pointing to the end of the vector
    typename std::vector<T>::iterator m_end;

    /// Current value of the sum
    float m_sum;

public:
    /// Constructs a moving average calculator with the given window size
    MovingAverage(int size) : m_data(size), m_windowSize(size),
        m_next(m_data.begin()), m_end(m_data.end()), m_sum() {}

    /// Adds a new entry to the internal storage
    void push_back(T value) {
        m_sum += (value - *m_next);
        *m_next = value;
        m_next++;
        if (m_next == m_end) {
            // Recalculate the sum here to avoid numeric drift
            m_next = m_data.begin();
            m_sum = std::accumulate(m_next, m_data.end(), 0);
        }
    }

    /// Returns the calculated average
    float value() const {
        return m_sum / m_windowSize;
    }
};

#endif

