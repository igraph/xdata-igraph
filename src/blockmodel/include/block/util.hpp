/* vim:set ts=4 sw=4 sts=4 et: */

#ifndef BLOCKMODEL_UTIL_HPP
#define BLOCKMODEL_UTIL_HPP

#include <algorithm>
#include <cmath>
#include <iterator>
#include <vector>

/// Calculates the value of the Akaike information criterion for some model
template <typename T>
double aic(const T& model) {
    return 2 * (model.getNumParameters() - model.getLogLikelihood());
}

/// Calculates the value of the Bayesian information criterion for some model
template <typename T>
double bic(const T& model) {
    return -2 * model.getLogLikelihood() +
		   model.getNumParameters() *
           std::log(model.getNumObservations());
}

/// Identity map
template <typename Key=int>
struct identity_map {
    Key operator[](const Key& x) {
        return x;
    }

    const Key& operator[](const Key& x) const {
        return x;
    }
};

/// Comparator that can be used to sort an index vector of a container
template <typename Map, typename Element=typename std::iterator_traits<Map>::value_type>
struct key_comparator {
    const Map& m_map;

    key_comparator(const Map& map) : m_map(map) {}

    bool operator()(Element e1, Element e2) {
        return m_map[e1] < m_map[e2];
    }
};

/// Comparator that can be used to sort pairs
template <int I>
struct pair_comparator {
    template <typename Pair>
    bool operator()(Pair p1, Pair p2);
};

/// Specialization of pair_comparator for sorting based on the first element
template <>
struct pair_comparator<1> {
    template <typename Pair>
    bool operator()(const Pair& p1, const Pair& p2) {
        return p1.first < p2.first;
    }
};

/// Specialization of pair_comparator for sorting based on the second element
template <>
struct pair_comparator<2> {
    template <typename Pair>
    bool operator()(const Pair& p1, const Pair& p2) {
        return p1.second < p2.second;
    }
};

/// Auxiliary class that maps the function call operator to the indexing operator
template <typename Container>
struct function_to_indexing {
public:
    typedef const typename std::iterator_traits<Container>::reference const_reference;
    typedef typename std::iterator_traits<Container>::difference_type size_type;
    typedef typename std::iterator_traits<Container>::reference       reference;

public:
    const Container& m_ref;

    function_to_indexing(const Container& ref) : m_ref(ref) {}

    reference operator()(size_type index) { return m_ref[index]; }
    const_reference operator()(size_type index) const { return m_ref[index]; }
};

/// Permuted view of a container
template <typename Container, typename Permutation>
class PermutedContainer {
public:
    typedef const typename std::iterator_traits<Container>::reference const_reference;
    typedef typename std::iterator_traits<Container>::difference_type size_type;
    typedef typename std::iterator_traits<Container>::reference       reference;
    typedef typename std::iterator_traits<Container>::value_type      value_type;

private:
    const Container& m_container;
    const Permutation& m_permutation;

public:
    PermutedContainer(const Container& container,
                      const Permutation& permutation)
        : m_container(container), m_permutation(permutation) {}

    reference operator[](size_type n) {
        return m_container[m_permutation[n]];
    }

    const_reference operator[](size_type n) const {
        return m_container[m_permutation[n]];
    }

    size_type size() const {
        return m_container.size();
    }
};

/// Calculates the rank vector of some elements in a container
template <typename InputIterator, typename OutputIterator>
void rank_vector(InputIterator begin, InputIterator end, OutputIterator out,
                 double tolerance = 0.0) {
    unsigned long i, j, n = end-begin;
    typedef typename std::iterator_traits<InputIterator>::value_type value_type;
    std::vector<long> ranks;
    std::vector<long> indices(end-begin);

    if (n == 0)
        return;

    i = 0;
    for (std::vector<long>::iterator it = indices.begin(); it != indices.end(); it++) {
        *it = i++;
    }
    std::sort(indices.begin(), indices.end(),
              key_comparator<InputIterator, long>(begin));

    unsigned long sumRanks = 0;
    unsigned long dupCount = 1;
    double rank;
    value_type curr, prev;
    PermutedContainer<InputIterator, std::vector<long> > sortedData(begin, indices);

    prev = sortedData[0];
    for (i = 1; i < n; i++) {
        curr = sortedData[i];

        if (std::fabs(curr - prev) <= tolerance) {
            dupCount++;
            sumRanks += i;
        } else {
            rank = (double)sumRanks / dupCount + 1;
            for (j = i - dupCount; j < i; j++)
                out[indices[j]] = rank;
            dupCount = 1;
            sumRanks = i;
        }

        prev = curr;
    }

    rank = (double)sumRanks / dupCount + 1;
    for (j = n - dupCount; j < n; j++)
        out[indices[j]] = rank;
}

#endif

