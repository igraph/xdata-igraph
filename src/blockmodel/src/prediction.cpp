/* vim:set ts=4 sw=4 sts=4 et: */

#include <block/prediction.h>
#include <igraph/cpp/matrix.h>
#include <igraph/cpp/vector.h>

using namespace igraph;

double AveragingPredictor::predictProbability(int v1, int v2) {
    if (v1 > v2)
        return m_counts(v2, v1) / m_numSamples;
    else
        return m_counts(v1, v2) / m_numSamples;
}

bool AveragingPredictor::takeSampleReal() {
    size_t n = m_pModel->getVertexCount();
    for (size_t v1 = 0; v1 < n; v1++)
        for (size_t v2 = v1+1; v2 < n; v2++)
            m_counts(v1, v2) += m_pModel->getEdgeProbability(v1, v2);

    return true;
}

