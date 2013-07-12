/* vim:set ts=4 sw=4 sts=4 et: */

#ifndef BLOCK_PREDICTION_H
#define BLOCK_PREDICTION_H

#include <block/blockmodel.h>
#include <igraph/cpp/matrix.h>
#include <igraph/cpp/vector.h>

/// Abstract predictor class for blockmodels
/**
 * Predictors take samples from the state of a model when their \c takeSample()
 * method is invoked, and then make predictions about the probability of a given
 * connection based on the previous samples taken. The prediction strategy is
 * entirely up to the concrete predictor implementation; it may use only the
 * last sample, or all the samples, or even none at all :)
 */
class Predictor {
protected:
    /// The model sampled by the predictor
    Blockmodel* m_pModel;

    /// The number of samples taken
    unsigned long m_numSamples;

public:
    /// Constructor
    explicit Predictor(Blockmodel* model) :
        m_pModel(model), m_numSamples(0) {}

    /// Virtual destructor that does nothing
    virtual ~Predictor() {}

    /// Returns the number of samples taken so far
    unsigned long getSampleCount() {
        return m_numSamples;
    }

    /// Predicts the probability of a connection
    virtual double predictProbability(int v1, int v2) = 0;

    /// Takes a sample using the current state of the model
    /**
     * This method also takes care of increasing the sample counter if the
     * current state was accepted for sampling.
     *
     * Returns \c true if the current state was used, \c false if it was
     * rejected.
     */
    bool takeSample() {
        bool result = takeSampleReal();
        if (result)
            m_numSamples++;
        return result;
    }

protected:
    /// Takes a sample using the current state of the model
    /**
     * Returns \c true if the current state was used, \c false if it was
     * rejected.
     */
    virtual bool takeSampleReal() = 0;
};

/// Predictor that averages over multiple samples
class AveragingPredictor : public Predictor {
private:
    /// A HUGE matrix that counts how many times a given edge was seen
    igraph::Matrix m_counts;

public:
    /// Constructor
    explicit AveragingPredictor(Blockmodel* model) :
        Predictor(model),
        m_counts(model->getVertexCount(), model->getVertexCount()) {
        m_counts.fill(0);
    }

    /// Predicts the probability of a connection
    virtual double predictProbability(int v1, int v2); 

protected:
    /// Takes a sample using the current state of the model
    /**
     * Returns \c true if the current state was used, \c false if it was
     * rejected.
     */
    virtual bool takeSampleReal();
};

#endif

