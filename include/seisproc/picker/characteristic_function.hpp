#pragma once

#include "../core/types.hpp"
#include "../core/waveform.hpp"

namespace seisproc {

/**
 * CharacteristicFunction - Various CF implementations for phase detection
 */
class CharacteristicFunction {
public:
    /**
     * Envelope (Hilbert transform based)
     * Good for detecting amplitude changes
     */
    static SampleVector envelope(const SampleVector& data);
    
    /**
     * Kurtosis - 4th moment, sensitive to impulsive arrivals
     * Computed in sliding window
     */
    static SampleVector kurtosis(const SampleVector& data, size_t window_size);
    
    /**
     * Modified Energy Ratio
     * CF = sum(x[i]^2) / sum(x[i-n]^2) for sliding window
     */
    static SampleVector energyRatio(const SampleVector& data, size_t window_size);
    
    /**
     * Polarization - For 3-component data
     * Detects linear particle motion (P-waves) vs elliptical (S-waves)
     */
    static SampleVector polarization(const SampleVector& z,
                                      const SampleVector& n,
                                      const SampleVector& e,
                                      size_t window_size);
    
    /**
     * Rectilinearity from polarization analysis
     * High for P-waves, low for S-waves
     */
    static SampleVector rectilinearity(const SampleVector& z,
                                         const SampleVector& n,
                                         const SampleVector& e,
                                         size_t window_size);
    
    /**
     * Carl-Sta-Trig function
     * Combines envelope and spectral information
     */
    static SampleVector carlStaTrig(const SampleVector& data, 
                                     double sample_rate,
                                     double sta_len, double lta_len);
    
    /**
     * Z-detector (cumulative sum)
     * Good for emergent arrivals
     */
    static SampleVector zDetector(const SampleVector& data, size_t window_size);
    
    /**
     * Skewness - 3rd moment
     * Sensitive to asymmetric arrivals
     */
    static SampleVector skewness(const SampleVector& data, size_t window_size);
    
    /**
     * Higher-order statistics CF
     * Combines multiple statistical measures
     */
    static SampleVector higherOrderStats(const SampleVector& data, 
                                          size_t window_size);
};

/**
 * Hilbert transform for envelope computation
 */
class HilbertTransform {
public:
    /**
     * Compute analytic signal using FFT
     */
    static SampleVector analyticSignal(const SampleVector& data);
    
    /**
     * Compute instantaneous phase
     */
    static SampleVector instantaneousPhase(const SampleVector& data);
    
    /**
     * Compute instantaneous frequency
     */
    static SampleVector instantaneousFrequency(const SampleVector& data, 
                                                 double sample_rate);
    
    /**
     * Compute envelope (magnitude of analytic signal)
     */
    static SampleVector envelope(const SampleVector& data);
};

} // namespace seisproc
