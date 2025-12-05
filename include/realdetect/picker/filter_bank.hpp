#pragma once

#include "../core/types.hpp"
#include <vector>
#include <complex>

namespace realdetect {

/**
 * IIRFilter - Infinite Impulse Response digital filter
 */
class IIRFilter {
public:
    IIRFilter() = default;
    
    // Create Butterworth bandpass filter
    static IIRFilter butterworth(int order, double low_freq, double high_freq,
                                  double sample_rate);
    
    // Create Butterworth lowpass filter
    static IIRFilter butterworthLowpass(int order, double cutoff, 
                                         double sample_rate);
    
    // Create Butterworth highpass filter
    static IIRFilter butterworthHighpass(int order, double cutoff,
                                          double sample_rate);
    
    // Apply filter to data (in-place)
    void apply(SampleVector& data) const;
    
    // Apply filter (copy)
    SampleVector filter(const SampleVector& data) const;
    
    // Zero-phase filtering (forward-backward)
    SampleVector filtfilt(const SampleVector& data) const;
    
    // Reset filter state
    void reset();
    
    // Filter coefficients access
    const std::vector<double>& b() const { return b_; }
    const std::vector<double>& a() const { return a_; }

private:
    std::vector<double> b_;  // Numerator coefficients
    std::vector<double> a_;  // Denominator coefficients
    mutable std::vector<double> z_;  // Filter state
    
    // Bilinear transform
    static void bilinear(const std::vector<std::complex<double>>& poles,
                         const std::vector<std::complex<double>>& zeros,
                         double k, double fs,
                         std::vector<double>& b, std::vector<double>& a);
};

/**
 * FilterBank - Collection of filters for multiband analysis
 */
class FilterBank {
public:
    FilterBank() = default;
    
    // Add bandpass filter
    void addBand(double low_freq, double high_freq, int order = 4);
    
    // Configure for sample rate
    void configure(double sample_rate);
    
    // Filter data through all bands
    std::vector<SampleVector> filter(const SampleVector& data) const;
    
    // Get combined envelope from all bands
    SampleVector combinedEnvelope(const SampleVector& data) const;
    
    // Band info
    size_t bandCount() const { return bands_.size(); }
    std::pair<double, double> bandFrequencies(size_t idx) const {
        if (idx < bands_.size()) return bands_[idx];
        return {0, 0};
    }
    
    // Clear all bands
    void clear();

private:
    double sample_rate_ = 100.0;
    std::vector<std::pair<double, double>> bands_;
    std::vector<IIRFilter> filters_;
    int order_ = 4;
};

/**
 * Simple FFT implementation for spectral analysis
 */
class FFT {
public:
    // Forward FFT
    static std::vector<std::complex<double>> forward(const SampleVector& data);
    
    // Inverse FFT
    static SampleVector inverse(const std::vector<std::complex<double>>& spectrum);
    
    // Power spectrum
    static SampleVector powerSpectrum(const SampleVector& data);
    
    // Amplitude spectrum
    static SampleVector amplitudeSpectrum(const SampleVector& data);
    
    // Frequency axis
    static SampleVector frequencies(size_t n, double sample_rate);
    
    // Dominant frequency
    static double dominantFrequency(const SampleVector& data, double sample_rate);

    // Cooley-Tukey FFT (public for use by other modules)
    static void fft(std::vector<std::complex<double>>& x, bool inverse = false);
    
    // Next power of 2
    static size_t nextPow2(size_t n);
};

} // namespace realdetect
