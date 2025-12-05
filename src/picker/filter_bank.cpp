#include "realdetect/picker/filter_bank.hpp"
#include <cmath>
#include <algorithm>

namespace realdetect {

// IIRFilter implementation

IIRFilter IIRFilter::butterworth(int order, double low_freq, double high_freq,
                                   double sample_rate) {
    IIRFilter filter;
    
    // Prewarp frequencies
    double nyq = sample_rate / 2.0;
    double wl = std::tan(M_PI * low_freq / sample_rate);
    double wh = std::tan(M_PI * high_freq / sample_rate);
    double bw = wh - wl;
    double w0 = std::sqrt(wl * wh);
    
    // Generate analog prototype poles
    std::vector<std::complex<double>> poles;
    for (int k = 0; k < order; k++) {
        double theta = M_PI * (2.0 * k + 1) / (2.0 * order) + M_PI / 2.0;
        poles.push_back(std::complex<double>(std::cos(theta), std::sin(theta)));
    }
    
    // Transform to bandpass
    std::vector<std::complex<double>> bp_poles;
    std::vector<std::complex<double>> bp_zeros;
    
    for (const auto& p : poles) {
        std::complex<double> scaled = p * bw / 2.0;
        std::complex<double> disc = scaled * scaled - w0 * w0;
        std::complex<double> sq = std::sqrt(disc);
        bp_poles.push_back(scaled + sq);
        bp_poles.push_back(scaled - sq);
    }
    
    // Zeros at origin for bandpass
    for (int i = 0; i < order; i++) {
        bp_zeros.push_back(0);
        bp_zeros.push_back(0);
    }
    
    // Gain at center frequency
    double k = std::pow(bw, order);
    
    // Bilinear transform
    bilinear(bp_poles, bp_zeros, k, sample_rate, filter.b_, filter.a_);
    
    // Normalize
    double sum_b = 0, sum_a = 0;
    for (auto x : filter.b_) sum_b += x;
    for (auto x : filter.a_) sum_a += x;
    
    // Initialize state
    filter.z_.resize(std::max(filter.b_.size(), filter.a_.size()) - 1, 0);
    
    return filter;
}

IIRFilter IIRFilter::butterworthLowpass(int order, double cutoff, double sample_rate) {
    IIRFilter filter;
    
    // Prewarp frequency
    double w = std::tan(M_PI * cutoff / sample_rate);
    
    // Generate poles
    std::vector<std::complex<double>> poles;
    for (int k = 0; k < order; k++) {
        double theta = M_PI * (2.0 * k + 1) / (2.0 * order) + M_PI / 2.0;
        poles.push_back(w * std::complex<double>(std::cos(theta), std::sin(theta)));
    }
    
    // No zeros for lowpass (all at -infinity)
    std::vector<std::complex<double>> zeros;
    
    double k = std::pow(w, order);
    bilinear(poles, zeros, k, sample_rate, filter.b_, filter.a_);
    filter.z_.resize(std::max(filter.b_.size(), filter.a_.size()) - 1, 0);
    
    return filter;
}

IIRFilter IIRFilter::butterworthHighpass(int order, double cutoff, double sample_rate) {
    IIRFilter filter;
    
    // Prewarp frequency
    double w = std::tan(M_PI * cutoff / sample_rate);
    
    // Generate poles (transformed from lowpass)
    std::vector<std::complex<double>> poles;
    for (int k = 0; k < order; k++) {
        double theta = M_PI * (2.0 * k + 1) / (2.0 * order) + M_PI / 2.0;
        std::complex<double> lp_pole(std::cos(theta), std::sin(theta));
        poles.push_back(w / lp_pole);
    }
    
    // Zeros at origin
    std::vector<std::complex<double>> zeros(order, 0);
    
    double k = 1.0;
    bilinear(poles, zeros, k, sample_rate, filter.b_, filter.a_);
    filter.z_.resize(std::max(filter.b_.size(), filter.a_.size()) - 1, 0);
    
    return filter;
}

void IIRFilter::bilinear(const std::vector<std::complex<double>>& poles,
                          const std::vector<std::complex<double>>& zeros,
                          double k, double fs,
                          std::vector<double>& b, std::vector<double>& a) {
    // Simplified bilinear transform
    // For a complete implementation, use polynomial multiplication
    
    size_t n = std::max(poles.size(), zeros.size());
    b.resize(n + 1, 0);
    a.resize(n + 1, 0);
    
    // Start with constant
    std::vector<std::complex<double>> num_poly = {k};
    std::vector<std::complex<double>> den_poly = {1.0};
    
    // Multiply in zeros: (1 - z * z^-1) for each zero
    for (const auto& z : zeros) {
        std::complex<double> z_d = (1.0 + z / fs) / (1.0 - z / fs);
        std::vector<std::complex<double>> factor = {-z_d, 1.0};
        
        std::vector<std::complex<double>> new_poly(num_poly.size() + 1, 0);
        for (size_t i = 0; i < num_poly.size(); i++) {
            for (size_t j = 0; j < factor.size(); j++) {
                new_poly[i + j] += num_poly[i] * factor[j];
            }
        }
        num_poly = new_poly;
    }
    
    // Multiply in poles
    for (const auto& p : poles) {
        std::complex<double> p_d = (1.0 + p / fs) / (1.0 - p / fs);
        std::vector<std::complex<double>> factor = {-p_d, 1.0};
        
        std::vector<std::complex<double>> new_poly(den_poly.size() + 1, 0);
        for (size_t i = 0; i < den_poly.size(); i++) {
            for (size_t j = 0; j < factor.size(); j++) {
                new_poly[i + j] += den_poly[i] * factor[j];
            }
        }
        den_poly = new_poly;
    }
    
    // Extract real parts (should be real for symmetric poles)
    b.resize(num_poly.size());
    a.resize(den_poly.size());
    
    for (size_t i = 0; i < num_poly.size(); i++) {
        b[i] = num_poly[i].real();
    }
    for (size_t i = 0; i < den_poly.size(); i++) {
        a[i] = den_poly[i].real();
    }
    
    // Normalize by a[0]
    double norm = a[0];
    for (auto& x : b) x /= norm;
    for (auto& x : a) x /= norm;
}

void IIRFilter::apply(SampleVector& data) const {
    SampleVector out = filter(data);
    data = out;
}

SampleVector IIRFilter::filter(const SampleVector& data) const {
    if (b_.empty() || a_.empty()) return data;
    
    size_t n = data.size();
    size_t nb = b_.size();
    size_t na = a_.size();
    
    SampleVector output(n, 0);
    std::vector<double> z = z_;  // Copy state
    
    for (size_t i = 0; i < n; i++) {
        // Direct Form II transposed
        output[i] = b_[0] * data[i] + (z.empty() ? 0 : z[0]);
        
        for (size_t j = 1; j < std::max(nb, na); j++) {
            double b_val = j < nb ? b_[j] : 0;
            double a_val = j < na ? a_[j] : 0;
            
            if (j - 1 < z.size()) {
                z[j - 1] = b_val * data[i] - a_val * output[i];
                if (j < z.size()) {
                    z[j - 1] += z[j];
                }
            }
        }
    }
    
    return output;
}

SampleVector IIRFilter::filtfilt(const SampleVector& data) const {
    // Forward filter
    SampleVector forward = filter(data);
    
    // Reverse
    std::reverse(forward.begin(), forward.end());
    
    // Backward filter
    SampleVector backward = filter(forward);
    
    // Reverse again
    std::reverse(backward.begin(), backward.end());
    
    return backward;
}

void IIRFilter::reset() {
    std::fill(z_.begin(), z_.end(), 0);
}

// FilterBank implementation

void FilterBank::addBand(double low_freq, double high_freq, int order) {
    bands_.push_back({low_freq, high_freq});
    order_ = order;
}

void FilterBank::configure(double sample_rate) {
    sample_rate_ = sample_rate;
    filters_.clear();
    
    for (const auto& band : bands_) {
        filters_.push_back(IIRFilter::butterworth(order_, band.first, band.second, 
                                                   sample_rate));
    }
}

std::vector<SampleVector> FilterBank::filter(const SampleVector& data) const {
    std::vector<SampleVector> result;
    
    for (const auto& filt : filters_) {
        result.push_back(filt.filtfilt(data));
    }
    
    return result;
}

SampleVector FilterBank::combinedEnvelope(const SampleVector& data) const {
    if (filters_.empty()) return data;
    
    auto filtered = filter(data);
    
    SampleVector combined(data.size(), 0);
    
    for (const auto& band : filtered) {
        // Compute envelope (simplified as absolute value)
        for (size_t i = 0; i < data.size(); i++) {
            combined[i] += band[i] * band[i];
        }
    }
    
    for (auto& x : combined) {
        x = std::sqrt(x);
    }
    
    return combined;
}

void FilterBank::clear() {
    bands_.clear();
    filters_.clear();
}

// FFT implementation

std::vector<std::complex<double>> FFT::forward(const SampleVector& data) {
    size_t n = nextPow2(data.size());
    std::vector<std::complex<double>> x(n);
    
    for (size_t i = 0; i < data.size(); i++) {
        x[i] = data[i];
    }
    for (size_t i = data.size(); i < n; i++) {
        x[i] = 0;
    }
    
    fft(x, false);
    return x;
}

SampleVector FFT::inverse(const std::vector<std::complex<double>>& spectrum) {
    std::vector<std::complex<double>> x = spectrum;
    fft(x, true);
    
    SampleVector result(x.size());
    for (size_t i = 0; i < x.size(); i++) {
        result[i] = x[i].real();
    }
    return result;
}

void FFT::fft(std::vector<std::complex<double>>& x, bool inverse) {
    size_t n = x.size();
    if (n <= 1) return;
    
    // Bit-reversal permutation
    for (size_t i = 1, j = 0; i < n; i++) {
        size_t bit = n >> 1;
        while (j & bit) {
            j ^= bit;
            bit >>= 1;
        }
        j ^= bit;
        if (i < j) std::swap(x[i], x[j]);
    }
    
    // Cooley-Tukey FFT
    for (size_t len = 2; len <= n; len *= 2) {
        double angle = 2 * M_PI / len * (inverse ? 1 : -1);
        std::complex<double> wlen(std::cos(angle), std::sin(angle));
        
        for (size_t i = 0; i < n; i += len) {
            std::complex<double> w = 1;
            for (size_t j = 0; j < len / 2; j++) {
                std::complex<double> u = x[i + j];
                std::complex<double> v = x[i + j + len/2] * w;
                x[i + j] = u + v;
                x[i + j + len/2] = u - v;
                w *= wlen;
            }
        }
    }
    
    if (inverse) {
        for (auto& val : x) val /= n;
    }
}

SampleVector FFT::powerSpectrum(const SampleVector& data) {
    auto spec = forward(data);
    SampleVector power(spec.size() / 2);
    
    for (size_t i = 0; i < power.size(); i++) {
        power[i] = std::norm(spec[i]);
    }
    
    return power;
}

SampleVector FFT::amplitudeSpectrum(const SampleVector& data) {
    auto spec = forward(data);
    SampleVector amp(spec.size() / 2);
    
    for (size_t i = 0; i < amp.size(); i++) {
        amp[i] = std::abs(spec[i]);
    }
    
    return amp;
}

SampleVector FFT::frequencies(size_t n, double sample_rate) {
    SampleVector freq(n / 2);
    double df = sample_rate / n;
    
    for (size_t i = 0; i < freq.size(); i++) {
        freq[i] = i * df;
    }
    
    return freq;
}

double FFT::dominantFrequency(const SampleVector& data, double sample_rate) {
    auto power = powerSpectrum(data);
    auto freq = frequencies(nextPow2(data.size()), sample_rate);
    
    size_t max_idx = 0;
    double max_power = 0;
    
    for (size_t i = 1; i < power.size(); i++) {
        if (power[i] > max_power) {
            max_power = power[i];
            max_idx = i;
        }
    }
    
    return max_idx < freq.size() ? freq[max_idx] : 0;
}

size_t FFT::nextPow2(size_t n) {
    size_t p = 1;
    while (p < n) p *= 2;
    return p;
}

} // namespace realdetect
