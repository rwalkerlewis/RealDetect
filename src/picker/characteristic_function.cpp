#include "seisproc/picker/characteristic_function.hpp"
#include "seisproc/picker/filter_bank.hpp"
#include <cmath>
#include <algorithm>
#include <numeric>

namespace seisproc {

SampleVector CharacteristicFunction::envelope(const SampleVector& data) {
    return HilbertTransform::envelope(data);
}

SampleVector CharacteristicFunction::kurtosis(const SampleVector& data, 
                                               size_t window_size) {
    size_t n = data.size();
    SampleVector result(n, 0);
    
    if (n < window_size) return result;
    
    for (size_t i = window_size; i < n; i++) {
        // Compute mean
        double sum = 0;
        for (size_t j = i - window_size; j < i; j++) {
            sum += data[j];
        }
        double mean = sum / window_size;
        
        // Compute 2nd and 4th moments
        double m2 = 0, m4 = 0;
        for (size_t j = i - window_size; j < i; j++) {
            double d = data[j] - mean;
            m2 += d * d;
            m4 += d * d * d * d;
        }
        m2 /= window_size;
        m4 /= window_size;
        
        // Kurtosis (excess)
        if (m2 > 1e-20) {
            result[i] = m4 / (m2 * m2) - 3.0;
        }
    }
    
    return result;
}

SampleVector CharacteristicFunction::energyRatio(const SampleVector& data,
                                                   size_t window_size) {
    size_t n = data.size();
    SampleVector result(n, 0);
    
    if (n < 2 * window_size) return result;
    
    for (size_t i = 2 * window_size; i < n; i++) {
        double pre_energy = 0, post_energy = 0;
        
        for (size_t j = i - 2 * window_size; j < i - window_size; j++) {
            pre_energy += data[j] * data[j];
        }
        for (size_t j = i - window_size; j < i; j++) {
            post_energy += data[j] * data[j];
        }
        
        if (pre_energy > 1e-20) {
            result[i] = post_energy / pre_energy;
        }
    }
    
    return result;
}

SampleVector CharacteristicFunction::polarization(const SampleVector& z,
                                                    const SampleVector& n,
                                                    const SampleVector& e,
                                                    size_t window_size) {
    // Polarization analysis using covariance matrix eigenvalue decomposition
    size_t len = std::min({z.size(), n.size(), e.size()});
    SampleVector result(len, 0);
    
    if (len < window_size) return result;
    
    for (size_t i = window_size; i < len; i++) {
        // Build covariance matrix
        double cov[3][3] = {{0}};
        
        for (size_t j = i - window_size; j < i; j++) {
            double vals[3] = {z[j], n[j], e[j]};
            for (int a = 0; a < 3; a++) {
                for (int b = 0; b < 3; b++) {
                    cov[a][b] += vals[a] * vals[b];
                }
            }
        }
        
        for (int a = 0; a < 3; a++) {
            for (int b = 0; b < 3; b++) {
                cov[a][b] /= window_size;
            }
        }
        
        // Compute eigenvalues (simple power iteration for largest)
        double trace = cov[0][0] + cov[1][1] + cov[2][2];
        
        // Simplified - use trace as polarization measure
        result[i] = std::sqrt(trace);
    }
    
    return result;
}

SampleVector CharacteristicFunction::rectilinearity(const SampleVector& z,
                                                      const SampleVector& n,
                                                      const SampleVector& e,
                                                      size_t window_size) {
    // Rectilinearity = 1 - (lambda2 + lambda3) / (2 * lambda1)
    // where lambda1 >= lambda2 >= lambda3 are eigenvalues
    
    size_t len = std::min({z.size(), n.size(), e.size()});
    SampleVector result(len, 0);
    
    if (len < window_size) return result;
    
    for (size_t i = window_size; i < len; i++) {
        // Build covariance matrix
        double cov[3][3] = {{0}};
        double mean[3] = {0};
        
        for (size_t j = i - window_size; j < i; j++) {
            mean[0] += z[j];
            mean[1] += n[j];
            mean[2] += e[j];
        }
        mean[0] /= window_size;
        mean[1] /= window_size;
        mean[2] /= window_size;
        
        for (size_t j = i - window_size; j < i; j++) {
            double vals[3] = {z[j] - mean[0], n[j] - mean[1], e[j] - mean[2]};
            for (int a = 0; a < 3; a++) {
                for (int b = 0; b < 3; b++) {
                    cov[a][b] += vals[a] * vals[b];
                }
            }
        }
        
        // Compute characteristic polynomial coefficients
        double I1 = cov[0][0] + cov[1][1] + cov[2][2];  // trace
        double I2 = cov[0][0]*cov[1][1] + cov[1][1]*cov[2][2] + cov[0][0]*cov[2][2]
                  - cov[0][1]*cov[0][1] - cov[1][2]*cov[1][2] - cov[0][2]*cov[0][2];
        
        // Approximate rectilinearity using I1 and I2
        if (I1 > 1e-20) {
            result[i] = 1.0 - 2.0 * std::sqrt(std::max(0.0, I2)) / I1;
            result[i] = std::max(0.0, std::min(1.0, result[i]));
        }
    }
    
    return result;
}

SampleVector CharacteristicFunction::carlStaTrig(const SampleVector& data,
                                                   double sample_rate,
                                                   double sta_len, double lta_len) {
    // Carl-Sta-Trig combines envelope and STA/LTA
    SampleVector env = envelope(data);
    
    size_t sta_samples = static_cast<size_t>(sta_len * sample_rate);
    size_t lta_samples = static_cast<size_t>(lta_len * sample_rate);
    
    size_t n = env.size();
    SampleVector result(n, 0);
    
    if (n < lta_samples) return result;
    
    double sta = 0, lta = 0;
    
    // Initialize
    for (size_t i = 0; i < lta_samples; i++) {
        lta += env[i];
    }
    lta /= lta_samples;
    
    for (size_t i = lta_samples - sta_samples; i < lta_samples; i++) {
        sta += env[i];
    }
    sta /= sta_samples;
    
    // Recursive computation
    for (size_t i = lta_samples; i < n; i++) {
        sta = sta + (env[i] - env[i - sta_samples]) / sta_samples;
        lta = lta + (env[i] - env[i - lta_samples]) / lta_samples;
        
        if (lta > 1e-20) {
            result[i] = (sta - lta) / lta;
        }
    }
    
    return result;
}

SampleVector CharacteristicFunction::zDetector(const SampleVector& data, 
                                                 size_t window_size) {
    // Z-detector: cumulative sum of squared samples normalized
    size_t n = data.size();
    SampleVector result(n, 0);
    
    if (n < window_size) return result;
    
    // Compute global mean and variance
    double global_mean = 0, global_var = 0;
    for (size_t i = 0; i < n; i++) {
        global_mean += data[i] * data[i];
    }
    global_mean /= n;
    
    for (size_t i = 0; i < n; i++) {
        double d = data[i] * data[i] - global_mean;
        global_var += d * d;
    }
    global_var = std::sqrt(global_var / n);
    
    if (global_var < 1e-20) return result;
    
    // Cumulative sum
    double cum_sum = 0;
    for (size_t i = 0; i < n; i++) {
        cum_sum += (data[i] * data[i] - global_mean) / global_var;
        result[i] = cum_sum;
    }
    
    return result;
}

SampleVector CharacteristicFunction::skewness(const SampleVector& data,
                                               size_t window_size) {
    size_t n = data.size();
    SampleVector result(n, 0);
    
    if (n < window_size) return result;
    
    for (size_t i = window_size; i < n; i++) {
        double sum = 0;
        for (size_t j = i - window_size; j < i; j++) {
            sum += data[j];
        }
        double mean = sum / window_size;
        
        double m2 = 0, m3 = 0;
        for (size_t j = i - window_size; j < i; j++) {
            double d = data[j] - mean;
            m2 += d * d;
            m3 += d * d * d;
        }
        m2 /= window_size;
        m3 /= window_size;
        
        double sigma = std::sqrt(m2);
        if (sigma > 1e-20) {
            result[i] = m3 / (sigma * sigma * sigma);
        }
    }
    
    return result;
}

SampleVector CharacteristicFunction::higherOrderStats(const SampleVector& data,
                                                        size_t window_size) {
    // Combine kurtosis and skewness
    SampleVector kurt = kurtosis(data, window_size);
    SampleVector skew = skewness(data, window_size);
    
    SampleVector result(data.size(), 0);
    for (size_t i = 0; i < data.size(); i++) {
        result[i] = std::sqrt(kurt[i] * kurt[i] + skew[i] * skew[i]);
    }
    
    return result;
}

// HilbertTransform implementation

SampleVector HilbertTransform::analyticSignal(const SampleVector& data) {
    // Compute Hilbert transform using FFT
    auto spectrum = FFT::forward(data);
    size_t n = spectrum.size();
    
    // Zero negative frequencies, double positive frequencies
    for (size_t i = 1; i < n / 2; i++) {
        spectrum[i] *= 2.0;
    }
    for (size_t i = n / 2 + 1; i < n; i++) {
        spectrum[i] = 0;
    }
    
    // Inverse FFT gives analytic signal
    return FFT::inverse(spectrum);
}

SampleVector HilbertTransform::envelope(const SampleVector& data) {
    auto spectrum = FFT::forward(data);
    size_t n = spectrum.size();
    
    // Create analytic signal spectrum
    for (size_t i = 1; i < n / 2; i++) {
        spectrum[i] *= 2.0;
    }
    for (size_t i = n / 2 + 1; i < n; i++) {
        spectrum[i] = 0;
    }
    
    // Get analytic signal (complex)
    std::vector<std::complex<double>> analytic = spectrum;
    
    // Inverse FFT
    FFT::fft(analytic, true);
    
    // Envelope is magnitude
    SampleVector env(data.size());
    for (size_t i = 0; i < data.size(); i++) {
        env[i] = std::abs(analytic[i]);
    }
    
    return env;
}

SampleVector HilbertTransform::instantaneousPhase(const SampleVector& data) {
    auto spectrum = FFT::forward(data);
    size_t n = spectrum.size();
    
    for (size_t i = 1; i < n / 2; i++) {
        spectrum[i] *= 2.0;
    }
    for (size_t i = n / 2 + 1; i < n; i++) {
        spectrum[i] = 0;
    }
    
    std::vector<std::complex<double>> analytic = spectrum;
    FFT::fft(analytic, true);
    
    SampleVector phase(data.size());
    for (size_t i = 0; i < data.size(); i++) {
        phase[i] = std::arg(analytic[i]);
    }
    
    return phase;
}

SampleVector HilbertTransform::instantaneousFrequency(const SampleVector& data,
                                                        double sample_rate) {
    SampleVector phase = instantaneousPhase(data);
    SampleVector freq(data.size(), 0);
    
    // Unwrap phase
    for (size_t i = 1; i < phase.size(); i++) {
        double diff = phase[i] - phase[i-1];
        while (diff > M_PI) diff -= 2 * M_PI;
        while (diff < -M_PI) diff += 2 * M_PI;
        freq[i] = diff * sample_rate / (2 * M_PI);
    }
    freq[0] = freq[1];
    
    return freq;
}

} // namespace seisproc
