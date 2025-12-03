#include "seisproc/picker/aic_picker.hpp"
#include "seisproc/picker/stalta_picker.hpp"
#include <cmath>
#include <algorithm>
#include <numeric>

namespace seisproc {

AICPicker::AICPicker()
    : pre_window_(2.0)
    , post_window_(1.0)
    , min_snr_(3.0)
{
}

std::vector<PickResult> AICPicker::pick(const Waveform& waveform) {
    std::vector<PickResult> picks;
    
    // Find candidates using energy ratio
    std::vector<size_t> candidates = findCandidates(waveform);
    
    const SampleVector& data = waveform.data();
    
    for (size_t candidate : candidates) {
        // Refine using AIC
        size_t window = static_cast<size_t>((pre_window_ + post_window_) * 
                                             waveform.sampleRate());
        size_t refined = refinePick(waveform, candidate, window);
        
        PickResult pick;
        pick.time = waveform.timeAt(refined);
        pick.phase_type = PhaseType::P;
        pick.sample_index = refined;
        
        // Estimate SNR
        size_t noise_samples = static_cast<size_t>(pre_window_ * waveform.sampleRate());
        size_t start = refined > noise_samples ? refined - noise_samples : 0;
        
        double noise_energy = 0;
        for (size_t i = start; i < refined; i++) {
            noise_energy += data[i] * data[i];
        }
        noise_energy /= (refined - start);
        
        double signal_energy = 0;
        size_t signal_end = std::min(refined + noise_samples, data.size());
        for (size_t i = refined; i < signal_end; i++) {
            signal_energy += data[i] * data[i];
        }
        signal_energy /= (signal_end - refined);
        
        if (noise_energy > 0) {
            pick.snr = std::sqrt(signal_energy / noise_energy);
        }
        
        if (pick.snr >= min_snr_) {
            picks.push_back(pick);
        }
    }
    
    return picks;
}

std::vector<size_t> AICPicker::findCandidates(const Waveform& waveform) const {
    std::vector<size_t> candidates;
    
    const SampleVector& data = waveform.data();
    size_t window = static_cast<size_t>(waveform.sampleRate());  // 1 second window
    
    if (data.size() < 2 * window) return candidates;
    
    // Compute energy ratio
    for (size_t i = window; i < data.size() - window; i += window/4) {
        double pre_energy = 0, post_energy = 0;
        
        for (size_t j = i - window; j < i; j++) {
            pre_energy += data[j] * data[j];
        }
        for (size_t j = i; j < i + window && j < data.size(); j++) {
            post_energy += data[j] * data[j];
        }
        
        if (pre_energy > 0 && post_energy / pre_energy > min_snr_ * min_snr_) {
            candidates.push_back(i);
        }
    }
    
    return candidates;
}

size_t AICPicker::refinePick(const Waveform& waveform, size_t approx_idx,
                              size_t window_size) {
    const SampleVector& data = waveform.data();
    
    size_t start = approx_idx > window_size/2 ? approx_idx - window_size/2 : 0;
    size_t end = std::min(approx_idx + window_size/2, data.size());
    
    if (end - start < 10) return approx_idx;
    
    // Extract window
    SampleVector window_data(data.begin() + start, data.begin() + end);
    
    // Compute AIC and find minimum
    SampleVector aic = computeAIC(window_data);
    size_t min_idx = findAICMinimum(aic, 5, aic.size() - 5);
    
    return start + min_idx;
}

SampleVector AICPicker::computeAIC(const SampleVector& data) {
    size_t n = data.size();
    SampleVector aic(n, 0);
    
    for (size_t k = 2; k < n - 2; k++) {
        // Variance before k
        double sum1 = 0, sum1_sq = 0;
        for (size_t i = 0; i < k; i++) {
            sum1 += data[i];
            sum1_sq += data[i] * data[i];
        }
        double var1 = (sum1_sq - sum1*sum1/k) / k;
        if (var1 < 1e-20) var1 = 1e-20;
        
        // Variance after k
        double sum2 = 0, sum2_sq = 0;
        for (size_t i = k; i < n; i++) {
            sum2 += data[i];
            sum2_sq += data[i] * data[i];
        }
        double var2 = (sum2_sq - sum2*sum2/(n-k)) / (n-k);
        if (var2 < 1e-20) var2 = 1e-20;
        
        aic[k] = k * std::log(var1) + (n - k - 1) * std::log(var2);
    }
    
    return aic;
}

size_t AICPicker::findAICMinimum(const SampleVector& aic, size_t start, size_t end) {
    if (end == 0 || end > aic.size()) end = aic.size();
    if (start >= end) return start;
    
    double min_val = aic[start];
    size_t min_idx = start;
    
    for (size_t i = start + 1; i < end; i++) {
        if (aic[i] < min_val) {
            min_val = aic[i];
            min_idx = i;
        }
    }
    
    return min_idx;
}

void AICPicker::setParameter(const std::string& name, double value) {
    if (name == "pre_window") pre_window_ = value;
    else if (name == "post_window") post_window_ = value;
    else if (name == "min_snr") min_snr_ = value;
}

double AICPicker::getParameter(const std::string& name) const {
    if (name == "pre_window") return pre_window_;
    if (name == "post_window") return post_window_;
    if (name == "min_snr") return min_snr_;
    return 0;
}

// ARPicker implementation

ARPicker::ARPicker()
    : ar_order_(6)
    , pre_window_(2.0)
    , post_window_(1.0)
{
}

std::vector<PickResult> ARPicker::pick(const Waveform& waveform) {
    // Use AR prediction error as characteristic function
    SampleVector cf = computePredictionError(waveform.data());
    
    // Apply STA/LTA to the CF
    STALTAPicker stalta;
    stalta.setSTALength(0.3);
    stalta.setLTALength(5.0);
    stalta.setTriggerRatio(5.0);
    
    Waveform cf_waveform(waveform.streamId(), waveform.sampleRate(), 
                          waveform.startTime());
    cf_waveform.data() = cf;
    
    return stalta.pick(cf_waveform);
}

std::vector<double> ARPicker::burgAR(const SampleVector& data) const {
    // Burg's method for AR coefficient estimation
    int n = data.size();
    int p = ar_order_;
    
    std::vector<double> a(p + 1, 0);
    a[0] = 1.0;
    
    std::vector<double> ef(n);
    std::vector<double> eb(n);
    
    for (int i = 0; i < n; i++) {
        ef[i] = eb[i] = data[i];
    }
    
    for (int m = 1; m <= p; m++) {
        // Compute reflection coefficient
        double num = 0, den = 0;
        for (int j = m; j < n; j++) {
            num += ef[j] * eb[j - 1];
            den += ef[j] * ef[j] + eb[j - 1] * eb[j - 1];
        }
        
        double k = -2.0 * num / den;
        
        // Update AR coefficients
        std::vector<double> a_new(p + 1, 0);
        a_new[0] = 1.0;
        for (int i = 1; i <= m; i++) {
            a_new[i] = a[i] + k * a[m - i];
        }
        a = a_new;
        
        // Update prediction errors
        std::vector<double> ef_new(n), eb_new(n);
        for (int j = m; j < n; j++) {
            ef_new[j] = ef[j] + k * eb[j - 1];
            eb_new[j] = eb[j - 1] + k * ef[j];
        }
        ef = ef_new;
        eb = eb_new;
    }
    
    return std::vector<double>(a.begin() + 1, a.end());
}

SampleVector ARPicker::computePredictionError(const SampleVector& data) const {
    int n = data.size();
    int p = ar_order_;
    
    if (n <= p) return SampleVector(n, 0);
    
    // Estimate AR coefficients from first portion
    size_t train_size = std::min(static_cast<size_t>(n / 4), size_t(1000));
    SampleVector train_data(data.begin(), data.begin() + train_size);
    std::vector<double> a = burgAR(train_data);
    
    // Compute prediction error
    SampleVector error(n, 0);
    
    for (int i = p; i < n; i++) {
        double pred = 0;
        for (int j = 0; j < p; j++) {
            pred -= a[j] * data[i - 1 - j];
        }
        error[i] = std::abs(data[i] - pred);
    }
    
    return error;
}

void ARPicker::setParameter(const std::string& name, double value) {
    if (name == "ar_order") ar_order_ = static_cast<int>(value);
    else if (name == "pre_window") pre_window_ = value;
    else if (name == "post_window") post_window_ = value;
}

double ARPicker::getParameter(const std::string& name) const {
    if (name == "ar_order") return ar_order_;
    if (name == "pre_window") return pre_window_;
    if (name == "post_window") return post_window_;
    return 0;
}

} // namespace seisproc
