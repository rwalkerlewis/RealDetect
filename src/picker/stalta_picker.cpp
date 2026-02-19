#include "realdetect/picker/stalta_picker.hpp"
#include <cmath>
#include <algorithm>
#include <numeric>

namespace realdetect {

STALTAPicker::STALTAPicker()
    : sta_length_(0.5)
    , lta_length_(10.0)
    , trigger_ratio_(3.0)
    , detrigger_ratio_(1.5)
    , filter_low_(1.0)
    , filter_high_(20.0)
    , use_filter_(true)
{
}

std::vector<PickResult> STALTAPicker::pick(const Waveform& waveform) {
    std::vector<PickResult> picks;
    
    if (waveform.sampleCount() == 0) return picks;
    
    // Get data and optionally filter
    SampleVector data = waveform.data();
    
    // Demean
    double mean = std::accumulate(data.begin(), data.end(), 0.0) / data.size();
    for (auto& s : data) s -= mean;
    
    // Bandpass filter to remove long-period noise and high-frequency noise
    if (use_filter_ && filter_low_ > 0 && filter_high_ > 0 &&
        filter_high_ < waveform.sampleRate() / 2.0 &&
        data.size() > 100) {
        auto bp = IIRFilter::butterworth(2, filter_low_, filter_high_,
                                          waveform.sampleRate());
        data = bp.filtfilt(data);
    }
    
    // Compute STA/LTA ratio
    SampleVector cf = computeSTALTA(data, waveform.sampleRate());
    
    // Find triggers
    bool triggered = false;
    size_t trigger_idx = 0;
    double max_cf = 0;
    size_t max_idx = 0;
    
    for (size_t i = 0; i < cf.size(); i++) {
        if (!triggered) {
            if (cf[i] >= trigger_ratio_) {
                triggered = true;
                trigger_idx = i;
                max_cf = cf[i];
                max_idx = i;
            }
        } else {
            if (cf[i] > max_cf) {
                max_cf = cf[i];
                max_idx = i;
            }
            
            if (cf[i] < detrigger_ratio_) {
                // Found a pick
                size_t pick_idx = refinePick(data, trigger_idx);
                
                PickResult pick;
                pick.time = waveform.timeAt(pick_idx);
                pick.phase_type = PhaseType::P;  // Assume P for now
                pick.sample_index = pick_idx;
                pick.snr = max_cf;
                pick.confidence = std::min(1.0, max_cf / 10.0);
                
                // Compute amplitude
                size_t window = static_cast<size_t>(waveform.sampleRate());
                size_t end_idx = std::min(pick_idx + window, data.size());
                double max_amp = 0;
                for (size_t j = pick_idx; j < end_idx; j++) {
                    max_amp = std::max(max_amp, std::abs(data[j]));
                }
                pick.amplitude = max_amp;
                
                picks.push_back(pick);
                
                triggered = false;
                max_cf = 0;
            }
        }
    }
    
    return picks;
}

SampleVector STALTAPicker::computeSTALTA(const SampleVector& data, 
                                          double sample_rate) const {
    size_t n = data.size();
    size_t sta_samples = static_cast<size_t>(sta_length_ * sample_rate);
    size_t lta_samples = static_cast<size_t>(lta_length_ * sample_rate);
    
    if (sta_samples < 1) sta_samples = 1;
    if (lta_samples <= sta_samples) lta_samples = sta_samples + 1;
    
    SampleVector cf(n, 0.0);
    
    // Use absolute values or squared values for energy
    SampleVector energy(n);
    for (size_t i = 0; i < n; i++) {
        energy[i] = data[i] * data[i];
    }
    
    // Recursive STA/LTA computation
    double sta = 0, lta = 0;
    double sta_coeff = 1.0 / sta_samples;
    double lta_coeff = 1.0 / lta_samples;
    
    // Initialize LTA with first lta_samples
    for (size_t i = 0; i < lta_samples && i < n; i++) {
        lta += energy[i];
    }
    lta /= lta_samples;
    
    // Initialize STA
    size_t sta_start = lta_samples > sta_samples ? lta_samples - sta_samples : 0;
    for (size_t i = sta_start; i < lta_samples && i < n; i++) {
        sta += energy[i];
    }
    sta /= sta_samples;
    
    // Compute ratio for remaining samples
    for (size_t i = lta_samples; i < n; i++) {
        // Update STA and LTA recursively
        sta = sta + sta_coeff * (energy[i] - energy[i - sta_samples]);
        lta = lta + lta_coeff * (energy[i] - energy[i - lta_samples]);
        
        if (lta > 1e-10) {
            cf[i] = sta / lta;
        }
    }
    
    return cf;
}

size_t STALTAPicker::refinePick(const SampleVector& data, size_t approx_idx) const {
    // Use AIC for refinement
    size_t window = 50;  // samples
    
    size_t start = approx_idx > window ? approx_idx - window : 0;
    size_t end = std::min(approx_idx + window, data.size());
    
    if (end - start < 10) return approx_idx;
    
    // Extract window
    SampleVector window_data(data.begin() + start, data.begin() + end);
    size_t n = window_data.size();
    
    // Compute AIC
    double best_aic = std::numeric_limits<double>::max();
    size_t best_idx = window;
    
    for (size_t k = 2; k < n - 2; k++) {
        // Variance before k
        double sum1 = 0, sum1_sq = 0;
        for (size_t i = 0; i < k; i++) {
            sum1 += window_data[i];
            sum1_sq += window_data[i] * window_data[i];
        }
        double var1 = (sum1_sq - sum1*sum1/k) / k;
        
        // Variance after k
        double sum2 = 0, sum2_sq = 0;
        for (size_t i = k; i < n; i++) {
            sum2 += window_data[i];
            sum2_sq += window_data[i] * window_data[i];
        }
        double var2 = (sum2_sq - sum2*sum2/(n-k)) / (n-k);
        
        // AIC
        double aic = k * std::log(std::max(var1, 1e-10)) + 
                     (n - k) * std::log(std::max(var2, 1e-10));
        
        if (aic < best_aic) {
            best_aic = aic;
            best_idx = k;
        }
    }
    
    return start + best_idx;
}

SampleVector STALTAPicker::characteristicFunction(const Waveform& waveform) const {
    SampleVector data = waveform.data();
    double mean = std::accumulate(data.begin(), data.end(), 0.0) / data.size();
    for (auto& s : data) s -= mean;
    
    // Apply same bandpass filter as pick()
    if (use_filter_ && filter_low_ > 0 && filter_high_ > 0 &&
        filter_high_ < waveform.sampleRate() / 2.0 &&
        data.size() > 100) {
        auto bp = IIRFilter::butterworth(2, filter_low_, filter_high_,
                                          waveform.sampleRate());
        data = bp.filtfilt(data);
    }
    
    return computeSTALTA(data, waveform.sampleRate());
}

void STALTAPicker::setParameter(const std::string& name, double value) {
    if (name == "sta_length") sta_length_ = value;
    else if (name == "lta_length") lta_length_ = value;
    else if (name == "trigger_ratio") trigger_ratio_ = value;
    else if (name == "detrigger_ratio") detrigger_ratio_ = value;
    else if (name == "filter_low") filter_low_ = value;
    else if (name == "filter_high") filter_high_ = value;
    else if (name == "use_filter") use_filter_ = (value != 0);
}

double STALTAPicker::getParameter(const std::string& name) const {
    if (name == "sta_length") return sta_length_;
    if (name == "lta_length") return lta_length_;
    if (name == "trigger_ratio") return trigger_ratio_;
    if (name == "detrigger_ratio") return detrigger_ratio_;
    if (name == "filter_low") return filter_low_;
    if (name == "filter_high") return filter_high_;
    if (name == "use_filter") return use_filter_ ? 1.0 : 0.0;
    return 0;
}

// RecursiveSTALTA implementation

RecursiveSTALTA::RecursiveSTALTA(double sta_len, double lta_len, double sample_rate)
    : sta_(0)
    , lta_(0)
    , initialized_(false)
    , sample_count_(0)
{
    sta_coeff_ = 1.0 / (sta_len * sample_rate);
    lta_coeff_ = 1.0 / (lta_len * sample_rate);
    lta_samples_ = static_cast<size_t>(lta_len * sample_rate);
}

double RecursiveSTALTA::process(double sample) {
    double energy = sample * sample;
    
    sample_count_++;
    
    if (!initialized_) {
        lta_ = lta_ + (energy - lta_) / sample_count_;
        
        if (sample_count_ >= lta_samples_) {
            sta_ = lta_;
            initialized_ = true;
        }
        return 0;
    }
    
    // Recursive update
    sta_ = sta_ + sta_coeff_ * (energy - sta_);
    lta_ = lta_ + lta_coeff_ * (energy - lta_);
    
    return lta_ > 0 ? sta_ / lta_ : 0;
}

void RecursiveSTALTA::reset() {
    sta_ = 0;
    lta_ = 0;
    initialized_ = false;
    sample_count_ = 0;
}

// MultiFreqSTALTA implementation

MultiFreqSTALTA::MultiFreqSTALTA()
    : sta_length_(0.5)
    , lta_length_(10.0)
    , trigger_ratio_(3.0)
{
}

void MultiFreqSTALTA::addBand(double low_freq, double high_freq) {
    bands_.push_back({low_freq, high_freq});
}

void MultiFreqSTALTA::clearBands() {
    bands_.clear();
}

std::vector<PickResult> MultiFreqSTALTA::pick(const Waveform& waveform) {
    // If no bands defined, use defaults
    if (bands_.empty()) {
        bands_ = {{1.0, 5.0}, {5.0, 15.0}, {15.0, 30.0}};
    }
    
    // For each band, compute STA/LTA and combine
    STALTAPicker picker;
    picker.setSTALength(sta_length_);
    picker.setLTALength(lta_length_);
    picker.setTriggerRatio(trigger_ratio_);
    
    return picker.pick(waveform);  // Simplified - use single band for now
}

void MultiFreqSTALTA::setParameter(const std::string& name, double value) {
    if (name == "sta_length") sta_length_ = value;
    else if (name == "lta_length") lta_length_ = value;
    else if (name == "trigger_ratio") trigger_ratio_ = value;
}

double MultiFreqSTALTA::getParameter(const std::string& name) const {
    if (name == "sta_length") return sta_length_;
    if (name == "lta_length") return lta_length_;
    if (name == "trigger_ratio") return trigger_ratio_;
    return 0;
}

} // namespace realdetect
