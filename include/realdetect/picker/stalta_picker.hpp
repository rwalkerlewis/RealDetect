#pragma once

#include "picker.hpp"
#include "filter_bank.hpp"

namespace realdetect {

/**
 * STALTAPicker - Classic Short-Term Average / Long-Term Average detector
 * 
 * The STA/LTA ratio is computed as:
 *   R(i) = STA(i) / LTA(i)
 * 
 * Where STA and LTA are computed using recursive formulas for efficiency.
 * Triggers when R exceeds threshold, detriggers when R falls below detrigger.
 */
class STALTAPicker : public BasePicker {
public:
    STALTAPicker();
    
    std::vector<PickResult> pick(const Waveform& waveform) override;
    std::string name() const override { return "STA/LTA"; }
    
    void setParameter(const std::string& name, double value) override;
    double getParameter(const std::string& name) const override;
    
    // Specific parameters
    void setSTALength(double seconds) { sta_length_ = seconds; }
    void setLTALength(double seconds) { lta_length_ = seconds; }
    void setTriggerRatio(double ratio) { trigger_ratio_ = ratio; }
    void setDetriggerRatio(double ratio) { detrigger_ratio_ = ratio; }
    
    double staLength() const { return sta_length_; }
    double ltaLength() const { return lta_length_; }
    double triggerRatio() const { return trigger_ratio_; }
    double detriggerRatio() const { return detrigger_ratio_; }
    
    // Get characteristic function (STA/LTA ratio)
    SampleVector characteristicFunction(const Waveform& waveform) const;

private:
    double sta_length_;       // STA window (seconds)
    double lta_length_;       // LTA window (seconds)
    double trigger_ratio_;    // Trigger threshold
    double detrigger_ratio_;  // Detrigger threshold
    
    // Optional bandpass filter
    double filter_low_;
    double filter_high_;
    bool use_filter_;
    
    // Compute STA/LTA
    SampleVector computeSTALTA(const SampleVector& data, double sample_rate) const;
    
    // Refine pick time using AIC
    size_t refinePick(const SampleVector& data, size_t approx_idx) const;
};

/**
 * RecursiveSTALTA - Optimized recursive implementation
 */
class RecursiveSTALTA {
public:
    RecursiveSTALTA(double sta_len, double lta_len, double sample_rate);
    
    // Process single sample, returns STA/LTA ratio
    double process(double sample);
    
    // Reset state
    void reset();
    
    // Current values
    double sta() const { return sta_; }
    double lta() const { return lta_; }
    double ratio() const { return lta_ > 0 ? sta_ / lta_ : 0; }

private:
    double sta_coeff_;
    double lta_coeff_;
    double sta_;
    double lta_;
    bool initialized_;
    size_t sample_count_;
    size_t lta_samples_;
};

/**
 * MultiFreqSTALTA - STA/LTA on multiple frequency bands
 */
class MultiFreqSTALTA : public BasePicker {
public:
    MultiFreqSTALTA();
    
    std::vector<PickResult> pick(const Waveform& waveform) override;
    std::string name() const override { return "MultiFreq-STA/LTA"; }
    
    void setParameter(const std::string& name, double value) override;
    double getParameter(const std::string& name) const override;
    
    // Add frequency band
    void addBand(double low_freq, double high_freq);
    void clearBands();

private:
    std::vector<std::pair<double, double>> bands_;
    double sta_length_;
    double lta_length_;
    double trigger_ratio_;
    FilterBank filter_bank_;
};

} // namespace realdetect
