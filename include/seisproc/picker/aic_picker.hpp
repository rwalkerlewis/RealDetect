#pragma once

#include "picker.hpp"

namespace seisproc {

/**
 * AICPicker - Akaike Information Criterion picker
 * 
 * The AIC function is computed as:
 *   AIC(k) = k * log(var(x[0:k])) + (N-k-1) * log(var(x[k:N]))
 * 
 * The pick time is at the minimum of the AIC function.
 * This is often used to refine picks from STA/LTA.
 */
class AICPicker : public BasePicker {
public:
    AICPicker();
    
    std::vector<PickResult> pick(const Waveform& waveform) override;
    std::string name() const override { return "AIC"; }
    
    void setParameter(const std::string& name, double value) override;
    double getParameter(const std::string& name) const override;
    
    // Refine a pick using AIC (returns sample index)
    size_t refinePick(const Waveform& waveform, size_t approx_idx, 
                      size_t window_size = 100);
    
    // Compute AIC function
    static SampleVector computeAIC(const SampleVector& data);
    
    // Find minimum in AIC (the pick)
    static size_t findAICMinimum(const SampleVector& aic, 
                                  size_t start = 0, size_t end = 0);

private:
    double pre_window_;   // Window before trigger (seconds)
    double post_window_;  // Window after trigger (seconds)
    double min_snr_;      // Minimum SNR
    
    // Pre-detection using energy change
    std::vector<size_t> findCandidates(const Waveform& waveform) const;
};

/**
 * AR-AIC Picker - Autoregressive AIC picker
 * 
 * Uses autoregressive prediction error as the characteristic function.
 * More sensitive to subtle phase arrivals.
 */
class ARPicker : public BasePicker {
public:
    ARPicker();
    
    std::vector<PickResult> pick(const Waveform& waveform) override;
    std::string name() const override { return "AR-AIC"; }
    
    void setParameter(const std::string& name, double value) override;
    double getParameter(const std::string& name) const override;
    
    // Set AR model order
    void setAROrder(int order) { ar_order_ = order; }
    int arOrder() const { return ar_order_; }

private:
    int ar_order_;        // AR model order (typically 2-10)
    double pre_window_;   // Analysis window before pick
    double post_window_;  // Analysis window after pick
    
    // Compute AR coefficients using Burg's method
    std::vector<double> burgAR(const SampleVector& data) const;
    
    // Compute prediction error
    SampleVector computePredictionError(const SampleVector& data) const;
};

} // namespace seisproc
