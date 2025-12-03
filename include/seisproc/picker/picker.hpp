#pragma once

#include "../core/types.hpp"
#include "../core/waveform.hpp"
#include "../core/event.hpp"
#include <vector>
#include <memory>
#include <functional>

namespace seisproc {

/**
 * PickResult - Result from a phase picker
 */
struct PickResult {
    TimePoint time;
    PhaseType phase_type;
    double snr;
    double confidence;
    double amplitude;
    double period;
    size_t sample_index;
    
    PickResult() : phase_type(PhaseType::Unknown), snr(0), 
                   confidence(0), amplitude(0), period(0), sample_index(0) {}
};

/**
 * BasePicker - Abstract base class for phase pickers
 */
class BasePicker {
public:
    virtual ~BasePicker() = default;
    
    // Pick phases from waveform
    virtual std::vector<PickResult> pick(const Waveform& waveform) = 0;
    
    // Get picker name
    virtual std::string name() const = 0;
    
    // Configuration
    virtual void setParameter(const std::string& name, double value) = 0;
    virtual double getParameter(const std::string& name) const = 0;
};

using PickerPtr = std::shared_ptr<BasePicker>;

/**
 * Callback for real-time picks
 */
using PickCallback = std::function<void(const StreamID&, const PickResult&)>;

/**
 * RealTimePicker - Continuous phase picking on streaming data
 */
class RealTimePicker {
public:
    RealTimePicker();
    
    // Set picker algorithm
    void setPicker(PickerPtr picker) { picker_ = picker; }
    
    // Process new data
    void process(const StreamID& stream_id, const Waveform& waveform);
    
    // Set callback for new picks
    void setPickCallback(PickCallback callback) { pick_callback_ = callback; }
    
    // Configuration
    void setMinimumGap(double seconds) { min_gap_ = seconds; }
    void setMinimumSNR(double snr) { min_snr_ = snr; }

private:
    PickerPtr picker_;
    PickCallback pick_callback_;
    double min_gap_;          // Minimum gap between picks on same channel
    double min_snr_;          // Minimum SNR for valid pick
    
    // Track last pick time per stream
    std::map<std::string, TimePoint> last_pick_time_;
};

} // namespace seisproc
