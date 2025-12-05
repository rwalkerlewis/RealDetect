#include "realdetect/core/waveform.hpp"
#include <cmath>
#include <algorithm>

namespace realdetect {

Waveform Waveform::resample(double new_rate) const {
    if (data_.empty() || new_rate <= 0 || new_rate == sample_rate_) {
        Waveform result = *this;
        return result;
    }
    
    Waveform result(stream_id_, new_rate, start_time_);
    
    double ratio = sample_rate_ / new_rate;
    size_t new_size = static_cast<size_t>(data_.size() / ratio);
    result.data_.resize(new_size);
    
    // Linear interpolation resampling
    for (size_t i = 0; i < new_size; i++) {
        double src_idx = i * ratio;
        size_t idx0 = static_cast<size_t>(src_idx);
        size_t idx1 = idx0 + 1;
        double frac = src_idx - idx0;
        
        if (idx1 >= data_.size()) {
            result.data_[i] = data_.back();
        } else {
            result.data_[i] = data_[idx0] * (1.0 - frac) + data_[idx1] * frac;
        }
    }
    
    return result;
}

} // namespace realdetect
