#pragma once

#include "types.hpp"
#include <memory>
#include <algorithm>
#include <numeric>
#include <cmath>

namespace realdetect {

/**
 * Waveform - Container for continuous seismic data
 */
class Waveform {
public:
    Waveform() : sample_rate_(0) {}
    
    Waveform(const StreamID& id, double sample_rate, TimePoint start_time)
        : stream_id_(id), sample_rate_(sample_rate), start_time_(start_time) {}
    
    // Accessors
    const StreamID& streamId() const { return stream_id_; }
    double sampleRate() const { return sample_rate_; }
    TimePoint startTime() const { return start_time_; }
    TimePoint endTime() const {
        if (data_.empty()) return start_time_;
        auto dur = std::chrono::microseconds(
            static_cast<int64_t>((data_.size() - 1) * 1e6 / sample_rate_));
        return start_time_ + dur;
    }
    
    size_t sampleCount() const { return data_.size(); }
    double duration() const { return data_.size() / sample_rate_; }
    
    // Data access
    const SampleVector& data() const { return data_; }
    SampleVector& data() { return data_; }
    
    Sample operator[](size_t idx) const { return data_[idx]; }
    Sample& operator[](size_t idx) { return data_[idx]; }
    
    // Append samples
    void append(Sample s) { data_.push_back(s); }
    void append(const SampleVector& samples) {
        data_.insert(data_.end(), samples.begin(), samples.end());
    }
    void append(const Sample* samples, size_t count) {
        data_.insert(data_.end(), samples, samples + count);
    }
    
    // Get sample at specific time
    Sample sampleAt(TimePoint t) const {
        auto dt = std::chrono::duration_cast<std::chrono::microseconds>(t - start_time_);
        double idx = dt.count() * 1e-6 * sample_rate_;
        size_t i = static_cast<size_t>(idx);
        if (i >= data_.size()) return 0;
        return data_[i];
    }
    
    // Get index for time
    int64_t indexAt(TimePoint t) const {
        auto dt = std::chrono::duration_cast<std::chrono::microseconds>(t - start_time_);
        return static_cast<int64_t>(dt.count() * 1e-6 * sample_rate_);
    }
    
    // Get time for index
    TimePoint timeAt(size_t idx) const {
        auto dur = std::chrono::microseconds(
            static_cast<int64_t>(idx * 1e6 / sample_rate_));
        return start_time_ + dur;
    }
    
    // Statistical operations
    Sample mean() const {
        if (data_.empty()) return 0;
        return std::accumulate(data_.begin(), data_.end(), 0.0) / data_.size();
    }
    
    Sample stddev() const {
        if (data_.size() < 2) return 0;
        Sample m = mean();
        Sample sum = 0;
        for (auto s : data_) sum += (s - m) * (s - m);
        return std::sqrt(sum / (data_.size() - 1));
    }
    
    Sample max() const {
        if (data_.empty()) return 0;
        return *std::max_element(data_.begin(), data_.end());
    }
    
    Sample min() const {
        if (data_.empty()) return 0;
        return *std::min_element(data_.begin(), data_.end());
    }
    
    Sample absMax() const {
        if (data_.empty()) return 0;
        Sample m = 0;
        for (auto s : data_) {
            Sample a = std::abs(s);
            if (a > m) m = a;
        }
        return m;
    }
    
    // Processing operations
    void demean() {
        Sample m = mean();
        for (auto& s : data_) s -= m;
    }
    
    void detrend() {
        if (data_.size() < 2) return;
        size_t n = data_.size();
        double sx = 0, sy = 0, sxy = 0, sxx = 0;
        for (size_t i = 0; i < n; i++) {
            sx += i;
            sy += data_[i];
            sxy += i * data_[i];
            sxx += i * i;
        }
        double slope = (n * sxy - sx * sy) / (n * sxx - sx * sx);
        double intercept = (sy - slope * sx) / n;
        for (size_t i = 0; i < n; i++) {
            data_[i] -= (slope * i + intercept);
        }
    }
    
    void normalize() {
        Sample m = absMax();
        if (m > 0) {
            for (auto& s : data_) s /= m;
        }
    }
    
    void taper(double fraction = 0.05) {
        if (data_.empty()) return;
        size_t taper_len = static_cast<size_t>(data_.size() * fraction);
        if (taper_len < 1) taper_len = 1;
        
        for (size_t i = 0; i < taper_len; i++) {
            double w = 0.5 * (1.0 - std::cos(M_PI * i / taper_len));
            data_[i] *= w;
            data_[data_.size() - 1 - i] *= w;
        }
    }
    
    // Slice waveform
    Waveform slice(size_t start_idx, size_t end_idx) const {
        Waveform result(stream_id_, sample_rate_, timeAt(start_idx));
        if (end_idx > data_.size()) end_idx = data_.size();
        if (start_idx < end_idx) {
            result.data_.assign(data_.begin() + start_idx, data_.begin() + end_idx);
        }
        return result;
    }
    
    Waveform slice(TimePoint start, TimePoint end) const {
        int64_t s = indexAt(start);
        int64_t e = indexAt(end);
        if (s < 0) s = 0;
        if (e < 0) e = 0;
        return slice(static_cast<size_t>(s), static_cast<size_t>(e));
    }
    
    // Resample to new rate
    Waveform resample(double new_rate) const;
    
    // Integrate
    void integrate() {
        if (data_.size() < 2) return;
        double dt = 1.0 / sample_rate_;
        for (size_t i = 1; i < data_.size(); i++) {
            data_[i] = data_[i-1] + data_[i] * dt;
        }
    }
    
    // Differentiate
    void differentiate() {
        if (data_.size() < 2) return;
        double dt = 1.0 / sample_rate_;
        for (size_t i = data_.size() - 1; i > 0; i--) {
            data_[i] = (data_[i] - data_[i-1]) / dt;
        }
        data_[0] = data_[1];
    }

private:
    StreamID stream_id_;
    double sample_rate_;
    TimePoint start_time_;
    SampleVector data_;
};

using WaveformPtr = std::shared_ptr<Waveform>;

} // namespace realdetect
