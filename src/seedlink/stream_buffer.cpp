#include "realdetect/seedlink/seedlink_client.hpp"
#include <algorithm>

namespace realdetect {

StreamBuffer::StreamBuffer(size_t max_seconds, double sample_rate)
    : sample_rate_(sample_rate)
    , max_samples_(static_cast<size_t>(max_seconds * sample_rate))
{
}

void StreamBuffer::append(const SeedLinkPacket& packet) {
    std::lock_guard<std::mutex> lock(mutex_);
    
    if (buffer_.empty()) {
        stream_id_ = packet.stream_id;
        sample_rate_ = packet.sample_rate;
        start_time_ = packet.start_time;
    }
    
    buffer_.insert(buffer_.end(), packet.samples.begin(), packet.samples.end());
    end_time_ = packet.start_time;
    
    // Update end time based on samples
    if (!packet.samples.empty()) {
        auto dur = std::chrono::microseconds(
            static_cast<int64_t>((packet.samples.size() - 1) * 1e6 / sample_rate_));
        end_time_ = packet.start_time + dur;
    }
    
    // Trim if buffer is too large
    if (buffer_.size() > max_samples_) {
        size_t trim = buffer_.size() - max_samples_;
        buffer_.erase(buffer_.begin(), buffer_.begin() + trim);
        
        auto dur = std::chrono::microseconds(
            static_cast<int64_t>(trim * 1e6 / sample_rate_));
        start_time_ += dur;
    }
}

void StreamBuffer::append(WaveformPtr waveform) {
    std::lock_guard<std::mutex> lock(mutex_);
    
    if (buffer_.empty()) {
        stream_id_ = waveform->streamId();
        sample_rate_ = waveform->sampleRate();
        start_time_ = waveform->startTime();
    }
    
    const auto& data = waveform->data();
    buffer_.insert(buffer_.end(), data.begin(), data.end());
    end_time_ = waveform->endTime();
    
    // Trim if needed
    if (buffer_.size() > max_samples_) {
        size_t trim = buffer_.size() - max_samples_;
        buffer_.erase(buffer_.begin(), buffer_.begin() + trim);
        
        auto dur = std::chrono::microseconds(
            static_cast<int64_t>(trim * 1e6 / sample_rate_));
        start_time_ += dur;
    }
}

WaveformPtr StreamBuffer::getWaveform(TimePoint start, TimePoint end) const {
    std::lock_guard<std::mutex> lock(mutex_);
    
    auto wf = std::make_shared<Waveform>(stream_id_, sample_rate_, start);
    
    if (buffer_.empty()) return wf;
    
    // Calculate indices
    auto start_offset = std::chrono::duration_cast<std::chrono::microseconds>(
        start - start_time_);
    auto end_offset = std::chrono::duration_cast<std::chrono::microseconds>(
        end - start_time_);
    
    int64_t start_idx = static_cast<int64_t>(start_offset.count() * 1e-6 * sample_rate_);
    int64_t end_idx = static_cast<int64_t>(end_offset.count() * 1e-6 * sample_rate_);
    
    if (start_idx < 0) start_idx = 0;
    if (end_idx > static_cast<int64_t>(buffer_.size())) 
        end_idx = buffer_.size();
    
    if (start_idx < end_idx) {
        wf->data().assign(buffer_.begin() + start_idx, buffer_.begin() + end_idx);
    }
    
    return wf;
}

WaveformPtr StreamBuffer::getLatest(double seconds) const {
    std::lock_guard<std::mutex> lock(mutex_);
    
    size_t samples = static_cast<size_t>(seconds * sample_rate_);
    if (samples > buffer_.size()) samples = buffer_.size();
    
    auto dur = std::chrono::microseconds(
        static_cast<int64_t>((buffer_.size() - samples) * 1e6 / sample_rate_));
    TimePoint start = start_time_ + dur;
    
    auto wf = std::make_shared<Waveform>(stream_id_, sample_rate_, start);
    
    if (samples > 0) {
        wf->data().assign(buffer_.end() - samples, buffer_.end());
    }
    
    return wf;
}

void StreamBuffer::clear() {
    std::lock_guard<std::mutex> lock(mutex_);
    buffer_.clear();
}

void StreamBuffer::trimTo(TimePoint earliest) {
    std::lock_guard<std::mutex> lock(mutex_);
    
    if (buffer_.empty()) return;
    
    auto offset = std::chrono::duration_cast<std::chrono::microseconds>(
        earliest - start_time_);
    int64_t idx = static_cast<int64_t>(offset.count() * 1e-6 * sample_rate_);
    
    if (idx > 0 && idx < static_cast<int64_t>(buffer_.size())) {
        buffer_.erase(buffer_.begin(), buffer_.begin() + idx);
        start_time_ = earliest;
    }
}

// MultiStreamBuffer implementation

MultiStreamBuffer::MultiStreamBuffer(size_t max_seconds)
    : max_seconds_(max_seconds)
{
}

void MultiStreamBuffer::addPacket(const SeedLinkPacket& packet) {
    std::lock_guard<std::mutex> lock(mutex_);
    
    std::string key = packet.stream_id.toString();
    
    auto it = buffers_.find(key);
    if (it == buffers_.end()) {
        buffers_[key] = std::make_unique<StreamBuffer>(max_seconds_, packet.sample_rate);
    }
    
    buffers_[key]->append(packet);
}

void MultiStreamBuffer::addWaveform(WaveformPtr waveform) {
    std::lock_guard<std::mutex> lock(mutex_);
    
    std::string key = waveform->streamId().toString();
    
    auto it = buffers_.find(key);
    if (it == buffers_.end()) {
        buffers_[key] = std::make_unique<StreamBuffer>(
            max_seconds_, waveform->sampleRate());
    }
    
    buffers_[key]->append(waveform);
}

StreamBuffer* MultiStreamBuffer::getBuffer(const StreamID& id) {
    std::lock_guard<std::mutex> lock(mutex_);
    
    auto it = buffers_.find(id.toString());
    return it != buffers_.end() ? it->second.get() : nullptr;
}

const StreamBuffer* MultiStreamBuffer::getBuffer(const StreamID& id) const {
    std::lock_guard<std::mutex> lock(mutex_);
    
    auto it = buffers_.find(id.toString());
    return it != buffers_.end() ? it->second.get() : nullptr;
}

std::vector<StreamID> MultiStreamBuffer::streams() const {
    std::lock_guard<std::mutex> lock(mutex_);
    
    std::vector<StreamID> result;
    for (const auto& [key, buf] : buffers_) {
        result.push_back(buf->streamId());
    }
    return result;
}

void MultiStreamBuffer::clear() {
    std::lock_guard<std::mutex> lock(mutex_);
    buffers_.clear();
}

void MultiStreamBuffer::trimAll(TimePoint earliest) {
    std::lock_guard<std::mutex> lock(mutex_);
    
    for (auto& [key, buf] : buffers_) {
        buf->trimTo(earliest);
    }
}

} // namespace realdetect
