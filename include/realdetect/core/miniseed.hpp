#pragma once

#include "types.hpp"
#include "waveform.hpp"
#include <cstring>
#include <cstdint>

namespace realdetect {

/**
 * MiniSEED Record - 512 or 4096 byte seismic data record
 */
class MiniSeedRecord {
public:
    static constexpr size_t HEADER_SIZE = 48;
    static constexpr size_t FIXED_HEADER_SIZE = 48;
    
    // Data encoding types
    enum class Encoding : uint8_t {
        ASCII = 0,
        INT16 = 1,
        INT24 = 2,
        INT32 = 3,
        FLOAT32 = 4,
        FLOAT64 = 5,
        STEIM1 = 10,
        STEIM2 = 11
    };
    
    MiniSeedRecord() : sample_rate_(0), sample_count_(0), 
                       record_length_(512), encoding_(Encoding::STEIM2) {}
    
    // Parse from raw data
    bool parse(const uint8_t* data, size_t length);
    
    // Serialize to raw data
    size_t serialize(uint8_t* buffer, size_t max_length) const;
    
    // Accessors
    const StreamID& streamId() const { return stream_id_; }
    TimePoint startTime() const { return start_time_; }
    double sampleRate() const { return sample_rate_; }
    size_t sampleCount() const { return sample_count_; }
    size_t recordLength() const { return record_length_; }
    Encoding encoding() const { return encoding_; }
    const SampleVector& samples() const { return samples_; }
    
    // Setters
    void setStreamId(const StreamID& id) { stream_id_ = id; }
    void setStartTime(TimePoint t) { start_time_ = t; }
    void setSampleRate(double rate) { sample_rate_ = rate; }
    void setSamples(const SampleVector& s) { 
        samples_ = s; 
        sample_count_ = s.size();
    }
    void setRecordLength(size_t len) { record_length_ = len; }
    void setEncoding(Encoding e) { encoding_ = e; }
    
    // Create waveform from record
    WaveformPtr toWaveform() const {
        auto wf = std::make_shared<Waveform>(stream_id_, sample_rate_, start_time_);
        wf->data() = samples_;
        return wf;
    }

private:
    StreamID stream_id_;
    TimePoint start_time_;
    double sample_rate_;
    size_t sample_count_;
    size_t record_length_;
    Encoding encoding_;
    SampleVector samples_;
    
    // Steim decompression
    bool decodeSteimData(const uint8_t* data, size_t length, int steim_level);
    bool decodeInt16Data(const uint8_t* data, size_t length);
    bool decodeInt32Data(const uint8_t* data, size_t length);
    bool decodeFloat32Data(const uint8_t* data, size_t length);
    
    // Time parsing (BTIME structure)
    static TimePoint parseTime(const uint8_t* btime);
};

/**
 * MiniSeedReader - Read MiniSEED files or streams
 */
class MiniSeedReader {
public:
    MiniSeedReader() = default;
    
    // Read from file
    bool open(const std::string& filename);
    
    // Read from memory
    bool parse(const uint8_t* data, size_t length);
    
    // Get all records
    const std::vector<MiniSeedRecord>& records() const { return records_; }
    
    // Convert to waveforms (merging continuous segments)
    std::vector<WaveformPtr> toWaveforms() const;

private:
    std::vector<MiniSeedRecord> records_;
};

} // namespace realdetect
