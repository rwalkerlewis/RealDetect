#include "seisproc/core/miniseed.hpp"
#include <fstream>
#include <cstring>
#include <iostream>
#include <algorithm>
#include <map>

namespace seisproc {

TimePoint MiniSeedRecord::parseTime(const uint8_t* btime) {
    // BTIME structure (10 bytes):
    // year (2), day of year (2), hour (1), min (1), sec (1), unused (1), 0.0001s (2)
    uint16_t year = (btime[0] << 8) | btime[1];
    uint16_t doy = (btime[2] << 8) | btime[3];
    uint8_t hour = btime[4];
    uint8_t min = btime[5];
    uint8_t sec = btime[6];
    uint16_t frac = (btime[8] << 8) | btime[9];
    
    // Convert to time_t
    std::tm tm = {};
    tm.tm_year = year - 1900;
    tm.tm_mday = 1;  // Day of year handling
    tm.tm_yday = doy - 1;
    tm.tm_hour = hour;
    tm.tm_min = min;
    tm.tm_sec = sec;
    
    // Calculate month and day from day of year
    int days_in_month[] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
    if ((year % 4 == 0 && year % 100 != 0) || year % 400 == 0) {
        days_in_month[1] = 29;  // Leap year
    }
    
    int remaining = doy;
    for (int m = 0; m < 12; m++) {
        if (remaining <= days_in_month[m]) {
            tm.tm_mon = m;
            tm.tm_mday = remaining;
            break;
        }
        remaining -= days_in_month[m];
    }
    
    time_t time = timegm(&tm);
    auto tp = std::chrono::system_clock::from_time_t(time);
    tp += std::chrono::microseconds(frac * 100);
    
    return tp;
}

bool MiniSeedRecord::parse(const uint8_t* data, size_t length) {
    if (length < HEADER_SIZE) return false;
    
    // Check for valid SEED sequence number
    for (int i = 0; i < 6; i++) {
        if (!std::isdigit(data[i]) && data[i] != ' ') return false;
    }
    
    // Data quality indicator (D, R, Q, M)
    char quality = data[6];
    if (quality != 'D' && quality != 'R' && quality != 'Q' && quality != 'M') {
        return false;
    }
    
    // Parse fixed header (48 bytes)
    // Station identifier
    stream_id_.station = std::string(reinterpret_cast<const char*>(data + 8), 5);
    stream_id_.station.erase(stream_id_.station.find_last_not_of(' ') + 1);
    
    stream_id_.location = std::string(reinterpret_cast<const char*>(data + 13), 2);
    stream_id_.channel = std::string(reinterpret_cast<const char*>(data + 15), 3);
    stream_id_.network = std::string(reinterpret_cast<const char*>(data + 18), 2);
    stream_id_.network.erase(stream_id_.network.find_last_not_of(' ') + 1);
    
    // Start time
    start_time_ = parseTime(data + 20);
    
    // Number of samples
    sample_count_ = (data[30] << 8) | data[31];
    
    // Sample rate factor and multiplier
    int16_t factor = (data[32] << 8) | data[33];
    int16_t multiplier = (data[34] << 8) | data[35];
    
    if (factor > 0 && multiplier > 0) {
        sample_rate_ = factor * multiplier;
    } else if (factor > 0 && multiplier < 0) {
        sample_rate_ = -factor / multiplier;
    } else if (factor < 0 && multiplier > 0) {
        sample_rate_ = -multiplier / factor;
    } else {
        sample_rate_ = multiplier / factor;
    }
    
    // Activity, I/O, data quality flags (bytes 36-38)
    
    // Number of blockettes
    uint8_t num_blockettes = data[39];
    
    // Time correction (not applied to start time)
    
    // Beginning of data
    uint16_t data_offset = (data[44] << 8) | data[45];
    
    // Beginning of blockettes
    uint16_t blockette_offset = (data[46] << 8) | data[47];
    
    // Parse blockettes to find encoding
    encoding_ = Encoding::STEIM2;  // Default
    
    if (num_blockettes > 0 && blockette_offset >= 48 && blockette_offset < length) {
        const uint8_t* bp = data + blockette_offset;
        
        while (bp < data + length - 4) {
            uint16_t type = (bp[0] << 8) | bp[1];
            uint16_t next = (bp[2] << 8) | bp[3];
            
            if (type == 1000) {  // Data Only SEED Blockette
                encoding_ = static_cast<Encoding>(bp[4]);
                record_length_ = 1 << bp[6];
                break;
            }
            
            if (next == 0 || next <= blockette_offset) break;
            bp = data + next;
        }
    }
    
    // Decode data
    if (data_offset > 0 && data_offset < length) {
        const uint8_t* dp = data + data_offset;
        size_t data_len = length - data_offset;
        
        switch (encoding_) {
            case Encoding::STEIM1:
            case Encoding::STEIM2:
                return decodeSteimData(dp, data_len, encoding_ == Encoding::STEIM1 ? 1 : 2);
            case Encoding::INT16:
                return decodeInt16Data(dp, data_len);
            case Encoding::INT32:
                return decodeInt32Data(dp, data_len);
            case Encoding::FLOAT32:
                return decodeFloat32Data(dp, data_len);
            default:
                std::cerr << "Unsupported encoding: " << static_cast<int>(encoding_) << std::endl;
                return false;
        }
    }
    
    return true;
}

bool MiniSeedRecord::decodeSteimData(const uint8_t* data, size_t length, int steim_level) {
    samples_.clear();
    samples_.reserve(sample_count_);
    
    if (length < 64) return false;
    
    // Each Steim frame is 64 bytes (16 32-bit words)
    // First word is control word, next 15 are data
    size_t num_frames = length / 64;
    
    int32_t x0 = 0;  // First sample (forward integration constant)
    int32_t xn = 0;  // Last sample (reverse integration constant)
    int32_t x = 0;   // Current sample
    
    for (size_t frame = 0; frame < num_frames && samples_.size() < sample_count_; frame++) {
        const uint8_t* fp = data + frame * 64;
        
        // Control word
        uint32_t ctrl = (fp[0] << 24) | (fp[1] << 16) | (fp[2] << 8) | fp[3];
        
        for (int word = 1; word < 16 && samples_.size() < sample_count_; word++) {
            const uint8_t* wp = fp + word * 4;
            uint32_t w = (wp[0] << 24) | (wp[1] << 16) | (wp[2] << 8) | wp[3];
            
            // Get nibble for this word
            int nibble = (ctrl >> (30 - (word - 1) * 2)) & 0x03;
            
            if (frame == 0 && word == 1) {
                // First data word is x0
                x0 = static_cast<int32_t>(w);
                x = x0;
                samples_.push_back(x);
                continue;
            }
            if (frame == 0 && word == 2) {
                // Second data word is xn (not used for decompression)
                xn = static_cast<int32_t>(w);
                continue;
            }
            
            if (nibble == 0) {
                // Non-data (special)
                continue;
            } else if (nibble == 1) {
                // Four 8-bit differences
                for (int i = 0; i < 4 && samples_.size() < sample_count_; i++) {
                    int8_t d = (w >> (24 - i * 8)) & 0xFF;
                    x += d;
                    samples_.push_back(x);
                }
            } else if (nibble == 2) {
                if (steim_level == 1) {
                    // Two 16-bit differences
                    for (int i = 0; i < 2 && samples_.size() < sample_count_; i++) {
                        int16_t d = (w >> (16 - i * 16)) & 0xFFFF;
                        x += d;
                        samples_.push_back(x);
                    }
                } else {
                    // Steim2: dnib determines format
                    int dnib = (w >> 30) & 0x03;
                    if (dnib == 1) {
                        // One 30-bit difference
                        int32_t d = (w & 0x3FFFFFFF);
                        if (d & 0x20000000) d |= 0xC0000000;  // Sign extend
                        x += d;
                        samples_.push_back(x);
                    } else if (dnib == 2) {
                        // Two 15-bit differences
                        for (int i = 0; i < 2 && samples_.size() < sample_count_; i++) {
                            int16_t d = (w >> (15 - i * 15)) & 0x7FFF;
                            if (d & 0x4000) d |= 0x8000;
                            x += d;
                            samples_.push_back(x);
                        }
                    } else if (dnib == 3) {
                        // Three 10-bit differences
                        for (int i = 0; i < 3 && samples_.size() < sample_count_; i++) {
                            int16_t d = (w >> (20 - i * 10)) & 0x3FF;
                            if (d & 0x200) d |= 0xFC00;
                            x += d;
                            samples_.push_back(x);
                        }
                    }
                }
            } else if (nibble == 3) {
                if (steim_level == 1) {
                    // One 32-bit difference
                    int32_t d = static_cast<int32_t>(w);
                    x += d;
                    samples_.push_back(x);
                } else {
                    // Steim2: dnib determines format
                    int dnib = (w >> 30) & 0x03;
                    if (dnib == 0) {
                        // Five 6-bit differences
                        for (int i = 0; i < 5 && samples_.size() < sample_count_; i++) {
                            int8_t d = (w >> (24 - i * 6)) & 0x3F;
                            if (d & 0x20) d |= 0xC0;
                            x += d;
                            samples_.push_back(x);
                        }
                    } else if (dnib == 1) {
                        // Six 5-bit differences
                        for (int i = 0; i < 6 && samples_.size() < sample_count_; i++) {
                            int8_t d = (w >> (25 - i * 5)) & 0x1F;
                            if (d & 0x10) d |= 0xE0;
                            x += d;
                            samples_.push_back(x);
                        }
                    } else if (dnib == 2) {
                        // Seven 4-bit differences
                        for (int i = 0; i < 7 && samples_.size() < sample_count_; i++) {
                            int8_t d = (w >> (24 - i * 4)) & 0x0F;
                            if (d & 0x08) d |= 0xF0;
                            x += d;
                            samples_.push_back(x);
                        }
                    }
                }
            }
        }
    }
    
    return samples_.size() == sample_count_;
}

bool MiniSeedRecord::decodeInt16Data(const uint8_t* data, size_t length) {
    samples_.clear();
    size_t count = std::min(sample_count_, length / 2);
    samples_.reserve(count);
    
    for (size_t i = 0; i < count; i++) {
        int16_t val = (data[i*2] << 8) | data[i*2 + 1];
        samples_.push_back(val);
    }
    return true;
}

bool MiniSeedRecord::decodeInt32Data(const uint8_t* data, size_t length) {
    samples_.clear();
    size_t count = std::min(sample_count_, length / 4);
    samples_.reserve(count);
    
    for (size_t i = 0; i < count; i++) {
        int32_t val = (data[i*4] << 24) | (data[i*4+1] << 16) | 
                      (data[i*4+2] << 8) | data[i*4+3];
        samples_.push_back(val);
    }
    return true;
}

bool MiniSeedRecord::decodeFloat32Data(const uint8_t* data, size_t length) {
    samples_.clear();
    size_t count = std::min(sample_count_, length / 4);
    samples_.reserve(count);
    
    for (size_t i = 0; i < count; i++) {
        uint32_t bits = (data[i*4] << 24) | (data[i*4+1] << 16) | 
                        (data[i*4+2] << 8) | data[i*4+3];
        float val;
        std::memcpy(&val, &bits, sizeof(float));
        samples_.push_back(val);
    }
    return true;
}

// MiniSeedReader implementation
bool MiniSeedReader::open(const std::string& filename) {
    std::ifstream file(filename, std::ios::binary);
    if (!file.is_open()) return false;
    
    // Get file size
    file.seekg(0, std::ios::end);
    size_t size = file.tellg();
    file.seekg(0, std::ios::beg);
    
    std::vector<uint8_t> buffer(size);
    file.read(reinterpret_cast<char*>(buffer.data()), size);
    
    return parse(buffer.data(), buffer.size());
}

bool MiniSeedReader::parse(const uint8_t* data, size_t length) {
    records_.clear();
    
    size_t offset = 0;
    while (offset + 48 < length) {
        MiniSeedRecord record;
        
        // Try different record lengths
        for (size_t rec_len : {512, 4096, 256, 1024, 2048, 8192}) {
            if (offset + rec_len > length) continue;
            
            if (record.parse(data + offset, rec_len)) {
                records_.push_back(record);
                offset += rec_len;
                break;
            }
        }
        
        // If parsing failed, try to find next record
        if (records_.empty() || offset == 0) {
            offset += 256;  // Skip ahead
        }
    }
    
    return !records_.empty();
}

std::vector<WaveformPtr> MiniSeedReader::toWaveforms() const {
    std::map<std::string, WaveformPtr> waveforms;
    
    for (const auto& rec : records_) {
        std::string key = rec.streamId().toString();
        
        auto it = waveforms.find(key);
        if (it == waveforms.end()) {
            waveforms[key] = rec.toWaveform();
        } else {
            // Merge if continuous
            it->second->append(rec.samples());
        }
    }
    
    std::vector<WaveformPtr> result;
    for (auto& [key, wf] : waveforms) {
        result.push_back(wf);
    }
    return result;
}

} // namespace seisproc
