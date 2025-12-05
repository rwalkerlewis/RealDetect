/**
 * Unit tests for MiniSEED parsing
 */

#include "test_framework.hpp"
#include "realdetect/core/miniseed.hpp"
#include "realdetect/core/waveform.hpp"
#include <cmath>
#include <fstream>

using namespace realdetect;
using namespace realdetect::test;

// ============================================================================
// MiniSeedRecord Tests
// ============================================================================

TEST(MiniSeedRecord, DefaultConstructor) {
    MiniSeedRecord record;
    
    ASSERT_EQ(record.recordLength(), 512u);
    ASSERT_EQ(record.sampleCount(), 0u);
    ASSERT_NEAR(record.sampleRate(), 0.0, 1e-10);
}

TEST(MiniSeedRecord, SetStreamId) {
    MiniSeedRecord record;
    
    StreamID id("IU", "ANMO", "00", "BHZ");
    record.setStreamId(id);
    
    ASSERT_EQ(record.streamId().network, "IU");
    ASSERT_EQ(record.streamId().station, "ANMO");
    ASSERT_EQ(record.streamId().location, "00");
    ASSERT_EQ(record.streamId().channel, "BHZ");
}

TEST(MiniSeedRecord, SetStartTime) {
    MiniSeedRecord record;
    
    auto now = std::chrono::system_clock::now();
    record.setStartTime(now);
    
    auto start = record.startTime();
    auto diff = std::chrono::duration_cast<std::chrono::microseconds>(now - start);
    
    ASSERT_LT(std::abs(diff.count()), 1000);  // Within 1ms
}

TEST(MiniSeedRecord, SetSampleRate) {
    MiniSeedRecord record;
    
    record.setSampleRate(100.0);
    ASSERT_NEAR(record.sampleRate(), 100.0, 1e-10);
    
    record.setSampleRate(40.0);
    ASSERT_NEAR(record.sampleRate(), 40.0, 1e-10);
}

TEST(MiniSeedRecord, SetSamples) {
    MiniSeedRecord record;
    
    SampleVector samples = {1.0, 2.0, 3.0, 4.0, 5.0};
    record.setSamples(samples);
    
    ASSERT_EQ(record.sampleCount(), 5u);
    ASSERT_EQ(record.samples().size(), 5u);
    ASSERT_NEAR(record.samples()[0], 1.0, 1e-10);
    ASSERT_NEAR(record.samples()[4], 5.0, 1e-10);
}

TEST(MiniSeedRecord, SetRecordLength) {
    MiniSeedRecord record;
    
    record.setRecordLength(4096);
    ASSERT_EQ(record.recordLength(), 4096u);
    
    record.setRecordLength(512);
    ASSERT_EQ(record.recordLength(), 512u);
}

TEST(MiniSeedRecord, SetEncoding) {
    MiniSeedRecord record;
    
    record.setEncoding(MiniSeedRecord::Encoding::STEIM1);
    ASSERT_TRUE(record.encoding() == MiniSeedRecord::Encoding::STEIM1);
    
    record.setEncoding(MiniSeedRecord::Encoding::INT32);
    ASSERT_TRUE(record.encoding() == MiniSeedRecord::Encoding::INT32);
}

TEST(MiniSeedRecord, ToWaveform) {
    MiniSeedRecord record;
    
    StreamID id("IU", "ANMO", "00", "BHZ");
    record.setStreamId(id);
    record.setSampleRate(100.0);
    
    auto now = std::chrono::system_clock::now();
    record.setStartTime(now);
    
    SampleVector samples;
    for (int i = 0; i < 100; i++) {
        samples.push_back(std::sin(2.0 * M_PI * 5.0 * i / 100.0));
    }
    record.setSamples(samples);
    
    auto wf = record.toWaveform();
    
    ASSERT_TRUE(wf != nullptr);
    ASSERT_EQ(wf->streamId().toString(), "IU.ANMO.00.BHZ");
    ASSERT_NEAR(wf->sampleRate(), 100.0, 1e-10);
    ASSERT_EQ(wf->sampleCount(), 100u);
}

// ============================================================================
// MiniSeedReader Tests
// ============================================================================

TEST(MiniSeedReader, DefaultConstructor) {
    MiniSeedReader reader;
    
    ASSERT_EQ(reader.records().size(), 0u);
}

TEST(MiniSeedReader, EmptyFile) {
    MiniSeedReader reader;
    
    // Try to open non-existent file
    bool result = reader.open("/tmp/nonexistent_file.mseed");
    
    // Should fail gracefully
    ASSERT_FALSE(result);
}

TEST(MiniSeedReader, ParseEmptyData) {
    MiniSeedReader reader;
    
    uint8_t empty_data[10] = {0};
    bool result = reader.parse(empty_data, 10);
    
    // Too small to be valid MiniSEED
    ASSERT_FALSE(result);
}

TEST(MiniSeedReader, ToWaveformsEmpty) {
    MiniSeedReader reader;
    
    auto waveforms = reader.toWaveforms();
    
    ASSERT_EQ(waveforms.size(), 0u);
}

// ============================================================================
// Encoding Types Tests
// ============================================================================

TEST(MiniSeedRecord, EncodingTypes) {
    // Test all encoding type values
    ASSERT_EQ(static_cast<uint8_t>(MiniSeedRecord::Encoding::ASCII), 0);
    ASSERT_EQ(static_cast<uint8_t>(MiniSeedRecord::Encoding::INT16), 1);
    ASSERT_EQ(static_cast<uint8_t>(MiniSeedRecord::Encoding::INT24), 2);
    ASSERT_EQ(static_cast<uint8_t>(MiniSeedRecord::Encoding::INT32), 3);
    ASSERT_EQ(static_cast<uint8_t>(MiniSeedRecord::Encoding::FLOAT32), 4);
    ASSERT_EQ(static_cast<uint8_t>(MiniSeedRecord::Encoding::FLOAT64), 5);
    ASSERT_EQ(static_cast<uint8_t>(MiniSeedRecord::Encoding::STEIM1), 10);
    ASSERT_EQ(static_cast<uint8_t>(MiniSeedRecord::Encoding::STEIM2), 11);
}

// ============================================================================
// Header Size Tests
// ============================================================================

TEST(MiniSeedRecord, HeaderConstants) {
    ASSERT_EQ(MiniSeedRecord::HEADER_SIZE, 48u);
    ASSERT_EQ(MiniSeedRecord::FIXED_HEADER_SIZE, 48u);
}

// ============================================================================
// Synthetic MiniSEED Data Tests
// ============================================================================

TEST(MiniSeedRecord, SyntheticData) {
    // Create synthetic waveform data
    MiniSeedRecord record;
    
    StreamID id("SY", "TEST", "00", "BHZ");
    record.setStreamId(id);
    record.setSampleRate(100.0);
    record.setStartTime(std::chrono::system_clock::now());
    
    // Generate 1 second of sine wave
    SampleVector samples;
    for (int i = 0; i < 100; i++) {
        samples.push_back(1000.0 * std::sin(2.0 * M_PI * 5.0 * i / 100.0));
    }
    record.setSamples(samples);
    
    // Convert to waveform
    auto wf = record.toWaveform();
    
    ASSERT_TRUE(wf != nullptr);
    ASSERT_EQ(wf->sampleCount(), 100u);
    
    // Check data integrity
    double max_diff = 0;
    for (size_t i = 0; i < samples.size(); i++) {
        double expected = samples[i];
        double actual = (*wf)[i];
        max_diff = std::max(max_diff, std::abs(expected - actual));
    }
    
    ASSERT_LT(max_diff, 1e-10);
}

TEST(MiniSeedRecord, MultiChannelData) {
    // Test multiple records from different channels
    std::vector<MiniSeedRecord> records;
    
    std::vector<std::string> channels = {"BHZ", "BHN", "BHE"};
    auto now = std::chrono::system_clock::now();
    
    for (const auto& chan : channels) {
        MiniSeedRecord record;
        StreamID id("IU", "ANMO", "00", chan);
        record.setStreamId(id);
        record.setSampleRate(40.0);
        record.setStartTime(now);
        
        SampleVector samples;
        for (int i = 0; i < 40; i++) {
            samples.push_back(std::sin(2.0 * M_PI * 1.0 * i / 40.0));
        }
        record.setSamples(samples);
        
        records.push_back(record);
    }
    
    ASSERT_EQ(records.size(), 3u);
    
    // Convert each to waveform
    for (size_t i = 0; i < records.size(); i++) {
        auto wf = records[i].toWaveform();
        ASSERT_TRUE(wf != nullptr);
        ASSERT_EQ(wf->streamId().channel, channels[i]);
    }
}

// ============================================================================
// Waveform Continuity Tests
// ============================================================================

TEST(MiniSeedRecord, ContinuousRecords) {
    // Test that multiple records can represent continuous data
    auto start_time = std::chrono::system_clock::now();
    
    std::vector<MiniSeedRecord> records;
    
    // Create 10 consecutive records
    for (int r = 0; r < 10; r++) {
        MiniSeedRecord record;
        StreamID id("XX", "TEST", "00", "BHZ");
        record.setStreamId(id);
        record.setSampleRate(100.0);
        
        // Each record starts where the previous ended
        auto record_time = start_time + 
            std::chrono::milliseconds(r * 1000);  // 1 second per record
        record.setStartTime(record_time);
        
        SampleVector samples;
        for (int i = 0; i < 100; i++) {
            double t = r + i / 100.0;  // Continuous time
            samples.push_back(std::sin(2.0 * M_PI * 2.0 * t));
        }
        record.setSamples(samples);
        
        records.push_back(record);
    }
    
    ASSERT_EQ(records.size(), 10u);
    
    // Check time continuity
    for (size_t i = 1; i < records.size(); i++) {
        auto prev_end = records[i-1].startTime() + 
            std::chrono::milliseconds(static_cast<int>(
                1000.0 * records[i-1].sampleCount() / records[i-1].sampleRate()));
        auto curr_start = records[i].startTime();
        
        auto gap = std::chrono::duration_cast<std::chrono::milliseconds>(
            curr_start - prev_end).count();
        
        // Allow 10ms gap tolerance
        ASSERT_LT(std::abs(gap), 10);
    }
}

// ============================================================================
// Sample Rate Tests
// ============================================================================

TEST(MiniSeedRecord, VariousSampleRates) {
    std::vector<double> sample_rates = {1.0, 20.0, 40.0, 50.0, 100.0, 200.0, 500.0, 1000.0};
    
    for (double rate : sample_rates) {
        MiniSeedRecord record;
        record.setSampleRate(rate);
        
        ASSERT_NEAR(record.sampleRate(), rate, 1e-10);
    }
}

// ============================================================================
// Large Data Tests
// ============================================================================

TEST(MiniSeedRecord, LargeWaveform) {
    MiniSeedRecord record;
    
    StreamID id("XX", "TEST", "00", "BHZ");
    record.setStreamId(id);
    record.setSampleRate(100.0);
    record.setStartTime(std::chrono::system_clock::now());
    
    // 1 hour of data at 100 Hz = 360,000 samples
    SampleVector samples;
    samples.reserve(360000);
    
    for (int i = 0; i < 360000; i++) {
        samples.push_back(std::sin(2.0 * M_PI * 1.0 * i / 100.0));
    }
    record.setSamples(samples);
    
    ASSERT_EQ(record.sampleCount(), 360000u);
    
    auto wf = record.toWaveform();
    ASSERT_TRUE(wf != nullptr);
    ASSERT_EQ(wf->sampleCount(), 360000u);
    
    // Duration should be 1 hour = 3600 seconds
    ASSERT_NEAR(wf->duration(), 3600.0, 0.1);
}
