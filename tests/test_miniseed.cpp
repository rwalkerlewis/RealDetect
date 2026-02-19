/**
 * Unit tests for MiniSEED parsing
 *
 * Tests MiniSEED record data structures and reading real EarthScope data.
 */

#include "test_framework.hpp"
#include "realdetect/core/miniseed.hpp"
#include "realdetect/core/waveform.hpp"
#include <cmath>
#include <fstream>

using namespace realdetect;
using namespace realdetect::test;

// ============================================================================
// MiniSeedRecord Data Structure Tests
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
    ASSERT_LT(std::abs(diff.count()), 1000);
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
    record.setStartTime(std::chrono::system_clock::now());
    SampleVector samples = {1.0, 2.0, 3.0, 4.0, 5.0};
    record.setSamples(samples);
    auto wf = record.toWaveform();
    ASSERT_TRUE(wf != nullptr);
    ASSERT_EQ(wf->streamId().toString(), "IU.ANMO.00.BHZ");
    ASSERT_NEAR(wf->sampleRate(), 100.0, 1e-10);
    ASSERT_EQ(wf->sampleCount(), 5u);
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
    bool result = reader.open("/tmp/nonexistent_file.mseed");
    ASSERT_FALSE(result);
}

TEST(MiniSeedReader, ParseEmptyData) {
    MiniSeedReader reader;
    uint8_t empty_data[10] = {0};
    bool result = reader.parse(empty_data, 10);
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
    ASSERT_EQ(static_cast<uint8_t>(MiniSeedRecord::Encoding::ASCII), 0);
    ASSERT_EQ(static_cast<uint8_t>(MiniSeedRecord::Encoding::INT16), 1);
    ASSERT_EQ(static_cast<uint8_t>(MiniSeedRecord::Encoding::INT24), 2);
    ASSERT_EQ(static_cast<uint8_t>(MiniSeedRecord::Encoding::INT32), 3);
    ASSERT_EQ(static_cast<uint8_t>(MiniSeedRecord::Encoding::FLOAT32), 4);
    ASSERT_EQ(static_cast<uint8_t>(MiniSeedRecord::Encoding::FLOAT64), 5);
    ASSERT_EQ(static_cast<uint8_t>(MiniSeedRecord::Encoding::STEIM1), 10);
    ASSERT_EQ(static_cast<uint8_t>(MiniSeedRecord::Encoding::STEIM2), 11);
}

TEST(MiniSeedRecord, HeaderConstants) {
    ASSERT_EQ(MiniSeedRecord::HEADER_SIZE, 48u);
    ASSERT_EQ(MiniSeedRecord::FIXED_HEADER_SIZE, 48u);
}

// ============================================================================
// Real EarthScope MiniSEED Data Tests
// ============================================================================

TEST(MiniSeedReader, ReadRealTestFixture) {
    // Load real MiniSEED data fetched from EarthScope
    MiniSeedReader reader;
    bool opened = reader.open("data/test/waveforms.mseed");
    if (!opened) opened = reader.open("../data/test/waveforms.mseed");

    if (!opened) {
        std::cout << "    [SKIP] No real test MiniSEED data available\n";
        return;
    }

    // Should have parsed real records
    ASSERT_GT(reader.records().size(), 0u);

    // Check record properties
    for (const auto& rec : reader.records()) {
        ASSERT_GT(rec.sampleRate(), 0);
        ASSERT_GT(rec.sampleCount(), 0u);
        ASSERT_FALSE(rec.streamId().station.empty());
    }
}

TEST(MiniSeedReader, RealDataToWaveforms) {
    MiniSeedReader reader;
    bool opened = reader.open("data/test/waveforms.mseed");
    if (!opened) opened = reader.open("../data/test/waveforms.mseed");

    if (!opened) {
        std::cout << "    [SKIP] No real test data\n";
        return;
    }

    auto waveforms = reader.toWaveforms();
    ASSERT_GT(waveforms.size(), 0u);

    for (const auto& wf : waveforms) {
        ASSERT_GT(wf->sampleCount(), 100u);  // Real data has substantial samples
        ASSERT_GT(wf->sampleRate(), 0);
        ASSERT_FALSE(wf->streamId().station.empty());

        // Real seismic data should have non-zero amplitudes
        ASSERT_GT(wf->absMax(), 0);
    }
}

TEST(MiniSeedReader, RealDataStationIdentification) {
    MiniSeedReader reader;
    bool opened = reader.open("data/test/waveforms.mseed");
    if (!opened) opened = reader.open("../data/test/waveforms.mseed");

    if (!opened) {
        std::cout << "    [SKIP] No real test data\n";
        return;
    }

    auto waveforms = reader.toWaveforms();

    // Should have CI network stations (Ridgecrest fixture)
    bool found_ci = false;
    for (const auto& wf : waveforms) {
        if (wf->streamId().network == "CI") found_ci = true;
    }
    ASSERT_TRUE(found_ci);
}

TEST(MiniSeedReader, ReadRealRidgecrestData) {
    MiniSeedReader reader;
    bool opened = reader.open("data/ridgecrest/waveforms.mseed");
    if (!opened) opened = reader.open("../data/ridgecrest/waveforms.mseed");

    if (!opened) {
        std::cout << "    [SKIP] No Ridgecrest MiniSEED data available\n";
        return;
    }

    auto waveforms = reader.toWaveforms();
    ASSERT_GE(waveforms.size(), 10u);  // Should have many stations

    // Verify BHZ channels present
    int bhz_count = 0;
    for (const auto& wf : waveforms) {
        if (wf->streamId().channel == "BHZ") bhz_count++;
    }
    ASSERT_GT(bhz_count, 5);

    std::cout << "    Loaded " << waveforms.size() << " waveforms ("
              << bhz_count << " BHZ) from real Ridgecrest data\n";
}

TEST(MiniSeedReader, VariousSampleRates) {
    std::vector<double> expected_rates = {1.0, 20.0, 40.0, 50.0, 100.0, 200.0, 500.0, 1000.0};
    MiniSeedRecord record;
    for (double rate : expected_rates) {
        record.setSampleRate(rate);
        ASSERT_NEAR(record.sampleRate(), rate, 1e-10);
    }
}
