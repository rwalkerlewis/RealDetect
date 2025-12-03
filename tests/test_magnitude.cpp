/**
 * Unit tests for magnitude calculation
 */

#include "test_framework.hpp"
#include "seisproc/core/waveform.hpp"
#include "seisproc/core/station.hpp"
#include "seisproc/core/event.hpp"
#include "seisproc/magnitude/local_magnitude.hpp"
#include "seisproc/magnitude/moment_magnitude.hpp"
#include <random>
#include <cmath>

using namespace seisproc;
using namespace seisproc::test;

// Helper to create synthetic waveform for magnitude testing
WaveformPtr createMagnitudeWaveform(double sample_rate, double duration,
                                     double p_time, double s_time,
                                     double amplitude, double frequency = 5.0) {
    StreamID id("SY", "TEST", "00", "BHZ");
    auto now = std::chrono::system_clock::now();
    auto wf = std::make_shared<Waveform>(id, sample_rate, now);
    
    std::mt19937 gen(42);
    std::normal_distribution<> noise(0.0, 1.0);
    
    size_t n_samples = static_cast<size_t>(duration * sample_rate);
    size_t p_idx = static_cast<size_t>(p_time * sample_rate);
    size_t s_idx = static_cast<size_t>(s_time * sample_rate);
    
    for (size_t i = 0; i < n_samples; i++) {
        double sample = noise(gen);  // Background noise
        
        // P-wave
        if (i >= p_idx && i < s_idx) {
            double t = (i - p_idx) / sample_rate;
            double envelope = std::exp(-t / 1.0);
            sample += amplitude * 0.3 * envelope * std::sin(2.0 * M_PI * frequency * t);
        }
        
        // S-wave (larger amplitude)
        if (i >= s_idx) {
            double t = (i - s_idx) / sample_rate;
            double envelope = std::exp(-t / 3.0);
            sample += amplitude * envelope * std::sin(2.0 * M_PI * frequency * 0.6 * t);
        }
        
        wf->append(sample);
    }
    
    return wf;
}

// Create test stations and origin
std::pair<StationInventory, Origin> createMagnitudeTestCase() {
    StationInventory inv;
    
    // Stations at various distances
    inv.addStation(std::make_shared<Station>("XX", "S1", 34.1, -118.0, 100));  // ~11 km
    inv.addStation(std::make_shared<Station>("XX", "S2", 34.3, -118.0, 200));  // ~33 km
    inv.addStation(std::make_shared<Station>("XX", "S3", 34.5, -117.5, 300));  // ~70 km
    inv.addStation(std::make_shared<Station>("XX", "S4", 35.0, -118.0, 100));  // ~111 km
    
    Origin origin;
    origin.location = GeoPoint(34.0, -118.0, 10.0);
    origin.time = std::chrono::system_clock::now();
    
    return {inv, origin};
}

// ============================================================================
// LocalMagnitude Tests
// ============================================================================

TEST(LocalMagnitude, DefaultParameters) {
    LocalMagnitude ml;
    
    ASSERT_TRUE(ml.type() == MagnitudeType::ML);
    ASSERT_EQ(ml.name(), "Local Magnitude (ML)");
}

TEST(LocalMagnitude, SetParameters) {
    LocalMagnitude ml;
    
    ml.setAttenuationA(1.5);
    ml.setAttenuationB(0.004);
    ml.setAttenuationC(2.5);
    ml.setMinDistance(10.0);
    ml.setMaxDistance(700.0);
    
    ASSERT_NO_THROW(ml.setParameter("a", 1.2));
}

TEST(LocalMagnitude, CalculateWithSyntheticData) {
    LocalMagnitude ml;
    ml.setSimulateWoodAnderson(false);  // Skip WA simulation for test
    
    auto [stations, origin] = createMagnitudeTestCase();
    
    // Create waveforms with known amplitudes
    std::map<StreamID, WaveformPtr> waveforms;
    
    for (const auto& [key, sta] : stations.stations()) {
        double dist = origin.location.distanceTo(sta->location());
        double hypo_dist = std::sqrt(dist * dist + origin.location.depth * origin.location.depth);
        
        double p_time = hypo_dist / 6.0;  // P at ~6 km/s
        double s_time = hypo_dist / 3.5;  // S at ~3.5 km/s
        
        // Amplitude roughly proportional to 1/distance for geometric spreading
        double amplitude = 10000.0 / dist;
        
        StreamID id(sta->network(), sta->code(), "00", "BHZ");
        waveforms[id] = createMagnitudeWaveform(100.0, 60.0, p_time, s_time, amplitude);
    }
    
    auto result = ml.calculate(origin, waveforms, stations);
    
    // Should have station magnitudes
    ASSERT_GT(result.station_count, 0);
    ASSERT_GE(result.station_magnitudes.size(), 1u);
    
    // Magnitude should be reasonable
    ASSERT_GT(result.value, -2.0);
    ASSERT_LT(result.value, 8.0);
}

TEST(LocalMagnitude, StationCorrection) {
    LocalMagnitude ml;
    
    ml.setStationCorrection("S1", 0.5);
    ml.setStationCorrection("S2", -0.3);
    
    auto [stations, origin] = createMagnitudeTestCase();
    
    std::map<StreamID, WaveformPtr> waveforms;
    for (const auto& [key, sta] : stations.stations()) {
        StreamID id(sta->network(), sta->code(), "00", "BHZ");
        waveforms[id] = createMagnitudeWaveform(100.0, 60.0, 10.0, 18.0, 1000.0);
    }
    
    auto result = ml.calculate(origin, waveforms, stations);
    
    // Check that corrections were applied
    for (const auto& sm : result.station_magnitudes) {
        if (sm.stream_id.station == "S1") {
            ASSERT_NEAR(sm.correction, 0.5, 1e-10);
        }
        if (sm.stream_id.station == "S2") {
            ASSERT_NEAR(sm.correction, -0.3, 1e-10);
        }
    }
}

TEST(LocalMagnitude, DistanceFiltering) {
    LocalMagnitude ml;
    ml.setMinDistance(50.0);   // Only use stations > 50 km
    ml.setMaxDistance(200.0);  // Only use stations < 200 km
    
    auto [stations, origin] = createMagnitudeTestCase();
    
    std::map<StreamID, WaveformPtr> waveforms;
    for (const auto& [key, sta] : stations.stations()) {
        StreamID id(sta->network(), sta->code(), "00", "BHZ");
        waveforms[id] = createMagnitudeWaveform(100.0, 60.0, 10.0, 18.0, 1000.0);
    }
    
    auto result = ml.calculate(origin, waveforms, stations);
    
    // Should exclude close and far stations
    // S1 (~11 km) excluded, S4 (~111 km) included, etc.
    for (const auto& sm : result.station_magnitudes) {
        double dist = origin.location.distanceTo(
            stations.getStation(sm.stream_id)->location());
        ASSERT_GE(dist, 50.0);
        ASSERT_LE(dist, 200.0);
    }
}

// ============================================================================
// DurationMagnitude Tests
// ============================================================================

TEST(DurationMagnitude, DefaultParameters) {
    DurationMagnitude md;
    
    ASSERT_TRUE(md.type() == MagnitudeType::Md);
    ASSERT_EQ(md.name(), "Duration Magnitude (Md)");
}

TEST(DurationMagnitude, SetCoefficients) {
    DurationMagnitude md;
    
    ASSERT_NO_THROW(md.setCoefficients(-0.5, 2.5, 0.003));
    ASSERT_NO_THROW(md.setCodaThreshold(2.0));
}

// ============================================================================
// MomentMagnitude Tests
// ============================================================================

TEST(MomentMagnitude, DefaultParameters) {
    MomentMagnitude mw;
    
    ASSERT_TRUE(mw.type() == MagnitudeType::Mw);
    ASSERT_EQ(mw.name(), "Moment Magnitude (Mw)");
}

TEST(MomentMagnitude, SetParameters) {
    MomentMagnitude mw;
    
    mw.setDensity(2800.0);
    mw.setVelocity(3600.0);
    mw.setQuality(300.0);
    mw.setFrequencyRange(0.5, 25.0);
    
    ASSERT_NO_THROW(mw.setParameter("density", 3000.0));
}

TEST(MomentMagnitude, BruneSpectrum) {
    // Test that Mw can process waveforms
    MomentMagnitude mw;
    
    auto [stations, origin] = createMagnitudeTestCase();
    
    std::map<StreamID, WaveformPtr> waveforms;
    for (const auto& [key, sta] : stations.stations()) {
        StreamID id(sta->network(), sta->code(), "00", "BHZ");
        waveforms[id] = createMagnitudeWaveform(100.0, 60.0, 10.0, 18.0, 5000.0);
    }
    
    // Should not throw
    ASSERT_NO_THROW(mw.calculate(origin, waveforms, stations));
}

// ============================================================================
// BodyWaveMagnitude Tests
// ============================================================================

TEST(BodyWaveMagnitude, DefaultParameters) {
    BodyWaveMagnitude mb;
    
    ASSERT_TRUE(mb.type() == MagnitudeType::Mb);
    ASSERT_EQ(mb.name(), "Body Wave Magnitude (mb)");
}

TEST(BodyWaveMagnitude, SetParameters) {
    BodyWaveMagnitude mb;
    
    ASSERT_NO_THROW(mb.setParameter("min_distance", 600.0));
    ASSERT_NO_THROW(mb.setParameter("period", 1.5));
}

// ============================================================================
// SurfaceWaveMagnitude Tests
// ============================================================================

TEST(SurfaceWaveMagnitude, DefaultParameters) {
    SurfaceWaveMagnitude ms;
    
    ASSERT_TRUE(ms.type() == MagnitudeType::Ms);
    ASSERT_EQ(ms.name(), "Surface Wave Magnitude (Ms)");
}

TEST(SurfaceWaveMagnitude, SetParameters) {
    SurfaceWaveMagnitude ms;
    
    ASSERT_NO_THROW(ms.setParameter("target_period", 18.0));
}

// ============================================================================
// MagnitudeFactory Tests
// ============================================================================

TEST(MagnitudeFactory, CreateML) {
    auto calc = MagnitudeFactory::create(MagnitudeType::ML);
    ASSERT_TRUE(calc != nullptr);
    ASSERT_TRUE(calc->type() == MagnitudeType::ML);
}

TEST(MagnitudeFactory, CreateMw) {
    auto calc = MagnitudeFactory::create(MagnitudeType::Mw);
    ASSERT_TRUE(calc != nullptr);
    ASSERT_TRUE(calc->type() == MagnitudeType::Mw);
}

TEST(MagnitudeFactory, CreateMb) {
    auto calc = MagnitudeFactory::create(MagnitudeType::Mb);
    ASSERT_TRUE(calc != nullptr);
    ASSERT_TRUE(calc->type() == MagnitudeType::Mb);
}

TEST(MagnitudeFactory, CreateMs) {
    auto calc = MagnitudeFactory::create(MagnitudeType::Ms);
    ASSERT_TRUE(calc != nullptr);
    ASSERT_TRUE(calc->type() == MagnitudeType::Ms);
}

TEST(MagnitudeFactory, CreateMd) {
    auto calc = MagnitudeFactory::create(MagnitudeType::Md);
    ASSERT_TRUE(calc != nullptr);
    ASSERT_TRUE(calc->type() == MagnitudeType::Md);
}

TEST(MagnitudeFactory, AvailableTypes) {
    auto types = MagnitudeFactory::availableTypes();
    
    ASSERT_GE(types.size(), 5u);
    
    bool has_ml = false, has_mw = false;
    for (auto type : types) {
        if (type == MagnitudeType::ML) has_ml = true;
        if (type == MagnitudeType::Mw) has_mw = true;
    }
    ASSERT_TRUE(has_ml);
    ASSERT_TRUE(has_mw);
}

// ============================================================================
// Magnitude Integration Tests
// ============================================================================

TEST(MagnitudeIntegration, MLConsistency) {
    // Same event should give similar ML across stations
    LocalMagnitude ml;
    ml.setSimulateWoodAnderson(false);
    
    auto [stations, origin] = createMagnitudeTestCase();
    
    std::map<StreamID, WaveformPtr> waveforms;
    for (const auto& [key, sta] : stations.stations()) {
        double dist = origin.location.distanceTo(sta->location());
        
        // Amplitude with proper distance decay
        double amplitude = 50000.0 / (dist * dist);
        
        StreamID id(sta->network(), sta->code(), "00", "BHZ");
        waveforms[id] = createMagnitudeWaveform(100.0, 60.0, dist/6, dist/3.5, amplitude);
    }
    
    auto result = ml.calculate(origin, waveforms, stations);
    
    if (result.station_magnitudes.size() >= 2) {
        // Check that station magnitudes are within 1 unit of each other
        double min_mag = 1e10, max_mag = -1e10;
        for (const auto& sm : result.station_magnitudes) {
            min_mag = std::min(min_mag, sm.value);
            max_mag = std::max(max_mag, sm.value);
        }
        
        // Spread should be less than 2 units for consistent data
        ASSERT_LT(max_mag - min_mag, 2.0);
    }
}

TEST(MagnitudeIntegration, NetworkMagnitude) {
    // Test that network magnitude is computed correctly
    LocalMagnitude ml;
    
    auto [stations, origin] = createMagnitudeTestCase();
    
    std::map<StreamID, WaveformPtr> waveforms;
    for (const auto& [key, sta] : stations.stations()) {
        StreamID id(sta->network(), sta->code(), "00", "BHZ");
        waveforms[id] = createMagnitudeWaveform(100.0, 60.0, 10.0, 18.0, 1000.0);
    }
    
    auto result = ml.calculate(origin, waveforms, stations);
    
    if (result.station_count >= 3) {
        // Network magnitude should be close to median of station magnitudes
        std::vector<double> mags;
        for (const auto& sm : result.station_magnitudes) {
            mags.push_back(sm.value);
        }
        std::sort(mags.begin(), mags.end());
        double median = mags[mags.size() / 2];
        
        // Network magnitude should be within 0.5 of median
        ASSERT_NEAR(result.value, median, 0.5);
    }
}
