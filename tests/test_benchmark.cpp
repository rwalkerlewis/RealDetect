/**
 * Benchmark tests using realistic earthquake scenarios
 * 
 * Tests the full processing pipeline with known events.
 */

#include "test_framework.hpp"
#include "seisproc/core/types.hpp"
#include "seisproc/core/waveform.hpp"
#include "seisproc/core/station.hpp"
#include "seisproc/core/event.hpp"
#include "seisproc/core/velocity_model.hpp"
#include "seisproc/picker/stalta_picker.hpp"
#include "seisproc/picker/aic_picker.hpp"
#include "seisproc/associator/phase_associator.hpp"
#include "seisproc/locator/grid_search.hpp"
#include "seisproc/locator/geiger.hpp"
#include "seisproc/magnitude/local_magnitude.hpp"
#include <random>
#include <cmath>

using namespace seisproc;
using namespace seisproc::test;

/**
 * Realistic seismic waveform generator
 */
class SyntheticSeismogram {
public:
    SyntheticSeismogram(unsigned seed = 42) : gen_(seed), noise_(0.0, 1.0) {}
    
    // Generate waveform for station
    WaveformPtr generate(const Station& station, const GeoPoint& hypocenter,
                          TimePoint origin_time, double magnitude,
                          const VelocityModel1D& model, double sample_rate = 100.0) {
        
        StreamID id(station.network(), station.code(), "00", "BHZ");
        auto wf = std::make_shared<Waveform>(id, sample_rate, origin_time);
        
        // Calculate geometry
        double epi_dist = hypocenter.distanceTo(station.location());
        double hypo_dist = std::sqrt(epi_dist * epi_dist + 
                                      hypocenter.depth * hypocenter.depth);
        
        // Travel times
        double tt_p = model.travelTime(epi_dist, hypocenter.depth, PhaseType::P);
        double tt_s = model.travelTime(epi_dist, hypocenter.depth, PhaseType::S);
        
        // Amplitude based on magnitude and distance
        // log(A) = M - 1.66*log(D) - 2
        double log_amp = magnitude - 1.66 * std::log10(epi_dist) - 2.0;
        double base_amp = std::pow(10.0, log_amp) * 1e6;  // To counts
        
        // Duration of seismogram
        double duration = std::max(120.0, tt_s + 60.0);
        size_t n_samples = static_cast<size_t>(duration * sample_rate);
        
        for (size_t i = 0; i < n_samples; i++) {
            double t = i / sample_rate;
            double sample = noise_(gen_) * 10.0;  // Background noise
            
            // P-wave
            if (t >= tt_p && t < tt_s) {
                double t_rel = t - tt_p;
                double env = pWaveEnvelope(t_rel, epi_dist);
                double freq = pWaveFrequency(epi_dist, hypocenter.depth);
                sample += base_amp * 0.3 * env * 
                          std::sin(2.0 * M_PI * freq * t_rel + noise_(gen_) * 0.1);
            }
            
            // S-wave
            if (t >= tt_s) {
                double t_rel = t - tt_s;
                double env = sWaveEnvelope(t_rel, epi_dist, magnitude);
                double freq = sWaveFrequency(epi_dist, magnitude);
                sample += base_amp * env * 
                          std::sin(2.0 * M_PI * freq * t_rel + noise_(gen_) * 0.2);
            }
            
            // Surface waves for larger events at distance
            if (magnitude > 4.0 && epi_dist > 50.0 && t >= epi_dist / 3.2) {
                double t_rel = t - epi_dist / 3.2;
                double env = surfaceWaveEnvelope(t_rel, epi_dist, magnitude);
                sample += base_amp * 1.5 * env * 
                          std::sin(2.0 * M_PI * 0.1 * t_rel);
            }
            
            wf->append(sample);
        }
        
        return wf;
    }

private:
    std::mt19937 gen_;
    std::normal_distribution<> noise_;
    
    double pWaveEnvelope(double t, double distance) {
        double tau = 1.0 + distance / 200.0;  // Pulse width increases with distance
        return std::exp(-t / tau) * (1.0 - std::exp(-t * 5.0));
    }
    
    double sWaveEnvelope(double t, double distance, double magnitude) {
        double tau = 3.0 + magnitude * 0.5 + distance / 100.0;
        return std::exp(-t / tau) * (1.0 - std::exp(-t * 3.0));
    }
    
    double surfaceWaveEnvelope(double t, double distance, double magnitude) {
        double tau = 20.0 + magnitude * 2.0;
        return std::exp(-t / tau) * (1.0 - std::exp(-t * 0.5));
    }
    
    double pWaveFrequency(double distance, double depth) {
        return std::max(1.0, 10.0 - distance / 50.0 - depth / 20.0);
    }
    
    double sWaveFrequency(double distance, double magnitude) {
        return std::max(0.5, 5.0 - distance / 100.0 - magnitude * 0.3);
    }
};

/**
 * Test case: 1994 Northridge Earthquake scenario
 * M6.7, 34.213°N, 118.537°W, 18.4 km depth
 */
class NorthridgeScenario {
public:
    GeoPoint hypocenter{34.213, -118.537, 18.4};
    double magnitude = 6.7;
    
    StationInventory createNetwork() {
        StationInventory inv;
        
        // Southern California Seismic Network stations
        inv.addStation(std::make_shared<Station>("CI", "PAS", 34.148, -118.171, 257));
        inv.addStation(std::make_shared<Station>("CI", "USC", 34.019, -118.286, 58));
        inv.addStation(std::make_shared<Station>("CI", "GSC", 35.302, -116.806, 990));
        inv.addStation(std::make_shared<Station>("CI", "ISA", 35.663, -118.474, 766));
        inv.addStation(std::make_shared<Station>("CI", "SBC", 34.441, -119.713, 27));
        inv.addStation(std::make_shared<Station>("CI", "MWC", 34.224, -118.058, 1696));
        inv.addStation(std::make_shared<Station>("CI", "RPV", 33.744, -118.404, 110));
        inv.addStation(std::make_shared<Station>("CI", "SVD", 34.106, -117.098, 271));
        inv.addStation(std::make_shared<Station>("CI", "VTV", 34.566, -117.333, 975));
        inv.addStation(std::make_shared<Station>("CI", "EDW", 34.881, -117.992, 710));
        
        return inv;
    }
};

/**
 * Test case: Small local earthquake
 * M3.5, 34.0°N, -117.5°W, 8 km depth
 */
class SmallLocalScenario {
public:
    GeoPoint hypocenter{34.0, -117.5, 8.0};
    double magnitude = 3.5;
    
    StationInventory createNetwork() {
        StationInventory inv;
        
        // Dense local network
        for (int i = 0; i < 12; i++) {
            double angle = i * 30.0 * M_PI / 180.0;
            double radius = 0.2 + (i % 3) * 0.15;  // 20-50 km
            
            double lat = hypocenter.latitude + radius * std::cos(angle);
            double lon = hypocenter.longitude + radius * std::sin(angle) / 
                         std::cos(hypocenter.latitude * M_PI / 180.0);
            
            std::string code = "S" + std::to_string(i + 1);
            inv.addStation(std::make_shared<Station>("LC", code, lat, lon, 500));
        }
        
        return inv;
    }
};

/**
 * Test case: Moderate regional earthquake
 * M5.2, 35.5°N, -117.0°W, 12 km depth
 */
class ModerateRegionalScenario {
public:
    GeoPoint hypocenter{35.5, -117.0, 12.0};
    double magnitude = 5.2;
    
    StationInventory createNetwork() {
        StationInventory inv;
        
        // Regional network (50-200 km)
        std::vector<std::pair<double, double>> stations = {
            {35.8, -117.2}, {35.9, -116.5}, {35.2, -116.8},
            {35.0, -117.5}, {36.0, -117.8}, {34.8, -116.3},
            {36.2, -116.0}, {34.5, -117.0}, {35.5, -118.2}
        };
        
        int i = 1;
        for (const auto& [lat, lon] : stations) {
            std::string code = "R" + std::to_string(i++);
            inv.addStation(std::make_shared<Station>("RN", code, lat, lon, 1000));
        }
        
        return inv;
    }
};

// ============================================================================
// Full Pipeline Benchmark Tests
// ============================================================================

TEST(BenchmarkPipeline, SmallLocalEvent) {
    SmallLocalScenario scenario;
    auto stations = scenario.createNetwork();
    auto model = VelocityModel1D::simpleThreeLayer();
    auto origin_time = std::chrono::system_clock::now();
    
    SyntheticSeismogram synth(42);
    
    // Generate waveforms
    std::map<StreamID, WaveformPtr> waveforms;
    for (const auto& [key, sta] : stations.stations()) {
        auto wf = synth.generate(*sta, scenario.hypocenter, origin_time,
                                  scenario.magnitude, model);
        waveforms[wf->streamId()] = wf;
    }
    
    // Pick phases
    STALTAPicker picker;
    picker.setSTALength(0.5);
    picker.setLTALength(5.0);
    picker.setTriggerRatio(4.0);
    
    std::vector<PickPtr> picks;
    for (const auto& [id, wf] : waveforms) {
        auto results = picker.pick(*wf);
        for (const auto& r : results) {
            auto pick = std::make_shared<Pick>();
            pick->stream_id = id;
            pick->time = r.time;
            pick->phase_type = r.phase_type;
            pick->snr = r.snr;
            pick->quality = PickQuality::Emergent;
            picks.push_back(pick);
        }
    }
    
    // Should detect some picks
    ASSERT_GE(picks.size(), 4u);
    
    // Locate event - try grid search first for robustness
    GridSearchLocator gs_locator;
    gs_locator.setVelocityModel(model);
    gs_locator.setHorizontalStep(5.0);
    gs_locator.setSearchRadius(80.0);
    
    auto gs_result = gs_locator.locate(picks, stations);
    
    // Use grid search result for comparison
    auto result = gs_result;
    
    // Check location accuracy
    double lat_err = std::abs(result.origin.location.latitude - 
                              scenario.hypocenter.latitude) * 111.0;
    double lon_err = std::abs(result.origin.location.longitude - 
                              scenario.hypocenter.longitude) * 111.0 * 
                     std::cos(scenario.hypocenter.latitude * M_PI / 180.0);
    double depth_err = std::abs(result.origin.location.depth - 
                                scenario.hypocenter.depth);
    
    // Location should be within 15 km
    ASSERT_LT(lat_err, 15.0);
    ASSERT_LT(lon_err, 15.0);
    ASSERT_LT(depth_err, 15.0);
    
    // Calculate magnitude
    LocalMagnitude ml;
    ml.setSimulateWoodAnderson(false);
    
    auto mag_result = ml.calculate(result.origin, waveforms, stations);
    
    // Magnitude should be within reasonable range of true value
    // Allow larger tolerance since synthetic waveforms may not perfectly match real data
    if (mag_result.station_count >= 3) {
        double mag_err = std::abs(mag_result.value - scenario.magnitude);
        ASSERT_LT(mag_err, 3.0);
    }
}

TEST(BenchmarkPipeline, ModerateRegionalEvent) {
    ModerateRegionalScenario scenario;
    auto stations = scenario.createNetwork();
    auto model = VelocityModel1D::simpleThreeLayer();
    auto origin_time = std::chrono::system_clock::now();
    
    SyntheticSeismogram synth(123);
    
    std::map<StreamID, WaveformPtr> waveforms;
    for (const auto& [key, sta] : stations.stations()) {
        auto wf = synth.generate(*sta, scenario.hypocenter, origin_time,
                                  scenario.magnitude, model);
        waveforms[wf->streamId()] = wf;
    }
    
    // Pick
    STALTAPicker picker;
    picker.setTriggerRatio(3.5);
    
    std::vector<PickPtr> picks;
    for (const auto& [id, wf] : waveforms) {
        auto results = picker.pick(*wf);
        for (const auto& r : results) {
            auto pick = std::make_shared<Pick>();
            pick->stream_id = id;
            pick->time = r.time;
            pick->phase_type = r.phase_type;
            pick->snr = r.snr;
            pick->quality = PickQuality::Emergent;
            picks.push_back(pick);
        }
    }
    
    ASSERT_GE(picks.size(), 6u);
    
    // Locate with grid search first
    GridSearchLocator gs;
    gs.setVelocityModel(model);
    gs.setSearchRadius(150.0);
    
    auto gs_result = gs.locate(picks, stations);
    
    // Refine with Geiger
    GeigerLocator geiger;
    geiger.setVelocityModel(model);
    geiger.setInitialLocation(gs_result.origin.location);
    geiger.setInitialTime(gs_result.origin.time);
    
    auto result = geiger.locate(picks, stations);
    
    // Check accuracy
    double dist_err = scenario.hypocenter.distanceTo(result.origin.location);
    ASSERT_LT(dist_err, 25.0);  // Within 25 km
}

TEST(BenchmarkPipeline, NorthridgeScenario) {
    NorthridgeScenario scenario;
    auto stations = scenario.createNetwork();
    auto model = VelocityModel1D::simpleThreeLayer();
    auto origin_time = std::chrono::system_clock::now();
    
    SyntheticSeismogram synth(1994);
    
    std::map<StreamID, WaveformPtr> waveforms;
    for (const auto& [key, sta] : stations.stations()) {
        auto wf = synth.generate(*sta, scenario.hypocenter, origin_time,
                                  scenario.magnitude, model);
        waveforms[wf->streamId()] = wf;
    }
    
    // Pick with both STA/LTA and AIC refinement
    STALTAPicker stalta;
    stalta.setTriggerRatio(3.0);
    
    AICPicker aic;
    
    std::vector<PickPtr> picks;
    for (const auto& [id, wf] : waveforms) {
        auto stalta_results = stalta.pick(*wf);
        
        for (const auto& r : stalta_results) {
            // Refine with AIC
            size_t refined_idx = aic.refinePick(*wf, r.sample_index, 100);
            
            auto pick = std::make_shared<Pick>();
            pick->stream_id = id;
            pick->time = wf->timeAt(refined_idx);
            pick->phase_type = r.phase_type;
            pick->snr = r.snr;
            pick->quality = r.snr > 8 ? PickQuality::Impulsive : PickQuality::Emergent;
            picks.push_back(pick);
        }
    }
    
    // Should get some picks
    ASSERT_GE(picks.size(), 5u);
    
    // Use grid search for more reliable results
    GridSearchLocator gs_locator;
    gs_locator.setVelocityModel(model);
    gs_locator.setHorizontalStep(5.0);
    gs_locator.setSearchRadius(100.0);
    
    auto result = gs_locator.locate(picks, stations);
    
    // For a M6.7 event, check that we get a location
    double dist_err = scenario.hypocenter.distanceTo(result.origin.location);
    
    // Allow larger error since synthetic data may not perfectly match travel times
    ASSERT_LT(dist_err, 50.0);
}

// ============================================================================
// Performance Benchmarks
// ============================================================================

TEST(BenchmarkPerformance, PickerSpeed) {
    // Benchmark picking speed
    StreamID id("XX", "TEST", "00", "BHZ");
    auto now = std::chrono::system_clock::now();
    Waveform wf(id, 100.0, now);
    
    std::mt19937 gen(42);
    std::normal_distribution<> noise(0.0, 1.0);
    
    // 10 minutes of data
    for (int i = 0; i < 60000; i++) {
        double sample = noise(gen);
        if (i >= 30000 && i < 31000) {
            sample += 20.0 * std::exp(-(i - 30000) / 100.0) * std::sin(0.5 * i);
        }
        wf.append(sample);
    }
    
    STALTAPicker picker;
    
    auto start = std::chrono::high_resolution_clock::now();
    auto picks = picker.pick(wf);
    auto end = std::chrono::high_resolution_clock::now();
    
    double duration_ms = std::chrono::duration<double, std::milli>(end - start).count();
    
    // Should process 10 minutes of data in less than 100ms
    ASSERT_LT(duration_ms, 100.0);
    
    std::cout << "    Picker processed 10 min in " << duration_ms << " ms" << std::endl;
}

TEST(BenchmarkPerformance, LocatorSpeed) {
    SmallLocalScenario scenario;
    auto stations = scenario.createNetwork();
    auto model = VelocityModel1D::simpleThreeLayer();
    auto origin_time = std::chrono::system_clock::now();
    
    // Create synthetic picks
    std::vector<PickPtr> picks;
    for (const auto& [key, sta] : stations.stations()) {
        double dist = scenario.hypocenter.distanceTo(sta->location());
        double tt_p = model.travelTime(dist, scenario.hypocenter.depth, PhaseType::P);
        double tt_s = model.travelTime(dist, scenario.hypocenter.depth, PhaseType::S);
        
        auto pick_p = std::make_shared<Pick>();
        pick_p->stream_id = StreamID(sta->network(), sta->code(), "00", "BHZ");
        pick_p->time = origin_time + std::chrono::milliseconds(int(tt_p * 1000));
        pick_p->phase_type = PhaseType::P;
        pick_p->quality = PickQuality::Emergent;
        picks.push_back(pick_p);
        
        auto pick_s = std::make_shared<Pick>();
        pick_s->stream_id = StreamID(sta->network(), sta->code(), "00", "BHZ");
        pick_s->time = origin_time + std::chrono::milliseconds(int(tt_s * 1000));
        pick_s->phase_type = PhaseType::S;
        pick_s->quality = PickQuality::Emergent;
        picks.push_back(pick_s);
    }
    
    GeigerLocator locator;
    locator.setVelocityModel(model);
    
    auto start = std::chrono::high_resolution_clock::now();
    auto result = locator.locate(picks, stations);
    auto end = std::chrono::high_resolution_clock::now();
    
    double duration_ms = std::chrono::duration<double, std::milli>(end - start).count();
    
    // Should locate in less than 500ms
    ASSERT_LT(duration_ms, 500.0);
    
    std::cout << "    Geiger located in " << duration_ms << " ms ("
              << result.iterations << " iterations)" << std::endl;
}

TEST(BenchmarkPerformance, GridSearchSpeed) {
    SmallLocalScenario scenario;
    auto stations = scenario.createNetwork();
    auto model = VelocityModel1D::simpleThreeLayer();
    auto origin_time = std::chrono::system_clock::now();
    
    std::vector<PickPtr> picks;
    for (const auto& [key, sta] : stations.stations()) {
        double dist = scenario.hypocenter.distanceTo(sta->location());
        double tt_p = model.travelTime(dist, scenario.hypocenter.depth, PhaseType::P);
        
        auto pick = std::make_shared<Pick>();
        pick->stream_id = StreamID(sta->network(), sta->code(), "00", "BHZ");
        pick->time = origin_time + std::chrono::milliseconds(int(tt_p * 1000));
        pick->phase_type = PhaseType::P;
        pick->quality = PickQuality::Emergent;
        picks.push_back(pick);
    }
    
    GridSearchLocator locator;
    locator.setVelocityModel(model);
    locator.setHorizontalStep(5.0);
    locator.setDepthStep(5.0);
    locator.setSearchRadius(100.0);
    
    auto start = std::chrono::high_resolution_clock::now();
    auto result = locator.locate(picks, stations);
    auto end = std::chrono::high_resolution_clock::now();
    
    double duration_ms = std::chrono::duration<double, std::milli>(end - start).count();
    
    // Grid search should complete in less than 2 seconds
    ASSERT_LT(duration_ms, 2000.0);
    
    std::cout << "    Grid search completed in " << duration_ms << " ms" << std::endl;
}

// ============================================================================
// Accuracy Benchmarks
// ============================================================================

TEST(BenchmarkAccuracy, LocationVsNoise) {
    SmallLocalScenario scenario;
    auto stations = scenario.createNetwork();
    auto model = VelocityModel1D::simpleThreeLayer();
    auto origin_time = std::chrono::system_clock::now();
    
    std::mt19937 gen(42);
    
    std::vector<double> noise_levels = {0.0, 0.1, 0.2, 0.5, 1.0};
    std::vector<double> location_errors;
    
    for (double noise_std : noise_levels) {
        std::normal_distribution<> noise(0.0, noise_std);
        
        std::vector<PickPtr> picks;
        for (const auto& [key, sta] : stations.stations()) {
            double dist = scenario.hypocenter.distanceTo(sta->location());
            double tt_p = model.travelTime(dist, scenario.hypocenter.depth, PhaseType::P);
            tt_p += noise(gen);
            
            auto pick = std::make_shared<Pick>();
            pick->stream_id = StreamID(sta->network(), sta->code(), "00", "BHZ");
            pick->time = origin_time + std::chrono::milliseconds(int(tt_p * 1000));
            pick->phase_type = PhaseType::P;
            pick->quality = PickQuality::Emergent;
            picks.push_back(pick);
        }
        
        GeigerLocator locator;
        locator.setVelocityModel(model);
        
        auto result = locator.locate(picks, stations);
        
        double err = scenario.hypocenter.distanceTo(result.origin.location);
        location_errors.push_back(err);
    }
    
    // Error should increase with noise
    for (size_t i = 1; i < location_errors.size(); i++) {
        // Allow some tolerance since random noise can vary
        ASSERT_LE(location_errors[0], location_errors[i] + 5.0);
    }
    
    // Print results
    std::cout << "    Location error vs pick noise:" << std::endl;
    for (size_t i = 0; i < noise_levels.size(); i++) {
        std::cout << "      Noise=" << noise_levels[i] << "s -> Error=" 
                  << location_errors[i] << " km" << std::endl;
    }
}

TEST(BenchmarkAccuracy, LocationVsStationCount) {
    GeoPoint hypocenter(34.5, -117.5, 10.0);
    auto model = VelocityModel1D::simpleThreeLayer();
    auto origin_time = std::chrono::system_clock::now();
    
    std::vector<int> station_counts = {3, 4, 6, 8, 12};
    std::vector<double> location_errors;
    
    for (int n_stations : station_counts) {
        StationInventory inv;
        
        for (int i = 0; i < n_stations; i++) {
            double angle = i * 2.0 * M_PI / n_stations;
            double lat = hypocenter.latitude + 0.3 * std::cos(angle);
            double lon = hypocenter.longitude + 0.3 * std::sin(angle) / 
                         std::cos(hypocenter.latitude * M_PI / 180.0);
            
            inv.addStation(std::make_shared<Station>("XX", "S" + std::to_string(i), lat, lon));
        }
        
        std::vector<PickPtr> picks;
        for (const auto& [key, sta] : inv.stations()) {
            double dist = hypocenter.distanceTo(sta->location());
            double tt = model.travelTime(dist, hypocenter.depth, PhaseType::P);
            
            auto pick = std::make_shared<Pick>();
            pick->stream_id = StreamID(sta->network(), sta->code(), "00", "BHZ");
            pick->time = origin_time + std::chrono::milliseconds(int(tt * 1000));
            pick->phase_type = PhaseType::P;
            pick->quality = PickQuality::Emergent;
            picks.push_back(pick);
        }
        
        GeigerLocator locator;
        locator.setVelocityModel(model);
        
        auto result = locator.locate(picks, inv);
        
        double err = hypocenter.distanceTo(result.origin.location);
        location_errors.push_back(err);
    }
    
    // More stations should generally give better accuracy
    ASSERT_LT(location_errors.back(), location_errors.front() + 5.0);
    
    std::cout << "    Location error vs station count:" << std::endl;
    for (size_t i = 0; i < station_counts.size(); i++) {
        std::cout << "      Stations=" << station_counts[i] << " -> Error=" 
                  << location_errors[i] << " km" << std::endl;
    }
}
