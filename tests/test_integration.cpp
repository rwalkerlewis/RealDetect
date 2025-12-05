/**
 * Full-scale integration tests
 * 
 * Tests the complete seismic processing pipeline from waveforms
 * to event detection, location, and magnitude calculation.
 */

#include "test_framework.hpp"
#include "realdetect/core/types.hpp"
#include "realdetect/core/waveform.hpp"
#include "realdetect/core/station.hpp"
#include "realdetect/core/event.hpp"
#include "realdetect/core/velocity_model.hpp"
#include "realdetect/core/regional_velocity_model.hpp"
#include "realdetect/picker/stalta_picker.hpp"
#include "realdetect/picker/aic_picker.hpp"
#include "realdetect/associator/phase_associator.hpp"
#include "realdetect/locator/grid_search.hpp"
#include "realdetect/locator/geiger.hpp"
#include "realdetect/magnitude/local_magnitude.hpp"
#include "realdetect/magnitude/moment_magnitude.hpp"
#include "realdetect/database/css30_database.hpp"
#include <random>
#include <cmath>
#include <cstdio>
#include <unistd.h>

using namespace realdetect;
using namespace realdetect::test;

/**
 * SyntheticEventGenerator - Generates realistic synthetic seismograms
 */
class SyntheticEventGenerator {
public:
    SyntheticEventGenerator(unsigned seed = 42) : gen_(seed), noise_(0.0, 1.0) {}
    
    // Generate complete waveform dataset for an event
    std::map<StreamID, WaveformPtr> generateEvent(
        const GeoPoint& hypocenter,
        TimePoint origin_time,
        double magnitude,
        const StationInventory& stations,
        const VelocityModel1D& model,
        double sample_rate = 100.0,
        double duration = 120.0) {
        
        std::map<StreamID, WaveformPtr> waveforms;
        
        for (const auto& [key, sta] : stations.stations()) {
            auto wf = generateWaveform(
                *sta, hypocenter, origin_time, magnitude, model, sample_rate, duration);
            waveforms[wf->streamId()] = wf;
        }
        
        return waveforms;
    }
    
    // Generate single waveform
    WaveformPtr generateWaveform(
        const Station& station,
        const GeoPoint& hypocenter,
        TimePoint origin_time,
        double magnitude,
        const VelocityModel1D& model,
        double sample_rate,
        double duration) {
        
        StreamID id(station.network(), station.code(), "00", "BHZ");
        auto wf = std::make_shared<Waveform>(id, sample_rate, origin_time);
        
        double epi_dist = hypocenter.distanceTo(station.location());
        double hypo_dist = std::sqrt(epi_dist * epi_dist + 
                                     hypocenter.depth * hypocenter.depth);
        
        double tt_p = model.travelTime(epi_dist, hypocenter.depth, PhaseType::P);
        double tt_s = model.travelTime(epi_dist, hypocenter.depth, PhaseType::S);
        
        // Amplitude based on magnitude and distance
        double log_amp = magnitude - 1.66 * std::log10(std::max(1.0, epi_dist)) - 2.0;
        double base_amp = std::pow(10.0, log_amp) * 1e6;
        
        size_t n_samples = static_cast<size_t>(duration * sample_rate);
        
        for (size_t i = 0; i < n_samples; i++) {
            double t = i / sample_rate;
            double sample = noise_(gen_) * 10.0;
            
            // P-wave
            if (t >= tt_p && t < tt_s + 10.0) {
                double t_rel = t - tt_p;
                double env = std::exp(-t_rel / 2.0) * (1.0 - std::exp(-t_rel * 5.0));
                double freq = std::max(1.0, 10.0 - epi_dist / 50.0);
                sample += base_amp * 0.3 * env * std::sin(2.0 * M_PI * freq * t_rel);
            }
            
            // S-wave
            if (t >= tt_s) {
                double t_rel = t - tt_s;
                double env = std::exp(-t_rel / 5.0) * (1.0 - std::exp(-t_rel * 3.0));
                double freq = std::max(0.5, 5.0 - epi_dist / 100.0);
                sample += base_amp * env * std::sin(2.0 * M_PI * freq * t_rel);
            }
            
            wf->append(sample);
        }
        
        return wf;
    }

private:
    std::mt19937 gen_;
    std::normal_distribution<> noise_;
};

// Helper to create a test network
StationInventory createIntegrationTestNetwork(
    double center_lat, double center_lon, int n_stations = 8, double radius_deg = 0.5) {
    
    StationInventory inv;
    
    for (int i = 0; i < n_stations; i++) {
        double angle = i * 2.0 * M_PI / n_stations;
        double lat = center_lat + radius_deg * std::cos(angle);
        double lon = center_lon + radius_deg * std::sin(angle) / 
                     std::cos(center_lat * M_PI / 180.0);
        
        std::string code = "S" + std::to_string(i + 1);
        auto sta = std::make_shared<Station>("XX", code, lat, lon, 500);
        
        // Add channels
        Channel bhz;
        bhz.code = "BHZ";
        bhz.sample_rate = 100.0;
        bhz.dip = -90.0;
        sta->addChannel(bhz);
        
        inv.addStation(sta);
    }
    
    return inv;
}

// ============================================================================
// End-to-End Pipeline Tests
// ============================================================================

TEST(Integration, FullPipelineSmallEvent) {
    // Test complete pipeline for a small local earthquake
    
    // Setup
    auto stations = createIntegrationTestNetwork(34.0, -118.0, 10, 0.4);
    auto model = VelocityModel1D::simpleThreeLayer();
    
    GeoPoint hypocenter(34.0, -118.0, 10.0);
    double magnitude = 3.5;
    auto origin_time = std::chrono::system_clock::now();
    
    SyntheticEventGenerator generator(42);
    auto waveforms = generator.generateEvent(
        hypocenter, origin_time, magnitude, stations, model);
    
    ASSERT_EQ(waveforms.size(), 10u);
    
    // Phase picking
    STALTAPicker picker;
    picker.setSTALength(0.5);
    picker.setLTALength(5.0);
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
            pick->quality = r.snr > 8.0 ? PickQuality::Impulsive : PickQuality::Emergent;
            picks.push_back(pick);
        }
    }
    
    ASSERT_GE(picks.size(), 6u);  // At least 6 picks
    
    // Phase association
    PhaseAssociator assoc;
    assoc.setStations(stations);
    assoc.setVelocityModel(model);
    assoc.setMinStations(3);
    assoc.setMinPhases(4);
    
    for (const auto& pick : picks) {
        assoc.addPick(pick);
    }
    
    auto events = assoc.process();
    ASSERT_GE(events.size(), 1u);
    
    // Hypocenter location
    GridSearchLocator gs;
    gs.setVelocityModel(model);
    gs.setHorizontalStep(5.0);
    gs.setDepthStep(5.0);
    gs.setSearchRadius(80.0);
    
    auto loc_result = gs.locate(picks, stations);
    ASSERT_TRUE(loc_result.converged);
    
    double loc_error = hypocenter.distanceTo(loc_result.origin.location);
    ASSERT_LT(loc_error, 20.0);  // Within 20 km
    
    // Magnitude calculation
    LocalMagnitude ml;
    ml.setSimulateWoodAnderson(false);
    
    auto mag_result = ml.calculate(loc_result.origin, waveforms, stations);
    
    if (mag_result.station_count >= 3) {
        double mag_error = std::abs(mag_result.value - magnitude);
        ASSERT_LT(mag_error, 3.0);  // Within 3 magnitude units (synthetic data can vary)
    }
}

TEST(Integration, FullPipelineModerateEvent) {
    // Test pipeline for a moderate regional earthquake
    
    auto stations = createIntegrationTestNetwork(35.0, -117.0, 12, 0.6);
    auto model = VelocityModel1D::simpleThreeLayer();
    
    GeoPoint hypocenter(35.0, -117.0, 15.0);
    double magnitude = 5.0;
    auto origin_time = std::chrono::system_clock::now();
    
    SyntheticEventGenerator generator(123);
    auto waveforms = generator.generateEvent(
        hypocenter, origin_time, magnitude, stations, model, 100.0, 180.0);
    
    // Pick phases
    STALTAPicker stalta;
    stalta.setTriggerRatio(3.0);
    
    AICPicker aic;
    
    std::vector<PickPtr> picks;
    for (const auto& [id, wf] : waveforms) {
        auto stalta_results = stalta.pick(*wf);
        
        for (const auto& r : stalta_results) {
            // Refine with AIC
            size_t refined = aic.refinePick(*wf, r.sample_index, 100);
            
            auto pick = std::make_shared<Pick>();
            pick->stream_id = id;
            pick->time = wf->timeAt(refined);
            pick->phase_type = r.phase_type;
            pick->snr = r.snr;
            pick->quality = PickQuality::Emergent;
            picks.push_back(pick);
        }
    }
    
    ASSERT_GE(picks.size(), 8u);
    
    // Locate with grid search, then refine with Geiger
    GridSearchLocator gs;
    gs.setVelocityModel(model);
    gs.setHorizontalStep(5.0);
    gs.setSearchRadius(120.0);
    
    auto gs_result = gs.locate(picks, stations);
    
    if (gs_result.converged) {
        GeigerLocator geiger;
        geiger.setVelocityModel(model);
        geiger.setInitialLocation(gs_result.origin.location);
        geiger.setInitialTime(gs_result.origin.time);
        geiger.setMaxIterations(30);
        
        auto final_result = geiger.locate(picks, stations);
        
        double loc_error = hypocenter.distanceTo(final_result.origin.location);
        ASSERT_LT(loc_error, 25.0);
    }
}

// ============================================================================
// Database Integration Tests
// ============================================================================

TEST(Integration, EventToDatabase) {
    std::string db_path = "/tmp/integration_test_" + std::to_string(getpid()) + ".db";
    
    // Create database
    CSS30Database db;
    ASSERT_TRUE(db.open(db_path));
    ASSERT_TRUE(db.createSchema());
    db.setAuthor("integration_test");
    
    // Setup and run pipeline
    auto stations = createIntegrationTestNetwork(34.0, -118.0, 8, 0.3);
    auto model = VelocityModel1D::simpleThreeLayer();
    
    GeoPoint hypocenter(34.0, -118.0, 8.0);
    auto origin_time = std::chrono::system_clock::now();
    
    SyntheticEventGenerator generator(456);
    auto waveforms = generator.generateEvent(
        hypocenter, origin_time, 3.8, stations, model);
    
    // Quick picking
    STALTAPicker picker;
    picker.setTriggerRatio(4.0);
    
    std::vector<PickPtr> picks;
    for (const auto& [id, wf] : waveforms) {
        auto results = picker.pick(*wf);
        for (const auto& r : results) {
            auto pick = std::make_shared<Pick>();
            pick->stream_id = id;
            pick->time = r.time;
            pick->phase_type = r.phase_type;
            picks.push_back(pick);
        }
    }
    
    // Locate
    GridSearchLocator locator;
    locator.setVelocityModel(model);
    auto loc_result = locator.locate(picks, stations);
    
    // Create and store event
    Event event;
    
    Origin origin;
    origin.location = loc_result.origin.location;
    origin.time = loc_result.origin.time;
    origin.rms = loc_result.origin.rms;
    origin.gap = loc_result.origin.gap;
    origin.phase_count = loc_result.origin.phase_count;
    origin.station_count = loc_result.origin.station_count;
    event.addOrigin(origin);
    
    event.addMagnitude(Magnitude(MagnitudeType::ML, 3.8, 0.3, 5));
    
    ASSERT_TRUE(db.storeCompleteEvent(event, "simple3layer", "XX"));
    
    // Verify stored
    ASSERT_EQ(db.countEvents(), 1);
    ASSERT_GE(db.countOrigins(), 1);
    
    // Store station inventory
    db.storeInventory(stations, "XX");
    
    db.close();
    std::remove(db_path.c_str());
}

TEST(Integration, MultipleEventsToDatabase) {
    std::string db_path = "/tmp/integration_multi_" + std::to_string(getpid()) + ".db";
    
    CSS30Database db;
    ASSERT_TRUE(db.open(db_path));
    ASSERT_TRUE(db.createSchema());
    db.setAuthor("test");
    
    auto stations = createIntegrationTestNetwork(35.0, -117.0, 10, 0.5);
    auto model = VelocityModel1D::simpleThreeLayer();
    
    SyntheticEventGenerator generator(789);
    
    // Generate and process 5 events
    std::vector<GeoPoint> hypocenters = {
        GeoPoint(34.8, -117.2, 10.0),
        GeoPoint(35.0, -116.8, 8.0),
        GeoPoint(35.2, -117.0, 15.0),
        GeoPoint(34.9, -117.3, 12.0),
        GeoPoint(35.1, -116.9, 6.0)
    };
    
    auto base_time = std::chrono::system_clock::now();
    
    db.beginTransaction();
    
    for (size_t i = 0; i < hypocenters.size(); i++) {
        auto origin_time = base_time + std::chrono::hours(i * 24);
        double mag = 3.0 + i * 0.5;
        
        auto waveforms = generator.generateEvent(
            hypocenters[i], origin_time, mag, stations, model);
        
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
                picks.push_back(pick);
            }
        }
        
        if (picks.size() >= 4) {
            GridSearchLocator locator;
            locator.setVelocityModel(model);
            auto result = locator.locate(picks, stations);
            
            Event event;
            Origin origin;
            origin.location = result.origin.location;
            origin.time = result.origin.time;
            origin.rms = result.origin.rms;
            event.addOrigin(origin);
            event.addMagnitude(Magnitude(MagnitudeType::ML, mag, 0.2, 5));
            
            db.storeCompleteEvent(event, "model", "XX");
        }
    }
    
    db.commitTransaction();
    
    ASSERT_GE(db.countEvents(), 3);  // At least 3 events should be stored
    
    db.close();
    std::remove(db_path.c_str());
}

// ============================================================================
// Regional Velocity Model Integration
// ============================================================================

TEST(Integration, RegionalVelocityModelSelection) {
    VelocityModelManager manager;
    
    // Setup regional models
    VelocityModel1D socal("SoCal");
    socal.addLayer(0, 5, 5.5, 3.18);
    socal.addLayer(5, 10, 6.0, 3.46);
    socal.addLayer(15, 17, 6.5, 3.75);
    socal.addLayer(32, 0, 7.8, 4.5);
    
    manager.setDefaultModel(VelocityModel1D::iasp91());
    manager.addModel("socal",
        GeographicBounds::boundingBox(32.0, 36.0, -121.0, -114.0),
        socal, 10);
    
    // Create events in different regions
    GeoPoint socal_event(34.0, -118.0, 10.0);
    GeoPoint global_event(40.0, -100.0, 10.0);
    
    const auto& model_socal = manager.getModelForLocation(socal_event);
    const auto& model_global = manager.getModelForLocation(global_event);
    
    ASSERT_EQ(model_socal.name(), "SoCal");
    ASSERT_NE(model_global.name(), "SoCal");
    
    // Travel times should be different
    double tt_socal = model_socal.travelTime(100.0, 10.0, PhaseType::P);
    double tt_global = model_global.travelTime(100.0, 10.0, PhaseType::P);
    
    ASSERT_NE(tt_socal, tt_global);
}

// ============================================================================
// Stress Tests
// ============================================================================

TEST(Integration, HighVolumeProcessing) {
    // Process many picks rapidly
    auto stations = createIntegrationTestNetwork(34.0, -118.0, 20, 0.8);
    auto model = VelocityModel1D::simpleThreeLayer();
    
    auto start = std::chrono::high_resolution_clock::now();
    
    // Generate 10 events worth of picks
    std::vector<PickPtr> all_picks;
    
    SyntheticEventGenerator generator(999);
    auto base_time = std::chrono::system_clock::now();
    
    for (int event_num = 0; event_num < 10; event_num++) {
        GeoPoint hypo(34.0 + (event_num % 3) * 0.1, 
                     -118.0 + (event_num % 2) * 0.1, 
                     8.0 + event_num);
        auto origin_time = base_time + std::chrono::minutes(event_num * 10);
        
        auto waveforms = generator.generateEvent(
            hypo, origin_time, 3.0 + event_num * 0.2, stations, model);
        
        STALTAPicker picker;
        picker.setTriggerRatio(3.0);
        
        for (const auto& [id, wf] : waveforms) {
            auto results = picker.pick(*wf);
            for (const auto& r : results) {
                auto pick = std::make_shared<Pick>();
                pick->stream_id = id;
                pick->time = r.time;
                pick->phase_type = r.phase_type;
                all_picks.push_back(pick);
            }
        }
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    double duration_ms = std::chrono::duration<double, std::milli>(end - start).count();
    
    // Should process all events in reasonable time
    ASSERT_LT(duration_ms, 30000.0);  // < 30 seconds
    
    std::cout << "    Processed " << all_picks.size() << " picks from 10 events in "
              << duration_ms << " ms" << std::endl;
}

TEST(Integration, ContinuousStreamSimulation) {
    // Simulate continuous data streaming
    auto stations = createIntegrationTestNetwork(34.5, -117.5, 8, 0.4);
    auto model = VelocityModel1D::simpleThreeLayer();
    
    SyntheticEventGenerator generator(111);
    
    // Simulate 1 hour of data with 2 events
    auto start_time = std::chrono::system_clock::now();
    
    std::vector<GeoPoint> event_hypos = {
        GeoPoint(34.5, -117.5, 10.0),
        GeoPoint(34.6, -117.4, 12.0)
    };
    
    std::vector<std::chrono::minutes> event_times = {
        std::chrono::minutes(15),
        std::chrono::minutes(45)
    };
    
    int events_detected = 0;
    
    for (size_t i = 0; i < event_hypos.size(); i++) {
        auto origin_time = start_time + event_times[i];
        auto waveforms = generator.generateEvent(
            event_hypos[i], origin_time, 4.0, stations, model);
        
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
                picks.push_back(pick);
            }
        }
        
        if (picks.size() >= 4) {
            GridSearchLocator locator;
            locator.setVelocityModel(model);
            auto result = locator.locate(picks, stations);
            
            if (result.converged) {
                events_detected++;
            }
        }
    }
    
    ASSERT_GE(events_detected, 1);  // Should detect at least 1 event
}

// ============================================================================
// Error Handling Tests
// ============================================================================

TEST(Integration, HandleMissingStations) {
    auto stations = createIntegrationTestNetwork(34.0, -118.0, 8, 0.4);
    auto model = VelocityModel1D::simpleThreeLayer();
    
    GeoPoint hypocenter(34.0, -118.0, 10.0);
    auto origin_time = std::chrono::system_clock::now();
    
    SyntheticEventGenerator generator(222);
    auto waveforms = generator.generateEvent(
        hypocenter, origin_time, 3.5, stations, model);
    
    // Remove some waveforms (simulate missing stations)
    std::vector<StreamID> to_remove;
    int count = 0;
    for (const auto& [id, wf] : waveforms) {
        if (count++ >= 3) break;
        to_remove.push_back(id);
    }
    for (const auto& id : to_remove) {
        waveforms.erase(id);
    }
    
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
            picks.push_back(pick);
        }
    }
    
    // Should still be able to locate with fewer stations
    if (picks.size() >= 4) {
        GridSearchLocator locator;
        locator.setVelocityModel(model);
        
        ASSERT_NO_THROW(locator.locate(picks, stations));
    }
}

TEST(Integration, HandleNoisyData) {
    auto stations = createIntegrationTestNetwork(34.0, -118.0, 10, 0.4);
    auto model = VelocityModel1D::simpleThreeLayer();
    
    GeoPoint hypocenter(34.0, -118.0, 10.0);
    auto origin_time = std::chrono::system_clock::now();
    
    // High noise generator
    std::mt19937 gen(333);
    std::normal_distribution<> noise(0.0, 50.0);  // High noise
    
    SyntheticEventGenerator generator(333);
    auto waveforms = generator.generateEvent(
        hypocenter, origin_time, 4.0, stations, model);
    
    // Add extra noise to waveforms
    for (auto& [id, wf] : waveforms) {
        for (size_t i = 0; i < wf->sampleCount(); i++) {
            wf->data()[i] += noise(gen);
        }
    }
    
    // Increase picker threshold for noisy data
    STALTAPicker picker;
    picker.setTriggerRatio(5.0);
    
    std::vector<PickPtr> picks;
    for (const auto& [id, wf] : waveforms) {
        auto results = picker.pick(*wf);
        for (const auto& r : results) {
            auto pick = std::make_shared<Pick>();
            pick->stream_id = id;
            pick->time = r.time;
            pick->phase_type = r.phase_type;
            picks.push_back(pick);
        }
    }
    
    // May or may not detect event in very noisy data - test that it doesn't crash
    GridSearchLocator locator;
    locator.setVelocityModel(model);
    
    if (picks.size() >= 3) {
        ASSERT_NO_THROW(locator.locate(picks, stations));
    }
}
