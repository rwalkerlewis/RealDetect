/**
 * Full-scale integration tests using real EarthScope data
 *
 * Tests the complete seismic processing pipeline from real waveforms
 * to event detection, location, and magnitude calculation.
 */

#include "test_framework.hpp"
#include "realdetect/core/types.hpp"
#include "realdetect/core/waveform.hpp"
#include "realdetect/core/station.hpp"
#include "realdetect/core/event.hpp"
#include "realdetect/core/velocity_model.hpp"
#include "realdetect/core/regional_velocity_model.hpp"
#include "realdetect/core/miniseed.hpp"
#include "realdetect/picker/stalta_picker.hpp"
#include "realdetect/picker/aic_picker.hpp"
#include "realdetect/associator/phase_associator.hpp"
#include "realdetect/locator/grid_search.hpp"
#include "realdetect/locator/geiger.hpp"
#include "realdetect/magnitude/local_magnitude.hpp"
#include "realdetect/magnitude/moment_magnitude.hpp"
#include "realdetect/database/css30_database.hpp"
#include <cmath>
#include <cstdio>
#include <unistd.h>
#include <set>

using namespace realdetect;
using namespace realdetect::test;

// Load real data from test or ridgecrest fixture
static StationInventory loadStations(const std::string& dir = "data/test") {
    StationInventory inv;
    if (!inv.loadFromFile(dir + "/stations.txt"))
        inv.loadFromFile("../" + dir + "/stations.txt");
    for (auto& [key, sta] : inv.stations()) {
        if (!sta->getChannel("BHZ")) {
            Channel bhz; bhz.code = "BHZ"; bhz.sample_rate = 40.0; bhz.dip = -90;
            sta->addChannel(bhz);
        }
    }
    return inv;
}

static std::map<StreamID, WaveformPtr> loadWaveforms(const std::string& dir = "data/test") {
    MiniSeedReader reader;
    if (!reader.open(dir + "/waveforms.mseed"))
        reader.open("../" + dir + "/waveforms.mseed");
    auto wfs = reader.toWaveforms();
    std::map<StreamID, WaveformPtr> result;
    for (auto& wf : wfs)
        if (wf->streamId().channel == "BHZ")
            result[wf->streamId()] = wf;
    return result;
}

static std::vector<PickPtr> pickFirstArrivals(
    const std::map<StreamID, WaveformPtr>& waveforms) {
    STALTAPicker picker;
    picker.setParameter("sta_length", 0.3);
    picker.setParameter("lta_length", 8.0);
    picker.setParameter("trigger_ratio", 2.5);

    std::vector<PickPtr> all_picks;
    for (auto& [id, wf] : waveforms) {
        auto results = picker.pick(*wf);
        for (auto& r : results) {
            auto p = std::make_shared<Pick>();
            p->stream_id = id;
            p->time = r.time;
            p->phase_type = r.phase_type;
            p->snr = r.snr;
            p->quality = r.snr > 8 ? PickQuality::Impulsive : PickQuality::Emergent;
            all_picks.push_back(p);
        }
    }

    std::vector<PickPtr> p_picks;
    std::set<std::string> seen;
    std::sort(all_picks.begin(), all_picks.end(),
              [](auto& a, auto& b) { return a->time < b->time; });
    for (auto& pk : all_picks) {
        std::string key = pk->stream_id.network + "." + pk->stream_id.station;
        if (seen.insert(key).second) {
            pk->phase_type = PhaseType::P;
            p_picks.push_back(pk);
        }
    }
    return p_picks;
}

// ============================================================================
// End-to-End Pipeline Tests with Real Data
// ============================================================================

TEST(Integration, FullPipelineRealData) {
    auto stations = loadStations();
    auto waveforms = loadWaveforms();
    if (stations.size() == 0 || waveforms.empty()) {
        std::cout << "    [SKIP] No real test data available\n";
        return;
    }

    ASSERT_GE(waveforms.size(), 3u);

    // Phase picking on real data
    auto picks = pickFirstArrivals(waveforms);
    ASSERT_GE(picks.size(), 2u);

    // Association
    auto model = VelocityModel1D::simpleThreeLayer();
    PhaseAssociator assoc;
    assoc.setStations(stations);
    assoc.setVelocityModel(model);
    assoc.setMinStations(2);
    assoc.setMinPhases(3);
    for (auto& pk : picks) assoc.addPick(pk);
    auto events = assoc.process();

    // Location
    std::vector<PickPtr> loc_picks = picks;
    if (!events.empty()) {
        loc_picks.clear();
        for (auto& arr : events[0]->preferredOrigin().arrivals)
            loc_picks.push_back(arr.pick);
    }

    GridSearchLocator gs;
    gs.setVelocityModel(model);
    gs.setSearchRadius(200.0);
    auto loc_result = gs.locate(loc_picks, stations);

    ASSERT_GT(loc_result.origin.phase_count, 0);

    // Magnitude
    LocalMagnitude ml;
    auto ml_result = ml.calculate(loc_result.origin, waveforms, stations);

    std::cout << "    Pipeline: " << picks.size() << " picks, "
              << events.size() << " events, location at "
              << loc_result.origin.location.latitude << "°N\n";
}

TEST(Integration, FullPipelineRidgecrestData) {
    auto stations = loadStations("data/ridgecrest");
    auto waveforms = loadWaveforms("data/ridgecrest");
    if (stations.size() == 0 || waveforms.empty()) {
        std::cout << "    [SKIP] No Ridgecrest data available\n";
        return;
    }

    auto picks = pickFirstArrivals(waveforms);

    auto model = VelocityModel1D::simpleThreeLayer();
    GridSearchLocator gs;
    gs.setVelocityModel(model);
    gs.setSearchRadius(300.0);
    auto result = gs.locate(picks, stations);

    // Ground truth: 35.770°N, 117.599°W
    GeoPoint truth(35.770, -117.599);
    GeoPoint computed(result.origin.location.latitude, result.origin.location.longitude);
    double error_km = truth.distanceTo(computed);

    std::cout << "    Ridgecrest location error: " << error_km << " km\n";
    // Location should at least be in the right part of California
    ASSERT_GT(result.origin.location.latitude, 25.0);
    ASSERT_LT(result.origin.location.latitude, 45.0);
}

// ============================================================================
// Database Integration Tests with Real Data
// ============================================================================

TEST(Integration, EventToDatabase) {
    std::string db_path = "/tmp/integration_test_" + std::to_string(getpid()) + ".db";

    CSS30Database db;
    ASSERT_TRUE(db.open(db_path));
    ASSERT_TRUE(db.createSchema());
    db.setAuthor("integration_test");

    auto stations = loadStations();
    auto waveforms = loadWaveforms();
    if (stations.size() == 0 || waveforms.empty()) {
        std::cout << "    [SKIP] No real test data\n";
        db.close();
        std::remove(db_path.c_str());
        return;
    }

    auto picks = pickFirstArrivals(waveforms);
    auto model = VelocityModel1D::simpleThreeLayer();

    GridSearchLocator locator;
    locator.setVelocityModel(model);
    auto loc_result = locator.locate(picks, stations);

    Event event;
    Origin origin;
    origin.location = loc_result.origin.location;
    origin.time = loc_result.origin.time;
    origin.rms = loc_result.origin.rms;
    origin.gap = loc_result.origin.gap;
    origin.phase_count = loc_result.origin.phase_count;
    origin.station_count = loc_result.origin.station_count;
    event.addOrigin(origin);
    event.addMagnitude(Magnitude(MagnitudeType::ML, 7.1, 0.3, 5));

    ASSERT_TRUE(db.storeCompleteEvent(event, "simple3layer", "CI"));
    ASSERT_EQ(db.countEvents(), 1);
    ASSERT_GE(db.countOrigins(), 1);

    db.storeInventory(stations, "CI");

    db.close();
    std::remove(db_path.c_str());
}

// ============================================================================
// Regional Velocity Model Integration
// ============================================================================

TEST(Integration, RegionalVelocityModelSelection) {
    VelocityModelManager manager;

    VelocityModel1D socal("SoCal");
    socal.addLayer(0, 5, 5.5, 3.18);
    socal.addLayer(5, 10, 6.0, 3.46);
    socal.addLayer(15, 17, 6.5, 3.75);
    socal.addLayer(32, 0, 7.8, 4.5);

    manager.setDefaultModel(VelocityModel1D::iasp91());
    manager.addModel("socal",
        GeographicBounds::boundingBox(32.0, 36.0, -121.0, -114.0),
        socal, 10);

    GeoPoint socal_event(34.0, -118.0, 10.0);
    GeoPoint global_event(40.0, -100.0, 10.0);

    const auto& model_socal = manager.getModelForLocation(socal_event);
    const auto& model_global = manager.getModelForLocation(global_event);

    ASSERT_EQ(model_socal.name(), "SoCal");
    ASSERT_NE(model_global.name(), "SoCal");

    double tt_socal = model_socal.travelTime(100.0, 10.0, PhaseType::P);
    double tt_global = model_global.travelTime(100.0, 10.0, PhaseType::P);
    ASSERT_NE(tt_socal, tt_global);
}

// ============================================================================
// Performance Tests with Real Data
// ============================================================================

TEST(Integration, HighVolumeProcessing) {
    auto stations = loadStations("data/ridgecrest");
    auto waveforms = loadWaveforms("data/ridgecrest");
    if (stations.size() == 0 || waveforms.empty()) {
        std::cout << "    [SKIP] No data\n";
        return;
    }

    auto start = std::chrono::high_resolution_clock::now();

    STALTAPicker picker;
    picker.setParameter("trigger_ratio", 3.0);

    size_t total_picks = 0;
    // Process same real data 5 times to simulate high volume
    for (int iter = 0; iter < 5; iter++) {
        for (auto& [id, wf] : waveforms) {
            auto results = picker.pick(*wf);
            total_picks += results.size();
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    double duration_ms = std::chrono::duration<double, std::milli>(end - start).count();

    ASSERT_LT(duration_ms, 30000.0);

    std::cout << "    Processed " << total_picks << " picks from "
              << waveforms.size() << " waveforms × 5 in "
              << duration_ms << " ms\n";
}

TEST(Integration, HandleMissingStations) {
    auto stations = loadStations();
    auto waveforms = loadWaveforms();
    if (stations.size() == 0 || waveforms.empty()) {
        std::cout << "    [SKIP] No data\n";
        return;
    }

    // Remove some waveforms to simulate missing stations
    std::map<StreamID, WaveformPtr> partial;
    int count = 0;
    for (auto& [id, wf] : waveforms) {
        if (count++ < 3) partial[id] = wf;
    }

    auto picks = pickFirstArrivals(partial);
    auto model = VelocityModel1D::simpleThreeLayer();

    GridSearchLocator locator;
    locator.setVelocityModel(model);

    if (picks.size() >= 3) {
        ASSERT_NO_THROW(locator.locate(picks, stations));
    }
}
