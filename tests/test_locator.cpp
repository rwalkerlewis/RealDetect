/**
 * Unit tests for hypocenter location algorithms
 *
 * Uses real seismic data from EarthScope for integration-style tests.
 */

#include "test_framework.hpp"
#include "realdetect/core/station.hpp"
#include "realdetect/core/event.hpp"
#include "realdetect/core/velocity_model.hpp"
#include "realdetect/core/miniseed.hpp"
#include "realdetect/locator/grid_search.hpp"
#include "realdetect/locator/geiger.hpp"
#include "realdetect/associator/phase_associator.hpp"
#include "realdetect/picker/stalta_picker.hpp"
#include <set>

using namespace realdetect;
using namespace realdetect::test;

// Ridgecrest ground truth
static const double TRUTH_LAT = 35.770;
static const double TRUTH_LON = -117.599;
static const double TRUTH_DEPTH = 8.0;

// Load real test fixtures
static StationInventory loadTestStations() {
    StationInventory inv;
    if (!inv.loadFromFile("data/test/stations.txt"))
        inv.loadFromFile("../data/test/stations.txt");
    for (auto& [key, sta] : inv.stations()) {
        if (!sta->getChannel("BHZ")) {
            Channel bhz; bhz.code = "BHZ"; bhz.sample_rate = 40.0; bhz.dip = -90;
            sta->addChannel(bhz);
        }
    }
    return inv;
}

static std::vector<WaveformPtr> loadTestWaveforms() {
    MiniSeedReader reader;
    if (!reader.open("data/test/waveforms.mseed"))
        reader.open("../data/test/waveforms.mseed");
    auto wfs = reader.toWaveforms();
    std::vector<WaveformPtr> bhz;
    for (auto& wf : wfs)
        if (wf->streamId().channel == "BHZ") bhz.push_back(wf);
    return bhz;
}

// Pick P-wave first arrivals from real data
static std::vector<PickPtr> pickRealData(const std::vector<WaveformPtr>& waveforms) {
    STALTAPicker picker;
    picker.setParameter("sta_length", 0.3);
    picker.setParameter("lta_length", 8.0);
    picker.setParameter("trigger_ratio", 2.5);

    std::vector<PickPtr> all_picks;
    for (auto& wf : waveforms) {
        auto results = picker.pick(*wf);
        for (auto& r : results) {
            auto p = std::make_shared<Pick>();
            p->stream_id = wf->streamId();
            p->time = r.time;
            p->phase_type = r.phase_type;
            p->snr = r.snr;
            p->quality = r.snr > 8 ? PickQuality::Impulsive : PickQuality::Emergent;
            all_picks.push_back(p);
        }
    }

    // Keep first arrival per station
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
// Travel Time Tests
// ============================================================================

TEST(TravelTimeTable, BasicTravelTime) {
    TravelTimeTable tt;
    auto model = VelocityModel1D::simpleThreeLayer();
    tt.initialize(model);

    double p_time = tt.getTime(100.0, 10.0, PhaseType::P);
    double s_time = tt.getTime(100.0, 10.0, PhaseType::S);
    ASSERT_GT(p_time, 0);
    ASSERT_GT(s_time, p_time);
}

TEST(TravelTimeTable, TravelTimeIncreasesWithDistance) {
    TravelTimeTable tt;
    auto model = VelocityModel1D::simpleThreeLayer();
    tt.initialize(model);

    double t1 = tt.getTime(50.0, 10.0, PhaseType::P);
    double t2 = tt.getTime(100.0, 10.0, PhaseType::P);
    double t3 = tt.getTime(200.0, 10.0, PhaseType::P);
    ASSERT_GT(t2, t1);
    ASSERT_GT(t3, t2);
}

TEST(TravelTimeTable, ZeroDistanceTime) {
    TravelTimeTable tt;
    auto model = VelocityModel1D::simpleThreeLayer();
    tt.initialize(model);
    double t = tt.getTime(0.0, 10.0, PhaseType::P);
    ASSERT_GT(t, 0);  // Still positive (depth travel time)
}

// ============================================================================
// Grid Search Locator Tests
// ============================================================================

TEST(GridSearchLocator, LocateRealData) {
    auto stations = loadTestStations();
    auto waveforms = loadTestWaveforms();
    if (stations.size() == 0 || waveforms.empty()) {
        std::cout << "    [SKIP] No real test data available\n";
        return;
    }

    auto picks = pickRealData(waveforms);
    if (picks.size() < 3) {
        std::cout << "    [SKIP] Insufficient picks from real data\n";
        return;
    }

    auto model = VelocityModel1D::simpleThreeLayer();
    GridSearchLocator gs;
    gs.setVelocityModel(model);
    gs.setHorizontalStep(5.0);
    gs.setSearchRadius(200.0);

    auto result = gs.locate(picks, stations);
    // Should produce a result (may not be perfect with real data)
    ASSERT_GT(result.origin.phase_count, 0);
    // Location should be on Earth
    ASSERT_GT(result.origin.location.latitude, -90.0);
    ASSERT_LT(result.origin.location.latitude, 90.0);
}

TEST(GridSearchLocator, SetParameters) {
    GridSearchLocator gs;
    gs.setHorizontalStep(2.0);
    gs.setDepthStep(3.0);
    gs.setSearchRadius(150.0);

    auto model = VelocityModel1D::simpleThreeLayer();
    gs.setVelocityModel(model);
    // Just verify no crash
    ASSERT_TRUE(true);
}

// ============================================================================
// Geiger Locator Tests
// ============================================================================

TEST(GeigerLocator, LocateRealData) {
    auto stations = loadTestStations();
    auto waveforms = loadTestWaveforms();
    if (stations.size() == 0 || waveforms.empty()) {
        std::cout << "    [SKIP] No real test data available\n";
        return;
    }

    auto picks = pickRealData(waveforms);
    if (picks.size() < 3) {
        std::cout << "    [SKIP] Insufficient picks from real data\n";
        return;
    }

    auto model = VelocityModel1D::simpleThreeLayer();

    // First do grid search for initial location
    GridSearchLocator gs;
    gs.setVelocityModel(model);
    gs.setSearchRadius(200.0);
    auto gs_result = gs.locate(picks, stations);

    // Then refine with Geiger
    GeigerLocator geiger;
    geiger.setVelocityModel(model);
    geiger.setInitialLocation(gs_result.origin.location);
    geiger.setInitialTime(gs_result.origin.time);
    geiger.setMaxIterations(30);
    geiger.setDampingFactor(0.4);

    auto result = geiger.locate(picks, stations);
    // Verify it produces a location
    ASSERT_GT(result.origin.location.latitude, 0.0);  // Northern hemisphere
}

TEST(GeigerLocator, ParameterSetting) {
    GeigerLocator geiger;
    geiger.setMaxIterations(50);
    geiger.setConvergenceThreshold(0.01);
    geiger.setDampingFactor(0.3);

    auto model = VelocityModel1D::simpleThreeLayer();
    geiger.setVelocityModel(model);
    ASSERT_TRUE(true);
}

// ============================================================================
// Locator Factory Tests
// ============================================================================

TEST(LocatorFactory, CreateGridSearch) {
    auto gs = LocatorFactory::create("grid");
    ASSERT_NE(gs, nullptr);
}

TEST(LocatorFactory, CreateGeiger) {
    auto geiger = LocatorFactory::create("geiger");
    ASSERT_NE(geiger, nullptr);
}

TEST(LocatorFactory, CreateOctTree) {
    auto oct = LocatorFactory::create("octree");
    ASSERT_NE(oct, nullptr);
}
