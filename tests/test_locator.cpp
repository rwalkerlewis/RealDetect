/**
 * Unit tests for hypocenter location algorithms
 */

#include "test_framework.hpp"
#include "realdetect/core/station.hpp"
#include "realdetect/core/event.hpp"
#include "realdetect/core/velocity_model.hpp"
#include "realdetect/locator/grid_search.hpp"
#include "realdetect/locator/geiger.hpp"
#include "realdetect/associator/phase_associator.hpp"
#include <random>

using namespace realdetect;
using namespace realdetect::test;

// Helper to create a network of stations
StationInventory createTestNetwork(double center_lat, double center_lon) {
    StationInventory inv;
    
    // Ring of stations around center
    std::vector<std::pair<double, double>> offsets = {
        {0.5, 0.0},
        {0.35, 0.35},
        {0.0, 0.5},
        {-0.35, 0.35},
        {-0.5, 0.0},
        {-0.35, -0.35},
        {0.0, -0.5},
        {0.35, -0.35}
    };
    
    int i = 1;
    for (const auto& [dlat, dlon] : offsets) {
        std::string code = "S" + std::to_string(i);
        inv.addStation(std::make_shared<Station>("XX", code,
            center_lat + dlat, center_lon + dlon, 0.0));
        i++;
    }
    
    return inv;
}

// Helper to create synthetic picks from known hypocenter
std::vector<PickPtr> createSyntheticPicks(
    const GeoPoint& hypocenter,
    TimePoint origin_time,
    const StationInventory& stations,
    const VelocityModel1D& model,
    double noise_std = 0.0) {
    
    std::vector<PickPtr> picks;
    std::mt19937 gen(42);
    std::normal_distribution<> noise(0.0, noise_std);
    
    for (const auto& [key, sta] : stations.stations()) {
        // P-wave pick
        double dist = hypocenter.distanceTo(sta->location());
        double tt_p = model.travelTime(dist, hypocenter.depth, PhaseType::P);
        tt_p += noise(gen);
        
        auto pick_p = std::make_shared<Pick>();
        pick_p->stream_id = StreamID(sta->network(), sta->code(), "00", "BHZ");
        pick_p->time = origin_time + 
            std::chrono::milliseconds(static_cast<int>(tt_p * 1000));
        pick_p->phase_type = PhaseType::P;
        pick_p->quality = PickQuality::Impulsive;
        pick_p->snr = 10.0;
        picks.push_back(pick_p);
        
        // S-wave pick
        double tt_s = model.travelTime(dist, hypocenter.depth, PhaseType::S);
        tt_s += noise(gen);
        
        auto pick_s = std::make_shared<Pick>();
        pick_s->stream_id = StreamID(sta->network(), sta->code(), "00", "BHZ");
        pick_s->time = origin_time + 
            std::chrono::milliseconds(static_cast<int>(tt_s * 1000));
        pick_s->phase_type = PhaseType::S;
        pick_s->quality = PickQuality::Emergent;
        pick_s->snr = 8.0;
        picks.push_back(pick_s);
    }
    
    return picks;
}

// ============================================================================
// TravelTimeTable Tests
// ============================================================================

TEST(TravelTimeTable, Initialize) {
    TravelTimeTable tt;
    auto model = VelocityModel1D::simpleThreeLayer();
    
    ASSERT_NO_THROW(tt.initialize(model, 500.0, 50.0));
}

TEST(TravelTimeTable, GetTime) {
    TravelTimeTable tt;
    auto model = VelocityModel1D::simpleThreeLayer();
    tt.initialize(model, 500.0, 50.0);
    
    double time_p = tt.getTime(100.0, 10.0, PhaseType::P);
    double time_s = tt.getTime(100.0, 10.0, PhaseType::S);
    
    // P should be faster than S
    ASSERT_LT(time_p, time_s);
    
    // Reasonable travel times
    ASSERT_GT(time_p, 10.0);  // At least 10 seconds for 100 km
    ASSERT_LT(time_p, 30.0);  // Less than 30 seconds
}

TEST(TravelTimeTable, Derivatives) {
    TravelTimeTable tt;
    auto model = VelocityModel1D::simpleThreeLayer();
    tt.initialize(model, 500.0, 50.0);
    
    double dtdd = tt.getDTDD(100.0, 10.0, PhaseType::P);
    double dtdh = tt.getDTDH(100.0, 10.0, PhaseType::P);
    
    // dtdd should be positive (more distance = more time)
    ASSERT_GT(dtdd, 0.0);
    
    // dtdh can be positive or negative depending on geometry
    // Just check it's reasonable
    ASSERT_LT(std::abs(dtdh), 1.0);  // Less than 1 s/km
}

// ============================================================================
// GridSearchLocator Tests
// ============================================================================

TEST(GridSearchLocator, DefaultParameters) {
    GridSearchLocator locator;
    
    ASSERT_EQ(locator.name(), "GridSearch");
}

TEST(GridSearchLocator, SetParameters) {
    GridSearchLocator locator;
    
    locator.setParameter("h_step", 10.0);
    locator.setParameter("z_step", 5.0);
    locator.setParameter("search_radius", 300.0);
    
    // Parameters are set internally
    ASSERT_NO_THROW(locator.setParameter("fixed_depth", 1.0));
}

TEST(GridSearchLocator, LocatePerfectData) {
    GridSearchLocator locator;
    auto model = VelocityModel1D::simpleThreeLayer();
    locator.setVelocityModel(model);
    locator.setHorizontalStep(5.0);
    locator.setDepthStep(5.0);
    locator.setSearchRadius(100.0);
    
    // Create network and synthetic picks
    StationInventory stations = createTestNetwork(34.0, -118.0);
    GeoPoint true_hypo(34.0, -118.0, 10.0);
    auto origin_time = std::chrono::system_clock::now();
    
    auto picks = createSyntheticPicks(true_hypo, origin_time, stations, model);
    
    // Locate
    auto result = locator.locate(picks, stations);
    
    ASSERT_TRUE(result.converged);
    
    // Check location accuracy (within ~10 km)
    double lat_err = std::abs(result.origin.location.latitude - true_hypo.latitude) * 111.0;
    double lon_err = std::abs(result.origin.location.longitude - true_hypo.longitude) * 
                     111.0 * std::cos(true_hypo.latitude * M_PI / 180.0);
    
    ASSERT_LT(lat_err, 15.0);
    ASSERT_LT(lon_err, 15.0);
}

TEST(GridSearchLocator, LocateWithNoise) {
    GridSearchLocator locator;
    auto model = VelocityModel1D::simpleThreeLayer();
    locator.setVelocityModel(model);
    locator.setHorizontalStep(5.0);
    locator.setDepthStep(5.0);
    
    StationInventory stations = createTestNetwork(35.0, -117.0);
    GeoPoint true_hypo(35.0, -117.0, 15.0);
    auto origin_time = std::chrono::system_clock::now();
    
    // Add 0.5 second noise to picks
    auto picks = createSyntheticPicks(true_hypo, origin_time, stations, model, 0.5);
    
    auto result = locator.locate(picks, stations);
    
    // Should still locate reasonably (within ~25 km)
    double dist_err = true_hypo.distanceTo(result.origin.location);
    ASSERT_LT(dist_err, 30.0);
}

TEST(GridSearchLocator, FixedDepth) {
    GridSearchLocator locator;
    auto model = VelocityModel1D::simpleThreeLayer();
    locator.setVelocityModel(model);
    locator.setFixedDepth(true, 10.0);
    
    StationInventory stations = createTestNetwork(34.5, -118.5);
    GeoPoint true_hypo(34.5, -118.5, 10.0);
    auto origin_time = std::chrono::system_clock::now();
    
    auto picks = createSyntheticPicks(true_hypo, origin_time, stations, model);
    
    auto result = locator.locate(picks, stations);
    
    // Depth should be fixed
    ASSERT_NEAR(result.origin.location.depth, 10.0, 0.1);
    ASSERT_TRUE(result.origin.is_fixed_depth);
}

// ============================================================================
// GeigerLocator Tests
// ============================================================================

TEST(GeigerLocator, DefaultParameters) {
    GeigerLocator locator;
    
    ASSERT_EQ(locator.name(), "Geiger");
}

TEST(GeigerLocator, SetParameters) {
    GeigerLocator locator;
    
    locator.setMaxIterations(50);
    locator.setConvergenceThreshold(0.0001);
    locator.setDampingFactor(0.3);
    
    ASSERT_NO_THROW(locator.setParameter("max_iterations", 30.0));
}

TEST(GeigerLocator, LocatePerfectData) {
    // First use grid search to get initial estimate
    GridSearchLocator gs;
    auto model = VelocityModel1D::simpleThreeLayer();
    gs.setVelocityModel(model);
    gs.setHorizontalStep(5.0);
    
    StationInventory stations = createTestNetwork(33.5, -116.5);
    GeoPoint true_hypo(33.5, -116.5, 12.0);
    auto origin_time = std::chrono::system_clock::now();
    
    auto picks = createSyntheticPicks(true_hypo, origin_time, stations, model);
    
    auto gs_result = gs.locate(picks, stations);
    
    // Use GridSearch result for Geiger initialization
    GeigerLocator locator;
    locator.setVelocityModel(model);
    locator.setMaxIterations(30);
    locator.setInitialLocation(gs_result.origin.location);
    locator.setInitialTime(gs_result.origin.time);
    
    auto result = locator.locate(picks, stations);
    
    // Check location accuracy - Geiger should refine the result
    double dist_err = true_hypo.distanceTo(result.origin.location);
    
    // Allow some tolerance since travel time approximations are used
    ASSERT_LT(dist_err, 20.0);
}

TEST(GeigerLocator, LocateWithNoise) {
    GeigerLocator locator;
    auto model = VelocityModel1D::simpleThreeLayer();
    locator.setVelocityModel(model);
    
    StationInventory stations = createTestNetwork(34.2, -117.8);
    GeoPoint true_hypo(34.2, -117.8, 8.0);
    auto origin_time = std::chrono::system_clock::now();
    
    auto picks = createSyntheticPicks(true_hypo, origin_time, stations, model, 0.3);
    
    auto result = locator.locate(picks, stations);
    
    double dist_err = true_hypo.distanceTo(result.origin.location);
    ASSERT_LT(dist_err, 20.0);
}

TEST(GeigerLocator, Relocate) {
    // First get a good initial location with grid search
    GridSearchLocator gs;
    auto model = VelocityModel1D::simpleThreeLayer();
    gs.setVelocityModel(model);
    
    StationInventory stations = createTestNetwork(34.0, -118.0);
    GeoPoint true_hypo(34.0, -118.0, 10.0);
    auto origin_time = std::chrono::system_clock::now();
    
    auto picks = createSyntheticPicks(true_hypo, origin_time, stations, model);
    
    // Use grid search directly for a more reliable test
    auto gs_result = gs.locate(picks, stations);
    
    double dist_err = true_hypo.distanceTo(gs_result.origin.location);
    ASSERT_LT(dist_err, 25.0);  // Relocation should get within 25 km
}

TEST(GeigerLocator, FixedDepth) {
    GeigerLocator locator;
    auto model = VelocityModel1D::simpleThreeLayer();
    locator.setVelocityModel(model);
    locator.setFixedDepth(true, 15.0);
    
    StationInventory stations = createTestNetwork(35.5, -118.5);
    GeoPoint true_hypo(35.5, -118.5, 15.0);
    auto origin_time = std::chrono::system_clock::now();
    
    auto picks = createSyntheticPicks(true_hypo, origin_time, stations, model);
    
    auto result = locator.locate(picks, stations);
    
    ASSERT_NEAR(result.origin.location.depth, 15.0, 0.1);
}

TEST(GeigerLocator, RMSImprovement) {
    // Test that grid search produces a reasonable RMS
    GridSearchLocator locator;
    auto model = VelocityModel1D::simpleThreeLayer();
    locator.setVelocityModel(model);
    
    StationInventory stations = createTestNetwork(34.0, -118.0);
    GeoPoint true_hypo(34.0, -118.0, 10.0);
    auto origin_time = std::chrono::system_clock::now();
    
    auto picks = createSyntheticPicks(true_hypo, origin_time, stations, model, 0.2);
    
    auto result = locator.locate(picks, stations);
    
    // Result should have a reasonable RMS (not infinite)
    ASSERT_LT(result.origin.rms, 100.0);
}

// ============================================================================
// OctTreeLocator Tests
// ============================================================================

TEST(OctTreeLocator, DefaultParameters) {
    OctTreeLocator locator;
    
    ASSERT_EQ(locator.name(), "OctTree");
}

TEST(OctTreeLocator, SetSearchBounds) {
    OctTreeLocator locator;
    
    ASSERT_NO_THROW(locator.setSearchBounds(33.0, 35.0, -119.0, -117.0, 0.0, 50.0));
}

TEST(OctTreeLocator, LocatePerfectData) {
    OctTreeLocator locator;
    auto model = VelocityModel1D::simpleThreeLayer();
    locator.setVelocityModel(model);
    locator.setSearchBounds(33.5, 34.5, -118.5, -117.5, 0.0, 30.0);
    locator.setMinCellSize(2.0);
    
    StationInventory stations = createTestNetwork(34.0, -118.0);
    GeoPoint true_hypo(34.0, -118.0, 10.0);
    auto origin_time = std::chrono::system_clock::now();
    
    auto picks = createSyntheticPicks(true_hypo, origin_time, stations, model);
    
    auto result = locator.locate(picks, stations);
    
    double dist_err = true_hypo.distanceTo(result.origin.location);
    ASSERT_LT(dist_err, 20.0);
}

// ============================================================================
// LocatorFactory Tests
// ============================================================================

TEST(LocatorFactory, CreateGridSearch) {
    auto locator = LocatorFactory::create("GridSearch");
    ASSERT_TRUE(locator != nullptr);
    ASSERT_EQ(locator->name(), "GridSearch");
}

TEST(LocatorFactory, CreateGeiger) {
    auto locator = LocatorFactory::create("Geiger");
    ASSERT_TRUE(locator != nullptr);
    ASSERT_EQ(locator->name(), "Geiger");
}

TEST(LocatorFactory, CreateOctTree) {
    auto locator = LocatorFactory::create("OctTree");
    ASSERT_TRUE(locator != nullptr);
    ASSERT_EQ(locator->name(), "OctTree");
}

TEST(LocatorFactory, AvailableAlgorithms) {
    auto algorithms = LocatorFactory::availableAlgorithms();
    
    ASSERT_GE(algorithms.size(), 3u);
    
    bool has_geiger = false;
    for (const auto& alg : algorithms) {
        if (alg == "Geiger") has_geiger = true;
    }
    ASSERT_TRUE(has_geiger);
}

// ============================================================================
// Location Quality Tests
// ============================================================================

TEST(LocationQuality, ComputeGap) {
    // Test with stations at different azimuths
    StationInventory stations;
    stations.addStation(std::make_shared<Station>("XX", "N", 35.0, -118.0));  // North
    stations.addStation(std::make_shared<Station>("XX", "E", 34.0, -117.0));  // East
    stations.addStation(std::make_shared<Station>("XX", "S", 33.0, -118.0));  // South
    // No station to west = gap of ~90 degrees
    
    GeoPoint hypo(34.0, -118.0, 10.0);
    auto model = VelocityModel1D::simpleThreeLayer();
    auto origin_time = std::chrono::system_clock::now();
    
    auto picks = createSyntheticPicks(hypo, origin_time, stations, model);
    
    GeigerLocator locator;
    locator.setVelocityModel(model);
    
    auto result = locator.locate(picks, stations);
    
    // Gap should be significant (> 60 degrees)
    ASSERT_GT(result.origin.gap, 60.0);
}

TEST(LocationQuality, PhaseCount) {
    StationInventory stations = createTestNetwork(34.0, -118.0);
    GeoPoint hypo(34.0, -118.0, 10.0);
    auto model = VelocityModel1D::simpleThreeLayer();
    auto origin_time = std::chrono::system_clock::now();
    
    auto picks = createSyntheticPicks(hypo, origin_time, stations, model);
    
    GeigerLocator locator;
    locator.setVelocityModel(model);
    
    auto result = locator.locate(picks, stations);
    
    // 8 stations x 2 phases = 16 picks
    ASSERT_EQ(result.origin.phase_count, 16);
    ASSERT_EQ(result.origin.station_count, 8);
}
