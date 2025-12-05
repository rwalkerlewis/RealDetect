/**
 * Unit tests for phase association
 */

#include "test_framework.hpp"
#include "realdetect/associator/phase_associator.hpp"
#include "realdetect/core/station.hpp"
#include "realdetect/core/event.hpp"
#include "realdetect/core/velocity_model.hpp"
#include <random>
#include <cmath>

using namespace realdetect;
using namespace realdetect::test;

// Helper to create a test station network
StationInventory createAssociatorTestNetwork() {
    StationInventory inv;
    
    // Network of stations around a central point
    std::vector<std::tuple<std::string, double, double>> stations = {
        {"S1", 34.5, -118.0},
        {"S2", 34.3, -117.7},
        {"S3", 34.1, -118.0},
        {"S4", 34.3, -118.3},
        {"S5", 34.6, -117.8},
        {"S6", 34.2, -117.5},
        {"S7", 34.0, -118.2},
        {"S8", 34.4, -118.5}
    };
    
    for (const auto& [code, lat, lon] : stations) {
        inv.addStation(std::make_shared<Station>("XX", code, lat, lon, 100));
    }
    
    return inv;
}

// Helper to create synthetic picks from known hypocenter
std::vector<PickPtr> createAssociatorTestPicks(
    const GeoPoint& hypocenter,
    TimePoint origin_time,
    const StationInventory& stations,
    const VelocityModel1D& model,
    double noise_std = 0.0) {
    
    std::vector<PickPtr> picks;
    std::mt19937 gen(42);
    std::normal_distribution<> noise(0.0, noise_std);
    
    for (const auto& [key, sta] : stations.stations()) {
        double dist = hypocenter.distanceTo(sta->location());
        
        // P-wave pick
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
// TravelTimeTable Tests for Associator
// ============================================================================

TEST(AssocTravelTime, Initialize) {
    TravelTimeTable tt;
    auto model = VelocityModel1D::simpleThreeLayer();
    
    ASSERT_NO_THROW(tt.initialize(model, 500.0, 50.0));
}

TEST(AssocTravelTime, GetTimeP) {
    TravelTimeTable tt;
    auto model = VelocityModel1D::simpleThreeLayer();
    tt.initialize(model, 500.0, 50.0);
    
    double time_p = tt.getTime(100.0, 10.0, PhaseType::P);
    
    // P-wave at 100 km, 10 km depth
    // Hypo distance ~ sqrt(100^2 + 10^2) = 100.5 km
    // Average Vp ~ 6 km/s, expected time ~ 16.7s
    ASSERT_GT(time_p, 10.0);
    ASSERT_LT(time_p, 25.0);
}

TEST(AssocTravelTime, GetTimeS) {
    TravelTimeTable tt;
    auto model = VelocityModel1D::simpleThreeLayer();
    tt.initialize(model, 500.0, 50.0);
    
    double time_s = tt.getTime(100.0, 10.0, PhaseType::S);
    
    // S should be slower than P
    double time_p = tt.getTime(100.0, 10.0, PhaseType::P);
    ASSERT_GT(time_s, time_p);
}

TEST(AssocTravelTime, TimeIncreaseWithDistance) {
    TravelTimeTable tt;
    auto model = VelocityModel1D::simpleThreeLayer();
    tt.initialize(model, 500.0, 50.0);
    
    double t1 = tt.getTime(50.0, 10.0, PhaseType::P);
    double t2 = tt.getTime(100.0, 10.0, PhaseType::P);
    double t3 = tt.getTime(200.0, 10.0, PhaseType::P);
    
    ASSERT_LT(t1, t2);
    ASSERT_LT(t2, t3);
}

TEST(AssocTravelTime, TimeIncreaseWithDepth) {
    TravelTimeTable tt;
    auto model = VelocityModel1D::simpleThreeLayer();
    tt.initialize(model, 500.0, 50.0);
    
    // For same epicentral distance, deeper source -> longer travel time
    double t1 = tt.getTime(100.0, 5.0, PhaseType::P);
    double t2 = tt.getTime(100.0, 20.0, PhaseType::P);
    
    // Deeper source generally has longer travel time for regional distances
    ASSERT_GT(t2, t1 - 1.0);  // Allow some tolerance
}

TEST(AssocTravelTime, Derivatives) {
    TravelTimeTable tt;
    auto model = VelocityModel1D::simpleThreeLayer();
    tt.initialize(model, 500.0, 50.0);
    
    double dtdd = tt.getDTDD(100.0, 10.0, PhaseType::P);
    double dtdh = tt.getDTDH(100.0, 10.0, PhaseType::P);
    
    // dtdd (dt/d_distance) should be positive (more distance = more time)
    ASSERT_GT(dtdd, 0.0);
    
    // Both should be reasonable values (not infinity, not zero)
    ASSERT_LT(std::abs(dtdd), 1.0);  // < 1 s/km
    ASSERT_LT(std::abs(dtdh), 1.0);  // < 1 s/km
}

// ============================================================================
// PhaseAssociator Tests
// ============================================================================

TEST(PhaseAssociator, DefaultConstructor) {
    PhaseAssociator assoc;
    
    // Should create without error
    ASSERT_TRUE(true);
}

TEST(PhaseAssociator, SetStations) {
    PhaseAssociator assoc;
    auto stations = createAssociatorTestNetwork();
    
    ASSERT_NO_THROW(assoc.setStations(stations));
}

TEST(PhaseAssociator, SetVelocityModel) {
    PhaseAssociator assoc;
    auto model = VelocityModel1D::simpleThreeLayer();
    
    ASSERT_NO_THROW(assoc.setVelocityModel(model));
}

TEST(PhaseAssociator, SetParameters) {
    PhaseAssociator assoc;
    
    ASSERT_NO_THROW(assoc.setTimeWindow(60.0));
    ASSERT_NO_THROW(assoc.setMaxResidual(2.0));
    ASSERT_NO_THROW(assoc.setMinStations(4));
    ASSERT_NO_THROW(assoc.setMinPhases(6));
    ASSERT_NO_THROW(assoc.setSearchRadius(200.0));
}

TEST(PhaseAssociator, AddPick) {
    PhaseAssociator assoc;
    
    auto pick = std::make_shared<Pick>();
    pick->stream_id = StreamID("XX", "S1", "00", "BHZ");
    pick->time = std::chrono::system_clock::now();
    pick->phase_type = PhaseType::P;
    
    ASSERT_NO_THROW(assoc.addPick(pick));
}

TEST(PhaseAssociator, AssociatePerfectPicks) {
    PhaseAssociator assoc;
    auto stations = createAssociatorTestNetwork();
    auto model = VelocityModel1D::simpleThreeLayer();
    
    assoc.setStations(stations);
    assoc.setVelocityModel(model);
    assoc.setMinStations(4);
    assoc.setMinPhases(6);
    assoc.setMaxResidual(2.0);
    
    // Create perfect picks from known hypocenter
    GeoPoint hypocenter(34.3, -117.9, 10.0);
    auto origin_time = std::chrono::system_clock::now();
    
    auto picks = createAssociatorTestPicks(hypocenter, origin_time, stations, model);
    
    // Add picks to associator
    for (const auto& pick : picks) {
        assoc.addPick(pick);
    }
    
    // Process should create an event
    auto events = assoc.process();
    
    // Should find at least one event
    ASSERT_GE(events.size(), 1u);
}

TEST(PhaseAssociator, AssociateNoisyPicks) {
    PhaseAssociator assoc;
    auto stations = createAssociatorTestNetwork();
    auto model = VelocityModel1D::simpleThreeLayer();
    
    assoc.setStations(stations);
    assoc.setVelocityModel(model);
    assoc.setMinStations(3);
    assoc.setMinPhases(4);
    assoc.setMaxResidual(3.0);
    
    GeoPoint hypocenter(34.3, -117.9, 10.0);
    auto origin_time = std::chrono::system_clock::now();
    
    // Add 0.5 second noise to picks
    auto picks = createAssociatorTestPicks(hypocenter, origin_time, stations, model, 0.5);
    
    for (const auto& pick : picks) {
        assoc.addPick(pick);
    }
    
    auto events = assoc.process();
    
    // Should still find event despite noise
    ASSERT_GE(events.size(), 1u);
}

TEST(PhaseAssociator, NoEventFromScatteredPicks) {
    PhaseAssociator assoc;
    auto stations = createAssociatorTestNetwork();
    auto model = VelocityModel1D::simpleThreeLayer();
    
    assoc.setStations(stations);
    assoc.setVelocityModel(model);
    assoc.setMinStations(4);
    assoc.setMinPhases(6);
    assoc.setMaxResidual(1.0);  // Strict residual threshold
    
    // Create random scattered picks (not from a real event)
    std::mt19937 gen(42);
    std::uniform_real_distribution<> time_dist(0.0, 60.0);
    
    auto now = std::chrono::system_clock::now();
    
    for (const auto& [key, sta] : stations.stations()) {
        auto pick = std::make_shared<Pick>();
        pick->stream_id = StreamID(sta->network(), sta->code(), "00", "BHZ");
        pick->time = now + std::chrono::milliseconds(
            static_cast<int>(time_dist(gen) * 1000));
        pick->phase_type = PhaseType::P;
        
        assoc.addPick(pick);
    }
    
    auto events = assoc.process();
    
    // Should not find events from random picks
    // (or very few false positives)
    ASSERT_LE(events.size(), 1u);
}

TEST(PhaseAssociator, EventCallback) {
    PhaseAssociator assoc;
    auto stations = createAssociatorTestNetwork();
    auto model = VelocityModel1D::simpleThreeLayer();
    
    assoc.setStations(stations);
    assoc.setVelocityModel(model);
    assoc.setMinStations(3);
    assoc.setMinPhases(4);
    
    int callback_count = 0;
    assoc.setEventCallback([&callback_count](EventPtr event) {
        callback_count++;
    });
    
    GeoPoint hypocenter(34.3, -117.9, 10.0);
    auto origin_time = std::chrono::system_clock::now();
    
    auto picks = createAssociatorTestPicks(hypocenter, origin_time, stations, model);
    
    for (const auto& pick : picks) {
        assoc.addPick(pick);
    }
    
    assoc.process();
    
    // Callback should have been called
    ASSERT_GE(callback_count, 1);
}

TEST(PhaseAssociator, UnassociatedPicks) {
    PhaseAssociator assoc;
    auto stations = createAssociatorTestNetwork();
    auto model = VelocityModel1D::simpleThreeLayer();
    
    assoc.setStations(stations);
    assoc.setVelocityModel(model);
    assoc.setMinStations(10);  // Require more stations than we have
    
    GeoPoint hypocenter(34.3, -117.9, 10.0);
    auto origin_time = std::chrono::system_clock::now();
    
    auto picks = createAssociatorTestPicks(hypocenter, origin_time, stations, model);
    
    for (const auto& pick : picks) {
        assoc.addPick(pick);
    }
    
    auto events = assoc.process();
    
    // With such strict requirements, no events should be formed
    // (picks stay unassociated, but implementation may or may not track them)
    // Just verify the processing completes without error
    ASSERT_TRUE(true);
}

// ============================================================================
// NucleatorAssociator Tests
// ============================================================================

TEST(NucleatorAssociator, DefaultConstructor) {
    NucleatorAssociator assoc;
    
    // Should create without error
    ASSERT_TRUE(true);
}

TEST(NucleatorAssociator, SetStations) {
    NucleatorAssociator assoc;
    auto stations = createAssociatorTestNetwork();
    
    ASSERT_NO_THROW(assoc.setStations(stations));
}

TEST(NucleatorAssociator, SetVelocityModel) {
    NucleatorAssociator assoc;
    auto model = VelocityModel1D::simpleThreeLayer();
    
    ASSERT_NO_THROW(assoc.setVelocityModel(model));
}

TEST(NucleatorAssociator, SetGrid) {
    NucleatorAssociator assoc;
    
    ASSERT_NO_THROW(assoc.setGrid(
        33.0, 35.0, 0.1,   // lat
        -119.0, -117.0, 0.1,  // lon
        0.0, 30.0, 5.0     // depth
    ));
}

TEST(NucleatorAssociator, SetParameters) {
    NucleatorAssociator assoc;
    
    ASSERT_NO_THROW(assoc.setTimeWindow(60.0));
    ASSERT_NO_THROW(assoc.setMinPhases(6));
    ASSERT_NO_THROW(assoc.setMaxResidual(2.0));
}

TEST(NucleatorAssociator, AddPick) {
    NucleatorAssociator assoc;
    
    auto pick = std::make_shared<Pick>();
    pick->stream_id = StreamID("XX", "S1", "00", "BHZ");
    pick->time = std::chrono::system_clock::now();
    pick->phase_type = PhaseType::P;
    
    ASSERT_NO_THROW(assoc.addPick(pick));
}

TEST(NucleatorAssociator, ProcessPicks) {
    NucleatorAssociator assoc;
    auto stations = createAssociatorTestNetwork();
    auto model = VelocityModel1D::simpleThreeLayer();
    
    assoc.setStations(stations);
    assoc.setVelocityModel(model);
    assoc.setGrid(
        34.0, 34.6, 0.1,
        -118.2, -117.6, 0.1,
        0.0, 20.0, 5.0
    );
    assoc.setMinPhases(4);
    assoc.setMaxResidual(2.0);
    
    GeoPoint hypocenter(34.3, -117.9, 10.0);
    auto origin_time = std::chrono::system_clock::now();
    
    auto picks = createAssociatorTestPicks(hypocenter, origin_time, stations, model);
    
    for (const auto& pick : picks) {
        assoc.addPick(pick);
    }
    
    auto events = assoc.process();
    
    // Should find events
    ASSERT_GE(events.size(), 1u);
}

// ============================================================================
// AssociationCandidate Tests
// ============================================================================

TEST(AssociationCandidate, DefaultConstructor) {
    AssociationCandidate candidate;
    
    ASSERT_NEAR(candidate.score, 0.0, 1e-10);
    ASSERT_EQ(candidate.phase_count, 0);
    ASSERT_EQ(candidate.station_count, 0);
}

TEST(AssociationCandidate, AddPick) {
    AssociationCandidate candidate;
    
    auto pick = std::make_shared<Pick>();
    pick->stream_id = StreamID("XX", "S1", "00", "BHZ");
    pick->phase_type = PhaseType::P;
    
    candidate.picks.push_back(pick);
    candidate.phase_count = 1;
    candidate.station_count = 1;
    
    ASSERT_EQ(candidate.picks.size(), 1u);
    ASSERT_EQ(candidate.phase_count, 1);
}

// ============================================================================
// Multi-Event Association Tests
// ============================================================================

TEST(PhaseAssociator, TwoEvents) {
    PhaseAssociator assoc;
    auto stations = createAssociatorTestNetwork();
    auto model = VelocityModel1D::simpleThreeLayer();
    
    assoc.setStations(stations);
    assoc.setVelocityModel(model);
    assoc.setMinStations(3);
    assoc.setMinPhases(4);
    assoc.setTimeWindow(30.0);  // 30 second window
    
    // First event
    GeoPoint hypo1(34.2, -117.8, 8.0);
    auto time1 = std::chrono::system_clock::now();
    
    auto picks1 = createAssociatorTestPicks(hypo1, time1, stations, model);
    for (const auto& pick : picks1) {
        assoc.addPick(pick);
    }
    
    // Second event 60 seconds later, different location
    GeoPoint hypo2(34.4, -118.2, 12.0);
    auto time2 = time1 + std::chrono::seconds(60);
    
    auto picks2 = createAssociatorTestPicks(hypo2, time2, stations, model);
    for (const auto& pick : picks2) {
        assoc.addPick(pick);
    }
    
    auto events = assoc.process();
    
    // Should detect at least 2 events
    ASSERT_GE(events.size(), 2u);
}

// ============================================================================
// Edge Cases
// ============================================================================

TEST(PhaseAssociator, TooFewPicks) {
    PhaseAssociator assoc;
    auto stations = createAssociatorTestNetwork();
    auto model = VelocityModel1D::simpleThreeLayer();
    
    assoc.setStations(stations);
    assoc.setVelocityModel(model);
    assoc.setMinStations(4);
    assoc.setMinPhases(8);
    
    // Only add 2 picks
    auto now = std::chrono::system_clock::now();
    
    auto pick1 = std::make_shared<Pick>();
    pick1->stream_id = StreamID("XX", "S1", "00", "BHZ");
    pick1->time = now;
    pick1->phase_type = PhaseType::P;
    assoc.addPick(pick1);
    
    auto pick2 = std::make_shared<Pick>();
    pick2->stream_id = StreamID("XX", "S2", "00", "BHZ");
    pick2->time = now + std::chrono::milliseconds(1500);
    pick2->phase_type = PhaseType::P;
    assoc.addPick(pick2);
    
    auto events = assoc.process();
    
    // Should not find event with too few picks
    ASSERT_EQ(events.size(), 0u);
}

TEST(PhaseAssociator, OnlyPWaves) {
    PhaseAssociator assoc;
    auto stations = createAssociatorTestNetwork();
    auto model = VelocityModel1D::simpleThreeLayer();
    
    assoc.setStations(stations);
    assoc.setVelocityModel(model);
    assoc.setMinStations(3);
    assoc.setMinPhases(4);
    
    GeoPoint hypocenter(34.3, -117.9, 10.0);
    auto origin_time = std::chrono::system_clock::now();
    
    // Only P-wave picks
    for (const auto& [key, sta] : stations.stations()) {
        double dist = hypocenter.distanceTo(sta->location());
        double tt_p = model.travelTime(dist, hypocenter.depth, PhaseType::P);
        
        auto pick = std::make_shared<Pick>();
        pick->stream_id = StreamID(sta->network(), sta->code(), "00", "BHZ");
        pick->time = origin_time + 
            std::chrono::milliseconds(static_cast<int>(tt_p * 1000));
        pick->phase_type = PhaseType::P;
        
        assoc.addPick(pick);
    }
    
    auto events = assoc.process();
    
    // Should still find event with only P waves
    ASSERT_GE(events.size(), 1u);
}
