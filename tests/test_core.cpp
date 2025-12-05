/**
 * Unit tests for core components
 */

#include "test_framework.hpp"
#include "realdetect/core/types.hpp"
#include "realdetect/core/waveform.hpp"
#include "realdetect/core/station.hpp"
#include "realdetect/core/event.hpp"
#include "realdetect/core/velocity_model.hpp"
#include "realdetect/core/config.hpp"

using namespace realdetect;
using namespace realdetect::test;

// ============================================================================
// GeoPoint Tests
// ============================================================================

TEST(GeoPoint, Constructor) {
    GeoPoint p(34.5, -118.0, 10.0);
    ASSERT_NEAR(p.latitude, 34.5, 1e-10);
    ASSERT_NEAR(p.longitude, -118.0, 1e-10);
    ASSERT_NEAR(p.depth, 10.0, 1e-10);
}

TEST(GeoPoint, DefaultConstructor) {
    GeoPoint p;
    ASSERT_NEAR(p.latitude, 0.0, 1e-10);
    ASSERT_NEAR(p.longitude, 0.0, 1e-10);
    ASSERT_NEAR(p.depth, 0.0, 1e-10);
}

TEST(GeoPoint, DistanceToSamePoint) {
    GeoPoint p1(34.5, -118.0);
    GeoPoint p2(34.5, -118.0);
    ASSERT_NEAR(p1.distanceTo(p2), 0.0, 1e-6);
}

TEST(GeoPoint, DistanceToNearby) {
    // Los Angeles to San Diego ~180 km
    GeoPoint la(34.05, -118.25);
    GeoPoint sd(32.72, -117.16);
    double dist = la.distanceTo(sd);
    ASSERT_NEAR(dist, 180.0, 20.0);  // Within 20 km
}

TEST(GeoPoint, DistanceLondon_NewYork) {
    // London to New York ~5570 km
    GeoPoint london(51.5, -0.12);
    GeoPoint nyc(40.71, -74.01);
    double dist = london.distanceTo(nyc);
    ASSERT_NEAR(dist, 5570.0, 50.0);  // Within 50 km
}

TEST(GeoPoint, AzimuthNorth) {
    GeoPoint p1(34.0, -118.0);
    GeoPoint p2(35.0, -118.0);  // Due north
    double az = p1.azimuthTo(p2);
    ASSERT_NEAR(az, 0.0, 1.0);  // Within 1 degree
}

TEST(GeoPoint, AzimuthEast) {
    GeoPoint p1(34.0, -118.0);
    GeoPoint p2(34.0, -117.0);  // Due east
    double az = p1.azimuthTo(p2);
    ASSERT_NEAR(az, 90.0, 1.0);
}

TEST(GeoPoint, AzimuthSouth) {
    GeoPoint p1(34.0, -118.0);
    GeoPoint p2(33.0, -118.0);  // Due south
    double az = p1.azimuthTo(p2);
    ASSERT_NEAR(az, 180.0, 1.0);
}

// ============================================================================
// StreamID Tests
// ============================================================================

TEST(StreamID, Constructor) {
    StreamID id("IU", "ANMO", "00", "BHZ");
    ASSERT_EQ(id.network, "IU");
    ASSERT_EQ(id.station, "ANMO");
    ASSERT_EQ(id.location, "00");
    ASSERT_EQ(id.channel, "BHZ");
}

TEST(StreamID, ToString) {
    StreamID id("IU", "ANMO", "00", "BHZ");
    ASSERT_EQ(id.toString(), "IU.ANMO.00.BHZ");
}

TEST(StreamID, Equality) {
    StreamID id1("IU", "ANMO", "00", "BHZ");
    StreamID id2("IU", "ANMO", "00", "BHZ");
    StreamID id3("IU", "ANMO", "00", "BHN");
    ASSERT_TRUE(id1 == id2);
    ASSERT_FALSE(id1 == id3);
}

TEST(StreamID, Comparison) {
    StreamID id1("AA", "STA1", "00", "BHZ");
    StreamID id2("AB", "STA1", "00", "BHZ");
    ASSERT_TRUE(id1 < id2);
}

// ============================================================================
// Waveform Tests
// ============================================================================

TEST(Waveform, Constructor) {
    StreamID id("IU", "ANMO", "00", "BHZ");
    auto now = std::chrono::system_clock::now();
    Waveform wf(id, 100.0, now);
    
    ASSERT_EQ(wf.streamId().toString(), "IU.ANMO.00.BHZ");
    ASSERT_NEAR(wf.sampleRate(), 100.0, 1e-10);
    ASSERT_EQ(wf.sampleCount(), 0u);
}

TEST(Waveform, AppendSamples) {
    StreamID id("IU", "ANMO", "00", "BHZ");
    auto now = std::chrono::system_clock::now();
    Waveform wf(id, 100.0, now);
    
    wf.append(1.0);
    wf.append(2.0);
    wf.append(3.0);
    
    ASSERT_EQ(wf.sampleCount(), 3u);
    ASSERT_NEAR(wf[0], 1.0, 1e-10);
    ASSERT_NEAR(wf[1], 2.0, 1e-10);
    ASSERT_NEAR(wf[2], 3.0, 1e-10);
}

TEST(Waveform, Mean) {
    StreamID id("IU", "ANMO", "00", "BHZ");
    auto now = std::chrono::system_clock::now();
    Waveform wf(id, 100.0, now);
    
    for (int i = 0; i < 100; i++) {
        wf.append(i);
    }
    
    // Mean of 0..99 is 49.5
    ASSERT_NEAR(wf.mean(), 49.5, 1e-10);
}

TEST(Waveform, Demean) {
    StreamID id("IU", "ANMO", "00", "BHZ");
    auto now = std::chrono::system_clock::now();
    Waveform wf(id, 100.0, now);
    
    wf.append(10.0);
    wf.append(20.0);
    wf.append(30.0);
    
    wf.demean();
    
    ASSERT_NEAR(wf.mean(), 0.0, 1e-10);
    ASSERT_NEAR(wf[0], -10.0, 1e-10);
    ASSERT_NEAR(wf[1], 0.0, 1e-10);
    ASSERT_NEAR(wf[2], 10.0, 1e-10);
}

TEST(Waveform, AbsMax) {
    StreamID id("IU", "ANMO", "00", "BHZ");
    auto now = std::chrono::system_clock::now();
    Waveform wf(id, 100.0, now);
    
    wf.append(-5.0);
    wf.append(10.0);
    wf.append(-15.0);
    wf.append(3.0);
    
    ASSERT_NEAR(wf.absMax(), 15.0, 1e-10);
}

TEST(Waveform, Normalize) {
    StreamID id("IU", "ANMO", "00", "BHZ");
    auto now = std::chrono::system_clock::now();
    Waveform wf(id, 100.0, now);
    
    wf.append(-5.0);
    wf.append(10.0);
    wf.append(-15.0);
    
    wf.normalize();
    
    ASSERT_NEAR(wf.absMax(), 1.0, 1e-10);
    ASSERT_NEAR(wf[2], -1.0, 1e-10);  // Was max magnitude
}

TEST(Waveform, Slice) {
    StreamID id("IU", "ANMO", "00", "BHZ");
    auto now = std::chrono::system_clock::now();
    Waveform wf(id, 100.0, now);
    
    for (int i = 0; i < 100; i++) {
        wf.append(static_cast<double>(i));
    }
    
    Waveform slice = wf.slice(10, 20);
    
    ASSERT_EQ(slice.sampleCount(), 10u);
    ASSERT_NEAR(slice[0], 10.0, 1e-10);
    ASSERT_NEAR(slice[9], 19.0, 1e-10);
}

TEST(Waveform, Duration) {
    StreamID id("IU", "ANMO", "00", "BHZ");
    auto now = std::chrono::system_clock::now();
    Waveform wf(id, 100.0, now);
    
    for (int i = 0; i < 1000; i++) {
        wf.append(0.0);
    }
    
    // 1000 samples at 100 Hz = 10 seconds
    ASSERT_NEAR(wf.duration(), 10.0, 1e-10);
}

// ============================================================================
// Station Tests
// ============================================================================

TEST(Station, Constructor) {
    Station sta("IU", "ANMO", 34.946, -106.457, 1850);
    
    ASSERT_EQ(sta.network(), "IU");
    ASSERT_EQ(sta.code(), "ANMO");
    ASSERT_NEAR(sta.latitude(), 34.946, 1e-10);
    ASSERT_NEAR(sta.longitude(), -106.457, 1e-10);
    ASSERT_NEAR(sta.elevation(), 1850, 1e-10);
}

TEST(Station, AddChannel) {
    Station sta("IU", "ANMO", 34.946, -106.457, 1850);
    
    Channel bhz;
    bhz.code = "BHZ";
    bhz.sample_rate = 40.0;
    bhz.dip = -90.0;
    sta.addChannel(bhz);
    
    const Channel* ch = sta.getChannel("BHZ");
    ASSERT_TRUE(ch != nullptr);
    ASSERT_EQ(ch->code, "BHZ");
    ASSERT_NEAR(ch->sample_rate, 40.0, 1e-10);
    ASSERT_TRUE(ch->isVertical());
}

TEST(Station, DistanceTo) {
    Station sta("CI", "PAS", 34.15, -118.17, 300);
    GeoPoint pt(35.15, -118.17);  // 1 degree north
    
    double dist = sta.distanceTo(pt);
    ASSERT_NEAR(dist, 111.0, 5.0);  // ~111 km per degree
}

TEST(StationInventory, AddAndGet) {
    StationInventory inv;
    
    auto sta1 = std::make_shared<Station>("IU", "ANMO", 34.946, -106.457);
    auto sta2 = std::make_shared<Station>("IU", "CCM", 38.056, -91.244);
    
    inv.addStation(sta1);
    inv.addStation(sta2);
    
    ASSERT_EQ(inv.size(), 2u);
    
    auto found = inv.getStation("IU", "ANMO");
    ASSERT_TRUE(found != nullptr);
    ASSERT_EQ(found->code(), "ANMO");
    
    auto not_found = inv.getStation("XX", "XXX");
    ASSERT_TRUE(not_found == nullptr);
}

TEST(StationInventory, StationsWithin) {
    StationInventory inv;
    
    // Stations around Los Angeles
    inv.addStation(std::make_shared<Station>("CI", "PAS", 34.15, -118.17));
    inv.addStation(std::make_shared<Station>("CI", "USC", 34.02, -118.29));
    inv.addStation(std::make_shared<Station>("CI", "FAR", 33.0, -117.0));  // Far away
    
    GeoPoint center(34.0, -118.0);
    auto nearby = inv.stationsWithin(center, 50.0);  // 50 km radius
    
    ASSERT_EQ(nearby.size(), 2u);  // PAS and USC, not FAR
}

// ============================================================================
// VelocityModel Tests
// ============================================================================

TEST(VelocityModel, SimpleThreeLayer) {
    auto model = VelocityModel1D::simpleThreeLayer();
    
    ASSERT_EQ(model.layerCount(), 3u);
    
    // Check surface velocity
    ASSERT_NEAR(model.vpAt(0.0), 5.5, 0.1);
    ASSERT_NEAR(model.vsAt(0.0), 3.18, 0.1);
    
    // Check mantle velocity
    ASSERT_NEAR(model.vpAt(40.0), 8.0, 0.1);
}

TEST(VelocityModel, IASP91) {
    auto model = VelocityModel1D::iasp91();
    
    ASSERT_GE(model.layerCount(), 4u);
    
    // Check Moho discontinuity
    ASSERT_GT(model.vpAt(36.0), model.vpAt(34.0));
}

TEST(VelocityModel, TravelTime) {
    auto model = VelocityModel1D::simpleThreeLayer();
    
    // P-wave at 100 km, 10 km depth
    double tt_p = model.travelTime(100.0, 10.0, PhaseType::P);
    
    // Should be roughly distance / velocity
    // Hypo distance = sqrt(100^2 + 10^2) ~ 100.5 km
    // Average Vp ~ 6 km/s
    // Expected tt ~ 16.7 s
    ASSERT_NEAR(tt_p, 16.7, 3.0);
    
    // S-wave should be slower
    double tt_s = model.travelTime(100.0, 10.0, PhaseType::S);
    ASSERT_GT(tt_s, tt_p);
}

TEST(VelocityModel, AddLayer) {
    VelocityModel1D model("Custom");
    
    model.addLayer(0.0, 10.0, 5.0, 2.9);
    model.addLayer(10.0, 0.0, 8.0, 4.5);
    
    ASSERT_EQ(model.layerCount(), 2u);
    ASSERT_NEAR(model.vpAt(5.0), 5.0, 1e-10);
    ASSERT_NEAR(model.vpAt(15.0), 8.0, 1e-10);
}

// ============================================================================
// Config Tests
// ============================================================================

TEST(Config, SetAndGet) {
    Config cfg;
    
    cfg.set("key1", "value1");
    cfg.set("key2", 42);
    cfg.set("key3", 3.14);
    cfg.set("key4", true);
    
    ASSERT_EQ(cfg.getString("key1"), "value1");
    ASSERT_EQ(cfg.getInt("key2"), 42);
    ASSERT_NEAR(cfg.getDouble("key3"), 3.14, 1e-10);
    ASSERT_TRUE(cfg.getBool("key4"));
}

TEST(Config, DefaultValues) {
    Config cfg;
    
    ASSERT_EQ(cfg.getString("missing", "default"), "default");
    ASSERT_EQ(cfg.getInt("missing", 99), 99);
    ASSERT_NEAR(cfg.getDouble("missing", 1.5), 1.5, 1e-10);
    ASSERT_FALSE(cfg.getBool("missing", false));
}

TEST(Config, HasKey) {
    Config cfg;
    
    cfg.set("exists", "yes");
    
    ASSERT_TRUE(cfg.has("exists"));
    ASSERT_FALSE(cfg.has("not_exists"));
}

// ============================================================================
// PhaseType Tests
// ============================================================================

TEST(PhaseType, ToString) {
    ASSERT_EQ(phaseTypeToString(PhaseType::P), "P");
    ASSERT_EQ(phaseTypeToString(PhaseType::S), "S");
    ASSERT_EQ(phaseTypeToString(PhaseType::Pn), "Pn");
    ASSERT_EQ(phaseTypeToString(PhaseType::Unknown), "?");
}

TEST(PhaseType, FromString) {
    ASSERT_TRUE(stringToPhaseType("P") == PhaseType::P);
    ASSERT_TRUE(stringToPhaseType("S") == PhaseType::S);
    ASSERT_TRUE(stringToPhaseType("Pn") == PhaseType::Pn);
    ASSERT_TRUE(stringToPhaseType("xxx") == PhaseType::Unknown);
}

// ============================================================================
// Pick Tests
// ============================================================================

TEST(Pick, Weight) {
    Pick pick;
    
    pick.quality = PickQuality::Impulsive;
    ASSERT_NEAR(pick.weight(), 0.0, 1e-10);
    
    pick.quality = PickQuality::Emergent;
    ASSERT_NEAR(pick.weight(), 1.0, 1e-10);
    
    pick.quality = PickQuality::Rejected;
    ASSERT_NEAR(pick.weight(), 4.0, 1e-10);
}

TEST(Pick, LocationWeight) {
    Pick pick;
    
    pick.quality = PickQuality::Impulsive;
    ASSERT_NEAR(pick.locationWeight(), 1.0, 1e-10);
    
    pick.quality = PickQuality::Rejected;
    ASSERT_NEAR(pick.locationWeight(), 0.0, 1e-10);
}

// ============================================================================
// Origin Tests
// ============================================================================

TEST(Origin, QualityCode) {
    Origin origin;
    
    // Quality A: gap < 90, phases >= 8, rms < 0.5
    origin.gap = 60;
    origin.phase_count = 10;
    origin.rms = 0.3;
    ASSERT_EQ(origin.qualityCode(), 'A');
    
    // Quality B
    origin.gap = 120;
    origin.phase_count = 7;
    origin.rms = 0.7;
    ASSERT_EQ(origin.qualityCode(), 'B');
    
    // Quality D
    origin.gap = 250;
    origin.phase_count = 3;
    origin.rms = 5.0;
    ASSERT_EQ(origin.qualityCode(), 'D');
}

// ============================================================================
// Event Tests
// ============================================================================

TEST(Event, CreateWithId) {
    Event event("test_event_001");
    ASSERT_EQ(event.id(), "test_event_001");
}

TEST(Event, AutoId) {
    Event event1;
    Event event2;
    
    ASSERT_NE(event1.id(), event2.id());
    ASSERT_FALSE(event1.id().empty());
}

TEST(Event, AddOrigin) {
    Event event;
    
    Origin origin1;
    origin1.location = GeoPoint(34.0, -118.0, 10.0);
    origin1.rms = 0.5;
    
    Origin origin2;
    origin2.location = GeoPoint(34.1, -118.1, 12.0);
    origin2.rms = 0.3;
    
    event.addOrigin(origin1);
    event.addOrigin(origin2);
    
    ASSERT_EQ(event.origins().size(), 2u);
    
    // Preferred origin is the last one
    ASSERT_NEAR(event.preferredOrigin().rms, 0.3, 1e-10);
}

TEST(Event, AddMagnitude) {
    Event event;
    
    event.addMagnitude(Magnitude(MagnitudeType::ML, 4.5, 0.2, 5));
    event.addMagnitude(Magnitude(MagnitudeType::Mw, 4.3, 0.1, 3));
    
    ASSERT_EQ(event.magnitudes().size(), 2u);
    
    // Mw is preferred over ML
    ASSERT_TRUE(event.preferredMagnitude().type == MagnitudeType::Mw);
}

TEST(Event, ConvenienceAccessors) {
    Event event;
    
    Origin origin;
    origin.location = GeoPoint(34.5, -117.5, 15.0);
    origin.time = std::chrono::system_clock::now();
    event.addOrigin(origin);
    
    event.addMagnitude(Magnitude(MagnitudeType::ML, 3.8));
    
    ASSERT_NEAR(event.latitude(), 34.5, 1e-10);
    ASSERT_NEAR(event.longitude(), -117.5, 1e-10);
    ASSERT_NEAR(event.depth(), 15.0, 1e-10);
    ASSERT_NEAR(event.magnitude(), 3.8, 1e-10);
}
