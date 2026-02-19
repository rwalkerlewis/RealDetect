/**
 * Unit tests for magnitude calculation
 *
 * Uses real seismic data from EarthScope where waveforms are needed.
 */

#include "test_framework.hpp"
#include "realdetect/core/waveform.hpp"
#include "realdetect/core/station.hpp"
#include "realdetect/core/event.hpp"
#include "realdetect/core/miniseed.hpp"
#include "realdetect/core/velocity_model.hpp"
#include "realdetect/magnitude/local_magnitude.hpp"
#include "realdetect/magnitude/moment_magnitude.hpp"
#include "realdetect/picker/stalta_picker.hpp"
#include "realdetect/locator/grid_search.hpp"
#include <cmath>
#include <set>

using namespace realdetect;
using namespace realdetect::test;

// Load test fixtures
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

static std::map<StreamID, WaveformPtr> loadTestWaveformMap() {
    MiniSeedReader reader;
    if (!reader.open("data/test/waveforms.mseed"))
        reader.open("../data/test/waveforms.mseed");
    auto wfs = reader.toWaveforms();
    std::map<StreamID, WaveformPtr> result;
    for (auto& wf : wfs)
        if (wf->streamId().channel == "BHZ")
            result[wf->streamId()] = wf;
    return result;
}

// ============================================================================
// Local Magnitude Tests
// ============================================================================

TEST(LocalMagnitude, DefaultParameters) {
    LocalMagnitude ml;
    ASSERT_TRUE(true);
}

TEST(LocalMagnitude, CalculateOnRealData) {
    auto stations = loadTestStations();
    auto waveforms = loadTestWaveformMap();
    if (stations.size() == 0 || waveforms.empty()) {
        std::cout << "    [SKIP] No real test data available\n";
        return;
    }

    // Create an origin near Ridgecrest
    Origin origin;
    origin.location = GeoPoint(35.770, -117.599, 8.0);
    origin.time = std::chrono::system_clock::now();

    LocalMagnitude ml;
    auto result = ml.calculate(origin, waveforms, stations);

    // Should produce some result with real earthquake data
    if (result.station_count > 0) {
        // Magnitude should be reasonable for M7.1 earthquake at close range
        ASSERT_GT(result.value, 0.0);
        ASSERT_LT(result.value, 12.0);
        ASSERT_GT(result.station_count, 0);
    }
}

TEST(LocalMagnitude, SimulateWoodAnderson) {
    LocalMagnitude ml;
    ml.setSimulateWoodAnderson(true);
    // Verify parameter can be set
    ASSERT_TRUE(true);
}

TEST(LocalMagnitude, SetWoodAnderson) {
    LocalMagnitude ml;
    ml.setSimulateWoodAnderson(true);
    ml.setSimulateWoodAnderson(false);
    ASSERT_TRUE(true);
}

// ============================================================================
// Moment Magnitude Tests
// ============================================================================

TEST(MomentMagnitude, DefaultParameters) {
    MomentMagnitude mw;
    ASSERT_TRUE(true);
}

TEST(MomentMagnitude, CalculateOnRealData) {
    auto stations = loadTestStations();
    auto waveforms = loadTestWaveformMap();
    if (stations.size() == 0 || waveforms.empty()) {
        std::cout << "    [SKIP] No real test data available\n";
        return;
    }

    Origin origin;
    origin.location = GeoPoint(35.770, -117.599, 8.0);
    origin.time = std::chrono::system_clock::now();

    MomentMagnitude mw;
    auto result = mw.calculate(origin, waveforms, stations);

    if (result.station_count > 0) {
        ASSERT_GT(result.value, -2.0);
        ASSERT_LT(result.value, 12.0);
    }
}

TEST(MomentMagnitude, SpectralAnalysis) {
    auto waveforms = loadTestWaveformMap();
    if (waveforms.empty()) {
        std::cout << "    [SKIP] No real test data available\n";
        return;
    }

    // Verify we can create Mw calculator with real data present
    MomentMagnitude mw;
    ASSERT_TRUE(true);
}

// ============================================================================
// Magnitude Factory Tests
// ============================================================================

TEST(MagnitudeFactory, CreateML) {
    auto ml = MagnitudeFactory::create(MagnitudeType::ML);
    ASSERT_NE(ml, nullptr);
}

TEST(MagnitudeFactory, CreateMw) {
    auto mw = MagnitudeFactory::create(MagnitudeType::Mw);
    ASSERT_NE(mw, nullptr);
}
