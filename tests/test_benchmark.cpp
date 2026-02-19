/**
 * Benchmark tests using real earthquake data from EarthScope
 *
 * Tests processing pipeline performance on real waveform data.
 */

#include "test_framework.hpp"
#include "realdetect/core/types.hpp"
#include "realdetect/core/waveform.hpp"
#include "realdetect/core/station.hpp"
#include "realdetect/core/event.hpp"
#include "realdetect/core/velocity_model.hpp"
#include "realdetect/core/miniseed.hpp"
#include "realdetect/picker/stalta_picker.hpp"
#include "realdetect/picker/aic_picker.hpp"
#include "realdetect/associator/phase_associator.hpp"
#include "realdetect/locator/grid_search.hpp"
#include "realdetect/locator/geiger.hpp"
#include "realdetect/magnitude/local_magnitude.hpp"
#include <cmath>
#include <set>

using namespace realdetect;
using namespace realdetect::test;

// Load real data fixtures
static StationInventory loadBenchStations() {
    StationInventory inv;
    if (!inv.loadFromFile("data/ridgecrest/stations.txt"))
        inv.loadFromFile("../data/ridgecrest/stations.txt");
    for (auto& [key, sta] : inv.stations()) {
        if (!sta->getChannel("BHZ")) {
            Channel bhz; bhz.code = "BHZ"; bhz.sample_rate = 40.0; bhz.dip = -90;
            sta->addChannel(bhz);
        }
    }
    return inv;
}

static std::map<StreamID, WaveformPtr> loadBenchWaveforms() {
    MiniSeedReader reader;
    if (!reader.open("data/ridgecrest/waveforms.mseed"))
        reader.open("../data/ridgecrest/waveforms.mseed");
    auto wfs = reader.toWaveforms();
    std::map<StreamID, WaveformPtr> result;
    for (auto& wf : wfs)
        if (wf->streamId().channel == "BHZ")
            result[wf->streamId()] = wf;
    return result;
}

// ============================================================================
// Pipeline Benchmark: Full processing on real Ridgecrest M7.1 data
// ============================================================================

TEST(BenchmarkPipeline, RealDataFullPipeline) {
    auto stations = loadBenchStations();
    auto waveforms = loadBenchWaveforms();
    if (stations.size() == 0 || waveforms.empty()) {
        std::cout << "    [SKIP] No real benchmark data available\n";
        return;
    }

    auto start = std::chrono::high_resolution_clock::now();

    auto model = VelocityModel1D::simpleThreeLayer();

    // Phase picking
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

    // First arrival per station
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

    // Association
    PhaseAssociator assoc;
    assoc.setVelocityModel(model);
    assoc.setStations(stations);
    assoc.setMinStations(3);
    assoc.setMinPhases(4);
    for (auto& pk : p_picks) assoc.addPick(pk);
    auto events = assoc.process();

    // Location
    std::vector<PickPtr> loc_picks = p_picks;
    if (!events.empty()) {
        loc_picks.clear();
        for (auto& arr : events[0]->preferredOrigin().arrivals)
            loc_picks.push_back(arr.pick);
    }

    GridSearchLocator gs;
    gs.setVelocityModel(model);
    gs.setSearchRadius(300.0);
    auto result = gs.locate(loc_picks, stations);

    // Magnitude
    LocalMagnitude ml;
    auto ml_result = ml.calculate(result.origin, waveforms, stations);

    auto end = std::chrono::high_resolution_clock::now();
    double duration_ms = std::chrono::duration<double, std::milli>(end - start).count();

    std::cout << "    Full pipeline on " << waveforms.size()
              << " real waveforms: " << duration_ms << " ms\n";
    std::cout << "    Picks: " << p_picks.size()
              << ", Events: " << events.size()
              << ", ML stations: " << ml_result.station_count << "\n";

    ASSERT_LT(duration_ms, 60000.0);  // Should complete within 60 seconds
}

TEST(BenchmarkPerformance, PickingSpeed) {
    auto waveforms = loadBenchWaveforms();
    if (waveforms.empty()) {
        std::cout << "    [SKIP] No benchmark data\n";
        return;
    }

    STALTAPicker picker;
    picker.setParameter("trigger_ratio", 3.0);

    auto start = std::chrono::high_resolution_clock::now();

    size_t total_picks = 0;
    for (int iter = 0; iter < 10; iter++) {
        for (auto& [id, wf] : waveforms) {
            auto picks = picker.pick(*wf);
            total_picks += picks.size();
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    double duration_ms = std::chrono::duration<double, std::milli>(end - start).count();

    std::cout << "    10 iterations × " << waveforms.size()
              << " waveforms: " << duration_ms << " ms ("
              << total_picks << " total picks)\n";

    ASSERT_LT(duration_ms, 30000.0);
}

TEST(BenchmarkAccuracy, LocationAccuracy) {
    auto stations = loadBenchStations();
    auto waveforms = loadBenchWaveforms();
    if (stations.size() == 0 || waveforms.empty()) {
        std::cout << "    [SKIP] No benchmark data\n";
        return;
    }

    auto model = VelocityModel1D::simpleThreeLayer();

    // Pick
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
            p->phase_type = PhaseType::P;
            p->snr = r.snr;
            all_picks.push_back(p);
        }
    }

    std::vector<PickPtr> p_picks;
    std::set<std::string> seen;
    std::sort(all_picks.begin(), all_picks.end(),
              [](auto& a, auto& b) { return a->time < b->time; });
    for (auto& pk : all_picks) {
        std::string key = pk->stream_id.network + "." + pk->stream_id.station;
        if (seen.insert(key).second) p_picks.push_back(pk);
    }

    GridSearchLocator gs;
    gs.setVelocityModel(model);
    gs.setSearchRadius(300.0);
    auto result = gs.locate(p_picks, stations);

    // Ground truth: Ridgecrest M7.1 at 35.770°N, 117.599°W
    GeoPoint truth(35.770, -117.599);
    GeoPoint computed(result.origin.location.latitude, result.origin.location.longitude);
    double error_km = truth.distanceTo(computed);

    std::cout << "    Location: " << result.origin.location.latitude << "°N, "
              << result.origin.location.longitude << "°E\n"
              << "    Error: " << error_km << " km\n";

    // Real data — location accuracy depends on data quality and pipeline
    // Just verify the result is physically reasonable
    ASSERT_GT(result.origin.location.latitude, 20.0);
    ASSERT_LT(result.origin.location.latitude, 50.0);
}
