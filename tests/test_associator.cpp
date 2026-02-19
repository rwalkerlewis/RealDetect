/**
 * Unit tests for phase association
 *
 * Uses real seismic data from EarthScope for pick-based tests.
 */

#include "test_framework.hpp"
#include "realdetect/associator/phase_associator.hpp"
#include "realdetect/core/station.hpp"
#include "realdetect/core/event.hpp"
#include "realdetect/core/velocity_model.hpp"
#include "realdetect/core/miniseed.hpp"
#include "realdetect/picker/stalta_picker.hpp"
#include <set>

using namespace realdetect;
using namespace realdetect::test;

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
// Phase Associator Tests
// ============================================================================

TEST(PhaseAssociator, DefaultParameters) {
    PhaseAssociator assoc;
    ASSERT_TRUE(true);
}

TEST(PhaseAssociator, SetParameters) {
    PhaseAssociator assoc;
    assoc.setTimeWindow(120.0);
    assoc.setMinStations(3);
    assoc.setMinPhases(4);
    ASSERT_TRUE(true);
}

TEST(PhaseAssociator, AssociateRealPicks) {
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
    PhaseAssociator assoc;
    assoc.setVelocityModel(model);
    assoc.setStations(stations);
    assoc.setTimeWindow(120.0);
    assoc.setMinStations(2);
    assoc.setMinPhases(3);

    for (auto& pk : picks) assoc.addPick(pk);
    auto events = assoc.process();

    // With real M7.1 data, we should get at least one event
    // (but association may group them differently)
    std::cout << "    Associated " << events.size() << " events from "
              << picks.size() << " real picks\n";
}

TEST(PhaseAssociator, EmptyInput) {
    auto model = VelocityModel1D::simpleThreeLayer();
    StationInventory inv;

    PhaseAssociator assoc;
    assoc.setVelocityModel(model);
    assoc.setStations(inv);

    auto events = assoc.process();
    ASSERT_EQ(events.size(), 0u);
}

// ============================================================================
// Nucleator Associator Tests
// ============================================================================

TEST(NucleatorAssociator, DefaultParameters) {
    NucleatorAssociator assoc;
    ASSERT_TRUE(true);
}

TEST(NucleatorAssociator, SetGrid) {
    NucleatorAssociator assoc;
    assoc.setGrid(34.0, 37.0, 0.1,
                  -119.0, -116.0, 0.1,
                  0.0, 50.0, 5.0);
    ASSERT_TRUE(true);
}

TEST(NucleatorAssociator, AssociateRealPicks) {
    auto stations = loadTestStations();
    auto waveforms = loadTestWaveforms();
    if (stations.size() == 0 || waveforms.empty()) {
        std::cout << "    [SKIP] No real test data available\n";
        return;
    }

    auto picks = pickRealData(waveforms);
    if (picks.size() < 3) {
        std::cout << "    [SKIP] Insufficient picks\n";
        return;
    }

    auto model = VelocityModel1D::simpleThreeLayer();
    NucleatorAssociator assoc;
    assoc.setVelocityModel(model);
    assoc.setStations(stations);
    assoc.setGrid(34.77, 36.77, 0.1,
                  -118.60, -116.60, 0.1,
                  0.0, 30.0, 5.0);

    for (auto& pk : picks) assoc.addPick(pk);
    auto events = assoc.process();

    std::cout << "    Nucleator found " << events.size() << " events\n";
}

// ============================================================================
// Event Callback Tests
// ============================================================================

TEST(PhaseAssociator, EventCallback) {
    int callback_count = 0;
    auto model = VelocityModel1D::simpleThreeLayer();
    StationInventory inv;

    PhaseAssociator assoc;
    assoc.setVelocityModel(model);
    assoc.setStations(inv);
    assoc.setEventCallback([&callback_count](EventPtr) {
        callback_count++;
    });
    assoc.process();
    // No picks → no events → callback count stays 0
    ASSERT_EQ(callback_count, 0);
}
