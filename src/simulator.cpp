/**
 * RealDetect Data Playback
 *
 * Plays back real seismic waveform data from MiniSEED files
 * through the full processing pipeline.
 *
 * Usage:
 *   realdetect_sim [data_dir]
 *   realdetect_sim data/ridgecrest
 *   realdetect_sim data/dprk
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <map>
#include <set>
#include <chrono>
#include <algorithm>

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

using namespace realdetect;

void runPlayback(const std::string& data_dir) {
    std::cout << "========================================\n"
              << "RealDetect — Real Data Playback\n"
              << "========================================\n\n";

    // Load station inventory
    std::string sta_file = data_dir + "/stations.txt";
    StationInventory stations;
    if (!stations.loadFromFile(sta_file)) {
        std::cerr << "ERROR: Cannot load stations from " << sta_file << "\n";
        return;
    }
    // Add BHZ channel
    for (auto& [key, sta] : stations.stations()) {
        if (!sta->getChannel("BHZ")) {
            Channel bhz;
            bhz.code = "BHZ";
            bhz.sample_rate = 40.0;
            bhz.dip = -90;
            sta->addChannel(bhz);
        }
    }
    std::cout << "Loaded " << stations.size() << " stations\n";

    // Load MiniSEED waveforms
    std::string mseed_file = data_dir + "/waveforms.mseed";
    MiniSeedReader reader;
    if (!reader.open(mseed_file)) {
        std::cerr << "ERROR: Cannot load MiniSEED from " << mseed_file << "\n";
        return;
    }
    auto waveforms_vec = reader.toWaveforms();

    std::map<StreamID, WaveformPtr> waveforms;
    for (auto& wf : waveforms_vec) {
        if (wf->streamId().channel == "BHZ")
            waveforms[wf->streamId()] = wf;
    }
    std::cout << "Loaded " << waveforms.size() << " waveforms from MiniSEED\n\n";

    if (waveforms.empty()) {
        std::cerr << "ERROR: No BHZ waveforms found\n";
        return;
    }

    // Use default velocity model
    VelocityModel1D velocity_model = VelocityModel1D::simpleThreeLayer();

    // Pick phases
    std::cout << "=== Phase Picking ===\n";
    STALTAPicker picker;
    picker.setParameter("sta_length", 0.5);
    picker.setParameter("lta_length", 10.0);
    picker.setParameter("trigger_ratio", 3.0);

    std::vector<PickPtr> all_picks;
    for (const auto& [id, wf] : waveforms) {
        auto picks = picker.pick(*wf);
        for (const auto& pr : picks) {
            auto pick = std::make_shared<Pick>();
            pick->stream_id = id;
            pick->time = pr.time;
            pick->phase_type = pr.phase_type;
            pick->snr = pr.snr;
            pick->amplitude = pr.amplitude;
            pick->is_automatic = true;
            pick->method = "STA/LTA";
            pick->quality = pr.snr > 8.0 ? PickQuality::Impulsive : PickQuality::Emergent;
            std::cout << "Pick: " << pick->stream_id.toString()
                      << " " << phaseTypeToString(pick->phase_type)
                      << " SNR=" << std::setprecision(1) << pick->snr << "\n";
            all_picks.push_back(pick);
        }
    }

    // Keep first arrival per station
    std::vector<PickPtr> p_picks;
    std::set<std::string> seen;
    std::sort(all_picks.begin(), all_picks.end(),
              [](auto& a, auto& b){ return a->time < b->time; });
    for (auto& pk : all_picks) {
        std::string key = pk->stream_id.network + "." + pk->stream_id.station;
        if (seen.insert(key).second) {
            pk->phase_type = PhaseType::P;
            p_picks.push_back(pk);
        }
    }
    std::cout << "\n" << p_picks.size() << " P-wave first arrivals\n";

    // Associate
    std::cout << "\n=== Event Association ===\n";
    PhaseAssociator assoc;
    assoc.setStations(stations);
    assoc.setVelocityModel(velocity_model);
    assoc.setMinStations(3);
    assoc.setMinPhases(4);
    for (const auto& pick : p_picks) assoc.addPick(pick);
    auto events = assoc.process();
    std::cout << "Found " << events.size() << " events\n";

    // Locate
    std::cout << "\n=== Event Location ===\n";
    std::vector<PickPtr> loc_picks = p_picks;
    if (!events.empty()) {
        auto event = events[0];
        loc_picks.clear();
        for (const auto& arr : event->preferredOrigin().arrivals)
            loc_picks.push_back(arr.pick);
    }

    GridSearchLocator gs;
    gs.setVelocityModel(velocity_model);
    gs.setParameter("grid_spacing", 5.0);
    gs.setParameter("search_radius", 300.0);
    auto result = gs.locate(loc_picks, stations);

    if (result.converged) {
        std::cout << "\nLocated event:\n"
                  << "  Location: " << std::setprecision(3)
                  << result.origin.location.latitude << "°N, "
                  << result.origin.location.longitude << "°E\n"
                  << "  Depth: " << std::setprecision(1) << result.origin.location.depth << " km\n"
                  << "  RMS: " << result.origin.rms << " s\n"
                  << "  Phases: " << result.origin.phase_count << "\n"
                  << "  Stations: " << result.origin.station_count << "\n"
                  << "  Quality: " << result.origin.qualityCode() << "\n";
    } else {
        std::cout << "Location failed to converge\n";
    }

    // Magnitude
    std::cout << "\n=== Magnitude Calculation ===\n";
    LocalMagnitude ml_calc;
    auto ml = ml_calc.calculate(result.origin, waveforms, stations);
    if (ml.station_count > 0) {
        std::cout << "  ML: " << std::setprecision(2) << ml.value
                  << " ± " << ml.uncertainty
                  << " (" << ml.station_count << " stations)\n";
    } else {
        std::cout << "  No magnitude calculated\n";
    }

    std::cout << "\n========================================\n"
              << "Playback complete\n"
              << "========================================\n";
}

void printUsage(const char* progname) {
    std::cout << "RealDetect Data Playback\n\n"
              << "Plays back real seismic data through the processing pipeline.\n\n"
              << "Usage: " << progname << " [data_dir]\n\n"
              << "Arguments:\n"
              << "  data_dir    Directory containing waveforms.mseed and stations.txt\n"
              << "              (default: data/ridgecrest)\n\n"
              << "Examples:\n"
              << "  " << progname << " data/ridgecrest\n"
              << "  " << progname << " data/dprk\n";
}

int main(int argc, char* argv[]) {
    if (argc > 1 && (std::string(argv[1]) == "-h" || std::string(argv[1]) == "--help")) {
        printUsage(argv[0]);
        return 0;
    }

    std::string data_dir = "data/ridgecrest";
    if (argc > 1) data_dir = argv[1];

    runPlayback(data_dir);
    return 0;
}
