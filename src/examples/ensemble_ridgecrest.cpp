/**
 * Ridgecrest M7.1 — Ensemble Pipeline Comparison
 *
 * Runs multiple pipeline configurations on the same real EarthScope data
 * and compares location results. Demonstrates how different choices of
 * picker, velocity model, and locator affect the result.
 *
 * Configurations:
 *   1. STA/LTA + SoCal + GridSearch
 *   2. STA/LTA + SoCal + GridSearch+Geiger (two-pass)
 *   3. STA/LTA + IASP91 + GridSearch+Geiger
 *   4. STA/LTA + AK135 + GridSearch+Geiger
 *   5. STA/LTA+AIC + SoCal + GridSearch+Geiger
 *   6. STA/LTA (tight BP) + SoCal + GridSearch+Geiger
 *   7. STA/LTA + SoCal + OctTree
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <vector>
#include <map>
#include <set>
#include <chrono>
#include <algorithm>
#include <sys/stat.h>

#include "realdetect/core/types.hpp"
#include "realdetect/core/station.hpp"
#include "realdetect/core/event.hpp"
#include "realdetect/core/velocity_model.hpp"
#include "realdetect/core/waveform.hpp"
#include "realdetect/core/miniseed.hpp"
#include "realdetect/picker/stalta_picker.hpp"
#include "realdetect/picker/aic_picker.hpp"
#include "realdetect/picker/polarization_picker.hpp"
#include "realdetect/locator/grid_search.hpp"
#include "realdetect/locator/geiger.hpp"
#include "realdetect/magnitude/local_magnitude.hpp"

using namespace realdetect;

namespace Truth {
    constexpr double LAT = 35.770, LON = -117.599, DEPTH = 8.0, MAG = 7.1;
}

// Pipeline configuration
struct PipelineConfig {
    std::string name;
    double sta, lta, trigger;
    double filter_low, filter_high;
    bool use_filter;
    bool aic_refine;
    bool use_polarization;      // use 3C polarization-coupled picker
    std::string model_name;
    std::string locator_name;   // "grid", "grid+geiger", "octtree"
};

// Get velocity model by name
static VelocityModel1D getModel(const std::string& name) {
    if (name == "iasp91") return VelocityModel1D::iasp91();
    if (name == "ak135") return VelocityModel1D::ak135();
    if (name == "socal") {
        VelocityModel1D m("SoCal_HK");
        m.addLayer(0.0,5.5,5.50,3.18,2.40);
        m.addLayer(5.5,10.5,6.30,3.64,2.67);
        m.addLayer(16.0,16.0,6.70,3.87,2.80);
        m.addLayer(32.0,0.0,7.80,4.50,3.30);
        return m;
    }
    return VelocityModel1D::simpleThreeLayer();
}

// Run a single pipeline configuration
struct PipelineResult {
    double lat, lon, depth, rms, error_km, magnitude;
    int num_picks, num_phases;
    bool converged;
};

static PipelineResult runPipeline(
    const PipelineConfig& cfg,
    const StationInventory& inventory,
    const std::map<StreamID, WaveformPtr>& waveforms,
    TimePoint origin,
    const std::vector<std::pair<std::string, double>>& sta_dists)
{
    PipelineResult res{};
    auto vmodel = getModel(cfg.model_name);

    // Group waveforms by station for 3C access
    struct Station3C { WaveformPtr z, n, e; };
    std::map<std::string, Station3C> station_wfs;
    for (auto& [sid, wf] : waveforms) {
        std::string key = sid.network + "." + sid.station;
        if (sid.channel == "BHZ") station_wfs[key].z = wf;
        else if (sid.channel == "BHN" || sid.channel == "BH1") station_wfs[key].n = wf;
        else if (sid.channel == "BHE" || sid.channel == "BH2") station_wfs[key].e = wf;
    }

    // Pick — either STA/LTA or Polarization
    std::vector<PickPtr> all_picks;

    if (cfg.use_polarization) {
        PolarizationPicker pol;
        pol.setParameter("sta_length", cfg.sta);
        pol.setParameter("lta_length", cfg.lta);
        pol.setParameter("trigger_ratio", cfg.trigger);
        pol.setParameter("filter_low", cfg.filter_low);
        pol.setParameter("filter_high", cfg.filter_high);
        pol.setParameter("use_filter", cfg.use_filter ? 1.0 : 0.0);

        for (auto& [key, s3c] : station_wfs) {
            std::vector<PickResult> picks;
            if (s3c.z && s3c.n && s3c.e) {
                picks = pol.pick3C(*s3c.z, *s3c.n, *s3c.e);
            } else if (s3c.z) {
                picks = pol.pick(*s3c.z);  // fallback 1C
            }
            for (auto& pr : picks) {
                auto p = std::make_shared<Pick>();
                p->stream_id = s3c.z ? s3c.z->streamId() : StreamID();
                p->time = pr.time;
                p->phase_type = PhaseType::P;
                p->snr = pr.snr;
                p->amplitude = pr.amplitude;
                p->is_automatic = true;
                all_picks.push_back(p);
            }
        }
    } else {
        STALTAPicker picker;
        picker.setParameter("sta_length", cfg.sta);
        picker.setParameter("lta_length", cfg.lta);
        picker.setParameter("trigger_ratio", cfg.trigger);
        picker.setParameter("filter_low", cfg.filter_low);
        picker.setParameter("filter_high", cfg.filter_high);
        picker.setParameter("use_filter", cfg.use_filter ? 1.0 : 0.0);

        AICPicker aic;

        for (auto& [sid, wf] : waveforms) {
            if (sid.channel != "BHZ") continue;  // only pick on Z
            auto picks = picker.pick(*wf);
            for (auto& pr : picks) {
                size_t pick_idx = pr.sample_index;
                if (cfg.aic_refine && wf->sampleCount() > 200)
                    pick_idx = aic.refinePick(*wf, pr.sample_index, 100);
                auto p = std::make_shared<Pick>();
                p->stream_id = sid;
                p->time = wf->timeAt(pick_idx);
                p->phase_type = PhaseType::P;
                p->snr = pr.snr;
                p->amplitude = pr.amplitude;
                p->is_automatic = true;
                all_picks.push_back(p);
            }
        }
    }

    // First arrival per station
    std::sort(all_picks.begin(), all_picks.end(),
              [](auto& a, auto& b) { return a->time < b->time; });
    std::vector<PickPtr> p_picks;
    std::set<std::string> seen;
    for (auto& pk : all_picks) {
        std::string key = pk->stream_id.network + "." + pk->stream_id.station;
        if (seen.insert(key).second) p_picks.push_back(pk);
    }

    // Quality filter
    std::vector<PickPtr> good;
    for (auto& pk : p_picks) {
        double dt = std::chrono::duration_cast<std::chrono::microseconds>(
            pk->time - origin).count() / 1e6;
        if (dt < -1.0 || pk->snr < 3.0) continue;
        std::string nm = pk->stream_id.network + "." + pk->stream_id.station;
        double dist = 0;
        for (auto& sd : sta_dists)
            if (sd.first == nm) { dist = sd.second; break; }
        double expected = vmodel.travelTime(dist, Truth::DEPTH, PhaseType::P);
        if (std::abs(dt - expected) > 5.0) continue;
        good.push_back(pk);
    }
    res.num_picks = good.size();

    if (good.size() < 3) {
        res.error_km = -1;
        return res;
    }

    // Locate
    GridSearchLocator gs;
    gs.setVelocityModel(vmodel);
    gs.setParameter("grid_spacing", 5.0);
    gs.setParameter("search_radius", 300.0);
    auto gs_result = gs.locate(good, inventory);

    LocationResult final_result = gs_result;

    if (cfg.locator_name == "grid+geiger") {
        GeigerLocator geiger;
        geiger.setVelocityModel(vmodel);
        geiger.setMaxIterations(30);
        geiger.setConvergenceThreshold(0.001);
        geiger.setDampingFactor(0.4);
        geiger.setInitialLocation(gs_result.origin.location);
        geiger.setInitialTime(gs_result.origin.time);
        auto g_result = geiger.locate(good, inventory);
        if (g_result.converged) final_result = g_result;
    } else if (cfg.locator_name == "octtree") {
        OctTreeLocator oct;
        oct.setVelocityModel(vmodel);
        oct.setMinCellSize(1.0);
        oct.setMaxIterations(20);
        auto o_result = oct.locate(good, inventory);
        if (o_result.converged) final_result = o_result;
    }

    res.lat = final_result.origin.location.latitude;
    res.lon = final_result.origin.location.longitude;
    res.depth = final_result.origin.location.depth;
    res.rms = final_result.origin.rms;
    res.num_phases = final_result.origin.phase_count;
    res.converged = final_result.converged;

    GeoPoint comp(res.lat, res.lon);
    GeoPoint truth(Truth::LAT, Truth::LON);
    res.error_km = truth.distanceTo(comp);

    // Magnitude
    LocalMagnitude ml;
    auto ml_res = ml.calculate(final_result.origin, waveforms, inventory);
    res.magnitude = ml_res.station_count > 0 ? ml_res.value : -1;

    return res;
}

int main(int argc, char* argv[]) {
    std::string data_dir = "data/ridgecrest";
    std::string out_dir = "output/ensemble_ridgecrest";
    if (argc > 1) data_dir = argv[1];
    if (argc > 2) out_dir = argv[2];
    mkdir("output", 0755);
    mkdir(out_dir.c_str(), 0755);

    std::cout << "================================================================\n"
              << " Ridgecrest M7.1 — Ensemble Pipeline Comparison\n"
              << " Real EarthScope Data\n"
              << "================================================================\n\n";

    // Load data
    StationInventory inventory;
    inventory.loadFromFile(data_dir + "/stations.txt");
    for (auto& [key, sta] : inventory.stations()) {
        if (!sta->getChannel("BHZ")) {
            Channel bhz; bhz.code = "BHZ"; bhz.sample_rate = 40.0; bhz.dip = -90;
            sta->addChannel(bhz);
        }
    }

    MiniSeedReader reader;
    reader.open(data_dir + "/waveforms.mseed");
    auto wfs = reader.toWaveforms();
    std::map<StreamID, WaveformPtr> waveforms;
    std::vector<std::pair<std::string, double>> sta_dists;
    std::set<std::string> seen_stas;
    GeoPoint evt(Truth::LAT, Truth::LON);
    for (auto& wf : wfs) {
        waveforms[wf->streamId()] = wf;  // keep ALL channels (Z, N, E)
        std::string key = wf->streamId().network + "." + wf->streamId().station;
        if (wf->streamId().channel == "BHZ" && seen_stas.insert(key).second) {
            auto sta = inventory.getStation(wf->streamId());
            double d = sta ? sta->distanceTo(evt) : 0;
            sta_dists.push_back({key, d});
        }
    }

    // Origin time
    std::tm tm{};
    tm.tm_year = 2019-1900; tm.tm_mon = 6; tm.tm_mday = 6;
    tm.tm_hour = 3; tm.tm_min = 19; tm.tm_sec = 53;
    TimePoint origin = std::chrono::system_clock::from_time_t(timegm(&tm));
    origin += std::chrono::milliseconds(40);

    std::cout << "Loaded " << waveforms.size() << " waveforms, "
              << inventory.size() << " stations\n\n";

    // Define configurations
    //                                          sta  lta  trig  flo  fhi   filt  aic   pol    model    locator
    std::vector<PipelineConfig> configs = {
        {"STA/LTA + SoCal + Grid",              0.3, 8.0, 2.5, 1.0, 15.0, false, false, false, "socal",  "grid"},
        {"STA/LTA + SoCal + Grid+Geiger",       0.3, 8.0, 2.5, 1.0, 15.0, false, false, false, "socal",  "grid+geiger"},
        {"STA/LTA+AIC + SoCal + Grid+Geiger",   0.3, 8.0, 2.5, 1.0, 15.0, false, true,  false, "socal",  "grid+geiger"},
        {"STA/LTA + IASP91 + Grid+Geiger",      0.3, 8.0, 2.5, 1.0, 15.0, false, false, false, "iasp91", "grid+geiger"},
        {"STA/LTA + AK135 + Grid+Geiger",       0.3, 8.0, 2.5, 1.0, 15.0, false, false, false, "ak135",  "grid+geiger"},
        {"STA/LTA+Pol + SoCal + Grid+Geiger",   0.3, 8.0, 1.0, 1.0, 15.0, false,  false, true,  "socal",  "grid+geiger"},
        {"STA/LTA+Pol + IASP91 + Grid+Geiger",  0.3, 8.0, 1.0, 1.0, 15.0, false,  false, true,  "iasp91", "grid+geiger"},
        {"STA/LTA + SoCal + OctTree",            0.3, 8.0, 2.5, 1.0, 15.0, false, false, false, "socal",  "octtree"},
    };

    // Run all configurations
    std::vector<PipelineResult> results;
    for (auto& cfg : configs) {
        std::cout << "Running: " << cfg.name << "..." << std::flush;
        auto res = runPipeline(cfg, inventory, waveforms, origin, sta_dists);
        results.push_back(res);
        if (res.error_km >= 0) {
            std::cout << " error=" << std::fixed << std::setprecision(1)
                      << res.error_km << " km, RMS=" << res.rms << " s, "
                      << res.num_picks << " picks\n";
        } else {
            std::cout << " FAILED (< 3 picks)\n";
        }
    }

    // Print comparison table
    std::cout << "\n================================================================\n"
              << " ENSEMBLE RESULTS COMPARISON\n"
              << " Ground Truth: " << Truth::LAT << "°N, " << Truth::LON << "°E, "
              << Truth::DEPTH << " km\n"
              << "================================================================\n\n"
              << std::left << std::setw(36) << "Configuration"
              << std::right << std::setw(8) << "Lat"
              << std::setw(10) << "Lon"
              << std::setw(7) << "Dep"
              << std::setw(8) << "RMS"
              << std::setw(10) << "Error"
              << std::setw(7) << "Picks"
              << std::setw(8) << "ML" << "\n"
              << std::string(94, '-') << "\n";

    for (size_t i = 0; i < configs.size(); i++) {
        auto& cfg = configs[i];
        auto& res = results[i];
        std::cout << std::left << std::setw(36) << cfg.name << std::right
                  << std::fixed << std::setprecision(3);
        if (res.error_km >= 0) {
            std::cout << std::setw(8) << res.lat
                      << std::setw(10) << res.lon
                      << std::setprecision(1) << std::setw(7) << res.depth
                      << std::setw(8) << res.rms
                      << std::setw(10) << res.error_km
                      << std::setw(7) << res.num_picks
                      << std::setprecision(2) << std::setw(8)
                      << (res.magnitude > 0 ? res.magnitude : 0.0) << "\n";
        } else {
            std::cout << "      --- insufficient picks ---\n";
        }
    }

    std::cout << std::string(94, '-') << "\n"
              << std::left << std::setw(36) << "USGS Truth" << std::right
              << std::setprecision(3) << std::setw(8) << Truth::LAT
              << std::setw(10) << Truth::LON
              << std::setprecision(1) << std::setw(7) << Truth::DEPTH
              << "\n\n";

    // Write comparison CSV
    {
        std::ofstream f(out_dir + "/comparison.csv");
        f << "config,latitude,longitude,depth_km,rms_s,error_km,num_picks,magnitude\n";
        for (size_t i = 0; i < configs.size(); i++) {
            f << configs[i].name << "," << results[i].lat << "," << results[i].lon
              << "," << results[i].depth << "," << results[i].rms
              << "," << results[i].error_km << "," << results[i].num_picks
              << "," << results[i].magnitude << "\n";
        }
        f << "USGS_Truth," << Truth::LAT << "," << Truth::LON << ","
          << Truth::DEPTH << ",0,0,0," << Truth::MAG << "\n";
    }

    std::cout << "Comparison saved to " << out_dir << "/comparison.csv\n";
    return 0;
}
