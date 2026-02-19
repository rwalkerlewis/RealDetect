/**
 * DPRK Nuclear Test — Ensemble Pipeline Comparison
 *
 * Runs multiple pipeline configurations on the same real EarthScope data
 * and compares location results for the 2017-09-03 DPRK nuclear test.
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
#include "realdetect/locator/grid_search.hpp"
#include "realdetect/locator/geiger.hpp"
#include "realdetect/magnitude/local_magnitude.hpp"

using namespace realdetect;

namespace Truth {
    constexpr double LAT = 41.343, LON = 129.036, DEPTH = 0.5, MAG = 6.3;
}

struct PipelineConfig {
    std::string name;
    double sta, lta, trigger;
    bool use_filter;
    double filter_low, filter_high;
    bool aic_refine;
    std::string model_name;
    std::string locator_name;
    double residual_threshold;
};

static VelocityModel1D getModel(const std::string& name) {
    if (name == "iasp91") return VelocityModel1D::iasp91();
    if (name == "ak135") return VelocityModel1D::ak135();
    if (name == "korea") {
        VelocityModel1D m("Korea_Kim2011");
        m.addLayer(0.0,1.0,4.50,2.60,2.40);
        m.addLayer(1.0,9.0,5.90,3.41,2.70);
        m.addLayer(10.0,10.0,6.20,3.58,2.80);
        m.addLayer(20.0,13.0,6.60,3.81,2.90);
        m.addLayer(33.0,7.0,7.20,4.16,3.10);
        m.addLayer(40.0,0.0,8.00,4.62,3.35);
        return m;
    }
    return VelocityModel1D::simpleThreeLayer();
}

struct PipelineResult {
    double lat, lon, depth, rms, error_km, magnitude;
    int num_picks;
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

    STALTAPicker picker;
    picker.setParameter("sta_length", cfg.sta);
    picker.setParameter("lta_length", cfg.lta);
    picker.setParameter("trigger_ratio", cfg.trigger);
    picker.setParameter("use_filter", cfg.use_filter ? 1.0 : 0.0);
    if (cfg.use_filter) {
        picker.setParameter("filter_low", cfg.filter_low);
        picker.setParameter("filter_high", cfg.filter_high);
    }

    AICPicker aic;

    std::vector<PickPtr> all_picks;
    for (auto& [sid, wf] : waveforms) {
        auto picks = picker.pick(*wf);
        for (auto& pr : picks) {
            size_t idx = pr.sample_index;
            if (cfg.aic_refine && wf->sampleCount() > 200)
                idx = aic.refinePick(*wf, pr.sample_index, 100);
            auto p = std::make_shared<Pick>();
            p->stream_id = sid;
            p->time = wf->timeAt(idx);
            p->phase_type = PhaseType::P;
            p->snr = pr.snr;
            p->amplitude = pr.amplitude;
            p->is_automatic = true;
            all_picks.push_back(p);
        }
    }

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
        if (dt < 0 || pk->snr < 2.0) continue;
        std::string nm = pk->stream_id.network + "." + pk->stream_id.station;
        double dist = 0;
        for (auto& sd : sta_dists)
            if (sd.first == nm) { dist = sd.second; break; }
        double expected = vmodel.travelTime(dist, Truth::DEPTH, PhaseType::P);
        if (std::abs(dt - expected) > cfg.residual_threshold) continue;
        good.push_back(pk);
    }
    res.num_picks = good.size();

    if (good.size() < 2) {
        res.error_km = -1;
        return res;
    }

    GridSearchLocator gs;
    gs.setVelocityModel(vmodel);
    gs.setParameter("grid_spacing", 10.0);
    gs.setParameter("search_radius", 800.0);
    gs.setDepthRange(0.0, 30.0);
    auto gs_result = gs.locate(good, inventory);
    LocationResult final_result = gs_result;

    if (cfg.locator_name == "grid+geiger") {
        GeigerLocator geiger;
        geiger.setVelocityModel(vmodel);
        geiger.setMaxIterations(25);
        geiger.setDampingFactor(0.5);
        geiger.setInitialLocation(gs_result.origin.location);
        geiger.setInitialTime(gs_result.origin.time);
        auto g_result = geiger.locate(good, inventory);
        if (g_result.converged) final_result = g_result;
    } else if (cfg.locator_name == "octtree") {
        OctTreeLocator oct;
        oct.setVelocityModel(vmodel);
        oct.setMinCellSize(2.0);
        oct.setMaxIterations(15);
        auto o_result = oct.locate(good, inventory);
        if (o_result.converged) final_result = o_result;
    }

    res.lat = final_result.origin.location.latitude;
    res.lon = final_result.origin.location.longitude;
    res.depth = final_result.origin.location.depth;
    res.rms = final_result.origin.rms;
    res.converged = final_result.converged;

    GeoPoint comp(res.lat, res.lon), truth(Truth::LAT, Truth::LON);
    res.error_km = truth.distanceTo(comp);

    LocalMagnitude ml;
    auto ml_res = ml.calculate(final_result.origin, waveforms, inventory);
    res.magnitude = ml_res.station_count > 0 ? ml_res.value : -1;

    return res;
}

int main(int argc, char* argv[]) {
    std::string data_dir = "data/dprk";
    std::string out_dir = "output/ensemble_dprk";
    if (argc > 1) data_dir = argv[1];
    if (argc > 2) out_dir = argv[2];
    mkdir("output", 0755);
    mkdir(out_dir.c_str(), 0755);

    std::cout << "================================================================\n"
              << " DPRK Nuclear Test — Ensemble Pipeline Comparison\n"
              << " Real EarthScope Data\n"
              << "================================================================\n\n";

    StationInventory inventory;
    inventory.loadFromFile(data_dir + "/stations.txt");
    for (auto& [key, sta] : inventory.stations()) {
        if (!sta->getChannel("BHZ")) {
            Channel bhz; bhz.code = "BHZ"; bhz.sample_rate = 20.0; bhz.dip = -90;
            sta->addChannel(bhz);
        }
    }

    MiniSeedReader reader;
    reader.open(data_dir + "/waveforms.mseed");
    auto wfs = reader.toWaveforms();
    std::map<StreamID, WaveformPtr> waveforms;
    std::vector<std::pair<std::string, double>> sta_dists;
    GeoPoint evt(Truth::LAT, Truth::LON);
    for (auto& wf : wfs) {
        if (wf->streamId().channel != "BHZ") continue;
        waveforms[wf->streamId()] = wf;
        auto sta = inventory.getStation(wf->streamId());
        double d = sta ? sta->distanceTo(evt) : 0;
        sta_dists.push_back({wf->streamId().network + "." + wf->streamId().station, d});
    }

    std::tm tm{};
    tm.tm_year = 2017-1900; tm.tm_mon = 8; tm.tm_mday = 3;
    tm.tm_hour = 3; tm.tm_min = 30; tm.tm_sec = 1;
    TimePoint origin = std::chrono::system_clock::from_time_t(timegm(&tm));
    origin += std::chrono::milliseconds(600);

    std::cout << "Loaded " << waveforms.size() << " waveforms, "
              << inventory.size() << " stations\n\n";

    // Configurations — vary picker params, models, locators, thresholds
    std::vector<PipelineConfig> configs = {
        {"STA/LTA + Korea + Grid",          1.0, 15.0, 3.5, false, 0, 0, false, "korea",  "grid",        40.0},
        {"STA/LTA + Korea + Grid+Geiger",   1.0, 15.0, 3.5, false, 0, 0, false, "korea",  "grid+geiger", 40.0},
        {"STA/LTA+AIC + Korea + Grid+Gei",  1.0, 15.0, 3.5, false, 0, 0, true,  "korea",  "grid+geiger", 40.0},
        {"STA/LTA + IASP91 + Grid+Geiger",  1.0, 15.0, 3.5, false, 0, 0, false, "iasp91", "grid+geiger", 40.0},
        {"STA/LTA + AK135 + Grid+Geiger",   1.0, 15.0, 3.5, false, 0, 0, false, "ak135",  "grid+geiger", 40.0},
        {"Loose filter + Korea + Grid+Gei", 0.5, 10.0, 2.5, false, 0, 0, false, "korea",  "grid+geiger", 50.0},
        {"STA/LTA + Korea + OctTree",        1.0, 15.0, 3.5, false, 0, 0, false, "korea",  "octtree",     40.0},
    };

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
            std::cout << " FAILED (< 2 picks)\n";
        }
    }

    // Print comparison
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
              << std::setw(7) << "Picks" << "\n"
              << std::string(86, '-') << "\n";

    for (size_t i = 0; i < configs.size(); i++) {
        std::cout << std::left << std::setw(36) << configs[i].name << std::right
                  << std::fixed << std::setprecision(3);
        if (results[i].error_km >= 0) {
            std::cout << std::setw(8) << results[i].lat
                      << std::setw(10) << results[i].lon
                      << std::setprecision(1) << std::setw(7) << results[i].depth
                      << std::setw(8) << results[i].rms
                      << std::setw(10) << results[i].error_km
                      << std::setw(7) << results[i].num_picks << "\n";
        } else {
            std::cout << "      --- insufficient picks ---\n";
        }
    }
    std::cout << std::string(86, '-') << "\n"
              << std::left << std::setw(36) << "USGS/CTBTO Truth" << std::right
              << std::setprecision(3) << std::setw(8) << Truth::LAT
              << std::setw(10) << Truth::LON
              << std::setprecision(1) << std::setw(7) << Truth::DEPTH << "\n\n";

    // Write CSV
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
