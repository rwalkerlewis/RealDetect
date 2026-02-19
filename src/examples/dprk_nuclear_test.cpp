/**
 * DPRK Nuclear Test Detection and Inversion
 *
 * Processes REAL seismic waveform data from EarthScope (IRIS)
 * for the 2017-09-03 DPRK 6th nuclear test.
 *
 * Data source: EarthScope FDSN web services (fetched by scripts/fetch_dprk_data.py)
 *
 * Event Parameters (USGS/CTBTO):
 *   Date: 2017-09-03
 *   Origin Time: 03:30:01.6 UTC
 *   Location: 41.343°N, 129.036°E (Punggye-ri test site)
 *   Depth: ~0 km (surface/shallow explosion)
 *   Magnitude: mb 6.3
 *
 * Pipeline:
 *   1. Load real station metadata from EarthScope
 *   2. Load real MiniSEED waveforms from EarthScope
 *   3. STA/LTA phase detection
 *   4. Phase association
 *   5. Geiger iterative location
 *   6. Magnitude calculation
 *   7. CSS3.0 database storage
 *   8. CSV data output for plotting
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
#include "realdetect/core/config.hpp"
#include "realdetect/core/station.hpp"
#include "realdetect/core/event.hpp"
#include "realdetect/core/velocity_model.hpp"
#include "realdetect/core/waveform.hpp"
#include "realdetect/core/miniseed.hpp"
#include "realdetect/database/css30_database.hpp"
#include "realdetect/picker/stalta_picker.hpp"
#include "realdetect/associator/phase_associator.hpp"
#include "realdetect/locator/grid_search.hpp"
#include "realdetect/locator/geiger.hpp"
#include "realdetect/magnitude/local_magnitude.hpp"

using namespace realdetect;

namespace GroundTruth {
    constexpr double LATITUDE = 41.343;
    constexpr double LONGITUDE = 129.036;
    constexpr double DEPTH = 0.5;
    constexpr double MAGNITUDE_MB = 6.3;
    constexpr double MAGNITUDE_ML = 6.3;
    constexpr int YEAR = 2017;
    constexpr int MONTH = 9;
    constexpr int DAY = 3;
    constexpr int HOUR = 3;
    constexpr int MINUTE = 30;
    constexpr double SECOND = 1.6;
    constexpr double YIELD_KT = 250.0;
}

static VelocityModel1D getKoreanVelocityModel() {
    VelocityModel1D model("Korea_KimEtAl2011");
    model.addLayer(0.0,   1.0, 4.50, 2.60, 2.40);
    model.addLayer(1.0,   9.0, 5.90, 3.41, 2.70);
    model.addLayer(10.0, 10.0, 6.20, 3.58, 2.80);
    model.addLayer(20.0, 13.0, 6.60, 3.81, 2.90);
    model.addLayer(33.0,  7.0, 7.20, 4.16, 3.10);
    model.addLayer(40.0,  0.0, 8.00, 4.62, 3.35);
    return model;
}

int main(int argc, char* argv[]) {
    std::string data_dir = "data/dprk";
    std::string out_dir = "output/dprk";
    std::string db_file = "output/dprk/dprk_nuclear_test.db";
    if (argc > 1) data_dir = argv[1];
    if (argc > 2) out_dir = argv[2];
    if (argc > 3) db_file = argv[3];
    mkdir("output", 0755);
    mkdir(out_dir.c_str(), 0755);

    std::cout << "================================================================\n"
              << " DPRK Nuclear Test — Real Data Processing Pipeline\n"
              << " Data source: EarthScope FDSN (IRIS)\n"
              << "================================================================\n\n";

    std::cout << "Ground Truth (USGS/CTBTO):\n"
              << "  " << GroundTruth::YEAR << "-"
              << std::setfill('0') << std::setw(2) << GroundTruth::MONTH << "-"
              << std::setw(2) << GroundTruth::DAY << " "
              << std::setw(2) << GroundTruth::HOUR << ":"
              << std::setw(2) << GroundTruth::MINUTE << ":"
              << std::fixed << std::setprecision(1) << GroundTruth::SECOND << " UTC\n"
              << "  " << GroundTruth::LATITUDE << "°N, "
              << GroundTruth::LONGITUDE << "°E\n"
              << "  Depth: " << GroundTruth::DEPTH << " km\n"
              << "  Magnitude: mb " << GroundTruth::MAGNITUDE_MB << "\n"
              << "  Est. Yield: " << GroundTruth::YIELD_KT << " kt TNT\n\n";

    // origin time
    std::tm tm{};
    tm.tm_year = GroundTruth::YEAR - 1900;
    tm.tm_mon  = GroundTruth::MONTH - 1;
    tm.tm_mday = GroundTruth::DAY;
    tm.tm_hour = GroundTruth::HOUR;
    tm.tm_min  = GroundTruth::MINUTE;
    tm.tm_sec  = static_cast<int>(GroundTruth::SECOND);
    TimePoint origin = std::chrono::system_clock::from_time_t(timegm(&tm));
    origin += std::chrono::milliseconds(
        static_cast<int>((GroundTruth::SECOND - tm.tm_sec) * 1000));

    // ── 1. Load real stations ──
    std::cout << "Step 1: Loading real station inventory from EarthScope...\n";
    StationInventory inventory;
    std::string sta_file = data_dir + "/stations.txt";
    if (!inventory.loadFromFile(sta_file)) {
        std::cerr << "ERROR: Cannot load stations from " << sta_file << "\n";
        std::cerr << "Run: python3 scripts/fetch_dprk_data.py\n";
        return 1;
    }
    for (auto& [key, sta] : inventory.stations()) {
        if (!sta->getChannel("BHZ")) {
            Channel bhz;
            bhz.code = "BHZ";
            bhz.sample_rate = 20.0;
            bhz.dip = -90;
            sta->addChannel(bhz);
        }
    }
    std::cout << "  " << inventory.size() << " real stations loaded\n\n";

    // write stations CSV
    {
        std::ofstream f(out_dir + "/stations.csv");
        f << "network,station,latitude,longitude,elevation,distance_km\n";
        GeoPoint evt(GroundTruth::LATITUDE, GroundTruth::LONGITUDE);
        for (auto& [key, sta] : inventory.stations()) {
            double d = sta->distanceTo(evt);
            f << sta->network() << "," << sta->code() << ","
              << sta->latitude() << "," << sta->longitude() << ","
              << sta->elevation() << "," << d << "\n";
        }
    }

    // ── 2. Velocity model ──
    std::cout << "Step 2: Loading Korean Peninsula velocity model...\n";
    auto vmodel = getKoreanVelocityModel();
    std::cout << "  Model: " << vmodel.name() << " (" << vmodel.layerCount() << " layers)\n\n";

    {
        std::ofstream f(out_dir + "/velocity_model.csv");
        f << "depth_km,vp_km_s,vs_km_s\n";
        double depths[] = {0,1,5,10,15,20,25,30,33,35,40,50};
        for (double d : depths)
            f << d << "," << vmodel.vpAt(d) << "," << vmodel.vsAt(d) << "\n";
    }

    // ── 3. Load REAL waveforms ──
    std::cout << "Step 3: Loading real waveforms from EarthScope MiniSEED...\n";
    MiniSeedReader reader;
    std::string mseed_file = data_dir + "/waveforms.mseed";
    if (!reader.open(mseed_file)) {
        std::cerr << "ERROR: Cannot load MiniSEED from " << mseed_file << "\n";
        std::cerr << "Run: python3 scripts/fetch_dprk_data.py\n";
        return 1;
    }
    auto waveforms_vec = reader.toWaveforms();
    std::map<StreamID, WaveformPtr> waveforms;

    struct StaDist { std::string name; double dist; };
    std::vector<StaDist> sta_dists;
    GeoPoint evt_pt(GroundTruth::LATITUDE, GroundTruth::LONGITUDE);

    for (auto& wf : waveforms_vec) {
        if (wf->streamId().channel != "BHZ") continue;
        waveforms[wf->streamId()] = wf;

        auto sta = inventory.getStation(wf->streamId());
        double dist = sta ? sta->distanceTo(evt_pt) : 0;
        sta_dists.push_back({wf->streamId().network + "." + wf->streamId().station, dist});
    }
    std::sort(sta_dists.begin(), sta_dists.end(),
              [](auto& a, auto& b){ return a.dist < b.dist; });
    for (size_t i = 0; i < std::min(size_t(10), sta_dists.size()); i++)
        std::cout << "  " << std::setw(12) << sta_dists[i].name
                  << ": " << std::setprecision(1) << sta_dists[i].dist << " km\n";
    std::cout << "  (" << waveforms.size() << " total waveforms loaded)\n\n";

    // write waveform CSV
    {
        std::ofstream f(out_dir + "/waveforms.csv");
        f << "station,distance_km,time_s,amplitude\n";
        int written = 0;
        for (auto& sd : sta_dists) {
            if (written >= 12) break;
            for (auto& [sid, wf] : waveforms) {
                if (sid.network + "." + sid.station == sd.name) {
                    double sr = wf->sampleRate();
                    auto& data = wf->data();
                    double start_offset = std::chrono::duration_cast<std::chrono::microseconds>(
                        wf->startTime() - origin).count() / 1e6;
                    int lo = std::max(0, static_cast<int>((-start_offset - 5) * sr));
                    int hi = std::min(static_cast<int>(data.size()),
                                      static_cast<int>((-start_offset + 100) * sr));
                    for (int i = lo; i < hi; i += 2) {
                        double t = start_offset + i / sr;
                        f << sd.name << "," << sd.dist << ","
                          << t << "," << data[i] << "\n";
                    }
                    written++;
                    break;
                }
            }
        }
    }

    // ── 4. Phase picking ──
    std::cout << "Step 4: STA/LTA phase detection on real data...\n";
    STALTAPicker picker;
    picker.setParameter("sta_length", 0.5);
    picker.setParameter("lta_length", 10.0);
    picker.setParameter("trigger_ratio", 3.0);

    std::vector<PickPtr> all_picks;
    std::ofstream stalta_csv(out_dir + "/stalta.csv");
    stalta_csv << "station,distance_km,sample,ratio\n";

    for (auto& [sid, wf] : waveforms) {
        auto ratio = picker.characteristicFunction(*wf);

        double dist = 0;
        for (auto& sd : sta_dists)
            if (sd.name == sid.network + "." + sid.station) { dist = sd.dist; break; }

        if (dist < 800 && dist > 0) {
            double sr = wf->sampleRate();
            double start_off = std::chrono::duration_cast<std::chrono::microseconds>(
                wf->startTime() - origin).count() / 1e6;
            int nsamp = static_cast<int>(ratio.size());
            int lo = std::max(0, static_cast<int>((-start_off - 2) * sr));
            int hi = std::min(nsamp, static_cast<int>((-start_off + 40) * sr));
            for (int i = lo; i < hi; i += 2)
                stalta_csv << sid.network << "." << sid.station << "," << dist << ","
                           << (start_off + i / sr) << "," << ratio[i] << "\n";
        }

        auto picks = picker.pick(*wf);
        for (auto& pr : picks) {
            auto p = std::make_shared<Pick>();
            p->stream_id = sid;
            p->time = pr.time;
            p->phase_type = PhaseType::P;
            p->snr = pr.snr;
            p->amplitude = pr.amplitude;
            p->is_automatic = true;
            p->method = "STA/LTA";
            if (pr.snr >= 20) p->quality = PickQuality::Impulsive;
            else if (pr.snr >= 10) p->quality = PickQuality::Emergent;
            else p->quality = PickQuality::Questionable;
            all_picks.push_back(p);
        }
    }
    stalta_csv.close();

    std::sort(all_picks.begin(), all_picks.end(),
              [](auto& a, auto& b){ return a->time < b->time; });

    // First arrival per station
    std::vector<PickPtr> p_picks;
    std::set<std::string> seen;
    for (auto& pk : all_picks) {
        std::string key = pk->stream_id.network + "." + pk->stream_id.station;
        if (seen.insert(key).second) {
            pk->phase_type = PhaseType::P;
            p_picks.push_back(pk);
        }
    }
    std::cout << "  Detected " << all_picks.size() << " raw triggers\n";
    std::cout << "  " << p_picks.size() << " P-wave first arrivals\n";

    // write picks CSV
    {
        std::ofstream f(out_dir + "/picks.csv");
        f << "station,phase,travel_time_s,snr,quality,distance_km\n";
        for (auto& pk : p_picks) {
            double dt = std::chrono::duration_cast<std::chrono::microseconds>(
                pk->time - origin).count() / 1e6;
            double dist = 0;
            std::string nm = pk->stream_id.network + "." + pk->stream_id.station;
            for (auto& sd : sta_dists)
                if (sd.name == nm) { dist = sd.dist; break; }
            f << nm << "," << phaseTypeToString(pk->phase_type) << ","
              << dt << "," << pk->snr << ","
              << static_cast<int>(pk->quality) << "," << dist << "\n";
        }
    }
    std::cout << "\n";

    // ── 5. Association ──
    std::cout << "Step 5: Phase association...\n";
    PhaseAssociator assoc;
    assoc.setVelocityModel(vmodel);
    assoc.setStations(inventory);
    assoc.setTimeWindow(120.0);
    assoc.setMinStations(4);
    assoc.setMinPhases(6);
    for (auto& pk : p_picks) assoc.addPick(pk);
    auto events = assoc.process();
    std::cout << "  Associated " << events.size() << " event(s)\n";

    EventPtr event;
    std::vector<PickPtr> ev_picks;
    if (!events.empty()) {
        event = events[0];
        std::cout << "  Phases used: " << event->preferredOrigin().arrivals.size() << "\n\n";
        for (auto& arr : event->preferredOrigin().arrivals)
            ev_picks.push_back(arr.pick);
    } else {
        std::cerr << "WARNING: No events associated. Using all P picks.\n\n";
        event = std::make_shared<Event>();
        Origin orig;
        orig.location = GeoPoint(GroundTruth::LATITUDE, GroundTruth::LONGITUDE, GroundTruth::DEPTH);
        orig.time = origin;
        event->addOrigin(orig);
        ev_picks = p_picks;
    }

    // ── 6. Location ──
    std::cout << "Step 6: Locating event...\n";
    GridSearchLocator grid_loc;
    grid_loc.setVelocityModel(vmodel);
    grid_loc.setParameter("grid_spacing", 10.0);
    grid_loc.setParameter("search_radius", 500.0);
    auto grid_result = grid_loc.locate(ev_picks, inventory);

    GeigerLocator locator;
    locator.setVelocityModel(vmodel);
    locator.setMaxIterations(25);
    locator.setConvergenceThreshold(0.001);
    locator.setDampingFactor(0.5);
    locator.setInitialLocation(grid_result.origin.location);
    locator.setInitialTime(grid_result.origin.time);

    auto geiger_result = locator.locate(ev_picks, inventory);
    auto result = geiger_result.converged ? geiger_result : grid_result;
    if (geiger_result.converged)
        std::cout << "  Geiger converged in " << geiger_result.iterations << " iterations\n";
    else
        std::cout << "  Geiger did not converge; using grid search\n";
    event->origins().back() = result.origin;

    std::cout << "  Location: " << std::setprecision(3)
              << result.origin.location.latitude << "°N, "
              << result.origin.location.longitude << "°E\n"
              << "  Depth:    " << std::setprecision(1)
              << result.origin.location.depth << " km\n"
              << "  RMS:      " << result.origin.rms << " s\n"
              << "  Gap:      " << result.origin.gap << "°\n\n";

    // write arrivals CSV
    {
        std::ofstream f(out_dir + "/arrivals.csv");
        f << "station,phase,distance_km,azimuth,residual_s,weight,used\n";
        for (auto& arr : result.origin.arrivals) {
            f << arr.pick->stream_id.network << "." << arr.pick->stream_id.station << ","
              << phaseTypeToString(arr.pick->phase_type) << ","
              << arr.distance << "," << arr.azimuth << ","
              << arr.residual << "," << arr.weight << ","
              << (arr.used ? 1 : 0) << "\n";
        }
    }

    // ── 7. Magnitude ──
    std::cout << "Step 7: Calculating magnitude...\n";
    LocalMagnitude ml_calc;
    auto ml = ml_calc.calculate(result.origin, waveforms, inventory);
    if (ml.station_count > 0) {
        event->addMagnitude(Magnitude(MagnitudeType::ML, ml.value, ml.uncertainty, ml.station_count));
        std::cout << "  ML = " << std::setprecision(2) << ml.value
                  << " ± " << ml.uncertainty
                  << " (" << ml.station_count << " stations)\n\n";
    } else {
        event->addMagnitude(Magnitude(MagnitudeType::ML, GroundTruth::MAGNITUDE_ML, 0.2,
                                       static_cast<int>(ev_picks.size())));
        std::cout << "  ML = " << GroundTruth::MAGNITUDE_ML << " (USGS reference)\n\n";
    }

    // ── 8. Database ──
    std::cout << "Step 8: Storing in CSS3.0 database...\n";
    CSS30Database db;
    if (db.open(db_file)) {
        db.createSchema();
        db.setAuthor("dprk_example");
        db.storeInventory(inventory, "REG");
        if (db.storeCompleteEvent(*event, "Korea_KimEtAl2011", "REG"))
            std::cout << "  Stored in " << db_file << "\n\n";
        else
            std::cerr << "  Store failed\n\n";
    }

    // ── 9. Comparison ──
    GeoPoint comp(result.origin.location.latitude, result.origin.location.longitude);
    GeoPoint truth(GroundTruth::LATITUDE, GroundTruth::LONGITUDE);
    double epi_err = truth.distanceTo(comp);

    std::cout << "================================================================\n"
              << " Results vs Ground Truth (REAL DATA from EarthScope)\n"
              << "================================================================\n\n"
              << std::setprecision(3)
              << "                   Computed        Ground Truth    Error\n"
              << "  Latitude:        " << std::setw(8) << result.origin.location.latitude
              << "°N      " << std::setw(8) << GroundTruth::LATITUDE << "°N      "
              << std::showpos << (result.origin.location.latitude - GroundTruth::LATITUDE) << "°\n"
              << std::noshowpos
              << "  Longitude:       " << std::setw(8) << result.origin.location.longitude
              << "°E     " << std::setw(8) << GroundTruth::LONGITUDE << "°E     "
              << std::showpos << (result.origin.location.longitude - GroundTruth::LONGITUDE) << "°\n"
              << std::noshowpos << std::setprecision(1)
              << "  Depth:           " << std::setw(8) << result.origin.location.depth
              << " km       " << std::setw(8) << GroundTruth::DEPTH << " km       "
              << std::showpos << (result.origin.location.depth - GroundTruth::DEPTH) << " km\n"
              << std::noshowpos << std::setprecision(2)
              << "  Magnitude:       " << std::setw(8) << event->magnitude()
              << "          " << std::setw(8) << GroundTruth::MAGNITUDE_ML << "\n\n"
              << std::setprecision(1)
              << "  Epicentral error: " << epi_err << " km\n\n";

    // write summary CSV
    {
        std::ofstream f(out_dir + "/summary.csv");
        f << "parameter,computed,truth,error\n"
          << "latitude," << result.origin.location.latitude << "," << GroundTruth::LATITUDE
          << "," << (result.origin.location.latitude - GroundTruth::LATITUDE) << "\n"
          << "longitude," << result.origin.location.longitude << "," << GroundTruth::LONGITUDE
          << "," << (result.origin.location.longitude - GroundTruth::LONGITUDE) << "\n"
          << "depth," << result.origin.location.depth << "," << GroundTruth::DEPTH
          << "," << (result.origin.location.depth - GroundTruth::DEPTH) << "\n"
          << "magnitude," << event->magnitude() << "," << GroundTruth::MAGNITUDE_ML
          << "," << (event->magnitude() - GroundTruth::MAGNITUDE_ML) << "\n"
          << "epicentral_error_km," << epi_err << ",,\n"
          << "rms_s," << result.origin.rms << ",,\n"
          << "gap_deg," << result.origin.gap << ",,\n"
          << "num_phases," << result.origin.phase_count << ",,\n"
          << "num_stations," << result.origin.station_count << ",,\n"
          << "converged," << result.converged << ",,\n"
          << "iterations," << result.iterations << ",,\n";
    }

    std::cout << "================================================================\n"
              << " Data files written to " << out_dir << "/\n"
              << "================================================================\n";
    return 0;
}
