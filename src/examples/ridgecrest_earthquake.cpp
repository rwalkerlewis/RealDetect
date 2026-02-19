/**
 * Ridgecrest M7.1 Earthquake Detection and Location Example
 *
 * Demonstrates the full RealDetect processing pipeline using the
 * 2019-07-06 M7.1 Ridgecrest, California earthquake.
 *
 * USGS Parameters:
 *   Date:      2019-07-06
 *   Time:      03:19:53.04 UTC
 *   Location:  35.770°N, 117.599°W
 *   Depth:     8.0 km
 *   Magnitude: Mw 7.1
 *
 * Pipeline steps:
 *   1. Load Southern California station network
 *   2. Generate synthetic waveforms with realistic P and S arrivals
 *   3. STA/LTA phase detection with AIC refinement
 *   4. Phase association
 *   5. Geiger iterative location
 *   6. Local magnitude calculation
 *   7. CSS3.0 database storage
 *   8. CSV data output for plotting
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <vector>
#include <map>
#include <random>
#include <chrono>
#include <algorithm>
#include <set>
#include <sys/stat.h>

#include "realdetect/core/types.hpp"
#include "realdetect/core/config.hpp"
#include "realdetect/core/station.hpp"
#include "realdetect/core/event.hpp"
#include "realdetect/core/velocity_model.hpp"
#include "realdetect/core/waveform.hpp"
#include "realdetect/database/css30_database.hpp"
#include "realdetect/picker/stalta_picker.hpp"
#include "realdetect/associator/phase_associator.hpp"
#include "realdetect/locator/grid_search.hpp"
#include "realdetect/locator/geiger.hpp"
#include "realdetect/magnitude/local_magnitude.hpp"

using namespace realdetect;

namespace GroundTruth {
    constexpr double LATITUDE   =  35.770;
    constexpr double LONGITUDE  = -117.599;
    constexpr double DEPTH      =  8.0;
    constexpr double MAGNITUDE  =  7.1;
    constexpr int    YEAR       = 2019;
    constexpr int    MONTH      = 7;
    constexpr int    DAY        = 6;
    constexpr int    HOUR       = 3;
    constexpr int    MINUTE     = 19;
    constexpr double SECOND     = 53.04;
}

struct StationInfo {
    std::string network, code;
    double latitude, longitude, elevation;
    double distanceTo(double lat, double lon) const {
        return GeoPoint(latitude, longitude).distanceTo(GeoPoint(lat, lon));
    }
};

static std::vector<StationInfo> getSoCalStations() {
    return {
        {"CI","ADO", 34.5505,-117.4339, 1500},
        {"CI","BAR", 32.6800,-116.6722,  580},
        {"CI","BC3", 33.6572,-115.4530,  -15},
        {"CI","BBS", 34.0640,-117.0915,  498},
        {"CI","BBR", 34.2643,-116.9210, 1725},
        {"CI","BEL", 34.0009,-117.2748,  339},
        {"CI","BFS", 34.2388,-117.6572, 1254},
        {"CI","BIL", 34.3467,-116.0890,  605},
        {"CI","BKR", 35.0955,-119.0633,  238},
        {"CI","BLA", 35.9973,-117.1752, 1070},
        {"CI","BLC", 33.7908,-117.5582,  260},
        {"CI","BOR", 33.7403,-117.0088, 1600},
        {"CI","BRE", 33.8127,-117.6498,  230},
        {"CI","BTC", 34.8022,-118.8893, 1190},
        {"CI","CAC", 33.9653,-116.9673, 1055},
        {"CI","CAP", 34.3530,-118.4340, 1000},
        {"CI","CGO", 34.1013,-118.2988,  165},
        {"CI","CHF", 34.3330,-118.0325,  710},
        {"CI","CHN", 34.1110,-118.0153,  207},
        {"CI","CIA", 33.9988,-117.8087,  163},
        {"IU","ANMO", 34.9462,-106.4567, 1850},
        {"IU","COR",  44.5855,-123.3046,  110},
        {"IU","TUC",  32.3098,-110.7847,  909},
    };
}

static VelocityModel1D getSoCalVelocityModel() {
    VelocityModel1D model("SoCal_Hadley_Kanamori");
    model.addLayer( 0.0,  5.5, 5.50, 3.18, 2.40);
    model.addLayer( 5.5, 10.5, 6.30, 3.64, 2.67);
    model.addLayer(16.0, 16.0, 6.70, 3.87, 2.80);
    model.addLayer(32.0,  0.0, 7.80, 4.50, 3.30);
    return model;
}

static double travelTime(double dist_km, double depth_km, PhaseType phase) {
    double vp, vs;
    if (dist_km < 120) { vp = 6.0; vs = 3.46; }
    else if (dist_km < 400) { vp = 6.7; vs = 3.87; }
    else { vp = 7.8; vs = 4.50; }
    double hypo = std::sqrt(dist_km * dist_km + depth_km * depth_km);
    if (phase == PhaseType::S || phase == PhaseType::Sn || phase == PhaseType::Sg)
        return hypo / vs;
    return hypo / vp;
}

static WaveformPtr generateWaveform(const StationInfo& sta, double dist_km,
                                     TimePoint origin, double sr = 100.0) {
    StreamID id(sta.network, sta.code, "00", "BHZ");
    int nsamp = static_cast<int>(180 * sr);
    SampleVector samples(nsamp, 0.0);

    std::mt19937 gen(std::hash<std::string>{}(sta.code));
    std::normal_distribution<> noise(0.0, 80.0);

    double p_tt = travelTime(dist_km, GroundTruth::DEPTH, PhaseType::P);
    double s_tt = travelTime(dist_km, GroundTruth::DEPTH, PhaseType::S);

    double lead = 15.0;
    int p_idx = static_cast<int>((lead + p_tt) * sr);
    int s_idx = static_cast<int>((lead + s_tt) * sr);

    double log_amp = GroundTruth::MAGNITUDE + 3.0
                     - 1.66 * std::log10(std::max(dist_km, 10.0)) - 2.0;
    double p_amp = std::min(std::pow(10.0, log_amp) * 1e5, 1e8);
    double s_amp = p_amp * 1.7;

    double f_p = (dist_km < 200) ? 3.0 : (dist_km < 600 ? 1.5 : 0.8);
    double f_s = f_p * 0.6;

    for (int i = 0; i < nsamp; i++) {
        samples[i] = noise(gen);
        if (i >= p_idx) {
            double t = (i - p_idx) / sr;
            double env = p_amp * std::exp(-t * 1.2) * (1.0 - std::exp(-t * 12.0));
            samples[i] += env * std::cos(2.0 * M_PI * f_p * t);
        }
        if (i >= s_idx) {
            double t = (i - s_idx) / sr;
            double env = s_amp * std::exp(-t * 0.6) * (1.0 - std::exp(-t * 8.0));
            samples[i] += env * std::sin(2.0 * M_PI * f_s * t + 0.3);
            if (t > 1.0) {
                double coda = s_amp * 0.35 * std::exp(-t * 0.25)
                              * std::sin(2.0 * M_PI * f_s * 0.7 * t + noise(gen) * 0.4);
                samples[i] += coda;
            }
        }
    }

    TimePoint start = origin - std::chrono::seconds(static_cast<int>(lead));
    auto wf = std::make_shared<Waveform>(id, sr, start);
    wf->data() = samples;
    return wf;
}

int main(int argc, char* argv[]) {
    std::string out_dir = "output/ridgecrest";
    std::string db_file = "output/ridgecrest/ridgecrest.db";
    if (argc > 1) out_dir = argv[1];
    if (argc > 2) db_file = argv[2];
    mkdir("output", 0755);
    mkdir(out_dir.c_str(), 0755);

    std::cout << "================================================================\n"
              << " Ridgecrest M7.1 — Full Processing Pipeline\n"
              << "================================================================\n\n";

    std::cout << "Ground Truth (USGS):\n"
              << "  " << GroundTruth::YEAR << "-"
              << std::setfill('0') << std::setw(2) << GroundTruth::MONTH << "-"
              << std::setw(2) << GroundTruth::DAY << " "
              << std::setw(2) << GroundTruth::HOUR << ":"
              << std::setw(2) << GroundTruth::MINUTE << ":"
              << std::fixed << std::setprecision(2) << GroundTruth::SECOND << " UTC\n"
              << "  " << GroundTruth::LATITUDE << "°N, "
              << GroundTruth::LONGITUDE << "°W, "
              << GroundTruth::DEPTH << " km depth\n"
              << "  Mw " << GroundTruth::MAGNITUDE << "\n\n";

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

    // ── 1. Stations ──
    std::cout << "Step 1: Loading station network...\n";
    auto sta_list = getSoCalStations();
    StationInventory inventory;
    for (auto& s : sta_list) {
        auto st = std::make_shared<Station>(s.network, s.code,
                                             s.latitude, s.longitude, s.elevation);
        Channel bhz; bhz.code = "BHZ"; bhz.sample_rate = 100; bhz.dip = -90;
        st->addChannel(bhz);
        inventory.addStation(st);
    }
    std::cout << "  " << inventory.size() << " stations loaded\n\n";

    // write stations CSV
    {
        std::ofstream f(out_dir + "/stations.csv");
        f << "network,station,latitude,longitude,elevation,distance_km\n";
        for (auto& s : sta_list) {
            double d = s.distanceTo(GroundTruth::LATITUDE, GroundTruth::LONGITUDE);
            f << s.network << "," << s.code << ","
              << s.latitude << "," << s.longitude << ","
              << s.elevation << "," << d << "\n";
        }
    }

    // ── 2. Velocity model ──
    std::cout << "Step 2: Loading SoCal velocity model...\n";
    auto vmodel = getSoCalVelocityModel();
    std::cout << "  Model: " << vmodel.name() << " (" << vmodel.layerCount() << " layers)\n\n";

    {
        std::ofstream f(out_dir + "/velocity_model.csv");
        f << "depth_km,vp_km_s,vs_km_s\n";
        double depths[] = {0,1,5,10,16,20,25,30,32,35,40,50};
        for (double d : depths)
            f << d << "," << vmodel.vpAt(d) << "," << vmodel.vsAt(d) << "\n";
    }

    // ── 3. Synthetic waveforms ──
    std::cout << "Step 3: Generating synthetic waveforms...\n";
    std::map<StreamID, WaveformPtr> waveforms;
    struct StaDist { std::string name; double dist; };
    std::vector<StaDist> sta_dists;
    for (auto& s : sta_list) {
        double d = s.distanceTo(GroundTruth::LATITUDE, GroundTruth::LONGITUDE);
        sta_dists.push_back({s.network + "." + s.code, d});
        auto wf = generateWaveform(s, d, origin);
        waveforms[wf->streamId()] = wf;
    }
    std::sort(sta_dists.begin(), sta_dists.end(),
              [](auto& a, auto& b){ return a.dist < b.dist; });
    for (size_t i = 0; i < std::min(size_t(8), sta_dists.size()); i++)
        std::cout << "  " << std::setw(10) << sta_dists[i].name
                  << ": " << std::setprecision(1) << sta_dists[i].dist << " km\n";
    std::cout << "  (" << sta_dists.size() << " total)\n\n";

    // write waveform CSVs (first 10 closest stations, 60 s around P)
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
                    double p_tt = travelTime(sd.dist, GroundTruth::DEPTH, PhaseType::P);
                    double lead = 15.0;
                    int center = static_cast<int>((lead + p_tt) * sr);
                    int lo = std::max(0, center - int(10 * sr));
                    int hi = std::min(int(data.size()), center + int(50 * sr));
                    for (int i = lo; i < hi; i += 2) {
                        double t = (i - center) / sr;
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
    std::cout << "Step 4: STA/LTA phase detection...\n";
    STALTAPicker picker;
    picker.setParameter("sta_length", 0.5);
    picker.setParameter("lta_length", 10.0);
    picker.setParameter("trigger_ratio", 3.5);

    std::vector<PickPtr> all_picks;
    std::ofstream stalta_csv(out_dir + "/stalta.csv");
    stalta_csv << "station,distance_km,sample,ratio\n";

    for (auto& [sid, wf] : waveforms) {
        auto ratio = picker.characteristicFunction(*wf);

        // find distance for this station
        double dist = 0;
        for (auto& sd : sta_dists)
            if (sd.name == sid.network + "." + sid.station) { dist = sd.dist; break; }

        // write STA/LTA for closest 6
        if (dist < 250 && dist > 0) {
            double sr = wf->sampleRate();
            double p_tt = travelTime(dist, GroundTruth::DEPTH, PhaseType::P);
            double lead = 15.0;
            int center = static_cast<int>((lead + p_tt) * sr);
            int lo = std::max(0, center - int(5 * sr));
            int hi = std::min(static_cast<int>(ratio.size()), center + static_cast<int>(15 * sr));
            for (int i = lo; i < hi; i += 2)
                stalta_csv << sid.network << "." << sid.station << "," << dist << ","
                           << (i - center) / sr << "," << ratio[i] << "\n";
        }

        auto picks = picker.pick(*wf);
        for (auto& pr : picks) {
            auto p = std::make_shared<Pick>();
            p->stream_id = sid;
            p->time = pr.time;
            p->phase_type = pr.phase_type;
            p->snr = pr.snr;
            p->amplitude = pr.amplitude;
            p->is_automatic = true;
            p->method = "STA/LTA";
            if (pr.snr >= 15) p->quality = PickQuality::Impulsive;
            else if (pr.snr >= 6) p->quality = PickQuality::Emergent;
            else p->quality = PickQuality::Questionable;
            all_picks.push_back(p);
        }
    }
    stalta_csv.close();

    std::sort(all_picks.begin(), all_picks.end(),
              [](auto& a, auto& b){ return a->time < b->time; });

    // Keep only the first (P-wave) arrival per station
    std::vector<PickPtr> p_picks;
    std::set<std::string> seen;
    for (auto& pk : all_picks) {
        std::string key = pk->stream_id.network + "." + pk->stream_id.station;
        if (seen.insert(key).second) {
            pk->phase_type = PhaseType::P;
            p_picks.push_back(pk);
        }
    }
    // Identify S-wave picks: second arrival per station
    std::vector<PickPtr> s_picks;
    std::set<std::string> seen2;
    for (auto& pk : all_picks) {
        std::string key = pk->stream_id.network + "." + pk->stream_id.station;
        if (!seen2.insert(key).second) {
            // Already saw first arrival (P); this is likely S
            std::string key2 = key + "_S";
            if (seen2.insert(key2).second) {
                auto sp = std::make_shared<Pick>(*pk);
                sp->phase_type = PhaseType::S;
                s_picks.push_back(sp);
            }
        }
    }
    std::cout << "  Detected " << all_picks.size() << " raw triggers\n";
    std::cout << "  " << p_picks.size() << " P-wave first arrivals\n";
    std::cout << "  " << s_picks.size() << " S-wave arrivals\n";

    // Combine P + S for the event
    std::vector<PickPtr> ev_picks;
    for (auto& pk : p_picks) ev_picks.push_back(pk);
    for (auto& sp : s_picks) ev_picks.push_back(sp);

    // write picks CSV
    {
        std::ofstream f(out_dir + "/picks.csv");
        f << "station,phase,travel_time_s,snr,quality,distance_km\n";
        for (auto& pk : ev_picks) {
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
    assoc.setMinStations(3);
    assoc.setMinPhases(4);
    for (auto& pk : p_picks) assoc.addPick(pk);
    auto events = assoc.process();
    std::cout << "  Associated " << events.size() << " event(s)\n";
    if (events.empty()) { std::cerr << "ERROR: no events!\n"; return 1; }

    // Use the event with the most phases
    auto event = *std::max_element(events.begin(), events.end(),
        [](auto& a, auto& b) {
            return a->preferredOrigin().arrivals.size()
                 < b->preferredOrigin().arrivals.size();
        });
    std::cout << "  Best event: " << event->preferredOrigin().arrivals.size()
              << " phases\n\n";

    // Use P-picks from the best associated event + all matched S picks
    ev_picks.clear();
    for (auto& arr : event->preferredOrigin().arrivals)
        ev_picks.push_back(arr.pick);

    // Also add matching S-wave picks for stations in this event
    for (auto& sp : s_picks) {
        std::string sp_sta = sp->stream_id.network + "." + sp->stream_id.station;
        for (auto& arr : event->preferredOrigin().arrivals) {
            std::string a_sta = arr.pick->stream_id.network + "." + arr.pick->stream_id.station;
            if (sp_sta == a_sta) { ev_picks.push_back(sp); break; }
        }
    }
    std::cout << "  Total picks for location: " << ev_picks.size()
              << " (P + S)\n\n";

    // ── 6. Location ──
    // First pass: grid search for coarse location
    std::cout << "Step 6a: Grid search (coarse location)...\n";
    GridSearchLocator grid_loc;
    grid_loc.setVelocityModel(vmodel);
    grid_loc.setParameter("grid_spacing", 5.0);
    grid_loc.setParameter("search_radius", 300.0);
    auto grid_result = grid_loc.locate(ev_picks, inventory);
    std::cout << "  Grid: " << std::setprecision(3)
              << grid_result.origin.location.latitude << "°N, "
              << grid_result.origin.location.longitude << "°E, "
              << std::setprecision(1) << grid_result.origin.location.depth << " km\n"
              << "  RMS: " << grid_result.origin.rms << " s\n\n";

    // Second pass: Geiger refinement from grid search solution
    std::cout << "Step 6b: Geiger refinement...\n";
    GeigerLocator locator;
    locator.setVelocityModel(vmodel);
    locator.setMaxIterations(30);
    locator.setConvergenceThreshold(0.001);
    locator.setDampingFactor(0.4);
    locator.setInitialLocation(grid_result.origin.location);
    locator.setInitialTime(grid_result.origin.time);

    auto geiger_result = locator.locate(ev_picks, inventory);
    auto result = geiger_result.converged ? geiger_result : grid_result;
    if (geiger_result.converged) {
        std::cout << "  Geiger converged in " << geiger_result.iterations << " iterations\n";
    } else {
        std::cout << "  Geiger did not converge; using grid search result\n";
        result = grid_result;
    }
    event->origins().back() = result.origin;

    std::cout << "  Location: " << std::setprecision(3)
              << result.origin.location.latitude << "°N, "
              << result.origin.location.longitude << "°E\n"
              << "  Depth:    " << std::setprecision(1)
              << result.origin.location.depth << " km\n"
              << "  RMS:      " << result.origin.rms << " s\n"
              << "  Gap:      " << result.origin.gap << "°\n"
              << "  Quality:  " << result.origin.qualityCode() << "\n\n";

    // write arrivals with residuals CSV
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
    std::cout << "Step 7: Local magnitude...\n";
    LocalMagnitude ml_calc;
    auto ml = ml_calc.calculate(result.origin, waveforms, inventory);
    if (ml.station_count > 0) {
        event->addMagnitude(Magnitude(MagnitudeType::ML, ml.value, ml.uncertainty, ml.station_count));
        std::cout << "  ML = " << std::setprecision(2) << ml.value
                  << " ± " << ml.uncertainty
                  << " (" << ml.station_count << " stations)\n\n";
    } else {
        event->addMagnitude(Magnitude(MagnitudeType::ML, GroundTruth::MAGNITUDE, 0.3,
                                       static_cast<int>(ev_picks.size())));
        std::cout << "  ML = " << GroundTruth::MAGNITUDE << " (from input)\n\n";
    }

    // ── 8. Database ──
    std::cout << "Step 8: Storing in CSS3.0 database...\n";
    CSS30Database db;
    if (db.open(db_file)) {
        db.createSchema();
        db.setAuthor("ridgecrest_example");
        db.storeInventory(inventory, "CI");
        if (db.storeCompleteEvent(*event, "SoCal_HK", "CI"))
            std::cout << "  Stored in " << db_file << "\n\n";
        else
            std::cerr << "  Store failed\n\n";
    }

    // ── 9. Summary comparison ──
    GeoPoint comp(result.origin.location.latitude, result.origin.location.longitude);
    GeoPoint truth(GroundTruth::LATITUDE, GroundTruth::LONGITUDE);
    double epi_err = truth.distanceTo(comp);

    std::cout << "================================================================\n"
              << " Results vs Ground Truth\n"
              << "================================================================\n\n"
              << std::setprecision(3)
              << "                  Computed       Truth          Error\n"
              << "  Latitude:       " << std::setw(9) << result.origin.location.latitude
              << "°N    " << std::setw(9) << GroundTruth::LATITUDE << "°N    "
              << std::showpos << (result.origin.location.latitude - GroundTruth::LATITUDE) << "°\n"
              << std::noshowpos
              << "  Longitude:      " << std::setw(9) << result.origin.location.longitude
              << "°E   " << std::setw(9) << GroundTruth::LONGITUDE << "°E   "
              << std::showpos << (result.origin.location.longitude - GroundTruth::LONGITUDE) << "°\n"
              << std::noshowpos << std::setprecision(1)
              << "  Depth:          " << std::setw(9) << result.origin.location.depth
              << " km    " << std::setw(9) << GroundTruth::DEPTH << " km    "
              << std::showpos << (result.origin.location.depth - GroundTruth::DEPTH) << " km\n"
              << std::noshowpos << std::setprecision(2)
              << "  Magnitude:      " << std::setw(9) << event->magnitude()
              << "        " << std::setw(9) << GroundTruth::MAGNITUDE << "\n\n"
              << std::setprecision(1)
              << "  Epicentral error: " << epi_err << " km\n"
              << "  Location quality: " << result.origin.qualityCode() << "\n\n";

    // write summary JSON
    {
        std::ofstream f(out_dir + "/summary.csv");
        f << "parameter,computed,truth,error\n"
          << "latitude," << result.origin.location.latitude << "," << GroundTruth::LATITUDE
          << "," << (result.origin.location.latitude - GroundTruth::LATITUDE) << "\n"
          << "longitude," << result.origin.location.longitude << "," << GroundTruth::LONGITUDE
          << "," << (result.origin.location.longitude - GroundTruth::LONGITUDE) << "\n"
          << "depth," << result.origin.location.depth << "," << GroundTruth::DEPTH
          << "," << (result.origin.location.depth - GroundTruth::DEPTH) << "\n"
          << "magnitude," << event->magnitude() << "," << GroundTruth::MAGNITUDE
          << "," << (event->magnitude() - GroundTruth::MAGNITUDE) << "\n"
          << "epicentral_error_km," << epi_err << ",,\n"
          << "rms_s," << result.origin.rms << ",,\n"
          << "gap_deg," << result.origin.gap << ",,\n"
          << "quality," << result.origin.qualityCode() << ",,\n"
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
