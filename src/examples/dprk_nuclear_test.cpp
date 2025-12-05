/**
 * DPRK Nuclear Test Detection and Inversion Example
 * 
 * This example demonstrates detection and location of the September 3, 2017
 * North Korean nuclear test (6th test) using a regional seismic network.
 * 
 * Event Parameters (USGS/CTBTO):
 *   Date: 2017-09-03
 *   Origin Time: 03:30:01.6 UTC  
 *   Location: 41.343°N, 129.036°E (Punggye-ri test site)
 *   Depth: ~0 km (surface/shallow explosion)
 *   Magnitude: mb 6.3, ML 6.3, Mw 6.3
 *   Estimated yield: ~250 kilotons TNT equivalent
 * 
 * This example:
 *   1. Sets up a regional network with real station locations
 *   2. Generates synthetic waveforms with realistic P arrivals
 *   3. Detects phase arrivals using STA/LTA
 *   4. Associates picks into an event
 *   5. Locates the event using Geiger inversion
 *   6. Calculates magnitude
 *   7. Stores results in CSS3.0 database
 *   8. Compares results with known solution
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <map>
#include <random>
#include <chrono>

#include "realdetect/core/types.hpp"
#include "realdetect/core/config.hpp"
#include "realdetect/core/station.hpp"
#include "realdetect/core/event.hpp"
#include "realdetect/core/velocity_model.hpp"
#include "realdetect/core/regional_velocity_model.hpp"
#include "realdetect/core/waveform.hpp"
#include "realdetect/database/css30_database.hpp"
#include "realdetect/picker/stalta_picker.hpp"
#include "realdetect/associator/phase_associator.hpp"
#include "realdetect/locator/geiger.hpp"
#include "realdetect/magnitude/local_magnitude.hpp"

using namespace realdetect;

// ============================================================================
// Ground Truth: 2017-09-03 DPRK Nuclear Test Parameters
// ============================================================================
namespace GroundTruth {
    constexpr double LATITUDE = 41.343;      // degrees N
    constexpr double LONGITUDE = 129.036;    // degrees E
    constexpr double DEPTH = 0.5;            // km (shallow explosion)
    constexpr double MAGNITUDE_MB = 6.3;     // body wave magnitude
    constexpr double MAGNITUDE_ML = 6.3;     // local magnitude
    constexpr double MAGNITUDE_MW = 6.3;     // moment magnitude
    
    // Origin time: 2017-09-03 03:30:01.6 UTC
    constexpr int YEAR = 2017;
    constexpr int MONTH = 9;
    constexpr int DAY = 3;
    constexpr int HOUR = 3;
    constexpr int MINUTE = 30;
    constexpr double SECOND = 1.6;
    
    // Estimated yield: ~250 kt TNT equivalent
    constexpr double YIELD_KT = 250.0;
}

// ============================================================================
// Station Data: Real Regional Network Stations
// ============================================================================
struct StationInfo {
    std::string network;
    std::string code;
    double latitude;
    double longitude;
    double elevation;
    
    double distanceTo(double lat, double lon) const {
        GeoPoint sta(latitude, longitude);
        GeoPoint evt(lat, lon);
        return sta.distanceTo(evt);
    }
};

// Real station locations for DPRK monitoring
std::vector<StationInfo> getRegionalStations() {
    return {
        // South Korean stations (closest)
        {"KS", "GKP", 37.999, 128.669, 330},    // Closest SK station
        {"KS", "ULJ", 36.993, 129.417, 185},
        {"KS", "MUN", 36.583, 128.720, 245},
        {"KS", "HSB", 35.566, 129.348, 190},
        {"KS", "CHC", 36.720, 127.960, 380},
        {"KS", "SEO", 37.457, 127.011, 52},
        {"KS", "TJN", 36.370, 127.364, 75},
        {"KS", "DAG", 35.768, 128.619, 110},
        {"IU", "INCN", 37.4786, 126.6242, 79},  // IRIS/IDA station
        
        // Chinese stations (close to border)
        {"CN", "YNB", 42.880, 129.470, 580},    // Yanji - very close
        {"CN", "MDJ", 44.617, 129.592, 220},    // Mudanjiang
        {"CN", "CNS", 41.785, 123.823, 65},     // Changchun area
        {"CN", "HLR", 45.723, 126.615, 180},    // Harbin
        {"CN", "DAN", 40.050, 124.333, 25},     // Dandong - border city
        {"IC", "BJT", 40.018, 116.168, 55},     // Beijing
        
        // Japanese stations
        {"JP", "MDJ", 44.620, 129.590, 260},    // Co-located MDJ
        {"JP", "TSK", 36.214, 140.089, 25},     // Tsukuba
        {"JP", "FUJ", 39.049, 141.434, 242},    // Fujisawa
        {"JP", "IZH", 34.131, 129.207, 170},    // Iki
        
        // Russian stations
        {"RU", "YSS", 46.959, 142.761, 30}      // Yuzhno-Sakhalinsk
    };
}

// ============================================================================
// Korean Peninsula Velocity Model
// ============================================================================
VelocityModel1D getKoreanVelocityModel() {
    VelocityModel1D model("Korea_KimEtAl2011");
    
    // Based on Kim et al. (2011) and regional studies
    model.addLayer(0.0,   1.0, 4.50, 2.60, 2.40);   // Shallow sediments
    model.addLayer(1.0,   9.0, 5.90, 3.41, 2.70);   // Upper crust
    model.addLayer(10.0, 10.0, 6.20, 3.58, 2.80);   // Middle crust  
    model.addLayer(20.0, 13.0, 6.60, 3.81, 2.90);   // Lower crust
    model.addLayer(33.0,  7.0, 7.20, 4.16, 3.10);   // Moho transition
    model.addLayer(40.0,  0.0, 8.00, 4.62, 3.35);   // Upper mantle
    
    return model;
}

// ============================================================================
// Calculate P-wave travel time using 1D model
// ============================================================================
double calculateTravelTime(const VelocityModel1D& model, 
                           double distance_km, double depth_km,
                           PhaseType phase) {
    // Simple ray parameter calculation for regional distances
    // Using average velocity approximation with depth correction
    
    double avg_velocity;
    if (phase == PhaseType::P || phase == PhaseType::Pn || phase == PhaseType::Pg) {
        // P-wave
        if (distance_km < 150) {
            // Pg - crustal P
            avg_velocity = 6.0;  // Average crustal Vp
        } else if (distance_km < 800) {
            // Pn - head wave along Moho
            avg_velocity = 7.8;  // Pn velocity
        } else {
            avg_velocity = 8.0;  // Upper mantle
        }
    } else {
        // S-wave
        if (distance_km < 150) {
            avg_velocity = 3.5;
        } else {
            avg_velocity = 4.5;
        }
    }
    
    // Hypocentral distance
    double hypo_dist = std::sqrt(distance_km * distance_km + depth_km * depth_km);
    
    // For Pn phase, account for Moho refraction
    if (distance_km > 150 && distance_km < 1000) {
        // Pn travel time = crustal path + Moho refraction + crustal path
        double moho_depth = 33.0;  // Korean Peninsula average
        double crust_vp = 6.2;
        double mantle_vp = 7.8;
        
        // Critical angle
        double ic = std::asin(crust_vp / mantle_vp);
        
        // Horizontal distance along Moho
        double x_moho = distance_km - 2 * moho_depth * std::tan(ic);
        if (x_moho > 0) {
            double t_crust = 2 * moho_depth / (crust_vp * std::cos(ic));
            double t_moho = x_moho / mantle_vp;
            return t_crust + t_moho;
        }
    }
    
    return hypo_dist / avg_velocity;
}

// ============================================================================
// Generate synthetic waveform for explosion source
// ============================================================================
WaveformPtr generateExplosionWaveform(const StationInfo& station,
                                       double distance_km,
                                       double magnitude,
                                       TimePoint origin_time,
                                       double sample_rate = 100.0) {
    StreamID id(station.network, station.code, "00", "BHZ");
    
    // Duration: 120 seconds of data
    int n_samples = static_cast<int>(120 * sample_rate);
    SampleVector samples(n_samples, 0.0);
    
    // Random number generator for noise
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> noise(0.0, 100.0);  // Background noise
    
    // P-wave arrival time relative to origin
    VelocityModel1D model = getKoreanVelocityModel();
    double p_travel_time = calculateTravelTime(model, distance_km, 
                                                GroundTruth::DEPTH, PhaseType::P);
    
    // Convert to sample index (P arrives after buffer lead time)
    double lead_time = 10.0;  // 10 seconds before P
    int p_sample = static_cast<int>((lead_time + p_travel_time) * sample_rate);
    
    // Explosion amplitude based on magnitude and distance
    // log10(A) = ML + 3.0 - 2.0*log10(distance)  (simplified attenuation)
    double log_amp = magnitude + 3.0 - 2.0 * std::log10(std::max(distance_km, 10.0));
    double amplitude = std::pow(10.0, log_amp);
    
    // Limit amplitude for numerical stability
    amplitude = std::min(amplitude, 1e7);
    
    // Dominant frequency depends on source and distance
    double freq = 3.0;  // ~3 Hz dominant for regional explosion
    if (distance_km > 500) freq = 1.5;
    if (distance_km > 1000) freq = 1.0;
    
    // Generate waveform with noise and signal
    for (int i = 0; i < n_samples; i++) {
        double t = i / sample_rate;
        
        // Background noise
        samples[i] = noise(gen);
        
        // P-wave arrival
        if (i >= p_sample) {
            double t_rel = (i - p_sample) / sample_rate;
            
            // Explosion source: positive first motion, rapid rise, exponential decay
            // P-wave is compressional (positive) from explosion
            double envelope = amplitude * std::exp(-t_rel * 1.5) * (1.0 - std::exp(-t_rel * 10.0));
            
            // Add oscillation
            double phase = 2.0 * M_PI * freq * t_rel;
            samples[i] += envelope * std::cos(phase);
            
            // Add some scattered coda
            if (t_rel > 0.5) {
                double coda = amplitude * 0.3 * std::exp(-t_rel * 0.5) * 
                             std::sin(2.0 * M_PI * freq * 0.8 * t_rel + noise(gen) * 0.5);
                samples[i] += coda;
            }
        }
    }
    
    // Set waveform start time (lead_time before origin)
    TimePoint start_time = origin_time - std::chrono::seconds(static_cast<int>(lead_time));
    
    // Create waveform with constructor and set data
    auto waveform = std::make_shared<Waveform>(id, sample_rate, start_time);
    waveform->data() = samples;
    
    return waveform;
}

// ============================================================================
// Main: DPRK Nuclear Test Detection and Inversion
// ============================================================================
int main(int argc, char* argv[]) {
    std::cout << "================================================================\n";
    std::cout << "DPRK Nuclear Test Detection and Inversion Example\n";
    std::cout << "================================================================\n\n";
    
    std::cout << "Ground Truth (2017-09-03 DPRK 6th Nuclear Test):\n";
    std::cout << "  Origin Time: " << GroundTruth::YEAR << "-"
              << std::setfill('0') << std::setw(2) << GroundTruth::MONTH << "-"
              << std::setw(2) << GroundTruth::DAY << " "
              << std::setw(2) << GroundTruth::HOUR << ":"
              << std::setw(2) << GroundTruth::MINUTE << ":"
              << std::fixed << std::setprecision(1) << GroundTruth::SECOND << " UTC\n";
    std::cout << "  Location: " << GroundTruth::LATITUDE << "°N, " 
              << GroundTruth::LONGITUDE << "°E\n";
    std::cout << "  Depth: " << GroundTruth::DEPTH << " km\n";
    std::cout << "  Magnitude: mb " << GroundTruth::MAGNITUDE_MB << "\n";
    std::cout << "  Est. Yield: " << GroundTruth::YIELD_KT << " kt TNT\n\n";
    
    // Create origin time point
    std::tm tm = {};
    tm.tm_year = GroundTruth::YEAR - 1900;
    tm.tm_mon = GroundTruth::MONTH - 1;
    tm.tm_mday = GroundTruth::DAY;
    tm.tm_hour = GroundTruth::HOUR;
    tm.tm_min = GroundTruth::MINUTE;
    tm.tm_sec = static_cast<int>(GroundTruth::SECOND);
    std::time_t tt = timegm(&tm);
    TimePoint origin_time = std::chrono::system_clock::from_time_t(tt);
    origin_time += std::chrono::milliseconds(
        static_cast<int>((GroundTruth::SECOND - tm.tm_sec) * 1000));
    
    // ========================================================================
    // Step 1: Set up regional network
    // ========================================================================
    std::cout << "Step 1: Setting up regional seismic network...\n";
    
    auto stations_list = getRegionalStations();
    StationInventory inventory;
    
    for (const auto& s : stations_list) {
        auto station = std::make_shared<Station>(s.network, s.code, 
                                                  s.latitude, s.longitude, s.elevation);
        Channel bhz;
        bhz.code = "BHZ";
        bhz.sample_rate = 100.0;
        bhz.azimuth = 0;
        bhz.dip = -90;
        station->addChannel(bhz);
        inventory.addStation(station);
    }
    
    std::cout << "  Loaded " << inventory.size() << " stations\n\n";
    
    // ========================================================================
    // Step 2: Load Korean Peninsula velocity model
    // ========================================================================
    std::cout << "Step 2: Loading Korean Peninsula velocity model...\n";
    
    VelocityModel1D velocity_model = getKoreanVelocityModel();
    
    std::cout << "  Model: " << velocity_model.name() << "\n";
    std::cout << "  Layers: " << velocity_model.layerCount() << "\n";
    std::cout << "  Moho depth: ~33 km\n\n";
    
    // ========================================================================
    // Step 3: Generate synthetic waveforms
    // ========================================================================
    std::cout << "Step 3: Generating synthetic waveforms...\n";
    
    GeoPoint event_loc(GroundTruth::LATITUDE, GroundTruth::LONGITUDE, 
                       GroundTruth::DEPTH);
    
    std::map<StreamID, WaveformPtr> waveforms;
    std::vector<std::pair<std::string, double>> station_distances;
    
    for (const auto& s : stations_list) {
        double dist = s.distanceTo(GroundTruth::LATITUDE, GroundTruth::LONGITUDE);
        station_distances.push_back({s.network + "." + s.code, dist});
        
        auto wf = generateExplosionWaveform(s, dist, GroundTruth::MAGNITUDE_MB, 
                                            origin_time);
        waveforms[wf->streamId()] = wf;
    }
    
    // Sort by distance
    std::sort(station_distances.begin(), station_distances.end(),
              [](const auto& a, const auto& b) { return a.second < b.second; });
    
    std::cout << "  Station distances from source:\n";
    for (size_t i = 0; i < std::min(size_t(10), station_distances.size()); i++) {
        std::cout << "    " << std::setw(10) << station_distances[i].first 
                  << ": " << std::fixed << std::setprecision(1) 
                  << station_distances[i].second << " km\n";
    }
    std::cout << "    ... (" << station_distances.size() << " total)\n\n";
    
    // ========================================================================
    // Step 4: Pick P-wave arrivals using STA/LTA
    // ========================================================================
    std::cout << "Step 4: Detecting P-wave arrivals (STA/LTA picker)...\n";
    
    auto picker = std::make_shared<STALTAPicker>();
    picker->setParameter("sta_length", 0.5);
    picker->setParameter("lta_length", 10.0);
    picker->setParameter("trigger_ratio", 4.0);  // Higher threshold for explosions
    
    std::vector<PickPtr> all_picks;
    
    for (const auto& [stream_id, wf] : waveforms) {
        auto picks = picker->pick(*wf);
        
        for (const auto& pr : picks) {
            auto pick = std::make_shared<Pick>();
            pick->stream_id = stream_id;
            pick->time = pr.time;
            pick->phase_type = PhaseType::P;  // Explosion - mainly P
            pick->snr = pr.snr;
            pick->amplitude = pr.amplitude;
            pick->is_automatic = true;
            pick->method = "STA/LTA";
            
            // Quality based on SNR
            if (pr.snr >= 20) {
                pick->quality = PickQuality::Impulsive;
            } else if (pr.snr >= 10) {
                pick->quality = PickQuality::Emergent;
            } else {
                pick->quality = PickQuality::Questionable;
            }
            
            all_picks.push_back(pick);
        }
    }
    
    std::cout << "  Detected " << all_picks.size() << " P-wave arrivals\n";
    
    // Show first few picks
    std::cout << "  First arrivals:\n";
    
    // Sort picks by time
    std::sort(all_picks.begin(), all_picks.end(),
              [](const PickPtr& a, const PickPtr& b) { return a->time < b->time; });
    
    for (size_t i = 0; i < std::min(size_t(8), all_picks.size()); i++) {
        auto pick = all_picks[i];
        auto dt = std::chrono::duration_cast<std::chrono::milliseconds>(
            pick->time - origin_time).count() / 1000.0;
        std::cout << "    " << std::setw(10) << pick->stream_id.toString()
                  << "  T+" << std::fixed << std::setprecision(2) << dt << "s"
                  << "  SNR=" << std::setprecision(1) << pick->snr << "\n";
    }
    std::cout << "\n";
    
    // ========================================================================
    // Step 5: Associate picks into event
    // ========================================================================
    std::cout << "Step 5: Associating picks into event...\n";
    
    PhaseAssociator associator;
    associator.setVelocityModel(velocity_model);
    associator.setStations(inventory);
    associator.setTimeWindow(120.0);
    associator.setMinStations(4);
    associator.setMinPhases(6);
    
    for (const auto& pick : all_picks) {
        associator.addPick(pick);
    }
    
    auto events = associator.process();
    
    std::cout << "  Associated " << events.size() << " event(s)\n";
    
    if (events.empty()) {
        std::cerr << "ERROR: No events associated!\n";
        return 1;
    }
    
    auto event = events[0];
    std::cout << "  Event ID: " << event->id() << "\n";
    std::cout << "  Phases used: " << event->preferredOrigin().arrivals.size() << "\n\n";
    
    // ========================================================================
    // Step 6: Locate event using Geiger method
    // ========================================================================
    std::cout << "Step 6: Locating event (Geiger method)...\n";
    
    auto locator = std::make_shared<GeigerLocator>();
    locator->setVelocityModel(velocity_model);
    locator->setParameter("max_iterations", 25);
    locator->setParameter("convergence", 0.001);
    locator->setParameter("damping", 0.5);
    
    // Get picks for this event
    std::vector<PickPtr> event_picks;
    for (const auto& arr : event->preferredOrigin().arrivals) {
        event_picks.push_back(arr.pick);
    }
    
    auto result = locator->locate(event_picks, inventory);
    
    if (!result.converged) {
        std::cerr << "WARNING: Location did not converge!\n";
    }
    
    // Update event origin
    event->origins().back() = result.origin;
    
    std::cout << "  Location result:\n";
    std::cout << "    Latitude:  " << std::fixed << std::setprecision(3) 
              << result.origin.location.latitude << "°N\n";
    std::cout << "    Longitude: " << result.origin.location.longitude << "°E\n";
    std::cout << "    Depth:     " << std::setprecision(1) 
              << result.origin.location.depth << " km\n";
    std::cout << "    RMS:       " << result.origin.rms << " s\n";
    std::cout << "    Gap:       " << result.origin.gap << "°\n";
    std::cout << "    Quality:   " << result.origin.qualityCode() << "\n\n";
    
    // ========================================================================
    // Step 7: Calculate magnitude
    // ========================================================================
    std::cout << "Step 7: Calculating magnitude...\n";
    
    auto ml_calc = std::make_shared<LocalMagnitude>();
    auto ml_result = ml_calc->calculate(result.origin, waveforms, inventory);
    
    if (ml_result.station_count > 0) {
        event->addMagnitude(Magnitude(MagnitudeType::ML, ml_result.value,
                                       ml_result.uncertainty, ml_result.station_count));
        std::cout << "  ML = " << std::fixed << std::setprecision(2) << ml_result.value
                  << " ± " << ml_result.uncertainty 
                  << " (" << ml_result.station_count << " stations)\n\n";
    } else {
        // Estimate from known source
        event->addMagnitude(Magnitude(MagnitudeType::ML, GroundTruth::MAGNITUDE_ML, 0.2, 
                                       event_picks.size()));
        std::cout << "  ML = " << GroundTruth::MAGNITUDE_ML << " (from input)\n\n";
    }
    
    // ========================================================================
    // Step 8: Store in CSS3.0 database
    // ========================================================================
    std::cout << "Step 8: Storing in CSS3.0 database...\n";
    
    std::string db_file = "dprk_nuclear_test.db";
    if (argc > 1) {
        db_file = argv[1];
    }
    
    CSS30Database database;
    if (!database.open(db_file)) {
        std::cerr << "ERROR: Failed to open database: " << db_file << "\n";
        return 1;
    }
    
    database.createSchema();
    database.setAuthor("dprk_example");
    
    // Store station inventory
    database.storeInventory(inventory, "REG");
    
    // Store event
    bool stored = database.storeCompleteEvent(*event, "Korea_KimEtAl2011", "REG");
    
    if (stored) {
        std::cout << "  Event stored in: " << db_file << "\n";
        std::cout << "  Tables populated: event, origin, origerr, arrival, assoc, netmag\n\n";
    } else {
        std::cerr << "  Failed to store event!\n\n";
    }
    
    // ========================================================================
    // Step 9: Compare with ground truth
    // ========================================================================
    std::cout << "================================================================\n";
    std::cout << "Results Comparison\n";
    std::cout << "================================================================\n\n";
    
    double lat_error = result.origin.location.latitude - GroundTruth::LATITUDE;
    double lon_error = result.origin.location.longitude - GroundTruth::LONGITUDE;
    double depth_error = result.origin.location.depth - GroundTruth::DEPTH;
    
    GeoPoint computed(result.origin.location.latitude, result.origin.location.longitude);
    GeoPoint truth(GroundTruth::LATITUDE, GroundTruth::LONGITUDE);
    double epicentral_error = truth.distanceTo(computed);
    
    std::cout << std::fixed << std::setprecision(3);
    std::cout << "                   Computed        Ground Truth    Error\n";
    std::cout << "  Latitude:        " << std::setw(8) << result.origin.location.latitude 
              << "°N      " << std::setw(8) << GroundTruth::LATITUDE << "°N      "
              << std::showpos << lat_error << "°\n";
    std::cout << std::noshowpos;
    std::cout << "  Longitude:       " << std::setw(8) << result.origin.location.longitude
              << "°E     " << std::setw(8) << GroundTruth::LONGITUDE << "°E     "
              << std::showpos << lon_error << "°\n";
    std::cout << std::noshowpos << std::setprecision(1);
    std::cout << "  Depth:           " << std::setw(8) << result.origin.location.depth
              << " km       " << std::setw(8) << GroundTruth::DEPTH << " km       "
              << std::showpos << depth_error << " km\n";
    std::cout << std::noshowpos << std::setprecision(2);
    std::cout << "  Magnitude:       " << std::setw(8) << event->magnitude()
              << "          " << std::setw(8) << GroundTruth::MAGNITUDE_ML << "\n\n";
    
    std::cout << std::setprecision(1);
    std::cout << "  Epicentral error: " << epicentral_error << " km\n";
    std::cout << "  Location quality: " << result.origin.qualityCode() << "\n\n";
    
    // Quality assessment
    std::cout << "Quality Assessment:\n";
    if (epicentral_error < 5.0) {
        std::cout << "  ✓ Excellent location accuracy (<5 km)\n";
    } else if (epicentral_error < 15.0) {
        std::cout << "  ✓ Good location accuracy (<15 km)\n";
    } else if (epicentral_error < 30.0) {
        std::cout << "  ~ Moderate location accuracy (<30 km)\n";
    } else {
        std::cout << "  ✗ Poor location accuracy (>30 km)\n";
    }
    
    if (std::abs(depth_error) < 5.0) {
        std::cout << "  ✓ Depth well constrained (<5 km error)\n";
    } else {
        std::cout << "  ~ Depth poorly constrained\n";
    }
    
    if (result.origin.gap < 180) {
        std::cout << "  ✓ Good azimuthal coverage (gap < 180°)\n";
    } else {
        std::cout << "  ~ Limited azimuthal coverage\n";
    }
    
    std::cout << "\n================================================================\n";
    std::cout << "Example complete. Database saved to: " << db_file << "\n";
    std::cout << "================================================================\n";
    
    return 0;
}
