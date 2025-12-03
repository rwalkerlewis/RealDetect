/**
 * SeisProc Simulator
 * 
 * Generates synthetic seismic waveforms for testing the processing system.
 * Creates realistic P and S wave arrivals with noise.
 */

#include <iostream>
#include <random>
#include <chrono>
#include <thread>
#include <cmath>
#include <fstream>
#include <iomanip>

#include "seisproc/core/types.hpp"
#include "seisproc/core/waveform.hpp"
#include "seisproc/core/station.hpp"
#include "seisproc/core/event.hpp"
#include "seisproc/core/velocity_model.hpp"
#include "seisproc/picker/stalta_picker.hpp"
#include "seisproc/picker/aic_picker.hpp"
#include "seisproc/associator/phase_associator.hpp"
#include "seisproc/locator/grid_search.hpp"
#include "seisproc/locator/geiger.hpp"
#include "seisproc/magnitude/local_magnitude.hpp"

using namespace seisproc;

class SeismicSimulator {
public:
    SeismicSimulator() : gen_(std::random_device{}()) {
        velocity_model_ = VelocityModel1D::simpleThreeLayer();
    }
    
    // Generate synthetic event
    void generateEvent(double lat, double lon, double depth, double magnitude) {
        event_lat_ = lat;
        event_lon_ = lon;
        event_depth_ = depth;
        event_mag_ = magnitude;
        event_time_ = std::chrono::system_clock::now();
        
        std::cout << "\n=== Generating synthetic event ===" << std::endl;
        std::cout << "Location: " << lat << "°N, " << lon << "°E, " << depth << " km" << std::endl;
        std::cout << "Magnitude: " << magnitude << std::endl;
    }
    
    // Generate waveform for a station
    WaveformPtr generateWaveform(const Station& station, double sample_rate = 100.0,
                                   double duration = 120.0) {
        GeoPoint event_loc(event_lat_, event_lon_, event_depth_);
        double distance = event_loc.distanceTo(station.location());
        
        // Calculate travel times
        double vp = 6.0, vs = 3.5;  // km/s
        double hypo_dist = std::sqrt(distance * distance + event_depth_ * event_depth_);
        double p_time = hypo_dist / vp;
        double s_time = hypo_dist / vs;
        
        // Create waveform
        StreamID id(station.network(), station.code(), "00", "BHZ");
        auto waveform = std::make_shared<Waveform>(id, sample_rate, event_time_);
        
        size_t n_samples = static_cast<size_t>(duration * sample_rate);
        SampleVector& data = waveform->data();
        data.resize(n_samples);
        
        // Generate noise background
        std::normal_distribution<> noise(0.0, 100.0);
        for (auto& s : data) {
            s = noise(gen_);
        }
        
        // Add P-wave arrival
        size_t p_idx = static_cast<size_t>(p_time * sample_rate);
        if (p_idx < n_samples) {
            addArrival(data, p_idx, sample_rate, 
                       getPAmplitude(distance, event_mag_),
                       5.0,   // Dominant frequency
                       2.0);  // Duration
        }
        
        // Add S-wave arrival
        size_t s_idx = static_cast<size_t>(s_time * sample_rate);
        if (s_idx < n_samples) {
            addArrival(data, s_idx, sample_rate,
                       getSAmplitude(distance, event_mag_),
                       2.5,   // Dominant frequency
                       5.0);  // Duration
        }
        
        // Add some surface waves for larger events
        if (event_mag_ > 4.0 && distance > 100) {
            double surface_time = distance / 3.0;  // ~3 km/s
            size_t surf_idx = static_cast<size_t>(surface_time * sample_rate);
            if (surf_idx < n_samples) {
                addArrival(data, surf_idx, sample_rate,
                           getSAmplitude(distance, event_mag_) * 1.5,
                           0.1,   // Low frequency
                           15.0); // Long duration
            }
        }
        
        std::cout << "Station " << station.code() 
                  << ": distance=" << distance << " km"
                  << ", P=" << p_time << "s, S=" << s_time << "s" << std::endl;
        
        return waveform;
    }
    
    // Create example station network
    StationInventory createTestNetwork() {
        StationInventory inventory;
        
        // Create a network around Southern California
        double center_lat = 34.0;
        double center_lon = -118.0;
        
        std::vector<std::pair<double, double>> offsets = {
            {0.0, 0.0},
            {0.5, 0.3},
            {-0.4, 0.6},
            {0.3, -0.5},
            {-0.6, -0.3},
            {0.8, 0.1},
            {-0.2, 0.9},
            {0.6, 0.7}
        };
        
        int i = 1;
        for (const auto& [dlat, dlon] : offsets) {
            std::string code = "S" + std::to_string(i);
            auto sta = std::make_shared<Station>("CI", code,
                center_lat + dlat, center_lon + dlon, 100.0);
            
            Channel bhz;
            bhz.code = "BHZ";
            bhz.sample_rate = 100.0;
            bhz.dip = -90.0;
            sta->addChannel(bhz);
            
            Channel bhn;
            bhn.code = "BHN";
            bhn.sample_rate = 100.0;
            bhn.azimuth = 0.0;
            bhn.dip = 0.0;
            sta->addChannel(bhn);
            
            Channel bhe;
            bhe.code = "BHE";
            bhe.sample_rate = 100.0;
            bhe.azimuth = 90.0;
            bhe.dip = 0.0;
            sta->addChannel(bhe);
            
            inventory.addStation(sta);
            i++;
        }
        
        return inventory;
    }

private:
    std::mt19937 gen_;
    VelocityModel1D velocity_model_;
    
    double event_lat_, event_lon_, event_depth_, event_mag_;
    TimePoint event_time_;
    
    void addArrival(SampleVector& data, size_t start_idx, double sample_rate,
                     double amplitude, double freq, double duration) {
        size_t n_samples = static_cast<size_t>(duration * sample_rate);
        double dt = 1.0 / sample_rate;
        
        for (size_t i = 0; i < n_samples && start_idx + i < data.size(); i++) {
            double t = i * dt;
            
            // Ricker wavelet envelope
            double sigma = 1.0 / (2.0 * M_PI * freq);
            double env = (1.0 - 0.5 * t / duration) * std::exp(-t * t / (4 * sigma * sigma));
            
            // Oscillation
            double osc = std::sin(2.0 * M_PI * freq * t);
            
            // Add some randomness
            std::normal_distribution<> noise(0.0, 0.1);
            double phase_noise = noise(gen_);
            
            data[start_idx + i] += amplitude * env * std::sin(2.0 * M_PI * freq * t + phase_noise);
        }
    }
    
    double getPAmplitude(double distance, double magnitude) {
        // Simplified amplitude-distance-magnitude relationship
        // log(A) = M - 1.66*log(D) - 2.0
        double log_amp = magnitude - 1.66 * std::log10(distance) - 2.0;
        return std::pow(10.0, log_amp) * 1e6;  // Convert to counts
    }
    
    double getSAmplitude(double distance, double magnitude) {
        // S-wave amplitude typically 1.5-2x P-wave
        return getPAmplitude(distance, magnitude) * 1.7;
    }
};

void runSimulation() {
    std::cout << "========================================" << std::endl;
    std::cout << "SeisProc Simulator - Testing Processing" << std::endl;
    std::cout << "========================================" << std::endl;
    
    SeismicSimulator sim;
    
    // Create test network
    StationInventory stations = sim.createTestNetwork();
    std::cout << "\nCreated network with " << stations.size() << " stations" << std::endl;
    
    // Generate a synthetic event
    sim.generateEvent(34.2, -117.8, 12.0, 4.5);
    
    // Generate waveforms for all stations
    std::vector<WaveformPtr> waveforms;
    for (const auto& [key, sta] : stations.stations()) {
        waveforms.push_back(sim.generateWaveform(*sta));
    }
    
    // Initialize processing components
    VelocityModel1D velocity_model = VelocityModel1D::simpleThreeLayer();
    
    STALTAPicker picker;
    picker.setParameter("sta_length", 0.5);
    picker.setParameter("lta_length", 10.0);
    picker.setParameter("trigger_ratio", 3.0);
    
    PhaseAssociator associator;
    associator.setStations(stations);
    associator.setVelocityModel(velocity_model);
    associator.setMinStations(3);
    associator.setMinPhases(4);
    
    GeigerLocator locator;
    locator.setVelocityModel(velocity_model);
    locator.setMaxIterations(20);
    
    // Process waveforms
    std::cout << "\n=== Phase Picking ===" << std::endl;
    
    std::vector<PickPtr> all_picks;
    
    for (const auto& waveform : waveforms) {
        auto picks = picker.pick(*waveform);
        
        for (const auto& pick_result : picks) {
            auto pick = std::make_shared<Pick>();
            pick->stream_id = waveform->streamId();
            pick->time = pick_result.time;
            pick->phase_type = pick_result.phase_type;
            pick->snr = pick_result.snr;
            pick->amplitude = pick_result.amplitude;
            pick->is_automatic = true;
            pick->method = "STA/LTA";
            pick->quality = PickQuality::Emergent;
            
            std::cout << "Pick: " << pick->stream_id.toString()
                      << " " << phaseTypeToString(pick->phase_type)
                      << " SNR=" << pick->snr << std::endl;
            
            all_picks.push_back(pick);
            associator.addPick(pick);
        }
    }
    
    // Associate picks into events
    std::cout << "\n=== Event Association ===" << std::endl;
    
    auto events = associator.process();
    std::cout << "Found " << events.size() << " events" << std::endl;
    
    // Locate events
    std::cout << "\n=== Event Location ===" << std::endl;
    
    for (auto& event : events) {
        std::vector<PickPtr> event_picks;
        for (const auto& arr : event->preferredOrigin().arrivals) {
            event_picks.push_back(arr.pick);
        }
        
        auto result = locator.locate(event_picks, stations);
        
        if (result.converged) {
            event->origins().back() = result.origin;
            
            std::cout << "\nLocated event:" << std::endl;
            std::cout << "  Location: " << result.origin.location.latitude << "°N, "
                      << result.origin.location.longitude << "°E" << std::endl;
            std::cout << "  Depth: " << result.origin.location.depth << " km" << std::endl;
            std::cout << "  RMS: " << result.origin.rms << " s" << std::endl;
            std::cout << "  Phases: " << result.origin.phase_count << std::endl;
            std::cout << "  Stations: " << result.origin.station_count << std::endl;
            std::cout << "  Quality: " << result.origin.qualityCode() << std::endl;
            std::cout << "  Iterations: " << result.iterations << std::endl;
            
            // Compare with true location
            std::cout << "\n  True location: 34.200°N, -117.800°E, 12.0 km" << std::endl;
            std::cout << "  Latitude error: " 
                      << std::abs(result.origin.location.latitude - 34.2) << "°" << std::endl;
            std::cout << "  Longitude error: "
                      << std::abs(result.origin.location.longitude - (-117.8)) << "°" << std::endl;
            std::cout << "  Depth error: "
                      << std::abs(result.origin.location.depth - 12.0) << " km" << std::endl;
        } else {
            std::cout << "Location failed to converge" << std::endl;
        }
        
        // Calculate magnitude
        std::cout << "\n=== Magnitude Calculation ===" << std::endl;
        
        std::map<StreamID, WaveformPtr> wf_map;
        for (const auto& wf : waveforms) {
            wf_map[wf->streamId()] = wf;
        }
        
        LocalMagnitude ml_calc;
        auto ml_result = ml_calc.calculate(result.origin, wf_map, stations);
        
        if (ml_result.station_count > 0) {
            event->addMagnitude(Magnitude(MagnitudeType::ML, ml_result.value,
                                           ml_result.uncertainty, ml_result.station_count));
            
            std::cout << "  ML: " << ml_result.value
                      << " ± " << ml_result.uncertainty
                      << " (" << ml_result.station_count << " stations)" << std::endl;
            std::cout << "  True magnitude: 4.5" << std::endl;
        }
        
        // Print event summary
        std::cout << "\n" << event->summary() << std::endl;
    }
    
    std::cout << "\n========================================" << std::endl;
    std::cout << "Simulation complete" << std::endl;
    std::cout << "========================================" << std::endl;
}

void printUsage(const char* progname) {
    std::cout << "SeisProc Simulator\n\n";
    std::cout << "Generates synthetic seismic events for testing.\n\n";
    std::cout << "Usage: " << progname << " [options]\n\n";
    std::cout << "Options:\n";
    std::cout << "  -m, --magnitude <M>    Event magnitude (default: 4.5)\n";
    std::cout << "  -d, --depth <km>       Event depth (default: 12.0)\n";
    std::cout << "  --lat <degrees>        Event latitude (default: 34.2)\n";
    std::cout << "  --lon <degrees>        Event longitude (default: -117.8)\n";
    std::cout << "  -h, --help             Show this help\n";
}

int main(int argc, char* argv[]) {
    // Parse command line
    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "-h" || arg == "--help") {
            printUsage(argv[0]);
            return 0;
        }
    }
    
    runSimulation();
    
    return 0;
}
