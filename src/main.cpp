/**
 * SeisProc - Real-time Seismic Event Processing System
 * 
 * Main application that:
 * 1. Connects to SeedLink servers for real-time data
 * 2. Picks phase arrivals using STA/LTA and AIC
 * 3. Associates picks into events
 * 4. Locates events using grid search + Geiger inversion
 * 5. Calculates magnitudes (ML, Mw)
 * 
 * Inspired by SeisComP architecture.
 */

#include <iostream>
#include <signal.h>
#include <atomic>
#include <thread>
#include <chrono>

#include "seisproc/core/types.hpp"
#include "seisproc/core/config.hpp"
#include "seisproc/core/station.hpp"
#include "seisproc/core/velocity_model.hpp"
#include "seisproc/seedlink/seedlink_client.hpp"
#include "seisproc/picker/stalta_picker.hpp"
#include "seisproc/picker/aic_picker.hpp"
#include "seisproc/associator/phase_associator.hpp"
#include "seisproc/locator/grid_search.hpp"
#include "seisproc/locator/geiger.hpp"
#include "seisproc/magnitude/local_magnitude.hpp"
#include "seisproc/magnitude/moment_magnitude.hpp"

using namespace seisproc;

// Global shutdown flag
std::atomic<bool> g_running(true);

void signalHandler(int signum) {
    std::cout << "\nReceived signal " << signum << ", shutting down..." << std::endl;
    g_running = false;
}

void printUsage(const char* progname) {
    std::cout << "Usage: " << progname << " [options]\n\n";
    std::cout << "Options:\n";
    std::cout << "  -c, --config <file>    Configuration file (default: seisproc.conf)\n";
    std::cout << "  -s, --server <host>    SeedLink server host (default: localhost)\n";
    std::cout << "  -p, --port <port>      SeedLink server port (default: 18000)\n";
    std::cout << "  -S, --stations <file>  Station inventory file\n";
    std::cout << "  -v, --velocity <file>  Velocity model file\n";
    std::cout << "  -h, --help             Show this help message\n";
    std::cout << "\n";
    std::cout << "Example:\n";
    std::cout << "  " << progname << " -s rtserve.iris.washington.edu -p 18000\n";
}

/**
 * SeisProc Application
 */
class SeisProc {
public:
    SeisProc() 
        : picker_(std::make_shared<STALTAPicker>())
        , locator_(std::make_shared<GeigerLocator>())
        , ml_calculator_(std::make_shared<LocalMagnitude>())
        , mw_calculator_(std::make_shared<MomentMagnitude>())
    {
        // Initialize with default velocity model
        velocity_model_ = VelocityModel1D::simpleThreeLayer();
        
        // Configure picker
        picker_->setParameter("sta_length", 0.5);
        picker_->setParameter("lta_length", 10.0);
        picker_->setParameter("trigger_ratio", 3.5);
        
        // Configure locator
        locator_->setVelocityModel(velocity_model_);
        
        // Configure associator
        associator_.setVelocityModel(velocity_model_);
        associator_.setTimeWindow(60.0);
        associator_.setMinStations(3);
        associator_.setMinPhases(4);
    }
    
    bool loadConfig(const std::string& filename) {
        if (!config_.loadFromFile(filename)) {
            std::cerr << "Failed to load config: " << filename << std::endl;
            return false;
        }
        
        // Apply configuration
        if (config_.has("picker.sta_length")) {
            picker_->setParameter("sta_length", config_.getDouble("picker.sta_length"));
        }
        if (config_.has("picker.lta_length")) {
            picker_->setParameter("lta_length", config_.getDouble("picker.lta_length"));
        }
        if (config_.has("picker.trigger_ratio")) {
            picker_->setParameter("trigger_ratio", config_.getDouble("picker.trigger_ratio"));
        }
        
        if (config_.has("locator.fixed_depth")) {
            locator_->setParameter("fixed_depth", config_.getDouble("locator.fixed_depth"));
        }
        
        if (config_.has("associator.min_stations")) {
            associator_.setMinStations(config_.getInt("associator.min_stations"));
        }
        if (config_.has("associator.min_phases")) {
            associator_.setMinPhases(config_.getInt("associator.min_phases"));
        }
        
        return true;
    }
    
    bool loadStations(const std::string& filename) {
        if (!stations_.loadFromFile(filename)) {
            std::cerr << "Failed to load stations: " << filename << std::endl;
            return false;
        }
        
        associator_.setStations(stations_);
        return true;
    }
    
    bool loadVelocityModel(const std::string& filename) {
        if (!velocity_model_.loadFromFile(filename)) {
            std::cerr << "Failed to load velocity model: " << filename << std::endl;
            return false;
        }
        
        locator_->setVelocityModel(velocity_model_);
        associator_.setVelocityModel(velocity_model_);
        return true;
    }
    
    bool connectToServer(const std::string& host, int port) {
        std::cout << "Connecting to SeedLink server " << host << ":" << port << std::endl;
        
        if (!client_.connect(host, port)) {
            std::cerr << "Failed to connect to server" << std::endl;
            return false;
        }
        
        // Set up callbacks
        client_.setWaveformCallback([this](WaveformPtr waveform) {
            processWaveform(waveform);
        });
        
        client_.setErrorCallback([](const std::string& error) {
            std::cerr << "SeedLink error: " << error << std::endl;
        });
        
        return true;
    }
    
    void addStream(const std::string& network, const std::string& station,
                   const std::string& selector = "BH?") {
        client_.addStream(network, station, selector);
    }
    
    bool start() {
        if (!client_.startStreaming()) {
            std::cerr << "Failed to start streaming" << std::endl;
            return false;
        }
        
        std::cout << "Started real-time processing" << std::endl;
        
        // Processing loop
        while (g_running) {
            // Process pending picks
            auto events = associator_.process();
            
            for (auto& event : events) {
                processEvent(event);
            }
            
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
        }
        
        client_.stopStreaming();
        return true;
    }
    
    void stop() {
        client_.disconnect();
    }
    
private:
    Config config_;
    StationInventory stations_;
    VelocityModel1D velocity_model_;
    SeedLinkClient client_;
    MultiStreamBuffer buffer_;
    
    std::shared_ptr<STALTAPicker> picker_;
    PhaseAssociator associator_;
    std::shared_ptr<GeigerLocator> locator_;
    std::shared_ptr<LocalMagnitude> ml_calculator_;
    std::shared_ptr<MomentMagnitude> mw_calculator_;
    
    std::mutex mutex_;
    std::vector<PickPtr> recent_picks_;
    
    void processWaveform(WaveformPtr waveform) {
        // Add to buffer
        buffer_.addWaveform(waveform);
        
        // Run picker
        auto picks = picker_->pick(*waveform);
        
        for (const auto& pick_result : picks) {
            auto pick = std::make_shared<Pick>();
            pick->stream_id = waveform->streamId();
            pick->time = pick_result.time;
            pick->phase_type = pick_result.phase_type;
            pick->snr = pick_result.snr;
            pick->amplitude = pick_result.amplitude;
            pick->is_automatic = true;
            pick->method = "STA/LTA";
            
            // Quality based on SNR
            if (pick_result.snr >= 10) {
                pick->quality = PickQuality::Impulsive;
            } else if (pick_result.snr >= 5) {
                pick->quality = PickQuality::Emergent;
            } else {
                pick->quality = PickQuality::Questionable;
            }
            
            std::cout << "Pick: " << pick->stream_id.toString() 
                      << " " << phaseTypeToString(pick->phase_type)
                      << " SNR=" << pick->snr << std::endl;
            
            // Add to associator
            associator_.addPick(pick);
            
            std::lock_guard<std::mutex> lock(mutex_);
            recent_picks_.push_back(pick);
        }
    }
    
    void processEvent(EventPtr event) {
        std::cout << "\n*** NEW EVENT DETECTED ***" << std::endl;
        
        // Get picks for this event
        std::vector<PickPtr> event_picks;
        for (const auto& arr : event->preferredOrigin().arrivals) {
            event_picks.push_back(arr.pick);
        }
        
        // Locate event
        auto location_result = locator_->locate(event_picks, stations_);
        
        if (location_result.converged) {
            event->origins().back() = location_result.origin;
            
            std::cout << "Location: " << location_result.origin.location.latitude << "° N, "
                      << location_result.origin.location.longitude << "° E, "
                      << location_result.origin.location.depth << " km depth" << std::endl;
            std::cout << "RMS: " << location_result.origin.rms << " s" << std::endl;
            std::cout << "Quality: " << location_result.origin.qualityCode() << std::endl;
            
            // Calculate magnitude
            std::map<StreamID, WaveformPtr> waveforms;
            
            for (const auto& pick : event_picks) {
                auto buf = buffer_.getBuffer(pick->stream_id);
                if (buf) {
                    // Get waveform around the event
                    TimePoint start = location_result.origin.time - std::chrono::seconds(10);
                    TimePoint end = location_result.origin.time + std::chrono::seconds(60);
                    waveforms[pick->stream_id] = buf->getWaveform(start, end);
                }
            }
            
            if (!waveforms.empty()) {
                // Calculate ML
                auto ml_result = ml_calculator_->calculate(
                    location_result.origin, waveforms, stations_);
                
                if (ml_result.station_count > 0) {
                    event->addMagnitude(Magnitude(MagnitudeType::ML, 
                                                   ml_result.value, 
                                                   ml_result.uncertainty,
                                                   ml_result.station_count));
                    
                    std::cout << "ML: " << ml_result.value 
                              << " ± " << ml_result.uncertainty 
                              << " (" << ml_result.station_count << " stations)" << std::endl;
                }
                
                // Calculate Mw for larger events
                if (ml_result.value > 3.0) {
                    auto mw_result = mw_calculator_->calculate(
                        location_result.origin, waveforms, stations_);
                    
                    if (mw_result.station_count > 0) {
                        event->addMagnitude(Magnitude(MagnitudeType::Mw,
                                                       mw_result.value,
                                                       mw_result.uncertainty,
                                                       mw_result.station_count));
                        
                        std::cout << "Mw: " << mw_result.value 
                                  << " ± " << mw_result.uncertainty << std::endl;
                    }
                }
            }
        } else {
            std::cout << "Location failed to converge" << std::endl;
        }
        
        std::cout << event->summary() << std::endl;
        std::cout << "****************************\n" << std::endl;
    }
};

int main(int argc, char* argv[]) {
    // Set up signal handlers
    signal(SIGINT, signalHandler);
    signal(SIGTERM, signalHandler);
    
    // Parse command line arguments
    std::string config_file = "seisproc.conf";
    std::string server_host = "localhost";
    int server_port = 18000;
    std::string stations_file;
    std::string velocity_file;
    
    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        
        if (arg == "-h" || arg == "--help") {
            printUsage(argv[0]);
            return 0;
        } else if ((arg == "-c" || arg == "--config") && i + 1 < argc) {
            config_file = argv[++i];
        } else if ((arg == "-s" || arg == "--server") && i + 1 < argc) {
            server_host = argv[++i];
        } else if ((arg == "-p" || arg == "--port") && i + 1 < argc) {
            server_port = std::stoi(argv[++i]);
        } else if ((arg == "-S" || arg == "--stations") && i + 1 < argc) {
            stations_file = argv[++i];
        } else if ((arg == "-v" || arg == "--velocity") && i + 1 < argc) {
            velocity_file = argv[++i];
        }
    }
    
    std::cout << "=======================================\n";
    std::cout << "SeisProc - Real-time Seismic Processing\n";
    std::cout << "=======================================\n\n";
    
    SeisProc app;
    
    // Load configuration
    app.loadConfig(config_file);
    
    // Load stations if provided
    if (!stations_file.empty()) {
        app.loadStations(stations_file);
    }
    
    // Load velocity model if provided
    if (!velocity_file.empty()) {
        app.loadVelocityModel(velocity_file);
    }
    
    // Connect to server
    if (!app.connectToServer(server_host, server_port)) {
        return 1;
    }
    
    // Add some default streams (example: IRIS network)
    app.addStream("IU", "ANMO", "BH?");
    app.addStream("IU", "CCM", "BH?");
    app.addStream("IU", "HRV", "BH?");
    app.addStream("IU", "SSPA", "BH?");
    app.addStream("US", "ACSO", "BH?");
    
    // Start processing
    if (!app.start()) {
        return 1;
    }
    
    app.stop();
    
    std::cout << "SeisProc shutdown complete" << std::endl;
    return 0;
}
