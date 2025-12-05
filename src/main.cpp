/**
 * RealDetect - Real-time Seismic Event Processing System
 * 
 * Main application that:
 * 1. Connects to SeedLink servers for real-time data OR plays back MiniSEED files
 * 2. Picks phase arrivals using STA/LTA, AIC, or ML detectors
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
#include <filesystem>
#include <algorithm>

#include "realdetect/core/types.hpp"
#include "realdetect/core/config.hpp"
#include "realdetect/core/station.hpp"
#include "realdetect/core/velocity_model.hpp"
#include "realdetect/core/regional_velocity_model.hpp"
#include "realdetect/core/miniseed.hpp"
#include "realdetect/database/css30_database.hpp"
#include "realdetect/seedlink/seedlink_client.hpp"
#include "realdetect/picker/stalta_picker.hpp"
#include "realdetect/picker/aic_picker.hpp"
#include "realdetect/picker/ml_picker.hpp"
#include "realdetect/associator/phase_associator.hpp"
#include "realdetect/locator/grid_search.hpp"
#include "realdetect/locator/geiger.hpp"
#include "realdetect/magnitude/local_magnitude.hpp"
#include "realdetect/magnitude/moment_magnitude.hpp"

using namespace realdetect;

// Global shutdown flag
std::atomic<bool> g_running(true);

void signalHandler(int signum) {
    std::cout << "\nReceived signal " << signum << ", shutting down..." << std::endl;
    g_running = false;
}

void printUsage(const char* progname) {
    std::cout << "Usage: " << progname << " [options]\n\n";
    std::cout << "Options:\n";
    std::cout << "  -c, --config <file>    Configuration file (default: realdetect.conf)\n";
    std::cout << "  -s, --server <host>    SeedLink server host (default: localhost)\n";
    std::cout << "  -p, --port <port>      SeedLink server port (default: 18000)\n";
    std::cout << "  -S, --stations <file>  Station inventory file\n";
    std::cout << "  -v, --velocity <file>  Velocity model file\n";
    std::cout << "  -V, --velocity-dir <dir> Directory with regional velocity models\n";
    std::cout << "  -d, --database <file>  CSS3.0 database file for output\n";
    std::cout << "  -m, --miniseed <file>  MiniSEED file for playback mode\n";
    std::cout << "  -M, --miniseed-dir <dir> Directory of MiniSEED files for playback\n";
    std::cout << "  --playback-speed <n>   Playback speed multiplier (default: 1.0, 0 = fast as possible)\n";
    std::cout << "  --picker <type>        Picker type: stalta, aic, ml (default: stalta)\n";
    std::cout << "  --ml-model <file>      ML model file for ML picker\n";
    std::cout << "  -h, --help             Show this help message\n";
    std::cout << "\n";
    std::cout << "Modes:\n";
    std::cout << "  Real-time:  Connect to SeedLink server for live data streaming\n";
    std::cout << "  Playback:   Process MiniSEED files (use -m or -M options)\n";
    std::cout << "\n";
    std::cout << "Picker Types:\n";
    std::cout << "  stalta      STA/LTA detector with AIC refinement (default)\n";
    std::cout << "  aic         AIC picker only\n";
    std::cout << "  ml          Machine learning detector (requires model file)\n";
    std::cout << "\n";
    std::cout << "Database Output:\n";
    std::cout << "  Events are stored in CSS3.0 schema format (SQLite database)\n";
    std::cout << "  Tables: event, origin, origerr, arrival, assoc, netmag, stamag\n";
    std::cout << "\n";
    std::cout << "Examples:\n";
    std::cout << "  # Real-time processing\n";
    std::cout << "  " << progname << " -s rtserve.iris.washington.edu -p 18000 -d catalog.db\n";
    std::cout << "\n";
    std::cout << "  # Playback mode with MiniSEED file\n";
    std::cout << "  " << progname << " -m data.mseed -d catalog.db\n";
    std::cout << "\n";
    std::cout << "  # ML picker with custom model\n";
    std::cout << "  " << progname << " --picker ml --ml-model phasenet.onnx -m data.mseed\n";
}

/**
 * PickerType - Available picker algorithms
 */
enum class PickerType {
    STALTA,
    AIC,
    ML
};

/**
 * RealDetect Application
 */
class RealDetect {
public:
    RealDetect() 
        : db_enabled_(false)
        , use_regional_models_(false)
        , playback_mode_(false)
        , playback_speed_(1.0)
        , picker_type_(PickerType::STALTA)
        , picker_(std::make_shared<STALTAPicker>())
        , locator_(std::make_shared<GeigerLocator>())
        , ml_calculator_(std::make_shared<LocalMagnitude>())
        , mw_calculator_(std::make_shared<MomentMagnitude>())
    {
        // Initialize with default velocity model
        velocity_model_ = VelocityModel1D::simpleThreeLayer();
        velocity_manager_.setDefaultModel(velocity_model_);
        
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
    
    void setPickerType(PickerType type, const std::string& ml_model_path = "") {
        picker_type_ = type;
        switch (type) {
            case PickerType::STALTA:
                picker_ = std::make_shared<STALTAPicker>();
                picker_->setParameter("sta_length", 0.5);
                picker_->setParameter("lta_length", 10.0);
                picker_->setParameter("trigger_ratio", 3.5);
                break;
            case PickerType::AIC:
                picker_ = std::make_shared<AICPicker>();
                break;
            case PickerType::ML:
                {
                    auto ml_picker = std::make_shared<MLPicker>();
                    if (!ml_model_path.empty()) {
                        ml_picker->loadModel(ml_model_path);
                    }
                    picker_ = ml_picker;
                }
                break;
        }
    }
    
    void setPlaybackMode(bool enabled) { playback_mode_ = enabled; }
    void setPlaybackSpeed(double speed) { playback_speed_ = speed; }
    
    bool loadConfig(const std::string& filename) {
        if (!config_.loadFromFile(filename)) {
            std::cerr << "Failed to load config: " << filename << std::endl;
            return false;
        }
        
        // Apply picker configuration
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
        
        // Load velocity model configuration
        if (config_.has("velocity_model.model")) {
            std::string model_name = config_.getString("velocity_model.model");
            if (model_name == "iasp91") {
                velocity_model_ = VelocityModel1D::iasp91();
            } else if (model_name == "ak135") {
                velocity_model_ = VelocityModel1D::ak135();
            } else if (model_name == "simple3layer") {
                velocity_model_ = VelocityModel1D::simpleThreeLayer();
            }
            velocity_manager_.setDefaultModel(velocity_model_);
            locator_->setVelocityModel(velocity_model_);
            associator_.setVelocityModel(velocity_model_);
        }
        
        // Regional velocity models
        if (config_.getBool("velocity_model.regional_models_enabled", false)) {
            use_regional_models_ = true;
            std::string models_dir = config_.getString("velocity_model.regional_models_dir", 
                                                        "config/velocity_models");
            if (!velocity_manager_.loadFromDirectory(models_dir)) {
                std::cerr << "Warning: Failed to load regional velocity models from " 
                          << models_dir << std::endl;
            } else {
                std::cout << "Loaded " << velocity_manager_.modelCount() 
                          << " regional velocity models" << std::endl;
            }
        }
        
        // Database configuration
        if (config_.getBool("database.enabled", false)) {
            std::string db_file = config_.getString("database.file", "realdetect_catalog.db");
            if (initDatabase(db_file)) {
                db_author_ = config_.getString("database.author", "realdetect");
                db_network_ = config_.getString("database.network", "XX");
                database_.setAuthor(db_author_);
            }
        }
        
        return true;
    }
    
    bool initDatabase(const std::string& filename) {
        if (!database_.open(filename)) {
            std::cerr << "Failed to open database: " << filename << std::endl;
            return false;
        }
        
        if (config_.getBool("database.auto_create_schema", true)) {
            if (!database_.createSchema()) {
                std::cerr << "Failed to create database schema" << std::endl;
                return false;
            }
        }
        
        db_enabled_ = true;
        std::cout << "CSS3.0 database opened: " << filename << std::endl;
        return true;
    }
    
    bool openDatabase(const std::string& filename) {
        if (!database_.open(filename)) {
            std::cerr << "Failed to open database: " << filename << std::endl;
            return false;
        }
        database_.createSchema();
        db_enabled_ = true;
        std::cout << "CSS3.0 database opened: " << filename << std::endl;
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
        
        velocity_manager_.setDefaultModel(velocity_model_);
        locator_->setVelocityModel(velocity_model_);
        associator_.setVelocityModel(velocity_model_);
        return true;
    }
    
    bool loadRegionalModels(const std::string& directory) {
        if (!velocity_manager_.loadFromDirectory(directory)) {
            std::cerr << "Failed to load regional velocity models from: " << directory << std::endl;
            return false;
        }
        use_regional_models_ = true;
        std::cout << "Loaded " << velocity_manager_.modelCount() 
                  << " regional velocity models" << std::endl;
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
        if (playback_mode_) {
            return startPlayback();
        }
        
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
    
    bool startPlayback() {
        std::cout << "Starting playback mode" << std::endl;
        
        if (playback_waveforms_.empty()) {
            std::cerr << "No waveforms loaded for playback" << std::endl;
            return false;
        }
        
        // Sort waveforms by start time
        std::vector<WaveformPtr> sorted_waveforms = playback_waveforms_;
        std::sort(sorted_waveforms.begin(), sorted_waveforms.end(),
            [](const WaveformPtr& a, const WaveformPtr& b) {
                return a->startTime() < b->startTime();
            });
        
        TimePoint playback_start = sorted_waveforms.front()->startTime();
        auto real_start = std::chrono::steady_clock::now();
        
        std::cout << "Playback starting from " << sorted_waveforms.size() 
                  << " waveforms" << std::endl;
        
        size_t processed = 0;
        for (const auto& waveform : sorted_waveforms) {
            if (!g_running) break;
            
            // Simulate timing if playback_speed > 0
            if (playback_speed_ > 0) {
                auto elapsed_data = std::chrono::duration_cast<std::chrono::milliseconds>(
                    waveform->startTime() - playback_start);
                auto target_elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
                    elapsed_data / playback_speed_);
                auto real_elapsed = std::chrono::steady_clock::now() - real_start;
                
                if (target_elapsed > real_elapsed) {
                    std::this_thread::sleep_for(target_elapsed - real_elapsed);
                }
            }
            
            // Process the waveform
            processWaveform(waveform);
            processed++;
            
            // Check for events periodically
            if (processed % 10 == 0) {
                auto events = associator_.process();
                for (auto& event : events) {
                    processEvent(event);
                }
            }
        }
        
        // Final event processing
        auto events = associator_.process();
        for (auto& event : events) {
            processEvent(event);
        }
        
        std::cout << "Playback complete. Processed " << processed << " waveforms." << std::endl;
        return true;
    }
    
    void stop() {
        client_.disconnect();
    }
    
    // Load MiniSEED file for playback
    bool loadMiniSeedFile(const std::string& filename) {
        MiniSeedReader reader;
        if (!reader.open(filename)) {
            std::cerr << "Failed to open MiniSEED file: " << filename << std::endl;
            return false;
        }
        
        auto waveforms = reader.toWaveforms();
        for (auto& wf : waveforms) {
            playback_waveforms_.push_back(wf);
        }
        
        std::cout << "Loaded " << waveforms.size() << " waveforms from " << filename << std::endl;
        playback_mode_ = true;
        return true;
    }
    
    // Load MiniSEED data from memory for playback
    bool loadMiniSeedData(const uint8_t* data, size_t length) {
        MiniSeedReader reader;
        if (!reader.parse(data, length)) {
            std::cerr << "Failed to parse MiniSEED data" << std::endl;
            return false;
        }
        
        auto waveforms = reader.toWaveforms();
        for (auto& wf : waveforms) {
            playback_waveforms_.push_back(wf);
        }
        
        std::cout << "Loaded " << waveforms.size() << " waveforms from memory block" << std::endl;
        playback_mode_ = true;
        return true;
    }
    
    // Load MiniSEED files from directory
    bool loadMiniSeedDirectory(const std::string& directory) {
        std::cout << "Loading MiniSEED files from directory: " << directory << std::endl;
        
        // Use filesystem to iterate directory
        size_t total_waveforms = 0;
        for (const auto& entry : std::filesystem::directory_iterator(directory)) {
            if (entry.is_regular_file()) {
                std::string ext = entry.path().extension().string();
                // Check for common miniseed extensions
                if (ext == ".mseed" || ext == ".miniseed" || ext == ".ms" || ext == ".seed") {
                    if (loadMiniSeedFile(entry.path().string())) {
                        total_waveforms += playback_waveforms_.size();
                    }
                }
            }
        }
        
        if (total_waveforms == 0) {
            std::cerr << "No MiniSEED files found in directory" << std::endl;
            return false;
        }
        
        playback_mode_ = true;
        return true;
    }
    
public:
    bool loadPlaybackFile(const std::string& filename) {
        return loadMiniSeedFile(filename);
    }
    
    bool loadPlaybackDirectory(const std::string& directory) {
        return loadMiniSeedDirectory(directory);
    }
    
    bool loadPlaybackData(const uint8_t* data, size_t length) {
        return loadMiniSeedData(data, length);
    }
    
private:
    Config config_;
    StationInventory stations_;
    VelocityModel1D velocity_model_;
    VelocityModelManager velocity_manager_;
    SeedLinkClient client_;
    MultiStreamBuffer buffer_;
    
    // CSS3.0 Database
    CSS30Database database_;
    bool db_enabled_;
    std::string db_author_;
    std::string db_network_;
    bool use_regional_models_;
    
    // Playback mode
    bool playback_mode_;
    double playback_speed_;
    std::vector<WaveformPtr> playback_waveforms_;
    
    // Picker
    PickerType picker_type_;
    std::shared_ptr<BasePicker> picker_;
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
            pick->method = picker_->name();
            
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
        
        // Select velocity model based on initial location estimate (if using regional models)
        std::string velocity_model_name = "default";
        if (use_regional_models_) {
            // Use initial event location to select velocity model
            const auto& initial_loc = event->preferredOrigin().location;
            velocity_model_name = velocity_manager_.getModelNameForLocation(initial_loc);
            const auto& regional_model = velocity_manager_.getModelForLocation(initial_loc);
            locator_->setVelocityModel(regional_model);
            std::cout << "Using velocity model: " << velocity_model_name << std::endl;
        }
        
        // Locate event
        auto location_result = locator_->locate(event_picks, stations_);
        
        if (location_result.converged) {
            // If using regional models, check if the final location needs a different model
            if (use_regional_models_) {
                std::string final_model = velocity_manager_.getModelNameForLocation(
                    location_result.origin.location);
                if (final_model != velocity_model_name) {
                    // Re-locate with the correct model
                    std::cout << "Relocating with velocity model: " << final_model << std::endl;
                    velocity_model_name = final_model;
                    const auto& regional_model = velocity_manager_.getModelForLocation(
                        location_result.origin.location);
                    locator_->setVelocityModel(regional_model);
                    location_result = locator_->locate(event_picks, stations_);
                }
            }
            
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
            
            // Store event in CSS3.0 database
            if (db_enabled_) {
                if (database_.storeCompleteEvent(*event, velocity_model_name, db_network_)) {
                    std::cout << "Event stored in database (evid: " << event->id() << ")" << std::endl;
                } else {
                    std::cerr << "Failed to store event in database" << std::endl;
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
    std::string config_file = "realdetect.conf";
    std::string server_host = "localhost";
    int server_port = 18000;
    std::string stations_file;
    std::string velocity_file;
    std::string velocity_dir;
    std::string database_file;
    std::string miniseed_file;
    std::string miniseed_dir;
    double playback_speed = 1.0;
    std::string picker_type = "stalta";
    std::string ml_model_file;
    
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
        } else if ((arg == "-V" || arg == "--velocity-dir") && i + 1 < argc) {
            velocity_dir = argv[++i];
        } else if ((arg == "-d" || arg == "--database") && i + 1 < argc) {
            database_file = argv[++i];
        } else if ((arg == "-m" || arg == "--miniseed") && i + 1 < argc) {
            miniseed_file = argv[++i];
        } else if ((arg == "-M" || arg == "--miniseed-dir") && i + 1 < argc) {
            miniseed_dir = argv[++i];
        } else if (arg == "--playback-speed" && i + 1 < argc) {
            playback_speed = std::stod(argv[++i]);
        } else if (arg == "--picker" && i + 1 < argc) {
            picker_type = argv[++i];
        } else if (arg == "--ml-model" && i + 1 < argc) {
            ml_model_file = argv[++i];
        }
    }
    
    std::cout << "==========================================\n";
    std::cout << "RealDetect - Real-time Seismic Processing\n";
    std::cout << "  with CSS3.0 Database Support\n";
    std::cout << "==========================================\n\n";
    
    RealDetect app;
    
    // Set picker type
    if (picker_type == "stalta") {
        app.setPickerType(PickerType::STALTA);
        std::cout << "Using STA/LTA picker" << std::endl;
    } else if (picker_type == "aic") {
        app.setPickerType(PickerType::AIC);
        std::cout << "Using AIC picker" << std::endl;
    } else if (picker_type == "ml") {
        app.setPickerType(PickerType::ML, ml_model_file);
        std::cout << "Using ML picker";
        if (!ml_model_file.empty()) {
            std::cout << " with model: " << ml_model_file;
        }
        std::cout << std::endl;
    } else {
        std::cerr << "Unknown picker type: " << picker_type << std::endl;
        return 1;
    }
    
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
    
    // Load regional velocity models if directory provided
    if (!velocity_dir.empty()) {
        app.loadRegionalModels(velocity_dir);
    }
    
    // Open database if provided via command line (overrides config)
    if (!database_file.empty()) {
        app.openDatabase(database_file);
    }
    
    // Check for playback mode
    bool playback_mode = !miniseed_file.empty() || !miniseed_dir.empty();
    
    if (playback_mode) {
        app.setPlaybackSpeed(playback_speed);
        
        if (!miniseed_file.empty()) {
            if (!app.loadPlaybackFile(miniseed_file)) {
                return 1;
            }
        }
        
        if (!miniseed_dir.empty()) {
            if (!app.loadPlaybackDirectory(miniseed_dir)) {
                return 1;
            }
        }
        
        std::cout << "Playback mode enabled (speed: " << playback_speed << "x)" << std::endl;
    } else {
        // Real-time mode - connect to server
        if (!app.connectToServer(server_host, server_port)) {
            return 1;
        }
        
        // Add some default streams (example: IRIS network)
        app.addStream("IU", "ANMO", "BH?");
        app.addStream("IU", "CCM", "BH?");
        app.addStream("IU", "HRV", "BH?");
        app.addStream("IU", "SSPA", "BH?");
        app.addStream("US", "ACSO", "BH?");
    }
    
    // Start processing
    if (!app.start()) {
        return 1;
    }
    
    app.stop();
    
    std::cout << "RealDetect shutdown complete" << std::endl;
    return 0;
}
