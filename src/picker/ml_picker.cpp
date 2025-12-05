#include "realdetect/picker/ml_picker.hpp"
#include <cmath>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <fstream>

namespace realdetect {

// ============================================================================
// BuiltinMLBackend Implementation
// ============================================================================

BuiltinMLBackend::BuiltinMLBackend() 
    : kurtosis_threshold_(5.0)
    , envelope_threshold_(3.0)
    , min_detection_gap_(0.5)
{
}

bool BuiltinMLBackend::load(const std::string& model_path) {
    // Load configuration from JSON file if provided
    if (model_path.empty()) {
        return true;  // Use defaults
    }
    
    std::ifstream file(model_path);
    if (!file.is_open()) {
        std::cerr << "BuiltinMLBackend: Could not open config file: " << model_path << std::endl;
        return false;
    }
    
    // Simple JSON-like parsing for configuration
    std::string line;
    while (std::getline(file, line)) {
        // Skip comments and empty lines
        if (line.empty() || line[0] == '#' || line[0] == '/') continue;
        
        auto pos = line.find(':');
        if (pos != std::string::npos) {
            std::string key = line.substr(0, pos);
            std::string value = line.substr(pos + 1);
            
            // Trim whitespace
            key.erase(0, key.find_first_not_of(" \t\""));
            key.erase(key.find_last_not_of(" \t\",") + 1);
            value.erase(0, value.find_first_not_of(" \t"));
            value.erase(value.find_last_not_of(" \t,") + 1);
            
            try {
                if (key == "kurtosis_threshold") {
                    kurtosis_threshold_ = std::stod(value);
                } else if (key == "envelope_threshold") {
                    envelope_threshold_ = std::stod(value);
                } else if (key == "min_detection_gap") {
                    min_detection_gap_ = std::stod(value);
                }
            } catch (...) {
                // Ignore parse errors
            }
        }
    }
    
    std::cout << "BuiltinMLBackend: Loaded configuration from " << model_path << std::endl;
    return true;
}

SampleVector BuiltinMLBackend::computeKurtosis(const SampleVector& data, size_t window_samples) const {
    SampleVector result(data.size(), 0.0);
    
    if (data.size() < window_samples) {
        return result;
    }
    
    for (size_t i = window_samples; i < data.size(); i++) {
        // Compute mean and variance
        double sum = 0, sum2 = 0;
        for (size_t j = i - window_samples; j < i; j++) {
            sum += data[j];
            sum2 += data[j] * data[j];
        }
        double mean = sum / window_samples;
        double var = sum2 / window_samples - mean * mean;
        double std = std::sqrt(std::max(var, 1e-10));
        
        // Compute kurtosis (4th moment)
        double m4 = 0;
        for (size_t j = i - window_samples; j < i; j++) {
            double z = (data[j] - mean) / std;
            m4 += z * z * z * z;
        }
        result[i] = m4 / window_samples - 3.0;  // Excess kurtosis
    }
    
    return result;
}

SampleVector BuiltinMLBackend::computeEnvelope(const SampleVector& data) const {
    // Simple envelope using absolute value and smoothing
    SampleVector envelope(data.size());
    
    // Compute absolute values
    for (size_t i = 0; i < data.size(); i++) {
        envelope[i] = std::abs(data[i]);
    }
    
    // Smooth with moving average (window of 10 samples)
    const size_t smooth_window = 10;
    SampleVector smoothed(data.size(), 0.0);
    
    double sum = 0;
    for (size_t i = 0; i < data.size(); i++) {
        sum += envelope[i];
        if (i >= smooth_window) {
            sum -= envelope[i - smooth_window];
            smoothed[i] = sum / smooth_window;
        } else {
            smoothed[i] = sum / (i + 1);
        }
    }
    
    return smoothed;
}

std::vector<std::tuple<double, PhaseType, double>> 
BuiltinMLBackend::predict(const SampleVector& samples, double sample_rate) {
    std::vector<std::tuple<double, PhaseType, double>> predictions;
    
    if (samples.empty() || sample_rate <= 0) {
        return predictions;
    }
    
    size_t window_samples = static_cast<size_t>(0.5 * sample_rate);  // 0.5s window
    
    // Compute characteristic functions
    SampleVector kurtosis = computeKurtosis(samples, window_samples);
    SampleVector envelope = computeEnvelope(samples);
    
    // Normalize envelope for threshold comparison
    double env_mean = 0, env_std = 0;
    for (size_t i = window_samples; i < envelope.size(); i++) {
        env_mean += envelope[i];
    }
    env_mean /= (envelope.size() - window_samples);
    
    for (size_t i = window_samples; i < envelope.size(); i++) {
        env_std += (envelope[i] - env_mean) * (envelope[i] - env_mean);
    }
    env_std = std::sqrt(env_std / (envelope.size() - window_samples));
    
    // Combined detection
    size_t min_gap_samples = static_cast<size_t>(min_detection_gap_ * sample_rate);
    size_t last_detection = 0;
    
    for (size_t i = window_samples; i < samples.size() - window_samples; i++) {
        // Skip if too close to last detection
        if (i < last_detection + min_gap_samples) {
            continue;
        }
        
        double kurtosis_score = kurtosis[i];
        double envelope_score = (envelope[i] - env_mean) / std::max(env_std, 1e-10);
        
        // Detect based on kurtosis (impulsive onset)
        if (kurtosis_score > kurtosis_threshold_) {
            double time_offset = static_cast<double>(i) / sample_rate;
            double probability = std::min(1.0, kurtosis_score / (kurtosis_threshold_ * 2));
            
            // P wave detection (assuming impulsive = P)
            predictions.emplace_back(time_offset, PhaseType::P, probability);
            last_detection = i;
        }
        // Detect based on envelope (amplitude increase)
        else if (envelope_score > envelope_threshold_) {
            double time_offset = static_cast<double>(i) / sample_rate;
            double probability = std::min(1.0, envelope_score / (envelope_threshold_ * 2));
            
            // Could be S wave or later arrival
            predictions.emplace_back(time_offset, PhaseType::S, probability * 0.8);
            last_detection = i;
        }
    }
    
    return predictions;
}

// ============================================================================
// MLPicker Implementation
// ============================================================================

MLPicker::MLPicker()
    : backend_(std::make_shared<BuiltinMLBackend>())
    , prob_threshold_(0.5)
    , min_gap_(0.5)
    , normalize_(true)
    , filter_low_(1.0)
    , filter_high_(20.0)
{
}

MLPicker::~MLPicker() = default;

std::string MLPicker::name() const {
    if (backend_) {
        return "ML/" + backend_->name();
    }
    return "ML";
}

bool MLPicker::loadModel(const std::string& model_path) {
    if (model_path.empty()) {
        std::cerr << "MLPicker: No model path provided, using built-in detector" << std::endl;
        backend_ = std::make_shared<BuiltinMLBackend>();
        return true;
    }
    
    // Detect model type from extension
    std::string ext;
    auto pos = model_path.rfind('.');
    if (pos != std::string::npos) {
        ext = model_path.substr(pos);
    }
    
    if (ext == ".json" || ext == ".conf" || ext == ".cfg") {
        // Built-in model with configuration
        auto builtin = std::make_shared<BuiltinMLBackend>();
        if (builtin->load(model_path)) {
            backend_ = builtin;
            std::cout << "MLPicker: Loaded built-in model configuration from " << model_path << std::endl;
            return true;
        }
        return false;
    }
    
#ifdef REALDETECT_USE_ONNX
    if (ext == ".onnx") {
        auto onnx = std::make_shared<ONNXBackend>();
        if (onnx->load(model_path)) {
            backend_ = onnx;
            std::cout << "MLPicker: Loaded ONNX model from " << model_path << std::endl;
            return true;
        }
        return false;
    }
#endif
    
    std::cerr << "MLPicker: Unsupported model format: " << ext << std::endl;
    std::cerr << "MLPicker: Falling back to built-in detector" << std::endl;
    
    // Try loading as configuration for built-in
    auto builtin = std::make_shared<BuiltinMLBackend>();
    builtin->load(model_path);
    backend_ = builtin;
    
    return true;
}

void MLPicker::setBackend(std::shared_ptr<MLModelBackend> backend) {
    if (backend) {
        backend_ = backend;
    }
}

std::string MLPicker::modelName() const {
    return backend_ ? backend_->name() : "None";
}

std::string MLPicker::modelVersion() const {
    return backend_ ? backend_->version() : "";
}

bool MLPicker::isModelLoaded() const {
    return backend_ && backend_->isLoaded();
}

void MLPicker::setParameter(const std::string& param_name, double value) {
    if (param_name == "probability_threshold" || param_name == "threshold") {
        prob_threshold_ = value;
    } else if (param_name == "min_gap" || param_name == "minimum_gap") {
        min_gap_ = value;
    } else if (param_name == "normalize") {
        normalize_ = (value != 0);
    } else if (param_name == "filter_low") {
        filter_low_ = value;
    } else if (param_name == "filter_high") {
        filter_high_ = value;
    }
}

double MLPicker::getParameter(const std::string& param_name) const {
    if (param_name == "probability_threshold" || param_name == "threshold") {
        return prob_threshold_;
    } else if (param_name == "min_gap" || param_name == "minimum_gap") {
        return min_gap_;
    } else if (param_name == "normalize") {
        return normalize_ ? 1.0 : 0.0;
    } else if (param_name == "filter_low") {
        return filter_low_;
    } else if (param_name == "filter_high") {
        return filter_high_;
    }
    return 0.0;
}

SampleVector MLPicker::preprocess(const Waveform& waveform) const {
    SampleVector data = waveform.data();
    
    if (data.empty()) {
        return data;
    }
    
    // Remove mean
    double mean = std::accumulate(data.begin(), data.end(), 0.0) / data.size();
    for (auto& sample : data) {
        sample -= mean;
    }
    
    // Normalize if enabled
    if (normalize_) {
        double max_abs = 0;
        for (const auto& sample : data) {
            max_abs = std::max(max_abs, std::abs(sample));
        }
        if (max_abs > 0) {
            for (auto& sample : data) {
                sample /= max_abs;
            }
        }
    }
    
    return data;
}

std::vector<PickResult> MLPicker::findPeaks(
    const std::vector<std::tuple<double, PhaseType, double>>& predictions,
    const Waveform& waveform) const 
{
    std::vector<PickResult> picks;
    
    TimePoint start_time = waveform.startTime();
    double sample_rate = waveform.sampleRate();
    
    for (const auto& [time_offset, phase_type, probability] : predictions) {
        if (probability < prob_threshold_) {
            continue;
        }
        
        PickResult pick;
        pick.time = start_time + std::chrono::microseconds(
            static_cast<int64_t>(time_offset * 1e6));
        pick.phase_type = phase_type;
        pick.confidence = probability;
        pick.snr = probability * 10;  // Approximate SNR from probability
        pick.sample_index = static_cast<size_t>(time_offset * sample_rate);
        
        // Estimate amplitude at pick
        if (pick.sample_index < waveform.data().size()) {
            pick.amplitude = std::abs(waveform.data()[pick.sample_index]);
        }
        
        picks.push_back(pick);
    }
    
    return picks;
}

std::vector<PickResult> MLPicker::pick(const Waveform& waveform) {
    if (!backend_ || !backend_->isLoaded()) {
        return {};
    }
    
    // Preprocess waveform
    SampleVector data = preprocess(waveform);
    
    if (data.empty()) {
        return {};
    }
    
    // Run ML model
    auto predictions = backend_->predict(data, waveform.sampleRate());
    
    // Convert predictions to picks
    return findPeaks(predictions, waveform);
}

} // namespace realdetect
