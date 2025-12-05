#pragma once

#include "picker.hpp"
#include <memory>
#include <string>
#include <vector>

namespace realdetect {

/**
 * MLModelBackend - Abstract interface for ML model inference backends
 * 
 * This allows different ML frameworks to be used (ONNX, TensorFlow, PyTorch, etc.)
 */
class MLModelBackend {
public:
    virtual ~MLModelBackend() = default;
    
    // Load model from file
    virtual bool load(const std::string& model_path) = 0;
    
    // Run inference on waveform data
    // Input: samples (continuous waveform data)
    // Output: vector of (time_offset, phase_type, probability) tuples
    virtual std::vector<std::tuple<double, PhaseType, double>> 
        predict(const SampleVector& samples, double sample_rate) = 0;
    
    // Get model info
    virtual std::string name() const = 0;
    virtual std::string version() const = 0;
    
    // Check if model is loaded
    virtual bool isLoaded() const = 0;
};

/**
 * BuiltinMLBackend - Simple built-in ML detector using characteristic functions
 * 
 * This provides a fallback ML-style detector when no external model is loaded.
 * It uses kurtosis and envelope-based detection with learned thresholds.
 */
class BuiltinMLBackend : public MLModelBackend {
public:
    BuiltinMLBackend();
    
    bool load(const std::string& model_path) override;
    
    std::vector<std::tuple<double, PhaseType, double>> 
        predict(const SampleVector& samples, double sample_rate) override;
    
    std::string name() const override { return "BuiltinML"; }
    std::string version() const override { return "1.0"; }
    bool isLoaded() const override { return true; }  // Always available

private:
    // Characteristic function parameters (can be loaded from file)
    double kurtosis_threshold_ = 5.0;
    double envelope_threshold_ = 3.0;
    double min_detection_gap_ = 0.5;  // seconds
    
    // Compute kurtosis characteristic function
    SampleVector computeKurtosis(const SampleVector& data, size_t window_samples) const;
    
    // Compute envelope characteristic function
    SampleVector computeEnvelope(const SampleVector& data) const;
};

#ifdef REALDETECT_USE_ONNX
/**
 * ONNXBackend - ONNX Runtime backend for ML models
 * 
 * Supports models like PhaseNet, EQTransformer, etc. exported to ONNX format.
 */
class ONNXBackend : public MLModelBackend {
public:
    ONNXBackend();
    ~ONNXBackend() override;
    
    bool load(const std::string& model_path) override;
    
    std::vector<std::tuple<double, PhaseType, double>> 
        predict(const SampleVector& samples, double sample_rate) override;
    
    std::string name() const override { return model_name_; }
    std::string version() const override { return model_version_; }
    bool isLoaded() const override { return session_ != nullptr; }

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
    void* session_ = nullptr;
    std::string model_name_;
    std::string model_version_;
};
#endif

/**
 * MLPicker - Machine Learning based phase picker
 * 
 * Uses neural network models to detect P and S wave arrivals.
 * Supports various backends:
 * - BuiltinML: Simple characteristic function based detector (default)
 * - ONNX: ONNX Runtime for models like PhaseNet, EQTransformer
 * - (Future: TensorFlow Lite, PyTorch, etc.)
 * 
 * Typical ML picking workflow:
 * 1. Preprocess waveform (normalize, filter)
 * 2. Run model inference to get phase probabilities over time
 * 3. Find peaks in probability functions
 * 4. Return picks at peak locations
 */
class MLPicker : public BasePicker {
public:
    MLPicker();
    ~MLPicker() override;
    
    std::vector<PickResult> pick(const Waveform& waveform) override;
    std::string name() const override;
    
    void setParameter(const std::string& name, double value) override;
    double getParameter(const std::string& name) const override;
    
    // Load ML model from file
    // Automatically detects model type from extension:
    //   .onnx -> ONNX backend
    //   .pb -> TensorFlow backend (if available)
    //   .pt -> PyTorch backend (if available)
    //   .json -> Built-in model configuration
    bool loadModel(const std::string& model_path);
    
    // Set probability threshold for detections
    void setProbabilityThreshold(double threshold) { prob_threshold_ = threshold; }
    double probabilityThreshold() const { return prob_threshold_; }
    
    // Set minimum gap between detections (seconds)
    void setMinimumGap(double seconds) { min_gap_ = seconds; }
    double minimumGap() const { return min_gap_; }
    
    // Get model information
    std::string modelName() const;
    std::string modelVersion() const;
    bool isModelLoaded() const;
    
    // Set backend explicitly
    void setBackend(std::shared_ptr<MLModelBackend> backend);

private:
    std::shared_ptr<MLModelBackend> backend_;
    double prob_threshold_;
    double min_gap_;
    
    // Preprocessing
    bool normalize_;
    double filter_low_;
    double filter_high_;
    
    // Preprocess waveform for ML inference
    SampleVector preprocess(const Waveform& waveform) const;
    
    // Post-process model output to find picks
    std::vector<PickResult> findPeaks(
        const std::vector<std::tuple<double, PhaseType, double>>& predictions,
        const Waveform& waveform) const;
};

} // namespace realdetect
