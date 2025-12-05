#include "realdetect/magnitude/local_magnitude.hpp"
#include <cmath>
#include <algorithm>
#include <numeric>

namespace realdetect {

LocalMagnitude::LocalMagnitude()
    : a_(1.0)
    , b_(0.00301)
    , c_(3.0)
    , min_distance_(5.0)
    , max_distance_(600.0)
    , simulate_wa_(true)
    , pre_p_window_(5.0)
    , s_window_(10.0)
{
    // Initialize Wood-Anderson simulation filter
    // WA: period=0.8s, damping=0.8, gain=2800
    // We approximate with a 2nd order bandpass filter
    wa_filter_ = IIRFilter::butterworth(2, 0.5, 10.0, 100.0);
}

void LocalMagnitude::setParameter(const std::string& name, double value) {
    if (name == "a") a_ = value;
    else if (name == "b") b_ = value;
    else if (name == "c") c_ = value;
    else if (name == "min_distance") min_distance_ = value;
    else if (name == "max_distance") max_distance_ = value;
    else if (name == "simulate_wa") simulate_wa_ = (value != 0);
}

MagnitudeResult LocalMagnitude::calculate(
    const Origin& origin,
    const std::map<StreamID, WaveformPtr>& waveforms,
    const StationInventory& stations) {
    
    MagnitudeResult result;
    result.type = MagnitudeType::ML;
    
    std::vector<double> magnitudes;
    
    for (const auto& [stream_id, waveform] : waveforms) {
        auto sta = stations.getStation(stream_id);
        if (!sta) continue;
        
        double distance = origin.location.distanceTo(sta->location());
        
        // Skip if outside distance range
        if (distance < min_distance_ || distance > max_distance_) continue;
        
        // Calculate station magnitude
        double amplitude;
        double station_ml = calculateStationMagnitude(waveform, origin, *sta, amplitude);
        
        if (std::isfinite(station_ml)) {
            StationMagnitude sta_mag;
            sta_mag.stream_id = stream_id;
            sta_mag.type = MagnitudeType::ML;
            sta_mag.value = station_ml;
            sta_mag.amplitude = amplitude;
            sta_mag.distance = distance;
            
            // Apply station correction if available
            auto corr_it = station_corrections_.find(sta->code());
            if (corr_it != station_corrections_.end()) {
                sta_mag.correction = corr_it->second;
                sta_mag.value += corr_it->second;
            }
            
            result.station_magnitudes.push_back(sta_mag);
            magnitudes.push_back(sta_mag.value);
        }
    }
    
    if (magnitudes.empty()) {
        result.value = 0;
        result.uncertainty = 0;
        result.station_count = 0;
        return result;
    }
    
    // Compute network magnitude (median is more robust)
    std::sort(magnitudes.begin(), magnitudes.end());
    result.value = magnitudes[magnitudes.size() / 2];  // Median
    result.station_count = magnitudes.size();
    
    // Compute uncertainty (standard deviation)
    double mean = std::accumulate(magnitudes.begin(), magnitudes.end(), 0.0) / 
                  magnitudes.size();
    double sq_sum = 0;
    for (double m : magnitudes) {
        sq_sum += (m - mean) * (m - mean);
    }
    result.uncertainty = std::sqrt(sq_sum / magnitudes.size());
    
    return result;
}

double LocalMagnitude::calculateStationMagnitude(
    const WaveformPtr& waveform,
    const Origin& origin,
    const Station& station,
    double& amplitude) {
    
    amplitude = 0;
    
    double distance = origin.location.distanceTo(station.location());
    
    // Find P and S arrival times
    // Approximate using velocity model
    double vp = 6.0;  // km/s
    double vs = 3.5;  // km/s
    double hypo_dist = std::sqrt(distance * distance + 
                                  origin.location.depth * origin.location.depth);
    
    double p_travel_time = hypo_dist / vp;
    double s_travel_time = hypo_dist / vs;
    
    TimePoint p_time = origin.time + 
        std::chrono::milliseconds(static_cast<int>(p_travel_time * 1000));
    TimePoint s_time = origin.time + 
        std::chrono::milliseconds(static_cast<int>(s_travel_time * 1000));
    
    // Measurement window: from just before S to s_window_ after S
    TimePoint start = s_time - std::chrono::milliseconds(500);
    TimePoint end = s_time + 
        std::chrono::milliseconds(static_cast<int>(s_window_ * 1000));
    
    // Get waveform slice
    auto segment_wf = waveform->slice(start, end);
    auto segment = std::make_shared<Waveform>(segment_wf);
    if (segment->sampleCount() < 10) return std::numeric_limits<double>::quiet_NaN();
    
    // Optionally simulate Wood-Anderson
    if (simulate_wa_) {
        segment = simulateWoodAnderson(segment);
    }
    
    // Find maximum amplitude
    amplitude = findMaxAmplitude(segment, segment->startTime(), segment->endTime());
    
    if (amplitude <= 0) return std::numeric_limits<double>::quiet_NaN();
    
    // Convert to mm if needed (assuming counts, need instrument response)
    // For now, assume 1 nm/count
    double amp_mm = amplitude * 1e-6;  // Convert nm to mm
    
    // Compute magnitude using attenuation formula
    // ML = log10(A) + a*log10(D/100) + b*(D-100) + c
    double ml = std::log10(amp_mm) + 
                a_ * std::log10(distance / 100.0) + 
                b_ * (distance - 100.0) + 
                c_;
    
    return ml;
}

WaveformPtr LocalMagnitude::simulateWoodAnderson(const WaveformPtr& waveform) const {
    // Wood-Anderson simulation
    // WA characteristics: T0=0.8s, h=0.8, V=2800
    // We approximate with bandpass filter and gain
    
    auto result = std::make_shared<Waveform>(*waveform);
    
    // Apply filter
    result->demean();
    result->detrend();
    
    SampleVector filtered = wa_filter_.filtfilt(result->data());
    result->data() = filtered;
    
    // Apply approximate WA gain (2800)
    for (auto& s : result->data()) {
        s *= 2800.0;
    }
    
    return result;
}

double LocalMagnitude::findMaxAmplitude(const WaveformPtr& waveform,
                                          TimePoint start, TimePoint end) const {
    const SampleVector& data = waveform->data();
    if (data.empty()) return 0;
    
    // Find absolute maximum in window
    double max_amp = 0;
    for (const auto& s : data) {
        max_amp = std::max(max_amp, std::abs(s));
    }
    
    return max_amp;
}

// DurationMagnitude implementation

DurationMagnitude::DurationMagnitude()
    : a_(-1.0)
    , b_(2.0)
    , c_(0.0035)
    , coda_threshold_(1.5)
    , noise_window_(5.0)
{
}

void DurationMagnitude::setParameter(const std::string& name, double value) {
    if (name == "a") a_ = value;
    else if (name == "b") b_ = value;
    else if (name == "c") c_ = value;
    else if (name == "coda_threshold") coda_threshold_ = value;
}

MagnitudeResult DurationMagnitude::calculate(
    const Origin& origin,
    const std::map<StreamID, WaveformPtr>& waveforms,
    const StationInventory& stations) {
    
    MagnitudeResult result;
    result.type = MagnitudeType::Md;
    
    std::vector<double> magnitudes;
    
    for (const auto& [stream_id, waveform] : waveforms) {
        auto sta = stations.getStation(stream_id);
        if (!sta) continue;
        
        double distance = origin.location.distanceTo(sta->location());
        
        // Estimate P time
        double vp = 6.0;
        double hypo_dist = std::sqrt(distance * distance + 
                                      origin.location.depth * origin.location.depth);
        double p_travel_time = hypo_dist / vp;
        TimePoint p_time = origin.time + 
            std::chrono::milliseconds(static_cast<int>(p_travel_time * 1000));
        
        // Find coda duration
        double duration = findCodaDuration(waveform, p_time);
        
        if (duration > 0) {
            // Md = a + b*log10(tau) + c*D
            double md = a_ + b_ * std::log10(duration) + c_ * distance;
            
            StationMagnitude sta_mag;
            sta_mag.stream_id = stream_id;
            sta_mag.type = MagnitudeType::Md;
            sta_mag.value = md;
            sta_mag.distance = distance;
            
            result.station_magnitudes.push_back(sta_mag);
            magnitudes.push_back(md);
        }
    }
    
    if (magnitudes.empty()) {
        return result;
    }
    
    // Network magnitude (median)
    std::sort(magnitudes.begin(), magnitudes.end());
    result.value = magnitudes[magnitudes.size() / 2];
    result.station_count = magnitudes.size();
    
    return result;
}

double DurationMagnitude::findCodaDuration(const WaveformPtr& waveform,
                                             TimePoint p_time) const {
    // Get waveform data
    const SampleVector& data = waveform->data();
    if (data.size() < 100) return 0;
    
    // Estimate noise level from pre-P window
    int64_t p_idx = waveform->indexAt(p_time);
    int noise_samples = static_cast<int>(noise_window_ * waveform->sampleRate());
    int noise_start = std::max(int64_t(0), p_idx - noise_samples);
    
    double noise_level = 0;
    for (int i = noise_start; i < p_idx && i < static_cast<int>(data.size()); i++) {
        noise_level += data[i] * data[i];
    }
    noise_level = std::sqrt(noise_level / std::max(1, static_cast<int>(p_idx - noise_start)));
    
    // Find where signal drops below threshold * noise
    double threshold = noise_level * coda_threshold_;
    
    // Use sliding window envelope
    int window = static_cast<int>(waveform->sampleRate());  // 1 second
    int end_idx = p_idx;
    
    for (int i = p_idx; i < static_cast<int>(data.size()) - window; i += window/2) {
        double env = 0;
        for (int j = i; j < i + window; j++) {
            env = std::max(env, std::abs(data[j]));
        }
        
        if (env < threshold) {
            end_idx = i;
            break;
        }
        end_idx = i;
    }
    
    // Duration in seconds
    return (end_idx - p_idx) / waveform->sampleRate();
}

} // namespace realdetect
