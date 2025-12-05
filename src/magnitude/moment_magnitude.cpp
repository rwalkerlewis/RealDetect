#include "realdetect/magnitude/moment_magnitude.hpp"
#include "realdetect/magnitude/local_magnitude.hpp"
#include "realdetect/picker/filter_bank.hpp"
#include <cmath>
#include <algorithm>
#include <numeric>

namespace realdetect {

MomentMagnitude::MomentMagnitude()
    : density_(2700.0)        // kg/m³
    , velocity_(3500.0)       // m/s (S-wave)
    , quality_factor_(200.0)
    , freq_min_(1.0)
    , freq_max_(20.0)
    , spectrum_window_(5.0)
    , min_distance_(10.0)
    , max_distance_(500.0)
{
}

void MomentMagnitude::setParameter(const std::string& name, double value) {
    if (name == "density") density_ = value;
    else if (name == "velocity") velocity_ = value;
    else if (name == "quality") quality_factor_ = value;
    else if (name == "freq_min") freq_min_ = value;
    else if (name == "freq_max") freq_max_ = value;
    else if (name == "spectrum_window") spectrum_window_ = value;
}

MagnitudeResult MomentMagnitude::calculate(
    const Origin& origin,
    const std::map<StreamID, WaveformPtr>& waveforms,
    const StationInventory& stations) {
    
    MagnitudeResult result;
    result.type = MagnitudeType::Mw;
    
    std::vector<double> magnitudes;
    
    for (const auto& [stream_id, waveform] : waveforms) {
        auto sta = stations.getStation(stream_id);
        if (!sta) continue;
        
        double distance = origin.location.distanceTo(sta->location());
        
        // Skip if outside distance range
        if (distance < min_distance_ || distance > max_distance_) continue;
        
        SpectralFit fit;
        double station_mw = calculateStationMw(waveform, origin, *sta, fit);
        
        if (std::isfinite(station_mw) && station_mw > 0 && station_mw < 10) {
            StationMagnitude sta_mag;
            sta_mag.stream_id = stream_id;
            sta_mag.type = MagnitudeType::Mw;
            sta_mag.value = station_mw;
            sta_mag.distance = distance;
            
            result.station_magnitudes.push_back(sta_mag);
            magnitudes.push_back(station_mw);
        }
    }
    
    if (magnitudes.empty()) {
        return result;
    }
    
    // Compute network magnitude (median)
    std::sort(magnitudes.begin(), magnitudes.end());
    result.value = magnitudes[magnitudes.size() / 2];
    result.station_count = magnitudes.size();
    
    // Uncertainty
    double mean = std::accumulate(magnitudes.begin(), magnitudes.end(), 0.0) /
                  magnitudes.size();
    double sq_sum = 0;
    for (double m : magnitudes) {
        sq_sum += (m - mean) * (m - mean);
    }
    result.uncertainty = std::sqrt(sq_sum / magnitudes.size());
    
    return result;
}

double MomentMagnitude::calculateStationMw(const WaveformPtr& waveform,
                                             const Origin& origin,
                                             const Station& station,
                                             SpectralFit& fit) {
    double distance = origin.location.distanceTo(station.location());
    double hypo_dist = std::sqrt(distance * distance + 
                                  origin.location.depth * origin.location.depth);
    
    // Estimate S arrival time
    double vs = velocity_ / 1000.0;  // km/s
    double s_travel_time = hypo_dist / vs;
    TimePoint s_time = origin.time + 
        std::chrono::milliseconds(static_cast<int>(s_travel_time * 1000));
    
    // Fit spectrum
    fit = fitSpectrum(waveform, s_time, spectrum_window_, hypo_dist * 1000);
    
    if (fit.omega_0 <= 0) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    
    // Calculate seismic moment
    // M0 = (4 * pi * rho * c³ * R * Omega_0) / (F * R_theta_phi)
    // F = free surface amplification = 2
    // R_theta_phi = radiation pattern ≈ 0.6 (average)
    
    double R = hypo_dist * 1000;  // meters
    double c = velocity_;         // m/s
    
    double F = 2.0;              // Free surface
    double radiation = 0.6;      // Average radiation pattern
    
    double M0 = (4.0 * M_PI * density_ * c * c * c * R * fit.omega_0) / 
                (F * radiation);
    
    // Convert to Mw
    // Mw = (2/3) * log10(M0) - 10.7  (M0 in dyne-cm)
    // Our M0 is in N-m, so we need to convert: 1 N-m = 10^7 dyne-cm
    double M0_dyne_cm = M0 * 1e7;
    double mw = (2.0/3.0) * std::log10(M0_dyne_cm) - 10.7;
    
    return mw;
}

MomentMagnitude::SpectralFit MomentMagnitude::fitSpectrum(
    const WaveformPtr& waveform,
    TimePoint start, double duration,
    double distance) const {
    
    SpectralFit result;
    result.omega_0 = 0;
    result.corner_freq = 0;
    result.misfit = 1e10;
    
    // Get segment
    TimePoint end = start + std::chrono::milliseconds(static_cast<int>(duration * 1000));
    auto segment_wf = waveform->slice(start, end);
    auto segment = std::make_shared<Waveform>(segment_wf);
    
    if (segment->sampleCount() < 64) {
        return result;
    }
    
    // Prepare data
    SampleVector data = segment->data();
    
    // Demean and taper
    double mean = std::accumulate(data.begin(), data.end(), 0.0) / data.size();
    for (auto& s : data) s -= mean;
    
    size_t taper_len = data.size() / 20;
    for (size_t i = 0; i < taper_len; i++) {
        double w = 0.5 * (1.0 - std::cos(M_PI * i / taper_len));
        data[i] *= w;
        data[data.size() - 1 - i] *= w;
    }
    
    // Compute amplitude spectrum
    SampleVector amp_spec = FFT::amplitudeSpectrum(data);
    SampleVector freqs = FFT::frequencies(FFT::nextPow2(data.size()), 
                                           segment->sampleRate());
    
    if (amp_spec.empty()) {
        return result;
    }
    
    // Apply attenuation correction
    for (size_t i = 0; i < amp_spec.size() && i < freqs.size(); i++) {
        amp_spec[i] *= attenuationCorrection(freqs[i], distance);
    }
    
    // Find low-frequency level (omega_0)
    // Average spectrum in low-frequency range
    double sum_low = 0;
    int count_low = 0;
    
    for (size_t i = 0; i < freqs.size(); i++) {
        if (freqs[i] >= freq_min_ && freqs[i] <= freq_min_ * 2) {
            sum_low += amp_spec[i];
            count_low++;
        }
    }
    
    if (count_low > 0) {
        result.omega_0 = sum_low / count_low;
    }
    
    // Estimate corner frequency by finding where spectrum drops by factor of 2
    double half_level = result.omega_0 / std::sqrt(2);
    
    for (size_t i = 0; i < freqs.size(); i++) {
        if (freqs[i] > freq_min_ * 2 && amp_spec[i] < half_level) {
            result.corner_freq = freqs[i];
            break;
        }
    }
    
    if (result.corner_freq == 0) {
        result.corner_freq = (freq_min_ + freq_max_) / 2;
    }
    
    // Compute misfit to Brune model
    double misfit = 0;
    int count = 0;
    
    for (size_t i = 0; i < freqs.size(); i++) {
        if (freqs[i] >= freq_min_ && freqs[i] <= freq_max_) {
            double model = bruneSpectrum(freqs[i], result.omega_0, result.corner_freq);
            double diff = std::log10(amp_spec[i] + 1e-10) - std::log10(model + 1e-10);
            misfit += diff * diff;
            count++;
        }
    }
    
    result.misfit = count > 0 ? std::sqrt(misfit / count) : 1e10;
    
    return result;
}

double MomentMagnitude::bruneSpectrum(double freq, double omega_0, 
                                        double corner_freq) const {
    // Brune (1970) source spectrum
    // Omega(f) = Omega_0 / (1 + (f/fc)^2)
    double ratio = freq / corner_freq;
    return omega_0 / (1.0 + ratio * ratio);
}

double MomentMagnitude::attenuationCorrection(double freq, double distance) const {
    // Q attenuation correction
    // A(f) = exp(pi * f * t / Q) where t = distance/velocity
    double travel_time = distance / velocity_;
    return std::exp(M_PI * freq * travel_time / quality_factor_);
}

// BodyWaveMagnitude implementation

BodyWaveMagnitude::BodyWaveMagnitude()
    : min_distance_(555.0)    // ~5 degrees
    , max_distance_(11100.0)  // ~100 degrees
    , period_(1.0)
{
}

void BodyWaveMagnitude::setParameter(const std::string& name, double value) {
    if (name == "min_distance") min_distance_ = value;
    else if (name == "max_distance") max_distance_ = value;
    else if (name == "period") period_ = value;
}

MagnitudeResult BodyWaveMagnitude::calculate(
    const Origin& origin,
    const std::map<StreamID, WaveformPtr>& waveforms,
    const StationInventory& stations) {
    
    MagnitudeResult result;
    result.type = MagnitudeType::Mb;
    
    std::vector<double> magnitudes;
    
    for (const auto& [stream_id, waveform] : waveforms) {
        // mb is measured on vertical component
        if (stream_id.channel.back() != 'Z') continue;
        
        auto sta = stations.getStation(stream_id);
        if (!sta) continue;
        
        double distance = origin.location.distanceTo(sta->location());
        
        if (distance < min_distance_ || distance > max_distance_) continue;
        
        double mb = calculateStationMb(waveform, origin, *sta);
        
        if (std::isfinite(mb)) {
            StationMagnitude sta_mag;
            sta_mag.stream_id = stream_id;
            sta_mag.type = MagnitudeType::Mb;
            sta_mag.value = mb;
            sta_mag.distance = distance;
            
            result.station_magnitudes.push_back(sta_mag);
            magnitudes.push_back(mb);
        }
    }
    
    if (!magnitudes.empty()) {
        std::sort(magnitudes.begin(), magnitudes.end());
        result.value = magnitudes[magnitudes.size() / 2];
        result.station_count = magnitudes.size();
    }
    
    return result;
}

double BodyWaveMagnitude::attenuation(double distance_deg, double depth) const {
    // Gutenberg-Richter attenuation function
    // Simplified approximation
    double Q = 0.0;
    
    if (distance_deg < 20) {
        Q = 5.9 + 0.05 * distance_deg;
    } else if (distance_deg < 100) {
        Q = 6.9 + 0.0015 * distance_deg;
    } else {
        Q = 7.0;
    }
    
    // Depth correction
    Q -= 0.005 * (depth - 33);
    
    return Q;
}

double BodyWaveMagnitude::calculateStationMb(const WaveformPtr& waveform,
                                               const Origin& origin,
                                               const Station& station) {
    double distance = origin.location.distanceTo(station.location());
    double distance_deg = distance / 111.195;
    
    // Estimate P arrival
    double p_velocity = 8.0;  // km/s average
    double p_time = distance / p_velocity;
    TimePoint p_arrival = origin.time + 
        std::chrono::milliseconds(static_cast<int>(p_time * 1000));
    
    // Get P-wave window (first few seconds)
    TimePoint end = p_arrival + std::chrono::seconds(5);
    auto segment_wf = waveform->slice(p_arrival, end);
    auto segment = std::make_shared<Waveform>(segment_wf);
    
    if (segment->sampleCount() < 10) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    
    // Find maximum amplitude and dominant period
    const SampleVector& data = segment->data();
    double max_amp = 0;
    
    for (const auto& s : data) {
        max_amp = std::max(max_amp, std::abs(s));
    }
    
    // Convert to microns (assuming nm input)
    double amp_microns = max_amp / 1000.0;
    
    // Q value
    double Q = attenuation(distance_deg, origin.location.depth);
    
    // mb = log10(A/T) + Q
    double mb = std::log10(amp_microns / period_) + Q;
    
    return mb;
}

// SurfaceWaveMagnitude implementation

SurfaceWaveMagnitude::SurfaceWaveMagnitude()
    : min_distance_(2200.0)   // ~20 degrees
    , max_distance_(17800.0)  // ~160 degrees
    , target_period_(20.0)
    , period_range_(5.0)
{
}

void SurfaceWaveMagnitude::setParameter(const std::string& name, double value) {
    if (name == "min_distance") min_distance_ = value;
    else if (name == "max_distance") max_distance_ = value;
    else if (name == "target_period") target_period_ = value;
}

MagnitudeResult SurfaceWaveMagnitude::calculate(
    const Origin& origin,
    const std::map<StreamID, WaveformPtr>& waveforms,
    const StationInventory& stations) {
    
    MagnitudeResult result;
    result.type = MagnitudeType::Ms;
    
    std::vector<double> magnitudes;
    
    for (const auto& [stream_id, waveform] : waveforms) {
        // Ms is measured on vertical component
        if (stream_id.channel.back() != 'Z') continue;
        
        auto sta = stations.getStation(stream_id);
        if (!sta) continue;
        
        double distance = origin.location.distanceTo(sta->location());
        
        if (distance < min_distance_ || distance > max_distance_) continue;
        
        double ms = calculateStationMs(waveform, origin, *sta);
        
        if (std::isfinite(ms)) {
            StationMagnitude sta_mag;
            sta_mag.stream_id = stream_id;
            sta_mag.type = MagnitudeType::Ms;
            sta_mag.value = ms;
            sta_mag.distance = distance;
            
            result.station_magnitudes.push_back(sta_mag);
            magnitudes.push_back(ms);
        }
    }
    
    if (!magnitudes.empty()) {
        std::sort(magnitudes.begin(), magnitudes.end());
        result.value = magnitudes[magnitudes.size() / 2];
        result.station_count = magnitudes.size();
    }
    
    return result;
}

double SurfaceWaveMagnitude::calculateStationMs(const WaveformPtr& waveform,
                                                   const Origin& origin,
                                                   const Station& station) {
    double distance = origin.location.distanceTo(station.location());
    double distance_deg = distance / 111.195;
    
    // Surface wave group velocity ~3.5 km/s
    double surface_velocity = 3.5;
    double surface_time = distance / surface_velocity;
    TimePoint surface_arrival = origin.time + 
        std::chrono::milliseconds(static_cast<int>(surface_time * 1000));
    
    // Get surface wave window (60 seconds around predicted arrival)
    TimePoint start = surface_arrival - std::chrono::seconds(30);
    TimePoint end = surface_arrival + std::chrono::seconds(60);
    auto segment_wf = waveform->slice(start, end);
    auto segment = std::make_shared<Waveform>(segment_wf);
    
    if (segment->sampleCount() < 100) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    
    // Bandpass filter around 20-second period
    double low_freq = 1.0 / (target_period_ + period_range_);
    double high_freq = 1.0 / (target_period_ - period_range_);
    
    IIRFilter bp = IIRFilter::butterworth(3, low_freq, high_freq, 
                                           segment->sampleRate());
    SampleVector filtered = bp.filtfilt(segment->data());
    
    // Find maximum amplitude
    double max_amp = 0;
    for (const auto& s : filtered) {
        max_amp = std::max(max_amp, std::abs(s));
    }
    
    // Convert to microns
    double amp_microns = max_amp / 1000.0;
    
    // Ms = log10(A/T) + 1.66*log10(D) + 3.3
    double ms = std::log10(amp_microns / target_period_) + 
                1.66 * std::log10(distance_deg) + 3.3;
    
    return ms;
}

// MagnitudeFactory

MagnitudeCalculatorPtr MagnitudeFactory::create(MagnitudeType type) {
    switch (type) {
        case MagnitudeType::ML:
            return std::make_shared<LocalMagnitude>();
        case MagnitudeType::Mw:
            return std::make_shared<MomentMagnitude>();
        case MagnitudeType::Mb:
            return std::make_shared<BodyWaveMagnitude>();
        case MagnitudeType::Ms:
            return std::make_shared<SurfaceWaveMagnitude>();
        case MagnitudeType::Md:
            return std::make_shared<DurationMagnitude>();
        default:
            return std::make_shared<LocalMagnitude>();
    }
}

std::vector<MagnitudeType> MagnitudeFactory::availableTypes() {
    return {MagnitudeType::ML, MagnitudeType::Mw, MagnitudeType::Mb, 
            MagnitudeType::Ms, MagnitudeType::Md};
}

} // namespace realdetect
