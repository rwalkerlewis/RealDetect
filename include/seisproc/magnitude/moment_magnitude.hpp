#pragma once

#include "magnitude.hpp"
#include "../picker/filter_bank.hpp"

namespace seisproc {

/**
 * MomentMagnitude - Mw from spectral analysis
 * 
 * Relationship:
 *   Mw = (2/3) * log10(M0) - 10.7
 * 
 * Where M0 is the seismic moment in dyne-cm.
 * 
 * M0 is estimated from the low-frequency spectral level:
 *   M0 = (4 * pi * rho * c^3 * R * Omega_0) / (F * R_theta_phi)
 * 
 * Where:
 *   rho = density
 *   c = wave velocity
 *   R = hypocentral distance
 *   Omega_0 = low-frequency spectral level
 *   F = free surface amplification (typically 2)
 *   R_theta_phi = radiation pattern
 */
class MomentMagnitude : public BaseMagnitude {
public:
    MomentMagnitude();
    
    MagnitudeResult calculate(
        const Origin& origin,
        const std::map<StreamID, WaveformPtr>& waveforms,
        const StationInventory& stations) override;
    
    MagnitudeType type() const override { return MagnitudeType::Mw; }
    std::string name() const override { return "Moment Magnitude (Mw)"; }
    
    void setParameter(const std::string& name, double value) override;
    
    // Medium parameters
    void setDensity(double rho) { density_ = rho; }  // kg/m³
    void setVelocity(double v) { velocity_ = v; }    // m/s
    void setQuality(double q) { quality_factor_ = q; } // Q factor
    
    // Analysis parameters
    void setFrequencyRange(double f_min, double f_max) {
        freq_min_ = f_min;
        freq_max_ = f_max;
    }
    void setSpectrumWindow(double seconds) { spectrum_window_ = seconds; }

private:
    // Medium parameters
    double density_;        // kg/m³
    double velocity_;       // m/s
    double quality_factor_; // Attenuation Q
    
    // Frequency range
    double freq_min_;
    double freq_max_;
    double spectrum_window_;
    
    // Distance limits
    double min_distance_;
    double max_distance_;
    
    // Spectral analysis result
    struct SpectralFit {
        double omega_0;      // Low-frequency level
        double corner_freq;  // Corner frequency
        double q_factor;     // Fitted Q
        double misfit;       // Fit quality
    };
    
    // Calculate spectrum and fit Brune model
    SpectralFit fitSpectrum(const WaveformPtr& waveform,
                             TimePoint start, double duration,
                             double distance) const;
    
    // Calculate single station Mw
    double calculateStationMw(const WaveformPtr& waveform,
                               const Origin& origin,
                               const Station& station,
                               SpectralFit& fit);
    
    // Brune spectrum model
    double bruneSpectrum(double freq, double omega_0, 
                          double corner_freq) const;
    
    // Attenuation correction
    double attenuationCorrection(double freq, double distance) const;
};

/**
 * BodyWaveMagnitude - mb from teleseismic P-waves
 */
class BodyWaveMagnitude : public BaseMagnitude {
public:
    BodyWaveMagnitude();
    
    MagnitudeResult calculate(
        const Origin& origin,
        const std::map<StreamID, WaveformPtr>& waveforms,
        const StationInventory& stations) override;
    
    MagnitudeType type() const override { return MagnitudeType::Mb; }
    std::string name() const override { return "Body Wave Magnitude (mb)"; }
    
    void setParameter(const std::string& name, double value) override;
    
    // mb formula: mb = log10(A/T) + Q(D,h)
    // where Q is the attenuation function
    
private:
    double min_distance_;  // Typically 5 degrees
    double max_distance_;  // Typically 100 degrees
    double period_;        // Measurement period (typically 1s)
    
    // Q function (Gutenberg-Richter attenuation)
    double attenuation(double distance_deg, double depth) const;
    
    double calculateStationMb(const WaveformPtr& waveform,
                               const Origin& origin,
                               const Station& station);
};

/**
 * SurfaceWaveMagnitude - Ms from Rayleigh waves
 */
class SurfaceWaveMagnitude : public BaseMagnitude {
public:
    SurfaceWaveMagnitude();
    
    MagnitudeResult calculate(
        const Origin& origin,
        const std::map<StreamID, WaveformPtr>& waveforms,
        const StationInventory& stations) override;
    
    MagnitudeType type() const override { return MagnitudeType::Ms; }
    std::string name() const override { return "Surface Wave Magnitude (Ms)"; }
    
    void setParameter(const std::string& name, double value) override;
    
    // Ms formula: Ms = log10(A/T) + 1.66*log10(D) + 3.3
    // Measured at 20-second period
    
private:
    double min_distance_;  // Typically 20 degrees
    double max_distance_;  // Typically 160 degrees
    double target_period_; // 20 seconds
    double period_range_;  // Allowed deviation
    
    double calculateStationMs(const WaveformPtr& waveform,
                               const Origin& origin,
                               const Station& station);
};

} // namespace seisproc
