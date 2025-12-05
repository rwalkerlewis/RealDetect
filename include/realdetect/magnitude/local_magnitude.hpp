#pragma once

#include "magnitude.hpp"
#include "../picker/filter_bank.hpp"

namespace realdetect {

/**
 * LocalMagnitude - ML (Richter local magnitude)
 * 
 * Original Richter formula:
 *   ML = log10(A) + 2.76*log10(D) - 2.48
 * 
 * Where:
 *   A = maximum trace amplitude on Wood-Anderson seismograph (mm)
 *   D = epicentral distance (km)
 * 
 * Modern version uses:
 *   ML = log10(A) + a*log10(D/100) + b*(D-100) + c
 * 
 * With station corrections and regional calibration.
 */
class LocalMagnitude : public BaseMagnitude {
public:
    LocalMagnitude();
    
    MagnitudeResult calculate(
        const Origin& origin,
        const std::map<StreamID, WaveformPtr>& waveforms,
        const StationInventory& stations) override;
    
    MagnitudeType type() const override { return MagnitudeType::ML; }
    std::string name() const override { return "Local Magnitude (ML)"; }
    
    void setParameter(const std::string& name, double value) override;
    
    // Attenuation coefficients (regional calibration)
    void setAttenuationA(double a) { a_ = a; }  // log distance term
    void setAttenuationB(double b) { b_ = b; }  // linear distance term
    void setAttenuationC(double c) { c_ = c; }  // constant
    
    // Distance limits
    void setMinDistance(double km) { min_distance_ = km; }
    void setMaxDistance(double km) { max_distance_ = km; }
    
    // Wood-Anderson simulation filter
    void setSimulateWoodAnderson(bool simulate) { simulate_wa_ = simulate; }
    
    // Station corrections
    void setStationCorrection(const std::string& station, double correction) {
        station_corrections_[station] = correction;
    }

private:
    // Attenuation model coefficients
    double a_, b_, c_;
    
    // Distance limits
    double min_distance_;
    double max_distance_;
    
    // Wood-Anderson simulation
    bool simulate_wa_;
    IIRFilter wa_filter_;
    
    // Station corrections
    std::map<std::string, double> station_corrections_;
    
    // Measurement window parameters
    double pre_p_window_;    // seconds before P
    double s_window_;        // seconds after S for measurement
    
    // Calculate single station magnitude
    double calculateStationMagnitude(
        const WaveformPtr& waveform,
        const Origin& origin,
        const Station& station,
        double& amplitude);
    
    // Simulate Wood-Anderson response
    WaveformPtr simulateWoodAnderson(const WaveformPtr& waveform) const;
    
    // Find maximum amplitude in window
    double findMaxAmplitude(const WaveformPtr& waveform,
                             TimePoint start, TimePoint end) const;
};

/**
 * DurationMagnitude - Md (coda duration magnitude)
 * 
 * Based on signal duration:
 *   Md = a + b*log10(D) + c*log10(tau)
 * 
 * Where:
 *   D = epicentral distance
 *   tau = coda duration (time from P to end of coda)
 */
class DurationMagnitude : public BaseMagnitude {
public:
    DurationMagnitude();
    
    MagnitudeResult calculate(
        const Origin& origin,
        const std::map<StreamID, WaveformPtr>& waveforms,
        const StationInventory& stations) override;
    
    MagnitudeType type() const override { return MagnitudeType::Md; }
    std::string name() const override { return "Duration Magnitude (Md)"; }
    
    void setParameter(const std::string& name, double value) override;
    
    // Calibration coefficients
    void setCoefficients(double a, double b, double c) {
        a_ = a; b_ = b; c_ = c;
    }
    
    // Coda threshold (factor above noise)
    void setCodaThreshold(double factor) { coda_threshold_ = factor; }

private:
    double a_, b_, c_;
    double coda_threshold_;
    double noise_window_;  // Window for noise estimation
    
    // Find coda duration
    double findCodaDuration(const WaveformPtr& waveform,
                             TimePoint p_time) const;
};

} // namespace realdetect
