#pragma once

#include "../core/types.hpp"
#include "../core/event.hpp"
#include "../core/station.hpp"
#include "../core/waveform.hpp"
#include <vector>
#include <memory>
#include <functional>

namespace realdetect {

/**
 * MagnitudeResult - Result from magnitude calculation
 */
struct MagnitudeResult {
    MagnitudeType type;
    double value;
    double uncertainty;
    int station_count;
    std::vector<StationMagnitude> station_magnitudes;
    
    MagnitudeResult() : type(MagnitudeType::Unknown), value(0), 
                        uncertainty(0), station_count(0) {}
};

/**
 * BaseMagnitude - Abstract base class for magnitude calculators
 */
class BaseMagnitude {
public:
    virtual ~BaseMagnitude() = default;
    
    // Calculate magnitude from waveforms
    virtual MagnitudeResult calculate(
        const Origin& origin,
        const std::map<StreamID, WaveformPtr>& waveforms,
        const StationInventory& stations) = 0;
    
    // Get magnitude type
    virtual MagnitudeType type() const = 0;
    virtual std::string name() const = 0;
    
    // Configuration
    virtual void setParameter(const std::string& name, double value) = 0;
};

using MagnitudeCalculatorPtr = std::shared_ptr<BaseMagnitude>;

/**
 * MagnitudeFactory - Create magnitude calculator instances
 */
class MagnitudeFactory {
public:
    static MagnitudeCalculatorPtr create(MagnitudeType type);
    static std::vector<MagnitudeType> availableTypes();
};

/**
 * WaveformProvider - Interface for getting waveforms
 */
class WaveformProvider {
public:
    virtual ~WaveformProvider() = default;
    
    virtual WaveformPtr getWaveform(const StreamID& id, 
                                     TimePoint start, TimePoint end) = 0;
    
    virtual std::vector<StreamID> availableStreams() const = 0;
};

} // namespace realdetect
