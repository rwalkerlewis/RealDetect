#pragma once

#include "types.hpp"
#include "station.hpp"
#include <vector>
#include <memory>
#include <optional>

namespace realdetect {

/**
 * Pick - A phase arrival pick
 */
struct Pick {
    StreamID stream_id;
    TimePoint time;
    PhaseType phase_type;
    PickQuality quality;
    
    double amplitude;       // Peak amplitude
    double period;          // Dominant period (seconds)
    double snr;             // Signal-to-noise ratio
    double uncertainty;     // Time uncertainty in seconds
    
    bool is_automatic;
    std::string method;     // Picking method used
    
    Pick() : phase_type(PhaseType::Unknown), quality(PickQuality::Questionable),
             amplitude(0), period(0), snr(0), uncertainty(0), is_automatic(true) {}
    
    // Weight from quality (0=best, 4=worst)
    double weight() const {
        switch (quality) {
            case PickQuality::Impulsive: return 0.0;
            case PickQuality::Emergent: return 1.0;
            case PickQuality::Questionable: return 2.0;
            case PickQuality::Poor: return 3.0;
            case PickQuality::Rejected: return 4.0;
        }
        return 4.0;
    }
    
    // Residual weight for location
    double locationWeight() const {
        double w = 1.0 - weight() / 4.0;
        return w * w;  // Square for more weight differentiation
    }
};

using PickPtr = std::shared_ptr<Pick>;

/**
 * Arrival - Pick associated with an event
 */
struct Arrival {
    PickPtr pick;
    double distance;        // Epicentral distance (km)
    double azimuth;         // Azimuth from event to station
    double takeoff_angle;   // Takeoff angle at source
    double residual;        // Observed - calculated time (seconds)
    double weight;          // Weight in location
    bool used;              // Used in location
    
    Arrival() : distance(0), azimuth(0), takeoff_angle(0),
                residual(0), weight(1.0), used(true) {}
};

/**
 * Magnitude - Single magnitude measurement
 */
struct Magnitude {
    MagnitudeType type;
    double value;
    double uncertainty;
    int station_count;
    std::string method;
    
    Magnitude() : type(MagnitudeType::Unknown), value(0), uncertainty(0), station_count(0) {}
    Magnitude(MagnitudeType t, double v, double u = 0, int cnt = 0)
        : type(t), value(v), uncertainty(u), station_count(cnt) {}
};

/**
 * StationMagnitude - Per-station magnitude
 */
struct StationMagnitude {
    StreamID stream_id;
    MagnitudeType type;
    double value;
    double amplitude;
    double distance;
    double correction;
    
    StationMagnitude() : type(MagnitudeType::Unknown), value(0),
                         amplitude(0), distance(0), correction(0) {}
};

/**
 * Origin - Event hypocenter solution
 */
struct Origin {
    TimePoint time;
    GeoPoint location;
    
    double latitude_error;   // km
    double longitude_error;  // km
    double depth_error;      // km
    double time_error;       // seconds
    
    double rms;              // RMS residual
    double gap;              // Azimuthal gap (degrees)
    int phase_count;         // Number of phases used
    int station_count;       // Number of stations
    
    std::string algorithm;   // Location algorithm used
    bool is_fixed_depth;
    
    std::vector<Arrival> arrivals;
    
    Origin() : latitude_error(0), longitude_error(0), depth_error(0),
               time_error(0), rms(0), gap(360), phase_count(0),
               station_count(0), is_fixed_depth(false) {}
    
    // Quality metrics
    bool isQualityA() const { return gap < 90 && phase_count >= 8 && rms < 0.5; }
    bool isQualityB() const { return gap < 135 && phase_count >= 6 && rms < 1.0; }
    bool isQualityC() const { return gap < 180 && phase_count >= 4 && rms < 2.0; }
    bool isQualityD() const { return gap < 270 && phase_count >= 3; }
    
    char qualityCode() const {
        if (isQualityA()) return 'A';
        if (isQualityB()) return 'B';
        if (isQualityC()) return 'C';
        if (isQualityD()) return 'D';
        return 'E';
    }
};

/**
 * Event - Complete seismic event
 */
class Event {
public:
    Event() : id_(generateId()) {}
    explicit Event(const std::string& id) : id_(id) {}
    
    // Identifiers
    const std::string& id() const { return id_; }
    void setId(const std::string& id) { id_ = id; }
    
    // Origins
    void addOrigin(const Origin& origin) { origins_.push_back(origin); }
    const std::vector<Origin>& origins() const { return origins_; }
    std::vector<Origin>& origins() { return origins_; }
    
    const Origin& preferredOrigin() const {
        return origins_.empty() ? dummy_origin_ : origins_.back();
    }
    
    // Magnitudes
    void addMagnitude(const Magnitude& mag) { magnitudes_.push_back(mag); }
    const std::vector<Magnitude>& magnitudes() const { return magnitudes_; }
    
    std::optional<Magnitude> getMagnitude(MagnitudeType type) const {
        for (const auto& m : magnitudes_) {
            if (m.type == type) return m;
        }
        return std::nullopt;
    }
    
    const Magnitude& preferredMagnitude() const {
        // Preference order: Mw > ML > Mb > Ms > Md
        for (auto type : {MagnitudeType::Mw, MagnitudeType::ML, 
                          MagnitudeType::Mb, MagnitudeType::Ms, MagnitudeType::Md}) {
            for (const auto& m : magnitudes_) {
                if (m.type == type) return m;
            }
        }
        return magnitudes_.empty() ? dummy_magnitude_ : magnitudes_.front();
    }
    
    // Station magnitudes
    void addStationMagnitude(const StationMagnitude& sm) { 
        station_magnitudes_.push_back(sm); 
    }
    const std::vector<StationMagnitude>& stationMagnitudes() const { 
        return station_magnitudes_; 
    }
    
    // Convenience accessors
    TimePoint time() const { return preferredOrigin().time; }
    const GeoPoint& location() const { return preferredOrigin().location; }
    double latitude() const { return preferredOrigin().location.latitude; }
    double longitude() const { return preferredOrigin().location.longitude; }
    double depth() const { return preferredOrigin().location.depth; }
    double magnitude() const { return preferredMagnitude().value; }
    MagnitudeType magnitudeType() const { return preferredMagnitude().type; }
    
    // Summary
    std::string summary() const;

private:
    std::string id_;
    std::vector<Origin> origins_;
    std::vector<Magnitude> magnitudes_;
    std::vector<StationMagnitude> station_magnitudes_;
    
    static Origin dummy_origin_;
    static Magnitude dummy_magnitude_;
    
    static std::string generateId();
};

using EventPtr = std::shared_ptr<Event>;

} // namespace realdetect
