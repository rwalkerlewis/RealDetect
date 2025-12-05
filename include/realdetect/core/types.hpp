#pragma once

#include <cstdint>
#include <string>
#include <vector>
#include <chrono>
#include <cmath>
#include <limits>

namespace realdetect {

// Time handling - using microsecond precision
using TimePoint = std::chrono::system_clock::time_point;
using Duration = std::chrono::microseconds;

// Sample data type
using Sample = double;
using SampleVector = std::vector<Sample>;

// Geographic coordinates
struct GeoPoint {
    double latitude;   // degrees, -90 to 90
    double longitude;  // degrees, -180 to 180
    double depth;      // km below surface (positive downward)
    
    GeoPoint() : latitude(0), longitude(0), depth(0) {}
    GeoPoint(double lat, double lon, double dep = 0) 
        : latitude(lat), longitude(lon), depth(dep) {}
    
    // Haversine distance in km
    double distanceTo(const GeoPoint& other) const {
        constexpr double R = 6371.0; // Earth radius in km
        double lat1 = latitude * M_PI / 180.0;
        double lat2 = other.latitude * M_PI / 180.0;
        double dLat = (other.latitude - latitude) * M_PI / 180.0;
        double dLon = (other.longitude - longitude) * M_PI / 180.0;
        
        double a = std::sin(dLat/2) * std::sin(dLat/2) +
                   std::cos(lat1) * std::cos(lat2) *
                   std::sin(dLon/2) * std::sin(dLon/2);
        double c = 2 * std::atan2(std::sqrt(a), std::sqrt(1-a));
        return R * c;
    }
    
    // Azimuth to another point in degrees
    double azimuthTo(const GeoPoint& other) const {
        double lat1 = latitude * M_PI / 180.0;
        double lat2 = other.latitude * M_PI / 180.0;
        double dLon = (other.longitude - longitude) * M_PI / 180.0;
        
        double y = std::sin(dLon) * std::cos(lat2);
        double x = std::cos(lat1) * std::sin(lat2) -
                   std::sin(lat1) * std::cos(lat2) * std::cos(dLon);
        double az = std::atan2(y, x) * 180.0 / M_PI;
        return std::fmod(az + 360.0, 360.0);
    }
};

// Stream identifier (SEED convention)
struct StreamID {
    std::string network;     // 2 char
    std::string station;     // 5 char max
    std::string location;    // 2 char
    std::string channel;     // 3 char (e.g., BHZ, HHN)
    
    StreamID() = default;
    StreamID(const std::string& net, const std::string& sta,
             const std::string& loc, const std::string& chan)
        : network(net), station(sta), location(loc), channel(chan) {}
    
    std::string toString() const {
        return network + "." + station + "." + location + "." + channel;
    }
    
    std::string seedId() const {
        return network + "_" + station + "_" + location + "_" + channel;
    }
    
    bool operator==(const StreamID& other) const {
        return network == other.network && station == other.station &&
               location == other.location && channel == other.channel;
    }
    
    bool operator<(const StreamID& other) const {
        if (network != other.network) return network < other.network;
        if (station != other.station) return station < other.station;
        if (location != other.location) return location < other.location;
        return channel < other.channel;
    }
};

// Phase types
enum class PhaseType {
    P,      // Primary/compressional
    S,      // Secondary/shear  
    Pn,     // P refracted at Moho
    Sn,     // S refracted at Moho
    Pg,     // P in crust
    Sg,     // S in crust
    Lg,     // Crustal guided wave
    Unknown
};

inline std::string phaseTypeToString(PhaseType pt) {
    switch (pt) {
        case PhaseType::P: return "P";
        case PhaseType::S: return "S";
        case PhaseType::Pn: return "Pn";
        case PhaseType::Sn: return "Sn";
        case PhaseType::Pg: return "Pg";
        case PhaseType::Sg: return "Sg";
        case PhaseType::Lg: return "Lg";
        default: return "?";
    }
}

inline PhaseType stringToPhaseType(const std::string& s) {
    if (s == "P") return PhaseType::P;
    if (s == "S") return PhaseType::S;
    if (s == "Pn") return PhaseType::Pn;
    if (s == "Sn") return PhaseType::Sn;
    if (s == "Pg") return PhaseType::Pg;
    if (s == "Sg") return PhaseType::Sg;
    if (s == "Lg") return PhaseType::Lg;
    return PhaseType::Unknown;
}

// Magnitude types
enum class MagnitudeType {
    ML,     // Local magnitude
    Mw,     // Moment magnitude
    Mb,     // Body wave magnitude
    Ms,     // Surface wave magnitude
    Md,     // Duration magnitude
    Unknown
};

inline std::string magnitudeTypeToString(MagnitudeType mt) {
    switch (mt) {
        case MagnitudeType::ML: return "ML";
        case MagnitudeType::Mw: return "Mw";
        case MagnitudeType::Mb: return "Mb";
        case MagnitudeType::Ms: return "Ms";
        case MagnitudeType::Md: return "Md";
        default: return "?";
    }
}

// Pick quality/weight
enum class PickQuality {
    Impulsive,  // 0 - best
    Emergent,   // 1
    Questionable, // 2
    Poor,       // 3
    Rejected    // 4 - worst
};

// Constants
namespace constants {
    constexpr double EARTH_RADIUS_KM = 6371.0;
    constexpr double DEG_TO_RAD = M_PI / 180.0;
    constexpr double RAD_TO_DEG = 180.0 / M_PI;
    constexpr double KM_PER_DEG = 111.195;  // Approximate
}

} // namespace realdetect
