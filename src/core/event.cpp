#include "realdetect/core/event.hpp"
#include <sstream>
#include <iomanip>
#include <random>
#include <chrono>

namespace realdetect {

// Static member definitions
Origin Event::dummy_origin_;
Magnitude Event::dummy_magnitude_;

std::string Event::generateId() {
    auto now = std::chrono::system_clock::now();
    auto time_t = std::chrono::system_clock::to_time_t(now);
    auto tm = *std::gmtime(&time_t);
    
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(1000, 9999);
    
    std::ostringstream oss;
    oss << "seisproc" 
        << std::put_time(&tm, "%Y%m%d%H%M%S")
        << "_" << dis(gen);
    return oss.str();
}

std::string Event::summary() const {
    std::ostringstream oss;
    
    const Origin& o = preferredOrigin();
    const Magnitude& m = preferredMagnitude();
    
    // Format time
    auto time_t = std::chrono::system_clock::to_time_t(o.time);
    auto tm = *std::gmtime(&time_t);
    
    oss << std::fixed << std::setprecision(3);
    oss << "Event: " << id_ << "\n";
    oss << "  Time: " << std::put_time(&tm, "%Y-%m-%d %H:%M:%S") << " UTC\n";
    oss << "  Location: " << o.location.latitude << "° N, " 
        << o.location.longitude << "° E\n";
    oss << "  Depth: " << o.location.depth << " km\n";
    
    if (m.type != MagnitudeType::Unknown) {
        oss << "  Magnitude: " << m.value << " " << magnitudeTypeToString(m.type) << "\n";
    }
    
    oss << "  Quality: " << o.qualityCode() << " (RMS=" << o.rms << "s, Gap=" 
        << o.gap << "°, Phases=" << o.phase_count << ")\n";
    
    return oss.str();
}

} // namespace realdetect
