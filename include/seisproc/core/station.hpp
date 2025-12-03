#pragma once

#include "types.hpp"
#include <map>
#include <memory>

namespace seisproc {

/**
 * Channel - Individual seismic channel information
 */
struct Channel {
    std::string code;           // e.g., BHZ, HHN
    double sample_rate;         // Hz
    double azimuth;             // degrees from north
    double dip;                 // degrees from horizontal (-90 to 90)
    double gain;                // counts per m/s
    std::string sensor_type;
    
    Channel() : sample_rate(0), azimuth(0), dip(0), gain(1.0) {}
    
    bool isVertical() const { return dip < -60 || dip > 60; }
    bool isHorizontal() const { return std::abs(dip) < 30; }
};

/**
 * Station - Seismic station with metadata
 */
class Station {
public:
    Station() : elevation_(0) {}
    
    Station(const std::string& network, const std::string& code,
            double lat, double lon, double elev = 0)
        : network_(network), code_(code),
          location_(lat, lon, -elev/1000.0), elevation_(elev) {}
    
    // Accessors
    const std::string& network() const { return network_; }
    const std::string& code() const { return code_; }
    const GeoPoint& location() const { return location_; }
    double latitude() const { return location_.latitude; }
    double longitude() const { return location_.longitude; }
    double elevation() const { return elevation_; }
    
    // Channel management
    void addChannel(const Channel& chan) { channels_[chan.code] = chan; }
    const Channel* getChannel(const std::string& code) const {
        auto it = channels_.find(code);
        return it != channels_.end() ? &it->second : nullptr;
    }
    const std::map<std::string, Channel>& channels() const { return channels_; }
    
    // Get stream ID for channel
    StreamID streamId(const std::string& location, const std::string& channel) const {
        return StreamID(network_, code_, location, channel);
    }
    
    // Distance to a point
    double distanceTo(const GeoPoint& point) const {
        return location_.distanceTo(point);
    }
    
    // Azimuth to a point
    double azimuthTo(const GeoPoint& point) const {
        return location_.azimuthTo(point);
    }

private:
    std::string network_;
    std::string code_;
    GeoPoint location_;
    double elevation_;  // meters above sea level
    std::map<std::string, Channel> channels_;
};

using StationPtr = std::shared_ptr<Station>;

/**
 * StationInventory - Collection of stations
 */
class StationInventory {
public:
    void addStation(StationPtr station) {
        std::string key = station->network() + "." + station->code();
        stations_[key] = station;
    }
    
    StationPtr getStation(const std::string& network, const std::string& code) const {
        std::string key = network + "." + code;
        auto it = stations_.find(key);
        return it != stations_.end() ? it->second : nullptr;
    }
    
    StationPtr getStation(const StreamID& id) const {
        return getStation(id.network, id.station);
    }
    
    const std::map<std::string, StationPtr>& stations() const { return stations_; }
    size_t size() const { return stations_.size(); }
    
    // Find stations within distance
    std::vector<StationPtr> stationsWithin(const GeoPoint& center, double max_dist_km) const {
        std::vector<StationPtr> result;
        for (const auto& [key, sta] : stations_) {
            if (sta->distanceTo(center) <= max_dist_km) {
                result.push_back(sta);
            }
        }
        return result;
    }
    
    // Load from file
    bool loadFromFile(const std::string& filename);
    bool saveToFile(const std::string& filename) const;

private:
    std::map<std::string, StationPtr> stations_;
};

} // namespace seisproc
