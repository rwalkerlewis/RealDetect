#include "seisproc/core/station.hpp"
#include <fstream>
#include <sstream>
#include <iostream>

namespace seisproc {

bool StationInventory::loadFromFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open station file: " << filename << std::endl;
        return false;
    }
    
    std::string line;
    int line_num = 0;
    
    while (std::getline(file, line)) {
        line_num++;
        
        // Skip comments and empty lines
        if (line.empty() || line[0] == '#') continue;
        
        std::istringstream iss(line);
        std::string network, code;
        double lat, lon, elev = 0;
        
        // Format: network station latitude longitude [elevation]
        if (!(iss >> network >> code >> lat >> lon)) {
            std::cerr << "Parse error at line " << line_num << std::endl;
            continue;
        }
        iss >> elev;  // Optional elevation
        
        auto station = std::make_shared<Station>(network, code, lat, lon, elev);
        
        // Try to read channel info (optional)
        std::string chan_code;
        double sample_rate, azimuth, dip;
        while (iss >> chan_code >> sample_rate >> azimuth >> dip) {
            Channel chan;
            chan.code = chan_code;
            chan.sample_rate = sample_rate;
            chan.azimuth = azimuth;
            chan.dip = dip;
            station->addChannel(chan);
        }
        
        addStation(station);
    }
    
    std::cout << "Loaded " << stations_.size() << " stations from " 
              << filename << std::endl;
    return true;
}

bool StationInventory::saveToFile(const std::string& filename) const {
    std::ofstream file(filename);
    if (!file.is_open()) {
        return false;
    }
    
    file << "# SeisProc Station Inventory\n";
    file << "# Network Station Latitude Longitude Elevation\n";
    
    for (const auto& [key, sta] : stations_) {
        file << sta->network() << " " 
             << sta->code() << " "
             << sta->latitude() << " "
             << sta->longitude() << " "
             << sta->elevation() << "\n";
    }
    
    return true;
}

} // namespace seisproc
