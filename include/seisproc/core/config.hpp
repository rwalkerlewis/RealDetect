#pragma once

#include <string>
#include <map>
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>

namespace seisproc {

/**
 * Config - Simple configuration parser
 */
class Config {
public:
    Config() = default;
    
    bool loadFromFile(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) return false;
        
        std::string line;
        std::string section;
        while (std::getline(file, line)) {
            // Trim whitespace
            line.erase(0, line.find_first_not_of(" \t\r\n"));
            line.erase(line.find_last_not_of(" \t\r\n") + 1);
            
            // Skip empty lines and comments
            if (line.empty() || line[0] == '#' || line[0] == ';') continue;
            
            // Section header
            if (line[0] == '[' && line.back() == ']') {
                section = line.substr(1, line.size() - 2);
                continue;
            }
            
            // Key-value pair
            auto pos = line.find('=');
            if (pos != std::string::npos) {
                std::string key = line.substr(0, pos);
                std::string value = line.substr(pos + 1);
                key.erase(key.find_last_not_of(" \t") + 1);
                value.erase(0, value.find_first_not_of(" \t"));
                
                std::string full_key = section.empty() ? key : section + "." + key;
                values_[full_key] = value;
            }
        }
        return true;
    }
    
    bool saveToFile(const std::string& filename) const {
        std::ofstream file(filename);
        if (!file.is_open()) return false;
        
        std::string current_section;
        for (const auto& [key, value] : values_) {
            auto pos = key.find('.');
            std::string section = pos != std::string::npos ? key.substr(0, pos) : "";
            std::string name = pos != std::string::npos ? key.substr(pos + 1) : key;
            
            if (section != current_section) {
                if (!current_section.empty()) file << "\n";
                if (!section.empty()) file << "[" << section << "]\n";
                current_section = section;
            }
            file << name << " = " << value << "\n";
        }
        return true;
    }
    
    // Getters
    std::string getString(const std::string& key, const std::string& default_val = "") const {
        auto it = values_.find(key);
        return it != values_.end() ? it->second : default_val;
    }
    
    int getInt(const std::string& key, int default_val = 0) const {
        auto it = values_.find(key);
        if (it == values_.end()) return default_val;
        try { return std::stoi(it->second); }
        catch (...) { return default_val; }
    }
    
    double getDouble(const std::string& key, double default_val = 0.0) const {
        auto it = values_.find(key);
        if (it == values_.end()) return default_val;
        try { return std::stod(it->second); }
        catch (...) { return default_val; }
    }
    
    bool getBool(const std::string& key, bool default_val = false) const {
        auto it = values_.find(key);
        if (it == values_.end()) return default_val;
        std::string v = it->second;
        std::transform(v.begin(), v.end(), v.begin(), ::tolower);
        return v == "true" || v == "yes" || v == "1" || v == "on";
    }
    
    std::vector<std::string> getStringList(const std::string& key) const {
        std::vector<std::string> result;
        auto it = values_.find(key);
        if (it != values_.end()) {
            std::stringstream ss(it->second);
            std::string item;
            while (std::getline(ss, item, ',')) {
                item.erase(0, item.find_first_not_of(" \t"));
                item.erase(item.find_last_not_of(" \t") + 1);
                if (!item.empty()) result.push_back(item);
            }
        }
        return result;
    }
    
    // Setters
    void set(const std::string& key, const std::string& value) {
        values_[key] = value;
    }
    
    void set(const std::string& key, int value) {
        values_[key] = std::to_string(value);
    }
    
    void set(const std::string& key, double value) {
        values_[key] = std::to_string(value);
    }
    
    void set(const std::string& key, bool value) {
        values_[key] = value ? "true" : "false";
    }
    
    bool has(const std::string& key) const {
        return values_.find(key) != values_.end();
    }
    
    const std::map<std::string, std::string>& all() const { return values_; }

private:
    std::map<std::string, std::string> values_;
};

} // namespace seisproc
