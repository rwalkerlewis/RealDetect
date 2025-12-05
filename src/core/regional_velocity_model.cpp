/**
 * Regional Velocity Model Manager Implementation
 */

#include "realdetect/core/regional_velocity_model.hpp"
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <filesystem>

namespace realdetect {

//=============================================================================
// GeographicBounds Implementation
//=============================================================================

GeographicBounds GeographicBounds::polygon(const std::vector<GeoPoint>& vertices) {
    GeographicBounds bounds;
    bounds.type_ = Type::Polygon;
    bounds.vertices_ = vertices;
    bounds.min_depth_ = 0;
    bounds.max_depth_ = 700;
    
    // Calculate bounding box
    bounds.min_lat_ = 90;
    bounds.max_lat_ = -90;
    bounds.min_lon_ = 180;
    bounds.max_lon_ = -180;
    
    for (const auto& v : vertices) {
        bounds.min_lat_ = std::min(bounds.min_lat_, v.latitude);
        bounds.max_lat_ = std::max(bounds.max_lat_, v.latitude);
        bounds.min_lon_ = std::min(bounds.min_lon_, v.longitude);
        bounds.max_lon_ = std::max(bounds.max_lon_, v.longitude);
    }
    
    return bounds;
}

GeographicBounds GeographicBounds::circle(const GeoPoint& center, double radius_km) {
    GeographicBounds bounds;
    bounds.type_ = Type::Circle;
    bounds.center_ = center;
    bounds.radius_ = radius_km;
    bounds.min_depth_ = 0;
    bounds.max_depth_ = 700;
    
    // Approximate bounding box (rough conversion)
    double deg_lat = radius_km / 111.0;
    double deg_lon = radius_km / (111.0 * std::cos(center.latitude * M_PI / 180.0));
    
    bounds.min_lat_ = center.latitude - deg_lat;
    bounds.max_lat_ = center.latitude + deg_lat;
    bounds.min_lon_ = center.longitude - deg_lon;
    bounds.max_lon_ = center.longitude + deg_lon;
    
    return bounds;
}

GeographicBounds GeographicBounds::boundingBox(double min_lat, double max_lat,
                                                 double min_lon, double max_lon,
                                                 double min_depth, double max_depth) {
    GeographicBounds bounds;
    bounds.type_ = Type::BoundingBox;
    bounds.min_lat_ = min_lat;
    bounds.max_lat_ = max_lat;
    bounds.min_lon_ = min_lon;
    bounds.max_lon_ = max_lon;
    bounds.min_depth_ = min_depth;
    bounds.max_depth_ = max_depth;
    
    return bounds;
}

bool GeographicBounds::contains(const GeoPoint& point) const {
    // Quick bounding box check
    if (point.latitude < min_lat_ || point.latitude > max_lat_ ||
        point.longitude < min_lon_ || point.longitude > max_lon_) {
        return false;
    }
    
    // Depth check
    if (point.depth < min_depth_ || point.depth > max_depth_) {
        return false;
    }
    
    switch (type_) {
        case Type::BoundingBox:
            return true;  // Already passed bounding box check
            
        case Type::Circle: {
            double dist = center_.distanceTo(point);
            return dist <= radius_;
        }
        
        case Type::Polygon:
            return pointInPolygon(point);
    }
    
    return false;
}

void GeographicBounds::getBoundingBox(double& min_lat, double& max_lat,
                                       double& min_lon, double& max_lon) const {
    min_lat = min_lat_;
    max_lat = max_lat_;
    min_lon = min_lon_;
    max_lon = max_lon_;
}

bool GeographicBounds::pointInPolygon(const GeoPoint& point) const {
    // Ray casting algorithm
    if (vertices_.size() < 3) return false;
    
    int n = vertices_.size();
    int count = 0;
    
    for (int i = 0, j = n - 1; i < n; j = i++) {
        double yi = vertices_[i].latitude;
        double yj = vertices_[j].latitude;
        double xi = vertices_[i].longitude;
        double xj = vertices_[j].longitude;
        
        if (((yi > point.latitude) != (yj > point.latitude)) &&
            (point.longitude < (xj - xi) * (point.latitude - yi) / (yj - yi) + xi)) {
            count++;
        }
    }
    
    return (count % 2) == 1;
}

std::string GeographicBounds::toString() const {
    std::ostringstream oss;
    
    switch (type_) {
        case Type::BoundingBox:
            oss << "box:" << min_lat_ << "," << max_lat_ << ","
                << min_lon_ << "," << max_lon_;
            break;
            
        case Type::Circle:
            oss << "circle:" << center_.latitude << "," << center_.longitude
                << "," << radius_;
            break;
            
        case Type::Polygon:
            oss << "polygon:";
            for (size_t i = 0; i < vertices_.size(); i++) {
                if (i > 0) oss << ";";
                oss << vertices_[i].latitude << "," << vertices_[i].longitude;
            }
            break;
    }
    
    return oss.str();
}

GeographicBounds GeographicBounds::fromString(const std::string& str) {
    size_t colon = str.find(':');
    if (colon == std::string::npos) {
        return boundingBox(-90, 90, -180, 180);  // Global default
    }
    
    std::string type = str.substr(0, colon);
    std::string spec = str.substr(colon + 1);
    
    if (type == "box") {
        double min_lat, max_lat, min_lon, max_lon;
        char c;
        std::istringstream iss(spec);
        iss >> min_lat >> c >> max_lat >> c >> min_lon >> c >> max_lon;
        return boundingBox(min_lat, max_lat, min_lon, max_lon);
    }
    else if (type == "circle") {
        double lat, lon, radius;
        char c;
        std::istringstream iss(spec);
        iss >> lat >> c >> lon >> c >> radius;
        return circle(GeoPoint(lat, lon), radius);
    }
    else if (type == "polygon") {
        std::vector<GeoPoint> vertices;
        std::istringstream iss(spec);
        std::string point_str;
        while (std::getline(iss, point_str, ';')) {
            double lat, lon;
            char c;
            std::istringstream pss(point_str);
            pss >> lat >> c >> lon;
            vertices.emplace_back(lat, lon);
        }
        if (!vertices.empty()) {
            return polygon(vertices);
        }
    }
    
    return boundingBox(-90, 90, -180, 180);
}

//=============================================================================
// VelocityModelManager Implementation
//=============================================================================

VelocityModelManager::VelocityModelManager() {
    // Initialize with IASP91 as default
    default_model_ = VelocityModel1D::iasp91();
}

void VelocityModelManager::addModel(const RegionalVelocityModel& model) {
    models_.push_back(model);
    sortByPriority();
}

void VelocityModelManager::addModel(const std::string& name, 
                                     const GeographicBounds& bounds,
                                     const VelocityModel1D& vmodel, 
                                     int priority) {
    RegionalVelocityModel rm;
    rm.name = name;
    rm.bounds = bounds;
    rm.model = vmodel;
    rm.priority = priority;
    addModel(rm);
}

bool VelocityModelManager::removeModel(const std::string& name) {
    auto it = std::remove_if(models_.begin(), models_.end(),
        [&name](const RegionalVelocityModel& m) { return m.name == name; });
    
    if (it != models_.end()) {
        models_.erase(it, models_.end());
        return true;
    }
    return false;
}

void VelocityModelManager::setDefaultModel(const VelocityModel1D& model) {
    default_model_ = model;
}

const VelocityModel1D& VelocityModelManager::getModelForLocation(const GeoPoint& location) const {
    // Models are sorted by priority, so first match wins
    for (const auto& rm : models_) {
        if (rm.contains(location)) {
            return rm.model;
        }
    }
    return default_model_;
}

const VelocityModel1D& VelocityModelManager::getModelForLocation(double lat, double lon, double depth) const {
    return getModelForLocation(GeoPoint(lat, lon, depth));
}

const VelocityModel1D* VelocityModelManager::getModelByName(const std::string& name) const {
    for (const auto& rm : models_) {
        if (rm.name == name) {
            return &rm.model;
        }
    }
    
    if (name == "default" || name == default_model_.name()) {
        return &default_model_;
    }
    
    return nullptr;
}

std::string VelocityModelManager::getModelNameForLocation(const GeoPoint& location) const {
    for (const auto& rm : models_) {
        if (rm.contains(location)) {
            return rm.name;
        }
    }
    return default_model_.name().empty() ? "default" : default_model_.name();
}

std::vector<std::string> VelocityModelManager::modelNames() const {
    std::vector<std::string> names;
    names.reserve(models_.size() + 1);
    
    for (const auto& rm : models_) {
        names.push_back(rm.name);
    }
    names.push_back(default_model_.name().empty() ? "default" : default_model_.name());
    
    return names;
}

bool VelocityModelManager::loadFromDirectory(const std::string& directory) {
    namespace fs = std::filesystem;
    
    if (!fs::exists(directory) || !fs::is_directory(directory)) {
        std::cerr << "Velocity model directory not found: " << directory << std::endl;
        return false;
    }
    
    // Look for config file first
    std::string config_file = directory + "/velocity_models.conf";
    if (fs::exists(config_file)) {
        return loadFromConfig(config_file);
    }
    
    // Otherwise, load all .vel, .mod, .txt files as velocity models
    int loaded = 0;
    for (const auto& entry : fs::directory_iterator(directory)) {
        if (!entry.is_regular_file()) continue;
        
        std::string ext = entry.path().extension().string();
        if (ext != ".vel" && ext != ".mod" && ext != ".txt") continue;
        
        VelocityModel1D model;
        if (VelocityModelLoader::load(entry.path().string(), model)) {
            // Use filename (without extension) as model name
            std::string name = entry.path().stem().string();
            model = VelocityModel1D(name);
            
            // Check for accompanying .bounds file
            std::string bounds_file = entry.path().string();
            bounds_file = bounds_file.substr(0, bounds_file.rfind('.')) + ".bounds";
            
            if (fs::exists(bounds_file)) {
                std::ifstream bf(bounds_file);
                std::string bounds_str;
                std::getline(bf, bounds_str);
                
                RegionalVelocityModel rm;
                rm.name = name;
                rm.model = model;
                rm.bounds = GeographicBounds::fromString(bounds_str);
                addModel(rm);
            } else {
                // Add as global model candidate
                if (loaded == 0) {
                    setDefaultModel(model);
                }
            }
            loaded++;
        }
    }
    
    std::cout << "Loaded " << loaded << " velocity models from " << directory << std::endl;
    return loaded > 0;
}

bool VelocityModelManager::loadFromConfig(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open velocity model config: " << filename << std::endl;
        return false;
    }
    
    std::string dir = filename.substr(0, filename.rfind('/'));
    if (dir.empty()) dir = ".";
    
    std::string line;
    std::string current_section;
    RegionalModelConfig config;
    
    auto processConfig = [this, &dir](const RegionalModelConfig& cfg) {
        if (cfg.name.empty() || cfg.model_file.empty()) return;
        
        std::string model_path = cfg.model_file;
        if (model_path[0] != '/') {
            model_path = dir + "/" + model_path;
        }
        
        VelocityModel1D model(cfg.name);
        if (!VelocityModelLoader::load(model_path, model)) {
            std::cerr << "Failed to load velocity model: " << model_path << std::endl;
            return;
        }
        
        if (cfg.bounds_type == "default" || cfg.bounds_spec.empty()) {
            // This is the default model
            setDefaultModel(model);
        } else {
            RegionalVelocityModel rm;
            rm.name = cfg.name;
            rm.description = cfg.description;
            rm.model = model;
            rm.priority = cfg.priority;
            rm.bounds = cfg.parseBounds();
            addModel(rm);
        }
    };
    
    while (std::getline(file, line)) {
        // Trim whitespace
        line.erase(0, line.find_first_not_of(" \t\r\n"));
        if (line.empty() || line[0] == '#' || line[0] == ';') continue;
        line.erase(line.find_last_not_of(" \t\r\n") + 1);
        
        // Section header
        if (line[0] == '[' && line.back() == ']') {
            if (!config.name.empty()) {
                processConfig(config);
            }
            current_section = line.substr(1, line.size() - 2);
            config = RegionalModelConfig();
            config.name = current_section;
            continue;
        }
        
        // Key=value pair
        auto pos = line.find('=');
        if (pos == std::string::npos) continue;
        
        std::string key = line.substr(0, pos);
        std::string value = line.substr(pos + 1);
        key.erase(key.find_last_not_of(" \t") + 1);
        value.erase(0, value.find_first_not_of(" \t"));
        
        if (key == "model_file") config.model_file = value;
        else if (key == "bounds_type") config.bounds_type = value;
        else if (key == "bounds") config.bounds_spec = value;
        else if (key == "priority") config.priority = std::stoi(value);
        else if (key == "description") config.description = value;
    }
    
    // Process last config
    if (!config.name.empty()) {
        processConfig(config);
    }
    
    std::cout << "Loaded " << models_.size() << " regional velocity models" << std::endl;
    return true;
}

bool VelocityModelManager::saveConfig(const std::string& filename) const {
    std::ofstream file(filename);
    if (!file.is_open()) return false;
    
    file << "# Regional Velocity Model Configuration\n";
    file << "# Generated by SeisProc\n\n";
    
    // Write default model
    file << "[default]\n";
    file << "model_file = " << default_model_.name() << ".vel\n";
    file << "bounds_type = default\n";
    file << "description = Default global velocity model\n\n";
    
    // Write regional models
    for (const auto& rm : models_) {
        file << "[" << rm.name << "]\n";
        file << "model_file = " << rm.name << ".vel\n";
        file << "bounds_type = " << (rm.bounds.type() == GeographicBounds::Type::BoundingBox ? "box" :
                                     rm.bounds.type() == GeographicBounds::Type::Circle ? "circle" : "polygon") << "\n";
        file << "bounds = " << rm.bounds.toString() << "\n";
        file << "priority = " << rm.priority << "\n";
        if (!rm.description.empty()) {
            file << "description = " << rm.description << "\n";
        }
        file << "\n";
    }
    
    return true;
}

std::vector<std::string> VelocityModelManager::findMatchingModels(const GeoPoint& location) const {
    std::vector<std::string> matches;
    
    for (const auto& rm : models_) {
        if (rm.contains(location)) {
            matches.push_back(rm.name);
        }
    }
    
    return matches;
}

void VelocityModelManager::sortByPriority() {
    std::sort(models_.begin(), models_.end(),
        [](const RegionalVelocityModel& a, const RegionalVelocityModel& b) {
            return a.priority > b.priority;  // Higher priority first
        });
}

//=============================================================================
// VelocityModelLoader Implementation
//=============================================================================

bool VelocityModelLoader::load(const std::string& filename, VelocityModel1D& model) {
    // Try to detect format from file extension or content
    std::string ext = filename.substr(filename.rfind('.') + 1);
    
    if (ext == "hyp" || ext == "hypo71") {
        return loadHypo71(filename, model);
    }
    else if (ext == "mod" || ext == "velest") {
        return loadVelest(filename, model);
    }
    else if (ext == "nll") {
        return loadNonLinLoc(filename, model);
    }
    
    // Default: SeisProc format
    return loadSeisProc(filename, model);
}

bool VelocityModelLoader::loadSeisProc(const std::string& filename, VelocityModel1D& model) {
    return model.loadFromFile(filename);
}

bool VelocityModelLoader::loadHypo71(const std::string& filename, VelocityModel1D& model) {
    std::ifstream file(filename);
    if (!file.is_open()) return false;
    
    std::string line;
    std::vector<std::pair<double, double>> layers;  // depth, velocity pairs
    double vpvs = 1.73;  // Default Vp/Vs ratio
    
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;
        
        std::istringstream iss(line);
        double val1, val2;
        
        // HYPO71 format: velocity at top of layer, then depth to top
        if (iss >> val1 >> val2) {
            // Check if this is VPVS ratio line (usually first line if single value)
            if (layers.empty() && !(iss >> val2)) {
                vpvs = val1;
                continue;
            }
            layers.emplace_back(val2, val1);  // depth, velocity
        }
    }
    
    // Convert to SeisProc format
    for (size_t i = 0; i < layers.size(); i++) {
        double depth = layers[i].first;
        double vp = layers[i].second;
        double vs = vp / vpvs;
        double thickness = (i + 1 < layers.size()) ? 
            layers[i + 1].first - depth : 0;  // 0 = halfspace
        
        model.addLayer(depth, thickness, vp, vs);
    }
    
    return !layers.empty();
}

bool VelocityModelLoader::loadVelest(const std::string& filename, VelocityModel1D& model) {
    std::ifstream file(filename);
    if (!file.is_open()) return false;
    
    std::string line;
    bool reading_p = false;
    bool reading_s = false;
    std::vector<std::pair<double, double>> p_layers;
    std::vector<double> s_velocities;
    
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        
        // VELEST format markers
        if (line.find("P-VELOCITY") != std::string::npos) {
            reading_p = true;
            reading_s = false;
            continue;
        }
        if (line.find("S-VELOCITY") != std::string::npos) {
            reading_p = false;
            reading_s = true;
            continue;
        }
        
        std::istringstream iss(line);
        double vel, depth;
        
        if (reading_p && iss >> vel >> depth) {
            p_layers.emplace_back(depth, vel);
        }
        else if (reading_s && iss >> vel) {
            s_velocities.push_back(vel);
        }
    }
    
    // Convert to SeisProc format
    for (size_t i = 0; i < p_layers.size(); i++) {
        double depth = p_layers[i].first;
        double vp = p_layers[i].second;
        double vs = (i < s_velocities.size()) ? s_velocities[i] : vp / 1.73;
        double thickness = (i + 1 < p_layers.size()) ?
            p_layers[i + 1].first - depth : 0;
        
        model.addLayer(depth, thickness, vp, vs);
    }
    
    return !p_layers.empty();
}

bool VelocityModelLoader::loadNonLinLoc(const std::string& filename, VelocityModel1D& model) {
    std::ifstream file(filename);
    if (!file.is_open()) return false;
    
    std::string line;
    
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;
        
        // NonLinLoc format: LAYER depth vp_top vp_grad vs_top vs_grad rho_top rho_grad
        if (line.substr(0, 5) != "LAYER") continue;
        
        std::istringstream iss(line.substr(5));
        double depth, vp_top, vp_grad, vs_top, vs_grad, rho_top, rho_grad;
        
        if (iss >> depth >> vp_top >> vp_grad >> vs_top >> vs_grad >> rho_top >> rho_grad) {
            // For simplicity, use top velocities (ignore gradients)
            // Thickness will be calculated when next layer is added
            model.addLayer(depth, 0, vp_top, vs_top, rho_top);
        }
    }
    
    return model.layerCount() > 0;
}

bool VelocityModelLoader::saveSeisProc(const std::string& filename, const VelocityModel1D& model) {
    std::ofstream file(filename);
    if (!file.is_open()) return false;
    
    file << "# SeisProc 1D Velocity Model\n";
    file << "# Model: " << model.name() << "\n";
    file << "# Format: Depth(km) Thickness(km) Vp(km/s) Vs(km/s) Density(g/cm³)\n\n";
    
    for (const auto& layer : model.layers()) {
        file << std::fixed << std::setprecision(2)
             << layer.top_depth << "\t"
             << layer.thickness << "\t"
             << layer.vp << "\t"
             << layer.vs << "\t"
             << layer.density << "\n";
    }
    
    return true;
}

//=============================================================================
// RegionalModelConfig Implementation
//=============================================================================

GeographicBounds RegionalModelConfig::parseBounds() const {
    if (bounds_type == "box") {
        // Format: min_lat,max_lat,min_lon,max_lon
        double min_lat, max_lat, min_lon, max_lon;
        char c;
        std::istringstream iss(bounds_spec);
        iss >> min_lat >> c >> max_lat >> c >> min_lon >> c >> max_lon;
        return GeographicBounds::boundingBox(min_lat, max_lat, min_lon, max_lon);
    }
    else if (bounds_type == "circle") {
        // Format: center_lat,center_lon,radius_km
        double lat, lon, radius;
        char c;
        std::istringstream iss(bounds_spec);
        iss >> lat >> c >> lon >> c >> radius;
        return GeographicBounds::circle(GeoPoint(lat, lon), radius);
    }
    else if (bounds_type == "polygon") {
        // Format: lat1,lon1;lat2,lon2;...
        std::vector<GeoPoint> vertices;
        std::istringstream iss(bounds_spec);
        std::string point_str;
        while (std::getline(iss, point_str, ';')) {
            double lat, lon;
            char c;
            std::istringstream pss(point_str);
            pss >> lat >> c >> lon;
            vertices.emplace_back(lat, lon);
        }
        return GeographicBounds::polygon(vertices);
    }
    
    // Default: global bounds
    return GeographicBounds::boundingBox(-90, 90, -180, 180);
}

} // namespace realdetect
