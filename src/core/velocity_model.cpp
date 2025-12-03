#include "seisproc/core/velocity_model.hpp"
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>

namespace seisproc {

double VelocityModel1D::travelTime(double distance, double depth, PhaseType phase) const {
    if (layers_.empty()) return 0;
    
    // Simple approximation: average velocity ray path
    // For more accurate times, use ray tracing
    double velocity = (phase == PhaseType::P || phase == PhaseType::Pn || 
                       phase == PhaseType::Pg) ? vpAt(depth/2) : vsAt(depth/2);
    
    // Hypocentral distance
    double hypo_dist = std::sqrt(distance * distance + depth * depth);
    
    return hypo_dist / velocity;
}

bool VelocityModel1D::loadFromFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open velocity model: " << filename << std::endl;
        return false;
    }
    
    layers_.clear();
    std::string line;
    
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;
        
        std::istringstream iss(line);
        double depth, thickness, vp, vs;
        double density = 2.7;
        
        if (!(iss >> depth >> thickness >> vp >> vs)) continue;
        iss >> density;
        
        addLayer(depth, thickness, vp, vs, density);
    }
    
    std::cout << "Loaded velocity model with " << layers_.size() 
              << " layers from " << filename << std::endl;
    return !layers_.empty();
}

VelocityModel1D VelocityModel1D::iasp91() {
    VelocityModel1D model("IASP91");
    
    // Simplified IASP91 crustal model
    model.addLayer(0.0, 20.0, 5.80, 3.36, 2.60);   // Upper crust
    model.addLayer(20.0, 15.0, 6.50, 3.75, 2.90);  // Lower crust
    model.addLayer(35.0, 77.5, 8.04, 4.47, 3.38);  // Uppermost mantle
    model.addLayer(112.5, 52.5, 8.05, 4.50, 3.38); // Upper mantle
    model.addLayer(165.0, 45.0, 8.17, 4.51, 3.39); // Transition
    model.addLayer(210.0, 0.0, 8.30, 4.52, 3.40);  // Upper mantle (halfspace)
    
    return model;
}

VelocityModel1D VelocityModel1D::ak135() {
    VelocityModel1D model("AK135");
    
    // Simplified AK135 model
    model.addLayer(0.0, 20.0, 5.80, 3.46, 2.60);
    model.addLayer(20.0, 15.0, 6.50, 3.85, 2.92);
    model.addLayer(35.0, 75.0, 8.04, 4.48, 3.32);
    model.addLayer(110.0, 100.0, 8.05, 4.50, 3.45);
    model.addLayer(210.0, 0.0, 8.30, 4.52, 3.50);
    
    return model;
}

VelocityModel1D VelocityModel1D::simpleThreeLayer() {
    VelocityModel1D model("Simple3Layer");
    
    // Simple 3-layer crustal model
    model.addLayer(0.0, 15.0, 5.50, 3.18, 2.60);   // Sediments + upper crust
    model.addLayer(15.0, 20.0, 6.30, 3.64, 2.85);  // Lower crust
    model.addLayer(35.0, 0.0, 8.00, 4.62, 3.30);   // Mantle halfspace
    
    return model;
}

} // namespace seisproc
