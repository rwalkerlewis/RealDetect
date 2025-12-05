#pragma once

#include "types.hpp"
#include <vector>
#include <memory>
#include <algorithm>

namespace realdetect {

/**
 * Layer in a 1D velocity model
 */
struct VelocityLayer {
    double top_depth;     // km
    double thickness;     // km (0 for halfspace)
    double vp;            // P velocity (km/s)
    double vs;            // S velocity (km/s)
    double density;       // g/cm³
    double qp;            // P-wave quality factor
    double qs;            // S-wave quality factor
    
    VelocityLayer() : top_depth(0), thickness(0), vp(6.0), vs(3.5),
                      density(2.7), qp(500), qs(250) {}
    
    VelocityLayer(double depth, double thick, double vp_, double vs_, 
                  double rho = 2.7, double qp_ = 500, double qs_ = 250)
        : top_depth(depth), thickness(thick), vp(vp_), vs(vs_),
          density(rho), qp(qp_), qs(qs_) {}
    
    double vpvs() const { return vs > 0 ? vp / vs : 1.73; }
    double poisson() const { 
        double r = vpvs();
        return (r*r - 2) / (2*(r*r - 1));
    }
};

/**
 * VelocityModel1D - 1D layered velocity model
 */
class VelocityModel1D {
public:
    VelocityModel1D() = default;
    VelocityModel1D(const std::string& name) : name_(name) {}
    
    // Add layers
    void addLayer(const VelocityLayer& layer) {
        layers_.push_back(layer);
        sortLayers();
    }
    
    void addLayer(double depth, double thickness, double vp, double vs,
                  double density = 2.7) {
        layers_.emplace_back(depth, thickness, vp, vs, density);
        sortLayers();
    }
    
    // Access
    const std::string& name() const { return name_; }
    const std::vector<VelocityLayer>& layers() const { return layers_; }
    size_t layerCount() const { return layers_.size(); }
    
    // Get velocity at depth
    double vpAt(double depth) const {
        for (size_t i = 0; i < layers_.size(); i++) {
            double bottom = (layers_[i].thickness > 0) ? 
                layers_[i].top_depth + layers_[i].thickness : 1e9;
            if (depth >= layers_[i].top_depth && depth < bottom) {
                return layers_[i].vp;
            }
        }
        return layers_.empty() ? 6.0 : layers_.back().vp;
    }
    
    double vsAt(double depth) const {
        for (size_t i = 0; i < layers_.size(); i++) {
            double bottom = (layers_[i].thickness > 0) ? 
                layers_[i].top_depth + layers_[i].thickness : 1e9;
            if (depth >= layers_[i].top_depth && depth < bottom) {
                return layers_[i].vs;
            }
        }
        return layers_.empty() ? 3.5 : layers_.back().vs;
    }
    
    // Get layer index at depth
    int layerIndexAt(double depth) const {
        for (size_t i = 0; i < layers_.size(); i++) {
            double bottom = (layers_[i].thickness > 0) ? 
                layers_[i].top_depth + layers_[i].thickness : 1e9;
            if (depth >= layers_[i].top_depth && depth < bottom) {
                return static_cast<int>(i);
            }
        }
        return layers_.size() - 1;
    }
    
    // Calculate travel time (simple linear interpolation)
    double travelTime(double distance, double depth, PhaseType phase) const;
    
    // Load from file
    bool loadFromFile(const std::string& filename);
    
    // Standard models
    static VelocityModel1D iasp91();
    static VelocityModel1D ak135();
    static VelocityModel1D simpleThreeLayer();

private:
    std::string name_;
    std::vector<VelocityLayer> layers_;
    
    void sortLayers() {
        std::sort(layers_.begin(), layers_.end(),
            [](const VelocityLayer& a, const VelocityLayer& b) {
                return a.top_depth < b.top_depth;
            });
    }
};

} // namespace realdetect
