#include "realdetect/core/velocity_model.hpp"
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>

namespace realdetect {

double VelocityModel1D::travelTime(double distance, double depth, PhaseType phase) const {
    if (layers_.empty()) return 0;
    
    bool is_p = (phase == PhaseType::P || phase == PhaseType::Pn || phase == PhaseType::Pg);
    
    // Build flat layer stack: each entry is (top_depth, bottom_depth, velocity)
    struct FlatLayer { double top, bot, vel; };
    std::vector<FlatLayer> stack;
    for (size_t i = 0; i < layers_.size(); i++) {
        double v = is_p ? layers_[i].vp : layers_[i].vs;
        if (v <= 0) continue;
        double bot = (layers_[i].thickness > 0) ?
            layers_[i].top_depth + layers_[i].thickness : 1e4;
        stack.push_back({layers_[i].top_depth, bot, v});
    }
    if (stack.empty()) return distance / 6.0;
    
    // --- 1) Direct ray: straight-line through average velocity ---
    //     (simple but used as upper-bound fallback)
    double avg_v = 0;
    double total_thick = 0;
    for (auto& L : stack) {
        double lo = L.top, hi = std::min(L.bot, std::max(depth, 1.0));
        if (lo >= hi) continue;
        avg_v += L.vel * (hi - lo);
        total_thick += (hi - lo);
    }
    if (total_thick > 0) avg_v /= total_thick;
    else avg_v = stack[0].vel;
    double hypo = std::sqrt(distance * distance + depth * depth);
    double t_direct = hypo / avg_v;
    
    // --- 2) Head-wave travel times at each velocity increase ---
    //     For each interface where velocity increases, compute the
    //     critically refracted (head wave) travel time.
    //     t_head = sum_i (2 * h_i / v_i * cos(ic_i)) + X_head / v_refractor
    //     where ic_i = arcsin(v_i / v_refractor)
    //     and X_head = distance - sum_i (2 * h_i * tan(ic_i))
    double t_best = t_direct;
    
    for (size_t ref = 1; ref < stack.size(); ref++) {
        double v_ref = stack[ref].vel;
        // Only consider interfaces where velocity increases
        if (v_ref <= stack[ref - 1].vel) continue;
        
        double interface_depth = stack[ref].top;
        // Source must be above the interface (or at it)
        if (depth > interface_depth + 1.0) continue;
        
        // Compute vertical travel time through each layer above interface
        double t_vertical = 0;
        double x_offset = 0;  // horizontal distance consumed by down+up legs
        bool valid = true;
        
        // Layers from source depth down to interface, then back up to surface
        // For simplicity, assume source at surface (depth~0) or correct for it
        for (size_t i = 0; i < ref; i++) {
            double v_i = stack[i].vel;
            if (v_i >= v_ref) { valid = false; break; } // Can't refract
            
            double ic = std::asin(v_i / v_ref);  // critical angle
            double h_i = stack[i].bot - stack[i].top;
            if (stack[i].bot > interface_depth)
                h_i = interface_depth - stack[i].top;
            if (h_i <= 0) continue;
            
            // Account for source depth being inside a layer
            double h_eff = h_i;
            if (stack[i].top < depth && depth < stack[i].bot) {
                h_eff = stack[i].bot - depth;
                if (stack[i].bot > interface_depth)
                    h_eff = interface_depth - depth;
            } else if (stack[i].bot <= depth) {
                // Layer is entirely above source — on the upgoing leg only
                // For head wave: down-leg uses source to interface,
                // up-leg uses interface to surface
                // Simplification: use full layer for the receiver-side leg
            }
            
            t_vertical += h_i / (v_i * std::cos(ic));  // receiver-side down+up leg
            x_offset += h_i * std::tan(ic);
        }
        
        if (!valid) continue;
        
        // Also add source-side vertical path (source to interface)
        double t_source = 0;
        double x_source = 0;
        for (size_t i = 0; i < ref; i++) {
            double v_i = stack[i].vel;
            if (v_i >= v_ref) { valid = false; break; }
            double ic = std::asin(v_i / v_ref);
            double layer_top = stack[i].top;
            double layer_bot = std::min(stack[i].bot, interface_depth);
            
            // Source-side: from source depth down to interface
            double h_src = 0;
            if (layer_bot > depth && layer_top < interface_depth) {
                h_src = layer_bot - std::max(layer_top, depth);
            }
            if (h_src > 0) {
                t_source += h_src / (v_i * std::cos(ic));
                x_source += h_src * std::tan(ic);
            }
        }
        if (!valid) continue;
        
        // Horizontal distance along the refractor
        double x_refractor = distance - x_offset - x_source;
        if (x_refractor < 0) continue;  // Too close for head wave
        
        double t_head = t_vertical + t_source + x_refractor / v_ref;
        if (t_head > 0 && t_head < t_best) {
            t_best = t_head;
        }
    }
    
    // --- 3) Layer-cake direct ray (better than straight-line average) ---
    //     Trace a ray through each layer using an average slowness approach
    if (distance > 0) {
        // Use the horizontal slowness p = sin(takeoff)/v_source
        // Try several takeoff angles and pick the one closest to target distance
        double best_tt_ray = t_direct;
        double src_v = is_p ? vpAt(depth) : vsAt(depth);
        
        // Binary search for the right takeoff angle
        double p_lo = 0.0001;
        double p_hi = 0.999 / src_v;  // max horizontal slowness
        
        for (int iter = 0; iter < 30; iter++) {
            double p_mid = (p_lo + p_hi) / 2.0;
            
            // Trace ray with this slowness: compute distance and time
            double X = 0, T = 0;
            bool ray_valid = true;
            
            // Downgoing from source to bottom of model (or turning point)
            for (size_t i = 0; i < stack.size(); i++) {
                double v = stack[i].vel;
                double layer_top = stack[i].top;
                double layer_bot = stack[i].bot;
                if (layer_bot > 200) layer_bot = 200; // limit depth
                
                if (layer_top >= 200) break;
                if (layer_bot <= depth) continue;
                
                double z_top = std::max(layer_top, depth);
                double z_bot = layer_bot;
                double dz = z_bot - z_top;
                if (dz <= 0) continue;
                
                double sin_val = p_mid * v;
                if (sin_val >= 1.0) { ray_valid = false; break; } // turning point
                
                double cos_val = std::sqrt(1.0 - sin_val * sin_val);
                X += dz * sin_val / cos_val;
                T += dz / (v * cos_val);
            }
            // Upgoing from depth back to surface
            for (size_t i = 0; i < stack.size(); i++) {
                double v = stack[i].vel;
                double layer_top = stack[i].top;
                double layer_bot = stack[i].bot;
                if (layer_bot > 200) layer_bot = 200;
                if (layer_top >= 200) break;
                
                double z_top = layer_top;
                double z_bot = std::min(layer_bot, depth);
                // For surface to source: the upgoing leg goes from interface to surface
                // Actually for source at depth, upgoing = surface layers
                z_bot = layer_bot;
                if (z_bot > depth) z_bot = depth;
                if (z_bot <= z_top) {
                    // Above source — this is the receiver-side upgoing leg
                    z_top = layer_top;
                    z_bot = std::min(layer_bot, depth);
                    if (z_bot <= z_top) continue;
                }
                double dz = z_bot - z_top;
                if (dz <= 0) continue;
                
                double sin_val = p_mid * v;
                if (sin_val >= 1.0) { ray_valid = false; break; }
                double cos_val = std::sqrt(1.0 - sin_val * sin_val);
                X += dz * sin_val / cos_val;
                T += dz / (v * cos_val);
            }
            
            if (!ray_valid) {
                p_hi = p_mid;
                continue;
            }
            
            if (X < distance) p_lo = p_mid;
            else p_hi = p_mid;
            
            if (std::abs(X - distance) < 0.5) { // within 0.5 km
                if (T > 0 && T < best_tt_ray) best_tt_ray = T;
                break;
            }
        }
        if (best_tt_ray > 0 && best_tt_ray < t_best) {
            t_best = best_tt_ray;
        }
    }
    
    return t_best;
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

} // namespace realdetect
