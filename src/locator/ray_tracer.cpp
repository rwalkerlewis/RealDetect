#include "realdetect/core/velocity_model.hpp"
#include <cmath>
#include <vector>

namespace realdetect {

/**
 * Simple 1D ray tracer for travel time computation
 */
class RayTracer1D {
public:
    struct RayPath {
        std::vector<double> depths;
        std::vector<double> distances;
        double travel_time;
        double takeoff_angle;
        double incidence_angle;
    };
    
    RayTracer1D(const VelocityModel1D& model) : model_(model) {}
    
    // Trace ray from source to receiver
    RayPath trace(double source_depth, double receiver_depth, 
                  double distance, PhaseType phase) {
        RayPath path;
        
        // Get velocity at source
        double v_source = (phase == PhaseType::P) ? 
            model_.vpAt(source_depth) : model_.vsAt(source_depth);
        
        // Simple straight-line approximation for now
        // Full implementation would use Snell's law and layer integration
        
        double hypo_dist = std::sqrt(distance * distance + 
            (source_depth - receiver_depth) * (source_depth - receiver_depth));
        
        path.travel_time = hypo_dist / v_source;
        path.takeoff_angle = std::atan2(distance, source_depth - receiver_depth) 
                             * 180.0 / M_PI;
        path.incidence_angle = 90.0 - path.takeoff_angle;
        
        // Simple path
        path.depths = {source_depth, receiver_depth};
        path.distances = {0, distance};
        
        return path;
    }
    
    // Compute travel time using Tau-P method (simplified)
    double travelTime(double distance, double source_depth, PhaseType phase) {
        // For small distances, use direct ray
        if (distance < 100) {
            double v = (phase == PhaseType::P) ? 
                model_.vpAt(source_depth / 2) : model_.vsAt(source_depth / 2);
            double hypo = std::sqrt(distance * distance + source_depth * source_depth);
            return hypo / v;
        }
        
        // For larger distances, include refraction at Moho
        const auto& layers = model_.layers();
        if (layers.empty()) {
            return distance / 6.0;  // Default
        }
        
        double tt = 0;
        double remaining_dist = distance;
        double current_depth = source_depth;
        
        // Trace down to each layer and back up
        for (const auto& layer : layers) {
            if (current_depth >= layer.top_depth + layer.thickness) continue;
            
            double v = (phase == PhaseType::P) ? layer.vp : layer.vs;
            double layer_dist = std::min(remaining_dist, layer.thickness * 2);
            
            tt += layer_dist / v;
            remaining_dist -= layer_dist;
            
            if (remaining_dist <= 0) break;
        }
        
        return tt;
    }

private:
    const VelocityModel1D& model_;
};

} // namespace realdetect
