#pragma once

#include "locator.hpp"

namespace seisproc {

/**
 * GridSearchLocator - Brute-force grid search for hypocenter
 * 
 * Searches over a 4D grid (lat, lon, depth, origin time) to find
 * the location that minimizes the travel time residuals.
 * 
 * Advantages:
 * - Always finds global minimum within search space
 * - No initial guess required
 * - Robust for sparse networks
 * 
 * Disadvantages:
 * - Computationally expensive for fine grids
 * - Resolution limited by grid spacing
 */
class GridSearchLocator : public BaseLocator {
public:
    GridSearchLocator();
    
    LocationResult locate(const std::vector<PickPtr>& picks,
                           const StationInventory& stations) override;
    
    LocationResult relocate(const EventPtr& event,
                             const StationInventory& stations) override;
    
    std::string name() const override { return "GridSearch"; }
    
    void setParameter(const std::string& name, double value) override;
    void setVelocityModel(const VelocityModel1D& model) override;
    
    // Grid configuration
    void setHorizontalStep(double km) { h_step_ = km; }
    void setDepthStep(double km) { z_step_ = km; }
    void setSearchRadius(double km) { search_radius_ = km; }
    void setDepthRange(double min_km, double max_km) { 
        depth_min_ = min_km; 
        depth_max_ = max_km; 
    }
    
    // Fixed depth option
    void setFixedDepth(bool fixed, double depth = 10.0) {
        fixed_depth_ = fixed;
        default_depth_ = depth;
    }
    
    // Multi-level refinement
    void setRefinementLevels(int levels) { refinement_levels_ = levels; }

private:
    VelocityModel1D velocity_model_;
    TravelTimeTable travel_times_;
    
    double h_step_;           // Horizontal grid step (km)
    double z_step_;           // Depth grid step (km)
    double search_radius_;    // Initial search radius (km)
    double depth_min_;        // Minimum depth (km)
    double depth_max_;        // Maximum depth (km)
    bool fixed_depth_;        // Fix depth to value
    double default_depth_;    // Default/fixed depth
    int refinement_levels_;   // Number of refinement iterations
    
    // Objective function (RMS residual)
    double computeRMS(double lat, double lon, double depth, TimePoint origin_time,
                      const std::vector<PickPtr>& picks,
                      const StationInventory& stations) const;
    
    // Get initial search center from picks
    GeoPoint estimateInitialLocation(const std::vector<PickPtr>& picks,
                                      const StationInventory& stations) const;
    
    // Estimate origin time for given location
    TimePoint estimateOriginTime(const GeoPoint& location,
                                  const std::vector<PickPtr>& picks,
                                  const StationInventory& stations) const;
};

/**
 * OctTreeLocator - Adaptive octree search
 * 
 * More efficient than uniform grid search by subdividing only
 * promising regions of the search space.
 */
class OctTreeLocator : public BaseLocator {
public:
    OctTreeLocator();
    
    LocationResult locate(const std::vector<PickPtr>& picks,
                           const StationInventory& stations) override;
    
    LocationResult relocate(const EventPtr& event,
                             const StationInventory& stations) override;
    
    std::string name() const override { return "OctTree"; }
    
    void setParameter(const std::string& name, double value) override;
    void setVelocityModel(const VelocityModel1D& model) override;
    
    // Configuration
    void setMinCellSize(double km) { min_cell_size_ = km; }
    void setMaxIterations(int iter) { max_iterations_ = iter; }
    void setSearchBounds(double lat_min, double lat_max,
                         double lon_min, double lon_max,
                         double depth_min, double depth_max);

private:
    VelocityModel1D velocity_model_;
    TravelTimeTable travel_times_;
    
    double min_cell_size_;
    int max_iterations_;
    
    // Search bounds
    double lat_min_, lat_max_;
    double lon_min_, lon_max_;
    double depth_min_, depth_max_;
    
    // Octree cell
    struct Cell {
        double lat_min, lat_max;
        double lon_min, lon_max;
        double depth_min, depth_max;
        double misfit;
        
        GeoPoint center() const {
            return GeoPoint((lat_min + lat_max) / 2,
                           (lon_min + lon_max) / 2,
                           (depth_min + depth_max) / 2);
        }
        
        double volume() const {
            return (lat_max - lat_min) * (lon_max - lon_min) * 
                   (depth_max - depth_min);
        }
    };
    
    // Subdivide cell into 8 children
    std::vector<Cell> subdivide(const Cell& parent) const;
    
    double computeMisfit(const Cell& cell,
                          const std::vector<PickPtr>& picks,
                          const StationInventory& stations);
};

} // namespace seisproc
