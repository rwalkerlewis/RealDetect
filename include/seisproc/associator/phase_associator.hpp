#pragma once

#include "../core/types.hpp"
#include "../core/event.hpp"
#include "../core/station.hpp"
#include "../core/velocity_model.hpp"
#include <vector>
#include <memory>
#include <functional>

namespace seisproc {

/**
 * AssociationCandidate - Potential event from clustered picks
 */
struct AssociationCandidate {
    TimePoint origin_time;
    GeoPoint location;
    std::vector<PickPtr> picks;
    double score;
    int phase_count;
    int station_count;
    
    AssociationCandidate() : score(0), phase_count(0), station_count(0) {}
};

/**
 * TravelTimeTable - Pre-computed travel times
 */
class TravelTimeTable {
public:
    TravelTimeTable() = default;
    
    // Initialize with velocity model
    void initialize(const VelocityModel1D& model, 
                    double max_distance = 1000.0,
                    double max_depth = 100.0);
    
    // Get travel time (seconds)
    double getTime(double distance_km, double depth_km, PhaseType phase) const;
    
    // Get derivative d(time)/d(distance)
    double getDTDD(double distance_km, double depth_km, PhaseType phase) const;
    
    // Get derivative d(time)/d(depth)
    double getDTDH(double distance_km, double depth_km, PhaseType phase) const;

private:
    std::vector<std::vector<double>> p_times_;  // [distance][depth]
    std::vector<std::vector<double>> s_times_;
    double dist_step_ = 1.0;   // km
    double depth_step_ = 1.0;  // km
    double max_dist_ = 1000.0;
    double max_depth_ = 100.0;
};

/**
 * PhaseAssociator - Associate picks into events
 */
class PhaseAssociator {
public:
    PhaseAssociator();
    
    // Set station inventory
    void setStations(const StationInventory& stations) { stations_ = &stations; }
    
    // Set velocity model
    void setVelocityModel(const VelocityModel1D& model);
    
    // Add a new pick
    void addPick(PickPtr pick);
    
    // Process picks and return new events
    std::vector<EventPtr> process();
    
    // Get unassociated picks
    std::vector<PickPtr> getUnassociatedPicks() const { return unassociated_picks_; }
    
    // Configuration
    void setTimeWindow(double seconds) { time_window_ = seconds; }
    void setMaxResidual(double seconds) { max_residual_ = seconds; }
    void setMinStations(int count) { min_stations_ = count; }
    void setMinPhases(int count) { min_phases_ = count; }
    void setSearchRadius(double km) { search_radius_ = km; }
    
    // Callbacks
    using EventCallback = std::function<void(EventPtr)>;
    void setEventCallback(EventCallback cb) { event_callback_ = cb; }

private:
    const StationInventory* stations_;
    VelocityModel1D velocity_model_;
    TravelTimeTable travel_times_;
    
    std::vector<PickPtr> pending_picks_;
    std::vector<PickPtr> unassociated_picks_;
    
    // Association parameters
    double time_window_;    // Max time span for association
    double max_residual_;   // Max allowed travel time residual
    int min_stations_;      // Minimum stations for event
    int min_phases_;        // Minimum phases for event
    double search_radius_;  // Search radius (km)
    
    EventCallback event_callback_;
    
    // Internal methods
    std::vector<AssociationCandidate> findCandidates();
    bool tryAssociate(const PickPtr& pick, AssociationCandidate& candidate);
    double computeResidual(const PickPtr& pick, const AssociationCandidate& candidate);
    EventPtr buildEvent(const AssociationCandidate& candidate);
};

/**
 * NucleatorAssociator - Grid-based nucleation approach
 * 
 * Places virtual sources on a grid and finds consistent pick clusters
 */
class NucleatorAssociator {
public:
    NucleatorAssociator();
    
    void setStations(const StationInventory& stations) { stations_ = &stations; }
    void setVelocityModel(const VelocityModel1D& model);
    
    // Define grid
    void setGrid(double lat_min, double lat_max, double lat_step,
                 double lon_min, double lon_max, double lon_step,
                 double depth_min, double depth_max, double depth_step);
    
    // Add picks
    void addPick(PickPtr pick);
    
    // Process and return events
    std::vector<EventPtr> process();
    
    // Configuration
    void setTimeWindow(double seconds) { time_window_ = seconds; }
    void setMinPhases(int count) { min_phases_ = count; }
    void setMaxResidual(double seconds) { max_residual_ = seconds; }

private:
    const StationInventory* stations_;
    VelocityModel1D velocity_model_;
    TravelTimeTable travel_times_;
    
    // Grid parameters
    double lat_min_, lat_max_, lat_step_;
    double lon_min_, lon_max_, lon_step_;
    double depth_min_, depth_max_, depth_step_;
    
    std::vector<PickPtr> pending_picks_;
    
    double time_window_;
    int min_phases_;
    double max_residual_;
    
    // Grid-based search
    struct GridPoint {
        GeoPoint location;
        int vote_count;
        std::vector<PickPtr> associated_picks;
    };
    
    std::vector<GridPoint> grid_points_;
    
    void initializeGrid();
    void voteForGrid(const PickPtr& pick);
    std::vector<EventPtr> extractEvents();
};

} // namespace seisproc
