#include "realdetect/locator/grid_search.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <set>

namespace realdetect {

GridSearchLocator::GridSearchLocator()
    : h_step_(5.0)
    , z_step_(5.0)
    , search_radius_(200.0)
    , depth_min_(0.0)
    , depth_max_(50.0)
    , fixed_depth_(false)
    , default_depth_(10.0)
    , refinement_levels_(3)
{
    velocity_model_ = VelocityModel1D::simpleThreeLayer();
    travel_times_.initialize(velocity_model_, 1000.0, 100.0);
}

void GridSearchLocator::setParameter(const std::string& name, double value) {
    if (name == "h_step") h_step_ = value;
    else if (name == "z_step") z_step_ = value;
    else if (name == "search_radius") search_radius_ = value;
    else if (name == "depth_min") depth_min_ = value;
    else if (name == "depth_max") depth_max_ = value;
    else if (name == "fixed_depth") fixed_depth_ = (value != 0);
    else if (name == "default_depth") default_depth_ = value;
}

void GridSearchLocator::setVelocityModel(const VelocityModel1D& model) {
    velocity_model_ = model;
    travel_times_.initialize(model, 1000.0, 100.0);
}

LocationResult GridSearchLocator::locate(const std::vector<PickPtr>& picks,
                                           const StationInventory& stations) {
    LocationResult result;
    result.algorithm = "GridSearch";
    
    if (picks.size() < 3) {
        result.converged = false;
        return result;
    }
    
    // Get initial estimate from station centroid
    GeoPoint initial = estimateInitialLocation(picks, stations);
    
    double current_radius = search_radius_;
    double current_h_step = h_step_;
    double current_z_step = z_step_;
    
    GeoPoint best_location = initial;
    TimePoint best_origin;
    double best_rms = std::numeric_limits<double>::max();
    
    // Multi-level refinement
    for (int level = 0; level < refinement_levels_; level++) {
        double lat_min = best_location.latitude - current_radius / 111.0;
        double lat_max = best_location.latitude + current_radius / 111.0;
        double lon_min = best_location.longitude - current_radius / 
                         (111.0 * std::cos(best_location.latitude * M_PI / 180.0));
        double lon_max = best_location.longitude + current_radius /
                         (111.0 * std::cos(best_location.latitude * M_PI / 180.0));
        
        double z_min = fixed_depth_ ? default_depth_ : depth_min_;
        double z_max = fixed_depth_ ? default_depth_ : depth_max_;
        
        // Grid search
        for (double lat = lat_min; lat <= lat_max; lat += current_h_step / 111.0) {
            for (double lon = lon_min; lon <= lon_max; 
                 lon += current_h_step / (111.0 * std::cos(lat * M_PI / 180.0))) {
                for (double depth = z_min; depth <= z_max; depth += current_z_step) {
                    
                    GeoPoint test_loc(lat, lon, depth);
                    TimePoint test_time = estimateOriginTime(test_loc, picks, stations);
                    
                    double rms = computeRMS(lat, lon, depth, test_time, picks, stations);
                    
                    if (rms < best_rms) {
                        best_rms = rms;
                        best_location = test_loc;
                        best_origin = test_time;
                    }
                }
            }
        }
        
        // Reduce search area for next level
        current_radius /= 3.0;
        current_h_step /= 3.0;
        current_z_step /= 2.0;
    }
    
    // Fill in result
    result.origin.location = best_location;
    result.origin.time = best_origin;
    result.origin.rms = best_rms;
    result.converged = (best_rms < 10.0);  // Reasonable threshold
    result.iterations = refinement_levels_;
    result.final_rms = best_rms;
    
    // Compute arrivals
    std::set<std::string> unique_stations;
    double min_az = 360, max_az = 0;
    
    for (const auto& pick : picks) {
        auto sta = stations.getStation(pick->stream_id);
        if (!sta) continue;
        
        Arrival arr;
        arr.pick = pick;
        arr.distance = best_location.distanceTo(sta->location());
        arr.azimuth = best_location.azimuthTo(sta->location());
        
        double expected_tt = travel_times_.getTime(arr.distance, best_location.depth,
                                                    pick->phase_type);
        auto obs_dt = std::chrono::duration_cast<std::chrono::milliseconds>(
            pick->time - best_origin);
        double obs_tt = obs_dt.count() / 1000.0;
        
        arr.residual = obs_tt - expected_tt;
        arr.weight = pick->locationWeight();
        arr.used = true;
        
        result.origin.arrivals.push_back(arr);
        unique_stations.insert(pick->stream_id.station);
        
        min_az = std::min(min_az, arr.azimuth);
        max_az = std::max(max_az, arr.azimuth);
    }
    
    result.origin.phase_count = picks.size();
    result.origin.station_count = unique_stations.size();
    result.origin.gap = 360.0 - (max_az - min_az);
    result.origin.algorithm = "GridSearch";
    result.origin.is_fixed_depth = fixed_depth_;
    
    // Estimate errors (rough approximation)
    result.origin.latitude_error = current_h_step;
    result.origin.longitude_error = current_h_step;
    result.origin.depth_error = fixed_depth_ ? 0 : current_z_step;
    result.origin.time_error = best_rms;
    
    return result;
}

LocationResult GridSearchLocator::relocate(const EventPtr& event,
                                            const StationInventory& stations) {
    std::vector<PickPtr> picks;
    for (const auto& arr : event->preferredOrigin().arrivals) {
        picks.push_back(arr.pick);
    }
    return locate(picks, stations);
}

double GridSearchLocator::computeRMS(double lat, double lon, double depth, 
                                       TimePoint origin_time,
                                       const std::vector<PickPtr>& picks,
                                       const StationInventory& stations) const {
    GeoPoint source(lat, lon, depth);
    double sum_sq = 0;
    double total_weight = 0;
    
    for (const auto& pick : picks) {
        auto sta = stations.getStation(pick->stream_id);
        if (!sta) continue;
        
        double dist = source.distanceTo(sta->location());
        double expected_tt = travel_times_.getTime(dist, depth, pick->phase_type);
        
        auto obs_dt = std::chrono::duration_cast<std::chrono::milliseconds>(
            pick->time - origin_time);
        double obs_tt = obs_dt.count() / 1000.0;
        
        double residual = obs_tt - expected_tt;
        double weight = pick->locationWeight();
        
        sum_sq += residual * residual * weight;
        total_weight += weight;
    }
    
    return total_weight > 0 ? std::sqrt(sum_sq / total_weight) : 1e10;
}

GeoPoint GridSearchLocator::estimateInitialLocation(const std::vector<PickPtr>& picks,
                                                      const StationInventory& stations) const {
    double sum_lat = 0, sum_lon = 0;
    int count = 0;
    
    for (const auto& pick : picks) {
        auto sta = stations.getStation(pick->stream_id);
        if (sta) {
            sum_lat += sta->latitude();
            sum_lon += sta->longitude();
            count++;
        }
    }
    
    if (count == 0) return GeoPoint(0, 0, default_depth_);
    
    return GeoPoint(sum_lat / count, sum_lon / count, default_depth_);
}

TimePoint GridSearchLocator::estimateOriginTime(const GeoPoint& location,
                                                   const std::vector<PickPtr>& picks,
                                                   const StationInventory& stations) const {
    // Find earliest P arrival and estimate origin time
    TimePoint earliest;
    double min_tt = std::numeric_limits<double>::max();
    bool found = false;
    
    for (const auto& pick : picks) {
        if (pick->phase_type != PhaseType::P && pick->phase_type != PhaseType::Pg) {
            continue;
        }
        
        auto sta = stations.getStation(pick->stream_id);
        if (!sta) continue;
        
        double dist = location.distanceTo(sta->location());
        double tt = travel_times_.getTime(dist, location.depth, PhaseType::P);
        
        TimePoint est_origin = pick->time - 
            std::chrono::milliseconds(static_cast<int>(tt * 1000));
        
        if (!found || tt < min_tt) {
            min_tt = tt;
            earliest = est_origin;
            found = true;
        }
    }
    
    return earliest;
}

// OctTreeLocator implementation

OctTreeLocator::OctTreeLocator()
    : min_cell_size_(1.0)
    , max_iterations_(1000)
    , lat_min_(-90), lat_max_(90)
    , lon_min_(-180), lon_max_(180)
    , depth_min_(0), depth_max_(100)
{
}

void OctTreeLocator::setParameter(const std::string& name, double value) {
    if (name == "min_cell_size") min_cell_size_ = value;
    else if (name == "max_iterations") max_iterations_ = static_cast<int>(value);
}

void OctTreeLocator::setVelocityModel(const VelocityModel1D& model) {
    velocity_model_ = model;
    travel_times_.initialize(model, 1000.0, 100.0);
}

void OctTreeLocator::setSearchBounds(double lat_min, double lat_max,
                                       double lon_min, double lon_max,
                                       double depth_min, double depth_max) {
    lat_min_ = lat_min; lat_max_ = lat_max;
    lon_min_ = lon_min; lon_max_ = lon_max;
    depth_min_ = depth_min; depth_max_ = depth_max;
}

LocationResult OctTreeLocator::locate(const std::vector<PickPtr>& picks,
                                        const StationInventory& stations) {
    LocationResult result;
    result.algorithm = "OctTree";
    
    if (picks.size() < 3) {
        result.converged = false;
        return result;
    }
    
    // Initialize with root cell
    Cell root;
    root.lat_min = lat_min_; root.lat_max = lat_max_;
    root.lon_min = lon_min_; root.lon_max = lon_max_;
    root.depth_min = depth_min_; root.depth_max = depth_max_;
    root.misfit = computeMisfit(root, picks, stations);
    
    std::vector<Cell> cells = {root};
    Cell best_cell = root;
    
    for (int iter = 0; iter < max_iterations_ && !cells.empty(); iter++) {
        // Find cell with lowest misfit
        auto it = std::min_element(cells.begin(), cells.end(),
            [](const Cell& a, const Cell& b) { return a.misfit < b.misfit; });
        
        Cell current = *it;
        cells.erase(it);
        
        if (current.misfit < best_cell.misfit) {
            best_cell = current;
        }
        
        // Check if cell is small enough
        double cell_size = std::max({
            (current.lat_max - current.lat_min) * 111.0,
            (current.lon_max - current.lon_min) * 111.0 * 
                std::cos((current.lat_min + current.lat_max) / 2 * M_PI / 180),
            current.depth_max - current.depth_min
        });
        
        if (cell_size < min_cell_size_) {
            continue;  // Don't subdivide further
        }
        
        // Subdivide
        auto children = subdivide(current);
        for (auto& child : children) {
            child.misfit = computeMisfit(child, picks, stations);
            cells.push_back(child);
        }
    }
    
    // Build result from best cell
    GeoPoint location = best_cell.center();
    GridSearchLocator gs;
    gs.setVelocityModel(velocity_model_);
    gs.setHorizontalStep(min_cell_size_ / 3);
    gs.setDepthStep(min_cell_size_ / 3);
    gs.setSearchRadius(min_cell_size_ * 2);
    
    result = gs.locate(picks, stations);
    result.algorithm = "OctTree";
    
    return result;
}

LocationResult OctTreeLocator::relocate(const EventPtr& event,
                                          const StationInventory& stations) {
    std::vector<PickPtr> picks;
    for (const auto& arr : event->preferredOrigin().arrivals) {
        picks.push_back(arr.pick);
    }
    return locate(picks, stations);
}

std::vector<OctTreeLocator::Cell> OctTreeLocator::subdivide(const Cell& parent) const {
    std::vector<Cell> children;
    
    double lat_mid = (parent.lat_min + parent.lat_max) / 2;
    double lon_mid = (parent.lon_min + parent.lon_max) / 2;
    double depth_mid = (parent.depth_min + parent.depth_max) / 2;
    
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            for (int k = 0; k < 2; k++) {
                Cell child;
                child.lat_min = (i == 0) ? parent.lat_min : lat_mid;
                child.lat_max = (i == 0) ? lat_mid : parent.lat_max;
                child.lon_min = (j == 0) ? parent.lon_min : lon_mid;
                child.lon_max = (j == 0) ? lon_mid : parent.lon_max;
                child.depth_min = (k == 0) ? parent.depth_min : depth_mid;
                child.depth_max = (k == 0) ? depth_mid : parent.depth_max;
                children.push_back(child);
            }
        }
    }
    
    return children;
}

double OctTreeLocator::computeMisfit(const Cell& cell,
                                        const std::vector<PickPtr>& picks,
                                        const StationInventory& stations) {
    GeoPoint center = cell.center();
    
    // Estimate origin time from P arrivals
    TimePoint origin;
    bool has_origin = false;
    
    for (const auto& pick : picks) {
        if (pick->phase_type == PhaseType::P || pick->phase_type == PhaseType::Pg) {
            auto sta = stations.getStation(pick->stream_id);
            if (!sta) continue;
            
            double dist = center.distanceTo(sta->location());
            double tt = travel_times_.getTime(dist, center.depth, PhaseType::P);
            origin = pick->time - std::chrono::milliseconds(static_cast<int>(tt * 1000));
            has_origin = true;
            break;
        }
    }
    
    if (!has_origin) return 1e10;
    
    // Compute RMS misfit
    double sum_sq = 0;
    int count = 0;
    
    for (const auto& pick : picks) {
        auto sta = stations.getStation(pick->stream_id);
        if (!sta) continue;
        
        double dist = center.distanceTo(sta->location());
        double expected_tt = travel_times_.getTime(dist, center.depth, pick->phase_type);
        
        auto obs_dt = std::chrono::duration_cast<std::chrono::milliseconds>(
            pick->time - origin);
        double obs_tt = obs_dt.count() / 1000.0;
        
        double residual = obs_tt - expected_tt;
        sum_sq += residual * residual;
        count++;
    }
    
    return count > 0 ? std::sqrt(sum_sq / count) : 1e10;
}

} // namespace realdetect
