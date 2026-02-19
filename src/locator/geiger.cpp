#include "realdetect/locator/geiger.hpp"
#include "realdetect/locator/grid_search.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <set>

namespace realdetect {

GeigerLocator::GeigerLocator()
    : max_iterations_(20)
    , convergence_thresh_(0.001)
    , damping_factor_(0.5)
    , fixed_depth_(false)
    , default_depth_(10.0)
{
    velocity_model_ = VelocityModel1D::simpleThreeLayer();
    travel_times_.initialize(velocity_model_, 1000.0, 100.0);
}

void GeigerLocator::setParameter(const std::string& name, double value) {
    if (name == "max_iterations") max_iterations_ = static_cast<int>(value);
    else if (name == "convergence") convergence_thresh_ = value;
    else if (name == "damping") damping_factor_ = value;
    else if (name == "fixed_depth") fixed_depth_ = (value != 0);
    else if (name == "default_depth") default_depth_ = value;
}

void GeigerLocator::setVelocityModel(const VelocityModel1D& model) {
    velocity_model_ = model;
    travel_times_.initialize(model, 3000.0, 100.0);
}

LocationResult GeigerLocator::locate(const std::vector<PickPtr>& picks,
                                       const StationInventory& stations) {
    LocationResult result;
    result.algorithm = "Geiger";
    
    if (picks.size() < 3) {
        result.converged = false;
        return result;
    }
    
    // Get initial location
    GeoPoint location = getInitialLocation(picks, stations);
    TimePoint origin_time = getInitialTime(location, picks, stations);
    
    int n_params = fixed_depth_ ? 3 : 4;  // x, y, [z], t
    int n_obs = picks.size();
    
    result.initial_rms = 0;
    
    // Iterative inversion
    bool converged = false;
    int iter = 0;
    double prev_rms = 1e10;
    
    for (; iter < max_iterations_ && !converged; iter++) {
        // Build linear system
        Eigen::MatrixXd G(n_obs, n_params);
        Eigen::VectorXd d(n_obs);
        Eigen::VectorXd w(n_obs);
        
        buildSystem(location, origin_time, picks, stations, G, d, w);
        
        // Compute RMS
        double rms = computeRMS(d, w);
        if (iter == 0) result.initial_rms = rms;
        
        // Check convergence
        if (std::abs(prev_rms - rms) < convergence_thresh_) {
            converged = true;
            break;
        }
        prev_rms = rms;
        
        // Solve for model update
        Eigen::VectorXd dm = solveSystem(G, d, w);
        
        // Apply update with damping
        double dx = dm(0) * damping_factor_;  // km
        double dy = dm(1) * damping_factor_;  // km
        double dt = dm(n_params - 1) * damping_factor_;  // seconds
        
        // Convert dx, dy to lat/lon
        double dlat = dy / 111.0;
        double dlon = dx / (111.0 * std::cos(location.latitude * M_PI / 180.0));
        
        location.latitude += dlat;
        location.longitude += dlon;
        
        if (!fixed_depth_) {
            double dz = dm(2) * damping_factor_;
            location.depth += dz;
            location.depth = std::max(0.0, std::min(100.0, location.depth));
        }
        
        origin_time += std::chrono::milliseconds(static_cast<int>(dt * 1000));
        
        // Check for large jumps (potential divergence)
        if (std::abs(dx) > 500 || std::abs(dy) > 500) {
            // Reset damping
            damping_factor_ *= 0.5;
        }
    }
    
    // Build final result
    result.converged = converged;
    result.iterations = iter;
    result.origin.location = location;
    result.origin.time = origin_time;
    result.origin.is_fixed_depth = fixed_depth_;
    result.origin.algorithm = "Geiger";
    
    // Compute final RMS and arrivals
    Eigen::MatrixXd G(n_obs, n_params);
    Eigen::VectorXd d(n_obs);
    Eigen::VectorXd w(n_obs);
    buildSystem(location, origin_time, picks, stations, G, d, w);
    
    result.final_rms = computeRMS(d, w);
    result.origin.rms = result.final_rms;
    
    // Compute uncertainties
    computeUncertainties(G, w, result.origin.rms, result.origin);
    
    // Add arrivals
    std::set<std::string> unique_stations;
    double min_az = 360, max_az = 0;
    
    for (size_t i = 0; i < picks.size(); i++) {
        const auto& pick = picks[i];
        auto sta = stations.getStation(pick->stream_id);
        if (!sta) continue;
        
        Arrival arr;
        arr.pick = pick;
        arr.distance = location.distanceTo(sta->location());
        arr.azimuth = location.azimuthTo(sta->location());
        arr.residual = d(i);
        arr.weight = w(i);
        arr.used = true;
        
        result.origin.arrivals.push_back(arr);
        unique_stations.insert(pick->stream_id.station);
        
        min_az = std::min(min_az, arr.azimuth);
        max_az = std::max(max_az, arr.azimuth);
    }
    
    result.origin.phase_count = picks.size();
    result.origin.station_count = unique_stations.size();
    result.origin.gap = 360.0 - (max_az - min_az);
    if (result.origin.gap < 0) result.origin.gap += 360;
    
    return result;
}

LocationResult GeigerLocator::relocate(const EventPtr& event,
                                         const StationInventory& stations) {
    // Set initial location from existing origin
    const Origin& prev = event->preferredOrigin();
    has_initial_location_ = true;
    initial_location_ = prev.location;
    has_initial_time_ = true;
    initial_time_ = prev.time;
    
    std::vector<PickPtr> picks;
    for (const auto& arr : prev.arrivals) {
        picks.push_back(arr.pick);
    }
    
    return locate(picks, stations);
}

GeigerLocator::TravelTimeDerivatives GeigerLocator::computeDerivatives(
    const GeoPoint& source, const GeoPoint& station, PhaseType phase) const {
    
    TravelTimeDerivatives result;
    
    double dist = source.distanceTo(station);
    double az = source.azimuthTo(station) * M_PI / 180.0;
    
    // Travel time
    result.time = travel_times_.getTime(dist, source.depth, phase);
    
    // Partial derivatives (numerical)
    double dd = 0.1;  // km perturbation
    
    double tt_plus_x = travel_times_.getTime(dist + dd * std::sin(az), 
                                               source.depth, phase);
    double tt_minus_x = travel_times_.getTime(dist - dd * std::sin(az),
                                                source.depth, phase);
    result.dt_dx = (tt_plus_x - tt_minus_x) / (2 * dd);
    
    double tt_plus_y = travel_times_.getTime(dist + dd * std::cos(az),
                                               source.depth, phase);
    double tt_minus_y = travel_times_.getTime(dist - dd * std::cos(az),
                                                source.depth, phase);
    result.dt_dy = (tt_plus_y - tt_minus_y) / (2 * dd);
    
    double dz = 0.5;  // km perturbation
    double tt_plus_z = travel_times_.getTime(dist, source.depth + dz, phase);
    double tt_minus_z = travel_times_.getTime(dist, 
                                                std::max(0.0, source.depth - dz), phase);
    result.dt_dz = (tt_plus_z - tt_minus_z) / (2 * dz);
    
    return result;
}

void GeigerLocator::buildSystem(const GeoPoint& location, TimePoint origin_time,
                                  const std::vector<PickPtr>& picks,
                                  const StationInventory& stations,
                                  Eigen::MatrixXd& G, Eigen::VectorXd& d,
                                  Eigen::VectorXd& w) {
    int n = picks.size();
    int m = fixed_depth_ ? 3 : 4;
    
    G.resize(n, m);
    d.resize(n);
    w.resize(n);
    
    for (int i = 0; i < n; i++) {
        const auto& pick = picks[i];
        auto sta = stations.getStation(pick->stream_id);
        
        if (!sta) {
            G.row(i).setZero();
            d(i) = 0;
            w(i) = 0;
            continue;
        }
        
        // Compute travel time and derivatives
        auto deriv = computeDerivatives(location, sta->location(), pick->phase_type);
        
        // Fill G matrix (Jacobian)
        G(i, 0) = deriv.dt_dx;  // dx (east)
        G(i, 1) = deriv.dt_dy;  // dy (north)
        
        if (!fixed_depth_) {
            G(i, 2) = deriv.dt_dz;  // dz (depth)
            G(i, 3) = -1.0;         // dt (origin time)
        } else {
            G(i, 2) = -1.0;         // dt (origin time)
        }
        
        // Compute residual
        auto obs_dt = std::chrono::duration_cast<std::chrono::microseconds>(
            pick->time - origin_time);
        double obs_tt = obs_dt.count() / 1e6;
        
        d(i) = obs_tt - deriv.time;
        
        // Weight
        w(i) = pick->locationWeight();
    }
}

Eigen::VectorXd GeigerLocator::solveSystem(const Eigen::MatrixXd& G,
                                             const Eigen::VectorXd& d,
                                             const Eigen::VectorXd& w) {
    // Weighted least squares: (G'WG + lambda*I)^-1 * G'Wd
    int m = G.cols();
    
    // Weight matrix (create diagonal matrix from weights vector)
    Eigen::MatrixXd W = Eigen::MatrixXd::asDiagonal(w);
    
    // Normal equations
    Eigen::MatrixXd GtWG = G.transpose() * W * G;
    
    // Add damping (Levenberg-Marquardt)
    Eigen::VectorXd diag = GtWG.diagonal();
    double mean_diag = 0;
    for (int i = 0; i < m; i++) mean_diag += diag(i);
    mean_diag /= m;
    double lambda = damping_factor_ * mean_diag;
    GtWG.addToDiagonal(Eigen::VectorXd::Constant(m, lambda));
    
    Eigen::VectorXd GtWd = G.transpose() * W * d;
    
    // Solve
    return GtWG.ldlt().solve(GtWd);
}

double GeigerLocator::computeRMS(const Eigen::VectorXd& residuals,
                                   const Eigen::VectorXd& weights) {
    double sum_sq = 0;
    double sum_w = 0;
    
    for (int i = 0; i < residuals.size(); i++) {
        if (weights(i) > 0) {
            sum_sq += residuals(i) * residuals(i) * weights(i);
            sum_w += weights(i);
        }
    }
    
    return sum_w > 0 ? std::sqrt(sum_sq / sum_w) : 0;
}

void GeigerLocator::computeUncertainties(const Eigen::MatrixXd& G,
                                           const Eigen::VectorXd& w,
                                           double rms,
                                           Origin& origin) {
    // Covariance matrix: sigma^2 * (G'WG)^-1
    Eigen::MatrixXd W = Eigen::MatrixXd::asDiagonal(w);
    Eigen::MatrixXd GtWG = G.transpose() * W * G;
    
    Eigen::MatrixXd cov;
    try {
        Eigen::MatrixXd inv = GtWG.inverse();
        double scale = rms * rms;
        cov = Eigen::MatrixXd(inv.rows(), inv.cols());
        for (int i = 0; i < inv.rows(); i++) {
            for (int j = 0; j < inv.cols(); j++) {
                cov(i, j) = scale * inv(i, j);
            }
        }
    } catch (...) {
        origin.latitude_error = 10.0;
        origin.longitude_error = 10.0;
        origin.depth_error = 10.0;
        origin.time_error = 1.0;
        return;
    }
    
    // Extract uncertainties
    origin.longitude_error = std::sqrt(std::abs(cov(0, 0)));  // km
    origin.latitude_error = std::sqrt(std::abs(cov(1, 1)));   // km
    
    if (!fixed_depth_) {
        origin.depth_error = std::sqrt(std::abs(cov(2, 2)));  // km
        origin.time_error = std::sqrt(std::abs(cov(3, 3)));   // seconds
    } else {
        origin.depth_error = 0;
        origin.time_error = std::sqrt(std::abs(cov(2, 2)));
    }
}

GeoPoint GeigerLocator::getInitialLocation(const std::vector<PickPtr>& picks,
                                             const StationInventory& stations) {
    if (has_initial_location_) {
        return initial_location_;
    }
    
    // Use station centroid as initial guess
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

TimePoint GeigerLocator::getInitialTime(const GeoPoint& location,
                                          const std::vector<PickPtr>& picks,
                                          const StationInventory& stations) {
    if (has_initial_time_) {
        return initial_time_;
    }
    
    // Find earliest P arrival and estimate origin time
    for (const auto& pick : picks) {
        if (pick->phase_type != PhaseType::P && pick->phase_type != PhaseType::Pg) {
            continue;
        }
        
        auto sta = stations.getStation(pick->stream_id);
        if (!sta) continue;
        
        double dist = location.distanceTo(sta->location());
        double tt = travel_times_.getTime(dist, location.depth, PhaseType::P);
        
        return pick->time - std::chrono::milliseconds(static_cast<int>(tt * 1000));
    }
    
    // Fallback: use first pick
    if (!picks.empty()) {
        return picks[0]->time - std::chrono::seconds(10);
    }
    
    return std::chrono::system_clock::now();
}

// NonLinLocLocator - placeholder implementation

NonLinLocLocator::NonLinLocLocator()
    : grid_spacing_(5.0)
    , search_radius_(200.0)
{
}

void NonLinLocLocator::setParameter(const std::string& name, double value) {
    if (name == "grid_spacing") grid_spacing_ = value;
    else if (name == "search_radius") search_radius_ = value;
}

void NonLinLocLocator::setVelocityModel(const VelocityModel1D& model) {
    velocity_model_ = model;
    travel_times_.initialize(model, 3000.0, 100.0);
}

LocationResult NonLinLocLocator::locate(const std::vector<PickPtr>& picks,
                                           const StationInventory& stations) {
    // Use OctTree + Geiger for now
    OctTreeLocator oct;
    oct.setVelocityModel(velocity_model_);
    LocationResult oct_result = oct.locate(picks, stations);
    
    // Refine with Geiger
    GeigerLocator geiger;
    geiger.setVelocityModel(velocity_model_);
    geiger.setInitialLocation(oct_result.origin.location);
    geiger.setInitialTime(oct_result.origin.time);
    
    LocationResult result = geiger.locate(picks, stations);
    result.algorithm = "NonLinLoc";
    
    return result;
}

LocationResult NonLinLocLocator::relocate(const EventPtr& event,
                                            const StationInventory& stations) {
    std::vector<PickPtr> picks;
    for (const auto& arr : event->preferredOrigin().arrivals) {
        picks.push_back(arr.pick);
    }
    return locate(picks, stations);
}

double NonLinLocLocator::computeLikelihood(const GeoPoint& point, 
                                             TimePoint origin_time,
                                             const std::vector<PickPtr>& picks,
                                             const StationInventory& stations) {
    double sum_sq = 0;
    
    for (const auto& pick : picks) {
        auto sta = stations.getStation(pick->stream_id);
        if (!sta) continue;
        
        double dist = point.distanceTo(sta->location());
        double expected_tt = travel_times_.getTime(dist, point.depth, pick->phase_type);
        
        auto obs_dt = std::chrono::duration_cast<std::chrono::milliseconds>(
            pick->time - origin_time);
        double obs_tt = obs_dt.count() / 1000.0;
        
        double residual = obs_tt - expected_tt;
        double sigma = pick->uncertainty > 0 ? pick->uncertainty : 0.5;
        
        sum_sq += (residual * residual) / (sigma * sigma);
    }
    
    return std::exp(-0.5 * sum_sq);
}

double NonLinLocLocator::edtLikelihood(const GeoPoint& point,
                                         const std::vector<PickPtr>& picks,
                                         const StationInventory& stations) {
    // Equal Differential Time (EDT) formulation
    double sum_sq = 0;
    int count = 0;
    
    for (size_t i = 0; i < picks.size(); i++) {
        auto sta_i = stations.getStation(picks[i]->stream_id);
        if (!sta_i) continue;
        
        double dist_i = point.distanceTo(sta_i->location());
        double tt_i = travel_times_.getTime(dist_i, point.depth, picks[i]->phase_type);
        
        for (size_t j = i + 1; j < picks.size(); j++) {
            auto sta_j = stations.getStation(picks[j]->stream_id);
            if (!sta_j) continue;
            
            double dist_j = point.distanceTo(sta_j->location());
            double tt_j = travel_times_.getTime(dist_j, point.depth, picks[j]->phase_type);
            
            // Observed differential time
            auto dt_obs = std::chrono::duration_cast<std::chrono::milliseconds>(
                picks[j]->time - picks[i]->time);
            double obs_diff = dt_obs.count() / 1000.0;
            
            // Calculated differential time
            double calc_diff = tt_j - tt_i;
            
            double residual = obs_diff - calc_diff;
            sum_sq += residual * residual;
            count++;
        }
    }
    
    if (count == 0) return 0;
    return std::exp(-0.5 * sum_sq / count);
}

// LocatorFactory

LocatorPtr LocatorFactory::create(const std::string& algorithm) {
    if (algorithm == "GridSearch" || algorithm == "grid") {
        return std::make_shared<GridSearchLocator>();
    } else if (algorithm == "OctTree" || algorithm == "octree") {
        return std::make_shared<OctTreeLocator>();
    } else if (algorithm == "Geiger" || algorithm == "geiger") {
        return std::make_shared<GeigerLocator>();
    } else if (algorithm == "NonLinLoc" || algorithm == "nll") {
        return std::make_shared<NonLinLocLocator>();
    }
    
    // Default to Geiger
    return std::make_shared<GeigerLocator>();
}

std::vector<std::string> LocatorFactory::availableAlgorithms() {
    return {"GridSearch", "OctTree", "Geiger", "NonLinLoc"};
}

} // namespace realdetect
