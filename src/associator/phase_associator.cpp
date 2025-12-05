#include "realdetect/associator/phase_associator.hpp"
#include <algorithm>
#include <set>
#include <cmath>
#include <iostream>

namespace realdetect {

// TravelTimeTable implementation

void TravelTimeTable::initialize(const VelocityModel1D& model,
                                   double max_distance, double max_depth) {
    max_dist_ = max_distance;
    max_depth_ = max_depth;
    
    size_t nd = static_cast<size_t>(max_distance / dist_step_) + 1;
    size_t nh = static_cast<size_t>(max_depth / depth_step_) + 1;
    
    p_times_.resize(nd, std::vector<double>(nh, 0));
    s_times_.resize(nd, std::vector<double>(nh, 0));
    
    for (size_t i = 0; i < nd; i++) {
        double dist = i * dist_step_;
        for (size_t j = 0; j < nh; j++) {
            double depth = j * depth_step_;
            p_times_[i][j] = model.travelTime(dist, depth, PhaseType::P);
            s_times_[i][j] = model.travelTime(dist, depth, PhaseType::S);
        }
    }
}

double TravelTimeTable::getTime(double distance_km, double depth_km, 
                                  PhaseType phase) const {
    if (p_times_.empty()) return 0;
    
    size_t di = static_cast<size_t>(distance_km / dist_step_);
    size_t hi = static_cast<size_t>(depth_km / depth_step_);
    
    di = std::min(di, p_times_.size() - 1);
    hi = std::min(hi, p_times_[0].size() - 1);
    
    bool is_p = (phase == PhaseType::P || phase == PhaseType::Pn || 
                 phase == PhaseType::Pg);
    
    return is_p ? p_times_[di][hi] : s_times_[di][hi];
}

double TravelTimeTable::getDTDD(double distance_km, double depth_km, 
                                  PhaseType phase) const {
    double t1 = getTime(distance_km, depth_km, phase);
    double t2 = getTime(distance_km + dist_step_, depth_km, phase);
    return (t2 - t1) / dist_step_;
}

double TravelTimeTable::getDTDH(double distance_km, double depth_km, 
                                  PhaseType phase) const {
    double t1 = getTime(distance_km, depth_km, phase);
    double t2 = getTime(distance_km, depth_km + depth_step_, phase);
    return (t2 - t1) / depth_step_;
}

// PhaseAssociator implementation

PhaseAssociator::PhaseAssociator()
    : stations_(nullptr)
    , time_window_(60.0)
    , max_residual_(3.0)
    , min_stations_(3)
    , min_phases_(4)
    , search_radius_(500.0)
{
}

void PhaseAssociator::setVelocityModel(const VelocityModel1D& model) {
    velocity_model_ = model;
    travel_times_.initialize(model, 1000.0, 100.0);
}

void PhaseAssociator::addPick(PickPtr pick) {
    pending_picks_.push_back(pick);
}

std::vector<EventPtr> PhaseAssociator::process() {
    std::vector<EventPtr> events;
    
    if (pending_picks_.size() < static_cast<size_t>(min_phases_)) {
        return events;
    }
    
    // Sort picks by time
    std::sort(pending_picks_.begin(), pending_picks_.end(),
              [](const PickPtr& a, const PickPtr& b) {
                  return a->time < b->time;
              });
    
    // Find candidates
    std::vector<AssociationCandidate> candidates = findCandidates();
    
    // Build events from candidates
    for (auto& candidate : candidates) {
        if (candidate.station_count >= min_stations_ &&
            candidate.phase_count >= min_phases_) {
            EventPtr event = buildEvent(candidate);
            if (event) {
                events.push_back(event);
                
                if (event_callback_) {
                    event_callback_(event);
                }
                
                // Remove associated picks from pending
                std::set<PickPtr> associated(candidate.picks.begin(), 
                                              candidate.picks.end());
                pending_picks_.erase(
                    std::remove_if(pending_picks_.begin(), pending_picks_.end(),
                                   [&associated](const PickPtr& p) {
                                       return associated.count(p) > 0;
                                   }),
                    pending_picks_.end());
            }
        }
    }
    
    // Move old unassociated picks
    auto now = std::chrono::system_clock::now();
    auto cutoff = now - std::chrono::seconds(static_cast<int>(time_window_ * 2));
    
    for (auto it = pending_picks_.begin(); it != pending_picks_.end(); ) {
        if ((*it)->time < cutoff) {
            unassociated_picks_.push_back(*it);
            it = pending_picks_.erase(it);
        } else {
            ++it;
        }
    }
    
    return events;
}

std::vector<AssociationCandidate> PhaseAssociator::findCandidates() {
    std::vector<AssociationCandidate> candidates;
    
    if (!stations_ || pending_picks_.empty()) return candidates;
    
    // Group picks by approximate origin time
    // Assume first P arrival is closest to origin
    
    for (size_t i = 0; i < pending_picks_.size(); i++) {
        const auto& first_pick = pending_picks_[i];
        
        // Skip if not P-wave
        if (first_pick->phase_type != PhaseType::P &&
            first_pick->phase_type != PhaseType::Pg) continue;
        
        // Get station location
        auto sta = stations_->getStation(first_pick->stream_id);
        if (!sta) continue;
        
        AssociationCandidate candidate;
        candidate.location = sta->location();
        candidate.location.depth = 10.0;  // Default depth
        candidate.picks.push_back(first_pick);
        
        // Estimate origin time (assume 8 km/s P velocity for now)
        double assumed_dist = 50.0;  // km
        double travel_time = assumed_dist / 8.0;
        candidate.origin_time = first_pick->time - 
            std::chrono::milliseconds(static_cast<int>(travel_time * 1000));
        
        // Find other picks that could belong to this event
        for (size_t j = i + 1; j < pending_picks_.size(); j++) {
            if (tryAssociate(pending_picks_[j], candidate)) {
                candidate.picks.push_back(pending_picks_[j]);
            }
        }
        
        // Count unique stations
        std::set<std::string> unique_stations;
        for (const auto& p : candidate.picks) {
            unique_stations.insert(p->stream_id.station);
        }
        candidate.station_count = unique_stations.size();
        candidate.phase_count = candidate.picks.size();
        
        if (candidate.station_count >= min_stations_) {
            candidates.push_back(candidate);
        }
    }
    
    return candidates;
}

bool PhaseAssociator::tryAssociate(const PickPtr& pick, 
                                    AssociationCandidate& candidate) {
    double residual = computeResidual(pick, candidate);
    return std::abs(residual) <= max_residual_;
}

double PhaseAssociator::computeResidual(const PickPtr& pick,
                                          const AssociationCandidate& candidate) {
    auto sta = stations_->getStation(pick->stream_id);
    if (!sta) return 1e10;
    
    // Distance from candidate location to station
    double dist = candidate.location.distanceTo(sta->location());
    
    // Expected travel time
    double tt = travel_times_.getTime(dist, candidate.location.depth, 
                                       pick->phase_type);
    
    // Observed travel time
    auto obs_dt = std::chrono::duration_cast<std::chrono::milliseconds>(
        pick->time - candidate.origin_time);
    double obs_tt = obs_dt.count() / 1000.0;
    
    return obs_tt - tt;
}

EventPtr PhaseAssociator::buildEvent(const AssociationCandidate& candidate) {
    auto event = std::make_shared<Event>();
    
    Origin origin;
    origin.time = candidate.origin_time;
    origin.location = candidate.location;
    origin.phase_count = candidate.phase_count;
    origin.station_count = candidate.station_count;
    
    // Create arrivals
    std::set<std::string> stations_used;
    double sum_sq_res = 0;
    double min_az = 360, max_az = 0;
    
    for (const auto& pick : candidate.picks) {
        auto sta = stations_->getStation(pick->stream_id);
        if (!sta) continue;
        
        Arrival arr;
        arr.pick = pick;
        arr.distance = candidate.location.distanceTo(sta->location());
        arr.azimuth = candidate.location.azimuthTo(sta->location());
        arr.residual = computeResidual(pick, candidate);
        arr.weight = pick->locationWeight();
        arr.used = true;
        
        sum_sq_res += arr.residual * arr.residual * arr.weight;
        origin.arrivals.push_back(arr);
        stations_used.insert(pick->stream_id.station);
        
        // Track azimuthal coverage
        min_az = std::min(min_az, arr.azimuth);
        max_az = std::max(max_az, arr.azimuth);
    }
    
    // Compute RMS
    double total_weight = 0;
    for (const auto& arr : origin.arrivals) {
        total_weight += arr.weight;
    }
    origin.rms = std::sqrt(sum_sq_res / total_weight);
    
    // Compute gap (simplified)
    origin.gap = 360.0 - (max_az - min_az);
    if (origin.gap < 0) origin.gap += 360;
    
    origin.algorithm = "Associator";
    
    event->addOrigin(origin);
    
    return event;
}

// NucleatorAssociator implementation

NucleatorAssociator::NucleatorAssociator()
    : stations_(nullptr)
    , time_window_(60.0)
    , min_phases_(4)
    , max_residual_(3.0)
{
}

void NucleatorAssociator::setVelocityModel(const VelocityModel1D& model) {
    velocity_model_ = model;
    travel_times_.initialize(model, 1000.0, 100.0);
}

void NucleatorAssociator::setGrid(double lat_min, double lat_max, double lat_step,
                                    double lon_min, double lon_max, double lon_step,
                                    double depth_min, double depth_max, double depth_step) {
    lat_min_ = lat_min; lat_max_ = lat_max; lat_step_ = lat_step;
    lon_min_ = lon_min; lon_max_ = lon_max; lon_step_ = lon_step;
    depth_min_ = depth_min; depth_max_ = depth_max; depth_step_ = depth_step;
    
    initializeGrid();
}

void NucleatorAssociator::initializeGrid() {
    grid_points_.clear();
    
    for (double lat = lat_min_; lat <= lat_max_; lat += lat_step_) {
        for (double lon = lon_min_; lon <= lon_max_; lon += lon_step_) {
            for (double dep = depth_min_; dep <= depth_max_; dep += depth_step_) {
                GridPoint gp;
                gp.location = GeoPoint(lat, lon, dep);
                gp.vote_count = 0;
                grid_points_.push_back(gp);
            }
        }
    }
}

void NucleatorAssociator::addPick(PickPtr pick) {
    pending_picks_.push_back(pick);
}

std::vector<EventPtr> NucleatorAssociator::process() {
    // Reset votes
    for (auto& gp : grid_points_) {
        gp.vote_count = 0;
        gp.associated_picks.clear();
    }
    
    // Vote for each pick
    for (const auto& pick : pending_picks_) {
        voteForGrid(pick);
    }
    
    return extractEvents();
}

void NucleatorAssociator::voteForGrid(const PickPtr& pick) {
    if (!stations_) return;
    
    auto sta = stations_->getStation(pick->stream_id);
    if (!sta) return;
    
    for (auto& gp : grid_points_) {
        double dist = gp.location.distanceTo(sta->location());
        double tt = travel_times_.getTime(dist, gp.location.depth, pick->phase_type);
        
        // Check if pick time is consistent
        auto pick_dt = std::chrono::duration_cast<std::chrono::milliseconds>(
            pick->time - std::chrono::system_clock::now()).count() / 1000.0;
        
        // This is simplified - would need proper origin time estimation
        gp.vote_count++;
        gp.associated_picks.push_back(pick);
    }
}

std::vector<EventPtr> NucleatorAssociator::extractEvents() {
    std::vector<EventPtr> events;
    
    // Find grid points with enough votes
    for (const auto& gp : grid_points_) {
        if (gp.vote_count >= min_phases_) {
            auto event = std::make_shared<Event>();
            
            Origin origin;
            origin.location = gp.location;
            origin.phase_count = gp.vote_count;
            
            // Count unique stations
            std::set<std::string> stations;
            for (const auto& p : gp.associated_picks) {
                stations.insert(p->stream_id.station);
            }
            origin.station_count = stations.size();
            
            origin.algorithm = "Nucleator";
            event->addOrigin(origin);
            
            events.push_back(event);
        }
    }
    
    return events;
}

} // namespace realdetect
