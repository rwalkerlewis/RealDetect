#pragma once

#include "locator.hpp"
#include <Eigen/Dense>

namespace realdetect {

/**
 * GeigerLocator - Iterative linearized inversion (Geiger's method)
 * 
 * The classic earthquake location algorithm using iterative least-squares.
 * 
 * We minimize the objective function:
 *   Phi = sum_i w_i * (t_obs_i - t_calc_i)^2
 * 
 * By Taylor expanding t_calc and solving the linear system:
 *   G * dm = d
 * 
 * Where:
 *   G = matrix of partial derivatives dt/d(lat,lon,depth,time)
 *   dm = model parameter updates
 *   d = residual vector (observed - calculated)
 */
class GeigerLocator : public BaseLocator {
public:
    GeigerLocator();
    
    LocationResult locate(const std::vector<PickPtr>& picks,
                           const StationInventory& stations) override;
    
    LocationResult relocate(const EventPtr& event,
                             const StationInventory& stations) override;
    
    std::string name() const override { return "Geiger"; }
    
    void setParameter(const std::string& name, double value) override;
    void setVelocityModel(const VelocityModel1D& model) override;
    
    // Configuration
    void setMaxIterations(int n) { max_iterations_ = n; }
    void setConvergenceThreshold(double thresh) { convergence_thresh_ = thresh; }
    void setDampingFactor(double factor) { damping_factor_ = factor; }
    void setFixedDepth(bool fixed, double depth = 10.0) {
        fixed_depth_ = fixed;
        default_depth_ = depth;
    }
    
    // Initial location estimate
    void setInitialLocation(const GeoPoint& loc) { initial_location_ = loc; }
    void setInitialTime(TimePoint t) { initial_time_ = t; }

private:
    VelocityModel1D velocity_model_;
    TravelTimeTable travel_times_;
    
    int max_iterations_;
    double convergence_thresh_;
    double damping_factor_;
    bool fixed_depth_;
    double default_depth_;
    
    GeoPoint initial_location_;
    TimePoint initial_time_;
    bool has_initial_location_ = false;
    bool has_initial_time_ = false;
    
    // Compute travel time and partial derivatives
    struct TravelTimeDerivatives {
        double time;      // Travel time
        double dt_dx;     // dt/dx (east-west)
        double dt_dy;     // dt/dy (north-south)
        double dt_dz;     // dt/dz (depth)
    };
    
    TravelTimeDerivatives computeDerivatives(
        const GeoPoint& source, const GeoPoint& station, PhaseType phase) const;
    
    // Build the linear system
    void buildSystem(const GeoPoint& location, TimePoint origin_time,
                     const std::vector<PickPtr>& picks,
                     const StationInventory& stations,
                     Eigen::MatrixXd& G, Eigen::VectorXd& d,
                     Eigen::VectorXd& w);
    
    // Solve the linear system with damping
    Eigen::VectorXd solveSystem(const Eigen::MatrixXd& G,
                                  const Eigen::VectorXd& d,
                                  const Eigen::VectorXd& w);
    
    // Compute RMS residual
    double computeRMS(const Eigen::VectorXd& residuals,
                       const Eigen::VectorXd& weights);
    
    // Compute location uncertainties from covariance matrix
    void computeUncertainties(const Eigen::MatrixXd& G,
                               const Eigen::VectorXd& w,
                               double rms,
                               Origin& origin);
    
    // Get initial estimate
    GeoPoint getInitialLocation(const std::vector<PickPtr>& picks,
                                 const StationInventory& stations);
    TimePoint getInitialTime(const GeoPoint& location,
                              const std::vector<PickPtr>& picks,
                              const StationInventory& stations);
};

/**
 * NonLinLocLocator - Probabilistic non-linear location
 * 
 * Uses equal differential time (EDT) formulation with
 * grid search + Oct-tree refinement.
 */
class NonLinLocLocator : public BaseLocator {
public:
    NonLinLocLocator();
    
    LocationResult locate(const std::vector<PickPtr>& picks,
                           const StationInventory& stations) override;
    
    LocationResult relocate(const EventPtr& event,
                             const StationInventory& stations) override;
    
    std::string name() const override { return "NonLinLoc"; }
    
    void setParameter(const std::string& name, double value) override;
    void setVelocityModel(const VelocityModel1D& model) override;
    
    // Probability density output
    struct PDFResult {
        GeoPoint maximum_likelihood;
        GeoPoint expectation;
        double covariance[3][3];
        std::vector<std::vector<std::vector<double>>> pdf_grid;
    };
    
    PDFResult getPDF() const { return pdf_result_; }

private:
    VelocityModel1D velocity_model_;
    TravelTimeTable travel_times_;
    
    double grid_spacing_;
    double search_radius_;
    
    PDFResult pdf_result_;
    
    // Compute likelihood at point
    double computeLikelihood(const GeoPoint& point, TimePoint origin_time,
                              const std::vector<PickPtr>& picks,
                              const StationInventory& stations);
    
    // EDT likelihood function
    double edtLikelihood(const GeoPoint& point,
                          const std::vector<PickPtr>& picks,
                          const StationInventory& stations);
};

} // namespace realdetect
