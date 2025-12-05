#pragma once

#include "../core/types.hpp"
#include "../core/event.hpp"
#include "../core/station.hpp"
#include "../core/velocity_model.hpp"
#include "../associator/phase_associator.hpp"
#include <vector>
#include <memory>

namespace realdetect {

/**
 * LocationResult - Result from location algorithm
 */
struct LocationResult {
    Origin origin;
    bool converged;
    int iterations;
    double initial_rms;
    double final_rms;
    std::string algorithm;
    
    LocationResult() : converged(false), iterations(0), 
                       initial_rms(0), final_rms(0) {}
};

/**
 * BaseLocator - Abstract base class for location algorithms
 */
class BaseLocator {
public:
    virtual ~BaseLocator() = default;
    
    // Locate event from picks
    virtual LocationResult locate(const std::vector<PickPtr>& picks,
                                   const StationInventory& stations) = 0;
    
    // Relocate existing event
    virtual LocationResult relocate(const EventPtr& event,
                                     const StationInventory& stations) = 0;
    
    // Get algorithm name
    virtual std::string name() const = 0;
    
    // Configuration
    virtual void setParameter(const std::string& name, double value) = 0;
    
    // Set velocity model
    virtual void setVelocityModel(const VelocityModel1D& model) = 0;
};

using LocatorPtr = std::shared_ptr<BaseLocator>;

/**
 * LocatorFactory - Create locator instances
 */
class LocatorFactory {
public:
    static LocatorPtr create(const std::string& algorithm);
    static std::vector<std::string> availableAlgorithms();
};

} // namespace realdetect
