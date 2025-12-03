#pragma once

/**
 * Regional Velocity Model Manager
 * 
 * Manages multiple 1D velocity models for different geographic regions.
 * Each region can be defined as:
 * - A polygon (list of lat/lon vertices)
 * - A circle (center point + radius)
 * - A bounding box (min/max lat/lon)
 * 
 * The manager selects the appropriate velocity model based on event location,
 * with support for priority ordering and fallback to a default global model.
 */

#include "velocity_model.hpp"
#include "types.hpp"
#include <vector>
#include <memory>
#include <map>
#include <optional>

namespace seisproc {

/**
 * GeographicBounds - Defines a geographic region
 */
class GeographicBounds {
public:
    enum class Type {
        Polygon,    // Closed polygon of lat/lon vertices
        Circle,     // Center point with radius
        BoundingBox // Simple rectangular region
    };
    
    // Create polygon bounds
    static GeographicBounds polygon(const std::vector<GeoPoint>& vertices);
    
    // Create circular bounds
    static GeographicBounds circle(const GeoPoint& center, double radius_km);
    
    // Create bounding box
    static GeographicBounds boundingBox(double min_lat, double max_lat,
                                         double min_lon, double max_lon,
                                         double min_depth = 0, double max_depth = 700);
    
    // Check if a point is inside the bounds
    bool contains(const GeoPoint& point) const;
    
    // Get bounding box (for quick rejection)
    void getBoundingBox(double& min_lat, double& max_lat,
                        double& min_lon, double& max_lon) const;
    
    // Accessors
    Type type() const { return type_; }
    const std::vector<GeoPoint>& vertices() const { return vertices_; }
    const GeoPoint& center() const { return center_; }
    double radius() const { return radius_; }
    double minDepth() const { return min_depth_; }
    double maxDepth() const { return max_depth_; }
    
    // Serialization
    std::string toString() const;
    static GeographicBounds fromString(const std::string& str);

private:
    Type type_;
    std::vector<GeoPoint> vertices_;
    GeoPoint center_;
    double radius_;         // km
    double min_lat_, max_lat_;
    double min_lon_, max_lon_;
    double min_depth_, max_depth_;
    
    // Polygon point-in-polygon test
    bool pointInPolygon(const GeoPoint& point) const;
};

/**
 * RegionalVelocityModel - A velocity model with geographic bounds
 */
struct RegionalVelocityModel {
    std::string name;           // Model identifier
    std::string description;    // Human-readable description
    GeographicBounds bounds;    // Geographic region
    VelocityModel1D model;      // The velocity model
    int priority;               // Higher priority = checked first (default 0)
    
    RegionalVelocityModel() : priority(0) {}
    
    bool contains(const GeoPoint& point) const {
        return bounds.contains(point);
    }
};

/**
 * VelocityModelManager - Manages regional velocity models
 */
class VelocityModelManager {
public:
    VelocityModelManager();
    
    // Add a regional model
    void addModel(const RegionalVelocityModel& model);
    void addModel(const std::string& name, const GeographicBounds& bounds,
                  const VelocityModel1D& vmodel, int priority = 0);
    
    // Remove a model by name
    bool removeModel(const std::string& name);
    
    // Set the default (global) model used when no regional model matches
    void setDefaultModel(const VelocityModel1D& model);
    const VelocityModel1D& defaultModel() const { return default_model_; }
    
    // Get the appropriate model for a location
    const VelocityModel1D& getModelForLocation(const GeoPoint& location) const;
    const VelocityModel1D& getModelForLocation(double lat, double lon, double depth = 0) const;
    
    // Get model by name
    const VelocityModel1D* getModelByName(const std::string& name) const;
    
    // Get the name of the model that would be used for a location
    std::string getModelNameForLocation(const GeoPoint& location) const;
    
    // List all models
    std::vector<std::string> modelNames() const;
    const std::vector<RegionalVelocityModel>& models() const { return models_; }
    
    // Load models from directory
    bool loadFromDirectory(const std::string& directory);
    
    // Load from configuration file
    bool loadFromConfig(const std::string& filename);
    
    // Save configuration
    bool saveConfig(const std::string& filename) const;
    
    // Statistics
    size_t modelCount() const { return models_.size(); }
    bool hasModels() const { return !models_.empty(); }
    
    // Diagnostic: Get all models that contain a point
    std::vector<std::string> findMatchingModels(const GeoPoint& location) const;

private:
    std::vector<RegionalVelocityModel> models_;
    VelocityModel1D default_model_;
    
    // Sort models by priority
    void sortByPriority();
};

/**
 * VelocityModelLoader - Load velocity models from various formats
 */
class VelocityModelLoader {
public:
    // Load model from file (auto-detect format)
    static bool load(const std::string& filename, VelocityModel1D& model);
    
    // Load from SeisProc format (depth thickness vp vs [density])
    static bool loadSeisProc(const std::string& filename, VelocityModel1D& model);
    
    // Load from HYPO71 format
    static bool loadHypo71(const std::string& filename, VelocityModel1D& model);
    
    // Load from VELEST format
    static bool loadVelest(const std::string& filename, VelocityModel1D& model);
    
    // Load from NonLinLoc format
    static bool loadNonLinLoc(const std::string& filename, VelocityModel1D& model);
    
    // Save to SeisProc format
    static bool saveSeisProc(const std::string& filename, const VelocityModel1D& model);
};

/**
 * RegionalModelConfig - Configuration for regional velocity models
 */
struct RegionalModelConfig {
    std::string name;
    std::string model_file;
    std::string bounds_type;   // "polygon", "circle", "box"
    std::string bounds_spec;   // Bounds specification string
    int priority;
    std::string description;
    
    // Parse bounds specification
    GeographicBounds parseBounds() const;
};

} // namespace seisproc
