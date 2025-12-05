/**
 * Unit tests for regional velocity model management
 */

#include "test_framework.hpp"
#include "realdetect/core/regional_velocity_model.hpp"
#include "realdetect/core/velocity_model.hpp"
#include <cmath>

using namespace realdetect;
using namespace realdetect::test;

// ============================================================================
// GeographicBounds Tests
// ============================================================================

TEST(GeographicBounds, BoundingBox) {
    auto bounds = GeographicBounds::boundingBox(33.0, 35.0, -119.0, -117.0);
    
    ASSERT_TRUE(bounds.type() == GeographicBounds::Type::BoundingBox);
    
    // Points inside
    ASSERT_TRUE(bounds.contains(GeoPoint(34.0, -118.0, 0.0)));
    ASSERT_TRUE(bounds.contains(GeoPoint(33.5, -118.5, 10.0)));
    ASSERT_TRUE(bounds.contains(GeoPoint(34.9, -117.1, 5.0)));
    
    // Points outside
    ASSERT_FALSE(bounds.contains(GeoPoint(32.0, -118.0, 0.0)));  // Too far south
    ASSERT_FALSE(bounds.contains(GeoPoint(36.0, -118.0, 0.0)));  // Too far north
    ASSERT_FALSE(bounds.contains(GeoPoint(34.0, -120.0, 0.0)));  // Too far west
    ASSERT_FALSE(bounds.contains(GeoPoint(34.0, -116.0, 0.0)));  // Too far east
}

TEST(GeographicBounds, BoundingBoxWithDepth) {
    auto bounds = GeographicBounds::boundingBox(33.0, 35.0, -119.0, -117.0, 0.0, 30.0);
    
    // Point within depth range
    ASSERT_TRUE(bounds.contains(GeoPoint(34.0, -118.0, 15.0)));
    
    // Point below depth range
    ASSERT_FALSE(bounds.contains(GeoPoint(34.0, -118.0, 50.0)));
}

TEST(GeographicBounds, Circle) {
    GeoPoint center(34.0, -118.0);
    double radius_km = 100.0;
    
    auto bounds = GeographicBounds::circle(center, radius_km);
    
    ASSERT_TRUE(bounds.type() == GeographicBounds::Type::Circle);
    ASSERT_NEAR(bounds.center().latitude, 34.0, 1e-10);
    ASSERT_NEAR(bounds.center().longitude, -118.0, 1e-10);
    ASSERT_NEAR(bounds.radius(), 100.0, 1e-10);
    
    // Point at center
    ASSERT_TRUE(bounds.contains(GeoPoint(34.0, -118.0, 0.0)));
    
    // Point ~50 km away (should be inside)
    // 0.5 degrees latitude ~ 55 km
    ASSERT_TRUE(bounds.contains(GeoPoint(34.4, -118.0, 0.0)));
    
    // Point ~200 km away (should be outside)
    ASSERT_FALSE(bounds.contains(GeoPoint(36.0, -118.0, 0.0)));
}

TEST(GeographicBounds, Polygon) {
    // Define a rectangular polygon
    std::vector<GeoPoint> vertices = {
        GeoPoint(33.0, -119.0),
        GeoPoint(33.0, -117.0),
        GeoPoint(35.0, -117.0),
        GeoPoint(35.0, -119.0)
    };
    
    auto bounds = GeographicBounds::polygon(vertices);
    
    ASSERT_TRUE(bounds.type() == GeographicBounds::Type::Polygon);
    ASSERT_EQ(bounds.vertices().size(), 4u);
    
    // Points inside
    ASSERT_TRUE(bounds.contains(GeoPoint(34.0, -118.0, 0.0)));
    
    // Points outside
    ASSERT_FALSE(bounds.contains(GeoPoint(32.0, -118.0, 0.0)));
    ASSERT_FALSE(bounds.contains(GeoPoint(36.0, -118.0, 0.0)));
}

TEST(GeographicBounds, GetBoundingBox) {
    auto bounds = GeographicBounds::circle(GeoPoint(34.0, -118.0), 100.0);
    
    double min_lat, max_lat, min_lon, max_lon;
    bounds.getBoundingBox(min_lat, max_lat, min_lon, max_lon);
    
    // Bounding box should enclose the circle
    ASSERT_LT(min_lat, 34.0);
    ASSERT_GT(max_lat, 34.0);
    ASSERT_LT(min_lon, -118.0);
    ASSERT_GT(max_lon, -118.0);
}

TEST(GeographicBounds, ToString) {
    auto bounds = GeographicBounds::boundingBox(33.0, 35.0, -119.0, -117.0);
    
    std::string str = bounds.toString();
    ASSERT_FALSE(str.empty());
}

TEST(GeographicBounds, FromString) {
    auto original = GeographicBounds::boundingBox(33.0, 35.0, -119.0, -117.0);
    std::string str = original.toString();
    
    auto parsed = GeographicBounds::fromString(str);
    
    // Should contain same point
    ASSERT_TRUE(parsed.contains(GeoPoint(34.0, -118.0, 0.0)));
}

// ============================================================================
// RegionalVelocityModel Tests
// ============================================================================

TEST(RegionalVelocityModel, DefaultConstructor) {
    RegionalVelocityModel model;
    
    ASSERT_TRUE(model.name.empty());
    ASSERT_EQ(model.priority, 0);
}

TEST(RegionalVelocityModel, Contains) {
    RegionalVelocityModel model;
    model.name = "test_model";
    model.bounds = GeographicBounds::boundingBox(33.0, 35.0, -119.0, -117.0);
    model.model = VelocityModel1D::simpleThreeLayer();
    model.priority = 10;
    
    ASSERT_TRUE(model.contains(GeoPoint(34.0, -118.0, 10.0)));
    ASSERT_FALSE(model.contains(GeoPoint(40.0, -118.0, 10.0)));
}

// ============================================================================
// VelocityModelManager Tests
// ============================================================================

TEST(VelocityModelManager, DefaultConstructor) {
    VelocityModelManager manager;
    
    ASSERT_EQ(manager.modelCount(), 0u);
    ASSERT_FALSE(manager.hasModels());
}

TEST(VelocityModelManager, SetDefaultModel) {
    VelocityModelManager manager;
    
    auto model = VelocityModel1D::simpleThreeLayer();
    manager.setDefaultModel(model);
    
    ASSERT_EQ(manager.defaultModel().layerCount(), model.layerCount());
}

TEST(VelocityModelManager, AddModel) {
    VelocityModelManager manager;
    
    RegionalVelocityModel rmodel;
    rmodel.name = "socal";
    rmodel.bounds = GeographicBounds::boundingBox(32.0, 36.0, -121.0, -114.0);
    rmodel.model = VelocityModel1D::simpleThreeLayer();
    rmodel.priority = 10;
    
    manager.addModel(rmodel);
    
    ASSERT_EQ(manager.modelCount(), 1u);
    ASSERT_TRUE(manager.hasModels());
}

TEST(VelocityModelManager, AddMultipleModels) {
    VelocityModelManager manager;
    
    // Add Southern California model
    manager.addModel("socal", 
        GeographicBounds::boundingBox(32.0, 36.0, -121.0, -114.0),
        VelocityModel1D::simpleThreeLayer(), 10);
    
    // Add Northern California model
    manager.addModel("nocal",
        GeographicBounds::boundingBox(36.0, 42.0, -125.0, -119.0),
        VelocityModel1D::iasp91(), 10);
    
    // Add Hawaii model
    manager.addModel("hawaii",
        GeographicBounds::circle(GeoPoint(19.5, -155.5), 200.0),
        VelocityModel1D::simpleThreeLayer(), 20);
    
    ASSERT_EQ(manager.modelCount(), 3u);
}

TEST(VelocityModelManager, RemoveModel) {
    VelocityModelManager manager;
    
    manager.addModel("test",
        GeographicBounds::boundingBox(33.0, 35.0, -119.0, -117.0),
        VelocityModel1D::simpleThreeLayer(), 0);
    
    ASSERT_EQ(manager.modelCount(), 1u);
    
    bool removed = manager.removeModel("test");
    ASSERT_TRUE(removed);
    ASSERT_EQ(manager.modelCount(), 0u);
    
    // Try to remove non-existent
    removed = manager.removeModel("nonexistent");
    ASSERT_FALSE(removed);
}

TEST(VelocityModelManager, GetModelForLocation) {
    VelocityModelManager manager;
    
    // Set default model
    manager.setDefaultModel(VelocityModel1D::iasp91());
    
    // Add regional model for SoCal
    VelocityModel1D socal_model("SoCal");
    socal_model.addLayer(0, 5, 5.5, 3.2);
    socal_model.addLayer(5, 10, 6.0, 3.5);
    socal_model.addLayer(15, 20, 6.5, 3.8);
    socal_model.addLayer(35, 0, 8.0, 4.6);
    
    manager.addModel("socal",
        GeographicBounds::boundingBox(32.0, 36.0, -121.0, -114.0),
        socal_model, 10);
    
    // Point in SoCal should get SoCal model
    const auto& model_socal = manager.getModelForLocation(GeoPoint(34.0, -118.0, 10.0));
    ASSERT_EQ(model_socal.name(), "SoCal");
    
    // Point outside should get default
    const auto& model_default = manager.getModelForLocation(GeoPoint(50.0, -100.0, 10.0));
    ASSERT_NE(model_default.name(), "SoCal");
}

TEST(VelocityModelManager, GetModelForLocationCoordinates) {
    VelocityModelManager manager;
    manager.setDefaultModel(VelocityModel1D::iasp91());
    
    manager.addModel("test",
        GeographicBounds::boundingBox(33.0, 35.0, -119.0, -117.0),
        VelocityModel1D::simpleThreeLayer(), 0);
    
    // Test using lat/lon/depth instead of GeoPoint
    const auto& model = manager.getModelForLocation(34.0, -118.0, 10.0);
    ASSERT_TRUE(model.layerCount() > 0);
}

TEST(VelocityModelManager, GetModelByName) {
    VelocityModelManager manager;
    
    manager.addModel("socal",
        GeographicBounds::boundingBox(32.0, 36.0, -121.0, -114.0),
        VelocityModel1D::simpleThreeLayer(), 10);
    
    const auto* model = manager.getModelByName("socal");
    ASSERT_TRUE(model != nullptr);
    
    const auto* null_model = manager.getModelByName("nonexistent");
    ASSERT_TRUE(null_model == nullptr);
}

TEST(VelocityModelManager, GetModelNameForLocation) {
    VelocityModelManager manager;
    manager.setDefaultModel(VelocityModel1D::iasp91());
    
    manager.addModel("socal",
        GeographicBounds::boundingBox(32.0, 36.0, -121.0, -114.0),
        VelocityModel1D::simpleThreeLayer(), 10);
    
    std::string name = manager.getModelNameForLocation(GeoPoint(34.0, -118.0, 10.0));
    ASSERT_EQ(name, "socal");
    
    std::string default_name = manager.getModelNameForLocation(GeoPoint(50.0, -100.0, 10.0));
    ASSERT_NE(default_name, "socal");
}

TEST(VelocityModelManager, ModelNames) {
    VelocityModelManager manager;
    
    manager.addModel("model_a",
        GeographicBounds::boundingBox(33.0, 35.0, -119.0, -117.0),
        VelocityModel1D::simpleThreeLayer(), 0);
    
    manager.addModel("model_b",
        GeographicBounds::circle(GeoPoint(40.0, -120.0), 100.0),
        VelocityModel1D::iasp91(), 0);
    
    auto names = manager.modelNames();
    
    // Should have at least the 2 models we added
    ASSERT_GE(names.size(), 2u);
    
    bool has_a = false, has_b = false;
    for (const auto& name : names) {
        if (name == "model_a") has_a = true;
        if (name == "model_b") has_b = true;
    }
    ASSERT_TRUE(has_a);
    ASSERT_TRUE(has_b);
}

TEST(VelocityModelManager, FindMatchingModels) {
    VelocityModelManager manager;
    
    // Overlapping regions
    manager.addModel("region1",
        GeographicBounds::boundingBox(33.0, 36.0, -120.0, -116.0),
        VelocityModel1D::simpleThreeLayer(), 0);
    
    manager.addModel("region2",
        GeographicBounds::boundingBox(34.0, 37.0, -119.0, -115.0),
        VelocityModel1D::iasp91(), 0);
    
    // Point in overlap
    auto matches = manager.findMatchingModels(GeoPoint(35.0, -117.0, 10.0));
    ASSERT_EQ(matches.size(), 2u);
    
    // Point in region1 only
    auto matches1 = manager.findMatchingModels(GeoPoint(33.5, -119.0, 10.0));
    ASSERT_GE(matches1.size(), 1u);
}

TEST(VelocityModelManager, Priority) {
    VelocityModelManager manager;
    manager.setDefaultModel(VelocityModel1D::iasp91());
    
    // Two overlapping regions with different priorities
    VelocityModel1D low_priority_model("LowPriority");
    low_priority_model.addLayer(0, 10, 5.0, 2.9);
    low_priority_model.addLayer(10, 0, 8.0, 4.5);
    
    VelocityModel1D high_priority_model("HighPriority");
    high_priority_model.addLayer(0, 10, 6.0, 3.5);
    high_priority_model.addLayer(10, 0, 8.5, 4.9);
    
    manager.addModel("low",
        GeographicBounds::boundingBox(33.0, 36.0, -120.0, -116.0),
        low_priority_model, 5);
    
    manager.addModel("high",
        GeographicBounds::boundingBox(34.0, 35.0, -119.0, -117.0),
        high_priority_model, 15);  // Higher priority
    
    // Point in both regions should get high priority model
    const auto& model = manager.getModelForLocation(GeoPoint(34.5, -118.0, 10.0));
    ASSERT_EQ(model.name(), "HighPriority");
}

// ============================================================================
// VelocityModelLoader Tests
// ============================================================================

TEST(VelocityModelLoader, LoadNonExistent) {
    VelocityModel1D model;
    
    bool result = VelocityModelLoader::load("/tmp/nonexistent_file.vel", model);
    ASSERT_FALSE(result);
}

// ============================================================================
// Integration Tests
// ============================================================================

TEST(VelocityModelManager, IntegrationWithLocator) {
    VelocityModelManager manager;
    
    // Set up default model
    manager.setDefaultModel(VelocityModel1D::iasp91());
    
    // Add regional models
    manager.addModel("socal",
        GeographicBounds::boundingBox(32.0, 36.0, -121.0, -114.0),
        VelocityModel1D::simpleThreeLayer(), 10);
    
    manager.addModel("hawaii",
        GeographicBounds::circle(GeoPoint(19.5, -155.5), 200.0),
        VelocityModel1D::simpleThreeLayer(), 20);
    
    // Simulate getting model for different event locations
    std::vector<GeoPoint> test_points = {
        GeoPoint(34.0, -118.0, 10.0),   // SoCal
        GeoPoint(19.8, -155.0, 5.0),    // Hawaii
        GeoPoint(40.0, -100.0, 15.0),   // Outside - should get default
    };
    
    std::vector<std::string> expected_models = {"socal", "hawaii", ""};
    
    for (size_t i = 0; i < test_points.size(); i++) {
        std::string name = manager.getModelNameForLocation(test_points[i]);
        if (i < 2) {
            ASSERT_EQ(name, expected_models[i]);
        } else {
            ASSERT_NE(name, "socal");
            ASSERT_NE(name, "hawaii");
        }
    }
}

TEST(VelocityModelManager, TravelTimeCalculation) {
    VelocityModelManager manager;
    
    // Different models should give different travel times
    VelocityModel1D fast_model("Fast");
    fast_model.addLayer(0, 35, 7.0, 4.0);  // Fast crust
    fast_model.addLayer(35, 0, 8.5, 4.9);
    
    VelocityModel1D slow_model("Slow");
    slow_model.addLayer(0, 35, 5.0, 2.9);  // Slow crust
    slow_model.addLayer(35, 0, 7.5, 4.3);
    
    manager.setDefaultModel(slow_model);
    manager.addModel("fast_region",
        GeographicBounds::boundingBox(34.0, 36.0, -119.0, -117.0),
        fast_model, 10);
    
    // Get models and calculate travel times
    GeoPoint in_fast(35.0, -118.0, 10.0);
    GeoPoint in_slow(40.0, -100.0, 10.0);
    
    const auto& model_fast = manager.getModelForLocation(in_fast);
    const auto& model_slow = manager.getModelForLocation(in_slow);
    
    double tt_fast = model_fast.travelTime(100.0, 10.0, PhaseType::P);
    double tt_slow = model_slow.travelTime(100.0, 10.0, PhaseType::P);
    
    // Fast model should give shorter travel time
    ASSERT_LT(tt_fast, tt_slow);
}
