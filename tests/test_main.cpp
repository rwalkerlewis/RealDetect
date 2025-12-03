/**
 * SeisProc Test Runner
 * 
 * Main entry point for running all unit tests.
 */

#include "test_framework.hpp"
#include <iostream>
#include <string>
#include <cstring>

using namespace seisproc::test;

// All test files are compiled together, tests are auto-registered

void printUsage(const char* prog) {
    std::cout << "Usage: " << prog << " [options] [suite_name]\n"
              << "\nOptions:\n"
              << "  -h, --help     Show this help message\n"
              << "  -l, --list     List all test suites\n"
              << "  -v, --verbose  Verbose output\n"
              << "\nIf suite_name is provided, only that suite will run.\n"
              << "Otherwise, all tests will run.\n";
}

int main(int argc, char* argv[]) {
    std::cout << R"(
╔══════════════════════════════════════════════════════════════╗
║                    SeisProc Test Suite                       ║
║              Real-time Seismic Processing Tests              ║
╚══════════════════════════════════════════════════════════════╝
)" << std::endl;

    std::string suite_filter;
    bool list_only = false;
    
    // Parse arguments
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            printUsage(argv[0]);
            return 0;
        } else if (strcmp(argv[i], "-l") == 0 || strcmp(argv[i], "--list") == 0) {
            list_only = true;
        } else if (strcmp(argv[i], "-v") == 0 || strcmp(argv[i], "--verbose") == 0) {
            // Verbose mode (could be used for more output)
        } else if (argv[i][0] != '-') {
            suite_filter = argv[i];
        }
    }
    
    if (list_only) {
        std::cout << "Available test suites:\n";
        std::cout << "  GeoPoint\n";
        std::cout << "  StreamID\n";
        std::cout << "  Waveform\n";
        std::cout << "  Station\n";
        std::cout << "  StationInventory\n";
        std::cout << "  VelocityModel\n";
        std::cout << "  Config\n";
        std::cout << "  PhaseType\n";
        std::cout << "  Pick\n";
        std::cout << "  Origin\n";
        std::cout << "  Event\n";
        std::cout << "  STALTAPicker\n";
        std::cout << "  RecursiveSTALTA\n";
        std::cout << "  AICPicker\n";
        std::cout << "  CharacteristicFunction\n";
        std::cout << "  IIRFilter\n";
        std::cout << "  FFT\n";
        std::cout << "  FilterBank\n";
        std::cout << "  TravelTimeTable\n";
        std::cout << "  GridSearchLocator\n";
        std::cout << "  GeigerLocator\n";
        std::cout << "  OctTreeLocator\n";
        std::cout << "  LocatorFactory\n";
        std::cout << "  LocationQuality\n";
        std::cout << "  LocalMagnitude\n";
        std::cout << "  DurationMagnitude\n";
        std::cout << "  MomentMagnitude\n";
        std::cout << "  BodyWaveMagnitude\n";
        std::cout << "  SurfaceWaveMagnitude\n";
        std::cout << "  MagnitudeFactory\n";
        std::cout << "  MagnitudeIntegration\n";
        std::cout << "  BenchmarkPipeline\n";
        std::cout << "  BenchmarkPerformance\n";
        std::cout << "  BenchmarkAccuracy\n";
        return 0;
    }
    
    std::vector<TestResult> results;
    
    if (suite_filter.empty()) {
        // Run all tests
        results = TestRegistry::instance().runAll();
    } else {
        // Run specific suite
        results = TestRegistry::instance().runSuite(suite_filter);
    }
    
    // Print summary
    printSummary(results);
    
    // Return error code if any tests failed
    for (const auto& r : results) {
        if (!r.passed) {
            return 1;
        }
    }
    
    return 0;
}
