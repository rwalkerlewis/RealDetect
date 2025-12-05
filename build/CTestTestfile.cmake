# CMake generated Testfile for 
# Source directory: /workspace
# Build directory: /workspace/build
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(unit_tests "/workspace/build/realdetect_tests")
set_tests_properties(unit_tests PROPERTIES  _BACKTRACE_TRIPLES "/workspace/CMakeLists.txt;181;add_test;/workspace/CMakeLists.txt;0;")
add_test(core_tests "/workspace/build/realdetect_tests" "GeoPoint" "StreamID" "Waveform" "Station" "Config")
set_tests_properties(core_tests PROPERTIES  _BACKTRACE_TRIPLES "/workspace/CMakeLists.txt;182;add_test;/workspace/CMakeLists.txt;0;")
add_test(picker_tests "/workspace/build/realdetect_tests" "STALTAPicker" "AICPicker" "CharacteristicFunction")
set_tests_properties(picker_tests PROPERTIES  _BACKTRACE_TRIPLES "/workspace/CMakeLists.txt;183;add_test;/workspace/CMakeLists.txt;0;")
add_test(locator_tests "/workspace/build/realdetect_tests" "TravelTimeTable" "GridSearchLocator" "GeigerLocator")
set_tests_properties(locator_tests PROPERTIES  _BACKTRACE_TRIPLES "/workspace/CMakeLists.txt;184;add_test;/workspace/CMakeLists.txt;0;")
add_test(magnitude_tests "/workspace/build/realdetect_tests" "LocalMagnitude" "MomentMagnitude")
set_tests_properties(magnitude_tests PROPERTIES  _BACKTRACE_TRIPLES "/workspace/CMakeLists.txt;185;add_test;/workspace/CMakeLists.txt;0;")
add_test(database_tests "/workspace/build/realdetect_tests" "CSS30Database" "CSS30Writer")
set_tests_properties(database_tests PROPERTIES  _BACKTRACE_TRIPLES "/workspace/CMakeLists.txt;186;add_test;/workspace/CMakeLists.txt;0;")
add_test(associator_tests "/workspace/build/realdetect_tests" "PhaseAssociator" "NucleatorAssociator")
set_tests_properties(associator_tests PROPERTIES  _BACKTRACE_TRIPLES "/workspace/CMakeLists.txt;187;add_test;/workspace/CMakeLists.txt;0;")
add_test(regional_velocity_tests "/workspace/build/realdetect_tests" "GeographicBounds" "RegionalVelocityModel" "VelocityModelManager")
set_tests_properties(regional_velocity_tests PROPERTIES  _BACKTRACE_TRIPLES "/workspace/CMakeLists.txt;188;add_test;/workspace/CMakeLists.txt;0;")
add_test(integration_tests "/workspace/build/realdetect_tests" "Integration")
set_tests_properties(integration_tests PROPERTIES  _BACKTRACE_TRIPLES "/workspace/CMakeLists.txt;189;add_test;/workspace/CMakeLists.txt;0;")
add_test(benchmark_tests "/workspace/build/realdetect_tests" "BenchmarkPipeline" "BenchmarkPerformance" "BenchmarkAccuracy")
set_tests_properties(benchmark_tests PROPERTIES  _BACKTRACE_TRIPLES "/workspace/CMakeLists.txt;190;add_test;/workspace/CMakeLists.txt;0;")
