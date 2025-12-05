cmake_minimum_required(VERSION 3.16)
project(RealDetect VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Build options
option(BUILD_TESTS "Build unit tests" ON)
option(USE_OPENMP "Enable OpenMP parallelization" ON)

# Find required packages
find_package(Threads REQUIRED)
find_package(Eigen3 3.3 QUIET)
find_package(SQLite3 QUIET)

if(NOT SQLite3_FOUND)
    # Try pkg-config as fallback
    find_package(PkgConfig)
    if(PKG_CONFIG_FOUND)
        pkg_check_modules(SQLITE3 sqlite3)
    endif()
    
    if(NOT SQLITE3_FOUND)
        message(STATUS "SQLite3 not found, will try to use system library directly")
        set(SQLITE3_LIBRARIES sqlite3)
    endif()
else()
    message(STATUS "Found SQLite3: ${SQLite3_LIBRARIES}")
endif()

if(NOT EIGEN3_FOUND)
    # If Eigen not found, use bundled header-only version
    message(STATUS "Eigen3 not found, using bundled headers")
    set(EIGEN3_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/include/eigen3)
else()
    message(STATUS "Found Eigen3: ${EIGEN3_INCLUDE_DIR}")
endif()

include_directories(${EIGEN3_INCLUDE_DIR})

if(USE_OPENMP)
    find_package(OpenMP)
    if(OpenMP_CXX_FOUND)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    endif()
endif()

# Compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -O3")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -DDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

# Include directories
include_directories(${CMAKE_SOURCE_DIR}/include)

# Source files
set(CORE_SOURCES
    src/core/miniseed.cpp
    src/core/station.cpp
    src/core/event.cpp
    src/core/waveform.cpp
    src/core/velocity_model.cpp
    src/core/regional_velocity_model.cpp
    src/core/config.cpp
)

set(DATABASE_SOURCES
    src/database/css30_database.cpp
)

set(SEEDLINK_SOURCES
    src/seedlink/seedlink_client.cpp
    src/seedlink/packet_parser.cpp
    src/seedlink/stream_buffer.cpp
)

set(PICKER_SOURCES
    src/picker/stalta_picker.cpp
    src/picker/aic_picker.cpp
    src/picker/ml_picker.cpp
    src/picker/characteristic_function.cpp
    src/picker/filter_bank.cpp
)

set(ASSOCIATOR_SOURCES
    src/associator/phase_associator.cpp
    src/associator/event_builder.cpp
    src/associator/travel_time.cpp
)

set(LOCATOR_SOURCES
    src/locator/grid_search.cpp
    src/locator/geiger.cpp
    src/locator/velocity_model_1d.cpp
    src/locator/ray_tracer.cpp
)

set(MAGNITUDE_SOURCES
    src/magnitude/local_magnitude.cpp
    src/magnitude/moment_magnitude.cpp
    src/magnitude/spectral_analysis.cpp
)

# Core library
add_library(realdetect_core STATIC
    ${CORE_SOURCES}
    ${DATABASE_SOURCES}
    ${SEEDLINK_SOURCES}
    ${PICKER_SOURCES}
    ${ASSOCIATOR_SOURCES}
    ${LOCATOR_SOURCES}
    ${MAGNITUDE_SOURCES}
)

# Link SQLite3
if(SQLite3_FOUND)
    target_link_libraries(realdetect_core ${SQLite3_LIBRARIES})
    target_include_directories(realdetect_core PRIVATE ${SQLite3_INCLUDE_DIRS})
elseif(SQLITE3_FOUND)
    target_link_libraries(realdetect_core ${SQLITE3_LIBRARIES})
    target_include_directories(realdetect_core PRIVATE ${SQLITE3_INCLUDE_DIRS})
else()
    target_link_libraries(realdetect_core sqlite3)
endif()

target_link_libraries(realdetect_core
    Threads::Threads
    ${CMAKE_DL_LIBS}
)

if(OpenMP_CXX_FOUND)
    target_link_libraries(realdetect_core OpenMP::OpenMP_CXX)
endif()

# Main executable
add_executable(realdetect src/main.cpp)
target_link_libraries(realdetect realdetect_core)

# Simulator for testing
add_executable(realdetect_sim src/simulator.cpp)
target_link_libraries(realdetect_sim realdetect_core)

# Install targets
install(TARGETS realdetect realdetect_sim
    RUNTIME DESTINATION bin
)

install(DIRECTORY include/
    DESTINATION include
)

# Copy config files
install(DIRECTORY config/
    DESTINATION share/realdetect/config
)

# Unit tests
if(BUILD_TESTS)
    enable_testing()
    
    set(TEST_SOURCES
        tests/test_main.cpp
        tests/test_core.cpp
        tests/test_picker.cpp
        tests/test_locator.cpp
        tests/test_magnitude.cpp
        tests/test_benchmark.cpp
    )
    
    add_executable(realdetect_tests ${TEST_SOURCES})
    target_include_directories(realdetect_tests PRIVATE ${CMAKE_SOURCE_DIR}/tests)
    target_link_libraries(realdetect_tests realdetect_core)
    
    # Add test target
    add_test(NAME unit_tests COMMAND realdetect_tests)
    
    # Add convenience targets for running tests
    add_custom_target(run_tests
        COMMAND realdetect_tests
        DEPENDS realdetect_tests
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        COMMENT "Running unit tests..."
    )
    
    # Run specific test suites
    add_custom_target(run_benchmark
        COMMAND realdetect_tests BenchmarkPipeline BenchmarkPerformance BenchmarkAccuracy
        DEPENDS realdetect_tests
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        COMMENT "Running benchmark tests..."
    )
endif()
