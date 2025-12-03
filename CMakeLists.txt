cmake_minimum_required(VERSION 3.16)
project(SeisProc VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Build options
option(BUILD_TESTS "Build unit tests" ON)
option(USE_OPENMP "Enable OpenMP parallelization" ON)

# Find required packages
find_package(Threads REQUIRED)
find_package(Eigen3 3.3 QUIET)

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
    src/core/config.cpp
)

set(SEEDLINK_SOURCES
    src/seedlink/seedlink_client.cpp
    src/seedlink/packet_parser.cpp
    src/seedlink/stream_buffer.cpp
)

set(PICKER_SOURCES
    src/picker/stalta_picker.cpp
    src/picker/aic_picker.cpp
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
add_library(seisproc_core STATIC
    ${CORE_SOURCES}
    ${SEEDLINK_SOURCES}
    ${PICKER_SOURCES}
    ${ASSOCIATOR_SOURCES}
    ${LOCATOR_SOURCES}
    ${MAGNITUDE_SOURCES}
)

target_link_libraries(seisproc_core
    Threads::Threads
    ${CMAKE_DL_LIBS}
)

if(OpenMP_CXX_FOUND)
    target_link_libraries(seisproc_core OpenMP::OpenMP_CXX)
endif()

# Main executable
add_executable(seisproc src/main.cpp)
target_link_libraries(seisproc seisproc_core)

# Simulator for testing
add_executable(seissim src/simulator.cpp)
target_link_libraries(seissim seisproc_core)

# Install targets
install(TARGETS seisproc seissim
    RUNTIME DESTINATION bin
)

install(DIRECTORY include/
    DESTINATION include
)

# Copy config files
install(DIRECTORY config/
    DESTINATION share/seisproc/config
)
