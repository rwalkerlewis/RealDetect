# SeisProc - Real-time Seismic Event Processing System

A comprehensive C++ seismic processing system inspired by [SeisComP](https://www.seiscomp.de/) for real-time earthquake detection, location, and magnitude calculation.

## Features

### Data Acquisition
- **SeedLink Client**: Real-time data streaming from SeedLink servers
- **MiniSEED Parser**: Native support for SEED data format with Steim1/Steim2 decompression
- **Multi-stream Buffering**: Efficient circular buffer management for multiple channels

### Phase Picking
- **STA/LTA Detector**: Classic Short-Term/Long-Term Average ratio trigger
- **AIC Picker**: Akaike Information Criterion for precise pick refinement
- **AR Picker**: Autoregressive prediction error method
- **Characteristic Functions**: 
  - Envelope (Hilbert transform)
  - Kurtosis (4th moment)
  - Energy ratio
  - Polarization analysis
- **Multi-frequency Analysis**: Parallel detection in multiple frequency bands
- **IIR Filtering**: Butterworth bandpass/highpass/lowpass filters

### Event Association
- **Phase Associator**: Travel-time based pick grouping
- **Nucleator**: Grid-based source nucleation for robust association
- **Travel-time Tables**: Pre-computed tables for efficient lookup

### Hypocenter Location
- **Grid Search**: Brute-force search with multi-level refinement
- **OctTree Search**: Adaptive octree subdivision for efficiency
- **Geiger's Method**: Iterative linearized inversion with damping
- **NonLinLoc-style**: Equal Differential Time (EDT) probability location
- **1D Velocity Models**: Support for layered velocity models (IASP91, AK135)

### Magnitude Calculation
- **ML (Local Magnitude)**: Wood-Anderson simulation, regional calibration
- **Mw (Moment Magnitude)**: Spectral analysis with Brune model fitting
- **Mb (Body Wave Magnitude)**: Teleseismic P-wave amplitude
- **Ms (Surface Wave Magnitude)**: 20-second Rayleigh wave
- **Md (Duration Magnitude)**: Coda duration method

## Building

### Requirements
- C++17 compatible compiler (GCC 7+, Clang 5+)
- CMake 3.16+
- Eigen3 (for linear algebra in Geiger inversion)
- pthreads

### Optional
- OpenMP (for parallel processing)

### Build Steps

```bash
# Create build directory
mkdir build && cd build

# Configure
cmake .. -DCMAKE_BUILD_TYPE=Release

# Build
make -j$(nproc)

# Install (optional)
sudo make install
```

## Usage

### Real-time Processing

```bash
# Connect to IRIS SeedLink and process data
./seisproc -s rtserve.iris.washington.edu -p 18000

# With custom configuration
./seisproc -c config/seisproc.conf -S config/stations.txt -v config/velocity_model.txt

# Connect to local server
./seisproc -s localhost -p 18000
```

### Simulator (for testing)

```bash
# Run synthetic event simulation
./seissim

# The simulator creates a test network and generates
# synthetic seismograms with realistic P and S arrivals
```

### Command Line Options

```
Usage: seisproc [options]

Options:
  -c, --config <file>    Configuration file (default: seisproc.conf)
  -s, --server <host>    SeedLink server host (default: localhost)
  -p, --port <port>      SeedLink server port (default: 18000)
  -S, --stations <file>  Station inventory file
  -v, --velocity <file>  Velocity model file
  -h, --help             Show this help message
```

## Configuration

See `config/seisproc.conf` for detailed configuration options:

```ini
[picker]
sta_length = 0.5          # STA window (seconds)
lta_length = 10.0         # LTA window (seconds)
trigger_ratio = 3.5       # STA/LTA trigger threshold

[locator]
algorithm = geiger        # Location algorithm (geiger, gridsearch, octtree)
fixed_depth = false       # Fix depth to default value
default_depth = 10.0      # Default/fixed depth (km)

[magnitude]
ml_enabled = true         # Calculate local magnitude
mw_enabled = true         # Calculate moment magnitude
```

## File Formats

### Station Inventory
```
# Network Station Latitude Longitude Elevation(m)
CI ADO    34.5505 -117.4339 1500
IU ANMO   34.9462 -106.4567 1850
```

### Velocity Model
```
# Depth(km) Thickness(km) Vp(km/s) Vs(km/s) Density(g/cm³)
0.0    5.0     5.50    3.18    2.60
5.0    10.0    6.00    3.46    2.75
15.0   17.0    6.50    3.75    2.90
32.0   0.0     7.80    4.50    3.30
```

## Architecture

```
seisproc/
├── include/seisproc/
│   ├── core/           # Core data structures
│   │   ├── types.hpp       # Basic types, GeoPoint, StreamID
│   │   ├── waveform.hpp    # Waveform container
│   │   ├── station.hpp     # Station/channel metadata
│   │   ├── event.hpp       # Event/Origin/Pick structures
│   │   ├── velocity_model.hpp
│   │   └── config.hpp
│   │
│   ├── seedlink/       # SeedLink data acquisition
│   │   └── seedlink_client.hpp
│   │
│   ├── picker/         # Phase detection
│   │   ├── picker.hpp
│   │   ├── stalta_picker.hpp
│   │   ├── aic_picker.hpp
│   │   ├── characteristic_function.hpp
│   │   └── filter_bank.hpp
│   │
│   ├── associator/     # Event association
│   │   └── phase_associator.hpp
│   │
│   ├── locator/        # Hypocenter location
│   │   ├── locator.hpp
│   │   ├── grid_search.hpp
│   │   └── geiger.hpp
│   │
│   └── magnitude/      # Magnitude calculation
│       ├── magnitude.hpp
│       ├── local_magnitude.hpp
│       └── moment_magnitude.hpp
│
└── src/                # Implementations
```

## Processing Pipeline

```
Data Flow:
                                                    
  SeedLink    →   Buffer   →   Picker   →   Associator   →   Locator   →   Magnitude
     │              │           │               │               │              │
  Real-time    Circular      STA/LTA        Travel-time      Geiger          ML
  MiniSEED     buffer        + AIC          clustering      inversion        Mw
```

## Algorithms

### STA/LTA Picker
The classic trigger algorithm computes the ratio of short-term to long-term energy:
```
R(t) = STA(t) / LTA(t)
Trigger when R > threshold
```

### Geiger Location
Iterative linearized inversion:
```
G · Δm = d
where G = ∂t/∂(x,y,z,t) Jacobian matrix
      d = observed - calculated residuals
```

### Local Magnitude
```
ML = log10(A) + a·log10(D/100) + b·(D-100) + c
where A = Wood-Anderson amplitude (mm)
      D = epicentral distance (km)
```

### Moment Magnitude  
```
Mw = (2/3)·log10(M0) - 10.7
where M0 = seismic moment from spectral analysis
```

## Performance

- Handles 100+ real-time streams
- <100ms latency from data arrival to pick
- ~1s for complete event processing
- Memory: ~50MB base + ~1MB per stream buffer

## References

- Withers et al. (1998). A comparison of select trigger algorithms for automated global seismic phase and event detection
- Lomax et al. (2000). Probabilistic earthquake location in 3D and layered models
- Brune (1970). Tectonic stress and the spectra of seismic shear waves
- Havskov & Ottemoller (2010). Routine Data Processing in Earthquake Seismology

## License

MIT License - see LICENSE file

## Contributing

Contributions welcome! Please submit issues and pull requests on GitHub.
