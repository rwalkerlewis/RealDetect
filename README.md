# RealDetect - Real-time Seismic Event Processing System

A comprehensive C++ seismic processing system inspired by [SeisComP](https://www.seiscomp.de/) for real-time earthquake detection, location, and magnitude calculation.

## Features

### Data Acquisition
- **SeedLink Client**: Real-time data streaming from SeedLink servers
- **MiniSEED Parser**: Native support for SEED data format with Steim1/Steim2 decompression
- **Multi-stream Buffering**: Efficient circular buffer management for multiple channels
- **Playback Mode**: Process MiniSEED files for offline analysis or testing

### Phase Picking
- **STA/LTA Detector**: Classic Short-Term/Long-Term Average ratio trigger
- **AIC Picker**: Akaike Information Criterion for precise pick refinement
- **ML Picker**: Machine learning based phase detection with pluggable backends
  - Built-in characteristic function detector
  - ONNX Runtime support for neural network models (PhaseNet, EQTransformer, etc.)
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

### Database Output (CSS3.0)
- **SQLite Database**: Event catalog storage in CSS3.0 standard schema
- **Complete Event Storage**: Events, origins, arrivals, associations, magnitudes
- **Station Metadata**: Site, sitechan, affiliation tables
- **Query Support**: Time, location, magnitude filtering
- **Flat File Export**: Export to CSS3.0 flat files for interoperability

### Regional Velocity Models
- **Geographic Regions**: Define 1D models for specific geographic areas
- **Region Types**: Bounding box, circle, or polygon boundaries
- **Priority System**: Automatic selection based on event location
- **Multiple Formats**: Support for HYPO71, VELEST, NonLinLoc formats
- **Fallback**: Global default model when no regional model matches

## Building

### Requirements
- C++17 compatible compiler (GCC 7+, Clang 5+)
- CMake 3.16+
- Eigen3 (for linear algebra in Geiger inversion)
- SQLite3 (for CSS3.0 database output)
- pthreads

### Optional
- OpenMP (for parallel processing)
- ONNX Runtime (for ML picker with neural network models)

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
./realdetect -s rtserve.iris.washington.edu -p 18000

# With custom configuration
./realdetect -c config/realdetect.conf -S config/stations.txt -v config/velocity_model.txt

# Connect to local server
./realdetect -s localhost -p 18000

# With CSS3.0 database output
./realdetect -s localhost -p 18000 -d catalog.db

# With regional velocity models
./realdetect -s localhost -p 18000 -V config/velocity_models -d catalog.db
```

### Playback Mode

Process archived MiniSEED data for offline analysis:

```bash
# Playback a single MiniSEED file
./realdetect -m data.mseed -d catalog.db

# Playback all files in a directory
./realdetect -M /path/to/miniseed/files -d catalog.db

# Fast playback (as fast as possible)
./realdetect -m data.mseed --playback-speed 0

# Slower playback for debugging (0.5x real-time)
./realdetect -m data.mseed --playback-speed 0.5
```

### ML Picker

Use machine learning for phase detection:

```bash
# Use built-in ML detector
./realdetect --picker ml -m data.mseed

# Use custom ML model (ONNX format)
./realdetect --picker ml --ml-model phasenet.onnx -m data.mseed

# Use ML picker with real-time data
./realdetect --picker ml -s localhost -p 18000
```

### Simulator (for testing)

```bash
# Run synthetic event simulation
./realdetect_sim

# The simulator creates a test network and generates
# synthetic seismograms with realistic P and S arrivals
```

### Command Line Options

```
Usage: realdetect [options]

Options:
  -c, --config <file>      Configuration file (default: realdetect.conf)
  -s, --server <host>      SeedLink server host (default: localhost)
  -p, --port <port>        SeedLink server port (default: 18000)
  -S, --stations <file>    Station inventory file
  -v, --velocity <file>    Velocity model file
  -V, --velocity-dir <dir> Directory with regional velocity models
  -d, --database <file>    CSS3.0 database file for output
  -m, --miniseed <file>    MiniSEED file for playback mode
  -M, --miniseed-dir <dir> Directory of MiniSEED files for playback
  --playback-speed <n>     Playback speed multiplier (default: 1.0, 0 = fast)
  --picker <type>          Picker type: stalta, aic, ml (default: stalta)
  --ml-model <file>        ML model file for ML picker
  -h, --help               Show this help message
```

## Configuration

See `config/realdetect.conf` for detailed configuration options:

```ini
[picker]
# Picker type: stalta, aic, ml
type = stalta
sta_length = 0.5          # STA window (seconds)
lta_length = 10.0         # LTA window (seconds)
trigger_ratio = 3.5       # STA/LTA trigger threshold

# ML picker settings
ml_model_path = path/to/model.onnx
ml_probability_threshold = 0.5

[playback]
enabled = false           # Enable playback mode
speed = 1.0              # Playback speed (0 = fast as possible)
file = data.mseed        # MiniSEED file to playback

[locator]
algorithm = geiger        # Location algorithm (geiger, gridsearch, octtree)
fixed_depth = false       # Fix depth to default value
default_depth = 10.0      # Default/fixed depth (km)

[velocity_model]
model = simple3layer      # Built-in model (iasp91, ak135, simple3layer)
regional_models_enabled = true
regional_models_dir = config/velocity_models

[database]
enabled = true            # Enable CSS3.0 database output
file = realdetect_catalog.db
author = realdetect
network = XX
auto_create_schema = true

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

### Regional Velocity Model Configuration
Place in `config/velocity_models/velocity_models.conf`:
```ini
[default]
model_file = iasp91.vel
bounds_type = default
description = IASP91 global reference model

[socal]
model_file = socal.vel
bounds_type = box
bounds = 32.0,36.0,-121.0,-114.0
priority = 10
description = Southern California velocity model

[hawaii]
model_file = hawaii.vel
bounds_type = circle
bounds = 19.5,-155.5,200.0
priority = 20
description = Hawaii volcanic region model

[cascadia]
model_file = cascadia.vel
bounds_type = polygon
bounds = 42.0,-125.0;42.0,-120.0;49.0,-120.0;49.0,-125.0
priority = 10
description = Cascadia subduction zone model
```

### CSS3.0 Database Schema
The database follows the CSS3.0 standard with the following core tables:
- **event**: Earthquake event identification
- **origin**: Hypocenter solutions (lat, lon, depth, time)
- **origerr**: Origin error ellipse and uncertainties
- **arrival**: Phase arrival picks
- **assoc**: Arrival-origin associations with residuals
- **netmag**: Network magnitudes
- **stamag**: Station magnitudes
- **site**: Station locations
- **sitechan**: Channel information
- **affiliation**: Network-station affiliations

## Architecture

```
realdetect/
├── include/realdetect/
│   ├── core/           # Core data structures
│   │   ├── types.hpp       # Basic types, GeoPoint, StreamID
│   │   ├── waveform.hpp    # Waveform container
│   │   ├── station.hpp     # Station/channel metadata
│   │   ├── event.hpp       # Event/Origin/Pick structures
│   │   ├── miniseed.hpp    # MiniSEED parser
│   │   ├── velocity_model.hpp
│   │   ├── regional_velocity_model.hpp  # Regional model manager
│   │   └── config.hpp
│   │
│   ├── database/       # CSS3.0 database support
│   │   ├── css30_schema.hpp   # CSS3.0 table definitions
│   │   └── css30_database.hpp # SQLite implementation
│   │
│   ├── seedlink/       # SeedLink data acquisition
│   │   └── seedlink_client.hpp
│   │
│   ├── picker/         # Phase detection
│   │   ├── picker.hpp
│   │   ├── stalta_picker.hpp
│   │   ├── aic_picker.hpp
│   │   ├── ml_picker.hpp       # ML-based picker
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
├── config/
│   ├── realdetect.conf         # Main configuration
│   ├── velocity_model.txt      # Default velocity model
│   └── velocity_models/        # Regional velocity models
│       ├── velocity_models.conf
│       ├── iasp91.vel
│       ├── socal.vel
│       └── ...
│
└── src/                # Implementations
```

## Processing Pipeline

```
Data Flow:
                                                                              ┌──────────┐
  SeedLink    →   Buffer   →   Picker   →   Associator   →   Locator   →   Magnitude   →   │  CSS3.0  │
     │              │           │               │               │              │            │ Database │
  Real-time    Circular      STA/LTA        Travel-time      Geiger          ML            └──────────┘
  MiniSEED     buffer        AIC/ML         clustering      inversion        Mw
     │                                           │
  Playback                                  Regional
    Mode                                 Velocity Model
                                           Selection
```

## Algorithms

### STA/LTA Picker
The classic trigger algorithm computes the ratio of short-term to long-term energy:
```
R(t) = STA(t) / LTA(t)
Trigger when R > threshold
```

### ML Picker
Machine learning based phase detection using neural networks:
```
1. Preprocess waveform (normalize, filter)
2. Run model inference → phase probabilities P(t), S(t)
3. Find peaks above threshold
4. Return picks at peak locations
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
- Ross et al. (2018). PhaseNet: A Deep-Neural-Network-Based Seismic Arrival Time Picking Method
- Mousavi et al. (2020). EQTransformer: An Attentive Deep-Learning Model for Simultaneous Earthquake Detection and Phase Picking

## License

MIT License - see LICENSE file

## Contributing

Contributions welcome! Please submit issues and pull requests on GitHub.
