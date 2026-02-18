# Multi-stage build for RealDetect
# Stage 1: Build the C++ application
FROM ubuntu:22.04 AS builder

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    libeigen3-dev \
    libsqlite3-dev \
    libgomp1 \
    pkg-config \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app
COPY CMakeLists.txt .
COPY include/ include/
COPY src/ src/
COPY tests/ tests/
COPY config/ config/

RUN mkdir build && cd build \
    && cmake .. -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTS=ON \
    && make -j$(nproc)

# Stage 2: Runtime image with both C++ app and Python web GUI
FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    libsqlite3-0 \
    libgomp1 \
    libeigen3-dev \
    python3 \
    python3-pip \
    && rm -rf /var/lib/apt/lists/*

RUN pip3 install --no-cache-dir \
    fastapi \
    uvicorn[standard] \
    aiosqlite \
    jinja2

WORKDIR /app

# Copy built binaries
COPY --from=builder /app/build/realdetect /app/bin/realdetect
COPY --from=builder /app/build/realdetect_sim /app/bin/realdetect_sim
COPY --from=builder /app/build/realdetect_tests /app/bin/realdetect_tests

# Copy config and web GUI
COPY config/ /app/config/
COPY gui/ /app/gui/

ENV PATH="/app/bin:${PATH}"

EXPOSE 8080

CMD ["python3", "/app/gui/server.py"]
