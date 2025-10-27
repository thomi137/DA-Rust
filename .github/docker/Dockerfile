FROM rust:1.90.0-slim

# Install system dependencies once and cache in the image
RUN apt-get update && apt-get install -y \
    libopenblas-dev \
    liblapack-dev \
    gfortran \
    pkg-config \
    libfontconfig1-dev \
    fontconfig \
    make \
    autoconf \
    file \
    git \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

RUN cargo install cargo-nextest --locked
