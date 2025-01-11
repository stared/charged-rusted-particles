#!/bin/bash

# Source cargo environment
source "$HOME/.cargo/env"

# Add wasm target
rustup target add wasm32-unknown-unknown

# Build the wasm package
wasm-pack build --target web

# Create build directory
mkdir -p build
cp index.html build/
cp -r pkg build/

# Find an available port starting from 8000
port=8000
while nc -z localhost $port 2>/dev/null; do
    port=$((port + 1))
    if [ $port -gt 8010 ]; then
        echo "No available ports found between 8000 and 8010"
        exit 1
    fi
done

echo "Starting server on port $port"
cd build
basic-http-server -a 127.0.0.1:$port . 