FROM ubuntu:24.04 as builder

RUN apt update && apt install -y \
    cmake \
    git \
    libboost-dev-all \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

RUN git clone https://github.com/twillis209/gps_cpp.git

RUN mkdir build && cd build && cmake .. && cmake --build .

# Stage 2: Create the final image
FROM ubuntu:24.04

WORKDIR /app

# Copy only the built binary from the builder stage
COPY --from=builder /app/build/your-binary-name .

CMD ["./your-binary-name"]
