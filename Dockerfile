FROM ubuntu:24.04 as builder

RUN apt update && apt install -y \
    cmake \
    git \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

RUN git clone -b docker --single-branch https://github.com/twillis209/gps_cpp.git

RUN mkdir gps_cpp/build

WORKDIR /app/gps_cpp/build

RUN cmake ..

RUN make

WORKDIR /app/gps_cpp/build/apps
