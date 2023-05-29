FROM ubuntu:22.04

# Install dependencies
RUN apt-get -q update && \
    DEBIAN_FRONTEND=noninteractive \
    TZ=Europe/Berlin \
    apt-get -q -y --no-install-recommends install \
    git \
    gcc \
    g++ \
    make \
    cmake \
    python3 \
    python3-dev \
    ca-certificates \
    google-perftools \
    libgoogle-perftools-dev \
    graphviz \
    ghostscript \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

COPY . /opt/svzerodplus

RUN mkdir -p /opt/svzerodplus-build/relwithdebinfo && \
    cd /opt/svzerodplus-build/relwithdebinfo && \
    cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo /opt/svzerodplus && \
    cmake --build . --target svzerodsolver

WORKDIR /opt
ENTRYPOINT ["/opt/svzerodplus/container/profiling/entrypoint.sh"]
