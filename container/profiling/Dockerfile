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

COPY . /opt/svzerod

RUN mkdir -p /opt/svzerod-build/relwithdebinfo && \
    cd /opt/svzerod-build/relwithdebinfo && \
    cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo /opt/svzerod && \
    cmake --build . --target svzerodsolver

WORKDIR /opt/data
RUN chmod +x /opt/svzerod/container/profiling/entrypoint.sh
ENTRYPOINT ["/opt/svzerod/container/profiling/entrypoint.sh"]
