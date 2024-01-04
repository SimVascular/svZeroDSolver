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

COPY . /opt/pysvzerod

RUN mkdir -p /opt/pysvzerod-build/relwithdebinfo && \
    cd /opt/pysvzerod-build/relwithdebinfo && \
    cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo /opt/pysvzerod && \
    cmake --build . --target svzerodsolver

WORKDIR /opt/data
RUN chmod +x /opt/pysvzerod/container/profiling/entrypoint.sh
ENTRYPOINT ["/opt/pysvzerod/container/profiling/entrypoint.sh"]
