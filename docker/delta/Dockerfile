FROM tensorflow/tensorflow:2.4.0-gpu@sha256:67dbafa3e7918a3efc67db49063aebdb282a0ebb1c124b7ca0db18207e6d7a22

RUN apt-get update && \
    apt-get install -y \
        git \
        libcurl4-openssl-dev\
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt

COPY requirements-delta.txt .

RUN pip install -U pip==21.1.2 && \
    pip install -r requirements-delta.txt && \
    rm requirements-delta.txt

RUN git clone https://github.com/aidenlab/straw.git && \
    cd straw && \
    git checkout 8f6175410a5c57645cc9d5116dd11db13106d72c && \
    pip install ./pybind11_python && \
    cd .. && \
    rm -rf straw

RUN git clone https://github.com/sa501428/deploy-delta.git && \
    cd deploy-delta && \
    git checkout v2.0-encode && \
    chmod +x Deploy.py

ENV PATH=/opt/deploy-delta:$PATH
