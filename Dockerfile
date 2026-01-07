# Dockerfile for Virall - Viral Assembler Pipeline
# Uses miniforge3 with mamba for fast, reliable package management

FROM condaforge/miniforge3:latest

LABEL maintainer="Virall Team"
LABEL description="Virall - Viral Assembly Pipeline"
LABEL version="1.0"

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive
ENV MAMBA_NO_BANNER=1
ENV CONDA_REMOTE_READ_TIMEOUT_SECS=300
ENV CONDA_REMOTE_CONNECT_TIMEOUT_SECS=60
ENV CONDA_REMOTE_MAX_RETRIES=5

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcrypt1 \
    build-essential \
    curl \
    wget \
    ca-certificates \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Create virall environment with Python 3.11
RUN mamba create -n virall -y python=3.11 && mamba clean -afy

# Activate virall environment for subsequent commands
SHELL ["mamba", "run", "-n", "virall", "/bin/bash", "-c"]

# Install Python dependencies
RUN mamba install -n virall -c conda-forge -y \
    "numpy>=2.0.0" \
    "pandas>=2.2.0" \
    "matplotlib>=3.8.0" \
    "seaborn>=0.13.0" \
    "plotly>=5.18.0" \
    "biopython>=1.80" \
    "scikit-learn>=1.4.0" \
    "click>=8.1.0" \
    "tqdm>=4.66.0" \
    "pyyaml>=6.0.0" \
    "loguru>=0.7.0" \
    "psutil>=5.9.0" \
    && mamba clean -afy

# Install bioinformatics tools
RUN mamba install -n virall -c bioconda -c conda-forge -y \
    samtools=1.22.1 \
    bwa=0.7.19 \
    minimap2=2.30 \
    spades=4.2.0 \
    flye=2.9.6 \
    fastp=1.0.1 \
    fastplong=0.4.1 \
    fastqc=0.12.1 \
    seqtk=1.3 \
    checkv=1.0.3 \
    bcftools=1.22 \
    pilon=1.24 \
    hmmer=3.4 \
    prodigal=2.6.3 \
    kaiju=1.10.1 \
    gsl=2.6 \
    && mamba clean -afy

# Set up directories for databases and work
ENV SOFTWARE_DIR=/opt/virall
ENV VOG_DB_DIR=${SOFTWARE_DIR}/databases/vog_db
ENV KAIJU_DB_DIR=${SOFTWARE_DIR}/databases/kaiju_db
ENV CHECKV_DB_DIR=${SOFTWARE_DIR}/databases/checkv_db
ENV CHECKV_DB=${CHECKV_DB_DIR}

RUN mkdir -p ${VOG_DB_DIR} ${KAIJU_DB_DIR} ${CHECKV_DB_DIR}

# Download VOG database
RUN cd ${VOG_DB_DIR} && \
    wget -q https://fileshare.lisc.univie.ac.at/vog/vog227/vog.hmm.tar.gz && \
    tar -xzf vog.hmm.tar.gz && \
    rm -f vog.hmm.tar.gz && \
    cat AllVOG.hmm > vog_all.hmm && \
    rm -f AllVOG.hmm && \
    hmmpress vog_all.hmm

# Download Kaiju database (viruses subset)
RUN cd ${KAIJU_DB_DIR} && \
    wget -q https://kaiju-idx.s3.eu-central-1.amazonaws.com/2024/kaiju_db_viruses_2024-10-04.tgz && \
    tar -xzf kaiju_db_viruses_2024-10-04.tgz && \
    rm -f kaiju_db_viruses_2024-10-04.tgz && \
    for f in kaiju_db_viruses*.fmi kaiju_db_viruses*nodes.dmp kaiju_db_viruses*names.dmp; do \
        [ -f "$f" ] && mv "$f" "$(echo $f | sed 's/kaiju_db_viruses[^.]*\./kaiju_db./g')"; \
    done 2>/dev/null || true

# Download CheckV database
RUN cd ${CHECKV_DB_DIR} && \
    curl -fL -o checkv-db.tar.gz https://portal.nersc.gov/CheckV/checkv-db-v1.5.tar.gz && \
    tar -xzf checkv-db.tar.gz && \
    rm -f checkv-db.tar.gz && \
    if [ -d "checkv-db-v1.5" ]; then mv checkv-db-v1.5/* . && rmdir checkv-db-v1.5; fi

# Create activation script for LD_LIBRARY_PATH
RUN mkdir -p /opt/conda/envs/virall/etc/conda/activate.d && \
    echo '#!/bin/bash' > /opt/conda/envs/virall/etc/conda/activate.d/virall_env.sh && \
    echo 'if [ -n "$CONDA_PREFIX" ] && [ -d "$CONDA_PREFIX/lib" ]; then' >> /opt/conda/envs/virall/etc/conda/activate.d/virall_env.sh && \
    echo '    export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"' >> /opt/conda/envs/virall/etc/conda/activate.d/virall_env.sh && \
    echo 'fi' >> /opt/conda/envs/virall/etc/conda/activate.d/virall_env.sh && \
    chmod +x /opt/conda/envs/virall/etc/conda/activate.d/virall_env.sh

# Set working directory
WORKDIR /data

# Set the default environment
ENV PATH="/opt/conda/envs/virall/bin:$PATH"
ENV LD_LIBRARY_PATH="/opt/conda/envs/virall/lib:${LD_LIBRARY_PATH:-}"

# Default command
CMD ["bash", "-c", "source /opt/conda/etc/profile.d/mamba.sh && mamba activate virall && bash"]
