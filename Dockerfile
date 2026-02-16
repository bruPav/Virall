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
# numpy 1.26.x is required for compatibility with both TensorFlow (<2.0) and numba (<2.4)
# This allows geNomad's neural network classification to work
RUN mamba install -n virall -c conda-forge -y \
    "numpy>=1.26.0,<2.0" \
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

# Install bioinformatics tools (versions relaxed so solver can resolve libzlib/gsl)
RUN mamba install -n virall -c bioconda -c conda-forge -y \
    samtools \
    bwa \
    minimap2 \
    spades \
    flye \
    fastp \
    fastplong \
    fastqc \
    seqtk \
    checkv \
    diamond \
    bcftools \
    pilon \
    hmmer \
    prodigal \
    kaiju \
    genomad \
    umi_tools \
    medaka \
    && mamba clean -afy

# Note: virall Python package is NOT installed in the container
# All Nextflow plotting functionality is self-contained in nextflow/bin/run_plots.py
# This simplifies the container and avoids version sync issues

# Set up directories for databases and work
ENV SOFTWARE_DIR=/opt/virall
ENV VOG_DB_DIR=${SOFTWARE_DIR}/databases/vog_db
ENV KAIJU_DB_DIR=${SOFTWARE_DIR}/databases/kaiju_db
ENV CHECKV_DB_DIR=${SOFTWARE_DIR}/databases/checkv_db
ENV CHECKV_DB=${CHECKV_DB_DIR}
ENV GENOMAD_DB_DIR=${SOFTWARE_DIR}/databases/genomad_db

RUN mkdir -p ${VOG_DB_DIR} ${KAIJU_DB_DIR} ${CHECKV_DB_DIR} ${GENOMAD_DB_DIR}

# Download VOG database (tarball may have AllVOG.hmm or individual VOG*.hmm files)
RUN cd ${VOG_DB_DIR} && \
    wget -q https://fileshare.lisc.univie.ac.at/vog/vog227/vog.hmm.tar.gz && \
    tar -xzf vog.hmm.tar.gz && \
    rm -f vog.hmm.tar.gz && \
    ( [ -f AllVOG.hmm ] && cat AllVOG.hmm > vog_all.hmm ) || find . -name '*.hmm' -exec cat {} + > vog_all.hmm && \
    hmmpress vog_all.hmm

# Download VOG annotation metadata (functional categories, virus-only flags)
RUN cd ${VOG_DB_DIR} && \
    wget -q https://fileshare.lisc.univie.ac.at/vog/vog227/vog.annotations.tsv.gz && gunzip vog.annotations.tsv.gz && \
    wget -q https://fileshare.lisc.univie.ac.at/vog/vog227/vog.virusonly.tsv.gz && gunzip vog.virusonly.tsv.gz && \
    wget -q https://fileshare.lisc.univie.ac.at/vog/vog227/vogdb.functional_categories.txt

# Download Kaiju database (viruses subset)
# Tarball may extract flat or into a subdir; flatten so .fmi and .dmp are in KAIJU_DB_DIR
RUN cd ${KAIJU_DB_DIR} && \
    wget -q https://kaiju-idx.s3.eu-central-1.amazonaws.com/2024/kaiju_db_viruses_2024-08-15.tgz -O kaiju_db_viruses.tgz && \
    tar -xzf kaiju_db_viruses.tgz && \
    rm -f kaiju_db_viruses.tgz && \
    subdir=$(find . -maxdepth 1 -type d ! -name . | head -1) && \
    if [ -n "$subdir" ]; then mv "$subdir"/* . 2>/dev/null; rmdir "$subdir" 2>/dev/null || true; fi && \
    for f in kaiju_db_viruses*.fmi kaiju_db_viruses*nodes.dmp kaiju_db_viruses*names.dmp; do \
        [ -f "$f" ] && mv "$f" "$(echo "$f" | sed 's/kaiju_db_viruses[^.]*\./kaiju_db./g; s/kaiju_db_viruses[^_]*_/kaiju_db_/g')"; \
    done 2>/dev/null || true && \
    test -n "$(find . -maxdepth 1 -name '*.fmi' | head -1)" || (echo "Kaiju DB build failed: no .fmi in ${KAIJU_DB_DIR}" && exit 1)

# Download CheckV database (v1.2+ tarball does not include .dmnd; build from .faa)
RUN cd ${CHECKV_DB_DIR} && \
    curl -fL -o checkv-db.tar.gz https://portal.nersc.gov/CheckV/checkv-db-v1.5.tar.gz && \
    tar -xzf checkv-db.tar.gz && \
    rm -f checkv-db.tar.gz && \
    if [ -d "checkv-db-v1.5" ]; then mv checkv-db-v1.5/* . && rmdir checkv-db-v1.5; fi && \
    if [ -f "genome_db/checkv_reps.faa" ] && [ ! -f "genome_db/checkv_reps.dmnd" ]; then \
        diamond makedb --in genome_db/checkv_reps.faa --db genome_db/checkv_reps; \
    fi && \
    test -f genome_db/checkv_reps.dmnd || (echo "CheckV DB build failed: genome_db/checkv_reps.dmnd missing" && exit 1)

# Download geNomad database (~4GB, for RNA viruses and eukaryotic DNA viruses)
# geNomad is used alongside CheckV: CheckV for phages, geNomad for other viruses
RUN genomad download-database ${GENOMAD_DB_DIR} && \
    test -d ${GENOMAD_DB_DIR}/genomad_db || test -f ${GENOMAD_DB_DIR}/virus_hallmark_annotation.txt || \
    (echo "geNomad DB download failed" && exit 1)

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
