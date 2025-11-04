#!/bin/bash
# HPC-Optimized Installation script for Virall
# Designed for HPC clusters with module systems, proxy requirements, and limited sudo access

set -e  # Exit on any error

echo "Virall - HPC Installation Script"
echo "=================================="
echo ""
echo "This script is optimized for HPC cluster environments."
echo "It will:"
echo "  - Skip system package installation (no sudo required)"
echo "  - Check for conda via modules"
echo "  - Use scratch space for large databases (if available)"
echo "  - Integrate with module systems"
echo "  - Support proxy configurations"
echo ""

# Detect HPC environment
HPC_MODE=true
if [ -n "$SLURM_JOB_ID" ] || [ -n "$PBS_JOBID" ] || [ -n "$LSB_JOBID" ]; then
    echo "Detected HPC batch job environment"
fi

# Helper function to check if a command is available
check_command() {
    command -v "$1" >/dev/null 2>&1
}

# Check for system modules (common on HPC systems)
echo "Checking for system modules (HPC environments)..."
HAS_MODULE_SYSTEM=false
if command -v module >/dev/null 2>&1; then
    HAS_MODULE_SYSTEM=true
    echo "Module system detected (likely HPC environment)"
    # Try to initialize module system if not already done
    if [ -z "$MODULEPATH" ]; then
        if [ -f /etc/profile.d/modules.sh ]; then
            source /etc/profile.d/modules.sh 2>/dev/null || true
        elif [ -f /usr/share/Modules/init/bash ]; then
            source /usr/share/Modules/init/bash 2>/dev/null || true
        fi
    fi
fi

# Skip system package installation on HPC (users typically don't have sudo)
echo ""
echo "Skipping system package installation (HPC mode)"
echo "  System packages are not required - all tools will be installed via conda"
echo "  If you need system packages, contact your HPC administrator"

# Set up proxy support for downloads
echo ""
echo "Checking proxy configuration..."
CURL_PROXY_FLAGS=""
WGET_PROXY_FLAGS=""
if [ -n "$HTTP_PROXY" ] || [ -n "$HTTPS_PROXY" ]; then
    echo "Detected proxy settings - will use for downloads"
    if [ -n "$HTTP_PROXY" ]; then
        CURL_PROXY_FLAGS="-x $HTTP_PROXY"
        WGET_PROXY_FLAGS="-e use_proxy=yes -e http_proxy=$HTTP_PROXY"
    fi
    # Configure conda to use proxies
    if [ -n "$HTTP_PROXY" ]; then
        echo "Configuring conda to use proxy settings..."
    fi
fi

# Download function with retry and timeout
download_with_retry() {
    local url=$1
    local output=$2
    local max_retries=${3:-3}
    local retry_count=0
    
    while [ $retry_count -lt $max_retries ]; do
        echo "Download attempt $((retry_count + 1)) of $max_retries..."
        if curl -fL --connect-timeout 30 --max-time 3600 \
               --retry 3 --retry-delay 5 \
               $CURL_PROXY_FLAGS -o "$output" "$url" 2>&1; then
            # Validate download
            if [ -s "$output" ]; then
                echo "Download successful"
                return 0
            else
                echo "Downloaded file is empty"
            fi
        else
            echo "Download failed"
        fi
        
        retry_count=$((retry_count + 1))
        if [ $retry_count -lt $max_retries ]; then
            local wait_time=$((retry_count * 10))
            echo "Retrying in $wait_time seconds..."
            sleep $wait_time
        fi
    done
    
    echo "Download failed after $max_retries attempts"
    return 1
}

# Find conda - check multiple locations
echo ""
echo "Locating conda..."
USE_CONDA=false
CONDA_FOUND_VIA=""

# First check if conda is in PATH
if command -v conda &> /dev/null; then
    USE_CONDA=true
    CONDA_FOUND_VIA="PATH"
    echo "Conda found in PATH"
# Check if conda is available via modules
elif [ "$HAS_MODULE_SYSTEM" = true ]; then
    # Try common conda module names
    for module_name in conda anaconda miniconda python/3 anaconda3; do
        if module avail "$module_name" &> /dev/null 2>&1; then
            echo "Found conda module: $module_name"
            echo "Loading module: $module_name"
            module load "$module_name" 2>/dev/null || true
            if command -v conda &> /dev/null; then
                USE_CONDA=true
                CONDA_FOUND_VIA="module:$module_name"
                echo "Conda loaded via module: $module_name"
                break
            fi
        fi
    done
fi

# Check common installation locations
if [ "$USE_CONDA" = false ]; then
    for conda_path in "$HOME/miniconda3" "$HOME/anaconda3" "$HOME/conda" \
                      "/usr/local/miniconda3" "/usr/local/anaconda3" \
                      "/opt/conda" "/opt/miniconda3"; do
        if [ -f "$conda_path/bin/conda" ]; then
            export PATH="$conda_path/bin:$PATH"
            USE_CONDA=true
            CONDA_FOUND_VIA="path:$conda_path"
            echo "Conda found at: $conda_path"
            break
        fi
    done
fi

# If still not found, provide instructions
if [ "$USE_CONDA" = false ]; then
    echo ""
    echo "ERROR: Conda not found!"
    echo ""
    echo "On HPC clusters, conda is typically available via modules. Try:"
    echo "  module avail conda"
    echo "  module avail anaconda"
    echo "  module load conda  # or module load anaconda"
    echo ""
    echo "Alternatively, if you have Miniconda installed:"
    echo "  export PATH=\"\$HOME/miniconda3/bin:\$PATH\""
    echo ""
    echo "After loading conda, re-run this script:"
    echo "  bash install_hpc.sh"
    exit 1
fi

# Verify conda is working
if ! conda --version >/dev/null 2>&1; then
    echo "ERROR: Conda found but not working properly"
    echo "Please check your conda installation"
    exit 1
fi

echo "Conda version: $(conda --version)"

# Configure conda for proxy if needed
if [ -n "$HTTP_PROXY" ] || [ -n "$HTTPS_PROXY" ]; then
    echo "Configuring conda proxy settings..."
    if [ -n "$HTTP_PROXY" ]; then
        conda config --set proxy_servers.http "$HTTP_PROXY" 2>/dev/null || true
    fi
    if [ -n "$HTTPS_PROXY" ]; then
        conda config --set proxy_servers.https "$HTTPS_PROXY" 2>/dev/null || true
    fi
fi

# Accept conda Terms of Service automatically
echo ""
echo "Accepting conda Terms of Service..."
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main 2>/dev/null || true
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r 2>/dev/null || true

# Function to activate conda environment properly
activate_conda_env() {
    local env_name=$1
    CONDA_BASE=$(conda info --base)
    
    # Initialize conda for this shell
    if [ -f "$CONDA_BASE/etc/profile.d/conda.sh" ]; then
        source "$CONDA_BASE/etc/profile.d/conda.sh"
    elif [ -f "$CONDA_BASE/etc/profile.d/conda.csh" ]; then
        source "$CONDA_BASE/etc/profile.d/conda.csh"
    fi
    
    # Activate environment
    conda activate "$env_name" || {
        echo "Failed to activate conda environment: $env_name"
        return 1
    }
    
    # Verify activation
    if [ "$CONDA_DEFAULT_ENV" != "$env_name" ]; then
        echo "Warning: Environment may not be activated correctly"
        echo "Current environment: $CONDA_DEFAULT_ENV"
    else
        echo "Conda environment '$env_name' activated successfully"
    fi
}

# Create conda environment for viral assembler
echo ""
echo "Creating conda environment 'virall'..."
if conda env list | grep -q "^virall "; then
    echo "Environment 'virall' already exists"
    read -p "Do you want to recreate it? (y/N): " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        echo "Removing existing environment..."
        conda env remove -n virall -y
        conda create -n virall -y python=3.11
    else
        echo "Using existing environment"
    fi
else
    conda create -n virall -y python=3.11
fi

# Activate environment
echo ""
echo "Activating conda environment..."
activate_conda_env virall

# Install Python dependencies with conda
echo ""
echo "Installing Python dependencies with conda..."
conda install -c conda-forge -y numpy pandas matplotlib seaborn plotly biopython scikit-learn click tqdm pyyaml loguru psutil

# Ensure mamba is available
echo ""
echo "Checking for mamba (faster dependency solver)..."
if ! command -v mamba &> /dev/null; then
    echo "Mamba not found. Installing mamba..."
    conda install -c conda-forge mamba -y
    echo "Mamba installed successfully"
else
    echo "Mamba found"
fi

# Check for existing MPI (cluster-provided)
check_existing_mpi() {
    if [ "$HAS_MODULE_SYSTEM" = true ]; then
        for mpi_module in openmpi intel-mpi mpich mvapich2; do
            if module avail "$mpi_module" &> /dev/null 2>&1; then
                echo "Cluster MPI module detected: $mpi_module"
                return 0
            fi
        done
    fi
    # Check if MPI libraries are in system paths
    if ldconfig -p 2>/dev/null | grep -q libmpi; then
        echo "System MPI libraries detected"
        return 0
    fi
    return 1
}

# Enhanced module tool checking
check_module_tool() {
    local tool=$1
    if [ "$HAS_MODULE_SYSTEM" = true ]; then
        # Try to find module (exact name first)
        if module avail "$tool" &> /dev/null 2>&1; then
            module load "$tool" 2>/dev/null || true
            if command -v "$tool" &> /dev/null; then
                return 0  # Found via module
            fi
        fi
        # Try variations
        for variant in "$tool" "${tool^^}" "${tool,,}"; do
            if module avail "$variant" &> /dev/null 2>&1; then
                module load "$variant" 2>/dev/null || true
                if command -v "$tool" &> /dev/null; then
                    return 0
                fi
            fi
        done
    fi
    return 1  # Not found via module
}

# Define all required bioinformatics tools
declare -A TOOLS_TO_INSTALL
TOOLS_TO_INSTALL=(
    ["samtools"]="samtools"
    ["bwa"]="bwa"
    ["minimap2"]="minimap2"
    ["spades"]="spades=4.2.0"
    ["flye"]="flye=2.9.6"
    ["fastp"]="fastp"
    ["fastplong"]="fastplong=0.4.1"
    ["fastqc"]="fastqc"
    ["checkv"]="checkv"
    ["bcftools"]="bcftools"
    ["pilon"]="pilon"
    ["hmmer"]="hmmer"
    ["prodigal"]="prodigal"
    ["kaiju"]="kaiju"
)

# Check which tools are already available
echo ""
echo "Checking for existing bioinformatics tools..."
MISSING_TOOLS=()
FOUND_TOOLS=()
MODULE_TOOLS=()

for tool in "${!TOOLS_TO_INSTALL[@]}"; do
    if check_module_tool "$tool"; then
        MODULE_TOOLS+=("$tool")
        echo "  $tool found via module system - will use module"
    elif check_command "$tool"; then
        FOUND_TOOLS+=("$tool")
        echo "  $tool found in PATH - skipping installation"
    else
        MISSING_TOOLS+=("$tool")
        echo "  $tool not found - will install via conda"
    fi
done

echo ""
echo "Summary:"
echo "  - ${#MODULE_TOOLS[@]} tools found via modules"
echo "  - ${#FOUND_TOOLS[@]} tools found in PATH"
echo "  - ${#MISSING_TOOLS[@]} tools need installation via conda"

# Provide module load instructions
if [ ${#MODULE_TOOLS[@]} -gt 0 ]; then
    echo ""
    echo "IMPORTANT: These tools are available via modules. Load them before using virall:"
    echo "  module load ${MODULE_TOOLS[*]}"
    echo ""
    echo "You may want to add this to your ~/.bashrc or job submission script:"
    echo "  module load ${MODULE_TOOLS[*]}"
fi

# Install missing bioinformatics tools using mamba
if [ ${#MISSING_TOOLS[@]} -gt 0 ]; then
    echo ""
    echo "Installing missing bioinformatics tools using mamba..."
    
    # Group tools by installation method for better dependency resolution
    MAPPING_TOOLS=()
    for tool in "${MISSING_TOOLS[@]}"; do
        if [[ "$tool" == "samtools" || "$tool" == "bwa" || "$tool" == "minimap2" ]]; then
            MAPPING_TOOLS+=("${TOOLS_TO_INSTALL[$tool]}")
        fi
    done
    
    ASSEMBLY_TOOLS=()
    for tool in "${MISSING_TOOLS[@]}"; do
        if [[ "$tool" == "spades" || "$tool" == "flye" ]]; then
            ASSEMBLY_TOOLS+=("${TOOLS_TO_INSTALL[$tool]}")
        fi
    done
    
    QC_TOOLS=()
    for tool in "${MISSING_TOOLS[@]}"; do
        if [[ "$tool" == "fastp" || "$tool" == "fastplong" || "$tool" == "fastqc" ]]; then
            QC_TOOLS+=("${TOOLS_TO_INSTALL[$tool]}")
        fi
    done
    
    OTHER_TOOLS=()
    for tool in "${MISSING_TOOLS[@]}"; do
        if [[ ! "$tool" == "samtools" && ! "$tool" == "bwa" && ! "$tool" == "minimap2" && \
              ! "$tool" == "spades" && ! "$tool" == "flye" && \
              ! "$tool" == "fastp" && ! "$tool" == "fastplong" && ! "$tool" == "fastqc" && \
              ! "$tool" == "kaiju" ]]; then
            OTHER_TOOLS+=("${TOOLS_TO_INSTALL[$tool]}")
        fi
    done
    
    # Install mapping tools
    if [ ${#MAPPING_TOOLS[@]} -gt 0 ]; then
        echo "Installing mapping tools: ${MAPPING_TOOLS[*]}..."
        mamba install -c bioconda -c conda-forge "${MAPPING_TOOLS[@]}" -y || {
            echo "Warning: Failed to install mapping tools together, trying individually..."
            for tool in "${MAPPING_TOOLS[@]}"; do
                mamba install -c bioconda -c conda-forge "$tool" -y || true
            done
        }
    fi
    
    # Install assembly tools
    if [ ${#ASSEMBLY_TOOLS[@]} -gt 0 ]; then
        echo "Installing assembly tools: ${ASSEMBLY_TOOLS[*]}..."
        for tool in "${ASSEMBLY_TOOLS[@]}"; do
            if [[ "$tool" == "spades=4.2.0" ]]; then
                # Check for existing MPI before installing SPAdes
                if check_existing_mpi; then
                    echo "Using cluster MPI - SPAdes should work with it"
                else
                    echo "No MPI detected - SPAdes will install OpenMPI via conda"
                fi
                
                mamba install -c conda-forge -c bioconda spades=4.2.0 -y || true
                
                # Ensure OpenMPI libraries are accessible
                echo "Verifying OpenMPI libraries are accessible..."
                if [ ! -f "$CONDA_PREFIX/lib/libmpi.so.40" ]; then
                    OPENMPI_LIB=$(find $CONDA_PREFIX -type f -name "libmpi.so.40*" 2>/dev/null | grep -E "libmpi\.so\.40$|libmpi\.so\.40\." | head -1)
                    
                    if [ -z "$OPENMPI_LIB" ]; then
                        CONDA_BASE=$(conda info --base)
                        OPENMPI_LIB=$(find "$CONDA_BASE/pkgs" -type f -name "libmpi.so.40*" 2>/dev/null | grep -E "libmpi\.so\.40$|libmpi\.so\.40\." | head -1)
                    fi
                    
                    if [ -n "$OPENMPI_LIB" ]; then
                        OPENMPI_LIB_DIR=$(dirname "$OPENMPI_LIB")
                        echo "Found OpenMPI library at: $OPENMPI_LIB"
                        echo "Linking OpenMPI libraries to environment lib directory..."
                        mkdir -p "$CONDA_PREFIX/lib"
                        for lib in "$OPENMPI_LIB_DIR"/*.so*; do
                            if [ -f "$lib" ] || [ -L "$lib" ]; then
                                lib_name=$(basename "$lib")
                                ln -sf "$lib" "$CONDA_PREFIX/lib/$lib_name" 2>/dev/null || true
                            fi
                        done
                        echo "OpenMPI libraries linked successfully"
                    else
                        echo "Warning: OpenMPI library not found - HMMER may have issues"
                    fi
                else
                    echo "OpenMPI library already in environment lib directory"
                fi
            elif [[ "$tool" == "flye=2.9.6" ]]; then
                mamba install -c bioconda -c conda-forge flye=2.9.6 -y || true
            fi
        done
    fi
    
    # Install QC tools
    if [ ${#QC_TOOLS[@]} -gt 0 ]; then
        echo "Installing QC tools: ${QC_TOOLS[*]}..."
        mamba install -c bioconda "${QC_TOOLS[@]}" -y || {
            echo "Warning: Failed to install QC tools together, trying individually..."
            for tool in "${QC_TOOLS[@]}"; do
                mamba install -c bioconda -c conda-forge "$tool" -y || true
            done
        }
    fi
    
    # Install other tools
    if [ ${#OTHER_TOOLS[@]} -gt 0 ]; then
        echo "Installing other tools: ${OTHER_TOOLS[*]}..."
        if [ -d "$CONDA_PREFIX/lib" ]; then
            export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$LD_LIBRARY_PATH"
        fi
        mamba install -c bioconda "${OTHER_TOOLS[@]}" -y || {
            echo "Warning: Failed to install other tools together, trying individually..."
            for tool in "${OTHER_TOOLS[@]}"; do
                mamba install -c bioconda "$tool" -y || true
            done
        }
    fi
    
    # Note about kaiju
    if [[ " ${MISSING_TOOLS[@]} " =~ " kaiju " ]]; then
        echo "Note: kaiju will be installed later during database setup"
    fi
else
    echo "All required bioinformatics tools are already available!"
fi

# Fix SPAdes PATH issue
echo ""
echo "Creating SPAdes symlink..."
SPADES_PATH=$(find $CONDA_PREFIX -name "spades.py" 2>/dev/null | head -1)
if [ -n "$SPADES_PATH" ]; then
    ln -sf "$SPADES_PATH" "$CONDA_PREFIX/bin/spades"
    echo "SPAdes symlink created: $SPADES_PATH -> $CONDA_PREFIX/bin/spades"
else
    SPADES_PATH=$(which spades.py 2>/dev/null)
    if [ -n "$SPADES_PATH" ]; then
        ln -sf "$SPADES_PATH" "$CONDA_PREFIX/bin/spades"
        echo "SPAdes symlink created: $SPADES_PATH -> $CONDA_PREFIX/bin/spades"
    else
        echo "Warning: SPAdes not found, may need manual symlink"
    fi
fi

# Verify critical tools
echo ""
echo "Verifying critical bioinformatics tools installation..."
VERIFICATION_FAILED=false

if check_command samtools; then
    if samtools --version >/dev/null 2>&1; then
        echo "  samtools is working correctly"
        samtools --version 2>&1 | head -1
    else
        echo "  samtools found but not working correctly"
        echo "    Attempting to fix by reinstalling samtools..."
        mamba install -c bioconda -c conda-forge samtools --force-reinstall -y 2>/dev/null || true
        if samtools --version >/dev/null 2>&1; then
            echo "  samtools fixed successfully"
        else
            echo "  Warning: samtools still has issues, but the pipeline has a fallback mode"
            VERIFICATION_FAILED=true
        fi
    fi
else
    echo "  Warning: samtools not found"
    VERIFICATION_FAILED=true
fi

for tool in bwa minimap2; do
    if check_command "$tool"; then
        if "$tool" --version >/dev/null 2>&1 || "$tool" -h >/dev/null 2>&1; then
            echo "  $tool is working correctly"
        else
            echo "  $tool found but may not work correctly"
        fi
    else
        echo "  Warning: $tool not found"
    fi
done

if [ "$VERIFICATION_FAILED" = true ]; then
    echo ""
    echo "Warning: Some tools may not be working correctly. The pipeline may have limited functionality."
fi

# Get the software installation directory
SOFTWARE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Determine database directory (use scratch space on HPC if available)
echo ""
echo "Determining database directory..."
DB_BASE_DIR=""
if [ -n "$VIRALL_DB_DIR" ]; then
    DB_BASE_DIR="$VIRALL_DB_DIR"
    echo "Using custom database directory: $DB_BASE_DIR"
elif [ -n "$SCRATCH" ] && [ -d "$SCRATCH" ]; then
    DB_BASE_DIR="$SCRATCH/virall_databases"
    echo "Using scratch space for databases: $DB_BASE_DIR"
    echo "  (This is recommended on HPC clusters to save home directory space)"
elif [ -n "$WORK" ] && [ -d "$WORK" ]; then
    DB_BASE_DIR="$WORK/virall_databases"
    echo "Using work directory for databases: $DB_BASE_DIR"
else
    DB_BASE_DIR="$SOFTWARE_DIR/databases"
    echo "Using software directory for databases: $DB_BASE_DIR"
    echo "  (Consider using scratch space if available: export VIRALL_DB_DIR=\$SCRATCH/virall_databases)"
fi

# Check available disk space
if command -v df >/dev/null 2>&1; then
    AVAILABLE_SPACE=$(df "$DB_BASE_DIR" 2>/dev/null | tail -1 | awk '{print $4}' || echo "0")
    REQUIRED_SPACE=5368709120  # 5GB in KB
    if [ "$AVAILABLE_SPACE" -lt "$REQUIRED_SPACE" ]; then
        echo "Warning: Insufficient disk space"
        echo "  Available: $((AVAILABLE_SPACE/1024/1024))GB, Required: ~5GB"
        echo "  Installation may fail during database downloads"
        read -p "Continue anyway? (y/N): " -n 1 -r
        echo
        if [[ ! $REPLY =~ ^[Yy]$ ]]; then
            echo "Installation cancelled"
            exit 1
        fi
    else
        echo "Disk space check passed: $((AVAILABLE_SPACE/1024/1024))GB available"
    fi
fi

mkdir -p "$DB_BASE_DIR"

# Install the virall package
echo ""
echo "Installing Virall package in development mode..."
cd "$SOFTWARE_DIR"
pip install -e . || {
    echo "Warning: Failed to install virall package"
    echo "  Continuing anyway - you can install it manually later"
}

# Set up VOG database
VOG_DB_DIR="$DB_BASE_DIR/vog_db"
echo ""
echo "Setting up VOG (Viral Orthologous Groups) database..."
echo "This will download ~1GB of viral protein data"
echo "Database location: $VOG_DB_DIR"
mkdir -p "$VOG_DB_DIR"

if [ -f "$VOG_DB_DIR/vog.hmm.h3m" ]; then
    echo "VOG database already exists and is ready to use"
else
    echo "Setting up VOG database (this may take 5-10 minutes)..."
    if [ -d "$CONDA_PREFIX/lib" ]; then
        export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$LD_LIBRARY_PATH"
    fi
    
    VOG_OUTPUT=$(python -c "
import os
import sys
from virall.core.vog_annotator import VOGAnnotator

conda_prefix = os.environ.get('CONDA_PREFIX', '')
if conda_prefix:
    lib_path = os.path.join(conda_prefix, 'lib')
    if os.path.isdir(lib_path):
        current_ld_path = os.environ.get('LD_LIBRARY_PATH', '')
        new_ld_path = f'{lib_path}:{current_ld_path}' if current_ld_path else lib_path
        os.environ['LD_LIBRARY_PATH'] = new_ld_path

try:
    annotator = VOGAnnotator()
    if annotator.setup_vog_database('$VOG_DB_DIR'):
        print('SUCCESS')
        sys.exit(0)
    else:
        print('FAILED')
        sys.exit(1)
except Exception as e:
    print(f'FAILED: {e}')
    sys.exit(1)
" 2>&1) || true
    
    VOG_EXIT_CODE=$?
    if [ $VOG_EXIT_CODE -eq 0 ] && echo "$VOG_OUTPUT" | grep -q "SUCCESS"; then
        if [ -f "$VOG_DB_DIR/vog.hmm.h3m" ]; then
            echo "VOG database setup completed successfully"
        else
            echo "Warning: VOG database setup reported success but database file not found"
        fi
    else
        echo "Warning: VOG database setup failed"
        echo "$VOG_OUTPUT" | grep -i error || echo "   Check the error messages above"
        
        if echo "$VOG_OUTPUT" | grep -q "libmpi.so"; then
            echo ""
            echo "Detected MPI library issue. Attempting to fix..."
            if [ -d "$CONDA_PREFIX/lib" ]; then
                export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$LD_LIBRARY_PATH"
                if [ -f "$VOG_DB_DIR/vog.hmm" ] && [ ! -f "$VOG_DB_DIR/vog.hmm.h3m" ]; then
                    echo "Pressing VOG database manually..."
                    if LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$LD_LIBRARY_PATH" hmmpress "$VOG_DB_DIR/vog.hmm" 2>&1; then
                        if [ -f "$VOG_DB_DIR/vog.hmm.h3m" ]; then
                            echo "VOG database pressed successfully!"
                        fi
                    fi
                fi
            fi
        fi
        
        if [ ! -f "$VOG_DB_DIR/vog.hmm.h3m" ]; then
            echo ""
            echo "You can run this later with:"
            echo "  export LD_LIBRARY_PATH=\"\$CONDA_PREFIX/lib:\$LD_LIBRARY_PATH\""
            echo "  python -c \"from virall.core.vog_annotator import VOGAnnotator; VOGAnnotator().setup_vog_database('$VOG_DB_DIR')\""
        fi
    fi
fi

# Install and set up Kaiju
echo ""
echo "Installing Kaiju for viral contig classification..."
if check_command kaiju; then
    echo "Kaiju already installed, skipping installation"
else
    echo "Installing Kaiju using mamba..."
    mamba install -c bioconda kaiju -y || {
        echo "Warning: Kaiju installation failed"
    }
fi

KAIJU_DB_DIR="$DB_BASE_DIR/kaiju_db"
if check_command kaiju; then
    echo ""
    echo "Setting up Kaiju viral database..."
    echo "Database location: $KAIJU_DB_DIR"
    mkdir -p "$KAIJU_DB_DIR"
    
    echo "Downloading Kaiju viral database (this may take a while)..."
    if kaiju-makedb -s viruses -d "$KAIJU_DB_DIR" 2>&1; then
        echo "Kaiju viral database setup completed"
        
        # Clean up temporary files
        echo "Cleaning up temporary database files..."
        cd "$KAIJU_DB_DIR"
        rm -f kaiju_db_viruses.bwt kaiju_db_viruses.sa
        echo "Cleanup completed"
        
        # Verify nodes.dmp exists
        if [ ! -f "$KAIJU_DB_DIR/nodes.dmp" ]; then
            echo "Warning: nodes.dmp not found. Downloading full taxonomy database..."
            if ! kaiju-makedb -d "$KAIJU_DB_DIR" 2>&1; then
                echo "Warning: Full taxonomy database download failed"
                echo "Trying manual download of nodes.dmp..."
                cd "$KAIJU_DB_DIR"
                if download_with_retry "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz" "taxdump.tar.gz"; then
                    echo "Extracting taxonomy files..."
                    tar -xzf taxdump.tar.gz
                    rm -f taxdump.tar.gz
                    echo "Manual taxonomy download completed"
                else
                    echo "Manual download failed. Kaiju may not work properly without nodes.dmp"
                fi
            fi
        fi
    else
        echo "Warning: Kaiju viral database setup failed"
        echo "You can run this later with: kaiju-makedb -s viruses -d $KAIJU_DB_DIR"
    fi
else
    echo "Warning: Kaiju not available - skipping database setup"
fi

# Set up CheckV database
CHECKV_DB_DIR="$DB_BASE_DIR/checkv_db"
echo ""
echo "Setting up CheckV database for viral contig quality assessment..."
echo "This will download ~3GB of viral genome data"
echo "Database location: $CHECKV_DB_DIR"
echo "Note: This step may take 10-30 minutes depending on your internet connection"
mkdir -p "$CHECKV_DB_DIR"

if [ -d "$CHECKV_DB_DIR" ] && [ "$(ls -A $CHECKV_DB_DIR 2>/dev/null)" ]; then
    echo "CheckV database already exists and is ready to use"
else
    echo "Setting up CheckV database (this may take 10-30 minutes)..."
    
    MAX_RETRIES=3
    RETRY_COUNT=0
    
    while [ $RETRY_COUNT -lt $MAX_RETRIES ]; do
        echo "Attempt $((RETRY_COUNT + 1)) of $MAX_RETRIES..."
        if checkv download_database "$CHECKV_DB_DIR" 2>&1; then
            echo "CheckV database setup completed successfully"
            break
        else
            RETRY_COUNT=$((RETRY_COUNT + 1))
            if [ $RETRY_COUNT -lt $MAX_RETRIES ]; then
                echo "CheckV database download failed (attempt $RETRY_COUNT/$MAX_RETRIES)"
                echo "Retrying in 60 seconds..."
                sleep 60
            else
                echo "Warning: CheckV database setup failed after $MAX_RETRIES attempts"
                echo "Trying alternative approach..."
                
                if checkv download_database "$CHECKV_DB_DIR" --force 2>&1; then
                    echo "CheckV database setup completed with alternative method"
                    break
                else
                    echo "Alternative method also failed"
                    echo "Trying manual download..."
                    
                    CHECKV_URL="https://portal.nersc.gov/CheckV/checkv-db-v1.5.tar.gz"
                    CHECKV_TAR="$CHECKV_DB_DIR/checkv-db.tar.gz"
                    
                    if download_with_retry "$CHECKV_URL" "$CHECKV_TAR" 3; then
                        if [ -s "$CHECKV_TAR" ] && tar -tzf "$CHECKV_TAR" >/dev/null 2>&1; then
                            echo "Extracting CheckV database..."
                            cd "$CHECKV_DB_DIR"
                            tar -xzf checkv-db.tar.gz && rm -f checkv-db.tar.gz
                            echo "CheckV database manually downloaded and extracted successfully"
                            break
                        else
                            echo "Downloaded file is invalid"
                            rm -f "$CHECKV_TAR"
                        fi
                    else
                        echo "Manual download failed. CheckV database will not be available."
                        echo "You can try downloading it later with:"
                        echo "  curl -L -o checkv-db.tar.gz https://portal.nersc.gov/CheckV/checkv-db-v1.5.tar.gz"
                        echo "  tar -xzf checkv-db.tar.gz"
                    fi
                fi
            fi
        fi
    done
fi

# Clean up
echo ""
echo "Cleaning up unnecessary files and folders..."
rm -rf "$SOFTWARE_DIR/build/" 2>/dev/null || true
rm -f "$SOFTWARE_DIR/nodes.dmp" "$SOFTWARE_DIR/merged.dmp" "$SOFTWARE_DIR/names.dmp" "$SOFTWARE_DIR/taxdump.tar.gz" 2>/dev/null || true

if [ -d "$SOFTWARE_DIR/viruses" ] && [ "$(ls -A $SOFTWARE_DIR/viruses 2>/dev/null)" ]; then
    echo "Moving viruses directory to databases/kaiju_db..."
    mkdir -p "$KAIJU_DB_DIR"
    mv "$SOFTWARE_DIR/viruses"/* "$KAIJU_DB_DIR/" 2>/dev/null || true
    rmdir "$SOFTWARE_DIR/viruses" 2>/dev/null || true
fi

for dir in data results logs models; do
    if [ -d "$SOFTWARE_DIR/$dir" ] && [ -z "$(ls -A $SOFTWARE_DIR/$dir 2>/dev/null)" ]; then
        rmdir "$SOFTWARE_DIR/$dir" 2>/dev/null || true
    fi
done

echo "Cleaning up completed"

# Test installation
echo ""
echo "Testing installation..."
python -c "from virall import ViralAssembler; print('Installation successful!')" || {
    echo "Warning: Installation test failed"
}

echo ""
echo "=================================="
echo "Installation completed!"
echo "=================================="
echo ""
echo "Next steps:"
echo "1. Activate the environment: conda activate virall"
if [ ${#MODULE_TOOLS[@]} -gt 0 ]; then
    echo "2. Load required modules: module load ${MODULE_TOOLS[*]}"
    echo "3. Prepare your sequencing reads (fastq)"
    echo "4. Run virall: virall --help"
else
    echo "2. Prepare your sequencing reads (fastq)"
    echo "3. Run virall: virall --help"
fi
echo ""
echo "Databases installed:"
echo "  - VOG database: $VOG_DB_DIR (for viral gene annotation)"
echo "  - Kaiju database: $KAIJU_DB_DIR (for viral contig classification)"
echo "  - CheckV database: $CHECKV_DB_DIR (for viral contig quality assessment)"
echo ""
if [ "$DB_BASE_DIR" != "$SOFTWARE_DIR/databases" ]; then
    echo "NOTE: Databases are stored in: $DB_BASE_DIR"
    echo "Make sure this location is accessible when running virall."
    echo "You may need to set database paths in your config file or environment."
fi
echo ""
if [ "$CONDA_FOUND_VIA" = *"module"* ]; then
    echo "IMPORTANT: Remember to load the conda module before using virall:"
    echo "  module load $(echo $CONDA_FOUND_VIA | cut -d: -f2)"
fi
echo ""
echo "For help, run: virall --help"
echo "Happy virus hunt!"

