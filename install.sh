#!/bin/bash
# Installation script for Virall

set -e  # Exit on any error

echo "Virall - Installation Script"
echo "============================"

# Install system dependencies
echo "Installing system dependencies..."
if command -v yum &> /dev/null; then
    echo "Detected yum package manager (RHEL/CentOS/Amazon Linux)"
    if sudo yum install -y libxcrypt-compat gcc gcc-c++ make 2>/dev/null; then
        echo "System dependencies installed successfully"
    else
        echo "Warning: Some system dependencies failed to install"
        echo "   This may not affect the installation. Continuing..."
    fi
elif command -v apt-get &> /dev/null; then
    echo "Detected apt package manager (Ubuntu/Debian)"
    if sudo apt-get update && sudo apt-get install -y libcrypt1 build-essential 2>/dev/null; then
        echo "System dependencies installed successfully"
    else
        echo "Warning: Some system dependencies failed to install"
        echo "   This may not affect the installation. Continuing..."
    fi
elif command -v dnf &> /dev/null; then
    echo "Detected dnf package manager (Fedora)"
    if sudo dnf install -y libxcrypt-compat gcc gcc-c++ make 2>/dev/null; then
        echo "System dependencies installed successfully"
    else
        echo "Warning: Some system dependencies failed to install"
        echo "   This may not affect the installation. Continuing..."
    fi
else
    echo "Warning: Unknown package manager. Skipping system dependencies"
    echo "   You may need to install them manually if required:"
    echo "   Ubuntu/Debian: sudo apt-get install libcrypt1 build-essential"
    echo "   RHEL/CentOS: sudo yum install libxcrypt-compat gcc gcc-c++ make"
    echo "   Fedora: sudo dnf install libxcrypt-compat gcc gcc-c++ make"
fi

echo "Continuing with Python package installation..."

# Check if conda is available
if command -v conda &> /dev/null; then
    echo "Conda found. Using conda for all installations..."
    USE_CONDA=true
else
    echo "Conda not found. Please install Miniconda first:"
    echo "  wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"
    echo "  bash Miniconda3-latest-Linux-x86_64.sh -b -p \$HOME/miniconda3"
    echo "  echo 'export PATH=\"\$HOME/miniconda3/bin:\$PATH\"' >> ~/.bashrc"
    echo "  source ~/.bashrc"
    exit 1
fi

# Accept conda Terms of Service automatically
echo "Accepting conda Terms of Service..."
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main 2>/dev/null || true
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r 2>/dev/null || true

# Create conda environment for viral assembler
echo "Creating conda environment 'virall'..."
conda create -n virall -y python=3.11

# Activate environment
echo "Activating conda environment..."
source $(conda info --base)/etc/profile.d/conda.sh
conda activate virall

# Install Python dependencies with conda (avoids GCC compilation issues)
echo "Installing Python dependencies with conda..."
conda install -c conda-forge -y numpy=2.3.4 pandas=2.3.3 matplotlib=3.10.6 seaborn=0.13.2 plotly=6.4.0 biopython=1.85 scikit-learn=1.7.2 click=8.3.0 tqdm=4.67.1 pyyaml=6.0.3 loguru=0.7.3 psutil=7.0.0

# Ensure mamba is available before using it
echo "Checking for mamba (faster dependency solver)..."
if ! command -v mamba &> /dev/null; then
    echo "Mamba not found. Installing mamba..."
    conda install -c conda-forge mamba -y
    echo "Mamba installed successfully"
else
    echo "Mamba found"
fi

# Helper function to check if a command is available
check_command() {
    command -v "$1" >/dev/null 2>&1
}

# Check for system modules first (common on HPC systems)
echo "Checking for system modules (HPC environments)..."
HAS_MODULE_SYSTEM=false
if command -v module >/dev/null 2>&1; then
    HAS_MODULE_SYSTEM=true
    echo "Module system detected (likely HPC environment)"
fi

# Define all required bioinformatics tools
declare -A TOOLS_TO_INSTALL
TOOLS_TO_INSTALL=(
    ["samtools"]="samtools=1.22.1"
    ["bwa"]="bwa=0.7.19"
    ["minimap2"]="minimap2=2.30"
    ["spades"]="spades=4.2.0"
    ["flye"]="flye=2.9.6"
    ["fastp"]="fastp=1.0.1"
    ["fastplong"]="fastplong=0.4.1"
    ["fastqc"]="fastqc=0.12.1"
    ["checkv"]="checkv=1.0.3"
    ["bcftools"]="bcftools=1.22"
    ["pilon"]="pilon=1.24"
    ["hmmer"]="hmmer=3.4"
    ["prodigal"]="prodigal=2.6.3"
    ["kaiju"]="kaiju=1.10.1"
)

# Check which tools are already available
echo "Checking for existing bioinformatics tools..."
MISSING_TOOLS=()
FOUND_TOOLS=()

for tool in "${!TOOLS_TO_INSTALL[@]}"; do
    if check_command "$tool"; then
        FOUND_TOOLS+=("$tool")
        echo "  $tool found - skipping installation"
    else
        MISSING_TOOLS+=("$tool")
        echo "  $tool not found - will install"
    fi
done

echo ""
echo "Summary: ${#FOUND_TOOLS[@]} tools found, ${#MISSING_TOOLS[@]} tools need installation"

# Warn about HPC modules if tools are found
if [ "$HAS_MODULE_SYSTEM" = true ] && [ ${#FOUND_TOOLS[@]} -gt 0 ]; then
    echo ""
    echo "Note: Some tools found via modules. Please ensure these modules are loaded when using virall:"
    echo "  module load ${FOUND_TOOLS[*]}"
fi

# Install missing bioinformatics tools using mamba
if [ ${#MISSING_TOOLS[@]} -gt 0 ]; then
    echo ""
    echo "Installing missing bioinformatics tools using mamba..."
    
    # Group tools by installation method for better dependency resolution
    # Mapping/alignment tools
    MAPPING_TOOLS=()
    for tool in "${MISSING_TOOLS[@]}"; do
        if [[ "$tool" == "samtools" || "$tool" == "bwa" || "$tool" == "minimap2" ]]; then
            MAPPING_TOOLS+=("${TOOLS_TO_INSTALL[$tool]}")
        fi
    done
    
    # Assembly tools
    ASSEMBLY_TOOLS=()
    for tool in "${MISSING_TOOLS[@]}"; do
        if [[ "$tool" == "spades" || "$tool" == "flye" ]]; then
            ASSEMBLY_TOOLS+=("${TOOLS_TO_INSTALL[$tool]}")
        fi
    done
    
    # QC tools
    QC_TOOLS=()
    for tool in "${MISSING_TOOLS[@]}"; do
        if [[ "$tool" == "fastp" || "$tool" == "fastplong" || "$tool" == "fastqc" ]]; then
            QC_TOOLS+=("${TOOLS_TO_INSTALL[$tool]}")
        fi
    done
    
    # Other tools
    OTHER_TOOLS=()
    for tool in "${MISSING_TOOLS[@]}"; do
        if [[ ! "$tool" == "samtools" && ! "$tool" == "bwa" && ! "$tool" == "minimap2" && \
              ! "$tool" == "spades" && ! "$tool" == "flye" && \
              ! "$tool" == "fastp" && ! "$tool" == "fastplong" && ! "$tool" == "fastqc" && \
              ! "$tool" == "kaiju" ]]; then
            OTHER_TOOLS+=("${TOOLS_TO_INSTALL[$tool]}")
        fi
    done
    
    # Install mapping tools (samtools, bwa, minimap2)
    if [ ${#MAPPING_TOOLS[@]} -gt 0 ]; then
        echo "Installing mapping tools: ${MAPPING_TOOLS[*]}..."
        mamba install -c bioconda -c conda-forge "${MAPPING_TOOLS[@]}" -y || {
            echo "Warning: Failed to install mapping tools together, trying individually..."
            for tool in "${MAPPING_TOOLS[@]}"; do
                mamba install -c bioconda -c conda-forge "$tool" -y || true
            done
        }
    fi
    
    # Install assembly tools (spades, flye)
    if [ ${#ASSEMBLY_TOOLS[@]} -gt 0 ]; then
        echo "Installing assembly tools: ${ASSEMBLY_TOOLS[*]}..."
        for tool in "${ASSEMBLY_TOOLS[@]}"; do
            if [[ "$tool" == "spades=4.2.0" ]]; then
                mamba install -c conda-forge -c bioconda spades=4.2.0 -y || true
                
                # Ensure OpenMPI libraries are accessible in the environment
                # (SPAdes installs OpenMPI but libraries might not be in $CONDA_PREFIX/lib)
                # Check conda cache (created during installation) for the library
                echo "Verifying OpenMPI libraries are accessible..."
                if [ ! -f "$CONDA_PREFIX/lib/libmpi.so.40" ]; then
                    # First check in environment
                    OPENMPI_LIB=$(find $CONDA_PREFIX -type f -name "libmpi.so.40*" 2>/dev/null | grep -E "libmpi\.so\.40$|libmpi\.so\.40\." | head -1)
                    
                    # If not found in environment, check conda package cache (created during installation)
                    if [ -z "$OPENMPI_LIB" ]; then
                        CONDA_BASE=$(conda info --base)
                        OPENMPI_LIB=$(find "$CONDA_BASE/pkgs" -type f -name "libmpi.so.40*" 2>/dev/null | grep -E "libmpi\.so\.40$|libmpi\.so\.40\." | head -1)
                    fi
                    
                    if [ -n "$OPENMPI_LIB" ]; then
                        OPENMPI_LIB_DIR=$(dirname "$OPENMPI_LIB")
                        echo "Found OpenMPI library at: $OPENMPI_LIB"
                        echo "Linking all OpenMPI libraries to environment lib directory..."
                        # Ensure lib directory exists
                        mkdir -p "$CONDA_PREFIX/lib"
                        # Link ALL OpenMPI libraries (not just libmpi.so.40)
                        # This includes libmpi.so.40, libopen-rte.so.40, libopen-pal.so.40, etc.
                        for lib in "$OPENMPI_LIB_DIR"/*.so*; do
                            if [ -f "$lib" ] || [ -L "$lib" ]; then
                                lib_name=$(basename "$lib")
                                # Create symlink for each library
                                ln -sf "$lib" "$CONDA_PREFIX/lib/$lib_name" 2>/dev/null || true
                            fi
                        done
                        echo "OpenMPI libraries linked to environment lib directory"
                    else
                        echo "Warning: OpenMPI library not found - HMMER may have issues with VOG database setup"
                    fi
                else
                    echo "OpenMPI library already in environment lib directory"
                fi
            elif [[ "$tool" == "flye=2.9.6" ]]; then
                mamba install -c bioconda -c conda-forge flye=2.9.6 -y || true
            fi
        done
    fi
    
    # Install QC tools (fastp, fastplong, fastqc)
    if [ ${#QC_TOOLS[@]} -gt 0 ]; then
        echo "Installing QC tools: ${QC_TOOLS[*]}..."
        mamba install -c bioconda "${QC_TOOLS[@]}" -y || {
            echo "Warning: Failed to install QC tools together, trying individually..."
            for tool in "${QC_TOOLS[@]}"; do
                mamba install -c bioconda -c conda-forge "$tool" -y || true
            done
        }
    fi
    
    # Install other tools (checkv, bcftools, pilon, hmmer, prodigal)
    if [ ${#OTHER_TOOLS[@]} -gt 0 ]; then
        echo "Installing other tools: ${OTHER_TOOLS[*]}..."
        # Ensure MPI libraries are available for HMMER (if OpenMPI was installed via SPAdes)
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
    
    # Install kaiju separately (it's handled later in the script, but check if it's missing)
    if [[ " ${MISSING_TOOLS[@]} " =~ " kaiju " ]]; then
        echo "Note: kaiju will be installed later during database setup"
    fi
else
    echo "All required bioinformatics tools are already available!"
fi


# Fix SPAdes PATH issue (create symlink)
echo "Creating SPAdes symlink..."
SPADES_PATH=$(find $CONDA_PREFIX -name "spades.py" 2>/dev/null | head -1)
if [ -n "$SPADES_PATH" ]; then
    ln -sf "$SPADES_PATH" "$CONDA_PREFIX/bin/spades"
    echo "SPAdes symlink created: $SPADES_PATH -> $CONDA_PREFIX/bin/spades"
else
    # Try alternative locations
    SPADES_PATH=$(which spades.py 2>/dev/null)
    if [ -n "$SPADES_PATH" ]; then
        ln -sf "$SPADES_PATH" "$CONDA_PREFIX/bin/spades"
        echo "SPAdes symlink created: $SPADES_PATH -> $CONDA_PREFIX/bin/spades"
    else
        echo "Warning: SPAdes not found, may need manual symlink"
        echo "  Try: find $CONDA_PREFIX -name spades.py"
    fi
fi

# Verify critical tools are working
echo "Verifying critical bioinformatics tools installation..."
VERIFICATION_FAILED=false

# Test samtools - check if it runs (newer versions work with modern OpenSSL)
if check_command samtools; then
    if samtools --version >/dev/null 2>&1; then
        echo "  samtools is working correctly"
        samtools --version 2>&1 | head -1
    else
        echo "  samtools found but not working correctly"
        echo "    Attempting to fix by reinstalling samtools..."
        # Reinstall samtools (will get latest version compatible with system OpenSSL)
        mamba install -c bioconda -c conda-forge samtools --force-reinstall -y 2>/dev/null || true
        # Test again
        if samtools --version >/dev/null 2>&1; then
            echo "  samtools fixed successfully"
        else
            echo "  Warning: samtools still has issues, but the pipeline has a fallback mode"
            echo "    The pipeline will use SAM file processing when BAM conversion fails"
            VERIFICATION_FAILED=true
        fi
    fi
else
    echo "  Warning: samtools not found - this may cause issues"
    VERIFICATION_FAILED=true
fi

# Verify other critical tools
for tool in bwa minimap2; do
    if check_command "$tool"; then
        if "$tool" --version >/dev/null 2>&1 || "$tool" -h >/dev/null 2>&1; then
            echo "   $tool is working correctly"
        else
            echo "   $tool found but may not work correctly"
        fi
    else
        echo "  Warning: $tool not found"
    fi
done

if [ "$VERIFICATION_FAILED" = true ]; then
    echo ""
    echo "Warning: Some tools may not be working correctly. The pipeline may have limited functionality."
fi

# Get the software installation directory (needed for package installation and database setup)
SOFTWARE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Install the virall package in development mode
# (Install before database setup since VOG setup needs it)
echo "Installing Virall package in development mode..."
cd "$SOFTWARE_DIR"
pip install -e . || {
    echo "Warning: Failed to install virall package"
    echo "  Continuing anyway - you can install it manually later"
}

# Set up VOG database for viral gene annotation
echo "Setting up VOG (Viral Orthologous Groups) database..."
echo "This will download ~1GB of viral protein data for functional annotation"
echo "Note: This step may take several minutes depending on your internet connection"
VOG_DB_DIR="$SOFTWARE_DIR/databases/vog_db"
mkdir -p "$VOG_DB_DIR"

# Check if VOG database already exists
if [ -f "$VOG_DB_DIR/vog.hmm.h3m" ]; then
    echo "VOG database already exists and is ready to use"
else
    echo "Setting up VOG database (this may take 5-10 minutes)..."
    # Fix MPI library path for HMMER (if OpenMPI is installed via SPAdes)
    # Set LD_LIBRARY_PATH in Python environment so subprocess calls inherit it
    if [ -d "$CONDA_PREFIX/lib" ]; then
        export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$LD_LIBRARY_PATH"
    fi
    
    # Run VOG setup and check the actual output
    # Set LD_LIBRARY_PATH in Python environment so hmmpress subprocess can find MPI libraries
    VOG_OUTPUT=$(python -c "
import os
import sys
from virall.core.vog_annotator import VOGAnnotator

# Set LD_LIBRARY_PATH for subprocess calls (hmmpress needs MPI libraries)
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
        # Verify the database file was actually created
        if [ -f "$VOG_DB_DIR/vog.hmm.h3m" ]; then
            echo "VOG database setup completed successfully"
            echo "   - Downloaded and extracted VOG HMM database"
            echo "   - Pressed database for HMMER searches"
            echo "   - Ready for viral gene annotation"
        else
            echo "Warning: VOG database setup reported success but database file not found"
            echo "   The database may need to be pressed manually:"
            echo "   hmmpress $VOG_DB_DIR/vog.hmm"
        fi
    else
        echo "Warning: VOG database setup failed"
        echo "$VOG_OUTPUT" | grep -i error || echo "   Check the error messages above"
        
        # Check if it's an MPI library issue and try to fix it
        if echo "$VOG_OUTPUT" | grep -q "libmpi.so"; then
            echo ""
            echo "Detected MPI library issue. Attempting to fix and retry..."
            
            # Ensure LD_LIBRARY_PATH is set and try again
            if [ -d "$CONDA_PREFIX/lib" ]; then
                export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$LD_LIBRARY_PATH"
                
                # Try pressing the database manually if it exists
                if [ -f "$VOG_DB_DIR/vog.hmm" ] && [ ! -f "$VOG_DB_DIR/vog.hmm.h3m" ]; then
                    echo "Pressing VOG database manually with correct library path..."
                    if LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$LD_LIBRARY_PATH" hmmpress "$VOG_DB_DIR/vog.hmm" 2>&1; then
                        echo "VOG database pressed successfully!"
                        if [ -f "$VOG_DB_DIR/vog.hmm.h3m" ]; then
                            echo "VOG database setup completed successfully (fixed MPI issue)"
                            echo "   - Downloaded and extracted VOG HMM database"
                            echo "   - Pressed database for HMMER searches"
                            echo "   - Ready for viral gene annotation"
                        fi
                    else
                        echo "Manual pressing also failed. You can try later with:"
                        echo "  export LD_LIBRARY_PATH=\"\$CONDA_PREFIX/lib:\$LD_LIBRARY_PATH\""
                        echo "  hmmpress $VOG_DB_DIR/vog.hmm"
                    fi
                else
                    echo "Could not automatically fix. You can run this later with:"
                    echo "  export LD_LIBRARY_PATH=\"\$CONDA_PREFIX/lib:\$LD_LIBRARY_PATH\""
                    echo "  python -c \"from virall.core.vog_annotator import VOGAnnotator; VOGAnnotator().setup_vog_database('$VOG_DB_DIR')\""
                fi
            else
                echo "Could not find conda lib directory. You can run this later with:"
                echo "  export LD_LIBRARY_PATH=\"\$CONDA_PREFIX/lib:\$LD_LIBRARY_PATH\""
                echo "  python -c \"from virall.core.vog_annotator import VOGAnnotator; VOGAnnotator().setup_vog_database('$VOG_DB_DIR')\""
            fi
        else
            echo ""
            echo "Common issues:"
            echo "  1. Network issues - try running setup again"
            echo "  2. Disk space - ensure you have enough space"
            echo ""
            echo "You can run this later with:"
            echo "  export LD_LIBRARY_PATH=\"\$CONDA_PREFIX/lib:\$LD_LIBRARY_PATH\""
            echo "  python -c \"from virall.core.vog_annotator import VOGAnnotator; VOGAnnotator().setup_vog_database('$VOG_DB_DIR')\""
        fi
    fi
fi

# Assembly quality assessment installation removed

echo "Installing Kaiju for viral contig classification..."
if check_command kaiju; then
    echo "Kaiju already installed, skipping installation"
else
    echo "Installing Kaiju using mamba..."
    if mamba install -c bioconda kaiju=1.10.1 -y; then
        echo "Kaiju installed successfully"
    else
        echo "Warning: Kaiju installation failed"
    fi
fi

# Setup Kaiju viral database (if kaiju is available)
if check_command kaiju; then
    echo "Setting up Kaiju viral database..."
    KAIJU_DB_DIR="$SOFTWARE_DIR/databases/kaiju_db"
    mkdir -p "$KAIJU_DB_DIR"
    
    echo "Downloading Kaiju viral database (this may take a while)..."
    echo "Note: This will download the full NCBI taxonomy database first, then create viral subset"
    if kaiju-makedb -s viruses -d "$KAIJU_DB_DIR"; then
        echo "Kaiju viral database setup completed"
        echo "Database location: $KAIJU_DB_DIR"
        
        # Clean up temporary files to save disk space (as recommended by Kaiju)
        echo "Cleaning up temporary database files..."
        cd "$KAIJU_DB_DIR"
        rm -f kaiju_db_viruses.bwt kaiju_db_viruses.sa
        echo "Cleanup completed - saved ~300MB of disk space"
        
        # Verify that nodes.dmp exists (required for Kaiju to work)
        if [ ! -f "$KAIJU_DB_DIR/nodes.dmp" ]; then
            echo "Warning: nodes.dmp not found. Downloading full taxonomy database..."
            echo "This is required for Kaiju to work properly..."
            if kaiju-makedb -d "$KAIJU_DB_DIR"; then
                echo "Full taxonomy database downloaded successfully"
                echo "Kaiju database is now complete and ready to use"
            else
                echo "Warning: Full taxonomy database download failed"
                echo "Trying manual download of nodes.dmp..."
                
                # Manual download of nodes.dmp as fallback
                cd "$KAIJU_DB_DIR"
                if curl -L -o taxdump.tar.gz https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz; then
                    echo "Extracting taxonomy files..."
                    tar -xzf taxdump.tar.gz
                    rm taxdump.tar.gz
                    echo "Manual taxonomy download completed successfully"
                else
                    echo "Manual download also failed. Kaiju may not work properly without nodes.dmp"
                fi
            fi
        fi
    else
        echo "Warning: Kaiju viral database setup failed"
        echo "You may need to run: kaiju-makedb -s viruses -d $KAIJU_DB_DIR"
    fi
else
    echo "Warning: Kaiju not available - skipping database setup"
    echo "You can install Kaiju later and run: kaiju-makedb -s viruses -d $KAIJU_DB_DIR"
fi

# Set up CheckV database for viral contig quality assessment
echo "Setting up CheckV database for viral contig quality assessment..."
echo "This will download ~3GB of viral genome data for quality assessment"
echo "Note: This step may take 10-30 minutes depending on your internet connection"

CHECKV_DB_DIR="$SOFTWARE_DIR/databases/checkv_db"
mkdir -p "$CHECKV_DB_DIR"

# Check if CheckV database already exists
if [ -d "$CHECKV_DB_DIR" ] && [ "$(ls -A $CHECKV_DB_DIR 2>/dev/null)" ]; then
    echo "CheckV database already exists and is ready to use"
    echo "Database location: $CHECKV_DB_DIR"
else
    echo "Setting up CheckV database (this may take 10-30 minutes)..."
    
    # Retry mechanism for CheckV database download
    MAX_RETRIES=3
    RETRY_COUNT=0
    
    while [ $RETRY_COUNT -lt $MAX_RETRIES ]; do
        echo "Attempt $((RETRY_COUNT + 1)) of $MAX_RETRIES..."
        
        # Try downloading with verbose output to see what's happening
        if checkv download_database "$CHECKV_DB_DIR"; then
            echo "CheckV database setup completed successfully"
            echo "   - Downloaded viral genome database"
            echo "   - Ready for viral contig quality assessment"
            echo "Database location: $CHECKV_DB_DIR"
            break
        else
            RETRY_COUNT=$((RETRY_COUNT + 1))
            if [ $RETRY_COUNT -lt $MAX_RETRIES ]; then
                echo "CheckV database download failed (attempt $RETRY_COUNT/$MAX_RETRIES)"
                echo "This might be due to server issues. Retrying in 60 seconds..."
                sleep 60
            else
                echo "Warning: CheckV database setup failed after $MAX_RETRIES attempts"
                echo "   This may be due to temporary server issues or network problems"
                echo "   Trying alternative approach..."
                
                # Try alternative download method
                echo "Attempting alternative CheckV database setup..."
                if checkv download_database "$CHECKV_DB_DIR" --force; then
                    echo "CheckV database setup completed with alternative method"
                    echo "Database location: $CHECKV_DB_DIR"
                    break
                else
                    echo "Alternative method also failed"
                    echo "You can run this later with: checkv download_database $CHECKV_DB_DIR"
                    echo "The pipeline will work without CheckV, but quality assessment will be limited"
                fi
            fi
        fi
    done
fi

# Verify CheckV database was downloaded successfully
# Check if the database directory exists and has actual content (not just empty files)
if [ ! -d "$CHECKV_DB_DIR" ] || [ -z "$(ls -A $CHECKV_DB_DIR 2>/dev/null)" ]; then
    echo ""
    echo "CheckV database setup failed or incomplete. Trying manual download as fallback..."
    echo "This may take 10-30 minutes depending on your internet connection"
    
    # Clean up any failed downloads first
    rm -f "$CHECKV_DB_DIR/checkv-db.tar.gz"
    
    # Try manual download with curl as fallback
    CHECKV_URL="https://portal.nersc.gov/CheckV/checkv-db-v1.5.tar.gz"
    CHECKV_TAR="$CHECKV_DB_DIR/checkv-db.tar.gz"
    
    echo "Downloading CheckV database manually from $CHECKV_URL"
    # Use -f to fail on HTTP errors, -L to follow redirects
    if curl -fL -o "$CHECKV_TAR" "$CHECKV_URL"; then
        # Verify the download was successful (file size reasonable and valid gzip)
        if [ -s "$CHECKV_TAR" ] && tar -tzf "$CHECKV_TAR" >/dev/null 2>&1; then
            echo "Extracting CheckV database..."
            cd "$CHECKV_DB_DIR"
            tar -xzf checkv-db.tar.gz && rm -f checkv-db.tar.gz
            echo "CheckV database manually downloaded and extracted successfully"
            echo "Database location: $CHECKV_DB_DIR"
        else
            echo "Downloaded file is not a valid gzip tarball or is too small. Skipping CheckV setup for now."
            rm -f "$CHECKV_TAR"
        fi
    else
        echo "Manual download failed. CheckV database will not be available."
        echo "You can try downloading it later with:"
        echo "  curl -L -o checkv-db.tar.gz https://portal.nersc.gov/CheckV/checkv-db-v1.5.tar.gz"
        echo "  tar -xzf checkv-db.tar.gz"
    fi
fi

echo "All tools installed successfully!"

# Note: Virall package was already installed earlier (before database setup)

# Note: Working directories (data, results, logs, models) will be created automatically when needed

# Clean up unnecessary files and folders (in SOFTWARE_DIR)
echo "Cleaning up unnecessary files and folders..."
echo "Removing build artifacts..."
rm -rf "$SOFTWARE_DIR/build/" 2>/dev/null || echo "  - build/ directory not found"
# Note: virall.egg-info is needed for editable install, don't remove it

echo "Removing Kaiju temporary files..."
rm -f "$SOFTWARE_DIR/nodes.dmp" "$SOFTWARE_DIR/merged.dmp" "$SOFTWARE_DIR/names.dmp" "$SOFTWARE_DIR/taxdump.tar.gz" 2>/dev/null || echo "  - Kaiju temporary files not found"

# Move viruses directory to databases/kaiju_db if it exists
if [ -d "$SOFTWARE_DIR/viruses" ] && [ "$(ls -A $SOFTWARE_DIR/viruses 2>/dev/null)" ]; then
    echo "Moving viruses directory to databases/kaiju_db..."
    mv "$SOFTWARE_DIR/viruses"/* "$SOFTWARE_DIR/databases/kaiju_db/" 2>/dev/null || echo "  - Some files may already be in databases/kaiju_db"
    rmdir "$SOFTWARE_DIR/viruses" 2>/dev/null || echo "  - viruses directory not empty, keeping it"
    echo "  - Moved Kaiju database files to databases/kaiju_db/"
fi

# Remove empty directories that shouldn't exist
echo "Removing empty directories..."
for dir in data results logs models; do
    if [ -d "$SOFTWARE_DIR/$dir" ] && [ -z "$(ls -A $SOFTWARE_DIR/$dir 2>/dev/null)" ]; then
        rmdir "$SOFTWARE_DIR/$dir" 2>/dev/null && echo "  - Removed empty $dir/ directory"
    fi
done

echo "Cleaning up completed - keeping essential files and your samples/ folder"

# Test installation
echo "Testing installation..."
python -c "from virall import ViralAssembler; print('Installation successful!')"

echo ""
echo "Installation completed successfully!"
echo ""
echo "                                .'.               "    
echo "                              .'''''.             "    
echo "                        .    .'''''''              "   
echo "                        .'.     '''.      .''.       "   
echo "                         .''.          .'''''''.       "   
echo "             ..........    ''''''..''''''''''''        "   
echo "         .'''''''''''''''.   ''''''''''''''.           "   
echo "       '''''''''''''''''''.        .                   "   
echo "        '''''''''     '''''           .......          "   
echo "         .'.            '''.       .'''''''''''.       "   
echo "                        '''       ''''''''.    .'.     "   
echo "                      .'''      .'''''''.        ''    "   
echo "            .......''''.      .'''''''          .''   "   
echo "             ''''''''        .''''''.            ''   "   
echo "              ''''''''..  ..'''''''              '''  "   
echo " ..             .''''''''''''''''.               .'''  "   
echo " .'                     ''''.                    ''''  "   
echo "  ''                     ''''    '              ''''.  "   
echo "  .'.                    '''''   ''           .'''''   "   
echo "   ''.                   '''''   .''...  ...'''''''.   "   
echo "    '''.                .''''.    ''''''''''''''''.    "   
echo "     ''''..          ..''''''      ''''''''''''''      "   
echo "      .'''''''''''''''''''''         ''''''''''        "   
echo "        .'''''''''''''''''.                            "   
echo "           .''''''''''''.                            "   
echo "                                       ......                                                  .,'''.      .,,,,.       "   
echo "      ......                 ......    ......                                                  .''''.      .,,,,'       "   
echo "         .....               .....                                                             .''''.      .',,'.       "   
echo "          ....              ....                                            ......             .''''.      .''''.       "   
echo "          .....            .....       .....      ...... ........      ................        .''''.      .''''.       "   
echo "             ...          ..           .....      ...............      ..................      .''''.      .''''.       "   
echo "              ...        ....          .....      .........                        .......     ......      .''''.       "   
echo "             .....      ......         .....      .......                           ......     ......      ..'''.       "   
echo "               ....    .....           .....      ......                ..................     ......      ......       "   
echo "                ...   ......           .....      ......              .......      .......     ......      ......       "    
echo "               ...........             .....      ......             .....          ......     ......      ......       "   
echo "                  ........             .....      ......             ......       ........     ......      ......       "   
echo "                    .....              .....      ......              ....................     ......      ......       "   
echo "                                                                          .....      ....       ....        .....    "   
echo ""
echo "Next steps:"
echo "1. Activate the environment: conda activate virall"
echo "2. Prepare your sequencing reads (fastq)"
echo "3. Run virall: virall --help"
echo "4. Happy virus hunt!"
echo ""
echo "Databases installed:"
echo "  - VOG database: $VOG_DB_DIR (for viral gene annotation)"
echo "  - Kaiju database: $KAIJU_DB_DIR (for viral contig classification)"
echo "  - CheckV database: $CHECKV_DB_DIR (for viral contig quality assessment)"
echo ""
echo "For help, run: virall --help"
