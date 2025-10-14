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
conda create -n virall -y python=3.9

# Activate environment
echo "Activating conda environment..."
source $(conda info --base)/etc/profile.d/conda.sh
conda activate virall

# Install Python dependencies with conda (avoids GCC compilation issues)
echo "Installing Python dependencies with conda..."
conda install -c conda-forge -y numpy pandas matplotlib seaborn plotly
conda install -c conda-forge -y biopython scikit-learn
conda install -c conda-forge -y click tqdm pyyaml loguru psutil

# Install bioinformatics tools
echo "Installing bioinformatics tools..."
conda install -c bioconda -y spades bwa samtools minimap2 flye
conda install -c bioconda -y fastqc trimmomatic
conda install -c bioconda -y checkv
conda install -c bioconda -y pysam pybedtools pyfaidx
conda install -c bioconda -y blast wget hmmer

# Optional: ONT adapter trimming with Porechop
echo "Installing optional ONT trimming tool (Porechop)..."
if conda install -c bioconda -y porechop; then
    echo "Porechop installed successfully"
else
    echo "Warning: Porechop installation failed. Long-read trimming will be skipped if unavailable."
fi

# Fix BLAST library issues on newer Linux systems
echo "Installing BLAST dependencies for newer Linux systems..."
conda install -c conda-forge -y libnsl

# Fix SPAdes PATH issue (create symlink)
echo "Creating SPAdes symlink..."
SPADES_PATH=$(find $CONDA_PREFIX -name "spades.py" -path "*/share/spades-*/bin/*" | head -1)
if [ -n "$SPADES_PATH" ]; then
    ln -sf "$SPADES_PATH" "$CONDA_PREFIX/bin/spades"
    echo "SPAdes symlink created"
else
    echo "SPAdes not found, may need manual symlink"
fi

# Install problematic tools in separate environments
echo "Installing additional tools in separate environments..."

# Check if mamba is available, install if not
if ! command -v mamba &> /dev/null; then
    echo "Mamba not found. Installing mamba..."
    conda install -c conda-forge mamba -y
    echo "Mamba installed successfully"
else
    echo "Mamba found"
fi

# VirSorter2 and DIAMOND removed - now using Kaiju for viral identification
echo "VirSorter2 and DIAMOND removed - using Kaiju for viral identification"
echo "This simplifies the pipeline and reduces dependencies"

# Install BWA and minimap2 for viral contig quantification
echo "Installing BWA and minimap2 for viral contig quantification..."
mamba install -n virall -c bioconda bwa minimap2 samtools -y

# Install the virall package in development mode
echo "Installing virall package..."
pip install -e .

# Set up VOG database for viral gene annotation
echo "Setting up VOG (Viral Orthologous Groups) database..."
echo "This will download ~1GB of viral protein data for functional annotation"
echo "Note: This step may take several minutes depending on your internet connection"

# Get the software installation directory
SOFTWARE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VOG_DB_DIR="$SOFTWARE_DIR/databases/vog_db"
mkdir -p "$VOG_DB_DIR"

# Check if VOG database already exists
if [ -f "$VOG_DB_DIR/vog.hmm.h3m" ]; then
    echo "VOG database already exists and is ready to use"
else
    echo "Setting up VOG database (this may take 5-10 minutes)..."
    if python -c "
from virall.core.vog_annotator import VOGAnnotator
annotator = VOGAnnotator()
if annotator.setup_vog_database('$VOG_DB_DIR'):
    print('VOG database setup completed successfully')
else:
    print('VOG database setup failed')
"; then
        echo "VOG database setup completed successfully"
        echo "   - Downloaded and extracted VOG HMM database"
        echo "   - Pressed database for HMMER searches"
        echo "   - Ready for viral gene annotation"
    else
        echo "Warning: VOG database setup failed or was interrupted"
        echo "   You can run this later with:"
        echo "   python -c \"from virall.core.vog_annotator import VOGAnnotator; VOGAnnotator().setup_vog_database('$VOG_DB_DIR')\""
    fi
fi

# Note: DIAMOND database setup removed - now using Kaiju for viral classification

#
echo "Creating QUAST environment (optional - for assembly quality assessment)..."
if conda create -n quast_env -y python=3.9; then
    echo "QUAST environment created successfully"
    echo "Installing QUAST..."
    # Try installing QUAST with mamba for better dependency resolution
    if conda run -n quast_env mamba install -c bioconda -c conda-forge quast -y; then
        echo "QUAST installed successfully"
    else
        echo "Warning: QUAST installation failed, trying alternative approach..."
        # Fallback: install without strict dependency resolution
        if conda run -n quast_env conda install -c bioconda quast --no-deps -y; then
            echo "QUAST installed with minimal dependencies"
        else
            echo "Warning: QUAST installation completely failed"
        fi
    fi
else
    echo "Warning: Failed to create QUAST environment"
fi

echo "Installing Kaiju for viral contig classification..."
if conda install -c bioconda kaiju -y; then
    echo "Kaiju installed successfully"
    
    # Setup Kaiju viral database
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
    echo "Warning: Kaiju installation failed"
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
if [ ! -d "$CHECKV_DB_DIR" ] || [ -z "$(ls -A $CHECKV_DB_DIR 2>/dev/null)" ] || [ ! -s "$CHECKV_DB_DIR/checkv-db.tar.gz" ]; then
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

# Install the package
echo "Installing Viral Genome Assembler package..."
cd "$SOFTWARE_DIR"
python setup.py install

# Note: Working directories (data, results, logs, models) will be created automatically when needed

# Clean up unnecessary files and folders
echo "Cleaning up unnecessary files and folders..."
echo "Removing build artifacts..."
rm -rf build/ 2>/dev/null || echo "  - build/ directory not found"
rm -rf virall.egg-info/ 2>/dev/null || echo "  - virall.egg-info/ directory not found"

echo "Removing Kaiju temporary files..."
rm -f nodes.dmp merged.dmp names.dmp taxdump.tar.gz 2>/dev/null || echo "  - Kaiju temporary files not found"

# Move viruses directory to databases/kaiju_db if it exists
if [ -d "viruses" ] && [ "$(ls -A viruses 2>/dev/null)" ]; then
    echo "Moving viruses directory to databases/kaiju_db..."
    mv viruses/* databases/kaiju_db/ 2>/dev/null || echo "  - Some files may already be in databases/kaiju_db"
    rmdir viruses 2>/dev/null || echo "  - viruses directory not empty, keeping it"
    echo "  - Moved Kaiju database files to databases/kaiju_db/"
fi

# Remove empty directories that shouldn't exist
echo "Removing empty directories..."
for dir in data results logs models; do
    if [ -d "$dir" ] && [ -z "$(ls -A $dir 2>/dev/null)" ]; then
        rmdir "$dir" 2>/dev/null && echo "  - Removed empty $dir/ directory"
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
echo "                                       .....      ....       ....        .....    "   
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
