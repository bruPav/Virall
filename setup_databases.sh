#!/bin/bash
# Setup Virall Databases
# This script downloads and sets up all databases for Virall
# Can be run on HPC after container is deployed

set -e

echo "Virall Database Setup Script"
echo "============================"
echo ""

# Default database location (can be overridden)
DB_BASE_DIR="${VIRALL_DB_DIR:-$PWD/databases}"
VOG_DB_DIR="$DB_BASE_DIR/vog_db"
KAIJU_DB_DIR="$DB_BASE_DIR/kaiju_db"
CHECKV_DB_DIR="$DB_BASE_DIR/checkv_db"

# Check if container is available
CONTAINER="${VIRALL_CONTAINER:-virall.sif}"
if [ ! -f "$CONTAINER" ]; then
    echo "ERROR: Container not found: $CONTAINER"
    echo "Please set VIRALL_CONTAINER environment variable or place virall.sif in current directory"
    exit 1
fi

echo "Using container: $CONTAINER"
echo "Database location: $DB_BASE_DIR"
echo ""

# Create database directories
mkdir -p "$VOG_DB_DIR"
mkdir -p "$KAIJU_DB_DIR"
mkdir -p "$CHECKV_DB_DIR"

# Check disk space
echo "Checking disk space..."
AVAILABLE_SPACE=$(df "$DB_BASE_DIR" 2>/dev/null | tail -1 | awk '{print $4}' || echo "0")
REQUIRED_SPACE=5368709120  # 5GB in KB
if [ "$AVAILABLE_SPACE" != "0" ] && [ "$AVAILABLE_SPACE" -lt "$REQUIRED_SPACE" ]; then
    echo "Warning: Insufficient disk space"
    echo "  Available: $((AVAILABLE_SPACE/1024/1024))GB, Required: ~5GB"
    read -p "Continue anyway? (y/N): " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        exit 1
    fi
fi

# Set up VOG database
echo ""
echo "Setting up VOG database (~1GB)..."
echo "This may take 5-10 minutes..."
if [ -f "$VOG_DB_DIR/vog.hmm.h3m" ]; then
    echo "VOG database already exists"
else
    singularity exec \
        --bind "$DB_BASE_DIR:/databases" \
        --env VOG_DB_DIR=/databases/vog_db \
        "$CONTAINER" \
        python -c "
import os
import sys
from virall.core.vog_annotator import VOGAnnotator

# Set LD_LIBRARY_PATH for subprocess calls
os.environ['LD_LIBRARY_PATH'] = '/opt/conda/lib:' + os.environ.get('LD_LIBRARY_PATH', '')

try:
    annotator = VOGAnnotator()
    if annotator.setup_vog_database('/databases/vog_db'):
        print('VOG database setup completed successfully')
        sys.exit(0)
    else:
        print('VOG database setup failed')
        sys.exit(1)
except Exception as e:
    print(f'VOG database setup error: {e}')
    sys.exit(1)
" || {
        echo "Warning: VOG database setup failed"
        echo "You can try manually later"
    }
fi

# Set up Kaiju database
echo ""
echo "Setting up Kaiju database (this may take 30-60 minutes)..."
if [ -f "$KAIJU_DB_DIR/kaiju_db_viruses.fmi" ]; then
    echo "Kaiju database already exists"
else
    singularity exec \
        --bind "$DB_BASE_DIR:/databases" \
        "$CONTAINER" \
        bash -c "cd /databases/kaiju_db && kaiju-makedb -s viruses -d /databases/kaiju_db" || {
        echo "Warning: Kaiju database setup failed"
        echo "Trying alternative method..."
        singularity exec \
            --bind "$DB_BASE_DIR:/databases" \
            "$CONTAINER" \
            bash -c "cd /databases/kaiju_db && kaiju-makedb -s viruses -d /databases/kaiju_db --force" || {
            echo "Kaiju database setup failed - you may need to install manually"
        }
    }
    
    # Clean up temporary files
    if [ -d "$KAIJU_DB_DIR" ]; then
        cd "$KAIJU_DB_DIR"
        rm -f kaiju_db_viruses.bwt kaiju_db_viruses.sa 2>/dev/null || true
    fi
    
    # Verify nodes.dmp exists
    if [ ! -f "$KAIJU_DB_DIR/nodes.dmp" ]; then
        echo "Downloading full Kaiju taxonomy database..."
        singularity exec \
            --bind "$DB_BASE_DIR:/databases" \
            "$CONTAINER" \
            bash -c "cd /databases/kaiju_db && kaiju-makedb -d /databases/kaiju_db" || {
            echo "Trying manual download of taxonomy files..."
            cd "$KAIJU_DB_DIR"
            curl -L -o taxdump.tar.gz https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz || \
            wget -O taxdump.tar.gz https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz || true
            if [ -f taxdump.tar.gz ]; then
                tar -xzf taxdump.tar.gz
                rm -f taxdump.tar.gz
            fi
        }
    fi
fi

# Set up CheckV database
echo ""
echo "Setting up CheckV database (~3GB)..."
echo "This may take 30-60 minutes..."
if [ -d "$CHECKV_DB_DIR" ] && [ "$(ls -A $CHECKV_DB_DIR 2>/dev/null)" ]; then
    echo "CheckV database already exists"
else
    MAX_RETRIES=3
    RETRY_COUNT=0
    
    while [ $RETRY_COUNT -lt $MAX_RETRIES ]; do
        echo "Attempt $((RETRY_COUNT + 1)) of $MAX_RETRIES..."
        if singularity exec \
            --bind "$DB_BASE_DIR:/databases" \
            "$CONTAINER" \
            checkv download_database /databases/checkv_db 2>&1; then
            echo "CheckV database setup completed successfully"
            break
        else
            RETRY_COUNT=$((RETRY_COUNT + 1))
            if [ $RETRY_COUNT -lt $MAX_RETRIES ]; then
                echo "CheckV download failed, retrying in 60 seconds..."
                sleep 60
            else
                echo "Trying alternative method..."
                if singularity exec \
                    --bind "$DB_BASE_DIR:/databases" \
                    "$CONTAINER" \
                    checkv download_database /databases/checkv_db --force 2>&1; then
                    echo "CheckV database setup completed with alternative method"
                    break
                else
                    echo "Trying manual download..."
                    CHECKV_URL="https://portal.nersc.gov/CheckV/checkv-db-v1.5.tar.gz"
                    CHECKV_TAR="$CHECKV_DB_DIR/checkv-db.tar.gz"
                    
                    curl -L -o "$CHECKV_TAR" "$CHECKV_URL" || \
                    wget -O "$CHECKV_TAR" "$CHECKV_URL" || true
                    
                    if [ -f "$CHECKV_TAR" ] && [ -s "$CHECKV_TAR" ]; then
                        cd "$CHECKV_DB_DIR"
                        tar -xzf checkv-db.tar.gz && rm -f checkv-db.tar.gz
                        echo "CheckV database manually downloaded and extracted"
                        break
                    else
                        echo "CheckV database setup failed after all attempts"
                        echo "You can try downloading manually later"
                    fi
                fi
            fi
        fi
    done
fi

echo ""
echo "=================================="
echo "Database setup completed!"
echo "=================================="
echo ""
echo "Databases installed at:"
echo "  - VOG: $VOG_DB_DIR"
echo "  - Kaiju: $KAIJU_DB_DIR"
echo "  - CheckV: $CHECKV_DB_DIR"
echo ""
echo "To use these databases with the container:"
echo "  singularity run \\"
echo "    --bind $DB_BASE_DIR:/opt/virall/databases \\"
echo "    $CONTAINER \\"
echo "    analyse --short-reads-1 reads.fq -o output/"
echo ""
echo "Or set environment variable:"
echo "  export VIRALL_DB_DIR=$DB_BASE_DIR"

