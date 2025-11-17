#!/bin/bash

# Compilation Script for COLOSS
# Compiles the code and optionally sets up GUI

set -e  # Exit on error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

# Function to print colored messages
print_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

print_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

print_step() {
    echo -e "${CYAN}==>${NC} $1"
}

# Get the directory where the script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR"

echo "=================================================="
echo "  COLOSS Compilation Script"
echo "=================================================="
echo ""

# Step 1: Detect and select Fortran compiler
print_step "Checking for Fortran compiler..."
echo ""

IFORT_FOUND=false
GFORTRAN_FOUND=false
SELECTED_COMPILER=""

# Check for ifort
if command -v ifort &> /dev/null; then
    IFORT_FOUND=true
    IFORT_VERSION=$(ifort --version 2>&1 | head -n 1)
    print_success "Intel Fortran (ifort) found: $IFORT_VERSION"
fi

# Check for gfortran
if command -v gfortran &> /dev/null; then
    GFORTRAN_FOUND=true
    GFORTRAN_VERSION=$(gfortran --version | head -n 1)
    print_success "GNU Fortran (gfortran) found: $GFORTRAN_VERSION"
fi

# Compiler selection logic
if [ "$IFORT_FOUND" = true ] && [ "$GFORTRAN_FOUND" = true ]; then
    # Both compilers found - ask user to choose
    echo ""
    print_warning "Both Intel Fortran and GNU Fortran are available!"
    echo ""
    echo "  1) Use ifort (Intel Fortran)"
    echo "  2) Use gfortran (GNU Fortran)"
    echo ""
    read -p "Select compiler [1-2]: " -n 1 -r COMPILER_CHOICE
    echo ""

    case $COMPILER_CHOICE in
        1)
            SELECTED_COMPILER="ifort"
            print_info "Using Intel Fortran (ifort)"
            ;;
        2)
            SELECTED_COMPILER="gfortran"
            print_info "Using GNU Fortran (gfortran)"
            ;;
        *)
            print_error "Invalid choice. Defaulting to gfortran."
            SELECTED_COMPILER="gfortran"
            ;;
    esac

elif [ "$IFORT_FOUND" = true ]; then
    # Only ifort found
    SELECTED_COMPILER="ifort"
    print_info "Using Intel Fortran (ifort)"

elif [ "$GFORTRAN_FOUND" = true ]; then
    # Only gfortran found
    SELECTED_COMPILER="gfortran"
    print_info "Using GNU Fortran (gfortran)"

else
    # No compiler found - need to install gfortran
    print_warning "No Fortran compiler found!"
    echo ""
    read -p "Do you want to install gfortran? (y/N): " -n 1 -r
    echo ""

    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        print_error "Fortran compiler is required to compile COLOSS. Exiting."
        exit 1
    fi

    print_info "Installing gfortran..."

    # Check for conda first (preferred method)
    if command -v conda &> /dev/null; then
        print_info "Installing gfortran via conda..."

        # Get conda base path and source it
        CONDA_BASE=$(conda info --base 2>/dev/null)
        if [ -n "$CONDA_BASE" ]; then
            source "$CONDA_BASE/etc/profile.d/conda.sh"
        fi

        conda install -y -c conda-forge gfortran
        print_success "gfortran installed via conda successfully!"

    # Check for package managers
    elif [[ "$OSTYPE" == "darwin"* ]]; then
        # macOS
        if command -v brew &> /dev/null; then
            print_info "Installing gfortran via Homebrew..."
            brew install gcc
            print_success "gfortran installed via Homebrew successfully!"
        else
            print_error "Homebrew not found. Please install Homebrew first:"
            echo "  /bin/bash -c \"\$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)\""
            exit 1
        fi

    elif [[ "$OSTYPE" == "linux-gnu"* ]]; then
        # Linux
        if command -v apt-get &> /dev/null; then
            print_info "Installing gfortran via apt-get..."
            sudo apt-get update && sudo apt-get install -y gfortran
            print_success "gfortran installed successfully!"
        elif command -v yum &> /dev/null; then
            print_info "Installing gfortran via yum..."
            sudo yum install -y gcc-gfortran
            print_success "gfortran installed successfully!"
        else
            print_error "No package manager found. Please install gfortran manually."
            exit 1
        fi

    else
        print_error "Unsupported operating system: $OSTYPE"
        exit 1
    fi

    # Verify installation
    if command -v gfortran &> /dev/null; then
        SELECTED_COMPILER="gfortran"
        GFORTRAN_VERSION=$(gfortran --version | head -n 1)
        print_success "gfortran is now available: $GFORTRAN_VERSION"
    else
        print_error "gfortran installation failed. Please install manually."
        exit 1
    fi
fi

echo ""

# Step 2: Check and configure LAPACK libraries
print_step "Checking for LAPACK libraries..."
LAPACK_FOUND=false
LAPACK_SOURCE=""
LAPACK_LIB=""

# Check for LAPACK in common locations
if [[ "$OSTYPE" == "darwin"* ]]; then
    # macOS: Check multiple sources in order of preference

    # 1. Check Homebrew installation
    if [ -d "/opt/homebrew/opt/lapack" ] && [ -f "/opt/homebrew/opt/lapack/lib/liblapack.dylib" ]; then
        LAPACK_FOUND=true
        LAPACK_SOURCE="Homebrew (/opt/homebrew/opt/lapack)"
        LAPACK_LIB="-L/opt/homebrew/opt/lapack/lib -llapack -lblas"
    elif [ -d "/usr/local/opt/lapack" ] && [ -f "/usr/local/opt/lapack/lib/liblapack.dylib" ]; then
        LAPACK_FOUND=true
        LAPACK_SOURCE="Homebrew (/usr/local/opt/lapack)"
        LAPACK_LIB="-L/usr/local/opt/lapack/lib -llapack -lblas"
    # 2. Check MacPorts installation
    elif [ -f "/opt/local/lib/liblapack.dylib" ]; then
        LAPACK_FOUND=true
        LAPACK_SOURCE="MacPorts (/opt/local/lib)"
        LAPACK_LIB="-L/opt/local/lib -llapack -lblas"
    # 3. Check Accelerate framework
    elif [ -f "/System/Library/Frameworks/Accelerate.framework/Accelerate" ] || \
         [ -f "/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/System/Library/Frameworks/Accelerate.framework/Accelerate" ]; then
        LAPACK_FOUND=true
        LAPACK_SOURCE="macOS Accelerate framework"
        LAPACK_LIB="-framework Accelerate"
    fi

    if [ "$LAPACK_FOUND" = true ]; then
        print_success "LAPACK found (via $LAPACK_SOURCE)"
    fi
else
    # Linux: check for liblapack
    if ldconfig -p 2>/dev/null | grep -q liblapack; then
        LAPACK_FOUND=true
        LAPACK_SOURCE="system ldconfig"
        LAPACK_LIB="-llapack -lblas"
        print_success "LAPACK found (via $LAPACK_SOURCE)"
    elif [ -f "/usr/lib/liblapack.so" ] || \
         [ -f "/usr/lib64/liblapack.so" ] || \
         [ -f "/usr/lib/x86_64-linux-gnu/liblapack.so" ]; then
        LAPACK_FOUND=true
        LAPACK_SOURCE="system libraries"
        LAPACK_LIB="-llapack -lblas"
        print_success "LAPACK found (via $LAPACK_SOURCE)"
    fi
fi

# If not found, offer to install via conda or package manager
if [ "$LAPACK_FOUND" = false ]; then
    print_warning "LAPACK libraries not found!"
    echo ""
    read -p "Do you want to install LAPACK? (y/N): " -n 1 -r
    echo ""
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        print_info "Installing LAPACK..."

        if [[ "$OSTYPE" == "darwin"* ]]; then
            # macOS
            if command -v brew &> /dev/null; then
                print_info "Installing LAPACK via Homebrew..."
                brew install lapack
                LAPACK_LIB="-L/opt/homebrew/opt/lapack/lib -llapack -lblas"
                print_success "LAPACK installed successfully!"
                LAPACK_FOUND=true
            else
                print_error "Homebrew not found. Please install LAPACK manually."
            fi
        else
            # Linux
            if command -v apt-get &> /dev/null; then
                print_info "Installing LAPACK via apt-get..."
                sudo apt-get update && sudo apt-get install -y liblapack-dev libblas-dev
                LAPACK_LIB="-llapack -lblas"
                print_success "LAPACK installed successfully!"
                LAPACK_FOUND=true
            elif command -v yum &> /dev/null; then
                print_info "Installing LAPACK via yum..."
                sudo yum install -y lapack-devel blas-devel
                LAPACK_LIB="-llapack -lblas"
                print_success "LAPACK installed successfully!"
                LAPACK_FOUND=true
            else
                print_error "No package manager found. Please install LAPACK manually."
            fi
        fi
    else
        print_warning "LAPACK is required for COLOSS. Using default flags (may not work)."
        LAPACK_LIB="-llapack -lblas"
    fi
fi

echo ""

# Step 3: Update Makefile with detected libraries
print_info "Configuring Makefile for your system..."

# Backup Makefile
cp Makefile Makefile.bak

# Update LAPACK library path in Makefile
if [[ "$OSTYPE" == "darwin"* ]]; then
    if [[ "$LAPACK_LIB" == *"Accelerate"* ]]; then
        sed -i.tmp "s|^LIB =.*|LIB = -framework Accelerate|" Makefile
    else
        sed -i.tmp "s|^LIB =.*|LIB = $LAPACK_LIB|" Makefile
    fi
else
    sed -i.tmp "s|^LIB =.*|LIB = $LAPACK_LIB|" Makefile
fi

# Update compiler in Makefile
sed -i.tmp2 "s|^FC =.*|FC = $SELECTED_COMPILER|" Makefile
sed -i.tmp3 "s|^F90 =.*|F90 = $SELECTED_COMPILER|" Makefile

# Clean up temporary files
rm -f Makefile.tmp Makefile.tmp2 Makefile.tmp3

print_success "Makefile configured successfully!"
echo ""

# Step 4: Clean previous build
print_step "Cleaning previous build..."
make clean 2>/dev/null || true
print_success "Clean completed!"
echo ""

# Step 5: Build adyo_v1_0 library
print_step "Building adyo_v1_0 library..."
cd adyo_v1_0
make clean 2>/dev/null || true
make
cd "$SCRIPT_DIR"
print_success "adyo_v1_0 library built successfully!"
echo ""

# Step 6: Compile COLOSS
print_step "Compiling COLOSS..."
echo ""

make

echo ""
echo "=================================================="
print_success "Compilation Completed Successfully!"
echo "=================================================="
echo ""
print_info "Executable is located at:"
echo ""
echo -e "  ${GREEN}✓${NC} COLOSS"
echo ""
print_info "To run COLOSS:"
echo -e "  ${CYAN}./COLOSS < test/6Li208Pb.in${NC}"
echo ""

# Step 7: Optional GUI setup
echo ""
print_step "Would you like to set up the GUI?"
echo ""
echo "  1) Skip GUI setup"
echo "  2) Set up GUI (creates conda environment)"
echo ""
read -p "Enter your choice [1-2]: " -n 1 -r GUI_CHOICE
echo ""
echo ""

if [ "$GUI_CHOICE" = "2" ]; then
    print_info "Setting up COLOSS GUI with conda environment..."
    echo ""

    # Check and install conda if needed
    print_info "Checking for conda installation..."
    if ! command -v conda &> /dev/null; then
        print_warning "conda not found!"
        echo ""
        read -p "Do you want to install Miniconda? (y/N): " -n 1 -r
        echo ""
        if [[ ! $REPLY =~ ^[Yy]$ ]]; then
            print_error "conda is required for GUI setup. Exiting."
            exit 1
        fi

        print_info "Installing Miniconda..."

        # Detect OS and architecture
        if [[ "$OSTYPE" == "darwin"* ]]; then
            # macOS
            if [[ $(uname -m) == "arm64" ]]; then
                MINICONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh"
            else
                MINICONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh"
            fi
        elif [[ "$OSTYPE" == "linux-gnu"* ]]; then
            # Linux
            MINICONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"
        else
            print_error "Unsupported operating system: $OSTYPE"
            exit 1
        fi

        # Download and install Miniconda
        TEMP_INSTALLER="/tmp/miniconda_installer.sh"
        print_info "Downloading Miniconda from $MINICONDA_URL..."

        # Check for download tools and use the first available one
        if command -v curl &> /dev/null; then
            curl -L -o "$TEMP_INSTALLER" "$MINICONDA_URL"
        elif command -v wget &> /dev/null; then
            wget -O "$TEMP_INSTALLER" "$MINICONDA_URL"
        else
            print_error "Neither curl nor wget found!"
            print_error "Please install one of them first:"
            echo "  macOS:  brew install curl"
            echo "  Ubuntu: sudo apt-get install curl"
            exit 1
        fi

        print_info "Running Miniconda installer..."
        bash "$TEMP_INSTALLER" -b -p "$HOME/miniconda3"
        rm "$TEMP_INSTALLER"

        # Initialize conda
        print_info "Initializing conda..."
        "$HOME/miniconda3/bin/conda" init bash
        "$HOME/miniconda3/bin/conda" init zsh 2>/dev/null || true

        # Source conda for current session
        source "$HOME/miniconda3/etc/profile.d/conda.sh"

        # Accept conda Terms of Service automatically
        print_info "Accepting conda Terms of Service..."
        conda config --set channel_priority flexible 2>/dev/null || true
        conda config --set auto_activate_base false 2>/dev/null || true

        print_success "Miniconda installed successfully!"
    else
        print_success "conda found: $(conda --version)"

        # Source conda to ensure it's available
        CONDA_BASE=$(conda info --base 2>/dev/null)
        if [ -n "$CONDA_BASE" ]; then
            source "$CONDA_BASE/etc/profile.d/conda.sh"
        fi
    fi
    echo ""

    # Create or update conda environment
    ENV_NAME="coloss_gui"
    print_info "Setting up conda environment: $ENV_NAME"

    if conda env list | grep -q "^$ENV_NAME "; then
        print_warning "Environment $ENV_NAME already exists"
        print_info "Removing and recreating environment for clean installation..."
        conda env remove -n $ENV_NAME -y
    fi

    print_info "Creating new conda environment with Python 3.10..."
    conda create -n $ENV_NAME python=3.10 -y

    print_success "Conda environment created"
    echo ""

    # Activate environment and install Python packages
    print_info "Installing Python dependencies in conda environment..."

    # Get conda base path
    CONDA_BASE=$(conda info --base)
    source "$CONDA_BASE/etc/profile.d/conda.sh"
    conda activate $ENV_NAME

    # Verify we're using the correct Python and pip
    print_info "Using Python: $(which python)"
    print_info "Python version: $(python --version)"

    # Upgrade pip first
    print_info "Upgrading pip..."
    python -m pip install --upgrade pip

    # Install requirements
    print_info "Installing requirements from coloss_gui/requirements.txt..."
    cd coloss_gui
    python -m pip install -r requirements.txt

    # Verify PySide6 installation
    print_info "Verifying PySide6 installation..."
    if python -c "import PySide6.QtWidgets; print('PySide6 version:', PySide6.__version__)" 2>/dev/null; then
        print_success "PySide6 installed and verified successfully"
    else
        print_warning "PySide6 installation verification failed!"
        print_info "Attempting to install via conda as fallback..."
        conda install -y -c conda-forge pyside6

        if python -c "import PySide6.QtWidgets" 2>/dev/null; then
            print_success "PySide6 installed via conda successfully"
        else
            print_error "Failed to install PySide6!"
            print_info "Please try manually:"
            echo "  conda activate coloss_gui"
            echo "  pip install PySide6"
            exit 1
        fi
    fi

    cd "$SCRIPT_DIR"

    print_success "Python dependencies installed successfully!"
    echo ""

    # Update run_coloss_gui.sh to use conda environment
    print_info "Updating GUI launcher script..."
    cat > run_coloss_gui.sh << 'EOF'
#!/bin/bash
# Launcher script for COLOSS GUI with conda environment

# Get the directory where this script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Source conda
CONDA_BASE=$(conda info --base 2>/dev/null)
if [ -z "$CONDA_BASE" ]; then
    # Try common conda locations
    if [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
        source "$HOME/miniconda3/etc/profile.d/conda.sh"
    elif [ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]; then
        source "$HOME/anaconda3/etc/profile.d/conda.sh"
    else
        echo "Error: conda not found!"
        echo "Please install conda or run compile.sh to set up the GUI"
        exit 1
    fi
else
    source "$CONDA_BASE/etc/profile.d/conda.sh"
fi

# Activate conda environment
conda activate coloss_gui

if [ $? -ne 0 ]; then
    echo "Error: Failed to activate coloss_gui environment"
    echo "Please run: ./compile.sh and choose option 2 to set up the GUI"
    exit 1
fi

# Change to the GUI directory
cd "$SCRIPT_DIR/coloss_gui"

# Run the GUI
echo "Starting COLOSS GUI..."
python main.py
EOF

    chmod +x run_coloss_gui.sh
    print_success "GUI launcher updated"

    echo ""
    echo "=================================================="
    print_success "GUI Setup Completed Successfully!"
    echo "=================================================="
    echo ""
    print_info "To run the COLOSS GUI:"
    echo ""
    echo -e "  ${CYAN}./run_coloss_gui.sh${NC}"
    echo ""
    echo "  or manually:"
    echo ""
    echo -e "  ${CYAN}conda activate coloss_gui${NC}"
    echo -e "  ${CYAN}cd coloss_gui && python main.py${NC}"
    echo ""
    print_warning "Note: You may need to restart your terminal or run:"
    echo -e "  ${CYAN}source ~/.bashrc${NC}  (or ~/.zshrc)"
    echo ""
fi

echo ""
print_success "All done! COLOSS is ready to use."
echo ""
