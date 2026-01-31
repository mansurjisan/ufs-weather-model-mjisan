#!/bin/bash
#==============================================================================
# Build GPU Regrid Shared Library with NVHPC on URSA
# This creates libgpu_regrid.so that can be linked from GCC-compiled code
# Submit: sbatch build_gpu_regrid_lib.sh
#==============================================================================
#SBATCH --partition=u1-h100
#SBATCH --job-name=build_gpuregrid
#SBATCH --output=build_gpuregrid_%j.out
#SBATCH --error=build_gpuregrid_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --gres=gpu:h100:1
#SBATCH --time=00:30:00
#SBATCH --account=gpu-nos-surge
#SBATCH --qos=gpu

set -e

echo "=============================================="
echo "Building GPU Regrid Shared Library with NVHPC"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Date:   $(date)"
echo "=============================================="

# Source directory (where this script lives)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd ${SCRIPT_DIR}

# Install location
INSTALL_DIR=/scratch5/purged/Mansur.Jisan/nvhpc_stack/gpu_regrid
mkdir -p ${INSTALL_DIR}/lib
mkdir -p ${INSTALL_DIR}/include

# Load modules
module purge
module load nvhpc/24.11
module load cuda/12.8.1

echo ""
echo "Compiler: $(nvfortran --version | head -1)"
echo "Source:   ${SCRIPT_DIR}"
echo "Install:  ${INSTALL_DIR}"
echo ""

#==============================================================================
# Compile Fortran module to shared library
#==============================================================================
echo "Compiling gpu_regrid_c_interface.F90..."

# Compile with OpenACC for GPU, create position-independent code for shared lib
nvfortran -c -fPIC -acc -gpu=cc90 -Minfo=accel \
    gpu_regrid_c_interface.F90 \
    -o gpu_regrid_c_interface.o

echo ""
echo "Creating shared library libgpu_regrid.so..."

# Create shared library
# Link OpenACC runtime and CUDA libraries
nvfortran -shared -fPIC -acc -gpu=cc90 \
    gpu_regrid_c_interface.o \
    -o libgpu_regrid.so

# Also create static library for flexibility
echo "Creating static library libgpu_regrid.a..."
ar rcs libgpu_regrid.a gpu_regrid_c_interface.o

#==============================================================================
# Install
#==============================================================================
echo ""
echo "Installing libraries and headers..."

cp libgpu_regrid.so ${INSTALL_DIR}/lib/
cp libgpu_regrid.a ${INSTALL_DIR}/lib/
cp gpu_regrid.h ${INSTALL_DIR}/include/

# Copy module file if created
if [ -f "gpu_regrid_c_interface.mod" ]; then
    cp gpu_regrid_c_interface.mod ${INSTALL_DIR}/include/
fi

# Clean up object files
rm -f *.o *.mod

#==============================================================================
# Create pkg-config file for easy integration
#==============================================================================
mkdir -p ${INSTALL_DIR}/lib/pkgconfig
cat > ${INSTALL_DIR}/lib/pkgconfig/gpu_regrid.pc << EOF
prefix=${INSTALL_DIR}
exec_prefix=\${prefix}
libdir=\${exec_prefix}/lib
includedir=\${prefix}/include

Name: gpu_regrid
Description: GPU-accelerated regridding library for ESMF/CMEPS
Version: 1.0.0
Libs: -L\${libdir} -lgpu_regrid -lnvf -lacchost -laccdevice
Cflags: -I\${includedir}
EOF

#==============================================================================
# Create environment setup script
#==============================================================================
cat > ${INSTALL_DIR}/setup_gpu_regrid.sh << EOF
#!/bin/bash
# Source this file to use gpu_regrid library
export GPU_REGRID_ROOT=${INSTALL_DIR}
export LD_LIBRARY_PATH=\${GPU_REGRID_ROOT}/lib:\${LD_LIBRARY_PATH}
export LIBRARY_PATH=\${GPU_REGRID_ROOT}/lib:\${LIBRARY_PATH}
export C_INCLUDE_PATH=\${GPU_REGRID_ROOT}/include:\${C_INCLUDE_PATH}
export PKG_CONFIG_PATH=\${GPU_REGRID_ROOT}/lib/pkgconfig:\${PKG_CONFIG_PATH}

# Also need nvhpc runtime libraries
module load nvhpc/24.11 2>/dev/null || true
EOF

echo ""
echo "=============================================="
echo "GPU Regrid Library Build Complete!"
echo ""
echo "Installed to: ${INSTALL_DIR}"
echo ""
echo "Libraries:"
ls -la ${INSTALL_DIR}/lib/libgpu_regrid.*
echo ""
echo "Headers:"
ls -la ${INSTALL_DIR}/include/
echo ""
echo "To use from GCC-compiled code:"
echo "  source ${INSTALL_DIR}/setup_gpu_regrid.sh"
echo "  gcc ... -I\${GPU_REGRID_ROOT}/include -L\${GPU_REGRID_ROOT}/lib -lgpu_regrid"
echo ""
echo "At runtime, ensure nvhpc libraries are available:"
echo "  module load nvhpc/24.11"
echo "  export LD_LIBRARY_PATH=\${GPU_REGRID_ROOT}/lib:\${LD_LIBRARY_PATH}"
echo "=============================================="
