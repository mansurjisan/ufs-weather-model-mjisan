/*
 * Test program for GPU Regrid C Interface
 * Compile with GCC, link against nvhpc-built libgpu_regrid.so
 *
 * Build:
 *   gcc -o test_c_interface test_c_interface.c \
 *       -I${GPU_REGRID_ROOT}/include \
 *       -L${GPU_REGRID_ROOT}/lib -lgpu_regrid \
 *       -Wl,-rpath,${GPU_REGRID_ROOT}/lib
 *
 * Run:
 *   module load nvhpc/24.11  # For runtime libraries
 *   ./test_c_interface
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gpu_regrid.h"

#define SRC_SIZE 1000
#define DST_SIZE 500
#define N_WEIGHTS 2000  /* Sparse: ~4 weights per destination point */
#define N_FIELDS 10

int main(int argc, char *argv[]) {
    int rc, map_id, available;
    int i, j;

    printf("==============================================\n");
    printf("GPU Regrid C Interface Test\n");
    printf("==============================================\n\n");

    /* Check GPU availability */
    gpu_regrid_check_gpu(&available);
    printf("GPU available: %s\n", available ? "YES" : "NO");

    /* Initialize */
    printf("Initializing GPU regrid module...\n");
    gpu_regrid_init(&rc);
    if (rc != 0) {
        fprintf(stderr, "ERROR: gpu_regrid_init failed with rc=%d\n", rc);
        return 1;
    }
    printf("  Init successful (rc=%d)\n", rc);

    /* Create test sparse matrix (COO format) */
    printf("\nCreating test weight matrix...\n");
    printf("  Source size:      %d\n", SRC_SIZE);
    printf("  Destination size: %d\n", DST_SIZE);
    printf("  Number of weights: %d\n", N_WEIGHTS);

    double *weights = (double*)malloc(N_WEIGHTS * sizeof(double));
    int *dst_indices = (int*)malloc(N_WEIGHTS * sizeof(int));
    int *src_indices = (int*)malloc(N_WEIGHTS * sizeof(int));

    /* Create a simple averaging pattern: each dst point averages ~4 src points */
    int weight_idx = 0;
    for (i = 0; i < DST_SIZE && weight_idx < N_WEIGHTS; i++) {
        /* Map dst[i] to approximately src[2*i] region */
        int src_base = (i * SRC_SIZE) / DST_SIZE;
        int n_contrib = (N_WEIGHTS / DST_SIZE);  /* Contributors per dst point */

        for (j = 0; j < n_contrib && weight_idx < N_WEIGHTS; j++) {
            int src_idx = (src_base + j) % SRC_SIZE;
            dst_indices[weight_idx] = i + 1;  /* 1-based for Fortran */
            src_indices[weight_idx] = src_idx + 1;
            weights[weight_idx] = 1.0 / n_contrib;  /* Equal weighting */
            weight_idx++;
        }
    }

    /* Store weights */
    printf("Storing weights on GPU...\n");
    gpu_regrid_store_weights(N_WEIGHTS, SRC_SIZE, DST_SIZE,
                             weights, dst_indices, src_indices,
                             &map_id, &rc);
    if (rc != 0) {
        fprintf(stderr, "ERROR: gpu_regrid_store_weights failed with rc=%d\n", rc);
        return 1;
    }
    printf("  Weights stored (map_id=%d, rc=%d)\n", map_id, rc);

    /* Create test source data */
    double *src_data = (double*)malloc(SRC_SIZE * sizeof(double));
    double *dst_data = (double*)malloc(DST_SIZE * sizeof(double));

    for (i = 0; i < SRC_SIZE; i++) {
        src_data[i] = sin(2.0 * 3.14159 * i / SRC_SIZE);
    }

    /* Apply single-field regridding */
    printf("\nApplying single-field regridding...\n");
    gpu_regrid_apply(map_id, SRC_SIZE, DST_SIZE, src_data, dst_data, &rc);
    if (rc != 0) {
        fprintf(stderr, "ERROR: gpu_regrid_apply failed with rc=%d\n", rc);
        return 1;
    }
    printf("  Single-field regrid successful (rc=%d)\n", rc);
    printf("  Sample output: dst[0]=%.6f, dst[%d]=%.6f\n",
           dst_data[0], DST_SIZE/2, dst_data[DST_SIZE/2]);

    /* Batch regridding test */
    printf("\nApplying batch regridding (%d fields)...\n", N_FIELDS);

    double *src_batch = (double*)malloc(SRC_SIZE * N_FIELDS * sizeof(double));
    double *dst_batch = (double*)malloc(DST_SIZE * N_FIELDS * sizeof(double));

    /* Column-major layout: src_batch[i + k*SRC_SIZE] = field k, point i */
    for (int k = 0; k < N_FIELDS; k++) {
        for (i = 0; i < SRC_SIZE; i++) {
            src_batch[i + k * SRC_SIZE] = sin(2.0 * 3.14159 * i / SRC_SIZE) * (k + 1);
        }
    }

    gpu_regrid_apply_batch(map_id, SRC_SIZE, DST_SIZE, N_FIELDS,
                           src_batch, dst_batch, &rc);
    if (rc != 0) {
        fprintf(stderr, "ERROR: gpu_regrid_apply_batch failed with rc=%d\n", rc);
        return 1;
    }
    printf("  Batch regrid successful (rc=%d)\n", rc);
    printf("  Sample output field 0: dst[0]=%.6f\n", dst_batch[0]);
    printf("  Sample output field 5: dst[0]=%.6f\n", dst_batch[5 * DST_SIZE]);

    /* Finalize */
    printf("\nFinalizing GPU regrid module...\n");
    gpu_regrid_finalize(&rc);
    if (rc != 0) {
        fprintf(stderr, "ERROR: gpu_regrid_finalize failed with rc=%d\n", rc);
        return 1;
    }
    printf("  Finalize successful (rc=%d)\n", rc);

    /* Cleanup */
    free(weights);
    free(dst_indices);
    free(src_indices);
    free(src_data);
    free(dst_data);
    free(src_batch);
    free(dst_batch);

    printf("\n==============================================\n");
    printf("All tests passed!\n");
    printf("==============================================\n");

    return 0;
}
