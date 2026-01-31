/*
 * GPU Regrid C Interface Header
 * Build library with nvhpc, link from GCC-compiled code
 */

#ifndef GPU_REGRID_H
#define GPU_REGRID_H

#ifdef __cplusplus
extern "C" {
#endif

/* Initialize GPU regrid module */
void gpu_regrid_init(int *rc);

/* Finalize GPU regrid module */
void gpu_regrid_finalize(int *rc);

/* Store regridding weights (COO format input, converted to CSR internally)
 * Returns map_id for use in apply functions
 */
void gpu_regrid_store_weights(
    int n_weights,
    int src_size,
    int dst_size,
    const double *weights,
    const int *dst_indices,
    const int *src_indices,
    int *map_id,
    int *rc
);

/* Apply regridding to a single field */
void gpu_regrid_apply(
    int map_id,
    int src_size,
    int dst_size,
    const double *src_data,
    double *dst_data,
    int *rc
);

/* Apply regridding to multiple fields (batch mode - more efficient) */
void gpu_regrid_apply_batch(
    int map_id,
    int src_size,
    int dst_size,
    int n_fields,
    const double *src_data,  /* src_size x n_fields, column-major */
    double *dst_data,        /* dst_size x n_fields, column-major */
    int *rc
);

/* Check if GPU is available */
void gpu_regrid_check_gpu(int *available);

#ifdef __cplusplus
}
#endif

#endif /* GPU_REGRID_H */
