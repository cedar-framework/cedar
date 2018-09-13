#ifndef CEDAR_CAPI_H
#define CEDAR_CAPI_H

#include <mpi.h>

#include <cedar/ctypes.h>

#ifdef __cplusplus
extern "C" {
#endif

#define CEDAR_SUCCESS     0
#define CEDAR_ERR_TOPO    1
#define CEDAR_ERR_KIND    2
#define CEDAR_ERR_HANDLE  3
#define CEDAR_ERR_MAT     4
#define CEDAR_ERR_VEC     5
#define CEDAR_ERR_DIM     6
#define CEDAR_ERR_STENCIL 7
#define CEDAR_ERR_SOLVER  8
#define CEDAR_ERR_OTHER   9

typedef int cedar_topo;
typedef int cedar_mat;
typedef int cedar_vec;
typedef int cedar_solver;

#define CEDAR_TOPO_NULL ((cedar_topo) 1)
#define CEDAR_MAT_NULL ((cedar_mat) 3)
#define CEDAR_SOLVER_NULL ((cedar_solver) 2)
#define CEDAR_VEC_NULL ((cedar_vec) 4)

/**
 * Create 2D distributed grid topology.
 *
 * @param[in]  comm MPI Communicator for distributed grid
 * @param[in]  ngx number of global grid points in first dimension
 * @param[in]  ngy number of global grid points in second dimension
 * @param[in]  lnx array of length nprocx of local grid points for each processor
 * @param[in]  lny array of length nprocy of local grid points for each processor
 * @param[in]  nprocx number of processors in first dimension
 * @param[in]  nprocy number of processors in second dimension
 * @param[out] newtopo new topology (handle)
 * @return
 *   - CEDAR_SUCCESS: no error
 */
int cedar_topo_create2d(MPI_Comm comm,
                        unsigned int ngx, unsigned int ngy,
                        unsigned int lnx[], unsigned int lny[],
                        int nprocx, int nprocy,
                        cedar_topo *newtopo);


/**
 * Create 3D distributed grid topology.
 *
 * @param[in]  comm MPI Communicator for distributed grid
 * @param[in]  ngx number of global grid points in first dimension
 * @param[in]  ngy number of global grid points in second dimension
 * @param[in]  ngz number of global grid points in third dimension
 * @param[in]  lnx array of length nprocx of local grid points for each processor
 * @param[in]  lny array of length nprocy of local grid points for each processor
 * @param[in]  lnz array of length nprocz of local grid points for each processor
 * @param[in]  nprocx number of processors in first dimension
 * @param[in]  nprocy number of processors in second dimension
 * @param[in]  nprocz number of processors in second dimension
 * @param[out] newtopo new topology (handle)
 * @return
 *   - CEDAR_SUCCESS: no error
 */
#ifdef ENABLE_3D
int cedar_topo_create3d(MPI_Comm comm,
                        unsigned int ngx, unsigned int ngy, unsigned int ngz,
                        unsigned int lnx[], unsigned int lny[], unsigned int lnz[],
                        int nprocx, int nprocy, int nprocz,
                        cedar_topo *newtopo);
#endif


/**
 * Free memory for distributed grid topology.
 *
 * @param topo topology to be destroyed
 * @return
 * - CEDAR_SUCCESS: no error
 * - CEDAR_ERR_TOPO: invalid topology object
 */
int cedar_topo_free(cedar_topo *topo);


/**
 * Create 2D vector (cedar grid function).
 *
 * @param[in] topo 2D distributed grid topology
 * @param[out] vec new vector (handle)
 * @return
 *   - CEDAR_SUCCESS: no error
 *   - CEDAR_ERR_TOPO: invalid topology object
 */
int cedar_vec_create2d(cedar_topo topo, cedar_vec *vec);


/**
 * Create 3D vector (cedar grid function).
 *
 * @param[in] topo 3D distributed grid topology
 * @param[out] vec new vector (handle)
 * @return
 *   - CEDAR_SUCCESS: no error
 *   - CEDAR_ERR_TOPO: invalid topology object
 */
#ifdef ENABLE_3D
int cedar_vec_create3d(cedar_topo topo, cedar_vec *vec);
#endif


/**
 * Free memory for cedar vector.
 *
 * @param vec vector to be destroyed
 * @return
 *   - CEDAR_SUCCESS: no error
 *   - CEDAR_ERR_VEC: invalid vector object
 */
int cedar_vec_free(cedar_vec *vec);

#include <cedar/2d/base_types.h>
typedef struct {
	unsigned int i;
	unsigned int j;
	bmg2_dir dir;
} cedar_coord_2d;
#ifdef ENABLE_3D
#include <cedar/3d/base_types.h>
typedef struct {
	unsigned int i;
	unsigned int j;
	unsigned int k;
	cdr3_dir dir;
} cedar_coord_3d;
#endif

typedef enum {
	CEDAR_STENCIL_FIVE_PT,
	CEDAR_STENCIL_NINE_PT
} cedar_stencil_2d;

typedef enum {
	CEDAR_STENCIL_SEVEN_PT,
	CEDAR_STENCIL_XXVII_PT
} cedar_stencil_3d;

/**
 * Create 2D cedar matrix.
 *
 * @param[in] topo 2D distributed grid topology
 * @param[out] mat new matrix
 * @return
 *   - CEDAR_SUCCESS: no error
 *   - CEDAR_ERR_TOPO: invalid topology object
 */
int cedar_mat_create2d(cedar_topo topo, cedar_stencil_2d sten_type, cedar_mat *mat);


/**
 * Create 3D cedar matrix.
 *
 * @param[in] topo 3D distributed grid topology
 * @param[out] mat new matrix
 * @return
 *   - CEDAR_SUCCESS: no error
 *   - CEDAR_ERR_TOPO: invalid topology object
 */
#ifdef ENABLE_3D
int cedar_mat_create3d(cedar_topo topo, cedar_stencil_3d sten_type, cedar_mat *mat);
#endif


/**
 * Set values in 2D cedar matrix.
 *
 * @param mat matrix to set values
 * @param nvals number of values to set
 * @param coords[nvals] array specifying stencil coordinates for values
 * @param vals[nvals] array of values to insert
 * @return
 *    - CEDAR_SUCCESS: no error
 *    - CEDAR_ERR_MAT: invalid matrix handle
 *    - CEDAR_ERR_DIM: invalid dimension
 *    - CEDAR_ERR_STENCIL: invalid stencil direction
 */
int cedar_mat_set2d(cedar_mat mat, unsigned int nvals, cedar_coord_2d coords[], cedar_real vals[]);


/**
 * Set values in 3D cedar matrix.
 *
 * @param mat matrix to set values
 * @param nvals number of values to set
 * @param coords[nvals] array specifying stencil coordinates for values
 * @param vals[nvals] array of values to insert
 * @return
 *    - CEDAR_SUCCESS: no error
 *    - CEDAR_ERR_MAT: invalid matrix handle
 *    - CEDAR_ERR_DIM: invalid dimension
 *    - CEDAR_ERR_STENCIL: invalid stencil direction
 */
#ifdef ENABLE_3D
int cedar_mat_set3d(cedar_mat mat, unsigned int nvals, cedar_coord_3d coords[], cedar_real vals[]);
#endif


/**
 * Computes y = mat * x
 */
int cedar_matvec(cedar_mat mat, cedar_vec x, cedar_vec y);


/**
 * Free memory for cedar matrix.
 *
 * @param mat matrix to be destroyed
 * @return
 *   - CEDAR_SUCCESS: no error
 *   - CEDAR_ERR_MAT: invalid mat object
 */
int cedar_mat_free(cedar_mat *mat);


/**
 * Create solver object and run cedar setup phase.
 *
 * @param mat fine grid matrix
 * @param solver solver object to create
 * @return
 *   - CEDAR_SUCCESS: no error
 *   - CEDAR_ERR_MAT: invalid mat handle
 */
int cedar_solver_create(cedar_mat mat, cedar_solver *solver);


/**
 * Run cedar solve.
 *
 * @param solver solver object
 * @param x solution vector
 * @param b right hand side
 * @return
 *   - CEDAR_SUCCESS: no error
 *   - CEDAR_ERR_VEC: invalid vector handle
 *   - CEDAR_ERR_SOLVER: invalid solver handle
 *   - CEDAR_ERR_DIM: dimension mismatch (between vectors and solver object)
 */
int cedar_solver_run(cedar_solver solver, cedar_vec x, cedar_vec b);


/**
* Free memory for cedar solver
*
* @param solver solver to destroy
*/
int cedar_solver_destroy(cedar_solver *solver);

#ifdef __cplusplus
}
#endif

#endif
