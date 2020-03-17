/**
 * @desc: Keccak 4-round pre-image attack implementation
 *  ╦╔═┌─┐┌─┐┌─┐┌─┐┬┌─  ╔═╗┌─┐┬ ┬  ┬┌─┐┬─┐
 *  ╠╩╗├┤ │  │  ├─┤├┴┐  ╚═╗│ ││ └┐┌┘├┤ ├┬┘
 *  ╩ ╩└─┘└─┘└─┘┴ ┴┴ ┴  ╚═╝└─┘┴─┘└┘ └─┘┴└─
 * @language: C++/CUDA
 */

#include <cuda_runtime.h>

#include <iostream>
#include <stdlib.h>
#include <inttypes.h>
#include <unistd.h>
#include <sys/time.h>
#include <sys/types.h>

#include "keccak.h"

/* ******************************************************
 * *********** Params and Global Variables **************
 * ******************************************************/

#define __DEBUG__

// parameters relates to linear constraints
#define LIN_CONST_VARNUM 800        // number of variables in linear constraints
#define LIN_CONST_EQNUM 707
#define LIN_ITER_EQNUM 63
#define LIN_CONST_SYSTEM_SIZE (801 * (LIN_CONST_EQNUM))
#define LIN_ITER_SYSTEM_SIZE (801 * (LIN_ITER_EQNUM))

// parameters relates to MQ system
#define MQ_VAR_NUM 31
#define MQ_EQ_NUM 48
#define MQ_XVAR_NUM ((((MQ_VAR_NUM) * ((MQ_VAR_NUM) + 1)) / 2) + MQ_VAR_NUM + 1)
#define MQ_SYSTEM_SIZE ((MQ_EQ_NUM) * (MQ_XVAR_NUM))

// parameters relates to base MQ system
#define BMQ_VAR_NUM 94
#define BMQ_EQ_NUM 48
#define BMQ_XVAR_NUM ((((BMQ_VAR_NUM) * ((BMQ_VAR_NUM) + 1)) / 2) + BMQ_VAR_NUM + 1)

// parameters relates to GPU configuration
#define GPU_THREADS_PER_BLOCK 64
#define GPU_BLOCK_PER_STREAM 512
#define GPU_TOTAL_THREADS ((GPU_STREAM_NUM) * (GPU_THREADS_PER_BLOCK) * (GPU_BLOCK_PER_STREAM))

/* ******************************************************
 * *************** Function Prototypes ******************
 * ******************************************************/
typedef struct KeccakSolver {

  bool verbose;                        /* decide whether to print out informative messages */

  /* iterative term to linearize last round terms */
  uint8_t iterate_terms[LIN_VART_EQNUM][LIN_ITER_EQNUM][LIN_CONST_VARNUM + 1];

  /* base mq system */
  uint8_t mq_polynomials[BMQ_EQ_NUM][BMQ_XVAR_NUM];

  uint32_t dev_id;                     /* specify which gpu device to use */

  uint8_t *host_mq_input_buffer;
  uint8_t *host_it_input_buffer;

  uint32_t *device_output_buffer;      /* if gpu threads get a result, it will be copied to output buffer */
  uint8_t *device_mq_input_buffer;     /* gpu buffer to store input mq systems */
  uint8_t *device_it_input_buffer;     /* gpu buffer to store input iterative constraints */
  uint8_t *device_constant_buffer;     /* gpu buffer to store linear constraints that will not alter during iteration
                                         * once gpu obtain a solution from the mq system, it will continue to solve the
                                         * remaining linear system and check whether the resulting hash is the goal */
} KeccakSolver;

void keccakSolverInit(KeccakSolver *keccakSolver);
void keccakSolverLoop(KeccakSolver *keccakSolver, uint64_t start, uint64_t end);
void keccakSolverExit(KeccakSolver *keccakSolver);

int main() {
    KeccakSolver keccakSolver;
    keccakSolver.dev_id = 0;
    keccakSolver.verbose = true;

    keccakSolverInit(&keccakSolver);
    keccakSolverLoop(&keccakSolver, 0, 0xFFFFF);
    keccakSolverExit(&keccakSolver);

    return 0;
}

/* ******************************************************
 * ************* Function Implementations ***************
 * ******************************************************/

__constance__ uint32_t file_mq_fvar_idx[BMQ_VAR_NUM] = {0};
__constance__ uint32_t file_reverse_mq_fvar_idx[LIN_CONST_VARNUM] = {0};
__constance__ uint32_t mq_fvar_idx[MQ_VAR_NUM] = {0};
__constance__ uint32_t reverse_mq_fvar_idx[MQ_VAR_NUM] = {0};

__host__
void
cuda_check(cudaError_t err, const char *file, const int line, bool fatal) {
    if (cudaSuccess != err) {
        fprintf(stderr, "[!] GPU error: %s:%d, code: %d, reason: %s\n",
                file, line, err, cudaGetErrorString(err));

        if (fatal) {
            fprintf(stderr, "[!] aborting...\n");
            exit(err);
        }
    }
}

#define CUDA_CHECK(call) do { \
    cuda_check((call), __FILE__, __LINE__, 1); \
} while(0)

/* function: loadConstantSystemToDevice
 * usage: read pre-processed const linear constraints and store them
 *      on gpu memory before iteration
 * arguments:
 *      keccakSolver: pointer to keccak solver
 */
void
loadConstantSystemToDevice(KeccakSolver *keccakSolver) {
    FILE *constFile = fopen("const.dat", "r");
    char line[LIN_CONST_VARNUM + 3];

    if (constFile == nullptr) {
        fprintf(stderr, "[!] error: cannot read const constraints file\n");
        fprintf(stderr, "[!] aborting...\n");
        exit(1);
    }

    uint32_t eq_idx, var_idx;
    if (keccakSolver->verbose) {
        printf("[+] loaded constant constraints constraints from file\n");
    }

    fclose(constFile);
}

/* function: loadConstantSystemToDevice
 * usage: read pre-processed iterative linear constraints and store them
 *      in keccak solver before iteration
 * arguments:
 *      keccakSolver: pointer to keccak solver
 */
void
loadVariateSystemToHost(KeccakSolver *keccakSolver) {
    uint32_t eq_idx, var_idx, term_idx;

    FILE *vartFile = fopen("vart.dat", "r");
    char line[LIN_CONST_VARNUM + 3];

    if (vartFile == nullptr) {
        fprintf(stderr, "[!] error: cannot read variate constraints file\n");
        fprintf(stderr, "[!] aborting...\n");
        exit(1);
    }

    if (keccakSolver->verbose) {
        printf("[+] loaded iterative constraints from file\n");
    }

    fclose(vartFile);

    FILE *vartTermFile = fopen("vart_terms.dat", "r");
    fclose(vartTermFile);

}

/* function: ctz
 * usage: return the number of trailing zeros in given number
 */
__device__ static __forceinline__
uint32_t
ctz(uint32_t c) {
    return __clz(__brev(c));
}

/* function: reducedRREF(linsys)
 * usage: perform Gaussian elimination to obtain an equivalent
 *       linear system in reduced row echelon form
 */
__device__ static __forceinline__
void
reducedRREF(uint8_t *linsys) {
    uint64_t i, j, k;
    const uint64_t col_num = LIN_CONST_VARNUM + 1;
    const uint64_t row_num = LIN_CONST_VARNUM;

    uint8_t tmp_row[col_num];
    memset(tmp_row, 0, sizeof(tmp_row));

    for (i = 0; i < LIN_CONST_VARNUM; i++) {
        // find the i-th pivot
        for (j = i; j < row_num; j++) {
            if (linsys[j * col_num + i]) {
                break;
            }
        }

        if (row_num == j) {
            continue;
        }

        memcpy(tmp_row, linsys + i * col_num, col_num);
        memcpy(linsys + i * col_num, linsys + j * col_num, col_num);
        memcpy(linsys + j * col_num, tmp_row, col_num);

        // for all the rows below pivot
        for (j = i + 1; j < row_num; ++j) {
            if (linsys[j * col_num + i]) { // sys[j][i]
                // subtract i-th row from the row
                for (k = 0; k < col_num; ++k) {
                    linsys[j * col_num + k] ^= linsys[i * col_num + k];
                }
            }
        }

        // for all the rows above pivot
        for (j = 0; j < i; ++j) {
            if (linsys[j * col_num + i]) { // sys[j][i]
                // subtract i-th row from the row
                for (k = 0; k < col_num; ++k) {
                    linsys[j * col_num + k] ^= linsys[i * col_num + k];
                }
            }
        }
    }
}

/* function: extractSolution(linsys, solution)
 * usage: extract a solution from a linear system in rref
 */
__device__ static __forceinline__
void
extractSolution(uint8_t *linsys, uint8_t *solution) {
    uint64_t eq_idx, var_idx;
    const uint64_t col_num = LIN_CONST_VARNUM + 1;
    const uint64_t row_num = LIN_CONST_VARNUM;

    for (var_idx = 0; var_idx < LIN_CONST_VARNUM; var_idx++) {
        for (eq_idx = var_idx; eq_idx < row_num; eq_idx++) {
            if (linsys[eq_idx * col_num + var_idx]) {
                break;
            }
        }
        if (row_num == eq_idx) {
            continue;
        }
        solution[var_idx] = linsys[eq_idx * col_num + col_num - 1];
    }

}

/* function: solveLinearSystem
 * usage: solve the given linear system using in-place Gaussian Jordan
 *        elimination
 * arguments:
 *        linsys: pointer to linear system
 *        solution: pointer to solution vector
 */
__device__ static __forceinline__
void
solveLinearSystem(uint8_t linsys[LIN_CONST_VARNUM * (LIN_CONST_VARNUM + 1)],
                  uint8_t *solution) {
    reducedRREF(linsys);
    extractSolution(linsys, solution);
}

__device__ static __forceinline__
void
diffEq(uint8_t *func, uint32_t idx, bool *result) {
    for (uint32_t i = 0; i < MQ_VAR_NUM + 1; i++)
        result[i] = 0x0;
    result[MQ_VAR_NUM] = func[MQ_XVAR_NUM - 1 - (MQ_VAR_NUM - idx)];
    uint32_t cursor = MQ_XVAR_NUM - MQ_VAR_NUM - 1;
    uint32_t i, j;
    uint32_t bound = (0 == idx) ? 1 : idx;
    for (i = MQ_VAR_NUM - 1; i >= bound; i--) {
        if (i == idx) {
            for (j = 1; j < i; j++) {
                result[i - j] ^= func[cursor - j];
            }
        } else {
            result[i] ^= func[cursor - (i - idx)];
        }
        cursor -= i + 1;
    }
}

/* function: findPartialDerivs
 * usage: find the partial derivatives of given mq system
 */
__device__ static __forceinline__
void
findPartialDerivs(uint8_t *mqs, bool derivs[MQ_EQ_NUM][MQ_VAR_NUM][MQ_VAR_NUM + 1]) {
    uint32_t eq_idx, var_idx;
    for (
        eq_idx = 0;
        eq_idx < MQ_EQ_NUM; eq_idx++) {
        for (
            var_idx = 0;
            var_idx < MQ_VAR_NUM; var_idx++) {
            diffEq(mqs + eq_idx * MQ_XVAR_NUM, var_idx, derivs[eq_idx][var_idx]);
        }
    }
}

/* function: kernelSolveMQSystems
 * usage: gpu kernel function to solve GPU_MQGROUP_NUM mq systems.
 *       in this implementation we apply Fast Exhaustive Search to find the
 *       solution of MQ system. once one thread finds a solution, it will continue
 *       to solve the remaining linear system and perform Keccak 4-round hash to
 *       test whether the result is a legal goal.
 * arguments:
 *      mq_input_buffer: pointer of buffer to store mq systems
 *      it_input_buffer: pointer of buffer to store linear constraints involved in iteration
 *      const_buffer: pointer of buffer to store constance constraints
 *      output_buffer: pointer of buffer to store potential preimages
 */
__global__
void
kernelIteration(uint64_t global_guessingbits_start,
                uint64_t global_guessingbits_end,
                uint8_t *mq_input_buffer, uint8_t *it_input_buffer,
                uint8_t *const_buffer, uint32_t *output_buffer) {
    uint32_t tid = threadIdx.x + blockIdx.x * blockDim.x;
    uint32_t output_buffer_offset = 25 * tid;

    bool derivs[MQ_EQ_NUM][MQ_VAR_NUM][MQ_VAR_NUM + 1];

    uint32_t pdiff2[MQ_EQ_NUM][MQ_VAR_NUM];
    uint32_t pdiff_eval[MQ_VAR_NUM];
    uint32_t mq_idx, var_idx, eq_idx, i, term;
    uint8_t *this_mq;
    uint8_t *this_it;
    uint8_t mq_solution[MQ_VAR_NUM][LIN_CONST_VARNUM + 1];
    uint8_t aggr_linsys[LIN_CONST_VARNUM * (LIN_CONST_VARNUM + 1)];
    uint8_t aggr_linsys_sol[LIN_CONST_VARNUM];

    uint32_t A[25];
    uint32_t initA[25];

    // fill coefficient in mq_solution
    memset(mq_solution, 0, sizeof(mq_solution));
    for (var_idx = 0; var_idx < MQ_VAR_NUM; var_idx++) {
        mq_solution[var_idx][780 + var_idx] = 1;
    }

    // main loop
    for (mq_idx = 0; mq_idx < GPU_MQGROUP_THREAD_SLICE; mq_idx++) {
        memset(pdiff2, 0, sizeof(pdiff2));
        this_mq = mq_input_buffer + mq_input_buffer_offset + mq_idx * MQ_SYSTEM_SIZE;
        this_it = it_input_buffer + it_input_buffer_offset + mq_idx * LIN_ITER_SYSTEM_SIZE;

        // set partial derivatives
        findPartialDerivs(this_mq, derivs);
        for (var_idx = 0; var_idx < MQ_VAR_NUM; var_idx++) {
            for (i = 0; i < MQ_VAR_NUM; i++) {
                for (eq_idx = 0; eq_idx < MQ_EQ_NUM; eq_idx++) {
                    term = (uint32_t) derivs[eq_idx][var_idx][i];
                    pdiff2[i][var_idx] |= term << eq_idx;
                }
            }
        }

        for (var_idx = 0; var_idx < MQ_VAR_NUM; var_idx++) {
            for (eq_idx = 0; eq_idx < MQ_EQ_NUM; eq_idx++) {
                if (0 == var_idx) {
                    term = (uint32_t) derivs[eq_idx][0][MQ_VAR_NUM];
                } else {
                    term = (uint32_t) derivs[eq_idx][var_idx][MQ_VAR_NUM] ^
                        derivs[eq_idx][var_idx][var_idx - 1];
                }
                pdiff_eval[var_idx] |= term << eq_idx;
            }
        }

        uint32_t func_eval = 0x0U;
        for (eq_idx = 0; eq_idx < MQ_EQ_NUM; eq_idx++) {
            term = (uint32_t) this_mq[eq_idx * MQ_XVAR_NUM + MQ_XVAR_NUM - 1];
            func_eval |= term << eq_idx;
        }

        // exhaustive search to find a solution fo mq system
        uint32_t count = 0;
        const uint32_t bound = (1U << MQ_VAR_NUM) - 1;
        uint32_t fp_idx, pre_fp_idx;
        while (func_eval && count < bound) {
            count++;
            fp_idx = ctz(count);

            if (count & (count - 1)) {
                pre_fp_idx = ctz(count ^ (1U << fp_idx));
                pdiff_eval[fp_idx] ^= pdiff2[fp_idx][pre_fp_idx];
            }

            func_eval ^= pdiff_eval[fp_idx];
        }

        if (!func_eval) {
            // found a solution for mq system
            for (var_idx = 0; var_idx < MQ_VAR_NUM; var_idx++) {
                mq_solution[var_idx][LIN_CONST_VARNUM] = ((count ^ count >> 1) >> var_idx) & 0x1U;
            }

            // solve the aggregate linear system
            memset(aggr_linsys, 0, sizeof(aggr_linsys));
            memset(aggr_linsys_sol, 0, sizeof(aggr_linsys_sol));

            memcpy(aggr_linsys, const_buffer, LIN_CONST_SYSTEM_SIZE * sizeof(uint8_t));
            memcpy(aggr_linsys + LIN_CONST_SYSTEM_SIZE, this_it, LIN_ITER_SYSTEM_SIZE * sizeof(uint8_t));
            memcpy(aggr_linsys + LIN_CONST_SYSTEM_SIZE + LIN_ITER_SYSTEM_SIZE, mq_solution, sizeof(mq_solution));
            solveLinearSystem(aggr_linsys, aggr_linsys_sol);

            memset(A, 0, sizeof(A));

            // set initial status
            for (i = 0; i < 32; i++) {
                A[0] |= (aggr_linsys_sol[i] << i);
                A[1] |= (aggr_linsys_sol[i + 32] << i);
                A[2] |= (aggr_linsys_sol[i + 64] << i);
                A[3] |= (aggr_linsys_sol[i + 96] << i);
                A[4] |= (aggr_linsys_sol[i + 128] << i);
                A[5] |= (aggr_linsys_sol[i + 160] << i);
                A[6] |= (aggr_linsys_sol[i + 192] << i);
                A[7] |= (aggr_linsys_sol[i + 224] << i);
                A[8] |= (aggr_linsys_sol[i + 256] << i);
                A[9] |= (aggr_linsys_sol[i + 288] << i);
                A[10] |= (aggr_linsys_sol[i + 320] << i);
                A[11] |= (aggr_linsys_sol[i + 352] << i);
                A[12] |= (aggr_linsys_sol[i + 384] << i);
                A[13] |= (aggr_linsys_sol[i + 416] << i);
                A[14] |= (aggr_linsys_sol[i + 448] << i);
                A[15] |= (aggr_linsys_sol[i + 480] << i);
                A[16] |= (aggr_linsys_sol[i + 512] << i);
                A[17] |= (aggr_linsys_sol[i + 544] << i);
                A[18] |= (aggr_linsys_sol[i + 576] << i);
                A[19] |= (aggr_linsys_sol[i + 608] << i);
                A[20] |= (aggr_linsys_sol[i + 640] << i);
                A[21] |= (aggr_linsys_sol[i + 672] << i);
                A[22] |= (aggr_linsys_sol[i + 704] << i);
                A[23] |= (aggr_linsys_sol[i + 736] << i);
                A[24] |= (aggr_linsys_sol[i + 768] << i);
            }
            memcpy(initA, A, sizeof(initA));

            if (verifyKeccakResult(A)) {
                // record the initial status to output buffer
                memcpy(output_buffer + output_buffer_offset, initA, sizeof(initA));
            }
        }
    }
}

/* function: keccakSolverInit
 * usage: initialize Keccak Solver. operations involve:
 *      memory allocation, load pre-processed files, etc.
 */
void
keccakSolverInit(KeccakSolver *keccakSolver) {
    // initialize gpu device
    CUDA_CHECK(cudaSetDevice(keccakSolver->dev_id));
    CUDA_CHECK(cudaDeviceSynchronize());
    CUDA_CHECK(cudaDeviceReset());

    // create cudaStreams
    for (uint32_t i = 0; i < GPU_STREAM_NUM; i++)
        CUDA_CHECK(cudaStreamCreate(&cudaStreams[i]));

    // allocate memory for device
    uint64_t kern_mq_input_malloc_size = GPU_MQGROUP_NUM * MQ_SYSTEM_SIZE * sizeof(uint8_t);
    uint64_t kern_it_input_malloc_size = GPU_MQGROUP_NUM * LIN_ITER_SYSTEM_SIZE * sizeof(uint8_t);
    uint64_t kern_const_malloc_size = LIN_CONST_SYSTEM_SIZE * sizeof(uint8_t);
    uint64_t kern_output_malloc_size = 25 * GPU_TOTAL_THREADS * sizeof(uint32_t);
    CUDA_CHECK(cudaMalloc(&keccakSolver->device_mq_input_buffer, kern_mq_input_malloc_size));
    CUDA_CHECK(cudaMalloc(&keccakSolver->device_it_input_buffer, kern_it_input_malloc_size));
    CUDA_CHECK(cudaMalloc(&keccakSolver->device_constant_buffer, kern_const_malloc_size));
    CUDA_CHECK(cudaMalloc(&keccakSolver->device_output_buffer, kern_output_malloc_size));

    // allocate memory for host
    CUDA_CHECK(cudaMallocHost(&keccakSolver->host_mq_input_buffer, kern_mq_input_malloc_size));
    CUDA_CHECK(cudaMallocHost(&keccakSolver->host_it_input_buffer, kern_it_input_malloc_size));

    loadConstantSystemToDevice(keccakSolver);
}

/* function: keccakSolverLoop
 * usage: main loop of Keccak Solver. iterating guessing bits from start to end
 */
void
keccakSolverLoop(KeccakSolver *keccakSolver, uint64_t start, uint64_t end) {
    dim3 threads_per_block(GPU_THREADS_PER_BLOCK);
    const uint32_t device_mq_buffer_ssize = (GPU_MQGROUP_NUM * MQ_SYSTEM_SIZE) / GPU_STREAM_NUM;
    const uint32_t device_it_buffer_ssize = (GPU_MQGROUP_NUM * LIN_ITER_SYSTEM_SIZE) / GPU_STREAM_NUM;

    // main loop
    kernelSolveMQSystems << <>>>;
    for (guessingBits = start; guessingBits < end; guessingBits += GPU_MQGROUP_NUM) {
        // wait gpu to finish the computation
        CUDA_CHECK(cudaDeviceSynchronize());
    }

    CUDA_CHECK(cudaDeviceSynchronize());
}

void
keccakSolverExit(KeccakSolver *keccakSolver) {
    // release registered memory
    CUDA_CHECK(cudaFree(&keccakSolver->device_mq_input_buffer));
    CUDA_CHECK(cudaFree(&keccakSolver->device_it_input_buffer));
    CUDA_CHECK(cudaFree(&keccakSolver->device_constant_buffer));
    CUDA_CHECK(cudaFree(&keccakSolver->device_output_buffer));
    CUDA_CHECK(cudaFree(&keccakSolver->host_mq_input_buffer));
    CUDA_CHECK(cudaFree(&keccakSolver->host_it_input_buffer));
}
