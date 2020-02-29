/**
 * @desc: Keccak 4-round pre-image attack implementation
 *  ╦╔═┌─┐┌─┐┌─┐┌─┐┬┌─  ╔═╗┌─┐┬ ┬  ┬┌─┐┬─┐
 *  ╠╩╗├┤ │  │  ├─┤├┴┐  ╚═╗│ ││ └┐┌┘├┤ ├┬┘
 *  ╩ ╩└─┘└─┘└─┘┴ ┴┴ ┴  ╚═╝└─┘┴─┘└┘ └─┘┴└─
 * @language: C++/CUDA
 */

#include <iostream>
#include <stdlib.h>
#include <inttypes.h>
#include <unistd.h>
#include <sys/time.h>
#include <sys/types.h>

#include <cuda_runtime.h>

#include "threadpool.h"

/* ******************************************************
 * *********** Params and Global Variables **************
 * ******************************************************/

#define __DEBUG__

// parameters relates to linear constraints
#define LIN_CONST_VARNUM 800        // number of variables in linear constraints
#define LIN_CONST_EQNUM 706
#define LIN_ITER_EQNUM 63
#define LIN_VART_EQNUM 12
#define LIN_CONST_SYSTEM_SIZE (801 * (LIN_CONST_EQNUM))
#define LIN_ITER_SYSTEM_SIZE (801 * (LIN_ITER_EQNUM + LIN_VART_EQNUM))

// parameters relates to MQ system
#define MQ_VAR_NUM 19
#define MQ_EQ_NUM 32
#define MQ_XVAR_NUM ((((MQ_VAR_NUM) * ((MQ_VAR_NUM) + 1)) / 2) + MQ_VAR_NUM + 1)
#define MQ_SYSTEM_SIZE ((MQ_EQ_NUM) * (MQ_XVAR_NUM))

// parameters relates to base MQ system
#define BMQ_VAR_NUM 94
#define BMQ_EQ_NUM 32
#define BMQ_XVAR_NUM ((((BMQ_VAR_NUM) * ((BMQ_VAR_NUM) + 1)) / 2) + BMQ_VAR_NUM + 1)

// parameters relates to Keccak analysis
#define ALPHA 0x44E72
#define BETA 0xBAA20F
#define DIFF_TOLERANCE 4
#define RC1 0x80000000
#define RC2 0x41010000
#define RC3 0x51010000
#define RC4 0x00010001

// parameters relates to GPU configuration
#define GPU_THREADS_PER_BLOCK 64
#define GPU_BLOCK_PER_STREAM 512
#define GPU_STREAM_NUM 8
#define GPU_TOTAL_THREADS ((GPU_STREAM_NUM) * (GPU_THREADS_PER_BLOCK) * (GPU_BLOCK_PER_STREAM))

// parameters relates to CPU threadpool
#define CPU_THREAD_NUM 20

#define GPU_MQGROUP_NUM 0x100000U       // 2 ** 24
#define GPU_MQGROUP_THREAD_SLICE ((GPU_MQGROUP_NUM) / (GPU_TOTAL_THREADS))

// use streams to parallel memory copy operation
cudaStream_t cudaStreams[GPU_STREAM_NUM];

/* ******************************************************
 * *************** Function Prototypes ******************
 * ******************************************************/

typedef struct KeccakSolver {

  bool verbose;                        /* decide whether print informative messages */

  /* linear constraints that will not alter during iteration */
  uint8_t const_constraints[LIN_CONST_EQNUM][LIN_CONST_VARNUM + 1];

  /* linear polynomials of which the constant part is being iterated*/
  uint8_t iterate_polynomials[LIN_ITER_EQNUM][LIN_CONST_VARNUM + 1];

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
    for (eq_idx = 0; eq_idx < LIN_CONST_EQNUM; eq_idx++) {
        fgets(line, sizeof(line), constFile);
        for (var_idx = 0; var_idx < LIN_CONST_VARNUM + 1; var_idx++) {
            keccakSolver->const_constraints[eq_idx][var_idx] = (uint8_t)(line[var_idx] - '0');
        }
    }

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

    for (eq_idx = 0; eq_idx < LIN_VART_EQNUM; eq_idx++) {
        fgets(line, sizeof(line), vartFile);
        for (var_idx = 0; var_idx < LIN_CONST_VARNUM + 1; var_idx++) {
            keccakSolver->iterate_polynomials[eq_idx][var_idx] = (uint8_t)(line[var_idx] - '0');
        }
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

/* function: ROR32
 * usage: perform right-rotate permutation
 */
__device__ __forceinline__
uint32_t
ROR32(uint32_t x, uint32_t a) {
    return (x >> a) | (x << (32 - a));
}

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

/* function: solveLinearSystem
 * usage: solve the given linear system using in-place Gaussian Jordan
 *        elimination
 * arguments:
 *        linsys: pointer to linear system
 *        solution: pointer to solution vector
 */
__device__ static __forceinline__
bool
solveLinearSystem(uint8_t linsys[LIN_CONST_VARNUM * (LIN_CONST_VARNUM + 1)], uint8_t *solution) {
    reducedRREF(linsys);
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
    for (eq_idx = 0; eq_idx < MQ_EQ_NUM; eq_idx++) {
        for (var_idx = 0; var_idx < MQ_VAR_NUM; var_idx++) {
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
kernelSolveMQSystems(uint8_t *mq_input_buffer, uint8_t *it_input_buffer,
                     uint8_t *const_buffer, uint32_t *output_buffer) {
    uint32_t tid = threadIdx.x + blockIdx.x * blockDim.x;
    uint32_t mq_input_buffer_offset = MQ_SYSTEM_SIZE * tid;
    uint32_t it_input_buffer_offset = LIN_ITER_SYSTEM_SIZE * tid;
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
    uint32_t tmpA[25];
    uint32_t C[5];
    uint32_t D[5];

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
            bool solvable = solveLinearSystem(aggr_linsys, aggr_linsys_sol);

            if (solvable) {
                // found a solution for initial status
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

                // Round 1 start here
                // Unrolled THETA operation
                C[0] = A[5];
                C[0] = C[0] ^ A[6];
                C[0] = C[0] ^ A[7];
                C[0] = C[0] ^ A[8];
                C[0] = C[0] ^ A[9];
                D[0] = ROR32(C[0], 1);

                C[1] = A[10];
                C[1] = C[1] ^ A[11];
                C[1] = C[1] ^ A[12];
                C[1] = C[1] ^ A[13];
                C[1] = C[1] ^ A[14];
                D[1] = ROR32(C[1], 1);

                C[2] = A[15];
                C[2] = C[2] ^ A[16];
                C[2] = C[2] ^ A[17];
                C[2] = C[2] ^ A[18];
                C[2] = C[2] ^ A[19];
                D[2] = ROR32(C[2], 1);

                C[3] = A[20];
                C[3] = C[3] ^ A[21];
                C[3] = C[3] ^ A[22];
                C[3] = C[3] ^ A[23];
                C[3] = C[3] ^ A[24];
                D[3] = ROR32(C[3], 1);

                C[4] = A[0];
                C[4] = C[4] ^ A[1];
                C[4] = C[4] ^ A[2];
                C[4] = C[4] ^ A[3];
                C[4] = C[4] ^ A[4];
                D[4] = ROR32(C[4], 1);

                A[0] = A[0] ^ C[3] ^ D[0];
                A[1] = A[1] ^ C[3] ^ D[0];
                A[2] = A[2] ^ C[3] ^ D[0];
                A[3] = A[3] ^ C[3] ^ D[0];
                A[4] = A[4] ^ C[3] ^ D[0];

                A[5] = A[5] ^ C[4] ^ D[1];
                A[6] = A[6] ^ C[4] ^ D[1];
                A[7] = A[7] ^ C[4] ^ D[1];
                A[8] = A[8] ^ C[4] ^ D[1];
                A[9] = A[9] ^ C[4] ^ D[1];

                A[10] = A[10] ^ C[0] ^ D[2];
                A[11] = A[11] ^ C[0] ^ D[2];
                A[12] = A[12] ^ C[0] ^ D[2];
                A[13] = A[13] ^ C[0] ^ D[2];
                A[14] = A[14] ^ C[0] ^ D[2];

                A[15] = A[15] ^ C[1] ^ D[3];
                A[16] = A[16] ^ C[1] ^ D[3];
                A[17] = A[17] ^ C[1] ^ D[3];
                A[18] = A[18] ^ C[1] ^ D[3];
                A[19] = A[19] ^ C[1] ^ D[3];

                A[20] = A[20] ^ C[2] ^ D[4];
                A[21] = A[21] ^ C[2] ^ D[4];
                A[22] = A[22] ^ C[2] ^ D[4];
                A[23] = A[23] ^ C[2] ^ D[4];
                A[24] = A[24] ^ C[2] ^ D[4];

                // Unrolled RHO + PI operation
                // tmpA[0] = A[0];
                tmpA[1] = A[1];
                tmpA[2] = A[2];
                tmpA[3] = A[3];
                tmpA[4] = A[4];
                tmpA[5] = A[5];
                tmpA[6] = A[6];
                tmpA[7] = A[7];
                tmpA[8] = A[8];
                tmpA[9] = A[9];
                tmpA[10] = A[10];
                tmpA[11] = A[11];
                tmpA[12] = A[12];
                tmpA[13] = A[13];
                tmpA[14] = A[14];
                tmpA[15] = A[15];
                tmpA[16] = A[16];
                tmpA[17] = A[17];
                tmpA[18] = A[18];
                tmpA[19] = A[19];
                tmpA[20] = A[20];
                tmpA[21] = A[21];
                tmpA[22] = A[22];
                tmpA[23] = A[23];
                tmpA[24] = A[24];

                // A[0] = tmpA[0];
                A[1] = ROR32(tmpA[15], 28);
                A[2] = ROR32(tmpA[5], 1);
                A[3] = ROR32(tmpA[20], 27);
                A[4] = ROR32(tmpA[10], 30);

                A[5] = ROR32(tmpA[6], 12);
                A[6] = ROR32(tmpA[21], 20);
                A[7] = ROR32(tmpA[11], 6);
                A[8] = ROR32(tmpA[1], 4);
                A[9] = ROR32(tmpA[16], 23);

                A[10] = ROR32(tmpA[12], 11);
                A[11] = ROR32(tmpA[2], 3);
                A[12] = ROR32(tmpA[17], 25);
                A[13] = ROR32(tmpA[7], 10);
                A[14] = ROR32(tmpA[22], 7);

                A[15] = ROR32(tmpA[18], 21);
                A[16] = ROR32(tmpA[8], 13);
                A[17] = ROR32(tmpA[23], 8);
                A[18] = ROR32(tmpA[13], 15);
                A[19] = ROR32(tmpA[3], 9);

                A[20] = ROR32(tmpA[24], 14);
                A[21] = ROR32(tmpA[14], 29);
                A[22] = ROR32(tmpA[4], 18);
                A[23] = ROR32(tmpA[19], 24);
                A[24] = ROR32(tmpA[9], 2);

                // Unrolled CHI operation
                C[0] = A[0] ^ ((~A[5]) & (A[10]));
                C[1] = A[5] ^ ((~A[10]) & (A[15]));
                C[2] = A[10] ^ ((~A[15]) & (A[20]));
                C[3] = A[15] ^ ((~A[20]) & (A[0]));
                C[4] = A[20] ^ ((~A[0]) & (A[5]));
                A[0] = C[0];
                A[5] = C[1];
                A[10] = C[2];
                A[15] = C[3];
                A[20] = C[4];

                C[0] = A[1] ^ ((~A[6]) & (A[11]));
                C[1] = A[6] ^ ((~A[11]) & (A[16]));
                C[2] = A[11] ^ ((~A[16]) & (A[21]));
                C[3] = A[16] ^ ((~A[21]) & (A[1]));
                C[4] = A[21] ^ ((~A[1]) & (A[6]));
                A[1] = C[0];
                A[6] = C[1];
                A[11] = C[2];
                A[16] = C[3];
                A[21] = C[4];

                C[0] = A[2] ^ ((~A[7]) & (A[12]));
                C[1] = A[7] ^ ((~A[12]) & (A[17]));
                C[2] = A[12] ^ ((~A[17]) & (A[22]));
                C[3] = A[17] ^ ((~A[22]) & (A[2]));
                C[4] = A[22] ^ ((~A[2]) & (A[7]));
                A[2] = C[0];
                A[7] = C[1];
                A[12] = C[2];
                A[17] = C[3];
                A[22] = C[4];

                C[0] = A[3] ^ ((~A[8]) & (A[13]));
                C[1] = A[8] ^ ((~A[13]) & (A[18]));
                C[2] = A[13] ^ ((~A[18]) & (A[23]));
                C[3] = A[18] ^ ((~A[23]) & (A[3]));
                C[4] = A[23] ^ ((~A[3]) & (A[8]));
                A[3] = C[0];
                A[8] = C[1];
                A[13] = C[2];
                A[18] = C[3];
                A[23] = C[4];

                C[0] = A[4] ^ ((~A[9]) & (A[14]));
                C[1] = A[9] ^ ((~A[14]) & (A[19]));
                C[2] = A[14] ^ ((~A[19]) & (A[24]));
                C[3] = A[19] ^ ((~A[24]) & (A[4]));
                C[4] = A[24] ^ ((~A[4]) & (A[9]));
                A[4] = C[0];
                A[9] = C[1];
                A[14] = C[2];
                A[19] = C[3];
                A[24] = C[4];

                A[0] = A[0] ^ RC1;
                // --- Round 1 end here ---

                // --- Round 2 start here ---
                // Unrolled THETA operation
                C[0] = A[5];
                C[0] = C[0] ^ A[6];
                C[0] = C[0] ^ A[7];
                C[0] = C[0] ^ A[8];
                C[0] = C[0] ^ A[9];
                D[0] = ROR32(C[0], 1);

                C[1] = A[10];
                C[1] = C[1] ^ A[11];
                C[1] = C[1] ^ A[12];
                C[1] = C[1] ^ A[13];
                C[1] = C[1] ^ A[14];
                D[1] = ROR32(C[1], 1);

                C[2] = A[15];
                C[2] = C[2] ^ A[16];
                C[2] = C[2] ^ A[17];
                C[2] = C[2] ^ A[18];
                C[2] = C[2] ^ A[19];
                D[2] = ROR32(C[2], 1);

                C[3] = A[20];
                C[3] = C[3] ^ A[21];
                C[3] = C[3] ^ A[22];
                C[3] = C[3] ^ A[23];
                C[3] = C[3] ^ A[24];
                D[3] = ROR32(C[3], 1);

                C[4] = A[0];
                C[4] = C[4] ^ A[1];
                C[4] = C[4] ^ A[2];
                C[4] = C[4] ^ A[3];
                C[4] = C[4] ^ A[4];
                D[4] = ROR32(C[4], 1);

                A[0] = A[0] ^ C[3] ^ D[0];
                A[1] = A[1] ^ C[3] ^ D[0];
                A[2] = A[2] ^ C[3] ^ D[0];
                A[3] = A[3] ^ C[3] ^ D[0];
                A[4] = A[4] ^ C[3] ^ D[0];

                A[5] = A[5] ^ C[4] ^ D[1];
                A[6] = A[6] ^ C[4] ^ D[1];
                A[7] = A[7] ^ C[4] ^ D[1];
                A[8] = A[8] ^ C[4] ^ D[1];
                A[9] = A[9] ^ C[4] ^ D[1];

                A[10] = A[10] ^ C[0] ^ D[2];
                A[11] = A[11] ^ C[0] ^ D[2];
                A[12] = A[12] ^ C[0] ^ D[2];
                A[13] = A[13] ^ C[0] ^ D[2];
                A[14] = A[14] ^ C[0] ^ D[2];

                A[15] = A[15] ^ C[1] ^ D[3];
                A[16] = A[16] ^ C[1] ^ D[3];
                A[17] = A[17] ^ C[1] ^ D[3];
                A[18] = A[18] ^ C[1] ^ D[3];
                A[19] = A[19] ^ C[1] ^ D[3];

                A[20] = A[20] ^ C[2] ^ D[4];
                A[21] = A[21] ^ C[2] ^ D[4];
                A[22] = A[22] ^ C[2] ^ D[4];
                A[23] = A[23] ^ C[2] ^ D[4];
                A[24] = A[24] ^ C[2] ^ D[4];

                // Unrolled RHO + PI operation
                // tmpA[0] = A[0];
                tmpA[1] = A[1];
                tmpA[2] = A[2];
                tmpA[3] = A[3];
                tmpA[4] = A[4];
                tmpA[5] = A[5];
                tmpA[6] = A[6];
                tmpA[7] = A[7];
                tmpA[8] = A[8];
                tmpA[9] = A[9];
                tmpA[10] = A[10];
                tmpA[11] = A[11];
                tmpA[12] = A[12];
                tmpA[13] = A[13];
                tmpA[14] = A[14];
                tmpA[15] = A[15];
                tmpA[16] = A[16];
                tmpA[17] = A[17];
                tmpA[18] = A[18];
                tmpA[19] = A[19];
                tmpA[20] = A[20];
                tmpA[21] = A[21];
                tmpA[22] = A[22];
                tmpA[23] = A[23];
                tmpA[24] = A[24];

                // A[0] = tmpA[0];
                A[1] = ROR32(tmpA[15], 28);
                A[2] = ROR32(tmpA[5], 1);
                A[3] = ROR32(tmpA[20], 27);
                A[4] = ROR32(tmpA[10], 30);

                A[5] = ROR32(tmpA[6], 12);
                A[6] = ROR32(tmpA[21], 20);
                A[7] = ROR32(tmpA[11], 6);
                A[8] = ROR32(tmpA[1], 4);
                A[9] = ROR32(tmpA[16], 23);

                A[10] = ROR32(tmpA[12], 11);
                A[11] = ROR32(tmpA[2], 3);
                A[12] = ROR32(tmpA[17], 25);
                A[13] = ROR32(tmpA[7], 10);
                A[14] = ROR32(tmpA[22], 7);

                A[15] = ROR32(tmpA[18], 21);
                A[16] = ROR32(tmpA[8], 13);
                A[17] = ROR32(tmpA[23], 8);
                A[18] = ROR32(tmpA[13], 15);
                A[19] = ROR32(tmpA[3], 9);

                A[20] = ROR32(tmpA[24], 14);
                A[21] = ROR32(tmpA[14], 29);
                A[22] = ROR32(tmpA[4], 18);
                A[23] = ROR32(tmpA[19], 24);
                A[24] = ROR32(tmpA[9], 2);

                // Unrolled CHI operation
                C[0] = A[0] ^ ((~A[5]) & (A[10]));
                C[1] = A[5] ^ ((~A[10]) & (A[15]));
                C[2] = A[10] ^ ((~A[15]) & (A[20]));
                C[3] = A[15] ^ ((~A[20]) & (A[0]));
                C[4] = A[20] ^ ((~A[0]) & (A[5]));
                A[0] = C[0];
                A[5] = C[1];
                A[10] = C[2];
                A[15] = C[3];
                A[20] = C[4];

                C[0] = A[1] ^ ((~A[6]) & (A[11]));
                C[1] = A[6] ^ ((~A[11]) & (A[16]));
                C[2] = A[11] ^ ((~A[16]) & (A[21]));
                C[3] = A[16] ^ ((~A[21]) & (A[1]));
                C[4] = A[21] ^ ((~A[1]) & (A[6]));
                A[1] = C[0];
                A[6] = C[1];
                A[11] = C[2];
                A[16] = C[3];
                A[21] = C[4];

                C[0] = A[2] ^ ((~A[7]) & (A[12]));
                C[1] = A[7] ^ ((~A[12]) & (A[17]));
                C[2] = A[12] ^ ((~A[17]) & (A[22]));
                C[3] = A[17] ^ ((~A[22]) & (A[2]));
                C[4] = A[22] ^ ((~A[2]) & (A[7]));
                A[2] = C[0];
                A[7] = C[1];
                A[12] = C[2];
                A[17] = C[3];
                A[22] = C[4];

                C[0] = A[3] ^ ((~A[8]) & (A[13]));
                C[1] = A[8] ^ ((~A[13]) & (A[18]));
                C[2] = A[13] ^ ((~A[18]) & (A[23]));
                C[3] = A[18] ^ ((~A[23]) & (A[3]));
                C[4] = A[23] ^ ((~A[3]) & (A[8]));
                A[3] = C[0];
                A[8] = C[1];
                A[13] = C[2];
                A[18] = C[3];
                A[23] = C[4];

                C[0] = A[4] ^ ((~A[9]) & (A[14]));
                C[1] = A[9] ^ ((~A[14]) & (A[19]));
                C[2] = A[14] ^ ((~A[19]) & (A[24]));
                C[3] = A[19] ^ ((~A[24]) & (A[4]));
                C[4] = A[24] ^ ((~A[4]) & (A[9]));
                A[4] = C[0];
                A[9] = C[1];
                A[14] = C[2];
                A[19] = C[3];
                A[24] = C[4];

                A[0] = A[0] ^ RC2;
                // --- Round 2 end here ---

                // --- Round 3 start here ---
                // Unrolled THETA operation
                C[0] = A[5];
                C[0] = C[0] ^ A[6];
                C[0] = C[0] ^ A[7];
                C[0] = C[0] ^ A[8];
                C[0] = C[0] ^ A[9];
                D[0] = ROR32(C[0], 1);

                C[1] = A[10];
                C[1] = C[1] ^ A[11];
                C[1] = C[1] ^ A[12];
                C[1] = C[1] ^ A[13];
                C[1] = C[1] ^ A[14];
                D[1] = ROR32(C[1], 1);

                C[2] = A[15];
                C[2] = C[2] ^ A[16];
                C[2] = C[2] ^ A[17];
                C[2] = C[2] ^ A[18];
                C[2] = C[2] ^ A[19];
                D[2] = ROR32(C[2], 1);

                C[3] = A[20];
                C[3] = C[3] ^ A[21];
                C[3] = C[3] ^ A[22];
                C[3] = C[3] ^ A[23];
                C[3] = C[3] ^ A[24];
                D[3] = ROR32(C[3], 1);

                C[4] = A[0];
                C[4] = C[4] ^ A[1];
                C[4] = C[4] ^ A[2];
                C[4] = C[4] ^ A[3];
                C[4] = C[4] ^ A[4];
                D[4] = ROR32(C[4], 1);

                A[0] = A[0] ^ C[3] ^ D[0];
                A[1] = A[1] ^ C[3] ^ D[0];
                A[2] = A[2] ^ C[3] ^ D[0];
                A[3] = A[3] ^ C[3] ^ D[0];
                A[4] = A[4] ^ C[3] ^ D[0];

                A[5] = A[5] ^ C[4] ^ D[1];
                A[6] = A[6] ^ C[4] ^ D[1];
                A[7] = A[7] ^ C[4] ^ D[1];
                A[8] = A[8] ^ C[4] ^ D[1];
                A[9] = A[9] ^ C[4] ^ D[1];

                A[10] = A[10] ^ C[0] ^ D[2];
                A[11] = A[11] ^ C[0] ^ D[2];
                A[12] = A[12] ^ C[0] ^ D[2];
                A[13] = A[13] ^ C[0] ^ D[2];
                A[14] = A[14] ^ C[0] ^ D[2];

                A[15] = A[15] ^ C[1] ^ D[3];
                A[16] = A[16] ^ C[1] ^ D[3];
                A[17] = A[17] ^ C[1] ^ D[3];
                A[18] = A[18] ^ C[1] ^ D[3];
                A[19] = A[19] ^ C[1] ^ D[3];

                A[20] = A[20] ^ C[2] ^ D[4];
                A[21] = A[21] ^ C[2] ^ D[4];
                A[22] = A[22] ^ C[2] ^ D[4];
                A[23] = A[23] ^ C[2] ^ D[4];
                A[24] = A[24] ^ C[2] ^ D[4];

                // Unrolled RHO + PI operation
                // tmpA[0] = A[0];
                tmpA[1] = A[1];
                tmpA[2] = A[2];
                tmpA[3] = A[3];
                tmpA[4] = A[4];
                tmpA[5] = A[5];
                tmpA[6] = A[6];
                tmpA[7] = A[7];
                tmpA[8] = A[8];
                tmpA[9] = A[9];
                tmpA[10] = A[10];
                tmpA[11] = A[11];
                tmpA[12] = A[12];
                tmpA[13] = A[13];
                tmpA[14] = A[14];
                tmpA[15] = A[15];
                tmpA[16] = A[16];
                tmpA[17] = A[17];
                tmpA[18] = A[18];
                tmpA[19] = A[19];
                tmpA[20] = A[20];
                tmpA[21] = A[21];
                tmpA[22] = A[22];
                tmpA[23] = A[23];
                tmpA[24] = A[24];

                // A[0] = tmpA[0];
                A[1] = ROR32(tmpA[15], 28);
                A[2] = ROR32(tmpA[5], 1);
                A[3] = ROR32(tmpA[20], 27);
                A[4] = ROR32(tmpA[10], 30);

                A[5] = ROR32(tmpA[6], 12);
                A[6] = ROR32(tmpA[21], 20);
                A[7] = ROR32(tmpA[11], 6);
                A[8] = ROR32(tmpA[1], 4);
                A[9] = ROR32(tmpA[16], 23);

                A[10] = ROR32(tmpA[12], 11);
                A[11] = ROR32(tmpA[2], 3);
                A[12] = ROR32(tmpA[17], 25);
                A[13] = ROR32(tmpA[7], 10);
                A[14] = ROR32(tmpA[22], 7);

                A[15] = ROR32(tmpA[18], 21);
                A[16] = ROR32(tmpA[8], 13);
                A[17] = ROR32(tmpA[23], 8);
                A[18] = ROR32(tmpA[13], 15);
                A[19] = ROR32(tmpA[3], 9);

                A[20] = ROR32(tmpA[24], 14);
                A[21] = ROR32(tmpA[14], 29);
                A[22] = ROR32(tmpA[4], 18);
                A[23] = ROR32(tmpA[19], 24);
                A[24] = ROR32(tmpA[9], 2);

                // Unrolled CHI operation
                C[0] = A[0] ^ ((~A[5]) & (A[10]));
                C[1] = A[5] ^ ((~A[10]) & (A[15]));
                C[2] = A[10] ^ ((~A[15]) & (A[20]));
                C[3] = A[15] ^ ((~A[20]) & (A[0]));
                C[4] = A[20] ^ ((~A[0]) & (A[5]));
                A[0] = C[0];
                A[5] = C[1];
                A[10] = C[2];
                A[15] = C[3];
                A[20] = C[4];

                C[0] = A[1] ^ ((~A[6]) & (A[11]));
                C[1] = A[6] ^ ((~A[11]) & (A[16]));
                C[2] = A[11] ^ ((~A[16]) & (A[21]));
                C[3] = A[16] ^ ((~A[21]) & (A[1]));
                C[4] = A[21] ^ ((~A[1]) & (A[6]));
                A[1] = C[0];
                A[6] = C[1];
                A[11] = C[2];
                A[16] = C[3];
                A[21] = C[4];

                C[0] = A[2] ^ ((~A[7]) & (A[12]));
                C[1] = A[7] ^ ((~A[12]) & (A[17]));
                C[2] = A[12] ^ ((~A[17]) & (A[22]));
                C[3] = A[17] ^ ((~A[22]) & (A[2]));
                C[4] = A[22] ^ ((~A[2]) & (A[7]));
                A[2] = C[0];
                A[7] = C[1];
                A[12] = C[2];
                A[17] = C[3];
                A[22] = C[4];

                C[0] = A[3] ^ ((~A[8]) & (A[13]));
                C[1] = A[8] ^ ((~A[13]) & (A[18]));
                C[2] = A[13] ^ ((~A[18]) & (A[23]));
                C[3] = A[18] ^ ((~A[23]) & (A[3]));
                C[4] = A[23] ^ ((~A[3]) & (A[8]));
                A[3] = C[0];
                A[8] = C[1];
                A[13] = C[2];
                A[18] = C[3];
                A[23] = C[4];

                C[0] = A[4] ^ ((~A[9]) & (A[14]));
                C[1] = A[9] ^ ((~A[14]) & (A[19]));
                C[2] = A[14] ^ ((~A[19]) & (A[24]));
                C[3] = A[19] ^ ((~A[24]) & (A[4]));
                C[4] = A[24] ^ ((~A[4]) & (A[9]));
                A[4] = C[0];
                A[9] = C[1];
                A[14] = C[2];
                A[19] = C[3];
                A[24] = C[4];

                A[0] = A[0] ^ RC3;
                // --- Round 3 end here ---

                // --- Round 4 start here ---
                // Unrolled THETA operation
                C[0] = A[5];
                C[0] = C[0] ^ A[6];
                C[0] = C[0] ^ A[7];
                C[0] = C[0] ^ A[8];
                C[0] = C[0] ^ A[9];
                D[0] = ROR32(C[0], 1);

                C[1] = A[10];
                C[1] = C[1] ^ A[11];
                C[1] = C[1] ^ A[12];
                C[1] = C[1] ^ A[13];
                C[1] = C[1] ^ A[14];
                D[1] = ROR32(C[1], 1);

                C[2] = A[15];
                C[2] = C[2] ^ A[16];
                C[2] = C[2] ^ A[17];
                C[2] = C[2] ^ A[18];
                C[2] = C[2] ^ A[19];
                D[2] = ROR32(C[2], 1);

                C[3] = A[20];
                C[3] = C[3] ^ A[21];
                C[3] = C[3] ^ A[22];
                C[3] = C[3] ^ A[23];
                C[3] = C[3] ^ A[24];
                D[3] = ROR32(C[3], 1);

                C[4] = A[0];
                C[4] = C[4] ^ A[1];
                C[4] = C[4] ^ A[2];
                C[4] = C[4] ^ A[3];
                C[4] = C[4] ^ A[4];
                D[4] = ROR32(C[4], 1);

                A[0] = A[0] ^ C[3] ^ D[0];
                A[1] = A[1] ^ C[3] ^ D[0];
                A[2] = A[2] ^ C[3] ^ D[0];
                A[3] = A[3] ^ C[3] ^ D[0];
                A[4] = A[4] ^ C[3] ^ D[0];

                A[5] = A[5] ^ C[4] ^ D[1];
                A[6] = A[6] ^ C[4] ^ D[1];
                A[7] = A[7] ^ C[4] ^ D[1];
                A[8] = A[8] ^ C[4] ^ D[1];
                A[9] = A[9] ^ C[4] ^ D[1];

                A[10] = A[10] ^ C[0] ^ D[2];
                A[11] = A[11] ^ C[0] ^ D[2];
                A[12] = A[12] ^ C[0] ^ D[2];
                A[13] = A[13] ^ C[0] ^ D[2];
                A[14] = A[14] ^ C[0] ^ D[2];

                A[15] = A[15] ^ C[1] ^ D[3];
                A[16] = A[16] ^ C[1] ^ D[3];
                A[17] = A[17] ^ C[1] ^ D[3];
                A[18] = A[18] ^ C[1] ^ D[3];
                A[19] = A[19] ^ C[1] ^ D[3];

                A[20] = A[20] ^ C[2] ^ D[4];
                A[21] = A[21] ^ C[2] ^ D[4];
                A[22] = A[22] ^ C[2] ^ D[4];
                A[23] = A[23] ^ C[2] ^ D[4];
                A[24] = A[24] ^ C[2] ^ D[4];

                // Unrolled RHO + PI operation
                // tmpA[0] = A[0];
                tmpA[1] = A[1];
                tmpA[2] = A[2];
                tmpA[3] = A[3];
                tmpA[4] = A[4];
                tmpA[5] = A[5];
                tmpA[6] = A[6];
                tmpA[7] = A[7];
                tmpA[8] = A[8];
                tmpA[9] = A[9];
                tmpA[10] = A[10];
                tmpA[11] = A[11];
                tmpA[12] = A[12];
                tmpA[13] = A[13];
                tmpA[14] = A[14];
                tmpA[15] = A[15];
                tmpA[16] = A[16];
                tmpA[17] = A[17];
                tmpA[18] = A[18];
                tmpA[19] = A[19];
                tmpA[20] = A[20];
                tmpA[21] = A[21];
                tmpA[22] = A[22];
                tmpA[23] = A[23];
                tmpA[24] = A[24];

                // A[0] = tmpA[0];
                A[1] = ROR32(tmpA[15], 28);
                A[2] = ROR32(tmpA[5], 1);
                A[3] = ROR32(tmpA[20], 27);
                A[4] = ROR32(tmpA[10], 30);

                A[5] = ROR32(tmpA[6], 12);
                A[6] = ROR32(tmpA[21], 20);
                A[7] = ROR32(tmpA[11], 6);
                A[8] = ROR32(tmpA[1], 4);
                A[9] = ROR32(tmpA[16], 23);

                A[10] = ROR32(tmpA[12], 11);
                A[11] = ROR32(tmpA[2], 3);
                A[12] = ROR32(tmpA[17], 25);
                A[13] = ROR32(tmpA[7], 10);
                A[14] = ROR32(tmpA[22], 7);

                A[15] = ROR32(tmpA[18], 21);
                A[16] = ROR32(tmpA[8], 13);
                A[17] = ROR32(tmpA[23], 8);
                A[18] = ROR32(tmpA[13], 15);
                A[19] = ROR32(tmpA[3], 9);

                A[20] = ROR32(tmpA[24], 14);
                A[21] = ROR32(tmpA[14], 29);
                A[22] = ROR32(tmpA[4], 18);
                A[23] = ROR32(tmpA[19], 24);
                A[24] = ROR32(tmpA[9], 2);

                // Unrolled CHI operation
                C[0] = A[0] ^ ((~A[5]) & (A[10]));
                C[1] = A[5] ^ ((~A[10]) & (A[15]));
                C[2] = A[10] ^ ((~A[15]) & (A[20]));
                C[3] = A[15] ^ ((~A[20]) & (A[0]));
                C[4] = A[20] ^ ((~A[0]) & (A[5]));
                A[0] = C[0];
                A[5] = C[1];
                A[10] = C[2];
                A[15] = C[3];
                A[20] = C[4];

                C[0] = A[1] ^ ((~A[6]) & (A[11]));
                C[1] = A[6] ^ ((~A[11]) & (A[16]));
                C[2] = A[11] ^ ((~A[16]) & (A[21]));
                C[3] = A[16] ^ ((~A[21]) & (A[1]));
                C[4] = A[21] ^ ((~A[1]) & (A[6]));
                A[1] = C[0];
                A[6] = C[1];
                A[11] = C[2];
                A[16] = C[3];
                A[21] = C[4];

                C[0] = A[2] ^ ((~A[7]) & (A[12]));
                C[1] = A[7] ^ ((~A[12]) & (A[17]));
                C[2] = A[12] ^ ((~A[17]) & (A[22]));
                C[3] = A[17] ^ ((~A[22]) & (A[2]));
                C[4] = A[22] ^ ((~A[2]) & (A[7]));
                A[2] = C[0];
                A[7] = C[1];
                A[12] = C[2];
                A[17] = C[3];
                A[22] = C[4];

                C[0] = A[3] ^ ((~A[8]) & (A[13]));
                C[1] = A[8] ^ ((~A[13]) & (A[18]));
                C[2] = A[13] ^ ((~A[18]) & (A[23]));
                C[3] = A[18] ^ ((~A[23]) & (A[3]));
                C[4] = A[23] ^ ((~A[3]) & (A[8]));
                A[3] = C[0];
                A[8] = C[1];
                A[13] = C[2];
                A[18] = C[3];
                A[23] = C[4];

                C[0] = A[4] ^ ((~A[9]) & (A[14]));
                C[1] = A[9] ^ ((~A[14]) & (A[19]));
                C[2] = A[14] ^ ((~A[19]) & (A[24]));
                C[3] = A[19] ^ ((~A[24]) & (A[4]));
                C[4] = A[24] ^ ((~A[4]) & (A[9]));
                A[4] = C[0];
                A[9] = C[1];
                A[14] = C[2];
                A[19] = C[3];
                A[24] = C[4];

                A[0] = A[0] ^ RC4;
                // --- Round 4 end here ---

                // check correctness of hash
                uint32_t hash0 = A[0];
                uint32_t hash1 = A[5];
                uint32_t hash2 = A[10] & 0xFFFF0000;
                uint32_t bit_diff = 0;
                bit_diff += __popc(hash0 ^ 0x751A16E5);
                bit_diff += __popc(hash1 ^ 0xE495E1E2);
                bit_diff += __popc(hash2 ^ 0xFF220000);

                if (bit_diff <= DIFF_TOLERANCE) {
                    // record the initial status to output buffer
                    memcpy(output_buffer + output_buffer_offset, initA, sizeof(initA));
                }
            }

        }
    }
}

typedef struct threadArgs {
  KeccakSolver *keccakSolver;
  uint64_t thread_gb;
} threadArgs;

/* function: threadUpdateMQSystems
 * usage: cpu threading function to update GPU_MQGROUP_NUM mq systems.
 *        the input arguments give the guessing bits in iteration, given of which this
 *        function can fill out the variate linear constraints and update the mq system.
 * arguments:
 *        threadArgs-keccakSolver: pointer to Keccak Solver
 *        threadArgs-thread_gb: guessing bits
 */
void
threadUpdateMQSystems(void *args) {

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

    threadpool_t *threadpool = threadpool_create(CPU_THREAD_NUM, GPU_MQGROUP_NUM, 0);
    uint64_t guessingBits;
    threadArgs thread_gb[GPU_MQGROUP_NUM];

    dim3 threads_per_block(GPU_THREADS_PER_BLOCK);
    const uint32_t device_mq_buffer_ssize = (GPU_MQGROUP_NUM * MQ_SYSTEM_SIZE) / GPU_STREAM_NUM;
    const uint32_t device_it_buffer_ssize = (GPU_MQGROUP_NUM * LIN_ITER_SYSTEM_SIZE) / GPU_STREAM_NUM;

    // main loop
    for (guessingBits = start; guessingBits < end; guessingBits += GPU_MQGROUP_NUM) {

        // use thread pool to update mq systems
        uint32_t taskId;
        for (taskId = 0; taskId < GPU_MQGROUP_NUM; taskId++) {
            thread_gb[taskId].keccakSolver = keccakSolver;
            thread_gb[taskId].thread_gb = guessingBits + taskId;
            threadpool_add(threadpool, threadUpdateMQSystems, (void *) &thread_gb[taskId], 0);
        }
        threadpool_join(threadpool, 0);

        // wait gpu to finish the computation
        CUDA_CHECK(cudaDeviceSynchronize());

        // launch kernel to solve MQ systems
        uint32_t streamId;
        uint32_t mq_buffer_soffset = 0;
        uint32_t it_buffer_soffset = 0;
        for (streamId = 0; streamId < GPU_STREAM_NUM; streamId++) {
            mq_buffer_soffset = device_mq_buffer_ssize * streamId;
            it_buffer_soffset = device_it_buffer_ssize * streamId;

            // copy mq system to device
            CUDA_CHECK(cudaMemcpyAsync(keccakSolver->device_mq_input_buffer + mq_buffer_soffset,
                                       keccakSolver->host_mq_input_buffer + mq_buffer_soffset,
                                       device_mq_buffer_ssize,
                                       cudaMemcpyHostToDevice,
                                       cudaStreams[streamId]));
            // copy iterative constraints to device
            CUDA_CHECK(cudaMemcpyAsync(keccakSolver->device_it_input_buffer + it_buffer_soffset,
                                       keccakSolver->host_it_input_buffer + it_buffer_soffset,
                                       device_it_buffer_ssize,
                                       cudaMemcpyHostToDevice,
                                       cudaStreams[streamId]));

            kernelSolveMQSystems << < GPU_BLOCK_PER_STREAM, threads_per_block, 0, cudaStreams[streamId] >> >
                (keccakSolver->device_mq_input_buffer,
                    keccakSolver->device_it_input_buffer,
                    keccakSolver->device_constant_buffer,
                    keccakSolver->device_output_buffer);
        }
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
