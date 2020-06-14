//
// Created by vitowu on 3/13/20.
//

#include "solver.h"

#include "threadpool.h"
#include "params.h"
#include "keccak.h"

__device__ static uint32_t __forceinline__
ctz(uint32_t c) {
    return __clz(__brev(c));
}

__device__ static void __forceinline__
diff_eq(uint8_t *func, uint32_t idx, uint8_t *result) {
    memset(result, 0x0, (MQ_VAR_NUM + 1) * sizeof(uint8_t));

    result[MQ_VAR_NUM] = func[MQ_XVAR_NUM - 1 - (MQ_VAR_NUM - idx)];

    uint32_t cursor = MQ_XVAR_NUM - 1 - MQ_VAR_NUM - 1;
    uint32_t i, j;
    uint32_t bound = (0 == idx) ? 1 : idx;
    for (i = MQ_VAR_NUM - 1; i >= bound; --i) {
        if (i == idx) {
            for (j = 1; j <= i; ++j) {
                result[i - j] ^= func[cursor - j];
            }
        } else {
            result[i] ^= func[cursor - (i - idx)];
        }
        cursor -= i + 1;
    }
}

__device__ static void __forceinline__
find_partial_derivs(uint8_t *mqsystem, uint8_t derivs[MQ_EQ_NUM][MQ_VAR_NUM][MQ_VAR_NUM + 1]) {
    uint32_t eq_idx, var_idx;
    for (eq_idx = 0; eq_idx < MQ_EQ_NUM; eq_idx++) {
        for (var_idx = 0; var_idx < MQ_VAR_NUM; var_idx++) {
            diff_eq(mqsystem + eq_idx * MQ_XVAR_NUM, var_idx, derivs[eq_idx][var_idx]);
        }
    }
}

__host__ __device__ static void
reduce_sys(uint8_t *mqsystem) {
    uint32_t eq_idx, var_idx, i, sqr_term_idx;
    for (eq_idx = 0; eq_idx < MQ_EQ_NUM; eq_idx++) {
        for (var_idx = 0; var_idx < MQ_VAR_NUM; var_idx++) {
            for (i = 0, sqr_term_idx = 0; i < var_idx; i++) {
                sqr_term_idx += i + 2;
            }
            mqsystem[eq_idx * MQ_XVAR_NUM + (MQ_XVAR_NUM - 2 - (MQ_VAR_NUM - 1 - var_idx))] ^=
                    mqsystem[eq_idx * MQ_XVAR_NUM + sqr_term_idx];
            mqsystem[eq_idx * MQ_XVAR_NUM + sqr_term_idx] = 0;
        }
    }
}

__device__ void __forceinline__
fast_exhaustive(uint8_t *mqsystem, uint8_t *solution) {
    uint64_t pdiff_eval[MQ_VAR_NUM];
    uint64_t func_eval = 0x0UL;
    uint32_t pre_fp_idx = 0;
    uint32_t count = 0;
    uint32_t fp_idx = 0;
    const uint32_t bound = (0x1U << MQ_VAR_NUM) - 1;
    uint64_t pdiff2[MQ_VAR_NUM][MQ_VAR_NUM];

    uint8_t derivs[MQ_EQ_NUM][MQ_VAR_NUM][MQ_VAR_NUM + 1];
    find_partial_derivs(mqsystem, derivs);

    uint32_t eq_idx, var_idx, i;
    uint64_t term;
    for (var_idx = 0; var_idx < MQ_VAR_NUM; var_idx++) {
        memset(pdiff2[var_idx], 0x0, sizeof(uint64_t) * MQ_VAR_NUM);
        for (i = 0; i < MQ_VAR_NUM; i++) {
            for (eq_idx = 0; eq_idx < MQ_EQ_NUM; eq_idx++) {
                term = derivs[eq_idx][var_idx][i];
                pdiff2[i][var_idx] |= term << eq_idx;
            }
        }
    }

    memset(pdiff_eval, 0x0, sizeof(uint64_t) * MQ_VAR_NUM);
    for (var_idx = 0; var_idx < MQ_VAR_NUM; var_idx++) {
        for (eq_idx = 0; eq_idx < MQ_EQ_NUM; eq_idx++) {
            if (var_idx == 0) {
                term = mqsystem[eq_idx * MQ_XVAR_NUM + deg2midx1(MQ_VAR_NUM, 0)];
            } else {
                term = mqsystem[eq_idx * MQ_XVAR_NUM + deg2midx1(MQ_VAR_NUM, var_idx)]
                       ^ derivs[eq_idx][var_idx][var_idx - 1];
            }
            pdiff_eval[var_idx] |= term << eq_idx;
        }
    }

    // brute forcing
    for (eq_idx = 0; eq_idx < MQ_EQ_NUM; eq_idx++) {
        term = mqsystem[eq_idx * MQ_XVAR_NUM + MQ_XVAR_NUM - 1];
        func_eval |= term << eq_idx;
    }

    while (func_eval && count < bound) {
        count++;
        fp_idx = ctz(count);

        if (count & (count - 1)) {
            pre_fp_idx = ctz(count ^ (0x1U << fp_idx));
            pdiff_eval[fp_idx] ^= pdiff2[fp_idx][pre_fp_idx];
        }

        func_eval ^= pdiff_eval[fp_idx];
    }

    if (!func_eval) {
        for (var_idx = 0; var_idx < MQ_VAR_NUM; var_idx++) {
            solution[var_idx] = (uint8_t) (((count ^ (count >> 1)) >> var_idx) & 0x1U);
        }
    } else {
        memset(solution, 0x0, MQ_VAR_NUM * sizeof(uint8_t));
    }
}

__host__ void
loadSystemsFromFile(KeccakSolver *keccakSolver) {
    char ch;
    uint32_t i, j;

    uint8_t constant_constr[LIN_CONST_EQNUM][801];
    uint8_t iterative_constr[LIN_ITER_EQNUM][801];
    uint8_t *mq_system[MQ_EQ_NUM];
    for (i = 0; i < MQ_EQ_NUM; i++)
        mq_system[i] = (uint8_t *) malloc(BMQ_XVAR_NUM * sizeof(uint8_t));
    uint8_t *append_system[AMQ_LIN_EQNUM];
    for (i = 0; i < AMQ_LIN_EQNUM; i++)
        append_system[i] = (uint8_t *) malloc(BMQ_XVAR_NUM * sizeof(uint8_t));

    PRINTF_STAMP("[+] Reading preprocessed systems from file\n");
    FILE *flin = fopen(keccakSolver->options.c_lin_analysis_file, "r");
    FILE *fmq = fopen(keccakSolver->options.mq_analysis_file, "r");
    FILE *fa = fopen(keccakSolver->options.a_lin_analysis_file, "r");
    FILE *fi = fopen(keccakSolver->options.i_lin_analysis_file, "r");

    if (flin == NULL) {
        EXIT_WITH_MSG("[!] cannot open constant linear constraint file\n");
    } else {
        for (i = 0; i < LIN_CONST_EQNUM; i++) {
            for (j = 0; j < 801; j++) {
                ch = fgetc(flin);
                constant_constr[i][j] = (uint8_t) (ch - '0');
            }
            fgetc(flin);
        }
    }
    PRINTF_STAMP("\t\tread constant linear constraints from file\n");

    if (fmq == NULL) {
        EXIT_WITH_MSG("[!] cannot open mq analysis file\n");
    } else {
        for (i = 0; i < MQ_EQ_NUM; i++) {
            for (j = 0; j < BMQ_XVAR_NUM; j++) {
                ch = fgetc(fmq);
                mq_system[i][j] = (uint8_t) (ch - '0');
            }
            fgetc(fmq);
        }
    }
    PRINTF_STAMP("\t\tread mq analysis from file\n");

    if (fi == NULL) {
        EXIT_WITH_MSG("[!] cannot open iterative linear constraint file\n");
    } else {
        for (i = 0; i < LIN_ITER_EQNUM; i++) {
            for (j = 0; j < 801; j++) {
                ch = fgetc(fi);
                iterative_constr[i][j] = (uint8_t) (ch - '0');
            }
            fgetc(fi);
        }
    }
    PRINTF_STAMP("\t\tread iterative linear constraints from file\n");

    if (fa == NULL) {
        EXIT_WITH_MSG("[!] cannot open append system file\n");
    } else {
        for (i = 0; i < AMQ_LIN_EQNUM; i++) {
            for (j = 0; j < BMQ_XVAR_NUM; j++) {
                ch = fgetc(fa);
                append_system[i][j] = (uint8_t) (ch - '0');
            }
            fgetc(fa);
        }
    }

    fclose(fmq);
    fclose(flin);
    fclose(fi);
    fclose(fa);
    PRINTF_STAMP("[+] Preprocessed systems have been load to program\n");
    PRINTF_STAMP("[+] Start reducing system to round 3\n");

    r3mqarg_t mqargs[MQ_EQ_NUM];
    r3aparg_t apargs[AMQ_LIN_EQNUM];

    extractRound3LinearDependency(&keccakSolver->mathSystem, constant_constr);
    reduceIterativeConstraints(&keccakSolver->mathSystem, iterative_constr);

    // reduceRound3AppendSystem(&keccakSolver->mathSystem, append_system);
    for (i = 0; i < AMQ_LIN_EQNUM; i++) {
        apargs[i].eq_idx = i;
        apargs[i].mathSystem = &keccakSolver->mathSystem;
        apargs[i].append_system = append_system;
        threadpool_add(keccakSolver->threadpool, reduceRound3AppendSystem, (void *) &apargs[i], 0);
    }

    // reduceRound3MQSystem(&keccakSolver->mathSystem, mq_system);
    for (i = 0; i < MQ_EQ_NUM; i++) {
        mqargs[i].eq_idx = i;
        mqargs[i].mathSystem = &keccakSolver->mathSystem;
        mqargs[i].mqsystem = mq_system;
        threadpool_add(keccakSolver->threadpool, reduceRound3MQSystem, (void *) &mqargs[i], 0);
    }
    threadpool_join(keccakSolver->threadpool, 0);

    uint8_t eval[800] = {0};
    for (i = 0; i < 800; i++) {
        if (keccakSolver->mathSystem.round3_lin2mq[i] == DEP_PLACEMENT) {
            for (j = 0; j < IMQ_VAR_NUM + 1; j++) {
                eval[i] ^= keccakSolver->mathSystem.round3_lin_dep[i][j];
            }
        } else {
            eval[i] = 1;
        }
    }

    // validate the reduction result
    for (uint32_t eq_idx = 0; eq_idx < LIN_ITER_EQNUM; eq_idx++) {
        uint32_t verifyResult1 = 0;
        uint32_t verifyResult2 = 0;
        for (i = 0; i < IMQ_VAR_NUM + 1; i++) {
            if (keccakSolver->mathSystem.round3_iter_system[eq_idx][i])
                verifyResult1 ^= 1;
        }
        for (i = 0; i < 800; i++) {
            if (iterative_constr[eq_idx][i])
                verifyResult2 ^= eval[i];
        }
        if (iterative_constr[eq_idx][800])
            verifyResult2 ^= 1;
        if (verifyResult1 != verifyResult2) {
            EXIT_WITH_MSG("reduced iterative constraints is not equal with the original ones, "
                          "please check the implementation.\n");
        }
    }

    for (uint32_t eq_idx = 0; eq_idx < MQ_EQ_NUM; eq_idx++) {
        uint32_t verifyResult1 = 0;
        uint32_t verifyResult2 = 0;
        for (i = 0; i < IMQ_VAR_NUM; i++) {
            for (j = i; j < IMQ_VAR_NUM; j++) {
                if (keccakSolver->mathSystem.round3_mq_system[eq_idx][deg2midx2(i, j)])
                    verifyResult1 ^= 1;
            }
            if (keccakSolver->mathSystem.round3_mq_system[eq_idx][deg2midx1(IMQ_VAR_NUM, i)])
                verifyResult1 ^= 1;
        }
        if (keccakSolver->mathSystem.round3_mq_system[eq_idx][IMQ_XVAR_NUM - 1])
            verifyResult1 ^= 1;

        for (i = 0; i < 800; i++) {
            for (j = i; j < 800; j++) {
                if (mq_system[eq_idx][deg2midx2(i, j)])
                    verifyResult2 ^= (eval[i] & eval[j]);
            }
            if (mq_system[eq_idx][deg2midx1(800, i)])
                verifyResult2 ^= eval[i];
        }
        if (mq_system[eq_idx][BMQ_XVAR_NUM - 1])
            verifyResult2 ^= 1;

        if (verifyResult1 != verifyResult2) {
            EXIT_WITH_MSG("reduced mq system is not equal with the original system, "
                          "please check the implementation.\n");
        }
    }

    for (uint32_t eq_idx = 0; eq_idx < AMQ_LIN_EQNUM; eq_idx++) {
        uint32_t verifyResult1 = 0;
        uint32_t verifyResult2 = 0;
        for (i = 0; i < IMQ_VAR_NUM; i++) {
            for (j = i; j < IMQ_VAR_NUM; j++) {
                if (keccakSolver->mathSystem.round3_append_system[eq_idx][deg2midx2(i, j)])
                    verifyResult1 ^= 1;
            }
            if (keccakSolver->mathSystem.round3_append_system[eq_idx][deg2midx1(IMQ_VAR_NUM, i)])
                verifyResult1 ^= 1;
        }
        if (keccakSolver->mathSystem.round3_append_system[eq_idx][IMQ_XVAR_NUM - 1])
            verifyResult1 ^= 1;

        for (i = 0; i < 800; i++) {
            for (j = i; j < 800; j++) {
                if (append_system[eq_idx][deg2midx2(i, j)])
                    verifyResult2 ^= (eval[i] & eval[j]);
            }
            if (append_system[eq_idx][deg2midx1(800, i)])
                verifyResult2 ^= eval[i];
        }
        if (append_system[eq_idx][BMQ_XVAR_NUM - 1])
            verifyResult2 ^= 1;

        if (verifyResult1 != verifyResult2) {
            EXIT_WITH_MSG("reduced append system is not equal with the original system, "
                          "please check the implementation.\n");
        }
    }

    uint32_t A[5][5];
    for (i = 0; i < 5; i++)
        memset(A[i], 0x0, 5 * sizeof(uint32_t));

    for (i = 0; i < 800; i++) {
        uint32_t idx_z = (uint32_t) (i % 32);
        uint32_t idx_y = (uint32_t) ((i / 32) % 5);
        uint32_t idx_x = (uint32_t) ((i / 32) / 5);

        A[idx_x][idx_y] |= (eval[i] << idx_z);
    }
    cpu_VerifyRound2(A);

    PRINTF_STAMP("[+] All reductions have finished, inner algebraic system are set to round 3\n");

    for (i = 0; i < MQ_EQ_NUM; i++)
        SFREE(mq_system[i]);
    for (i = 0; i < AMQ_LIN_EQNUM; i++)
        SFREE(append_system[i]);
}

__global__ void
kernelLoop(uint8_t *device_output_buffer,
           uint8_t *device_mq_buffer) {
    uint64_t thread_id = (threadIdx.x + blockIdx.x * blockDim.x);
    uint8_t *kern_mq_buffer = device_mq_buffer + thread_id * (MQ_SYSTEM_SIZE);
    uint8_t *kern_output_buffer = device_output_buffer + thread_id * MQ_VAR_NUM;
    fast_exhaustive(kern_mq_buffer, kern_output_buffer);
}

void
readMessageFromFile(FILE *message_file,
                    uint32_t **message) {
    char line[256];
    PRINTF_STAMP("Read message from file: \n");
    while (fgets(line, sizeof(line), message_file)) {
        uint32_t i, j, val;
        sscanf(line, "A[%d][%d]: 0x%08x", &i, &j, &val);
        message[i][j] = val;
        PRINTF_STAMP("\tA[%d][%d]: 0x%08x\n", i, j, message[i][j]);
    }
}

__host__ void
keccakSolverInit(KeccakSolver *keccakSolver, int argc, char **argv) {
    PRINTF_STAMP("[+] Initialize Keccak Solver\n");
    PRINTF_STAMP("[+] Parse command line options\n");
    options_init(&keccakSolver->options);
    options_parse(&keccakSolver->options, argc, argv);
    PRINTF_STAMP("[+] Command line options have been parse\n");

    PRINTF_STAMP("[+] Initial math system\n");
    initMathSystem(&keccakSolver->mathSystem);

    keccakSolver->threadpool = threadpool_create(keccakSolver->options.cpu_thread_num, 0x2000, 0);

    const size_t device_mbufer_size = CHUNK_SIZE * MQ_SYSTEM_SIZE * sizeof(uint8_t);
    const size_t device_obufer_size = CHUNK_SIZE * MQ_VAR_NUM * sizeof(uint8_t);
    PRINTF_STAMP("[+] Allocate buffer on GPU\n");
    PRINTF_STAMP("\t\tdevice_mbuffer_size = 0x%lx bytes\n", device_mbufer_size);
    PRINTF_STAMP("\t\tdevice_obuffer_size = 0x%lx bytes\n", device_obufer_size);

    CUDA_CHECK(cudaSetDevice(keccakSolver->options.dev_id));
    CUDA_CHECK(cudaDeviceSynchronize());
    CUDA_CHECK(cudaDeviceReset());
    CUDA_CHECK(cudaMalloc(&keccakSolver->device_mq_buffer, device_mbufer_size));
    CUDA_CHECK(cudaMalloc(&keccakSolver->device_output_buffer, device_obufer_size));

    loadSystemsFromFile(keccakSolver);

#ifdef TEST_PRE
    FILE *message_file = fopen(keccakSolver->options.test_message_file, "r");
    if (message_file == nullptr) {
        PRINTF_ERR_STAMP("cannot open testMessage file, exit.");
        exit(1);
    }
    keccakSolver->test_message = (uint32_t **) malloc(5 * sizeof(uint32_t *));
    for (uint32_t i = 0; i < 5; i++) {
        keccakSolver->test_message[i] = (uint32_t *) malloc(5 * sizeof(uint32_t));
        memset(keccakSolver->test_message[i], 0, 5 * sizeof(uint32_t));
    }
    readMessageFromFile(message_file, keccakSolver->test_message);
    fclose(message_file);
#endif
}

__host__ void
threadUpdateMQSystem(void *arg) {
    tmqarg_t *args = (tmqarg_t *) arg;
    uint64_t guessingBits = args->guessingBits;
    uint8_t *mqbuffer = args->mqbuffer;
    uint32_t *mq2lin = args->mq2lin;
    uint32_t *lin2mq = args->lin2mq;
    uint8_t *lindep = args->lindep;
    MathSystem *mathSystem = args->mathSystem;
    guessingBitsToMqSystem(mathSystem, guessingBits, mqbuffer, mq2lin, lin2mq, lindep);
    reduce_sys(mqbuffer);
}

__host__ static inline bool
verify_sol(uint8_t *solution, uint8_t *sys, const uint64_t eq_num,
           const uint64_t var_num, const uint64_t term_num,
           const uint64_t start) {
    for (uint64_t i = 0; i < eq_num; i++) { // for each equation
        uint64_t res = 0;

        for (uint64_t mul_1 = 0; mul_1 < var_num; mul_1++) {
            for (uint64_t mul_2 = mul_1; mul_2 < var_num; mul_2++) {
                if (sys[i * term_num + deg2midx2(mul_1, mul_2)] == 1) {
                    res ^= (solution[mul_1] & solution[mul_2]);
                }
            }
        }

        for (uint64_t var_idx = 0; var_idx < var_num; var_idx++) {
            res ^= (solution[var_idx] & sys[i * term_num + deg2midx1(var_num, var_idx)]);
        }

        res ^= sys[i * term_num + term_num - 1];
        if (res == 1) { // the equation is evaluated to 1
            return false;
        }
    }

    return true;
}

__host__ uint64_t
precalculate_round3_value(KeccakSolver *keccakSolver) {
    uint64_t ret = 0;
    uint64_t eq_val = 0;
    uint64_t i, j;
    uint8_t **round3_iter_system = keccakSolver->mathSystem.round3_iter_system;
    uint32_t *round3_mq2lin = keccakSolver->mathSystem.round3_mq2lin;
    uint32_t **message = keccakSolver->test_message;

    for (i = 0; i < LIN_ITER_EQNUM; i++) {
        eq_val = 0;
        for (j = 0; j < IMQ_VAR_NUM; j++) {
            uint32_t var_idx = round3_mq2lin[j];
            uint32_t idx_x = (uint32_t) ((var_idx / 32) / 5);
            uint32_t idx_y = (uint32_t) ((var_idx / 32) % 5);
            uint32_t idx_z = (uint32_t) (var_idx % 32);
            eq_val ^= (round3_iter_system[i][j] & ((message[idx_x][idx_y] >> idx_z) & 0x1U));
        }
        eq_val ^= round3_iter_system[i][IMQ_VAR_NUM];
        ret |= (eq_val << i);
    }
    return ret;
}

__host__ bool
validate_iter_round3_value(uint64_t candidate, uint64_t precal_val) {
    // accept candidates with difference only on lowest 3 bits
    uint64_t filter_out = ~0x7UL;
    return (candidate & filter_out) == (precal_val & filter_out);
}

__host__ void
threadCheckResult(void *arg) {
    tcheckarg_t *args = (tcheckarg_t *) arg;
    uint8_t *result_buffer = args->result_buffer;
    uint32_t *lin2mq = args->lin2mq;
    uint8_t *lindep = args->lindep;
    uint32_t *minDiff = args->minDiff;
    uint8_t *mq_buffer = args->mqbuffer;

    uint32_t i, j, idx_x, idx_y, idx_z;
    bool result_found = false;
    for (i = 0; i < MQ_VAR_NUM; i++) {
        if (result_buffer[i] != 0x0) {
            result_found = true;
            break;
        }
    }

    if (result_found) {
        if (!verify_sol(result_buffer, mq_buffer, MQ_EQ_NUM, MQ_VAR_NUM, MQ_XVAR_NUM, 0))
            PRINTF_ERR_STAMP("mq system solve error\n");
        uint32_t A[5][5];
        for (i = 0; i < 5; i++)
            memset(A[i], 0x0, 5 * sizeof(uint32_t));
        for (i = 0; i < 800; i++) {
            idx_z = (uint32_t) (i % 32);
            idx_y = (uint32_t) ((i / 32) % 5);
            idx_x = (uint32_t) ((i / 32) / 5);
            uint32_t val = 0;

            if (lin2mq[i] == DEP_PLACEMENT) {
                for (j = 0; j < MQ_VAR_NUM; j++) {
                    if (result_buffer[j])
                        val ^= (lindep[i * (MQ_VAR_NUM + 1) + j]);
                }
                val ^= (lindep[i * (MQ_VAR_NUM + 1) + MQ_VAR_NUM]);
                A[idx_x][idx_y] |= (val << idx_z);
            } else {
                if (result_buffer[lin2mq[i]])
                    A[idx_x][idx_y] |= (1 << idx_z);
            }
        }

//        printStatus(A);

        if (cpu_VerifyKeccakResult(A, minDiff)) {
            PRINTF_STAMP("\t\tone thread found valid preimage: ");
            printStatus(A);
            *args->preimage_found = true;
        }
    }
}


__host__ void
keccakSolverLoop(KeccakSolver *keccakSolver) {
#ifdef TEST_PRE
    uint64_t precal_val = precalculate_round3_value(keccakSolver);
    uint8_t mqsystem[MQ_SYSTEM_SIZE] = {0};
    uint8_t result[MQ_VAR_NUM] = {0};
    uint8_t lin_dep[800 * (MQ_VAR_NUM + 1)] = {0};
    uint32_t mq2lin[MQ_VAR_NUM] = {0};
    uint32_t lin2mq[800] = {0};
    guessingBitsToMqSystem(&(keccakSolver->mathSystem), precal_val, mqsystem, mq2lin, lin2mq, lin_dep);
    reduce_sys(mqsystem);

    // a fake solve of mq system: use known message to substitute mq variables
    for (uint32_t i = 0; i < MQ_VAR_NUM; i++) {
        uint32_t idx_z = (uint32_t) (mq2lin[i] % 32);
        uint32_t idx_y = (uint32_t) ((mq2lin[i] / 32) % 5);
        uint32_t idx_x = (uint32_t) ((mq2lin[i] / 32) / 5);
        result[i] = (uint8_t) ((keccakSolver->test_message[idx_x][idx_y] >> idx_z) & 1);
    }

    uint32_t A[5][5];
    for (uint32_t i = 0; i < 5; i++)
        memset(A[i], 0x0, 5 * sizeof(uint32_t));
    for (uint32_t i = 0; i < 800; i++) {
        uint32_t idx_z = (uint32_t) (i % 32);
        uint32_t idx_y = (uint32_t) ((i / 32) % 5);
        uint32_t idx_x = (uint32_t) ((i / 32) / 5);
        uint32_t val = 0;

        if (lin2mq[i] == DEP_PLACEMENT) {
            for (uint32_t j = 0; j < MQ_VAR_NUM; j++) {
                if (result[j])
                    val ^= (lin_dep[i * (MQ_VAR_NUM + 1) + j]);
            }
            val ^= (lin_dep[i * (MQ_VAR_NUM + 1) + MQ_VAR_NUM]);
            A[idx_x][idx_y] |= (val << idx_z);
        } else {
            if (result[lin2mq[i]])
                A[idx_x][idx_y] |= (1 << idx_z);
        }
    }

    printStatus(A);

    if (!verify_sol(result, mqsystem, MQ_EQ_NUM, MQ_VAR_NUM, MQ_XVAR_NUM, 0))
        PRINTF_ERR_STAMP("[t] Testing: mq system solve error.\n");
    else {
        PRINTF_STAMP("[t] Testing: mq system solve correct!\n");
    }
#else
    CUDA_CHECK(cudaDeviceSynchronize());
    dim3 tpb(GPU_THREADS_PER_BLOCK);

    /* set the smallest searching space of guessing bits */
    uint64_t search_interval = keccakSolver->options.gbend - keccakSolver->options.gbstart;
    search_interval = search_interval > CHUNK_SIZE ? search_interval : CHUNK_SIZE;
    uint64_t gbstart = keccakSolver->options.gbstart;
    uint64_t gbend = gbstart + search_interval;

    /* temporary memory buffer for copy data back from kernel */
    size_t mbuffer_size = (CHUNK_SIZE * MQ_SYSTEM_SIZE) * sizeof(uint8_t);
    size_t rbuffer_size = (CHUNK_SIZE * MQ_VAR_NUM) * sizeof(uint8_t);

    /* temporary memory buffer for validation */
    size_t lindep_buffer_size = (CHUNK_SIZE * 800 * (MQ_VAR_NUM + 1)) * sizeof(uint8_t);
    size_t mq2lin_buffer_size = (CHUNK_SIZE * MQ_VAR_NUM) * sizeof(uint32_t);
    size_t lin2mq_buffer_size = (CHUNK_SIZE * 800) * sizeof(uint32_t);

    uint8_t *mqbuffer = (uint8_t *) malloc(mbuffer_size);
    uint8_t *result_buffer = (uint8_t *) malloc(rbuffer_size);
    uint8_t *lin_dep_buffer = (uint8_t *) malloc(lindep_buffer_size);
    uint32_t *mq2lin_buffer = (uint32_t *) malloc(mq2lin_buffer_size);
    uint32_t *lin2mq_buffer = (uint32_t *) malloc(lin2mq_buffer_size);

    bool preimage_found = false;

        /* exhaustively search between gbstart and gbend */
        for (uint64_t gb = gbstart; gb < gbend; gb += CHUNK_SIZE) {

        CUDA_CHECK(cudaMemset(keccakSolver->device_mq_buffer, 0, mbuffer_size));
        CUDA_CHECK(cudaMemset(keccakSolver->device_output_buffer, 0, rbuffer_size));
        PRINTF_STAMP("[+] Init new guessing bits, starts at: 0x%lx, ends at: 0x%lx\n", gb, gb + CHUNK_SIZE);
        PRINTF_STAMP("\t\tupdating mq system\n");

        /* update MQ System */
        memset(mqbuffer, 0, CHUNK_SIZE * MQ_SYSTEM_SIZE);
        memset(result_buffer, 0, rbuffer_size);
        memset(mqbuffer, 0, mbuffer_size);
        memset(lin_dep_buffer, 0, lindep_buffer_size);
        memset(mq2lin_buffer, 0, mq2lin_buffer_size);
        memset(lin2mq_buffer, 0, lin2mq_buffer_size);

        tmqarg_t args[CHUNK_SIZE];
        uint64_t this_gb;
        for (this_gb = gb; this_gb < gb + CHUNK_SIZE; this_gb++) {
            uint64_t thread_id = this_gb - gb;
            args[thread_id].guessingBits = this_gb;
            args[thread_id].mqbuffer = mqbuffer + thread_id * MQ_SYSTEM_SIZE;
            args[thread_id].lindep = lin_dep_buffer + thread_id * (800 * (MQ_VAR_NUM + 1));
            args[thread_id].mq2lin = mq2lin_buffer + thread_id * (MQ_VAR_NUM);
            args[thread_id].lin2mq = lin2mq_buffer + thread_id * 800;
            args[thread_id].mathSystem = &keccakSolver->mathSystem;
            threadpool_add(keccakSolver->threadpool, threadUpdateMQSystem, (void *) &args[thread_id], 0);
        }
        threadpool_join(keccakSolver->threadpool, 0);

        PRINTF_STAMP("\t\tlaunch GPU kernel to solve mq systems\n");
        /* use gpu to solve mq system*/
        CUDA_CHECK(cudaDeviceSynchronize());
        CUDA_CHECK(cudaMemcpy(keccakSolver->device_mq_buffer, mqbuffer, mbuffer_size, cudaMemcpyHostToDevice));
        kernelLoop << < GPU_BLOCK_NUM, tpb >> > (keccakSolver->device_output_buffer, keccakSolver->device_mq_buffer);
        CUDA_CHECK(cudaDeviceSynchronize());

        PRINTF_STAMP("\t\tGPU kernel finished\n");

        /* check output result */
        PRINTF_STAMP("\t\tverifying result...\n");
        CUDA_CHECK(cudaMemcpy(result_buffer, keccakSolver->device_output_buffer, rbuffer_size, cudaMemcpyDeviceToHost));
        tcheckarg_t check_args[CHUNK_SIZE];
        uint32_t minDiff[CHUNK_SIZE];
        memset(minDiff, 80, CHUNK_SIZE * sizeof(uint32_t));
        for (uint64_t thread_id = 0; thread_id < CHUNK_SIZE; thread_id++) {
            check_args[thread_id].preimage_found = &preimage_found;
            check_args[thread_id].result_buffer = result_buffer + thread_id * MQ_VAR_NUM;
            check_args[thread_id].lindep = lin_dep_buffer + thread_id * (800 * (MQ_VAR_NUM + 1));
            check_args[thread_id].mq2lin = mq2lin_buffer + thread_id * (MQ_VAR_NUM);
            check_args[thread_id].lin2mq = lin2mq_buffer + thread_id * 800;
            check_args[thread_id].mqbuffer = mqbuffer + thread_id * MQ_SYSTEM_SIZE;
            check_args[thread_id].minDiff = &minDiff[thread_id];
            threadpool_add(keccakSolver->threadpool, threadCheckResult, (void *) &check_args[thread_id], 0);
        }
        threadpool_join(keccakSolver->threadpool, 0);

        /* check flag */
        if (preimage_found) {
            PRINTF_STAMP("[-] Valid pre-image found\n");
            PRINTF_STAMP("[-] Terminating...\n");
            break;
        } else {
            PRINTF_STAMP("[+] No valid pre-image found\n");
            uint32_t m = 80;
            for (uint32_t i = 0; i < CHUNK_SIZE; i++) {
                if (minDiff[i] < m)
                    m = minDiff[i];
            }
            PRINTF_STAMP("\t\tminimal differ last round: %d\n", m);
        }
    }
    free(mqbuffer);
    free(result_buffer);
    free(lin_dep_buffer);
    free(mq2lin_buffer);
    free(lin2mq_buffer);
#endif
}

__host__ void
keccakSolverExit(KeccakSolver *keccakSolver) {
    options_free(&keccakSolver->options);
    freeMathSystem(&keccakSolver->mathSystem);
    threadpool_destroy(keccakSolver->threadpool, 0);
    CUDA_CHECK(cudaDeviceSynchronize());
    CUDA_CHECK(cudaFree(keccakSolver->device_mq_buffer));
    CUDA_CHECK(cudaFree(keccakSolver->device_output_buffer));
    EXIT_WITH_MSG("[-] Keccak solver exit\n");

#ifdef TEST_PRE
    for (uint32_t i = 0; i < 5; i++) {
        SFREE(keccakSolver->test_message[i]);
    }
    SFREE(keccakSolver->test_message);
#endif
}
