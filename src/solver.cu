//
// Created by vitowu on 3/13/20.
//

#include "solver.h"

#include "mq.h"
#include "threadpool.h"
#include "params.h"

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
    for (i = 0; i < MQ_EQ_NUM; i++)
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

    PRINTF_STAMP("[+] All reductions have finished, inner algebraic system are set to round 3\n");

    for (i = 0; i < MQ_EQ_NUM; i++)
        SFREE(mq_system[i]);
    for (i = 0; i < AMQ_LIN_EQNUM; i++)
        SFREE(append_system[i]);
}

__global__ void
kernelLoop(uint32_t *device_output_buffer,
           uint8_t *device_mq_buffer) {
    uint64_t thread_id = (threadIdx.x + blockIdx.x * blockDim.x);

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

    keccakSolver->threadpool = threadpool_create(keccakSolver->options.cpu_thread_num, CHUNK_SIZE, 0);

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
}

__host__ void
threadCheckResult(void *arg) {
    tcheckarg_t *args = (tcheckarg_t *) arg;
    *args->preimage_found = true;
}

__host__ void
keccakSolverLoop(KeccakSolver *keccakSolver) {
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
    size_t lindep_buffer_size = (CHUNK_SIZE * 800 * (MQ_VAR_NUM + 1)) * sizeof(uint32_t);
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
        PRINTF_STAMP("[+] Init new guessing bits, starts at: 0x%lx\n", gb);
        PRINTF_STAMP("\t\tupdating mq system\n");

        /* update MQ System */
        memset(mqbuffer, 0, CHUNK_SIZE * MQ_SYSTEM_SIZE);
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
        PRINTF_STAMP("\t\tverifying result...\n");
        /* check output result */
        CUDA_CHECK(cudaMemcpy(result_buffer, keccakSolver->device_output_buffer, rbuffer_size, cudaMemcpyDeviceToHost));
        tcheckarg_t check_args[CHUNK_SIZE];
        for (uint64_t thread_id = 0; thread_id < CHUNK_SIZE; thread_id++) {
            check_args[thread_id].preimage_found = &preimage_found;
            check_args[thread_id].result_buffer = result_buffer + thread_id * MQ_VAR_NUM;
            check_args[thread_id].lindep = lin_dep_buffer + thread_id * (800 * (MQ_VAR_NUM + 1));
            check_args[thread_id].mq2lin = mq2lin_buffer + thread_id * (MQ_VAR_NUM);
            check_args[thread_id].lin2mq = lin2mq_buffer + thread_id * 800;
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
            PRINTF_STAMP("[+] Start a new chunk.\n");
        }
    }
    free(mqbuffer);
    free(result_buffer);
    free(lin_dep_buffer);
    free(mq2lin_buffer);
    free(lin2mq_buffer);
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
}
