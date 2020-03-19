//
// Created by vitowu on 3/13/20.
//

#include "solver.h"

#include "threadpool.h"
#include "params.h"


// linear dependency passed to thread
__constant__ uint8_t lin_dep[800][MQ_VAR_NUM + 1];

// mq index (upto 31) to linear index (upto 800)
__constant__ uint32_t mqidx2lin[MQ_VAR_NUM];

__host__ void
loadSystemsFromFile(KeccakSolver *keccakSolver) {
    char ch;
    uint32_t i, j;

    uint8_t constant_constr[LIN_ITER_EQNUM][801];
    uint8_t iterative_constr[LIN_CONST_EQNUM][801];
    uint8_t *mq_system[MQ_EQ_NUM];
    for (i = 0; i < MQ_EQ_NUM; i++)
        mq_system[i] = (uint8_t *) malloc(BMQ_XVAR_NUM * sizeof(uint8_t));
    uint8_t *append_system[AMQ_LIN_EQNUM];
    for (i = 0; i < MQ_EQ_NUM; i++)
        append_system[i] = (uint8_t *) malloc(BMQ_XVAR_NUM * sizeof(uint8_t));

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
    PRINTF_STAMP("read constant linear constraints from file\n");

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
    PRINTF_STAMP("read mq analysis from file\n");

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
    PRINTF_STAMP("read iterative linear constraints from file\n");

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

    extractRound3LinearDependency(&keccakSolver->mathSystem, constant_constr);
    reduceIterativeConstraints(&keccakSolver->mathSystem, iterative_constr);
    reduceRound3AppendSystem(&keccakSolver->mathSystem, append_system);
    reduceRound3MQSystem(&keccakSolver->mathSystem, mq_system);
    PRINTF_STAMP("all reductions have finished\n");

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
    options_init(&keccakSolver->options);
    options_parse(&keccakSolver->options, argc, argv);

    initMathSystem(&keccakSolver->mathSystem);

    keccakSolver->threadpool = threadpool_create(keccakSolver->options.cpu_thread_num, CHUNK_SIZE, 0);

    const size_t device_mbufer_size = CHUNK_SIZE * MQ_SYSTEM_SIZE * sizeof(uint8_t);
    const size_t device_obufer_size = CHUNK_SIZE * MQ_VAR_NUM * sizeof(uint8_t);
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
    uint8_t *lindep = args->lindep;
    MathSystem *mathSystem = args->mathSystem;
    guessingBitsToMqSystem(mathSystem, guessingBits, mqbuffer, mq2lin, lindep);
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

    uint8_t *mqbuffer = (uint8_t *) malloc(mbuffer_size);
    uint8_t *result_buffer = (uint8_t *) malloc(rbuffer_size);
    uint8_t *lin_dep_buffer = (uint8_t *) malloc(lindep_buffer_size);
    uint32_t *mq2lin_buffer = (uint32_t *) malloc(mq2lin_buffer_size);

    bool preimage_found = false;

    /* exhaustively search between gbstart and gbend */
    for (uint64_t gb = gbstart; gb < gbend; gb += CHUNK_SIZE) {

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
            args[thread_id].mathSystem = &keccakSolver->mathSystem;
            threadpool_add(keccakSolver->threadpool, threadUpdateMQSystem, (void *) &args[thread_id], 0);
        }
        threadpool_join(keccakSolver->threadpool, 0);


        /* use gpu to solve mq system*/
        CUDA_CHECK(cudaDeviceSynchronize());
        CUDA_CHECK(cudaMemcpy(keccakSolver->device_mq_buffer, mqbuffer, mbuffer_size, cudaMemcpyHostToDevice));
        kernelLoop << < GPU_BLOCK_NUM, tpb >> > (keccakSolver->device_output_buffer, keccakSolver->device_mq_buffer);
        CUDA_CHECK(cudaDeviceSynchronize());


        /* check output result */
        CUDA_CHECK(cudaMemcpy(result_buffer, keccakSolver->device_output_buffer, rbuffer_size, cudaMemcpyDeviceToHost));
        tcheckarg_t check_args[CHUNK_SIZE];
        for (uint64_t thread_id = 0; thread_id < CHUNK_SIZE; thread_id++) {
            check_args[thread_id].preimage_found = &preimage_found;
            check_args[thread_id].result_buffer = result_buffer + thread_id * MQ_VAR_NUM;
            check_args[thread_id].lindep = lin_dep_buffer + thread_id * (800 * (MQ_VAR_NUM + 1));
            check_args[thread_id].mq2lin = mq2lin_buffer + thread_id * (MQ_VAR_NUM);
            threadpool_add(keccakSolver->threadpool, threadCheckResult, (void *) &check_args[thread_id], 0);
        }
        threadpool_join(keccakSolver->threadpool, 0);

        /* check flag */
        if (preimage_found) {
            PRINTF_STAMP("one valid pre-image found!\n");
            break;
        }
    }
    free(mqbuffer);
    free(result_buffer);
    free(lin_dep_buffer);
    free(mq2lin_buffer);
}

__host__ void
keccakSolverExit(KeccakSolver *keccakSolver) {
    options_free(&keccakSolver->options);
    freeMathSystem(&keccakSolver->mathSystem);
    threadpool_destroy(keccakSolver->threadpool, 0);
    CUDA_CHECK(cudaDeviceSynchronize());
    CUDA_CHECK(cudaFree(keccakSolver->device_mq_buffer));
    CUDA_CHECK(cudaFree(keccakSolver->device_output_buffer));
    EXIT_WITH_MSG("keccak solver exit\n");
}