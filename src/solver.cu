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

    // load mq system
    FILE *fmq = fopen(keccakSolver->options.mq_analysis_file, "r");
    uint32_t *file_mq_system = (uint32_t *) malloc(MQ_EQ_NUM * BMQ_XVAR_NUM * sizeof(uint32_t));

    if (fmq == NULL) {
        EXIT_WITH_MSG("[!] cannot open mq analysis file\n");
    } else {
        for (i = 0; i < MQ_EQ_NUM; i++) {
            for (j = 0; j < BMQ_XVAR_NUM; j++) {
                ch = fgetc(fmq);
                file_mq_system[i * BMQ_XVAR_NUM + j] = (uint8_t) (ch - '0');
            }
            fgetc(fmq);
        }
    }
    fclose(fmq);
    PRINTF_STAMP("[+] read mq analysis from file\n");

    free(file_mq_system);
//    // load linear system
//    FILE *flin = fopen(keccakSolver->options.c_lin_analysis_file, "r");
//    uint8_t local_c_constr[LIN_CONST_SYSTEM_SIZE];
//    if (flin == NULL) {
//        EXIT_WITH_MSG("[!] cannot open constant linear constraint file\n");
//    } else {
//        for (i = 0; i < LIN_CONST_EQNUM; i++) {
//            for (j = 0; j < 801; j++) {
//                ch = fgetc(flin);
//                local_c_constr[i * 801 + j] = (uint8_t) (ch - '0');
//            }
//            fgetc(flin);
//        }
//    }
//    fclose(flin);
//    PRINTF_STAMP("[+] read constant linear constraints from file\n");
//
//    FILE *fi = fopen(keccakSolver->options.i_lin_analysis_file, "r");
//    uint8_t local_i_constr[LIN_ITER_SYSTEM_SIZE];
//    if (fi == NULL) {
//        EXIT_WITH_MSG("[!] cannot open iterative linear constraint file\n");
//    } else {
//        for (i = 0; i < LIN_ITER_EQNUM; i++) {
//            for (j = 0; j < 801; j++) {
//                ch = fgetc(fi);
//                local_i_constr[i * 801 + j] = (uint8_t) (ch - '0');
//            }
//            fgetc(fi);
//        }
//    }
//    fclose(fi);
//    PRINTF_STAMP("[+] read iterative linear constraints from file\n");
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

    CUDA_CHECK(cudaSetDevice(keccakSolver->options.dev_id));
    CUDA_CHECK(cudaDeviceSynchronize());
    CUDA_CHECK(cudaDeviceReset());

    keccakSolver->threadpool = threadpool_create(keccakSolver->options.cpu_thread_num, CHUNK_SIZE, 0);

    const size_t device_mbufer_size = CHUNK_SIZE * MQ_SYSTEM_SIZE * sizeof(uint8_t);
    const size_t device_obufer_size = CHUNK_SIZE * 25 * sizeof(uint32_t);
    CUDA_CHECK(cudaMalloc(&keccakSolver->device_mq_buffer, device_mbufer_size));
    CUDA_CHECK(cudaMalloc(&keccakSolver->device_output_buffer, device_obufer_size));

    loadSystemsFromFile(keccakSolver);
}

typedef struct tmqarg_t {
  uint64_t guessingBits;
  uint64_t thread_id;
  uint8_t *mqbuffer;
} tmqarg_t;

typedef struct tcheckarg_t {
  uint32_t *result_buffer;
  bool *preimage_found;
} tcheckarg_t;

__host__ void
threadUpdateMQSystem(void *arg) {

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

    uint64_t search_interval = keccakSolver->options.gbend - keccakSolver->options.gbstart;
    search_interval = search_interval > CHUNK_SIZE ? search_interval : CHUNK_SIZE;
    uint64_t gbstart = keccakSolver->options.gbstart;
    uint64_t gbend = gbstart + search_interval;

    size_t mbuffer_size = (CHUNK_SIZE * MQ_SYSTEM_SIZE) * sizeof(uint8_t);
    uint8_t *mqbuffer = (uint8_t *) malloc(mbuffer_size);
    size_t rbuffer_size = (CHUNK_SIZE * 25) * sizeof(uint32_t);
    uint32_t *result_buffer = (uint32_t *) malloc(rbuffer_size);
    bool preimage_found = false;

    for (uint64_t gb = gbstart; gb < gbend; gb += CHUNK_SIZE) {
        /* check output result */
        CUDA_CHECK(cudaMemcpy(result_buffer, keccakSolver->device_output_buffer, rbuffer_size, cudaMemcpyDeviceToHost));
        tcheckarg_t check_args[CHUNK_SIZE];
        for (uint64_t thread_id = 0; thread_id < CHUNK_SIZE; thread_id++) {
            check_args[thread_id].preimage_found = &preimage_found;
            check_args[thread_id].result_buffer = result_buffer + thread_id * 25;
            threadpool_add(keccakSolver->threadpool, threadCheckResult, (void *) &check_args[thread_id], 0);
        }
        threadpool_join(keccakSolver->threadpool, 0);

        if (preimage_found) {
            PRINTF_STAMP("one valid pre-image found!\n");
            break;
        }
        /* update MQ System */
        memset(mqbuffer, 0, CHUNK_SIZE * MQ_SYSTEM_SIZE);
        tmqarg_t args[CHUNK_SIZE];
        uint64_t this_gb;
        for (this_gb = gb; this_gb < gb + CHUNK_SIZE; this_gb++) {
            uint64_t thread_id = this_gb - gb;
            args[thread_id].thread_id = thread_id;
            args[thread_id].guessingBits = this_gb;
            args[thread_id].mqbuffer = mqbuffer + thread_id * MQ_SYSTEM_SIZE;
            threadpool_add(keccakSolver->threadpool, threadUpdateMQSystem, (void *) &args[thread_id], 0);
        }
        threadpool_join(keccakSolver->threadpool, 0);

        /* use gpu to solve mq system*/
        CUDA_CHECK(cudaDeviceSynchronize());
        CUDA_CHECK(cudaMemcpy(keccakSolver->device_mq_buffer, mqbuffer, mbuffer_size, cudaMemcpyHostToDevice));
        kernelLoop << < GPU_BLOCK_NUM, tpb >> > (keccakSolver->device_output_buffer, keccakSolver->device_mq_buffer);
    }
    free(mqbuffer);
}

__host__ void
keccakSolverExit(KeccakSolver *keccakSolver) {
    options_free(&keccakSolver->options);
    threadpool_destroy(keccakSolver->threadpool, 0);
    CUDA_CHECK(cudaDeviceSynchronize());
    CUDA_CHECK(cudaFree(keccakSolver->device_mq_buffer));
    CUDA_CHECK(cudaFree(keccakSolver->device_output_buffer));
    EXIT_WITH_MSG("keccak solver exit\n");
}