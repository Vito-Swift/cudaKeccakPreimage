//
// Created by vitowu on 3/13/20.
//

#include "solver.h"

#define KECCAK_VNUM 800

#define LIN_CONST_VNUM 800
#define LIN_CONST_XVNUM 801
#define LIN_CONST_EQNUM 707
#define LIN_ITER_EQNUM 63
#define LIN_CONST_SYSTEM_SIZE ((LIN_CONST_XVNUM) * (LIN_CONST_EQNUM))
#define LIN_ITER_SYSTEM_SIZE ((LIN_CONST_XVNUM) * (LIN_ITER_EQNUM))

#define MQ_VAR_NUM 31
#define MQ_EQ_NUM 48
#define MQ_XVAR_NUM ((((MQ_VAR_NUM) * ((MQ_VAR_NUM) + 1)) / 2) + MQ_VAR_NUM + 1)
#define MQ_SYSTEM_SIZE ((MQ_EQ_NUM) * (MQ_XVAR_NUM))

#define BMQ_VAR_NUM 94
#define BMQ_EQ_NUM 48
#define BMQ_XVAR_NUM ((((BMQ_VAR_NUM) * ((BMQ_VAR_NUM) + 1)) / 2) + BMQ_VAR_NUM + 1)
#define BMQ_SYSTEM_SIZE ((BMQ_EQ_NUM) * (BMQ_XVAR_NUM))

#define GPU_THREADS_PER_BLOCK 64
#define GPU_BLOCK_NUM 64
#define GPU_TOTAL_THREADS ((GPU_THREADS_PER_BLOCK) * (GPU_BLOCK_NUM))
#define CHUNK_SIZE 0xFFFFFFF

__constant__ uint32_t file_mq_fvar_idx[BMQ_VAR_NUM];
__constant__ uint32_t file_reverse_mq_fvar_idx[KECCAK_VNUM];
__constant__ uint32_t mq_fvar_idx[MQ_VAR_NUM];
__constant__ uint32_t reverse_mq_fvar_idx[KECCAK_VNUM];
__constant__ uint8_t iterative_constraints[LIN_ITER_SYSTEM_SIZE];

__host__ void
loadSystemsFromFile(KeccakSolver *keccakSolver) {
    char buffer[BMQ_XVAR_NUM + 3];
    char ch;
    uint32_t i, j;

    // load mq system
    FILE *fmq = fopen(keccakSolver->options.mq_analysis_file, "r");
    uint32_t local_file_mq_system[BMQ_SYSTEM_SIZE];
    uint32_t local_file_mq_fvar_idx[BMQ_VAR_NUM];
    uint32_t local_reverse_mq_fvar_idx[KECCAK_VNUM];
    memset(local_file_mq_fvar_idx, 0, BMQ_VAR_NUM);
    memset(local_reverse_mq_fvar_idx, 0xFFFFFFFF, KECCAK_VNUM);

    if (fmq == NULL) {
        EXIT_WITH_MSG("[!] cannot open mq analysis file\n");
    } else {
        uint32_t left_idx;
        uint32_t right_idx;
        for (i = 0; i < BMQ_VAR_NUM; i++) {
            fgets(buffer, '\n', fmq);
            strtok(buffer, "\n");
            sscanf(buffer, "%d %d", &left_idx, &right_idx);
            local_file_mq_fvar_idx[left_idx] = right_idx;
            local_reverse_mq_fvar_idx[right_idx] = left_idx;
        }
        for (i = 0; i < BMQ_EQ_NUM; i++) {
            for (j = 0; j < BMQ_XVAR_NUM; j++) {
                ch = fgetc(fmq);
                local_file_mq_system[i * BMQ_XVAR_NUM + j] = (uint8_t) (ch - '0');
            }
            fgetc(fmq);
        }
    }
    fclose(fmq);
    PRINTF_STAMP("[+] read mq analysis from file\n");
    PRINTF_STAMP("\t\tequation_number: %d\n", BMQ_EQ_NUM);
    PRINTF_STAMP("\t\tvar_number: %d\n", BMQ_VAR_NUM);
    PRINTF_STAMP("\t\txvar_number: %d\n", BMQ_XVAR_NUM);

    PRINTF_STAMP("[+] copy mq system to device...\n");
    CUDA_CHECK(cudaMemcpyToSymbol(file_mq_fvar_idx, local_file_mq_fvar_idx, BMQ_VAR_NUM * sizeof(uint32_t)));
    CUDA_CHECK(cudaMemcpyToSymbol(file_reverse_mq_fvar_idx, local_reverse_mq_fvar_idx, KECCAK_VNUM * sizeof(uint32_t)));
    CUDA_CHECK(cudaMemcpy(keccakSolver->device_mq_buffer,
                          local_file_mq_system,
                          BMQ_SYSTEM_SIZE * sizeof(uint8_t),
                          cudaMemcpyHostToDevice));

    // load linear system
    FILE *flin = fopen(keccakSolver->options.c_lin_analysis_file, "r");
    uint8_t local_c_constr[LIN_CONST_SYSTEM_SIZE];
    if (flin == NULL) {
        EXIT_WITH_MSG("[!] cannot open constant linear constraint file\n");
    } else {
        for (i = 0; i < LIN_CONST_EQNUM; i++) {
            for (j = 0; j < LIN_CONST_XVNUM; j++) {
                ch = fgetc(flin);
                local_c_constr[i * LIN_CONST_XVNUM + j] = (uint8_t) (ch - '0');
            }
            fgetc(flin);
        }
    }
    fclose(flin);
    PRINTF_STAMP("[+] read constant linear constraints from file\n");
    PRINTF_STAMP("\t\tequation_number: %d\n", LIN_CONST_EQNUM);
    PRINTF_STAMP("\t\tvar_number: %d\n", LIN_CONST_VNUM);
    PRINTF_STAMP("\t\txvar_number: %d\n", LIN_CONST_XVNUM);

    PRINTF_STAMP("[+] copy linear system to device...\n");
    CUDA_CHECK(cudaMemcpy(keccakSolver->device_c_constr_buffer,
                          local_c_constr,
                          LIN_CONST_SYSTEM_SIZE * sizeof(uint8_t),
                          cudaMemcpyHostToDevice));

    FILE *fi = fopen(keccakSolver->options.i_lin_analysis_file, "r");
    uint8_t local_i_constr[LIN_ITER_SYSTEM_SIZE];
    if (fi == NULL) {
        EXIT_WITH_MSG("[!] cannot open iterative linear constraint file\n");
    } else {
        for (i = 0; i < LIN_ITER_EQNUM; i++) {
            for (j = 0; j < LIN_CONST_XVNUM; j++) {
                ch = fgetc(fi);
                local_i_constr[i * LIN_CONST_XVNUM + j] = (uint8_t) (ch - '0');
            }
            fgetc(fi);
        }
    }
    fclose(fi);
    PRINTF_STAMP("[+] read iterative linear constraints from file\n");
    PRINTF_STAMP("\t\tequation_number: %d\n", LIN_ITER_EQNUM);
    PRINTF_STAMP("\t\tvar_number: %d\n", LIN_CONST_VNUM);
    PRINTF_STAMP("\t\txvar_number: %d\n", LIN_CONST_XVNUM);

    PRINTF_STAMP("[+] copy iterative linear system to device...\n");
    CUDA_CHECK(cudaMemcpyToSymbol(iterative_constraints, local_i_constr, LIN_ITER_SYSTEM_SIZE * sizeof(uint8_t)));

    CUDA_CHECK(cudaDeviceSynchronize());
}

__global__ void
kernelLoop(uint32_t *kern_output_buffer,
           uint8_t *linear_system_buffer,
           uint8_t *constance_constrs,
           uint8_t *bmq_system,
           uint64_t gbchunk_start,
           uint64_t gbchunk_end) {
    uint64_t thread_gbstart;
    uint64_t thread_gbend;
    uint64_t thread_id = (threadIdx.x + blockIdx.x * blockDim.x);
    uint64_t thread_interval = (gbchunk_end - gbchunk_start) / GPU_TOTAL_THREADS;

    thread_gbstart = thread_id * thread_interval;
    thread_gbend = (thread_id + 1) * thread_interval;
    uint64_t gb;

    uint8_t *kern_linear_system_buffer = linear_system_buffer +
        (LIN_ITER_SYSTEM_SIZE + LIN_CONST_SYSTEM_SIZE) * thread_id;
    // main loop
    for (gb = thread_gbstart; gb < thread_gbend; gb++) {
        memcpy(kern_linear_system_buffer, linear_system_buffer, LIN_CONST_SYSTEM_SIZE * sizeof(uint8_t));
        memcpy(kern_linear_system_buffer + LIN_CONST_SYSTEM_SIZE,
               iterative_constraints,
               LIN_ITER_SYSTEM_SIZE * sizeof(uint8_t));
        if (thread_id == 1)
            printf("hi\n");
    }
}

__host__ void
keccakSolverInit(KeccakSolver *keccakSolver, int argc, char **argv) {
    options_init(&keccakSolver->options);
    options_parse(&keccakSolver->options, argc, argv);

    CUDA_CHECK(cudaSetDevice(keccakSolver->options.dev_id));
    CUDA_CHECK(cudaDeviceSynchronize());
    CUDA_CHECK(cudaDeviceReset());
    CUDA_CHECK(cudaMalloc(&keccakSolver->device_output_buffer, 25 * GPU_TOTAL_THREADS * sizeof(uint32_t)));
    CUDA_CHECK(cudaMalloc(&keccakSolver->device_mq_buffer, BMQ_SYSTEM_SIZE * sizeof(uint8_t)));
    CUDA_CHECK(cudaMalloc(&keccakSolver->device_c_constr_buffer, LIN_CONST_SYSTEM_SIZE * sizeof(uint8_t)));
    CUDA_CHECK(cudaMalloc(&keccakSolver->device_linsys_buffer,
                          (LIN_CONST_SYSTEM_SIZE + LIN_ITER_SYSTEM_SIZE) * GPU_TOTAL_THREADS * sizeof(uint8_t)));
    loadSystemsFromFile(keccakSolver);
}

__host__ void
keccakSolverLoop(KeccakSolver *keccakSolver) {
    CUDA_CHECK(cudaDeviceSynchronize());
    dim3 tpb(GPU_THREADS_PER_BLOCK);

    uint64_t search_interval = keccakSolver->options.gbend - keccakSolver->options.gbstart;
    if (search_interval > CHUNK_SIZE) {
        for (uint64_t i = keccakSolver->options.gbstart; i < keccakSolver->options.gbend; i += CHUNK_SIZE) {
            kernelLoop << < GPU_BLOCK_NUM, tpb >> >
                (keccakSolver->device_output_buffer,
                    keccakSolver->device_linsys_buffer,
                    keccakSolver->device_c_constr_buffer,
                    keccakSolver->device_mq_buffer,
                    i, i + CHUNK_SIZE);
        }
    } else {
        kernelLoop << < GPU_BLOCK_NUM, tpb >> >
            (keccakSolver->device_output_buffer,
                keccakSolver->device_linsys_buffer,
                keccakSolver->device_c_constr_buffer,
                keccakSolver->device_mq_buffer,
                keccakSolver->options.gbstart,
                keccakSolver->options.gbend);
    }
    CUDA_CHECK(cudaDeviceSynchronize());
}

__host__ void
keccakSolverExit(KeccakSolver *keccakSolver) {
    options_free(&keccakSolver->options);
    CUDA_CHECK(cudaDeviceSynchronize());
    CUDA_CHECK(cudaFree(keccakSolver->device_output_buffer));
    CUDA_CHECK(cudaFree(keccakSolver->device_mq_buffer));
    CUDA_CHECK(cudaFree(keccakSolver->device_c_constr_buffer));
    CUDA_CHECK(cudaFree(keccakSolver->device_linsys_buffer));
    EXIT_WITH_MSG("keccak solver exit\n");
}