//
// Created by vitowu on 3/18/20.
//

#ifndef KECCAKSOLVER_PARAMS_H
#define KECCAKSOLVER_PARAMS_H

#define LIN_CONST_EQNUM 707
#define BMQ_XVAR_NUM 321201


#ifdef TEST_PRE
#define CHUNK_SIZE (0x100UL)
#define MQ_VAR_NUM 30UL
#define MQ_EQ_NUM 38UL
#define LIN_ITER_EQNUM 54UL
#define IMQ_VAR_NUM 94UL
#define AMQ_VAR_NUM 40UL
#define AMQ_LIN_EQNUM 10UL
#else
#define CHUNK_SIZE (0x20000UL)
#define MQ_VAR_NUM 31UL
#define MQ_EQ_NUM 38UL
#define LIN_ITER_EQNUM 53UL
#define IMQ_VAR_NUM 94UL
#define AMQ_VAR_NUM 41UL
#define AMQ_LIN_EQNUM 10UL
#endif

#define LIN_CONST_SYSTEM_SIZE (801 * (LIN_CONST_EQNUM))
#define LIN_ITER_SYSTEM_SIZE (801* (LIN_ITER_EQNUM))
#define IMQ_XVAR_NUM ((((IMQ_VAR_NUM) * ((IMQ_VAR_NUM) + 1)) / 2) + IMQ_VAR_NUM + 1)
#define AMQ_XVAR_NUM ((((AMQ_VAR_NUM) * ((AMQ_VAR_NUM) + 1)) / 2) + AMQ_VAR_NUM + 1)
#define MQ_XVAR_NUM ((((MQ_VAR_NUM) * ((MQ_VAR_NUM) + 1)) / 2) + MQ_VAR_NUM + 1)
#define MQ_SYSTEM_SIZE ((MQ_EQ_NUM) * (MQ_XVAR_NUM))

#define GPU_THREADS_PER_BLOCK (32)
#define GPU_BLOCK_NUM ((CHUNK_SIZE) / (GPU_THREADS_PER_BLOCK))

#endif //KECCAKSOLVER_PARAMS_H
