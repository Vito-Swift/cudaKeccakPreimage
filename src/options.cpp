//
// Created by vitowu on 3/13/20.
//

#include "options.h"
#include "utils.h"

static void
print_usage(char *prg_name) {
    PRINTF_STAMP("\n"
                 "Usage: %s [OPTIONS]\n"
                 "\n"
                 "\n"
                 "\n", prg_name);
}

/* copy parsed option into a char* pointer, take at most 1024 chars */
static inline void
copy_opt(char **str, char *optarg) {
    if (NULL == ((*str) = strndup(optarg, 1024))) {
        PRINTF_ERR("[!] invalid input file\n");
    }
}

void
options_init(Options *options) {
#ifdef TEST_PRE
    char variate_filename[] = "testBench/test_variate_analysis.dat";
    char append_filename[] = "testBench/test_append_analysis.dat";
    char lin_filename[] = "testBench/test_lin_analysis.dat";
    char mq_filename[] = "testBench/test_mq_analysis.dat";
    char message_file[] = "testBench/testMessage.txt";
    copy_opt(&options->test_message_file, message_file);
#else
    char variate_filename[] = "variate_analysis.dat";
    char append_filename[] = "append_analysis.dat";
    char lin_filename[] = "lin_analysis.dat";
    char mq_filename[] = "mq_analysis.dat";
#endif
    copy_opt(&options->i_lin_analysis_file, variate_filename);
    copy_opt(&options->a_lin_analysis_file, append_filename);
    copy_opt(&options->c_lin_analysis_file, lin_filename);
    copy_opt(&options->mq_analysis_file, mq_filename);

    options->dev_id = 0;
    options->cpu_thread_num = 1;
    options->gbstart = 0;
    options->gbend = 0;
}

void
options_free(Options *options) {
    if (options->c_lin_analysis_file)
        free(options->c_lin_analysis_file);

    if (options->mq_analysis_file)
        free(options->mq_analysis_file);

    if (options->i_lin_analysis_file)
        free(options->i_lin_analysis_file);

    if (options->test_message_file)
        free(options->test_message_file);
}

/* for getopt_long */
#define OP_CLIN_FILE    1
#define OP_MQ_FILE      2
#define OP_ILIN_FILE    3
#define OP_ALIN_FILE    4
#define OP_DEV_ID       5
#define OP_THREAD_NUM   6
#define OP_GB_START     7
#define OP_GB_END       8
#define OP_TEST_MSG_F   9

static struct option keccak_long_opts[] = {
        {"clin_file", required_argument, 0, OP_CLIN_FILE},
        {"mq_file",   required_argument, 0, OP_MQ_FILE},
        {"ilin_file", required_argument, 0, OP_ILIN_FILE},
        {"dev_id",    required_argument, 0, OP_DEV_ID},
        {"t",         required_argument, 0, OP_THREAD_NUM},
        {"alin_file", required_argument, 0, OP_ALIN_FILE},
        {"gb_start",  optional_argument, 0, OP_GB_START},
        {"gb_end",    optional_argument, 0, OP_GB_END},
        {"tm_file",   optional_argument, 0, OP_TEST_MSG_F},
        {"help", 0,                      0, 'h'},
        {0,      0,                      0, 0}
};

void
options_parse(Options *options, int argc, char **argv) {
    int c, opt_idx;
    if (argc == 1) {
        print_usage(argv[0]);
        exit(0);
    }
    while (-1 != (c = getopt_long(argc, argv, "h", keccak_long_opts, &opt_idx))) {
        switch (c) {
            case 0:
                /* If opts option set a flag, don't do anything */
                if (keccak_long_opts[opt_idx].flag == 0) {
                    PRINTF_STAMP("\toption %s: %s\n", keccak_long_opts[opt_idx].name,
                                 optarg ? optarg : "null");
                }
                break;

            case 'h':
                print_usage(argv[0]);
                options_free(options);
                exit(0);

            case OP_CLIN_FILE:
                copy_opt(&options->c_lin_analysis_file, optarg);
                PRINTF_STAMP("\t\toption constant linear analysis file: %s\n", options->c_lin_analysis_file);
                break;

            case OP_ALIN_FILE:
                copy_opt(&options->a_lin_analysis_file, optarg);
                PRINTF_STAMP("\t\toption append linear analysis file: %s\n", options->a_lin_analysis_file);
                break;

            case OP_ILIN_FILE:
                copy_opt(&options->i_lin_analysis_file, optarg);
                PRINTF_STAMP("\t\toption iterative linear analysis file: %s\n", options->i_lin_analysis_file);
                break;

            case OP_MQ_FILE:
                copy_opt(&options->mq_analysis_file, optarg);
                PRINTF_STAMP("\t\toption mq analysis file: %s\n", options->mq_analysis_file);
                break;

            case OP_GB_START:
                options->gbstart = strtoul(optarg, NULL, 0);
                PRINTF_STAMP("\t\toption start guessing bits: 0x%lx\n", options->gbstart);
                break;

            case OP_GB_END:
                options->gbend = strtoul(optarg, NULL, 0);
                PRINTF_STAMP("\t\toption end gussing bits: 0x%lx\n", options->gbend);
                break;

            case OP_DEV_ID:
                options->dev_id = (uint32_t) atoi(optarg);
                PRINTF_STAMP("\t\toption dev id: %d\n", options->dev_id);
                break;

            case OP_THREAD_NUM:
                options->cpu_thread_num = (uint64_t) strtoul(optarg, NULL, 0);
                PRINTF_STAMP("\t\toptions cpu thread num: %ld\n", options->cpu_thread_num);
                break;

            case OP_TEST_MSG_F:
                copy_opt(&options->test_message_file, optarg);
                PRINTF_STAMP("\t\ttest mode, option test message file: %s\n", options->test_message_file);
                break;

            case '?':
                /* getopt_long already printed an error message */
                break;

            default:
                EXIT_WITH_MSG("[!] unknown error, exit...\n");
        }
    }
}