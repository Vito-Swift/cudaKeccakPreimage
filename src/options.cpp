//
// Created by vitowu on 3/13/20.
//

#include "options.h"

static void
print_usage(char *prg_name) {
    printf("\n"
           "Usage: %s [OPTIONS]\n"
           "\n"
           "\n"
           "\n", prg_name);
}

void
options_init(Options *options) {
    options->c_lin_analysis_file = nullptr;
    options->i_lin_analysis_file = nullptr;
    options->mq_analysis_file = nullptr;
    options->dev_id = 0;
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
}

/* for getopt_long */
#define OP_CLIN_FILE    1
#define OP_MQ_FILE      2
#define OP_ILIN_FILE    3
#define OP_DEV_ID       4
#define OP_GB_START     5
#define OP_GB_END       6

static struct option keccak_long_opts[] = {
    {"clin_file", required_argument, 0, OP_CLIN_FILE},
    {"mq_file", required_argument, 0, OP_MQ_FILE},
    {"ilin_file", required_argument, 0, OP_ILIN_FILE},
    {"dev_id", required_argument, 0, OP_DEV_ID},
    {"gb_start", required_argument, 0, OP_GB_START},
    {"gb_end", required_argument, 0, OP_GB_END},
    {"help", 0, 0, 'h'},
    {0, 0, 0, 0}
};

/* copy parsed option into a char* pointer, take at most 1024 chars */
static inline void
copy_opt(char **str, char *optarg) {
    if (NULL == ((*str) = strndup(optarg, 1024))) {
        fprintf(stderr, "[!] invalid input file\n");
    }
}

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
                    printf("option %s: %s\n", keccak_long_opts[opt_idx].name,
                           optarg ? optarg : "null");
                }
                break;

            case 'h':print_usage(argv[0]);
                options_free(options);
                exit(0);

            case OP_CLIN_FILE:copy_opt(&options->c_lin_analysis_file, optarg);
                printf("option constant linear analysis file: %s\n", options->c_lin_analysis_file);
                break;

            case OP_ILIN_FILE:copy_opt(&options->i_lin_analysis_file, optarg);
                printf("option iterative linear analysis file: %s\n", options->i_lin_analysis_file);
                break;

            case OP_MQ_FILE:copy_opt(&options->mq_analysis_file, optarg);
                printf("option mq analysis file: %s\n", options->mq_analysis_file);
                break;

            case OP_GB_START: options->gbstart = strtoul(optarg, NULL, 0);
                printf("option start guessing bits: 0x%lx\n", options->gbstart);
                break;

            case OP_GB_END: options->gbend = strtoul(optarg, NULL, 0);
                printf("option end gussing bits: 0x%lx\n", options->gbend);
                break;

            case OP_DEV_ID:options->dev_id = (uint32_t) atoi(optarg);
                printf("option dev id: %d\n", options->dev_id);
                break;

            case '?':
                /* getopt_long already printed an error message */
                break;

            default:fprintf(stderr, "[!] unknown error, exit...\n");
        }
    }
}