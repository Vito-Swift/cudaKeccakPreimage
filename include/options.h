//
// Created by vitowu on 3/13/20.
//

#ifndef KECCAKSOLVER_OPTIONS_H
#define KECCAKSOLVER_OPTIONS_H

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <time.h>
#include <stdbool.h>

typedef struct {
  char *c_lin_analysis_file;
  char *i_lin_analysis_file;
  char *a_lin_analysis_file;
  char *mq_analysis_file;
  uint64_t cpu_thread_num;
  uint32_t dev_id;
  uint64_t gbstart;
  uint64_t gbend;
} Options;

void
options_init(Options *options);

void
options_free(Options *options);

void
options_parse(Options *options, int argc, char **argv);

#endif //KECCAKSOLVER_OPTIONS_H
