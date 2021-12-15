//
// Created by MA Chenhao on 6/12/2021.
//

#ifndef LDSCVX_ARGS_H
#define LDSCVX_ARGS_H

#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <getopt.h>

#define USAGE_TXT							   \
    "usage: \n"                                \
    "\t[-g input directed graph file]\n"       \
    "\t[-t repeated iterations]\n"             \
    "\t[-k topk]\n"

class Args {
public:
    char *address{};
    char ds_address[50];
    int NT = 100;
    int topk = 5;
    ~Args();
    void usage(char *msg, int exit_status);
    void parse_args(int argc, char *argv[]);
};


#endif //LDSCVX_ARGS_H
