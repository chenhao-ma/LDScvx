//
// Created by MA Chenhao on 6/12/2021.
//

#include "Args.h"

Args::~Args() {
    delete address;
}

void Args::usage(char *msg, int exit_status) {
    fprintf(exit_status == 0 ? stdout : stderr, "%s", USAGE_TXT);

    if (msg) {
        fprintf(exit_status == 0 ? stdout : stderr, "\n%s\n", msg);
    }
    exit(exit_status);
}

void Args::parse_args(int argc, char **argv) {
    int c;
    opterr = 0;
    if (argc < 2) {
        usage(nullptr, 0);
    }

    while ((c = getopt(argc, argv, "g:t:k:")) != -1) {
        switch (c) {
            case 'g': {
                printf("%s\n", optarg);
                address = strdup(optarg);
                break;
            }

            case 't': {
                sscanf(optarg, "%d", &NT);
                break;
            }

            case 'k': {
                sscanf(optarg, "%d", &topk);
            }

            default:
                break;
        }
    }
    char dataset[5] = "";
    strncat(dataset, address, 2);
    printf("%s\n", dataset);
    sprintf(ds_address, "output/%s-%d.txt", dataset, topk);

    printf("%s\n", ds_address);
}
