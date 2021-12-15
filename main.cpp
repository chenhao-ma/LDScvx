#include <iostream>
#include "Args.h"
#include "Graph.h"

int main(int argc, char *argv[]) {
    setbuf(stdout, NULL);
    clock_t begin = clock();

    Args *args = new Args();
    args->parse_args(argc, argv);
    FILE* d_file = fopen(args->address, "r");
    Graph g = Graph(d_file, args->NT, args->topk);

    clock_t io_end = clock();

    g.findLDS();
    g.output(args->ds_address);

    clock_t end = clock();
    double io_secs = double(io_end - begin) / CLOCKS_PER_SEC;
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    printf("io time: %.4f, total time: %.4f\n", io_secs, elapsed_secs);
    return 0;
}
