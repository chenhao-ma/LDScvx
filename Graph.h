//
// Created by MA Chenhao on 29/11/2021.
//

#ifndef LDSCVX_GRAPH_H
#define LDSCVX_GRAPH_H

#include <vector>
#include <stack>
#include <algorithm>
#include <queue>

#include "FlowNetwork.h"

using namespace std;

class Graph {
public:
    Graph(FILE *pFile, int NT, int topk);
    void findLDS();

    void output(char *ds_address);

private:
    int n;
    unsigned long m;
    int NT;
    int CT;
    int topk;
    int nsg;
    int max_d;
    int num_verify;
    bool check_first = false;
    bool fn_baseline = false;
    vector<vector<int>> adj;
    vector<bool> selected;
    vector<bool> active;
    vector<int> slt_nodes;
    vector<int> slt_edges;
    stack<vector<int>> stk_nodes;
    stack<vector<int>> stk_edges;
    stack<int> stk_CT;
    vector<vector<int>> ldses;
    vector<double> lds_rho;
    vector<int> lds_num;
    vector<int> deg;
    vector<int> pos;
    vector<int> sg;
    vector<pair<int, int>> edges;
    vector<double> alpha;
    vector<double> r;
    vector<double> rho_u;
    vector<double> rho_gu;
    vector<double> rho_l;
    vector<double> val;
    vector<int> nag;
    vector<int> fa;
    vector<pair<int, int>> cmpt;
    vector<int> veri_vtx;

    void frank_wolfe();
    void pava();
    void check_sg();
    void pruning();
    void compute_core();
    void prune_by_core();
    bool verify_LDS(vector<int> & nodes, double g);
    int find_fa(int x);
    void connected_components();
};


#endif //LDSCVX_GRAPH_H
