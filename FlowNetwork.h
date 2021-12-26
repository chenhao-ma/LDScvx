//
// Created by MA Chenhao on 3/12/2021.
//

#ifndef LDSCVX_FLOWNETWORK_H
#define LDSCVX_FLOWNETWORK_H

#include <vector>
#include <queue>
#include "EdgeFN.h"
#include <cmath>
#include <map>

using namespace std;


class FlowNetwork {
public:
    int n;
    vector <vector <EdgeFN>> adj;
    vector <double> excess;
    vector <int> dist, count;
    vector <bool> active;
    vector <vector <int>> B;
    vector<int> nums;
    map<int, int> mapping;
    std::vector<int> ori_id;
    int b;
    double m;

    FlowNetwork(vector<pair<int, int>> edges, double g, bool LDSvalidate = false);

    void add_edge(int from, int to, double cap);

    void enqueue (int v);

    void push (EdgeFN &e);

    void gap (int k);

    void relabel (int v);

    void discharge (int v);

    double get_maxflow(int s, int t, bool need_initial = true);

    double get_mincut(int s, int t, std::vector<int> &S, bool need_initial = true);


};



#endif //LDSCVX_FLOWNETWORK_H
