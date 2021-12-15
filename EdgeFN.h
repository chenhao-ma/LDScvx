//
// Created by MA Chenhao on 3/12/2021.
//

#ifndef LDSCVX_EDGEFN_H
#define LDSCVX_EDGEFN_H


class EdgeFN {
public:
    int from, to, index;
    double cap, flow;
    EdgeFN(int from, int to, double cap, double flow, double index);
};


#endif //LDSCVX_EDGEFN_H
