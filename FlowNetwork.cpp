//
// Created by MA Chenhao on 3/12/2021.
//

#include "FlowNetwork.h"

FlowNetwork::FlowNetwork(vector<pair<int, int>> edges, double g) {
    n = 0;
    for (auto & edge : edges) {
        if (mapping.find(edge.first) == mapping.end()) {
            ++n;
            mapping.insert(make_pair(edge.first, n));
        }
        if (mapping.find(edge.second) == mapping.end()) {
            ++n;
            mapping.insert(make_pair(edge.second, n));
        }
    }
    n += 2;
    adj.resize(n);

    for (auto &edge : edges) {
        int u = mapping.find(edge.second)->second;
        int v = mapping.find(edge.first)->second;
        add_edge(u, v, 1);
        add_edge(v, u, 1);
    }

    for (int u = 1; u < n - 1; u++) {
        add_edge(0, u, adj[u].size() / 4.0);
    }

    for (int u = 1; u < n - 1; u++) {
        add_edge(u, n - 1, g);
    }
}

void FlowNetwork::add_edge(int from, int to, double cap) {
    adj[from].push_back(EdgeFN(from, to, cap, 0, adj[to].size()));
    if (from == to) {
        adj[from].back().index++;
    }
    adj[to].push_back(EdgeFN(to, from, 0, 0, adj[from].size() - 1));
}

void FlowNetwork::enqueue (int v) {
    if (!active[v] && excess[v] > 0 && dist[v] < n) {
        active[v] = true;
        B[dist[v]].push_back(v);
        b = std::max(b, dist[v]);
    }
}

void FlowNetwork::push (EdgeFN &e) {
    double amt = std::min(excess[e.from], e.cap - e.flow);
    if (dist[e.from] == dist[e.to] + 1 && amt > 0) {
        e.flow += amt;
        adj[e.to][e.index].flow -= amt;
        excess[e.to] += amt;
        excess[e.from] -= amt;
        enqueue(e.to);
    }
}

void FlowNetwork::gap (int k) {
    for (int v = 0; v < n; v++) {
        if (dist[v] >= k) {
            count[dist[v]]--;
            dist[v] = std::max(dist[v], n);
            count[dist[v]]++;
            enqueue(v);
        }
    }
}

void FlowNetwork::relabel(int v) {
    count[dist[v]]--;
    dist[v] = n;
    for (auto e : adj[v]) if (e.cap - e.flow > 0) {
            dist[v] = std::min(dist[v], dist[e.to] + 1);
        }
    count[dist[v]]++;
    enqueue(v);
}

void FlowNetwork::discharge(int v) {
    for (auto &e : adj[v]) {
        if (excess[v] > 0) {
            push(e);
        } else {
            break;
        }
    }

    if (excess[v] > 0) {
        if (count[dist[v]] == 1) {
            gap(dist[v]);
        } else {
            relabel(v);
        }
    }
}

double FlowNetwork::get_maxflow(int s, int t, bool need_initial) {
    if (need_initial) {
        dist = std::vector<int>(n, 0);
        excess = std::vector<double>(n, 0);
        count = std::vector<int>(n + 1, 0);
        active = std::vector<bool>(n, false);
        B = std::vector<std::vector<int>>(n);
        b = 0;

        for (auto &e: adj[s]) {
            excess[s] += e.cap;
        }

        count[0] = n;
        enqueue(s);
        active[t] = true;
    }

    while (b >= 0) {
        if (!B[b].empty()) {
            int v = B[b].back();
            B[b].pop_back();
            active[v] = false;
            discharge(v);
        } else {
            b--;
        }
    }

    return excess[t];
}