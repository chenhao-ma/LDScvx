//
// Created by MA Chenhao on 29/11/2021.
//
#include "Graph.h"

#define DELTA 1e-5

Graph::Graph(FILE *file, int NT, int topk) {
    printf("Graph construction\n");
    fscanf(file, "%d%lu", &n, &m);

    adj.resize(n);
    r.resize(n);
    deg.resize(n);
    selected.resize(n, true);
    active.resize(n, true);
    lds_num.resize(n, -1);
    num_verify = 0;
    veri_vtx.resize(n, 0);
    rho_gu.resize(n);
    rho_u.resize(n);
    rho_l.resize(n);
    val.resize(n);
    nag.resize(n);
    sg.resize(n);
    fa.resize(n);
    slt_nodes.clear();
    slt_edges.clear();

//    printf("ok1\n");
    for (int i = 0; i < m; i++) {
        int u, v;
        fscanf(file, "%d%d", &u, &v);
        edges.emplace_back(u, v);
        adj[u].push_back(v);
        adj[v].push_back(u);
    }
    alpha.resize(m);
//    printf("ok2\n");

    this->NT = NT;
    this->topk = topk;
}

void Graph::frank_wolfe() {
    printf("#node %lu #edge %lu\n", slt_nodes.size(), slt_edges.size());
    if (CT == 0) {
        for (auto &edge: slt_edges) {
            alpha[edge] = 0.5;
        }
    }

    for (auto & node : slt_nodes) {
        r[node] = 0;
    }
    for (auto & e : slt_edges) {
        r[edges[e].first] += alpha[e];
        r[edges[e].second] += 1 - alpha[e];
    }


    for (int t = CT + 1; t <= CT + NT; t++) {
        double gamma_t = 2.0 / (t + 2);
        for (auto & e : slt_edges) {
//            alpha[e] *= 1 - gamma_t;
            if (r[edges[e].first] < r[edges[e].second]) {
                r[edges[e].first] += gamma_t * (1 - alpha[e]);
                r[edges[e].second] -= gamma_t * (1 - alpha[e]);
                alpha[e] = alpha[e] * (1 - gamma_t) + gamma_t;
            }
            else if (r[edges[e].first] > r[edges[e].second]) {
                r[edges[e].first] -= gamma_t * alpha[e];
                r[edges[e].second] += gamma_t * alpha[e];
                alpha[e] = alpha[e] * (1 - gamma_t);
            }
        }
//        for (auto & node : slt_nodes) {
//            r[node] = 0;
//        }
//        for (auto & e : slt_edges) {
//            r[edges[e].first] += alpha[e];
//            r[edges[e].second] += 1 - alpha[e];
//        }
    }
//    for (int i = 0; i < n; i++) {
//        printf("%.4f ", r[i]);
//    }
//    printf("\n");
}

/* Pool Adjacent Violators Algorithm. Try to decompose the graph into Stable Groups.
 * deg reused here!
 */
void Graph::pava() {
    sort(slt_nodes.begin(), slt_nodes.end(), [this](int a, int b)->bool {
        return r[a] > r[b];
    });
    for (int i = 0; i < slt_nodes.size(); i++) {
//        if (i < 10) printf("%d %.4f\n", slt_nodes[i], r[slt_nodes[i]]);
        deg[slt_nodes[i]] = i;
    }
    vector<int> ne(slt_nodes.size(), 0);
    for (auto e : slt_edges) {
        int u = deg[edges[e].first];
        int v = deg[edges[e].second];
        ne[(u > v)? u:v]++;
    }
//    for (int i = 0; i < n; i++) {
//        printf("%d ", ne[i]);
//    }
//    printf("\n");

    nag[0] = 1;
    val[0] = ne[0];
    nsg = 0;
    for (int i = 1; i < slt_nodes.size(); i++) {
        nsg += 1;
        val[nsg] = ne[i];
        nag[nsg] = 1;
        while ((nsg > 0) && (val[nsg] > val[nsg - 1] - 1e-5)) {
            val[nsg - 1] = (nag[nsg] * val[nsg] + nag[nsg - 1] * val[nsg - 1])/(nag[nsg]+nag[nsg - 1]);
            nag[nsg - 1] += nag[nsg];
            nsg--;
        }
    }
    nsg++;

    printf("nsg %d\n", nsg);

    int cur = 0;
    for (int i = 0; i < nsg; i++) {
        for (int j = cur; j < cur + nag[i]; j++) {
            sg[slt_nodes[j]] = i;
        }
        cur += nag[i];
    }
}

void Graph::check_sg() {
    if (nsg <= 1) return;
    sort(slt_edges.begin(), slt_edges.end(), [this](int a, int b)->bool {
        int a1 = min(sg[edges[a].first], sg[edges[a].second]);
        int a2 = max(sg[edges[a].first], sg[edges[a].second]);
        int b1 = min(sg[edges[b].first], sg[edges[b].second]);
        int b2 = max(sg[edges[b].first], sg[edges[b].second]);
        if (a1 != b1)
            return a1 < b1;
        else
            return a2 < b2;
    });
    vector<bool> valid(nsg);
    vector<int> bin(nsg + 1, 0);
    for (int i = 0; i < slt_edges.size(); i++) {
        int a = min(sg[edges[slt_edges[i]].first], sg[edges[slt_edges[i]].second]);
        ++bin[a];
    }
    int s = 0;
    for (int i = 0; i <= nsg; i++) {
        int tmp = s;
        s += bin[i];
        bin[i] = tmp;
    }

    vector<double> max_r(nsg);
    int cur = 0;
    for (int i = 0; i < nsg; i++) {
        max_r[i] = r[slt_nodes[cur]];
//        double tmp_min = r[slt_nodes[cur]];
        for (int j = cur + 1; j < cur + nag[i]; j++) {
            max_r[i] = max(max_r[i], r[slt_nodes[j]]);
//            tmp_min = min(tmp_min, r[slt_nodes[j]]);
        }
//        printf("%d %d %.4f %.4f\n", i, cur, max_r[i], tmp_min);
        cur += nag[i];
    }

    double min_r = r[slt_nodes[0]];
    cur = 0;
    vector<int> cpt(bin);
    clock_t start = clock();
    for (int i = 0; i < nsg - 1; i++) {
//        printf("checking %d\n", i);
        for (int j = cur; j < cur + nag[i]; j++) {
            min_r = min(min_r, r[slt_nodes[j]]);
        }
        double min_t = min_r;
        vector<double> max_tmp(max_r.begin() + i + 1, max_r.end());
        double max_t = *max_element(begin(max_tmp), end(max_tmp));
//        printf("%.4f %.4f\n", min_t, max_t);
        queue<tuple<int, int, double>> q;
        valid[i] = true;
//        if (i % 3 != 0) {
//            valid[i] = false;
//            continue;
//        }
        for (int j = 0; j <= i; j++) {
            while (cpt[j] < bin[j + 1]) {
                int e = slt_edges[cpt[j]];
                if (max(sg[edges[e].first], sg[edges[e].second]) > i) {
                    q.push(make_tuple(j, cpt[j], alpha[e]));
                    int u = edges[e].first;
                    int v = edges[e].second;
                    if (sg[u] > sg[v]) {
                        r[u] += 1 - alpha[e];
                        r[v] -= 1 - alpha[e];
                        max_r[sg[u]] = max(max_r[sg[u]], r[u]);
                        max_t = max(max_t, r[u]);
                        min_t = min(min_t, r[v]);
                        alpha[e] = 1;
                    } else {
                        r[u] -= alpha[e];
                        r[v] += alpha[e];
                        min_t = min(min_t, r[u]);
                        max_r[sg[v]] = max(max_r[sg[v]], r[v]);
                        max_t = max(max_t, r[v]);
                        alpha[e] = 0;
                    }
                    if (min_t <= max_t) {
                        valid[i] = false;
                        break;
                    }
                }
                ++cpt[j];
            }
            if (!valid[i]) break;
        }
        if (valid[i]) {
            min_r = min(min_t, min_r);
        } else {
            for (int j = i + 1; j < max_r.size(); j++) {
                max_r[j] = max_tmp[j - i - 1];
            }
            while (!q.empty()) {
                auto tp = q.front(); q.pop();
                cpt[get<0>(tp)] = min(cpt[get<0>(tp)], get<1>(tp));
                int e = slt_edges[get<1>(tp)];
                int u = edges[e].first, v = edges[e].second;
                r[u] += get<2>(tp) - alpha[e];
                r[v] -= get<2>(tp) - alpha[e];
                alpha[e] = get<2>(tp);
            }
        }
        clock_t cur_t = clock();
//        printf("used %.4f\n", (double) (cur_t - start) / CLOCKS_PER_SEC);
    }

    check_first = (nsg == 1) || valid[0];
    //merge stable groups
    if (nsg > 1) {
        vector<int> n_nag;
        for (int i = 0; i < nsg; i++) {
            if (i == 0 || valid[i - 1]) {
                n_nag.push_back(nag[i]);
            } else {
                n_nag[n_nag.size() - 1] += nag[i];
            }
        }
        nag = n_nag;

        cur = 0;
        nsg = nag.size();
        double minr = m;
        for (int i = 0; i < nsg; i++) {
//            printf("nsg %d %d ", i, nag[i]);
            double tmp_r = m;
            for (int j = cur; j < cur + nag[i]; j++) {
                sg[slt_nodes[j]] = i;
                tmp_r = min(tmp_r, r[slt_nodes[j]]);
                minr = min(minr, r[slt_nodes[j]]);
            }
//            printf("%.4f %.4f\n", minr, tmp_r);
            cur += nag[i];
        }
    }
    printf("updated nsg %d\n", nsg);

}

//void Graph::check_sg() {
//    vector<bool> valid(nsg);
//    queue<pair<int, double>> q;
//    int cur = 0;
//    for (int i = 0; i < nsg - 1; i++){
//        cur += nag[i];
//        for (auto e : slt_edges) {
//            int u = edges[e].first, v = edges[e].second;
//            if (min(sg[u], sg[v]) <= i && max(sg[u], sg[v]) > i) {
//                q.push(make_pair(e, alpha[e]));
//                if (sg[u] > sg[v]) {
//                    r[u] += 1 - alpha[e];
//                    r[v] -= 1 - alpha[e];
//                    alpha[e] = 1;
//                } else if (sg[u] < sg[v]) {
//                    r[u] -= alpha[e];
//                    r[v] += alpha[e];
//                    alpha[e] = 0;
//                }
//            }
//        }
//        double min_r = r[slt_nodes[0]];
//        for (int j = 0; j < cur; j++) {
//            min_r = min(min_r, r[slt_nodes[j]]);
//        }
//        double max_r = r[slt_nodes[cur]];
//        for (int j = cur; j < slt_nodes.size(); j++) {
//            max_r = max(max_r, r[slt_nodes[j]]);
//        }
//        printf("%d %.4f %.4f\n", i, min_r, max_r);
//        if (min_r > max_r) {
//            valid[i] = true;
//        } else {
//            valid[i] = false;
//            while (!q.empty()) {
//                auto tmp = q.front(); q.pop();
//                int u = edges[tmp.first].first, v = edges[tmp.first].second;
//                r[u] += tmp.second - alpha[tmp.first];
//                r[v] -= tmp.second - alpha[tmp.first];
//                alpha[tmp.first] = tmp.second;
//            }
//        }
//    }
//
//    check_first = (nsg == 1) || valid[0];
//    //merge stable groups
//    if (nsg > 1) {
//        vector<int> n_nag;
//        for (int i = 0; i < nsg; i++) {
//            if (i == 0 || valid[i - 1]) {
//                n_nag.push_back(nag[i]);
//            } else {
//                n_nag[n_nag.size() - 1] += nag[i];
//            }
//        }
//        nag = n_nag;
//
//        cur = 0;
//        nsg = nag.size();
//        double minr = m;
//        for (int i = 0; i < nsg; i++) {
//            printf("nsg %d %d ", i, nag[i]);
//            double tmp_r = m;
//            for (int j = cur; j < cur + nag[i]; j++) {
//                sg[slt_nodes[j]] = i;
//                tmp_r = min(tmp_r, r[slt_nodes[j]]);
//                minr = min(minr, r[slt_nodes[j]]);
//            }
//            printf("%.4f %.4f\n", minr, tmp_r);
//            cur += nag[i];
//        }
//    }
//    printf("updated nsg %d\n", nsg);
//
//}

//void Graph::check_sg() {
//    sort(slt_edges.begin(), slt_edges.end(), [this](int a, int b)->bool {
//        int a1 = min(sg[edges[a].first], sg[edges[a].second]);
//        int a2 = max(sg[edges[a].first], sg[edges[a].second]);
//        int b1 = min(sg[edges[b].first], sg[edges[b].second]);
//        int b2 = max(sg[edges[b].first], sg[edges[b].second]);
//        if (a1 != b1)
//            return a1 < b1;
//        else
//            return a2 < b2;
//    });
//
//    vector<vector<int>::iterator> its;
//    auto it = slt_edges.begin();
//
//    vector<double> max_r(nsg);
//    int cur = 0;
//    for (int i = 0; i < nsg; i++) {
//        its.push_back(it);
//        while (it != slt_edges.end() && min(sg[edges[*it].first], sg[edges[*it].second]) <= i) it++;
//
//        max_r[i] = r[slt_nodes[cur]];
//        double tmp_min = r[slt_nodes[cur]];
//        for (int j = cur + 1; j < cur + nag[i]; j++) {
//            max_r[i] = max(max_r[i], r[slt_nodes[j]]);
//            tmp_min = min(tmp_min, r[slt_nodes[j]]);
//        }
////        printf("%d %d %.4f %.4f\n", i, cur, max_r[i], tmp_min);
//        cur += nag[i];
//    }
//
//    vector<bool> valid(nsg);
//
//    double min_r = r[slt_nodes[0]];
//    cur = 0;
//
//    for (int i = 0; i < nsg - 1; i++) {
////        printf("i %d\n", i);
//        printf("%d %d %.4f %.4f\n", i, cur, min_r, max_r[i + 1]);
//        for (int j = cur; j < cur + nag[i]; j++) {
//            min_r = min(min_r, r[slt_nodes[j]]);
//        }
//        printf("%d %d %.4f %.4f\n", i, cur, min_r, max_r[i + 1]);
//        cur += nag[i];
//        vector<double> tmp_max_r(max_r.begin() + i + 1, max_r.end());
//        double mmax = *max_element(begin(tmp_max_r), end(tmp_max_r));
//        double mmin = min_r;
//        stack<int> s_alpha;
//        stack<int> s_pos;
//        valid[i] = true;
//        for (int j = 0; j <= i; j++) {
//            while (its[j] != slt_edges.end() && min(sg[edges[*its[j]].first], sg[edges[*its[j]].second]) <= j) {
//                if (max(sg[edges[*its[j]].first], sg[edges[*its[j]].second]) > i) {
//                    s_alpha.push(alpha[*its[j]]);
//                    s_pos.push(j);
//                    int u = edges[*its[j]].first, v = edges[*its[j]].second;
//                    if (sg[u] > sg[v]) {
//                        r[u] += 1 - alpha[*its[j]];
//                        max_r[sg[u]] = max(max_r[sg[u]], r[u]);
//                        mmax = max(mmax, r[u]);
//                        r[v] -= 1 - alpha[*its[j]];
//                        mmin = min(mmin, r[v]);
//                        alpha[*its[j]] = 1;
//                    } else if (sg[u] < sg[v]) {
//                        r[u] -= alpha[*its[j]];
//                        mmin = min(mmin, r[u]);
//                        r[v] += alpha[*its[j]];
//                        max_r[sg[v]] = max(max_r[sg[v]], r[v]);
//                        mmax = max(mmax, r[v]);
//                        alpha[*its[j]] = 0;
//                    }
//
//                }
//                its[j]++;
//                if (mmin <= mmax) {
//                    valid[i] = false;
//                    break;
//                }
//            }
//            if (!valid[i]) break;
//        }
//        printf("%d ", i);
//        printf(valid[i] ? "true\n" : "false\n");
//        if (valid[i]) {
//            min_r = min(min_r, mmin);
//        } else {
//            for (int j = i + 1; j < max_r.size(); j++) {
//                max_r[j] = tmp_max_r[j - i - 1];
//            }
//            while (!s_alpha.empty()) {
//                double val = s_alpha.top(); s_alpha.pop();
//                auto j = s_pos.top(); s_pos.pop();
//                its[j]--;
//                int u = edges[*its[j]].first, v = edges[*its[j]].second;
//                r[u] += val - alpha[*its[j]];
//                r[v] -= val - alpha[*its[j]];
//                alpha[*its[j]] = val;
//            }
//        }
//        printf("%d %.4f %.4f %.4f %.4f\n", i, min_r, mmin, mmax, max_r[i + 1]);
//    }
//
//    check_first = (nsg == 1) || valid[0];
//    //merge stable groups
//    if (nsg > 1) {
//        vector<int> n_nag;
//        for (int i = 0; i < nsg; i++) {
//            if (i == 0 || valid[i - 1]) {
//                n_nag.push_back(nag[i]);
//            } else {
//                n_nag[n_nag.size() - 1] += nag[i];
//            }
//        }
//        nag = n_nag;
//
//        cur = 0;
//        nsg = nag.size();
//        double minr = m;
//        for (int i = 0; i < nsg; i++) {
//            printf("nsg %d %d ", i, nag[i]);
//            double tmp_r = m;
//            for (int j = cur; j < cur + nag[i]; j++) {
//                sg[slt_nodes[j]] = i;
//                tmp_r = min(tmp_r, r[slt_nodes[j]]);
//                minr = min(minr, r[slt_nodes[j]]);
//            }
//            printf("%.4f %.4f\n", minr, tmp_r);
//            cur += nag[i];
//        }
//    }
//    printf("updated nsg %d\n", nsg);
//}

//void Graph::check_sg() {
//    vector<double> min_r(nsg);
//    vector<double> max_r(nsg);
//    int cur = 0;
////    for (int i = 0; i < nsg; i++) {
////        min_r[i] = r[slt_nodes[cur]];
////        max_r[i] = r[slt_nodes[cur]];
////        for (int j = cur + 1; j < cur + nag[i]; j++) {
////            min_r[i] = min(min_r[i], r[slt_nodes[j]]);
////            max_r[i] = max(max_r[i], r[slt_nodes[j]]);
////        }
////        cur += nag[i];
////        printf("i %d min %.4f max %.4f\n", i, min_r[i], max_r[i]);
////        if (i > 0) min_r[i] = min(min_r[i], min_r[i - 1]);
////    }
//
//    for (auto e : slt_edges) {
//        int u = edges[e].first;
//        int v = edges[e].second;
//        if (sg[u] > sg[v]) {
//            r[u] += 1 - alpha[e];
//            r[v] -= 1 - alpha[e];
//            alpha[e] = 1;
//        } else if (sg[u] < sg[v]) {
//            r[u] -= alpha[e];
//            r[v] += alpha[e];
//            alpha[e] = 0;
//        }
//    }
//
////    vector<double> min_r(nsg);
////    vector<double> max_r(nsg);
//    cur = 0;
//    for (int i = 0; i < nsg; i++) {
//        min_r[i] = r[slt_nodes[cur]];
//        max_r[i] = r[slt_nodes[cur]];
//        for (int j = cur + 1; j < cur + nag[i]; j++) {
//            min_r[i] = min(min_r[i], r[slt_nodes[j]]);
//            max_r[i] = max(max_r[i], r[slt_nodes[j]]);
//        }
//        cur += nag[i];
//        printf("i %d min %.4f max %.4f\n", i, min_r[i], max_r[i]);
//        if (i > 0) min_r[i] = min(min_r[i], min_r[i - 1]);
//    }
////    for (int i = 0; i < n; i++) {
////        printf("%.4f ", r[i]);
////    }
////    printf("\n");
//
//    for (int i = nsg - 2; i >= 0; i--) {
//        max_r[i] = max(max_r[i], max_r[i + 1]);
//    }
//
//    for (int i = 0; i < nsg; i++) {
//        printf("min %.4f max %.4f\n", min_r[i], max_r[i]);
//    }
//
//    vector<bool> valid(nsg);
//    for (int i = 0; i < nsg - 1; i++) {
//        valid[i] = min_r[i] > max_r[i + 1];
//    }
//    check_first = (nsg == 1) || valid[0];
//
//    //merge stable groups
//    if (nsg > 1) {
//        vector<int> n_nag;
//        for (int i = 0; i < nsg; i++) {
//            if (i == 0 || valid[i - 1]) {
//                n_nag.push_back(nag[i]);
//            } else {
//                n_nag[n_nag.size() - 1] += nag[i];
//            }
//        }
//        nag = n_nag;
//
//        cur = 0;
//        nsg = nag.size();
//        for (int i = 0; i < nsg; i++) {
//            for (int j = cur; j < cur + nag[i]; j++) {
//                sg[slt_nodes[j]] = i;
//            }
//            cur += nag[i];
//        }
//    }
//    printf("updated nsg %d\n", nsg);
//}

void Graph::pruning() {
    vector<double> min_r(nsg);
    vector<double> max_r(nsg);

    int cur = 0;
    for (int i = 0; i < nsg; i++) {
        min_r[i] = r[slt_nodes[cur]];
        max_r[i] = r[slt_nodes[cur]];
        for (int j = cur + 1; j < cur + nag[i]; j++) {
            min_r[i] = min(min_r[i], r[slt_nodes[j]]);
            max_r[i] = max(max_r[i], r[slt_nodes[j]]);
        }
        for (int j = cur; j < cur + nag[i]; j++) {
            selected[slt_nodes[j]] = true;
            deg[slt_nodes[j]] = 0;
            rho_u[slt_nodes[j]] = min(rho_u[slt_nodes[j]], max_r[i]);
            rho_l[slt_nodes[j]] = max(rho_l[slt_nodes[j]], min_r[i]);
            if (rho_u[slt_nodes[j]] < rho_l[slt_nodes[j]]) {
                selected[slt_nodes[j]] = false;
            }
            if (active[slt_nodes[j]]) {
                rho_gu[slt_nodes[j]] = min(rho_gu[slt_nodes[j]], max_r[i]);
            }
        }
        cur += nag[i];
    }

    for (auto e : slt_edges) {
        int u = edges[e].first;
        int v = edges[e].second;
        if (sg[u] > sg[v]) {
            selected[u] = false;
        }
        if (sg[u] < sg[v]) {
            selected[v] = false;
        }
    }

//    for (auto u : slt_nodes) {
//        if (selected[u]) {
//            for (auto v : adj[u]) {
//                if (rho_l[v] > rho_u[u]) {
//                    selected[u] = false;
//                }
//            }
//        }
//    }

    for (auto e : slt_edges) {
        int u = edges[e].first;
        int v = edges[e].second;
        if (selected[u] && selected[v]) {
            ++deg[u];
            ++deg[v];
        }
    }

    queue<int> q;
    for (auto u : slt_nodes) {
        if (selected[u] && deg[u] < rho_l[u]) {
            selected[u] = false;
            q.push(u);
        }
    }

    while (!q.empty()) {
        int u = q.front(); q.pop();
        for (auto v : adj[u]) {
            if (selected[v]) {
                if (--deg[v] < rho_l[v]) {
                    selected[v] = false;
                    q.push(v);
                }
            }
        }
    }
    vector<int> tmp_nodes;
    vector<bool> inactive_sg(nsg, false);
    for (auto u : slt_nodes) {
        if (selected[u])
            tmp_nodes.push_back(u);
        else
            inactive_sg[sg[u]] = true;
    }
    slt_nodes = tmp_nodes;

    for (auto u : slt_nodes) {
        if (inactive_sg[sg[u]]) {
            active[u] = false;
        }
    }

    vector<int> tmp_edges;
    for (auto e : slt_edges) {
        if (selected[edges[e].first] && selected[edges[e].second]) {
            tmp_edges.push_back(e);
        }
    }
    slt_edges = tmp_edges;
    sort(slt_edges.begin(), slt_edges.end(), [this](int a, int b)->bool {
        return sg[edges[a].first] < sg[edges[b].first];
    });

    while (slt_nodes.size() > 0 && sg[slt_nodes.back()] > 0) {
        int n_nodes = slt_nodes.size() - 1;
        while (n_nodes > 0 && sg[slt_nodes[n_nodes]] == sg[slt_nodes.back()]) {
            --n_nodes;
        }
        n_nodes++;
        vector<int> t_nodes(slt_nodes.begin() + n_nodes, slt_nodes.end());
        stk_nodes.push(t_nodes);
        slt_nodes.resize(n_nodes);

        if (sg[t_nodes[0]] > sg[edges[slt_edges.back()].first]) {
            stk_nodes.pop();
            continue;
        }

        int n_edges = slt_edges.size() - 1;
        while (n_edges > 0 && sg[edges[slt_edges[n_edges]].first] == sg[edges[slt_edges.back()].first]) {
            --n_edges;
        }
        n_edges++;
        vector<int> t_edges(slt_edges.begin() + n_edges, slt_edges.end());
        stk_edges.push(t_edges);
        slt_edges.resize(n_edges);
        printf("sg %d %d %lu %lu\n", sg[t_nodes[0]], sg[edges[t_edges[0]].first], t_nodes.size(), t_edges.size());
        stk_CT.push(CT);
    }
}

void Graph::findLDS() {
    compute_core();
    prune_by_core();
    CT = 0;
    while (!slt_nodes.empty() || !stk_nodes.empty()) {
        clock_t start = clock();
        if (slt_nodes.empty()) {
            slt_nodes = stk_nodes.top(); stk_nodes.pop();
            slt_edges = stk_edges.top(); stk_edges.pop();
            CT = stk_CT.top(); stk_CT.pop();
        }
        frank_wolfe();
        clock_t t_fw = clock();
        printf("fw time: %.4f\n", double(t_fw - start) / CLOCKS_PER_SEC);
        CT += NT;
        pava();
        clock_t t_pv = clock();
        printf("pava time: %.4f\n", double(t_pv - t_fw) / CLOCKS_PER_SEC);
        check_sg();
        clock_t t_check_sg = clock();
        printf("check sg time: %.4f\n", double(t_check_sg - t_pv) / CLOCKS_PER_SEC);
        pruning();
        clock_t t_prune = clock();
        printf("pruning time: %.4f\n", double(t_prune - t_check_sg) / CLOCKS_PER_SEC);
        if (!slt_nodes.empty() && check_first) {
            double g = (double) slt_edges.size() / slt_nodes.size();
            vector<pair<int, int>> tmp_edges;
            for (auto e : slt_edges) {
                tmp_edges.push_back(edges[e]);
            }
            //TODO check by edge number? and connected component
            FlowNetwork fn = FlowNetwork(tmp_edges, g);
            double max_flow = fn.get_maxflow(0, fn.n-1);
            if (abs(max_flow - slt_edges.size()) <= 1e-3) {
                printf("ldses candidate: #nodes %lu #edges %lu\n", slt_nodes.size(), slt_edges.size());
                if (verify_LDS(g)) {
                    connected_components();
                    int cur_u = 0;
                    for (auto pr : cmpt) {
                        vector<int> tmp_nodes(slt_nodes.begin() + cur_u, slt_nodes.begin() + pr.first);
                        ldses.push_back(tmp_nodes);
                        lds_rho.push_back(g);
                        cur_u = pr.first;
                    }
//                    ldses.push_back(slt_nodes);
                    if (ldses.size() >= topk)
                        break;
                }
                slt_nodes.clear();
                slt_edges.clear();
            }
        }
        clock_t t_verify_LDS = clock();
        printf("verifyLDS time: %.4f\n", double(t_verify_LDS - t_prune) / CLOCKS_PER_SEC);
    }
}

void Graph::compute_core() {
    max_d = 0;
    for (int i = 0; i < n; i++) {
        deg[i] = adj[i].size();
        max_d = max(deg[i], max_d);
    }
    printf("max_d %d\n", max_d);
    vector<int> bin(max_d + 1, 0);
    pos.resize(n);
    vector<int> vert(n);
    for (int i = 0; i < n; i++) {
        ++bin[deg[i]];
    }
    int start = 0;
    for (int d = 0; d <= max_d; d++) {
        int num = bin[d];
        bin[d] = start;
        start += num;
    }
    for (int i = 0; i < n; i++) {
        pos[i] = bin[deg[i]]++;
        vert[pos[i]] = i;
    }
    for (int d = max_d; d > 0; d--) {
        bin[d] = bin[d - 1];
    }
    bin[0] = 0;
    for (int i = 0; i < n; i++) {
        int u = vert[i];
        for (auto v : adj[u]) {
            if (deg[v] > deg[u]) {
                int dv = deg[v], pv = pos[v];
                int pw = bin[dv], w = vert[pw];
                if (v != w) {
                    pos[v] = pw; vert[pv] = w;
                    pos[w] = pv; vert[pw] = v;
                }
                ++bin[dv]; --deg[v];
            }
        }
    }
}

void Graph::prune_by_core() {
    for (int u = 0; u < n; u++) {
        rho_l[u] = deg[u] / 2.0;
        rho_u[u] = deg[u];
        rho_gu[u] = deg[u];
    }
    vector<double> tmp_rho_l(n);
    queue<int> q;
    for (int u = 0; u < n; u++) {
        tmp_rho_l[u] = rho_l[u];
        for (auto v : adj[u]) {
            tmp_rho_l[u] = max(tmp_rho_l[u], rho_l[v]);
        }
        if (tmp_rho_l[u] > rho_u[u] + DELTA){
            selected[u] = false;
            q.push(u);
        }
    }

    vector<int> vtx;
    while (!q.empty()) {
        int u = q.front(); q.pop();
        vtx.push_back(u);
        for (auto v : adj[u]) {
            if (selected[v] && pos[u] > pos[v]) {
                --deg[v];
                if (tmp_rho_l[v] > deg[v] + DELTA) {
                    selected[v] = false;
                    q.push(v);
                }
            }
        }
    }

    for (int i = 0; i < n; i++) {
        if (selected[i]) slt_nodes.push_back(i);
    }
    for (int i = 0; i < m; i++) {
        if (selected[edges[i].first] && selected[edges[i].second]) {
            slt_edges.push_back(i);
        }
    }
    printf("#slt node %lu #slt edge %lu\n", slt_nodes.size(), slt_edges.size());

    sort(vtx.begin(), vtx.end(), [this](int a, int b)->bool {
            return rho_gu[a] > rho_gu[b];
    });
    for (auto u : vtx) {
        if (active[u]) {
            active[u] = false;
            double threshold = rho_gu[u] + DELTA;
            q.push(u);
            while (!q.empty()) {
                int v = q.front(); q.pop();
                for (auto w : adj[v]) {
                    if (active[w] && rho_l[w] <= threshold) {
                        active[w] = false;
                        q.push(w);
                    }
                }
            }
        }
    }

    int cnt = 0;
    for (int i = 0; i < n; i++) {
        if (active[i]) ++cnt;
    }
    printf("#active nodes %d\n", cnt);
}

bool Graph::verify_LDS(double g) {
    for (auto u : slt_nodes) {
        lds_num[u] = ldses.size();
    }
    bool flag = true;

    ++num_verify;
    queue<int> q;
    vector<pair<int, int>> tmp_edges;
    for (auto u : slt_nodes) {
        if (veri_vtx[u] != num_verify) {
            q.push(u);
            veri_vtx[u] = num_verify;
        }
        while (!q.empty()) {
            int v = q.front(); q.pop();
            for (auto w : adj[v]) {
//                if (slt_nodes.size() == 13 && slt_edges.size() == 68) printf("%d %d %.4f %.4f\n", v, w, rho_gu[w], rho_l[w]);
                if (rho_gu[w] >= g) {
                    if (veri_vtx[w] != num_verify) {
                        if (lds_num[w] != -1 && lds_num[w] < lds_num[u]) {
                            flag = false;
                        }
//                        veri_vtx[w] = num_verify;
//                        q.push(w);
                        if (rho_l[w] <= g) {
                            veri_vtx[w] = num_verify;
                            q.push(w);
                        } else {
                            flag = false;
                            tmp_edges.emplace_back(v, v);
                        }
                    }
                    if (v < w && rho_l[w] <= g)
                        tmp_edges.emplace_back(v, w);
                }
            }
        }
    }
    printf("size of tmp edges %lu\n", tmp_edges.size());
    if (flag) return flag;

    flag = true;
    FlowNetwork fn = FlowNetwork(tmp_edges, g - 1.0 / n / n, true);
    vector<int> tmp_nodes;
    fn.get_mincut(0, fn.n - 1, tmp_nodes);
    ++num_verify;
    printf("size of tmp nodes %lu\n", tmp_nodes.size());
    for (auto u : tmp_nodes) {
//        printf("%d ", u);
        veri_vtx[u] = num_verify;
    }
    printf("\n");
    for (auto u : slt_nodes) {
        for (auto v : adj[u]) {
            if (lds_num[v] != lds_num[u] && veri_vtx[v] == num_verify) {
                flag = false;
                break;
            }
        }
        if (!flag) break;
    }

    if (flag) return flag;
    printf("validate failed\n");
    //TODO add further verification

    for (auto u : slt_nodes) {
        lds_num[u] = -1;
    }
    return false;
}

void Graph::output(char *ds_address) {
    printf("num of LDS %lu\n", ldses.size());
    FILE* dsFile = fopen(ds_address, "w");
    for (int i = 0; i < ldses.size(); i++) {

        fprintf(dsFile, "%lu %.4f\n", ldses[i].size(), lds_rho[i]);
        for (auto & u : ldses[i]) {
            fprintf(dsFile, "%d ", u);
        }
        fprintf(dsFile, "\n");
    }

    fclose(dsFile);
}

int Graph::find_fa(int x) {
    if (x != fa[x])
        fa[x] = find_fa(fa[x]);
    return fa[x];
//    return (x == fa[x]) ? x : (fa[x] = find_fa(fa[x]));
}

void Graph::connected_components() {
    for (auto u : slt_nodes) {
        fa[u] = u;
    }

    for (auto e : slt_edges) {
//        if (find_fa(edges[e].first) != find_fa(edges[e].second))
        int x = find_fa(edges[e].first);
        int y = find_fa(edges[e].second);
//        printf("%d %d %d %d\n", edges[e].first, edges[e].second, x, y);
        fa[x] = y;
    }
    for (auto u : slt_nodes) {
        find_fa(u);
    }

    sort(slt_nodes.begin(), slt_nodes.end(), [this](int a, int b)->bool {
        return fa[a] > fa[b];
    });
    sort(slt_edges.begin(), slt_edges.end(), [this](int a, int b)->bool {
        return fa[edges[a].first] > fa[edges[b].first];
    });

//    printf("fa ");
//    for (int i = 0; i < n; i++) {
//        printf("%d ", find_fa(i));
//    }
//    printf("\n");

    cmpt.clear();
    int x = 1, y = 1;
    while (x < slt_nodes.size() && y < slt_edges.size()) {
        while (x < slt_nodes.size() && fa[slt_nodes[x]] == fa[slt_nodes[x - 1]]) x++;
        while (y < slt_edges.size() && fa[edges[slt_edges[y]].first] == fa[edges[slt_edges[y - 1]].first]) y++;
        cmpt.emplace_back(x, y);
        ++x;
        ++y;
    }
//    for (auto pr : cmpt) {
//        printf("%d %d\n", pr.first, pr.second);
//    }
}








