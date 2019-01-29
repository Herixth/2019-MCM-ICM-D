#include "Graph_def.h"


// ----------- Vertice ------------
Vertice::Vertice() { __vis = false; }

Vertice::~Vertice() { }

Vertice::Vertice(std::string lab, int type, int maxcup, int num, int I, int II, int III, \
    std::map<Label, double>& pri_c) {
    __lab = lab, __type = type, __maxcup = maxcup, __num = num;
    __priors[__I] = I, __priors[__II] = II, __priors[__III] = III;
    __delta = 0;
    __vis = false;
    if (type == EXTPORT) {
        __maxcup = INF;
        __num = 0, __priors[__I] = 0;
        __priors[__II] = 0, __priors[__III] = 0;
    }
}

std::string Vertice::getLab() const {
    return __lab;
}

int Vertice::getPrNum(Label label) {
    return __priors[label];
}

int Vertice::getNum() const {
    return __num;
}

int Vertice::getDelta() const {
    return __delta;
}

int Vertice::getMaxcup() const {
    return __maxcup;
}

// delta > 0 输出, 更新__delta
// delta < 0 输入, 不更新__delta
void Vertice::update(Label label, int delta) {
    __priors[label] -= delta;
    __num -= delta;

    __delta += (delta > 0) * delta;
}

bool Vertice::operator < (const Vertice& rhs) const {
    return __num < rhs.getNum();
}

bool Vertice::operator == (const Vertice& rhs) const {
    return __num == rhs.getNum();
}

bool Vertice::has_visted() const {
    return __vis;
}

void Vertice::reset_vis(bool val) {
    __vis = val;
}

void Vertice::reset_delta(int delta) {
    __delta = delta;
}

bool Vertice::is_ext() const {
    return __type;
}

void Vertice::clear_num() {
    __num = 0;
}

// ----------- Edge ------------
Edge::Edge() { __vis = false; }

Edge::~Edge() { }

Edge::Edge(int fir, int sec, int pass, int flow, int maxcup) {
    fir_ver = fir, sec_ver = sec;
    passage_type = pass; __flow = flow;
    __maxcup = maxcup; __gama = 0;
    __vis = false;
}

int Edge::getFirver() const {
    return fir_ver;
}

int Edge::getSecver() const {
    return sec_ver;
}

int Edge::getFlow() const {
    return __flow;
}

int Edge::getMaxcup() const {
    return __maxcup;
}

double Edge::getGama() const {
    return __gama;
}

void Edge::reset_gama(double gama) {
    __gama = gama;
}

bool Edge::has_visted() const {
    return __vis;
}

void Edge::reset_vis(bool vis) {
    __vis = vis;
}

int Edge::getType() const {
    return passage_type;
}

ESC_graph::ESC_graph() { }

ESC_graph::~ESC_graph() { 
    inFile.close();
    outFile_TOT.close();
    outFile_I.close();
    outFile_II.close();
    outFile_III.close();
    outFile_ver.close();
    outFile_clock.close();
    outFile_ord_to_ver.close();
    outFile_result.close();
}

ESC_graph::ESC_graph(const char* infilename) {
    inFile.open(infilename);
    outFile_TOT.open("total.txt");
    outFile_I.open("att_I.txt");
    outFile_II.open("att_II.txt");
    outFile_III.open("att_III.txt");
    outFile_ver.open("ver.txt");
    outFile_clock.open("clock.txt");
    outFile_ord_to_ver.open("Conversion.csv");
    outFile_result.open("result.csv");
    if (!inFile.good())
        exit(-1);
}

/**
 * @brief 
 *      construct Graph G(V, E) with abstract Louvre condition
 *      input-file format:
 *          No.1: N M unit_t
 *          For vertices from 0 to N - 1:
 *              No.*: label  type  maxcup num priors(in pri_c's order)
 *          For edges from 0 to M - 1:
 *              No.*: first_vertice  second_vertice  type  flow  maxcup  
 */
void ESC_graph::init() {
    // pri_c     remain
    pri_c[__I]    = 0.6,   remain[__I]   = 0;
    pri_c[__II]   = 0.3,   remain[__II]  = 0;
    pri_c[__III]  = 0.1,   remain[__III] = 0;
    time_cost[__I]   = 0;
    time_cost[__II]  = 0;
    time_cost[__III] = 0;
    time_cost[__TOT] = 0;  remain[__TOT] = 0;
    inFile >> N >> M >> unit_t;
    // vertice
    int type = 0, maxcup = 0, num = 0;
    std::string lab;
    int I = 0, II = 0, III = 0;
    for (int inc = 0; inc < N; inc ++) {
        int __ = 0;
        inFile >> lab >> type >> __;
#ifdef __RELEASE__
        num = __;
        maxcup = int(std::ceil(_maxcup * num));
        if (!maxcup)
            maxcup = stair_maxcup;
        III = int(_III * num);
        II = int(_II * num);
        I = num - III - II;
#else
        maxcup = __;
        inFile >> num >> I >> II >> III;
#endif // __RELEASE__
        if (type != EXTPORT) {
            remain[__TOT] += num, remain[__I] += I;
            remain[__III] += III, remain[__II] += II;
        }
        if (type == EXTPORT)
            ext_ord.push_back(vertices.size());
        vertices.push_back(Vertice(lab, type, maxcup, num, I, II, III, pri_c));
    }
    __sort(ext_ord);
    // edge
    int pass = 0, flow = 0;
    std::string fir, sec;
    for (int inc = 0; inc < M; inc ++) {
        inFile >> fir >> sec >> pass;
#ifdef __RELEASE__
        if (pass == NARROW) flow = flow_NARROW;
        else if (pass == WIDE) flow = flow_WIDE;
        else if (pass == STAIR) flow = flow_STAIR;
        flow = (rand() % flow) + 1;
        maxcup = flow;
#else
        inFile >> flow >> maxcup;
#endif // __RELEASE__
        if (type == DISABLE) {
            M--;
            continue;
        }
        int fir_ord = __find(fir);
        int sec_ord = __find(sec);
        int edge_ord = edges.size();

        index[fir_ord].push_back(edge_ord);
        index[sec_ord].push_back(edge_ord);

        edges.push_back(Edge(fir_ord, sec_ord, pass, flow, maxcup));
    }
}

bool ESC_graph::comp(const int fir, const int sec) {
    return vertices[fir].getNum() >= vertices[sec].getNum();
}


inline int ESC_graph::__find(std::string& label) {
    for (int inc = 0; inc < N; inc++)
        if (vertices[inc].getLab() == label)
            return inc;
    return NULL;
}

inline int ESC_graph::get_adj(int curr, int idx) {
    int adj = edges[index[curr][idx]].getFirver();
    if (adj == curr)
        adj = edges[index[curr][idx]].getSecver();
    return adj;
}

inline void ESC_graph::__clf_ver(int curr) {
    AdjUE.clear(), AdjNE.clear();
    VECTOR(iter, int, index[curr]) {
        int adj = get_adj(curr, iter - index[curr].begin());
        if (vertices[adj].has_visted()) {
            AdjUE.push_back(iter - index[curr].begin());
        }
        else {
            AdjNE.push_back(iter - index[curr].begin());
        }
    }
}

void ESC_graph::__sort(std::vector<int> &obj) {
    int size = obj.size();
    for (int inc = 1; inc < size; inc++)
        for (int snc = inc; snc > 0 && comp(obj[snc], obj[snc - 1]); snc--)
            std::swap(obj[snc], obj[snc - 1]);
}

void ESC_graph::__calc_gama(int curr) {
    int radix = 0;
    VECTOR(iter, int, AdjNE) {
        radix += std::min(unit_t * edges[index[curr][*iter]].getFlow(),
            vertices[get_adj(curr, *iter)].getNum());
    }
    // set gama
    VECTOR(iter, int, AdjNE) {
        double _gama = 0.0;
        if (radix)
            _gama = 1.0 * std::min(unit_t * edges[index[curr][*iter]].getFlow(), \
                vertices[get_adj(curr, *iter)].getNum()) / radix;
        else
            _gama = 1.0 / AdjNE.size();
        

        edges[index[curr][*iter]].reset_gama(_gama);
    }
}

void ESC_graph::__update_num(int curr) {
    if (vertices[curr].is_ext())
        return;
    // real tot num
    // radix
    double radix = 0;
    bool tied_to_ext = false;
    int ext_num = 0;
    VECTOR(iter, int, AdjUE) {
        int adj = get_adj(curr, *iter);
        radix += vertices[adj].getDelta() + vertices[adj].getNum();
        if (vertices[adj].is_ext()) {
            tied_to_ext = true;
            ext_num++;
        }
    }
    std::vector<int> limit_A;
    std::vector<int> limit_B;
    std::vector<int> limit_C;

    int pass_cnt[3] = { 0, 0, 0 };
    VECTOR(iter, int, AdjUE) {
        pass_cnt[edges[index[curr][*iter]].getType()] ++;
    }
    int empty_radix = pass_cnt[2];
    int level = STAIR;
    if (!empty_radix) {
        empty_radix = pass_cnt[1];
        level = WIDE;
    }
    if (!empty_radix) {
        empty_radix = pass_cnt[0];
        level = NARROW;
    }

    VECTOR(iter, int, AdjUE) {
        int adj = get_adj(curr, *iter);
        int edge_ord = index[curr][*iter];
        if (tied_to_ext) {
            if (vertices[adj].is_ext()) {
                limit_A.push_back(vertices[curr].getNum() * 1.0 / ext_num);
            }
            else {
                limit_A.push_back(0);
            }
        }
        else {
            if (radix) {
                double rate = 1.0 * (vertices[adj].getDelta() + vertices[adj].getNum()) / radix;
                limit_A.push_back(rate * vertices[curr].getNum());
            }
            else {
                if (edges[edge_ord].getType() == level)
                    limit_A.push_back(vertices[curr].getNum() * (1.0 / empty_radix));
                else
                    limit_A.push_back(0);
            }
        }
        limit_B.push_back(std::max(0.0, edges[edge_ord].getGama() * (vertices[adj].getMaxcup() - vertices[adj].getNum())));
        limit_C.push_back(edges[edge_ord].getFlow() * unit_t);
    }
    
    int sz = AdjUE.size();
    int sum = 0;
    for (int inc = 0; inc < sz - 1; inc++) {
        sum += limit_A[inc];
    }
    limit_A[sz - 1] = vertices[curr].getNum() - sum;
    std::vector<int> limit_tot;
    for (int inc = 0; inc < sz; inc++) {
        limit_tot.push_back(std::min(limit_C[inc], std::min(limit_A[inc], limit_B[inc])));
    }
    
    // tot num
    radix = vertices[curr].getNum();
    if (!radix)
        return;
    
    for (int inc = 0; inc < sz; inc++) {
        int re = std::min(limit_tot[inc], vertices[curr].getNum());
        
        int delta = 0;
        // __I  high priority
        int gos = int(std::ceil(limit_tot[inc] * 1.0 * vertices[curr].getPrNum(__I) / radix));
        if (gos < 0)
            gos = 0;
        else if (gos > re)
            gos = re;
        else
            gos = std::min(gos, re);
        if (gos > vertices[curr].getPrNum(__I))
            gos = vertices[curr].getPrNum(__I);
        if (vertices[get_adj(curr, AdjUE[inc])].is_ext()) {
            remain[__I] -= gos;
            remain[__TOT] -= gos;
        }
        vertices[curr].update(__I, gos);
        vertices[get_adj(curr, AdjUE[inc])].update(__I, -gos);
        re -= gos;

        // __II middle priority
        gos = int(std::ceil(limit_tot[inc] * 1.0 * vertices[curr].getPrNum(__II) / radix));
        if (gos < 0)
            gos = 0;
        else if (gos > re)
            gos = re;
        else
            gos = std::min(gos, re);
        if (gos > vertices[curr].getPrNum(__II))
            gos = vertices[curr].getPrNum(__II);
        re -= gos;
        if (vertices[get_adj(curr, AdjUE[inc])].is_ext()) {
            remain[__II] -= gos;
            remain[__TOT] -= gos;
        }
        vertices[curr].update(__II, gos);
        vertices[get_adj(curr, AdjUE[inc])].update(__II, -gos);

        // __III low priority
        gos = re;
        if (gos < 0)
            gos = 0;
        else if (gos > re)
            gos = re;
        else
            gos = std::min(gos, re);
        if (gos > vertices[curr].getPrNum(__III))
            gos = vertices[curr].getPrNum(__III);
        if (vertices[get_adj(curr, AdjUE[inc])].is_ext()) {
            remain[__III] -= gos;
            remain[__TOT] -= gos;
        }
        vertices[curr].update(__III, gos);
        vertices[get_adj(curr, AdjUE[inc])].update(__III, -gos);
    }
    if (!vertices[curr].getNum()) {
        order.push_back(vertices[curr].getLab());
    }
}

void ESC_graph::__update_ver() {
    // load ext_ord
    std::queue<int> Qcont;
    VECTOR(iter, int, ext_ord) {
        Qcont.push(*iter);
    }
    
    // opti-BFS
    while (!Qcont.empty()) {
        int Q_curr = Qcont.front();
        Qcont.pop();
        // calc_gama
        __clf_ver(Q_curr);
        __calc_gama(Q_curr);

        // update_num
        __update_num(Q_curr);

        vertices[Q_curr].reset_vis(true);

        // bfs
        AdjUE.clear();
        VECTOR(iter, int, AdjNE) {
            AdjUE.push_back(get_adj(Q_curr, *iter));
        }
        __sort(AdjUE);
        VECTOR(iter, int, AdjUE) {
            Qcont.push(*iter);
        }
    }

}

void ESC_graph::__update_edge() {

}

void ESC_graph::simulate() {
    int cnt = 0;
    // init output file
    for (int inc = 0; inc < N; inc++) {
        if (vertices[inc].is_ext())
            continue;
        outFile_ver << cnt << std::endl;
        outFile_ord_to_ver << cnt << "," << vertices[inc].getLab() << std::endl;
        cnt++;
    }
    cnt = 0;
    while (true) {
        // judge
        if (remain[__I] <= 0 && !time_cost[__I])
            time_cost[__I]   = cnt * unit_t;
        if (remain[__II] <= 0 && !time_cost[__II])
            time_cost[__II]  = cnt * unit_t;
        if (remain[__III] <= 0 && !time_cost[__III])
            time_cost[__III] = cnt * unit_t;

        if (!remain[__TOT]) {
            time_cost[__TOT] = cnt * unit_t;
            break;
        }

        // update vertices
        __update_ver();

        // update edges
        __update_edge();

        // record infor
            // clock
        outFile_clock << (++cnt) * unit_t << std::endl;
            // total
        VECTOR(iter, Vertice, vertices) {
            if ((*iter).is_ext())
                continue;
            outFile_TOT << (*iter).getNum() << " ";
        }
        outFile_TOT << std::endl;
            // I
        VECTOR(iter, Vertice, vertices) {
            if ((*iter).is_ext())
                continue;
            outFile_I << (*iter).getPrNum(__I) << " ";
        }
        outFile_I << std::endl;
            // II
        VECTOR(iter, Vertice, vertices) {
            if ((*iter).is_ext())
                continue;
            outFile_II << (*iter).getPrNum(__II) << " ";
        }
        outFile_II << std::endl;
            // III
        VECTOR(iter, Vertice, vertices) {
            if ((*iter).is_ext())
                continue;
            outFile_III << (*iter).getPrNum(__III) << " ";
        }
        outFile_III << std::endl;
        // reset vertices edges
        VECTOR(iter, Vertice, vertices) {
            (*iter).reset_delta();
            (*iter).reset_vis(false);
            //if ((*iter).is_ext())
            //    (*iter).clear_num();
        }
        
        VECTOR(iter, Edge, edges) {
            (*iter).reset_gama(0.0);
            (*iter).reset_vis(false);
        }
    }

    outFile_result << "I_t,II_t,III_t,TOTAL_t" << std::endl;
    outFile_result << time_cost[__I] << "," << time_cost[__II] << \
       "," << time_cost[__III] << "," << time_cost[__TOT] << std::endl;
    std::cout << "total: " << time_cost[__TOT] << std::endl;
}
