#pragma once
#ifndef GRAPH_DEF
#include <iostream>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>
#include <queue>
#include <cstdlib>
#include <ctime>
#include <map>

#define __RELEASE__
#define __DEBUG__
// 结点类别
#define EXTPORT  0x00000001  // 出口
#define INTERNAL 0x00000000  // 内部
// 人群分类
#define __I      0x00000000
#define __II     0x00000001
#define __III    0x00000002
#define __TOT    0x00000003
// 通道种类
#define NARROW   0x00000000
#define WIDE     0x00000001
#define STAIR    0x00000002
#define DISABLE  0x00000003

#define INF      0x0fffffff

#define VECTOR(Iter, Type, Container) \
for (std::vector<Type>::iterator Iter = Container.begin();\
Iter != Container.end(); Iter ++)

#ifdef __RELEASE__
const int flow_WIDE = 6;
const int flow_STAIR = 4;
const int flow_NARROW = 1;
const double _I = 0.3;
const double _II = 0.2;
const double _III = 0.5;
const double _maxcup   = 1.2;
const int stair_maxcup = 30;
#endif

typedef int Label;
const int maxN = int(1e4);


class Prior;
class Vertice;
class Edge;
class ESC_graph;


class Vertice {
public:
    Vertice();
    Vertice(std::string, int, int, int, int, int, int,
        std::map < Label, double>& pri_c);
    ~Vertice();

    // prior
    int    getPrNum(Label);
    // -----
    
    std::string getLab() const;
    int getNum() const;
    int getDelta() const;
    int getMaxcup() const;
    
    void update(Label, int);

    bool operator < (const Vertice&) const;
    bool operator == (const Vertice&) const;
    
    bool is_ext() const;
    bool has_visted() const;
    void reset_vis(bool);
    void reset_delta(int = 0);
    void clear_num();
private:
    std::map<Label, int> __priors;
    std::string __lab;
    int __num;
    int __delta;
    int __maxcup;
    Label __type;

    bool __vis;
};

class Edge {
public:
    Edge();
    Edge(int, int, int, int, int);
    ~Edge();

    int getFirver() const;
    int getSecver() const;
    int getFlow() const;
    int getMaxcup() const;
    int getType() const;
    double getGama() const;

    bool has_visted() const;
    void reset_vis(bool);
    void reset_gama(double);
private:
    int fir_ver;
    int sec_ver;

    int passage_type;
    int __flow;
    int __maxcup;
    double __gama;

    bool __vis;
};

class ESC_graph {
public:
    ESC_graph();
    ESC_graph(const char*);
    ~ESC_graph();
    
    void init();
    void simulate();

private:
    bool comp(const int, const int);
    void __sort(std::vector<int> &);
    void __update_ver();
    void __update_edge();

    inline int get_adj(int, int);
    inline int __find(std::string&);
    inline void __clf_ver(int);

    void __calc_gama(int);
    void __update_num(int);
    std::map<Label, double> pri_c;
    std::map<Label, int> time_cost;
    std::map<Label, int> remain;

    std::ifstream inFile;
    std::ofstream outFile_TOT;
    std::ofstream outFile_I;
    std::ofstream outFile_II;
    std::ofstream outFile_III;
    std::ofstream outFile_ver;
    std::ofstream outFile_clock;
    std::ofstream outFile_ord_to_ver;
    std::ofstream outFile_result;

    std::vector<Vertice> vertices;
    std::vector<Edge> edges;
    std::vector<int> index[maxN];
    std::vector<int> ext_ord;

    std::vector<int> AdjUE;
    std::vector<int> AdjNE;
    std::vector<std::string> order;

    int N; // vertice num
    int M; // edge num

    int unit_t;
};

#endif // !GRAPH_DEF
