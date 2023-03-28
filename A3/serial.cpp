#include <iostream>
#include <fstream>
#include <vector>
#include <map>

#include <algorithm>

using namespace std;

struct Node {
    int id;
    int degree;
    vector<int> neighbors;
};

struct edge {
    int u, v;

    edge(int x,int y){
        if(x < y){
            u = x;
            v = y;
        }else{
            u = y;
            v = x;
        }
    }

    bool operator<(const edge& other) const {
        if (u == other.u) {
            return v < other.v;
        }
        return u < other.u;
    }

    // Define the comparison function
//    bool operator<(const edge& other) const {
//        if (u == other.u && v == other.v) {
//            return false;
//        }
//        else if (u == other.v && v == other.u) {
//            return false;
//        }
//        else if (u < other.u || (u == other.u && v < other.v)) {
//            return true;
//        }
//        else {
//            return false;
//        }
//    }
};


int main() {
    ifstream infile("test4.gra", ios::binary);
    
    // Read n and m
    int n, m;
    infile.read(reinterpret_cast<char*>(&n), sizeof(n));
    infile.read(reinterpret_cast<char*>(&m), sizeof(m));
    
    // Read node information
    vector<Node> nodes(n);
    for (int i = 0; i < n; i++) {
        infile.read(reinterpret_cast<char*>(&nodes[i].id), sizeof(nodes[i].id));
        infile.read(reinterpret_cast<char*>(&nodes[i].degree), sizeof(nodes[i].degree));
        nodes[i].neighbors.resize(nodes[i].degree);
        infile.read(reinterpret_cast<char*>(nodes[i].neighbors.data()), sizeof(int) * nodes[i].degree);
    }

    infile.close();


    map<edge,int> sups;
    map<edge,int> settled;


    for(auto n: nodes){
        for(auto v: n.neighbors){
                sups[{n.id,v}] = 0;
        }
    }

    for(auto &v: sups){
        for(int i = 0; i < n; i++){
            if(sups.count({v.first.u,i}) && sups.count({v.first.v,i})){
                sups[v.first]++;
            }
        }
    }

    auto sups_copy = sups;

//    for(auto n: nodes){
//        for(auto v1: n.neighbors){
//            for(auto v2: n.neighbors){
//                if(v1 != v2){
//                    if(sups.count({v1,v2})){
//                        sups[{v1,v2}] ++;
////                        sups[{n.id,v1}] ++;
////                        sups[{n.id,v2}] ++;
//                    }
//                }
//            }
//        }
//    }

    int k = 3;

    while(sups.size() > 0){
        vector<edge> F_k;

        for(auto x: sups){
            if(x.second < k - 2){
                F_k.push_back(x.first);
            }
        }
        while (F_k.size() > 0){
            for(auto e: F_k){
                for(int i = 0; i < n; i++){
                    if(sups.count({e.u,i}) && sups.count({e.v,i})){
                        sups[{e.u,i}] --;
                        sups[{e.v,i}] --;
                    }
                }

                settled[e] = k - 1;
                sups.erase(e);
            }
            F_k.clear();
            for(auto x: sups){
                if(x.second < k - 2){
                    F_k.push_back(x.first);
                }
            }
        }
        k++;
    }

    for(auto x: settled){
        cout << "Edge (" << x.first.u<<" , " << x.first.v<<")  : Truss no. = " << x.second << " Support = " << sups_copy[x.first] << endl;

//        cout << "Truss no. of Edge (" << x.first.u<<" , " << x.first.v<<") : " << x.second << endl;
    }


    return 0;
}
