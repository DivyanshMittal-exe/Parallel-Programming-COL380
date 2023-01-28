#include <unistd.h>
#include <iostream>
#include <fstream>
#include <chrono>
#include <algorithm>
#include "matrify.h"

#include <vector>
#include <map>
#include <set>
#include <omp.h>

using namespace std;

struct Chunk {
    int x, y;
    
    vector<vector<int>> d_og;
    vector<vector<int>>* d = &d_og;

    Chunk(const Chunk& other):x(other.x), y(other.y), d(other.d){}

    Chunk(int x = 0, int  y = 0, int m = 1) : x(x), y(y), d_og(m, std::vector<int>(m)) {
        d = &d_og;
    }

};

bool compare_by_row(const Chunk &a, const Chunk &b) {
    if (a.x == b.x) {
        return a.y < b.y;
    }
    return a.x < b.x;
}

bool compare_by_col(const Chunk &a, const Chunk &b) {
    if (a.y == b.y) {
        return a.x < b.x;
    }
    return a.y < b.y;
}

int main(int argc, char *argv[])
{

    int n,m,k;

    ifstream input("input.bin", ios::binary);


    input.read((char*)&n, 4);
    input.read((char*)&m, 4);
    input.read((char*)&k, 4);

    vector<Chunk> chunks(2*k);
    set<int> indices;
//    map<int,int> indices;
//    vector<Chunk> chunks_by_col(k);

    for (int i = 0; i < k; i++) {
        int x, y;
        input.read((char*)&x, 4);
        input.read((char*)&y, 4);

            indices.insert(x);
        indices.insert(y);

        chunks[i] = Chunk(x,y,m);
        for(int j = 0; j < m; j ++){
            for(int k = 0; k < m;k++){
                input.read((char*)&(chunks[i].d_og[j][k] ), 1);
            }
        }
        chunks[k + i] = chunks[i];
        chunks[k+i].x = chunks[i].y;
        chunks[k+i].y = chunks[i].x;
    }
    sort(chunks.begin(), chunks.end(), compare_by_row);


    #pragma omp parallel for
    for (auto it = indices.begin(); it != indices.end(); ++it) {
        int ind = *it;
        auto it = std::lower_bound(chunks.begin(), chunks.end(), chunk{x, 0, {}});


    }
//    sort(chunks_by_row.begin(), chunks_by_row.end(), compare_by_row);
//    sort(chunks_by_col.begin(), chunks_by_col.end(), compare_by_col);



    return 0;
}
