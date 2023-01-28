#include <unistd.h>
#include <iostream>
#include <fstream>
#include <chrono>
#include "matrify.h"

#include <vector>

using namespace std;

struct chunck {
    int x, y;
//    vector<vector<int>> d_og;

    vector<vector<int>>& d;

    chunk(const chunck& other): d(other.d){}

    chunk(int x = 0, int  y = 0, int m = 1) : x(x), y(y), d(m, std::vector<int>(m)) {}

//    chunck(int x = 0, int y = 0, int m = 1) {
//        i = x;
//        j = y;
//        d_og.resize(m);
//        for (int i = 0; i < m; i++) {
//            d_og[i].resize(m);
//        }
//        d = d_og;
//    }
};

bool CompareChunk(const Chunk &a, const Chunk &b) {
    if (a.x == b.x) {
        return a.y < b.y;
    }
    return a.x < b.x;
}

bool CompareChunk(const Chunk &a, const Chunk &b) {
    if (a.x == b.x) {l
        return a.y < b.y;
    }
    return a.x < b.x;
}

int main(int argc, char *argv[])
{

    int n,m,k;

    ifstream input("input.bin", ios::binary);

    // read n, m, and k from the file
    input.read((char*)&n, 4);
    input.read((char*)&m, 4);
    input.read((char*)&k, 4);

    vector<chunk> chunks(k);

    for (int i = 0; i < k; i++) {
        int x, y;
        input.read((char*)&x, 4);
        input.read((char*)&y, 4);

        chunks[i] = chunk(x,y,m);

        for(int j = 0; j < m; j ++){
            for(int k = 0; k < m;k++){
                input.read((char*)&(chunks[i].d[j][k] ), 1);
            }
        }

    }

    return 0;
}
