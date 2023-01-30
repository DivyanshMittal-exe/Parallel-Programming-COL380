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

#include "library.hpp"

using namespace std;

struct Chunk {
    int x, y;
    
    vector<vector<int>> d_og;
    vector<vector<int>>*d = &d_og;

    Chunk(const Chunk& other):x(other.x), y(other.y), d(other.d){}

    Chunk(int x = 0, int  y = 0, int m = 1) : x(x), y(y), d_og(m, std::vector<int>(m,0)), d(&d_og) {}

    bool operator<(const Chunk& other) const {
        if (x != other.x) {
            return x < other.x;
        } else {
            return y < other.y;
        }
    }

    Chunk operator*(const Chunk& other) const {
        assert(y == other.x);
        assert(d->size() == other.d->size());

        int m = d->size();
        Chunk result(x, other.y, m);
        for (int i = 0; i < m; i++) {
            for (int k = 0; k < m; k++) {
                for (int j = 0; j < m; j++) {
                    (*(result.d))[i][j] = Outer((*(result.d))[i][j], Inner((*d)[i][k], (*(other.d))[k][j]));
                }
            }
        }
        return result;
    }

    Chunk operator+(const Chunk& other) const {
        assert(d->size() == other.d->size());
        assert(x == other.x && y == other.y);

        int m = d->size();
        Chunk result(x, y, m);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) {
                (*(result.d))[i][j] = Outer((*d)[i][j], (*(other.d))[i][j]);
            }
        }
        return result;
    }

    Chunk& operator+=(const Chunk& other) {
        assert(d->size() == other.d->size());
        assert(x == other.x && y == other.y);

        int m = d->size();
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) {
                (*d)[i][j] = Outer((*d)[i][j], (*(other.d))[i][j]);
            }
        }
        return *this;
    }

};



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
            for(int t = 0; t < m;t++){
                input.read((char*)&(chunks[i].d_og[j][t] ), 1);
            }
        }
        chunks[k + i] = chunks[i];
        chunks[k+i].x = chunks[i].y;
        chunks[k+i].y = chunks[i].x;
    }
    sort(chunks.begin(), chunks.end());

    vector<map<int,Chunk>> f_output(indices.size());

    int indices_index = 0;

    #pragma omp parallel for  schedule(dynamic, 1) shared(indices_index)
    for (auto it = indices.begin(); it != indices.end(); ++it) {
        int ind = *it;
        int ind_index = -1;
        #pragma omp atomic capture
            ind_index = indices_index++;
        auto it_to_row = std::lower_bound(chunks.begin(), chunks.end(), Chunk{ind, 0, {}});
        for(;it_to_row != chunks.end() && it_to_row->x == ind; it_to_row++){
            auto it_to_col = std::lower_bound(chunks.begin(), chunks.end(), Chunk{it_to_row->y, 0, {}});

            for(;it_to_col!= chunks.end() && it_to_col->x == it_to_row->y && it_to_col->y <= ind; it_to_col++){
                if(f_output[ind_index].count(it_to_col->y) == 0){
                    f_output[ind_index][it_to_col->y] = Chunk(ind,it_to_col->y,m);
                }
                f_output[ind_index][it_to_col->y] += (*it_to_row)*(*it_to_col);


            }
        }

    }

    int final_chunk_count = 0;
    #pragma omp parallel for reduction(+:final_chunk_count)
        for (int i = 0; i < f_output.size(); i++) {
            final_chunk_count += f_output[i].size();
        }

    ofstream output("output.bin", ios::binary);


    output.write((char*)&n, 4);
    output.write((char*)&m, 4);
    output.write((char*)&final_chunk_count, 4);

    for(const auto int_chunk_map: f_output){
        for (const auto &[k_val, chunk_val] : int_chunk_map) {
            output.write((char*)&(chunk_val->x), 4);
            output.write((char*)&(chunk_val->y), 4);

            for(int j = 0; j < m; j ++){
                for(int t = 0; t < m;t++){
                    output.write((char*)&((chunk_val->d)[j][t]) + 2, 2);
                }
            }
        }
    }

//    sort(chunks_by_row.begin(), chunks_by_row.end(), compare_by_row);
//    sort(chunks_by_col.begin(), chunks_by_col.end(), compare_by_col);



    return 0;
}
