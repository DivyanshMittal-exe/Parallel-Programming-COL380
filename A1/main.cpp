//#include <unistd.h>
#include <iostream>
#include <fstream>
#include <chrono>
#include <algorithm>

#include <vector>
#include <map>
#include <set>
#include <omp.h>

#include "library.hpp"

using namespace std;

template <typename T>
struct Chunk {
    int x, y;

    T* d;


    Chunk(const Chunk<T>& other):x(other.x), y(other.y), d(other.d){}

    Chunk(int x = 0, int  y = 0, int m = 1) : x(x), y(y) {
        d = (T*)std::malloc(sizeof(T)*m*m);
    }

    bool operator<(const Chunk<T>& other) const {
        if (x != other.x) {
            return x < other.x;
        } else {
            return y < other.y;
        }
    }

    void clr(int m){
        for (int i = 0; i < m*m; ++i) {
            d[i] = 0;
        }
    }


};

void mult_add(int chunk_size,Chunk<int> & result, Chunk<unsigned char > & r1,Chunk<unsigned char > & r2 ){
    assert(result.x == r1.x);
    assert(result.y == r2.y);
    assert(r1.y == r2.x);
    int i,j,k;

    if(r1.x <= r1.y && r2.x <= r2.y){
        for (i = 0; i < chunk_size; ++i) {
            for (j = 0; j < chunk_size; ++j) {
                for (k = 0; k < chunk_size; ++k) {
                        result.d[i*chunk_size + k] += r1.d[i*chunk_size + j] * r2.d[j*chunk_size + k];
                    }
                }
            }

    }else if(r1.y <= r1.x && r2.x <= r2.y){

        for (j = 0; j < chunk_size; ++j) {
            for (i = 0; i < chunk_size; ++i) {
                for (k = 0; k < chunk_size; ++k) {
                    result.d[i*chunk_size + k] += r1.d[j*chunk_size + i] * r2.d[j*chunk_size + k];
                }
            }
        }

    }else if(r1.x <= r1.y && r2.y <= r2.x){
        for (i = 0; i < chunk_size; ++i) {
            for (k = 0; k < chunk_size; ++k) {
                for (j = 0; j < chunk_size; ++j) {

                    result.d[i*chunk_size + k] += r1.d[i*chunk_size + j] * r2.d[k*chunk_size + j];
                }
            }
        }
    }else{
        cout << "Hmmmmmm";
    }

//    int (*array_ptr)[chunk_size] = (int(*)[chunk_size])result.d;

//    int (arr)[chunk_size][chunk_size] = *reinterpret_cast<int (*)[chunk_size][chunk_size]>(result.d);




//
//    if (r1.x > r1.y) {
//        if (r2.x > r2.y) {
//            for (i = 0; i < chunk_size; ++i) {
//                for (k = 0; k < chunk_size; ++k) {
//                    for (j = 0; j < chunk_size; ++j) {
//                        result.d[k][i] += r1.d[j][i] * r2.d[k][j];
//                    }
//                }
//            }
//        } else {
//            for (i = 0; i < chunk_size; ++i) {
//                for (k = 0; k < chunk_size; ++k) {
//                    for (j = 0; j < chunk_size; ++j) {
//                        result.d[k][i] += r1.d[j][i] * r2.d[j][k];
//                    }
//                }
//            }
//        }
//    } else {
//        if (r2.x > r2.y) {
//            for (i = 0; i < chunk_size; ++i) {
//                for (k = 0; k < chunk_size; ++k) {
//                    for (j = 0; j < chunk_size; ++j) {
//                        result.d[k][i] += r1.d[i][j] * r2.d[k][j];
//                    }
//                }
//            }
//        } else {
//            for (i = 0; i < chunk_size; ++i) {
//                for (k = 0; k < chunk_size; ++k) {
//                    for (j = 0; j < chunk_size; ++j) {
//                        result.d[k][i] += r1.d[i][j] * r2.d[j][k];
//                    }
//                }
//            }
//        }
//    }





}



int main(int argc, char *argv[])
{

    int n,m,chunk_count;

    ifstream input(argv[1], ios::binary);


    input.read((char*)&n, 4);
    input.read((char*)&m, 4);
    input.read((char*)&chunk_count, 4);

    vector<Chunk<unsigned char>> chunks(2*chunk_count);
    set<int> indices;
//    map<int,int> indices;
//    vector<Chunk> chunks_by_col(k);

    for (int i = 0; i < chunk_count; i++) {
        int x = 0, y = 0;
        input.read((char*)&x, 4);
        input.read((char*)&y, 4);

        indices.insert(x);
        indices.insert(y);

        chunks[i] = Chunk<unsigned char>(x,y,m);

        input.read((char*)(chunks[i].d), m*m);


    }

//    cout << indices.size();

    for (int i = 0; i < chunk_count; ++i) {
        chunks[chunk_count + i] = chunks[i];
        chunks[chunk_count+i].x = chunks[i].y;
        chunks[chunk_count+i].y = chunks[i].x;
    }


    sort(chunks.begin(), chunks.end());

    vector<map<int,Chunk<int>>> f_output(indices.size());

//    int indices_index = 0;
    vector<int> indices_vec(indices.begin(), indices.end());


//#pragma omp parallel for  schedule(dynamic, 1)
    for (int ind_index = 0; ind_index < indices_vec.size(); ind_index ++) {
        int i_index = indices_vec[ind_index];

        auto row_element = std::lower_bound(chunks.begin(), chunks.end(), Chunk<unsigned char>{i_index, 0, {}});
        for(;row_element != chunks.end() && row_element->x == i_index; row_element++){
            int j_index = row_element->y;
            auto col_element = std::lower_bound(chunks.begin(), chunks.end(), Chunk<unsigned char>{j_index, i_index, {}});

            for(; col_element != chunks.end() && col_element->x == j_index; col_element++){
                int k_index = col_element->y;


                if(f_output[ind_index].count(k_index) == 0){
                    f_output[ind_index][k_index] = Chunk<int>(i_index, k_index, m);
                    f_output[ind_index][k_index].clr(m);
                }

                mult_add(m,f_output[ind_index][k_index],*row_element,*col_element);

//                f_output[ind_index][col_element->y] += (*row_element) * (*col_element);


            }
        }

    }

    int final_chunk_count = 0;
//    #pragma omp parallel for reduction(+:final_chunk_count)
        for (int i = 0; i < f_output.size(); i++) {
            final_chunk_count += f_output[i].size();
        }

    ofstream output(argv[2], ios::binary);


    output.write((char*)&n, 4);
    output.write((char*)&m, 4);
    output.write((char*)&final_chunk_count, 4);
    const int max_val = 0xFFFF;



    for(const auto int_chunk_map: f_output){
        for (const auto &[k_val, chunk_val] : int_chunk_map) {
            output.write((char*)&(chunk_val.x), 4);
            output.write((char*)&(chunk_val.y), 4);

            for (int i = 0; i < m*m; ++i) {
                if(chunk_val.d[i] > max_val){
                        output.write((char*)&(max_val), 2);
                    }else{
                        output.write((char*)&(chunk_val.d[i]), 2);
                    }
            }
//            for(int j = 0; j < m; j ++){
//                for(int t = 0; t < m;t++){
//                    if((*(chunk_val.d))[j][t] > max_val){
//                        output.write((char*)&(max_val) + 2, 2);
//                    }else{
//                        output.write((char*)&((*(chunk_val.d))[j][t]) + 2, 2);
//                    }
//
//                }
//            }
        }
    }

//    sort(chunks_by_row.begin(), chunks_by_row.end(), compare_by_row);
//    sort(chunks_by_col.begin(), chunks_by_col.end(), compare_by_col);



    return 0;
}
