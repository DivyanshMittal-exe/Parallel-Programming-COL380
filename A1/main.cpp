#include <iostream>
#include <algorithm>

#include <vector>
#include <map>
#include <set>
#include <omp.h>

#include <chrono>

#include "library.hpp"

//#define  time_test


using namespace std;

template<typename T>
struct Chunk {
    int x, y;

    T *d;


    Chunk(const Chunk<T> &other) : x(other.x), y(other.y), d(other.d) {}

    Chunk(int x = 0, int y = 0, int m = 1) : x(x), y(y) {
        d = (T *) std::malloc(sizeof(T) * m * m);
    }

    bool operator<(const Chunk<T> &other) const {
        if (x != other.x) {
            return x < other.x;
        } else {
            return y < other.y;
        }
    }

    void clr(int m) {
//        cout << "clearing" << sizeof(d) << endl;
        memset(d, 0, sizeof(T)*m*m);

    }

    const void print() const{
        cout << x << " " << y;

        int m = 5;
        for(int i = 0; i < m; ++i){
            cout << "\n";
            for (int j = 0; j < m; ++j) {
                cout << d[i*m + j] << " ";
            }
        }
        cout << "\n";
    }

};

void mult_add(int chunk_size, Chunk<unsigned int> &result, const Chunk<unsigned char> &r1, const Chunk<unsigned char> &r2) {
    assert(result.x == r1.x);
    assert(result.y == r2.y);
    assert(r1.y == r2.x);
    int i, j, k;




    if (r1.x <= r1.y && r2.x <= r2.y) {
        for (i = 0; i < chunk_size; ++i) {
            for (j = 0; j < chunk_size; ++j) {
//#pragma  omp simd
                for (k = 0; k < chunk_size; ++k) {
                    result.d[i * chunk_size + k] = Outer(result.d[i * chunk_size + k],Inner((int)r1.d[i * chunk_size + j] , (int)r2.d[j * chunk_size + k])) ;
                }
            }
        }

    } else if (r1.y <= r1.x && r2.x <= r2.y) {

        for (j = 0; j < chunk_size; ++j) {
            for (i = 0; i < chunk_size; ++i) {
//#pragma  omp simd
                for (k = 0; k < chunk_size; ++k) {
                    result.d[i * chunk_size + k] = Outer(result.d[i * chunk_size + k],Inner((int)r1.d[j * chunk_size + i] , (int)r2.d[j * chunk_size + k])) ;
                }
            }
        }

    } else if (r1.x <= r1.y && r2.y <= r2.x) {
        for (i = 0; i < chunk_size; ++i) {
            for (k = 0; k < chunk_size; ++k) {
//#pragma  omp simd
                for (j = 0; j < chunk_size; ++j) {

                    result.d[i * chunk_size + k] =Outer(result.d[i * chunk_size + k],  Inner((int)r1.d[i * chunk_size + j] ,(int)r2.d[k * chunk_size + j]));
                }
            }
        }
    } else {
        cout << "Hmmmmmm";
    }


}


int main(int argc, char *argv[]) {

#ifdef time_test
    auto start = chrono::high_resolution_clock::now();
#endif


    int n, m, chunk_count;

    ifstream input(argv[1], ios::binary);


    input.read((char *) &n, 4);
    input.read((char *) &m, 4);
    input.read((char *) &chunk_count, 4);

    vector<Chunk<unsigned char>> chunks(chunk_count);

    int i_i_pairs = 0;

    set<int> indices;


    for (int i = 0; i < chunk_count; i++) {
        int x = 0, y = 0;
        input.read((char *) &x, 4);
        input.read((char *) &y, 4);

        indices.insert(x);
        indices.insert(y);

        chunks[i] = Chunk<unsigned char>(x, y, m);



        input.read((char *) (chunks[i].d), m * m);

        if(x == y){
            i_i_pairs ++;
        }


    }

    for (int i = 0; i < chunk_count; ++i) {

        if(chunks[i].x != chunks[i].y){
            Chunk<unsigned char> new_chunk = chunks[i];
            new_chunk.x = chunks[i].y;
            new_chunk.y = chunks[i].x;
            chunks.push_back(new_chunk);
        }
    }

//
//#pragma  omp parallel
//    {
//#pragma  omp single
//        {
//            for (int i = 0; i < chunk_count; ++i) {
//
//#pragma  task
//                {
//                    chunks[chunk_count + i] = chunks[i];
//                    chunks[chunk_count + i].x = chunks[i].y;
//                    chunks[chunk_count + i].y = chunks[i].x;
//
//                }
//
//            }
//#pragma  omp taskwait
//        }
//
//    }
//


    sort(chunks.begin(), chunks.end());

    vector<map<int, Chunk<unsigned int>>> f_output(indices.size());


    vector<int> indices_vec(indices.begin(), indices.end());


#pragma  omp parallel
    {
#pragma  omp single
        {


            for (int ind_index = 0; ind_index < indices_vec.size(); ind_index++) {
#pragma  omp task shared(f_output, chunks, indices_vec, m) firstprivate(ind_index)
                {
                    int i_index = indices_vec[ind_index];

                    auto row_element = std::lower_bound(chunks.begin(), chunks.end(),
                                                        Chunk<unsigned char>{i_index, 0, {}});
                    for (; row_element != chunks.end() && row_element->x == i_index; row_element++) {
                        int j_index = row_element->y;
                        auto col_element = std::lower_bound(chunks.begin(), chunks.end(),
                                                            Chunk<unsigned char>{j_index, i_index, {}});

                        for (; col_element != chunks.end() && col_element->x == j_index; col_element++) {
                            int k_index = col_element->y;


                            if (f_output[ind_index].count(k_index) == 0) {
                                f_output[ind_index][k_index] = Chunk<unsigned int>(i_index, k_index, m);
                                f_output[ind_index][k_index].clr(m);
                            }

                            mult_add(m, f_output[ind_index][k_index], *row_element, *col_element);


                        }
                    }

                }


            }
#pragma  omp taskwait
        }


    }


    int final_chunk_count = 0;
#pragma  omp parallel for reduction(+:final_chunk_count)
    for (int i = 0; i < f_output.size(); i++) {
        final_chunk_count += f_output[i].size();
    }

    vector<vector<Chunk<unsigned short>>> f_output_us(f_output.size());
    for (int i = 0; i < f_output.size(); i++) {
        f_output_us[i].resize(f_output[i].size());
    }

    const int max_val = 0xFFFF;
    Chunk<unsigned int> zero_chunk(-1,-1,m);
    zero_chunk.clr(m);

#pragma  omp parallel
    {
#pragma  omp single
        {
            for (int i = 0; i < f_output.size(); ++i) {
#pragma  omp task firstprivate(i) shared(final_chunk_count,f_output,m,zero_chunk,f_output_us) default(none)
                {
                    int j = 0;
                    for (auto it = f_output[i].begin(); it != f_output[i].end();) {
                        auto const &[k_val, chunk_val] = *it;
                        if (memcmp(chunk_val.d, zero_chunk.d, m * m * sizeof(int)) == 0) {
                            it = f_output[i].erase(it);
#pragma  omp atomic
                            final_chunk_count--;
                        } else {
                            ++it;
                            f_output_us[i][j] = Chunk<unsigned short>(chunk_val.x, chunk_val.y, m);
                            for (int k = 0; k < m; ++k) {
                                for (int l = 0; l < m; ++l) {
//                                     f_output_us[i][j].d[m * k + l] =  chunk_val.d[m * k + l];

                                    if(chunk_val.d[m * k + l] < max_val){
                                        f_output_us[i][j].d[m * k + l] =  chunk_val.d[m * k + l];
                                    }else{
                                        f_output_us[i][j].d[m * k + l] = max_val;
                                    }



                                }
                            }
                            ++j;
                        }
                    }
                }
            }
#pragma  omp taskwait
        }
    }



    ofstream output(argv[2], ios::binary);


    output.write((char *) &n, 4);
    output.write((char *) &m, 4);
    output.write((char *) &final_chunk_count, 4);


//    #pragma  omp parallel for
    for (int i = 0; i < f_output.size(); ++i) {
        for (int j = 0; j < f_output[i].size(); ++j) {
            const auto &f_ch_to_p = f_output_us[i][j];
//        #pragma  omp critical
            {
                output.write((char *) &(f_ch_to_p.x), 4);
                output.write((char *) &(f_ch_to_p.y), 4);
                output.write((char *) (f_ch_to_p.d), 2 * m * m);
            }

        }
    }

#ifdef time_test

    output.close();
  auto stop = chrono::high_resolution_clock::now();
  auto duration = chrono::duration_cast<std::chrono::milliseconds>(stop - start);

    cout << duration.count() ;
#endif

    return 0;
}
