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
    int x, y, m;

    T* d;


    Chunk(const Chunk<T>& other):x(other.x), y(other.y), m(other.m), d(other.d){}

    Chunk(int x = 0, int  y = 0, int m = 1) : x(x), y(y), m(m) {
        d = (T*)std::malloc(sizeof(T)*m*m);
    }

    bool operator<(const Chunk<T>& other) const {
        if (x != other.x) {
            return x < other.x;
        } else {
            return y < other.y;
        }
    }

    bool operator==(const Chunk<T>& other) const {
        if (x != other.x || y != other.y || m != other.m) {
            return false;
        }

//        if(x == 8 && y == 1121){
//            for (int i = 0; i <m; i++) {
//                for (int j = 0; j < m; j++) {
//                    cout << d[m*i + j] << " ";
//                }
//            }
//            cout << endl;
//
//        }

        for (int i = 0; i <m; i++) {
            for (int j = 0; j < m; j++) {
                if (d[m*i + j] != other.d[m*i + j]) {
                    return false;
                }
            }
        }
        return true;
    }

    const void print() const{
        cout << x << " " << y;


        for(int i = 0; i < m; ++i){
            cout << "\n";
            for (int j = 0; j < m; ++j) {
                cout << d[i*m + j] << " ";
            }
        }
        cout << "\n";
    }




};



int main(int argc, char *argv[])
{

    int n,m,k;

    ifstream input(argv[1], fstream::binary);


    input.read((char*)&n, 4);
    input.read((char*)&m, 4);
    input.read((char*)&k, 4);

    if (argc == 2){

        cout << k << endl;

        return  0;
    }


    int n2,m2,k2;

    ifstream input2(argv[2], ios::binary);


    input2.read((char*)&n2, 4);
    input2.read((char*)&m2, 4);
    input2.read((char*)&k2, 4);

    //IMP

//    k2 = k = 210;
//    k2 = 210;
//    k = 210;
//    cout << n << " " << m << " " << k << endl;
//    cout << n << " " << m << " " << k2 << endl;



    cout << "Starting F \n";

    vector<Chunk<unsigned short >> chunks(k);

    for (int i = 0; i < k; i++) {
        int x = 0, y = 0;
        input.read((char*)&x, 4);
        input.read((char*)&y, 4);



        chunks[i] = Chunk<unsigned short>(x,y,m);

        input.read((char *)(chunks[i].d), 2 * m *  m);


    }

    for(auto chunk:chunks){
        if(chunk.x == chunk.y){
            cout << "Hm " << chunk.x << endl;
        }

//        cout << chunk.x << " " << chunk.y << endl;
    }

    sort(chunks.begin(), chunks.end());
////
//    for (int i = 0; i < chunks.size(); ++i) {
//        chunks[i].print();
//    }


    cout << "Starting S \n";

    vector<Chunk<unsigned short >> chunks2(k2);

    for (int i = 0; i < k2; i++) {
        int x, y;
        input2.read((char*)&x, 4);
        input2.read((char*)&y, 4);


        chunks2[i] = Chunk<unsigned short >(x,y,m2);
        input2.read((char *)(chunks2[i].d), 2 * m *  m);

//        for(int j = 0; j < m2; j ++){
//            for(int t = 0; t < m2;t++){
//                input2.read((char*)&(chunks2[i].d_og[j][t] ), 2);
//            }
//        }

//        chunks2[i].print();

    }

    sort(chunks2.begin(), chunks2.end());

    for(auto chunk:chunks2){
        if(chunk.x == chunk.y){
            cout << "Hm " << chunk.x << endl;
        }

//        cout << chunk.x << " " << chunk.y << endl;
    }

    sort(chunks.begin(), chunks.end());

//    for (int i = 0; i < chunks2.size(); ++i) {
//        chunks2[i].print();
//    }

    cout << endl;
    cout << endl;
    cout << endl;
    cout << endl;
    cout << endl;
    cout << endl;

    int t = 0;
    bool equal = true;
    if (n != n2 || m != m2 || k != k2) {
        equal = false;
    } else if (chunks.size() != chunks2.size()) {
        equal = false;
    } else {
        for (int i = 0; i < chunks2.size(); i++) {
            if (!(chunks[i] == chunks2[i])) {

                cout << i << " f\n";
                chunks[i].print();

                cout << "s\n";
                chunks2[i].print();

                equal = false;
                ++t;
                break;
            }
        }
    }

    cout << t;
    if (equal) {
        cout << "equal" << endl;
    } else {
        cout << "not equal" << endl;
    }





    return 0;
}

