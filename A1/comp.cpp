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

struct Chunk {
    int x , y ;

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

    bool operator==(const Chunk& other) const {
        if (x != other.x || y != other.y || d->size() != other.d->size()) {

            return false;
        }
        auto m = d->size();
        for (int i = 0; i <m; i++) {
            for (int j = 0; j < m; j++) {
                if (d_og[i][j] != d_og[i][j]) {
                    return false;
                }
            }
        }
        return true;
    }

    void print(){
        cout << x << " " << y;
        for(auto x: d_og){
            cout << '\n';
            for(auto y: x){
                cout << y << " ";
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

    cout << n << " " << m << " " << k << endl;

//    return 0;

    k = 210;


    cout << "Starting F \n";

    vector<Chunk> chunks(k);

    for (int i = 0; i < k; i++) {
        int x = 0, y = 0;
        input.read((char*)&x, 4);
        input.read((char*)&y, 4);


        chunks[i] = Chunk(x,y,m);
        for(int j = 0; j < m; j ++){
            for(int t = 0; t < m;t++){
                input.read((char*)&(chunks[i].d_og[j][t] ), 2);
            }
        }

        chunks[i].print();


    }
    sort(chunks.begin(), chunks.end());

    int n2,m2,k2;

    ifstream input2(argv[2], ios::binary);


    input2.read((char*)&n2, 4);
    input2.read((char*)&m2, 4);
    input2.read((char*)&k2, 4);

    k2 = k;

    cout << "Starting S \n";

    vector<Chunk> chunks2(k2);

    for (int i = 0; i < k2; i++) {
        int x, y;
        input2.read((char*)&x, 4);
        input2.read((char*)&y, 4);


        chunks2[i] = Chunk(x,y,m2);
        for(int j = 0; j < m2; j ++){
            for(int t = 0; t < m2;t++){
                input2.read((char*)&(chunks2[i].d_og[j][t] ), 2);
            }
        }

        chunks2[i].print();

    }

    cout << endl;
    cout << endl;
    cout << endl;
    cout << endl;
    cout << endl;
    cout << endl;


    bool equal = true;
//    if (n != n2 || m != m2 || k != k2) {
//        equal = false;
//    } else if (chunks.size() != chunks2.size()) {
//        equal = false;
//    } else {
        for (int i = 0; i < chunks2.size(); i++) {
            if (!(chunks[i] == chunks2[i])) {

                cout << i << " f\n";
                chunks[i].print();

                cout << "s\n";
                chunks2[i].print();

                equal = false;
                break;
            }
        }
//    }

    if (equal) {
        cout << "equal" << endl;
    } else {
        cout << "not equal" << endl;
    }





    return 0;
}

