#include <iostream>
#include <algorithm>
#include <fstream>

#include <vector>
#include <map>
#include <set>

#include <random>
#include <vector>

using namespace std;

int main(int argc, char *argv[]){

    std::random_device rd;
    std::mt19937 gen(rd());

//    cout << argv[1];


    int shift_by = std::stoi(argv[1]);


    int num_chunk = 1<<shift_by;
    int m = 10;
    int n = 10*num_chunk;

    int max_r = n/m;

    uniform_int_distribution<> index_dis(0, max_r - 1);
    uniform_int_distribution<> chunk_dis(0, 255);

    set<pair<int,int>> indices;

    std::string str = "test_case_";
    str += argv[1];

    ofstream output(str, ios::binary);

    output.write((char *) &n, 4);
    output.write((char *) &m, 4);
    output.write((char *) &num_chunk, 4);

    int chunk_p = 0;

    while (chunk_p < num_chunk){
        int a = index_dis(gen);
        int b = index_dis(gen);

        int i = min(a,b);
        int j = max(a,b);

        if(!indices.count({i,j})){
            indices.insert({i,j});
            output.write((char *) &i, 4);
            output.write((char *) &j, 4);

            for (int k = 0; k < m; ++k) {
                for (int l = 0; l < m; ++l) {
                    char v =  chunk_dis(gen);
                    output.write((char *)&v,1);
                }
            }

            ++chunk_p;
        }


    }





    return 0;
}


