#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>

using namespace std;

int main(int argc, char* argv[]) {

    string filename1 = "output.txt";
    string filename2 = "v_me.txt";


    // open the file for reading
    ifstream file1(filename1);

    map<pair<int, int>, int> truss_numbers1;

    string line1;
    while (getline(file1, line1)) {

        int u, v, t;
        sscanf(line1.c_str(), "Truss no. of Edge (%d , %d) : %d", &u, &v, &t);
        truss_numbers1[{u, v}] = t;
        truss_numbers1[{v, u}] = t;

//        cout << u << v << t << endl;
    }

    ifstream file2(filename2);

    map<pair<int, int>, int> truss_numbers2;

    string line2;
    while (getline(file2, line2)) {

        int u, v, t;
        sscanf(line2.c_str(), "Truss no. of Edge (%d , %d) : %d", &u, &v, &t);
        truss_numbers2[{u, v}] = t;
        cout << u << " " << v << " " << t << endl;
    }

    // print the map
    for (const auto& entry : truss_numbers1) {

        if(truss_numbers2[{entry.first.first, entry.first.second}] != entry.second){
            cout << "[" << entry.first.first << ", " << entry.first.second << "] = " << entry.second << "| " <<  truss_numbers2[{entry.first.first, entry.first.second}]  << endl;
        }

    }


    return 0;
}
