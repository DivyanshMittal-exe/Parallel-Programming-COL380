#include <iostream>
#include <fstream>
#include <vector>

using namespace std;
vector<vector<int>> read_graph_file(string filename) {
    ifstream in(filename, ios::in);
    int num_nodes, num_edges;
    in >> num_nodes >> num_edges;
    vector<vector<int>> graph(num_nodes);
    for (int i = 0; i < num_nodes; i++) {
        int node_i, deg_i;
        in >> node_i >> deg_i;
        for (int j = 0; j < deg_i; j++) {
            int neighbor_j;
            in >> neighbor_j;
            graph[node_i].push_back(neighbor_j);
        }
    }
    in.close();
    return graph;
}

int main(int argc, char *argv[]) {
    if (argc != 2 && argc != 3) {
        return 1;
    }

//    vector<vector<int>> graph = {{1,2,3,4,5},{0,2,3,4,6},{0,1,3},{0,1,2,4},{0,1,3},{0,6},{1,5}};
//    vector<vector<int>> graph = {{1},{0,2,3},{1,3},{1,2,4,5}};
//    vector<vector<int>> graph = {{1},{0,2},{1,3,4},{2,4},{2,3,5,6},{4,6,7,8},{4,5,7,8},{5,6,8},{5,6,7}};
    vector<vector<int>> graph = {{1,2,4},{0,2,4},{1,0,3,4},{2,4},{0,1,2,3}};

    string filename;
    if(argc == 2){
         filename = argv[1];

    }else{
        graph = read_graph_file(argv[1]);
        filename = argv[2];
    }

    int num_nodes = graph.size();
    int num_edges = 0;

    for (int i = 0; i < num_nodes; i++) {
        num_edges += graph[i].size();
    }

    num_edges/=2;

    ofstream out(filename + ".gra", ios::binary | ios::out);

    string mdat_name = filename + ".dat";
    ofstream mdat(mdat_name.c_str(), ios::binary | ios::out);

    int offset = 0;
    out.write(reinterpret_cast<const char*>(&num_nodes), sizeof(int));
    offset+=4;
    out.write(reinterpret_cast<const char*>(&num_edges), sizeof(int));
    offset+=4;

    for (int i = 0; i < num_nodes; i++) {
        mdat.write(reinterpret_cast<const char*>(&offset), sizeof(int));
        int node_i = i;
        int deg_i = graph[i].size();

        out.write(reinterpret_cast<const char*>(&node_i), sizeof(int));
        offset+=4;
        out.write(reinterpret_cast<const char*>(&deg_i), sizeof(int));
        offset+=4;

        for (int j = 0; j < deg_i; j++) {
            int neighbor_j = graph[i][j];
            int neighbor_j_size = graph[neighbor_j].size();
            out.write(reinterpret_cast<const char*>(&neighbor_j), sizeof(int));
            offset+=4;
//            out.write(reinterpret_cast<const char*>(&neighbor_j_size), sizeof(int));
//            offset+=4;

        }
    }
    out.close();
    return 0;
}
