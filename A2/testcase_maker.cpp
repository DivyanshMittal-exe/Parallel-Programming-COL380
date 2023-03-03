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
        cerr << "Usage: " << argv[0] << " <output file>" << endl;
        return 1;
    }

    vector<vector<int>> graph = {{1, 2}, {0, 2}, {0, 1}};

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
    ofstream out(filename, ios::binary | ios::out);
    out.write(reinterpret_cast<const char*>(&num_nodes), sizeof(int));
    out.write(reinterpret_cast<const char*>(&num_edges), sizeof(int));
    for (int i = 0; i < num_nodes; i++) {
        int node_i = i;
        int deg_i = graph[i].size();
        out.write(reinterpret_cast<const char*>(&node_i), sizeof(int));
        out.write(reinterpret_cast<const char*>(&deg_i), sizeof(int));
        for (int j = 0; j < deg_i; j++) {
            int neighbor_j = graph[i][j];
            int neighbor_j_size = graph[neighbor_j].size();
            out.write(reinterpret_cast<const char*>(&neighbor_j), sizeof(int));
            out.write(reinterpret_cast<const char*>(&neighbor_j_size), sizeof(int));
        }
    }
    out.close();
    return 0;
}
