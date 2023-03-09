#include <iostream>
#include <algorithm>

#include <vector>
#include <mpi/mpi.h>
#include <set>
#include <climits>
#include <fstream>
#include <string.h>

#define DEBUG_MODE 1

using namespace std;

const int INF = INT_MAX;

struct Tri_Struct {
    int u, v, w;

    Tri_Struct(int u, int v, int w) {
        int vertices[] = {u, v, w};
        std::sort(vertices, vertices + 3);
        this->u = vertices[0];
        this->v = vertices[1];
        this->w = vertices[2];
    }

    bool operator<(const Tri_Struct& other) const {
        if (u != other.u) {
            return u < other.u;
        } else if (v != other.v) {
            return v < other.v;
        } else {
            return w < other.w;
        }
    }

    bool operator==(const Tri_Struct& other) const {
        return ((u == other.u) && (v == other.v) && (w == other.w));
    }
};


struct query_struct{
    int u;
    int v;
};

struct node_neighbour {
    int node_val;
    int node_deg;
};


struct node {
    int node_val;
    int node_deg;
//    node_neighbour* neighbours;
    vector<node_neighbour> neighbours;

    node(){}
//    node(int val, int deg) : node_val(val), node_deg(deg){
//        neighbours = (node_neighbour*)calloc(node_deg,sizeof(node_neighbour));
//    }
    node(int val, int deg) : node_val(val), node_deg(deg), neighbours(deg) {
    }
};



int main(int argc, char *argv[]) {


#if DEBUG_MODE
    cout << "HELLO" << endl;
#endif


    MPI_Init(&argc, &argv);




    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);


    char inp_f_name[] = "test_case/test0/test-input-0.gra";
    char mdt_f_name[] = "test_case/test0/test-header-0.dat";


#if DEBUG_MODE
    cout << "Initialised" << endl;
#endif


    MPI_File input_data;
    MPI_File_open(MPI_COMM_WORLD, inp_f_name, MPI_MODE_RDONLY, MPI_INFO_NULL, &input_data);

    MPI_File meta_data;
    MPI_File_open(MPI_COMM_WORLD, mdt_f_name, MPI_MODE_RDONLY, MPI_INFO_NULL, &meta_data);

    int n, m;
    MPI_File_read_at_all(input_data, 0, &n, 1, MPI_INT, MPI_STATUS_IGNORE);
    MPI_File_read_at_all(input_data, 4, &m, 1, MPI_INT, MPI_STATUS_IGNORE);


#if DEBUG_MODE
    cout << "Read n,m " << n << " " << m << endl;
#endif



    int all_offsets[n];
    MPI_File_read_at_all(meta_data, 0, all_offsets, n, MPI_INT , MPI_STATUS_IGNORE);

    #if DEBUG_MODE
        cout << "Read offsets " << all_offsets[rank] << endl;
    #endif


    map<int,int> node_vertex_map;

    int my_num_nodes = (n/size) + (rank < n%size);

    vector<node> graph;


    for (int i = rank; i < n; i+= size) {
        int node_val,node_deg;

        MPI_File_read_at_all(input_data, all_offsets[i], &node_val, 1, MPI_INT , MPI_STATUS_IGNORE);
        MPI_File_read_at_all(input_data, all_offsets[i] + 4, &node_deg, 1, MPI_INT , MPI_STATUS_IGNORE);


        auto new_node = node(node_val,node_deg);
        MPI_File_read_at_all(input_data, all_offsets[i] + 8, new_node.neighbours.data(), 2*node_deg, MPI_INT , MPI_STATUS_IGNORE);

        graph.push_back(new_node);

    }

#if DEBUG_MODE
    cout << "Read graph " << graph.size() << endl;
#endif

    map<pair<int,int>,vector<int>> edge_to_third_node_map;

    for(const auto &node_g: graph){
        for(const auto &node_g_n_1: node_g.neighbours){
            for(const auto &node_g_n_2: node_g.neighbours){
                if(node_g_n_1.node_val < node_g_n_2.node_val){
                    edge_to_third_node_map[{node_g_n_1.node_val,node_g_n_2.node_val}].push_back(node_g.node_val);
                }
//                    queries_for_tau[node_mapping[node_g_n_1]].push_back({node_g_n_1,node_g_n_2});
            }
        }
    }

    map<int,int> node_mapping;


    map<pair<int,int>,int> tau_hat_e;
    set<Tri_Struct> tau_hat_tri;

    vector<node_neighbour> temp;

    for(int i = 0; i < n; i++){

        int n_val_temp,n_deg_temp;
        MPI_File_read_at_all(input_data, all_offsets[i], &n_val_temp, 1, MPI_INT , MPI_STATUS_IGNORE);

        MPI_File_read_at_all(input_data, all_offsets[i] + 4, &n_deg_temp, 1, MPI_INT , MPI_STATUS_IGNORE);

        node_mapping[n_val_temp] = i%size;

        temp.resize(n_deg_temp);
        MPI_File_read_at_all(input_data, all_offsets[i] + 8, temp.data(), 2*n_deg_temp, MPI_INT , MPI_STATUS_IGNORE);

        for(const auto &n_node: temp){
            if(n_val_temp < n_node.node_val){
                if(edge_to_third_node_map.count({n_val_temp,n_node.node_val})){
                    for(auto x: edge_to_third_node_map[{n_val_temp,n_node.node_val}]){
                        tau_hat_e[{x,n_val_temp}]++;
                        tau_hat_e[{x,n_node.node_val}]++;
                        tau_hat_tri.insert({x,n_val_temp,n_node.node_val});
                    }
                }
            }
        }
    }

//     for(auto x: tau_hat_e){
//         auto a = x.first;
//         cout << a.first  << " " << a.second <<" " << x.second << endl;
//     }

    MPI_Finalize();

    return 0;
}
