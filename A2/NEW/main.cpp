#include <iostream>
#include <algorithm>
#include <cassert>
#include <vector>
#include <mpi/mpi.h>
#include <set>
#include <climits>
#include <fstream>
#include <string.h>

#define DEBUG_MODE 1
#define EXP_MODE 0

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

    bool operator==(const query_struct& other) const {
        return ((u == other.u) && (v == other.v));
    }
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

    MPI_Datatype MPI_QUERY_STRUCT;

    int block_lengths[2] = {1, 1};
    MPI_Aint displacements[2] = {0, sizeof(int)};
    MPI_Datatype types[2] = {MPI_INT, MPI_INT};

    MPI_Type_create_struct(2, block_lengths, displacements, types, &MPI_QUERY_STRUCT);
    MPI_Type_commit(&MPI_QUERY_STRUCT);



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

//    vector<node> graph;

    map<int,node> graph;

    for (int i = rank; i < n; i+= size) {
        int node_val,node_deg;

        MPI_File_read_at_all(input_data, all_offsets[i], &node_val, 1, MPI_INT , MPI_STATUS_IGNORE);
        MPI_File_read_at_all(input_data, all_offsets[i] + 4, &node_deg, 1, MPI_INT , MPI_STATUS_IGNORE);


        auto new_node = node(node_val,node_deg);
        MPI_File_read_at_all(input_data, all_offsets[i] + 8, new_node.neighbours.data(), 2*node_deg, MPI_INT , MPI_STATUS_IGNORE);

//        graph.push_back(new_node);

        graph[node_val] = new_node;
    }

#if DEBUG_MODE
    cout << "Read graph " << graph.size() << endl;
#endif

    map<pair<int,int>,vector<int>> edge_to_third_node_map;

    map<pair<int,int>,pair<int,char>> tau_hat_e;

    for(const auto &node_g: graph){
        for(const auto &node_g_n_1: node_g.second.neighbours){
            tau_hat_e[{node_g.second.node_val,node_g_n_1.node_val}] = {2,0};
            for(const auto &node_g_n_2: node_g.second.neighbours){
                if(node_g_n_1.node_val < node_g_n_2.node_val){
                    edge_to_third_node_map[{node_g_n_1.node_val,node_g_n_2.node_val}].push_back(node_g.second.node_val);
                }
//                    queries_for_tau[node_mapping[node_g_n_1]].push_back({node_g_n_1,node_g_n_2});
            }
        }
    }

    map<int,int> node_mapping;


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
                        tau_hat_e[{x,n_val_temp}].first++;
                        tau_hat_e[{x,n_node.node_val}].first++;
                        tau_hat_tri.insert({x,n_val_temp,n_node.node_val});
                    }
                }
            }
        }
    }

#if DEBUG_MODE
    for(const auto & tau_it: tau_hat_e){
        cout << tau_it.first.first << " " tau_it.first.second << " " << tau_it.second.first << " " << tau_it.second.second << endl;
    }
    for(const auto& tri: tau_hat_tri){
        cout << tri.u << " " << tri.v << " " << tri.w;
    }
//    cout << "Read graph " << graph.size() << endl;
#endif

//    int MPI_Alltoallv(const void *sendbuf, const int *sendcounts,
//                      const int *sdispls, MPI_Datatype sendtype, void *recvbuf,
//                      const int *recvcounts, const int *rdispls, MPI_Datatype recvtype,
//                      MPI_Comm comm)

#if EXP_MODE


    vector<int> query_recvcounts(size,0);
    for(const auto &node_g: graph){
        for(const auto &node_g_n_1: node_g.second.neighbours){
            query_recvcounts[node_mapping[node_g_n_1.node_val]]++;
        }
    }
    vector<int> query_rdispls(size,0);
    for(int i = 1; i < size; i++){
        query_rdispls[i] = query_rdispls[i-1] + query_recvcounts[i-1];
    }

    vector<query_struct> query_recv_buf(query_rdispls[size-1] + query_recvcounts[size-1]);





    vector<vector<query_struct>> decrement_queries(size);

//    auto k_min_it = min_element(tau_hat_e.begin(), tau_hat_e.end(),[](const auto& a, const auto& b) { return a.second < b.second; });

    auto k_min_it = std::min_element(tau_hat_e.begin(), tau_hat_e.end(),
                                    [](const auto& a, const auto& b) {
                                        if (a.second.second == 0 && b.second.second == 0) {
                                            return a.second.first < b.second.first;
                                        } else if (a.second.second == 0) {
                                            return true;
                                        } else {
                                            return false;
                                        }
                                    });

    assert((k_min_it->second).second == 0);

    int k_min_loc = (k_min_it->second).first;


    int k_min_global;
    MPI_Allreduce(&k_min_loc, &k_min_global, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

    for (auto const& pr : tau_hat_e) {
        auto const& edge = pr.first;
        auto const& tau_val = pr.second;


        if(tau_val.first == k_min_global && tau_val.second == 0){

            auto node_to_enum = graph[edge.first];
            for(const auto &neighb: node_to_enum.neighbours){
                // If a triangle exists
                if(neighb.node_val != edge.second && tau_hat_tri.count({edge.first,edge.second,neighb.node_val})){
                    auto &under_consider = tau_hat_e[{edge.first,neighb.node_val}];
//                     Case where both are min
                    if(under_consider.first == k_min_global && under_consider.second == 0){
                        // Only 1 of them sends
                        if(edge.second < neighb.node_val){
                            decrement_queries[node_mapping[edge.second]].push_back({edge.second,neighb.node_val});
                            decrement_queries[node_mapping[neighb.node_val]].push_back({neighb.node_val,edge.second});

                        }
                    }else{
                        decrement_queries[node_mapping[edge.second]].push_back({edge.second,neighb.node_val});
                        decrement_queries[node_mapping[neighb.node_val]].push_back({neighb.node_val,edge.second});

                    }
                }
            }
        }
    }


    for(auto &x:decrement_queries){
        x.push_back({-1,-1});
    }

    for (auto & pr : tau_hat_e) {
        auto const& edge = pr.first;
        auto & tau_val = pr.second;
        if(tau_val.first == k_min_global && tau_val.second == 0){
            tau_val.second = 1;
        }
    }




//    int MPI_Alltoallv(const void *sendbuf, const int *sendcounts,
//                      const int *sdispls, MPI_Datatype sendtype, void *recvbuf,
//                      const int *recvcounts, const int *rdispls, MPI_Datatype recvtype,
//                      MPI_Comm comm)

    vector<int> query_sendcounts(size,0);
    vector<int> query_sdispls(size,0);



    for (int i = 1; i < size; ++i) {
        query_sendcounts[i-1] = decrement_queries[i-1].size();
        query_sdispls[i] = query_sdispls[i-1] + query_sendcounts[i-1];
    }
    query_sendcounts[size-1] = decrement_queries[size-1].size();


    vector<query_struct> all_together;
    all_together.reserve(query_sendcounts[size-1]  + decrement_queries[size-1].size() );
    for (auto& dec_q_node: decrement_queries){
        move(dec_q_node.begin(), dec_q_node.end(), std::back_inserter(all_together));
    }


    MPI_Alltoallv(all_together.data(),query_sendcounts.data(),
                  query_sdispls.data(),MPI_QUERY_STRUCT,query_recv_buf.data(),
                  query_recvcounts.data(),query_rdispls.data(),MPI_QUERY_STRUCT,
                  MPI_COMM_WORLD);


    for(int i = 0;i < size;i++){
        for (int j = query_rdispls[i]; j < query_rdispls[i] + query_recvcounts[i]; ++j) {
            if(query_recv_buf[j].u == -1){
                break;
            }

            auto &edge_to_pot_dec = tau_hat_e[{query_recv_buf[j].u,query_recv_buf[j].v}];

            if(edge_to_pot_dec.second == 0){
                edge_to_pot_dec.first--;
            }
        }
    }

#endif

    MPI_Finalize();

    return 0;
}
