#include <iostream>
#include <algorithm>

#include <vector>
#include <mpi/mpi.h>

#include <fstream>
#include <string.h>

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
        return ((u == other.u) && (v == other.v) && (w == other.w))
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

    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);


    MPI_File input_data;
    MPI_File_open(MPI_COMM_WORLD, argv[2], MPI_MODE_RDONLY, MPI_INFO_NULL, &input_data);

    MPI_File meta_data;
    MPI_File_open(MPI_COMM_WORLD, argv[3], MPI_MODE_RDONLY, MPI_INFO_NULL, &meta_data);

    int n, m;
    MPI_File_read_at_all(input_data, 0, &n, 1, MPI_INT, MPI_STATUS_IGNORE);
    MPI_File_read_at_all(input_data, 4, &m, 1, MPI_INT, MPI_STATUS_IGNORE);

    int avg_work_chunk = n/size + (n%size != 0);
    int my_chunk_start = 4 * rank * avg_work_chunk;

    int my_chunk_offsets[avg_work_chunk];
    memset(my_chunk_offsets,0,4*avg_work_chunk);
    MPI_File_read_at_all(meta_data, my_chunk_start, my_chunk_offsets, avg_work_chunk, MPI_INT , MPI_STATUS_IGNORE);
    

    vector<node> graph;

    for(int i  = 0; i < avg_work_chunk; i++){
        if(my_chunk_offsets[i] == 0)
            break;
        
        int node_val,node_deg;
        MPI_File_read_at_all(input_data, my_chunk_offsets[i], &node_val, 1, MPI_INT , MPI_STATUS_IGNORE);
        MPI_File_read_at_all(input_data, my_chunk_offsets[i] + 4, &node_deg, 1, MPI_INT , MPI_STATUS_IGNORE);


        auto new_node = node(node_deg,node_val);
        MPI_File_read_at_all(input_data, my_chunk_offsets[i] + 8, new_node.neighbours.data(), 2*node_deg, MPI_INT , MPI_STATUS_IGNORE);

        graph.push_back(new_node);

    }
    // Uptill here graph reading for each node



    MPI_Datatype MPI_NODE_VAL_ONLY_ARRAY_MINE;

    //   I think this is smart, but might be something dumb as well. Lol, can break IDK
    MPI_Type_vector(avg_work_chunk, 4, sizeof (node), MPI_CHAR, &MPI_NODE_VAL_ONLY_ARRAY_MINE);
    MPI_Type_commit(&MPI_NODE_VAL_ONLY_ARRAY_MINE);


    int vertex_node[size*avg_work_chunk];
    memset(vertex_node,0,4*size*avg_work_chunk);


    // Now this line helps gather info who has which node
    MPI_Allgather(graph.data(), 1, MPI_NODE_VAL_ONLY_ARRAY_MINE, vertex_node, avg_work_chunk, MPI_INT, MPI_COMM_WORLD);



    // Mapping of nodes
    map<int,int> node_mapping;


    // Here iteration only till n
    for (int i = 0; i < n; ++i) {
        node_mapping[vertex_node[i]] = i/avg_work_chunk;
        // 0 to avg_work_chunk-1 => 0, avg_work_chunk to 2*avg_work_chunk-1 => 1... Seems correct
    }


    // Now I start making Tau Hat
    map<pair<int,int>,int> tau_hat_e;


    // Here initialisation by 2, as tau_hat = supp + 2
    for(const auto &node_g: graph){
        for(const auto &node_g_n: node_g.neighbours){

            tau_hat_e[{node_g.node_val,node_g_n.node_val}] = 2;

        }
    }



    // Tau Hat queies stored here
    vector<vector<query_struct>> queries_for_tau(size);

    // Here aggregate the total number of queries to be made, maybe have a struct of queries or something
    for(const auto &node_g: graph){
        for(const auto &node_g_n_1: node_g.neighbours){
            for(const auto &node_g_n_2: node_g.neighbours){
                if(node_g_n_1!=node_g_n_2)
                    queries_for_tau[node_mapping[node_g_n_1]].push_back({node_g_n_1,node_g_n_2});
            }
        }
    }



    // Scrap this, gatherv takes care of it, as recv_counts[i] = 0, if nothing to send

    // Remember Last query is {-1,-1} so that everyone sends atleast 1 query
//    for(auto &ind: queries_for_tau){
//        ind.push_back({-1,-1});
//    }

    // Gathers how much data each guy will send
    vector<int> recv_counts(size);

    // Each node does it so gather for this node only if i == rank. Recv parameters only relevant then
    for (int i = 0; i < size; ++i) {
        int send_size = queries_for_tau[i].size();
        MPI_Gather(&send_size, 1, MPI_INT, recvcounts.data(), 1, MPI_INT, i, MPI_COMM_WORLD);
    }

    vector<int> displs(size);
    displs[0] = 0;
    for (int i = 1; i < size; i++) {
        displs[i] = displs[i-1] + recvcounts[i-1];
    }

    int total_queries = displs[size-1] + recvcounts[size-1];
    vector<query_struct> recv_buf(total_queries);

    // Query Struct Datatype
    MPI_Datatype MPI_QUERY_STRUCT;

    int block_lengths[2] = {1, 1};
    MPI_Aint displacements[2] = {0, sizeof(int)};
    MPI_Datatype types[2] = {MPI_INT, MPI_INT};

    MPI_Type_create_struct(2, block_lengths, displacements, types, &MPI_QUERY_STRUCT);
    MPI_Type_commit(&MPI_QUERY_STRUCT);

    // MPI Gathering all queries


//    for (int i = 0; i < size; ++i) {
//        MPI_Gatherv(queries_for_tau[i].data(), queries_for_tau[i].size(), MPI_QUERY_STRUCT, recv_buf.data(), recv_counts.data(), displs.data(), MPI_QUERY_STRUCT, i, MPI_COMM_WORLD);
//    }

    // Check if this Allgatherv is correct
    MPI_Allgatherv(queries_for_tau[rank].data(), queries_for_tau[rank].size(), MPI_QUERY_STRUCT,
                   recv_buf.data(), recv_counts.data(), displs.data(), MPI_QUERY_STRUCT, MPI_COMM_WORLD);



    // Here we get final tau_hat = supp + 2, as 2 was initialised in the start
    vector<char> ret_for_tri(total_queries,0);


    for(int i = 0; i < total_queries; i++){
        int u = recv_buf[i].u;
        int v = recv_buf[i].v;
        if(tau_hat_e.count({u,v})){
            tau_hat_e[{u,v}]++;
            ret_for_tri[i] = 1;
        }
    }


    // Now the code to return existance of edge back
    vector<int> displs_for_tri_recv(size + 1);
    vector<int> count_for_tri_recv(size);
    displs_for_tri_recv[0] = 0;
    for(int i = 1; i < size+1; i++){

        count_for_tri_recv[i-1] = queries_for_tau[i-1].size();
        displs_for_tri_recv[i] = displs_for_tri_recv[i-1] + queries_for_tau[i-1].size();
    }

    vector<char> tri_recv_buf(displs_for_tri_recv[size]);




    for (int i = 0; i < size; ++i) {
        MPI_Scatterv(ret_for_tri.data(),recv_counts.data(), displs.data(),MPI_CHAR, &tri_recv_buf[displs_for_tri_recv[i]], queries_for_tau[i].size(), MPI_CHAR, i, MPI_COMM_WORLD);
    }

    // Check this as well
    // Needs obv changes
//    MPI_Alltoallv(ret_for_tri.data(), count_for_tri_recv.data(), displs_for_tri_recv.data(), MPI_CHAR,
//                  tri_recv_buf.data(), recv_counts.data(), displs_for_tri_recv.data(), MPI_CHAR,
//                  MPI_COMM_WORLD);


    map<Tri_Struct,int> tau_hat_tri;

    int q_it = 0;

    for(const auto &node_g: graph){
        for(const auto &node_g_n_1: node_g.neighbours){
            for(const auto &node_g_n_2: node_g.neighbours){
                if(node_g_n_1!=node_g_n_2) {
                    if (tri_recv_buf[q_it] == 1) {
                        tau_hat_tri[Tri_Struct(node_g.node_val,node_g_n_1.node_val,node_g_n_2.node_val)] = INF;
                    }
                    q_it ++;
                }
            }
        }
    }


    map<pair<int,int>,int> g_e;
    map<pair<pair<int,int>,int>,int> h_e_j;


    for (const auto& [key, value] : tau_hat_e) {
        g_e[key] = value-2;
        for (int i = 0; i < value; ++i) {
            h_e_j[{key,i}] = 0;
         }
    }

    // Intialisation ends here

    auto k_min_it = min_element(tau_hat_e.begin(), tau_hat_e.end(),[](const auto& a, const auto& b) { return a.second < b.second; });
    auto k_max_it = max_element(tau_hat_e.begin(), tau_hat_e.end(),[](const auto& a, const auto& b) { return a.second < b.second; });
    int k_min_loc = k_min_it->second
    int k_max_loc = k_max_it->second

    int k_min_global;
    int k_max_global;
    MPI_Allreduce(&k_min_loc, &k_min_global, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&k_max_loc, &k_max_global, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);


    return 0;
}
