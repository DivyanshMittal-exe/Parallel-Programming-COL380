#include <iostream>
#include <algorithm>
#include <cassert>
#include <vector>
#include <mpi/mpi.h>
#include <set>
#include <climits>
#include <fstream>
#include <string.h>
#include <queue>
#include <map>
#include <chrono>

#define DEBUG_MODE 1
#define EXP_MODE 1
#define VERBOSE_ONE_MODE_CODE 1

#if DEBUG_MODE
#define DEBUG_STAT //cout << "Here by " << rank << endl;
#else
#define DEBUG_STAT ;
#endif

using namespace std;

const int INF = INT_MAX;


struct HASH_FOR_PAIR {
    template <class T1, class T2>
    size_t operator()(const pair<T1, T2>& p) const
    {
        return p.first >= p.second ? p.first * p.first + p.first + p.second : p.second * p.second + p.first;

    }
};

struct Tri_Struct {
    int u, v, w;

    Tri_Struct(int u, int v, int w) {
        int vertices[] = {u, v, w};
        sort(vertices, vertices + 3);
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



struct tau_edge{
    int u,v,tau;
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
};


struct node {
    int node_val;
    int node_deg;
    vector<node_neighbour> neighbours;

    node(){}

    node(int val, int deg) : node_val(val), node_deg(deg), neighbours(deg) {
    }
};



int main(int argc, char *argv[]) {



    auto start = chrono::steady_clock::now();


    MPI_Init(&argc, &argv);
    int task_id = -1;
    string input_path;
    string header_path;
    string output_path;
    int verbose = -1;
    int start_k = -1;
    int end_k = -1;
    int p = -1;

    for (int i = 1; i < argc; i++) {
        string arg = argv[i];

        if (arg.find("--taskid=") != string::npos) {
            task_id = stoi(arg.substr(arg.find("=") + 1));
        } else if (arg.find("--inputpath=") != string::npos) {
            input_path = arg.substr(arg.find("=") + 1);
        } else if (arg.find("--headerpath=") != string::npos) {
            header_path = arg.substr(arg.find("=") + 1);
        } else if (arg.find("--outputpath=") != string::npos) {
            output_path = arg.substr(arg.find("=") + 1);
        } else if (arg.find("--verbose=") != string::npos) {
            verbose = stoi(arg.substr(arg.find("=") + 1));
        } else if (arg.find("--startk=") != string::npos) {
            start_k = stoi(arg.substr(arg.find("=") + 1));
        } else if (arg.find("--endk=") != string::npos) {
            end_k = stoi(arg.substr(arg.find("=") + 1));
        } else if (arg.find("--p=") != string::npos) {
            p = stoi(arg.substr(arg.find("=") + 1));
        }
    }



    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    MPI_Datatype MPI_QUERY_STRUCT;
    MPI_Datatype type[2] = { MPI_INT, MPI_INT };
    int blocklen[2] = { 1, 1 };
    MPI_Aint disp[2];
    disp[0] = offsetof(query_struct, u);
    disp[1] = offsetof(query_struct, v);

    MPI_Type_create_struct(2, blocklen, disp, type, &MPI_QUERY_STRUCT);
    MPI_Type_commit(&MPI_QUERY_STRUCT);


    MPI_Datatype MPI_TAU_EDGE;
    int blocklengths[3] = { 1, 1, 1 };
    MPI_Datatype types[3] = { MPI_INT, MPI_INT, MPI_INT };
    MPI_Aint offsets[3];
    offsets[0] = offsetof(struct tau_edge, u);
    offsets[1] = offsetof(struct tau_edge, v);
    offsets[2] = offsetof(struct tau_edge, tau);
    MPI_Type_create_struct(3, blocklengths, offsets, types, &MPI_TAU_EDGE);
    MPI_Type_commit(&MPI_TAU_EDGE);




    MPI_File input_data;
    MPI_File_open(MPI_COMM_WORLD, input_path.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &input_data);

    MPI_File meta_data;
    MPI_File_open(MPI_COMM_WORLD, header_path.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &meta_data);

    int n, m;
    MPI_File_read_at_all(input_data, 0, &n, 1, MPI_INT, MPI_STATUS_IGNORE);
    MPI_File_read_at_all(input_data, 4, &m, 1, MPI_INT, MPI_STATUS_IGNORE);






    int all_offsets[n];
    MPI_File_read_at_all(meta_data, 0, all_offsets, n, MPI_INT , MPI_STATUS_IGNORE);


    vector<node> graph;

    for (int i = rank; i < n; i+= size) {
        int node_val,node_deg;

        MPI_File_read_at(input_data, all_offsets[i], &node_val, 1, MPI_INT , MPI_STATUS_IGNORE);
        MPI_File_read_at(input_data, all_offsets[i] + 4, &node_deg, 1, MPI_INT , MPI_STATUS_IGNORE);


        auto new_node = node(node_val,node_deg);
        MPI_File_read_at(input_data, all_offsets[i] + 8, new_node.neighbours.data(), node_deg, MPI_INT , MPI_STATUS_IGNORE);

        graph.push_back(new_node);
    }



    map<pair<int,int>,pair<int,int>> tau_hat_e;




    for(const auto &node_g: graph){
        for(const auto &node_g_n_1: node_g.neighbours){
            tau_hat_e[{node_g.node_val,node_g_n_1.node_val}] = {0,0};
        }
    }





//    map<pair<int,int>,int> tau_hat_tri;

    vector<vector<int>> tau_hat_tri(n);



    for(int i = 0; i < n; i++){


        int n_val_temp,n_deg_temp;

        MPI_File_read_at_all(input_data, all_offsets[i], &n_val_temp, 1, MPI_INT , MPI_STATUS_IGNORE);

        MPI_File_read_at_all(input_data, all_offsets[i] + 4, &n_deg_temp, 1, MPI_INT , MPI_STATUS_IGNORE);



        vector<node_neighbour> temp(n_deg_temp);
        MPI_File_read_at_all(input_data, all_offsets[i] + 8, temp.data(), n_deg_temp, MPI_INT , MPI_STATUS_IGNORE);



        for(const auto &n_node_1: temp){
            if(n_node_1.node_val % size == rank){
                auto &other = temp;
                auto &me = graph[n_node_1.node_val/size].neighbours;

                auto it_1 = me.begin();
                auto it_2 = other.begin();

                while(it_1 != me.end() && it_2 != other.end()){
                    if(it_1->node_val == it_2->node_val){
                        tau_hat_e[{n_node_1.node_val,n_val_temp}].first++;

                        tau_hat_tri[n_val_temp].push_back(it_2->node_val);

//                        tau_hat_tri[{n_val_temp,it_2->node_val}] = 0;
                        it_2++;
                        it_1++;
                    }else if(it_1->node_val < it_2->node_val){
                        it_1++;
                    }else{
                        it_2++;
                    }
                }
            }
        }


    }


    int k = 3;

    while (true) {


        int my_edges = 0;
        for (auto it = tau_hat_e.begin(); it != tau_hat_e.end(); ++it) {
            const auto &a = it->second;

            if (a.second == 0) {
                my_edges += 1;
                break;
            }
        }


        int all_edges;
        MPI_Allreduce(&my_edges, &all_edges, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

        if (all_edges == 0) {
            break;
        }


        while (true) {
            vector <vector<query_struct>> decrement_queries(size);

            int f_k_count = 0;
            for (auto const &pr: tau_hat_e) {
                auto const &edge = pr.first;
                auto const &tau_val = pr.second;
                if (tau_val.first < k - 2 && tau_val.second == 0) {
                    f_k_count++;
                    break;
                }
            }


            int f_k_count_global;
            MPI_Allreduce(&f_k_count, &f_k_count_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

            if (f_k_count_global == 0) {
                break;
            }


            for (auto const &pr: tau_hat_e) {
                auto const &edge = pr.first;
                auto const &tau_val = pr.second;


                if (tau_val.first < k - 2 && tau_val.second == 0) {

                    auto node_to_enum = graph[edge.first / size];
                    for (const auto &neighb: node_to_enum.neighbours) {

                        if (neighb.node_val != edge.second &&
                            std::binary_search(tau_hat_tri[edge.second].begin(), tau_hat_tri[edge.second].end(),
                                               neighb.node_val)
                                ) {
                            auto &under_consider = tau_hat_e[{edge.first, neighb.node_val}];
//                     Case where both are min
                            if (under_consider.first < k - 2 && under_consider.second == 0) {
                                // Only 1 of them sends
                                if (edge.second < neighb.node_val) {
                                    decrement_queries[edge.second % size].push_back({edge.second, neighb.node_val});
                                    decrement_queries[neighb.node_val % size].push_back(
                                            {neighb.node_val, edge.second});

                                }
                            } else {
                                if (under_consider.second == 0) {
                                    decrement_queries[edge.second % size].push_back({edge.second, neighb.node_val});
                                    decrement_queries[neighb.node_val % size].push_back(
                                            {neighb.node_val, edge.second});

                                }

                            }
                        }
                    }
                }
            }


            for (auto &pr: tau_hat_e) {
                auto const &edge = pr.first;
                auto &tau_val = pr.second;
                if (tau_val.first < k - 2 && tau_val.second == 0) {
                    tau_val.first = k - 1;
                    tau_val.second = 1;
                }
            }


            vector<int> query_sendcounts(size, 0);
            vector<int> query_sdispls(size, 0);


            for (int i = 1; i < size; ++i) {
                query_sendcounts[i - 1] = decrement_queries[i - 1].size();
                query_sdispls[i] = query_sdispls[i - 1] + query_sendcounts[i - 1];
            }
            query_sendcounts[size - 1] = decrement_queries[size - 1].size();


            vector <query_struct> all_together;
            all_together.reserve(query_sdispls[size - 1] + query_sendcounts[size - 1]);
            for (auto &dec_q_node: decrement_queries) {
                move(dec_q_node.begin(), dec_q_node.end(), back_inserter(all_together));
            }


            vector<int> query_recvcounts(size, 0);

            MPI_Alltoall(query_sendcounts.data(), 1, MPI_INT, query_recvcounts.data(), 1, MPI_INT, MPI_COMM_WORLD);

            DEBUG_STAT


            vector<int> query_rdispls(size, 0);

            vector <query_struct> query_recv_buf(query_rdispls[size - 1] + query_recvcounts[size - 1]);

            for (int i = 1; i < size; i++) {
                query_rdispls[i] = query_rdispls[i - 1] + query_recvcounts[i - 1];
            }


            MPI_Alltoallv(all_together.data(), query_sendcounts.data(),
                          query_sdispls.data(), MPI_QUERY_STRUCT, query_recv_buf.data(),
                          query_recvcounts.data(), query_rdispls.data(), MPI_QUERY_STRUCT,
                          MPI_COMM_WORLD);


            for (auto q: query_recv_buf) {

                auto &edge_to_pot_dec = tau_hat_e[{q.u, q.v}];

                if (edge_to_pot_dec.second == 0) {
                    edge_to_pot_dec.first--;
                }
            }
        }

            k++;
        }




    if(verbose == 0){
        int max_k_loc = 0;
        int max_k_glob;
        for(const auto &tau_it: tau_hat_e){
            max_k_loc = max(max_k_loc,tau_it.second.first);
        }


        MPI_Reduce(&max_k_loc,&max_k_glob,1,MPI_INT,MPI_MAX,0,MPI_COMM_WORLD);

        max_k_glob -= 2;

        if(rank == 0) {
            string buff = "";
            for(int i = start_k; i <= end_k; i++){
                if(i > max_k_glob){
                    buff += "0\n";
                }else{
                    buff += "1\n";
                }
            }

            MPI_File outfile;
            MPI_File_open(MPI_COMM_SELF, output_path.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &outfile);
            MPI_File_write(outfile, buff.c_str(), buff.size(), MPI_CHAR, MPI_STATUS_IGNORE);
            MPI_File_close(&outfile);
        }

    }else {


        vector <tau_edge> settled_edges;

        for (const auto &tau_it: tau_hat_e) {
            auto const &edge = tau_it.first;
            auto const &tau_val = tau_it.second;

            if (edge.first < edge.second) {
                auto new_tau_edge = tau_edge{edge.first, edge.second, tau_val.first};
                settled_edges.push_back(new_tau_edge);
            }

        }

        int my_settle_size = settled_edges.size();
        vector<int> settled_edge_recv_cnt(size, 0);

        MPI_Allgather(&my_settle_size, 1, MPI_INT, settled_edge_recv_cnt.data(), 1, MPI_INT, MPI_COMM_WORLD);

        vector<int> settled_edge_rdispls(size, 0);

        for (int i = 1; i < size; i++) {
            settled_edge_rdispls[i] = settled_edge_rdispls[i - 1] + settled_edge_recv_cnt[i - 1];
        }


        vector <tau_edge> all_edges(m);


        MPI_Allgatherv(settled_edges.data(), settled_edges.size(), MPI_TAU_EDGE,
                       all_edges.data(), settled_edge_recv_cnt.data(), settled_edge_rdispls.data(),
                       MPI_TAU_EDGE, MPI_COMM_WORLD);


        vector<vector<pair<int,int>>> final_graph(n);

        for(const auto& tedge: all_edges){
            final_graph[tedge.u].push_back({tedge.v,tedge.tau});
            final_graph[tedge.v].push_back({tedge.u,tedge.tau});
        }


        int diff_in_k = end_k - start_k + 1;
        int num_it = diff_in_k/size + (diff_in_k%size != 0);

        int i_start_from = start_k + rank - size;

        MPI_File outfile;
        MPI_File_open(MPI_COMM_WORLD, output_path.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &outfile);

        for(int i = 0; i < num_it; i++){
            i_start_from += size;
            string buff = "";
            if(i_start_from <= end_k){
                int theta = 0;
                vector<int> con_comp(n,-1);
                for(int i = 0; i < n; i++){
                    if(con_comp[i] == -1){
                        deque<int> bfs_q;
                        int flag = 0;
                        for(auto n: final_graph[i]){
                            if(n.second >= i_start_from + 2){
                                flag = 1;
                                break;
                            }
                        }
                        if(!flag)
                            continue;
                        bfs_q.push_back(i);
                        while (!bfs_q.empty()){
                            auto e = bfs_q.front();
                            bfs_q.pop_front();
                            con_comp[e] = theta;
                            for(auto node: final_graph[e]){
                                if(node.second >= i_start_from + 2 && con_comp[node.first] == -1){
                                    bfs_q.push_back(node.first);
                                }
                            }
                        }
                        theta++;
                    }
                }
                if(theta != 0)
                {
                    buff+="1\n";
                }
                buff += to_string(theta);
                buff += "\n";
                for (int j = 0; j < theta ; ++j) {
                    for (int l = 0; l < con_comp.size(); ++l) {
                        if(con_comp[l] == j){
                            buff += to_string(l);
                            buff += " ";
                        }
                    }
                    buff[buff.size()-1] = '\n';
                }
            }

            MPI_File_write_ordered(outfile,buff.c_str(),buff.size(),MPI_CHAR,MPI_STATUS_IGNORE);
            MPI_Barrier(MPI_COMM_WORLD);
        }

        MPI_File_close(&outfile);


    }
    MPI_Finalize();

    auto end = chrono::steady_clock::now();
    cout << "Elapsed time in milliseconds: "
         << chrono::duration_cast<chrono::milliseconds>(end - start).count()
         << " ms" << endl;

    return 0;
}