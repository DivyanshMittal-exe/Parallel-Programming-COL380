#include <iostream>
#include <algorithm>

#include <vector>
#include <mpi/mpi.h>

#include <fstream>
#include <string.h>

using namespace std;


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

    MPI_Datatype MPI_NODE_VAL_ONLY_ARRAY_MINE;

    //   I think this is smart, but might be something dumb as well. Lol, can break IDK
    MPI_Type_vector(avg_work_chunk, 4, sizeof (node), MPI_CHAR, &MPI_NODE_VAL_ONLY_ARRAY_MINE);
    MPI_Type_commit(&MPI_NODE_VAL_ONLY_ARRAY_MINE);


    int vertex_node[size*avg_work_chunk];
    memset(vertex_node,0,4*size*avg_work_chunk);


    MPI_Allgather(graph.data(), 1, MPI_NODE_VAL_ONLY_ARRAY_MINE, vertex_node, avg_work_chunk, MPI_INT, MPI_COMM_WORLD);
    
    map<int,int> node_mapping;


    // Here iteration only till n
    for (int i = 0; i < n; ++i) {
        node_mapping[vertex_node[i]] = i/avg_work_chunk;
        // 0 to avg_work_chunk-1 => 0, avg_work_chunk to 2*avg_work_chunk-1 => 1... Seems correct
    }

    map<pair<int,int>,int> tau_hat_e;

    for(const auto &node_g: graph){
        for(const auto &node_g_n: node_g.neighbours){

            tau_hat_e[{node_g.node_val,node_g_n.node_val}] = 0;

        }
    }

    // Here aggregate the total number of queries to be made, maybe have a struct of queries or something

    return 0;
}
