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
    node_neighbour* neighbours;
//    vector<node_neighbour> neighbours;

    node(){}
    node(int val, int deg) : node_val(val), node_deg(deg){
        neighbours = (node_neighbour*)calloc(node_deg,sizeof(node_neighbour));
    }
//    node(int val, int deg) : node_val(val), node_deg(deg), neighbours(deg) {
//    }
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
        MPI_File_read_at_all(input_data, my_chunk_offsets[i] + 8, new_node.neighbours, 2*node_deg, MPI_INT , MPI_STATUS_IGNORE);

        graph.push_back(new_node);

    }



    // ifstream input(argv[1], ios::binary);


    // input.read((char *) &n, 4);
    // input.read((char *) &m, 4);

    // vector<node> graph(n);

    // for (int i = 0; i < n; ++i) {
    //     int node_val = 0,node_deg =0;
    //     input.read((char *) &node_val, 4);
    //     input.read((char *) &node_deg, 4);
    //     graph[i] = node(node_deg,node_val);

    //     input.read((char *) &(graph[i].neighbours), 4*2*node_deg);

    // }


    return 0;
}
