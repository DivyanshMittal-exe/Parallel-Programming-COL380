#include<bits/stdc++.h>
#include<mpi.h>

using namespace std;

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    MPI_Aint lb = 0, extent = world_size*sizeof(int);
    MPI_Datatype etype = MPI_INT, filetype, contig;
    MPI_Offset disp;

    MPI_Type_contiguous(1, MPI_INT, &contig);
    MPI_Type_create_resized(contig, lb, extent, &filetype);
    MPI_Type_commit(&filetype);

    disp = my_rank*sizeof(int);

    MPI_File input_data;
    MPI_File_open(MPI_COMM_WORLD, "test0/test-input-0.gra", MPI_MODE_RDONLY, MPI_INFO_NULL, &input_data);
    // MPI_File_open(MPI_COMM_WORLD, "/tmp/COL380/A2/test0/test-input-0.gra", MPI_MODE_RDONLY, MPI_INFO_NULL, &input_data);

    MPI_File meta_data;
    MPI_File_open(MPI_COMM_WORLD, "test0/test-header-0.dat", MPI_MODE_RDONLY, MPI_INFO_NULL, &meta_data);
    // MPI_File_open(MPI_COMM_WORLD, "/tmp/COL380/A2/test0/test-header-0.dat", MPI_MODE_RDONLY, MPI_INFO_NULL, &meta_data);

    MPI_File_set_view(meta_data, disp, etype, filetype, "native", MPI_INFO_NULL);

    int n, m;
    MPI_File_read_at_all(input_data, 0, &n, 1, MPI_INT, MPI_STATUS_IGNORE);
    MPI_File_read_at_all(input_data, 4, &m, 1, MPI_INT, MPI_STATUS_IGNORE);

    int my_num_nodes = (n/world_size) + (my_rank < n%world_size);

    int my_offsets[my_num_nodes];

    MPI_File_read_at_all(meta_data, 0, my_offsets, my_num_nodes, MPI_INT , MPI_STATUS_IGNORE);

    int my_vertex_ids[my_num_nodes];
    int my_vertex_deg[my_num_nodes];

    // neighbour id, neighbour degree, neighbour triangles
    vector<vector<int>> neighbour_data(my_num_nodes);

    // contains the coordinates in the neighbour data given a neighbour's value
    vector<map<int,int>> neighbour_coordinates(my_num_nodes);

    for(int i=0; i<my_num_nodes; i++)
    {
        MPI_File_read_at(input_data, my_offsets[i], &(my_vertex_ids[i]), 1, MPI_INT , MPI_STATUS_IGNORE);
        MPI_File_read_at(input_data, my_offsets[i]+4, &(my_vertex_deg[i]), 1, MPI_INT , MPI_STATUS_IGNORE);

        int neighbour_start = my_offsets[i]+8;

        for(int j=0; j<my_vertex_deg[i]; j++)
        {
            int neighbour;
            MPI_File_read_at(input_data, neighbour_start+4*j, &neighbour, 1, MPI_INT , MPI_STATUS_IGNORE);

            int neighbour_offset;
            MPI_File_read_at(meta_data, neighbour*4, &neighbour_offset, 1, MPI_INT, MPI_STATUS_IGNORE);

            int neighbour_degree;
            MPI_File_read_at(input_data, neighbour_offset+4, &neighbour_degree, 1, MPI_INT, MPI_STATUS_IGNORE);

            if(neighbour_degree > my_vertex_deg[i] || ((neighbour_degree == my_vertex_deg[i]) && (neighbour > my_vertex_ids[i])))
            {

                neighbour_data[i].push_back(neighbour);
                neighbour_data[i].push_back(neighbour_degree);
                neighbour_data[i].push_back(0);

                neighbour_coordinates[i].insert({neighbour, neighbour_data[i].size()-3});
            }
        }
    }

    MPI_File_close(&input_data);
    MPI_File_close(&meta_data);


    // own_node, node1 (with owner1), node2
    vector<vector<int>> queries(world_size);

    for(int i=0; i<my_num_nodes; i++)
    {
        int id = my_vertex_ids[i];
        int degree = my_vertex_deg[i];
        if(degree <2)
        {
            continue;
        }
        for(int j=0; j<degree; j++)
        {
            int deg_j = neighbour_data[i][3*j+1];
            int id_j = neighbour_data[i][3*j];

            for(int k=j+1; k<degree; k++)
            {
                int deg_k = neighbour_data[i][3*k+1];
                int id_k = neighbour_data[i][3*k];
                
                int smaller;
                int larger;

                if(deg_k>deg_j || ((deg_k == deg_j) && (id_k > id_j)))
                {
                    smaller = id_j;
                    larger = id_k;
                }
                else
                {
                    smaller = id_k;
                    larger = id_j;
                }

                int owner = smaller%world_size;

                queries[owner].push_back(id);
                queries[owner].push_back(smaller);
                queries[owner].push_back(larger);
            }
        }
    }

    int queries_sizes[world_size];

    for(int i=0; i<world_size; i++)
    {
        queries_sizes[i] = queries[i].size();
    }

    int query_size_recv[world_size];

    MPI_Alltoall(queries_sizes,1,MPI_INT, query_size_recv, 1, MPI_INT, MPI_COMM_WORLD);

    vector<vector<int>> queries_recv(world_size);

    for(int i=0; i<world_size; i++)
    {
        for(int j=0; j<query_size_recv[i]; j++)
        {
            if(j == my_rank)
            {
                continue;
            }
            queries_recv[i].push_back(0);
        }
    }

    for(int i=0; i<world_size; i++)
    {
        if(my_rank == i)
        {
            for(int j=0; j<world_size; j++)
            {
                if(j == i)
                {
                    continue;
                }
                MPI_Send(&(queries[j][0]), queries_sizes[j], MPI_INT, j, i, MPI_COMM_WORLD);
            }
        }
        else
        {
            MPI_Recv(&(queries_recv[i][0]), query_size_recv[i], MPI_INT, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

    }

    if(my_rank == 5)
    {
        int a = neighbour_data[0][12];
        cout<<a<<endl;
        
    }

    MPI_Finalize();
    return(0);
}