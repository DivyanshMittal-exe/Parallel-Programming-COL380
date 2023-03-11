#include<bits/stdc++.h>
#include<mpi.h>

using namespace std;

struct vertex
{
    int node_val;
    int node_deg;
    map<int,int> edge_triangles;
    vector<pair<int,int>> neighbour_degrees;
};

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    MPI_File input_data;
    MPI_File_open(MPI_COMM_WORLD, "test0/test-input-0.gra", MPI_MODE_RDONLY, MPI_INFO_NULL, &input_data);
    // MPI_File_open(MPI_COMM_WORLD, "/tmp/COL380/A2/test0/test-input-0.gra", MPI_MODE_RDONLY, MPI_INFO_NULL, &input_data);

    MPI_Aint lb = 0, extent = world_size*sizeof(int);
    MPI_Datatype etype = MPI_INT, filetype, contig;
    MPI_Offset disp;

    MPI_Type_contiguous(1, MPI_INT, &contig);

    MPI_Type_create_resized(contig, lb, extent, &filetype);

    MPI_Type_commit(&filetype);

    disp = my_rank*sizeof(int);

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

    vector<vertex> my_vertices;
    
    for(int i=0; i<my_num_nodes; i++)
    {
        
        struct vertex new_vertex;

        MPI_File_read_at(input_data, my_offsets[i], &(new_vertex.node_val), 1, MPI_INT , MPI_STATUS_IGNORE);
        MPI_File_read_at(input_data, my_offsets[i]+4, &(new_vertex.node_deg), 1, MPI_INT , MPI_STATUS_IGNORE);

        int neighbour_start = my_offsets[i]+8;

        for(int j=0; j<new_vertex.node_deg; j++)
        {
            int neighbour;
            MPI_File_read_at(input_data, neighbour_start+4*j, &neighbour, 1, MPI_INT , MPI_STATUS_IGNORE);

            int neighbour_offset;
            MPI_File_read_at(meta_data, neighbour*4, &neighbour_offset, 1, MPI_INT, MPI_STATUS_IGNORE);

            int neighbour_degree;
            MPI_File_read_at(input_data, neighbour_offset+4, &neighbour_degree, 1, MPI_INT, MPI_STATUS_IGNORE);

            if((neighbour_degree>new_vertex.node_deg) || (neighbour_degree == new_vertex.node_deg && neighbour>new_vertex.node_val))
            {
                new_vertex.edge_triangles.insert({neighbour,0});
                pair<int,int> new_pair;
                new_pair.first = neighbour;
                new_pair.second = neighbour_degree;
                new_vertex.neighbour_degrees.push_back(new_pair);
            }
        }

        my_vertices.push_back(new_vertex);
    }

    MPI_File_close(&meta_data);
    MPI_File_close(&input_data);

    vector<vector<pair<int,int>>> queries(world_size);

    for(int i=0; i<my_vertices.size(); i++)
    {
        struct vertex curr = my_vertices[i];
        
        if(curr.neighbour_degrees.size()>1)
        {
            for(int j=0; j<curr.neighbour_degrees.size(); j++)
            {
                for(int k=j+1; k<curr.neighbour_degrees.size(); k++)
                {
                    int smaller, larger;
                    if((curr.neighbour_degrees[j].second>curr.neighbour_degrees[k].second) || (curr.neighbour_degrees[j].second == curr.neighbour_degrees[k].second && curr.neighbour_degrees[j].first > curr.neighbour_degrees[k].second))
                    {
                        smaller = k;
                        larger = j;
                    }
                    else
                    {
                        smaller = j;
                        larger = k;
                    }

                    pair<int,int> temp;
                    temp.first = smaller;
                    temp.second = larger;
                    queries[smaller%world_size].push_back(temp);
                }
            }
        }
    }

    int query_lengths[world_size];

    for(int i=0; i<world_size; i++)
    {
        query_lengths[i] = queries[i].size();
    }

    int length_recv[world_size];

    MPI_Alltoall(query_lengths,1,MPI_INT, length_recv, 1, MPI_INT, MPI_COMM_WORLD);

    MPI_Finalize();

    return(0);
}