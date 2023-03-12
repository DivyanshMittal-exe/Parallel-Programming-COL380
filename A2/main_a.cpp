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
    MPI_File_open(MPI_COMM_WORLD, "test1.gra", MPI_MODE_RDONLY, MPI_INFO_NULL, &input_data);

    MPI_File meta_data;
    MPI_File_open(MPI_COMM_WORLD, "test1.dat", MPI_MODE_RDONLY, MPI_INFO_NULL, &meta_data);

    MPI_File_set_view(meta_data, disp, etype, filetype, "native", MPI_INFO_NULL);

    int n, m;

    MPI_File_read_at_all(input_data, 0, &n, 1, MPI_INT, MPI_STATUS_IGNORE);
    MPI_File_read_at_all(input_data, 4, &m, 1, MPI_INT, MPI_STATUS_IGNORE);

    int my_num_nodes = (n/world_size) + (my_rank < n%world_size);

    int my_offsets[my_num_nodes];

    MPI_File_read_at_all(meta_data, 0, my_offsets, my_num_nodes, MPI_INT , MPI_STATUS_IGNORE);

    int my_vertex_ids[my_num_nodes];
    int my_vertex_deg[my_num_nodes];

    int prefc[my_num_nodes][n];

    for(int i=0; i<my_num_nodes; i++)
    {
        for(int j=0; j<n; j++)
        {
            prefc[my_num_nodes][n] = 0;
        }
    }

    // neighbour id, neighbour triangles, qr count, ktruss value, present
    vector<vector<int>> neighbour_data(my_num_nodes);

    // contains the coordinates in the neighbour data given a neighbour's value
    vector<map<int,int>> neighbour_coordinates(my_num_nodes);

    MPI_File_close(&meta_data);

    MPI_File_open(MPI_COMM_WORLD, "test1.dat", MPI_MODE_RDONLY, MPI_INFO_NULL, &meta_data);

    int all_degrees[n];

    for(int i=0; i<n; i++)
    {
        int offset;
        MPI_File_read_at(meta_data, i*4, &offset, 1, MPI_INT, MPI_STATUS_IGNORE);
        MPI_File_read_at(input_data, offset+4, &(all_degrees[i]), 1, MPI_INT , MPI_STATUS_IGNORE);
    }

    for(int i=0; i<my_num_nodes; i++)
    {
        MPI_File_read_at(input_data, my_offsets[i], &(my_vertex_ids[i]), 1, MPI_INT , MPI_STATUS_IGNORE);
        MPI_File_read_at(input_data, my_offsets[i]+4, &(my_vertex_deg[i]), 1, MPI_INT , MPI_STATUS_IGNORE);

        int neighbour_start = my_offsets[i]+8;

        for(int j=0; j<my_vertex_deg[i]; j++)
        {
            int neighbour;
            MPI_File_read_at(input_data, neighbour_start+4*j, &neighbour, 1, MPI_INT , MPI_STATUS_IGNORE);

            if(all_degrees[neighbour] > my_vertex_deg[i] || ((all_degrees[neighbour] == my_vertex_deg[i]) && (neighbour > my_vertex_ids[i])))
            {
                neighbour_data[i].push_back(neighbour);
                neighbour_data[i].push_back(0);
                neighbour_data[i].push_back(0);
                neighbour_data[i].push_back(0);
                neighbour_data[i].push_back(1);
                neighbour_coordinates[i].insert({neighbour, neighbour_data[i].size()-5});
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
        int degree = neighbour_data[i].size()/5;
        if(degree <2)
        {
            continue;
        }
        for(int j=0; j<degree; j++)
        {
            int id_j = neighbour_data[i][5*j];
            int deg_j = all_degrees[id_j];

            for(int k=j+1; k<degree; k++)
            {
                int id_k = neighbour_data[i][5*k];
                int deg_k = all_degrees[id_k];
                
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
        if(i == my_rank)
        {
            continue;
        }
        for(int j=0; j<query_size_recv[i]; j++)
        {
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

    // if(my_rank == 5)
    // {
    //     cout<<queries[4][queries[4].size()-1]<<" 5"<<endl;
    // }
    // if(my_rank == 4)
    // {
    //     cout<<queries_recv[5][queries_recv[5].size()-1]<<" 4"<<endl;
    // }

    vector<vector<int>> query_response(world_size);

    for(int i=0; i<world_size; i++)
    {
        if(i == my_rank)
        {
            continue;
        }
        for(int j=0; j<queries_recv[i].size()/3; j++)
        {
            int id = queries_recv[i][3*j];
            int node1 = queries_recv[i][3*j+1];
            int node2 = queries_recv[i][3*j+2];

            int node1_index = node1/world_size;

            if(neighbour_coordinates[node1_index].find(node2) != neighbour_coordinates[node1_index].end())
            {
                int coord = neighbour_coordinates[node1_index][node2];
                neighbour_data[node1_index][coord+1]++;
                neighbour_data[node1_index][coord+2]++;

                prefc[node1_index][id]++;

                query_response[i].push_back(id);
                query_response[i].push_back(node1);
                query_response[i].push_back(node2);
            }
        }
    }   

    int query_resp_lengths[world_size];

    for(int i=0; i<world_size; i++)
    {
        query_resp_lengths[i] = query_response[i].size();
    }

    int query_resp_lengths_recv[world_size];

    MPI_Alltoall(query_resp_lengths,1,MPI_INT, query_resp_lengths_recv, 1, MPI_INT, MPI_COMM_WORLD);

    // if(my_rank == 4)
    // {
    //     cout<<query_resp_lengths[3]<<" 4"<<endl;
    // }
    // if(my_rank == 3)
    // {
    //     cout<<query_resp_lengths_recv[4]<<" 3"<<endl;
    // }

    vector<vector<int>> query_resp_recv(world_size);
    
    for(int i=0; i<world_size; i++)
    {
        if(i == my_rank)
        {
            continue;
        }
        for(int j=0; j<query_resp_lengths_recv[i]; j++)
        {
            query_resp_recv[i].push_back(0);
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
                MPI_Send(&(query_response[j][0]), query_resp_lengths[j], MPI_INT, j, i, MPI_COMM_WORLD);
            }
        }
        else
        {
            MPI_Recv(&(query_resp_recv[i][0]), query_resp_lengths_recv[i], MPI_INT, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

    }

    // if(my_rank == 0)
    // {
    //     cout<<query_response[1][5]<<" 0"<<endl;
    // }
    // if(my_rank == 1)
    // {
    //     cout<<query_resp_recv[0][5]<<" 1"<<endl;
    // }

    for(int i=0; i<world_size; i++)
    {
        if(i == my_rank)
        {
            continue;
        }
        for(int j=0; j<query_resp_recv[i].size()/3; j++)
        {
            int id = query_resp_recv[i][3*j];
            int node1 = query_resp_recv[i][3*j+1];
            int node2 = query_resp_recv[i][3*j+2];

            int node_index = id/world_size;

            int coord1 = neighbour_coordinates[node_index][node1];
            int coord2 = neighbour_coordinates[node_index][node2];
            neighbour_data[node_index][coord1+1]++;
            neighbour_data[node_index][coord2+1]++;
        }
    }
    
    for(int i=0; i<queries[my_rank].size()/3; i++)
    {
        int id = queries[my_rank][3*i];
        int node1 = queries[my_rank][3*i+1];
        int node2 = queries[my_rank][3*i+2];

        int id_index = id/world_size;
        
        int node1_index = node1/world_size;

        int coord_id_node1 = neighbour_coordinates[id_index][node1];
        int coord_id_node2 = neighbour_coordinates[id_index][node2];

        if(neighbour_coordinates[node1_index].find(node2) != neighbour_coordinates[node1_index].end())
        {
            int coord_node1_node2 = neighbour_coordinates[node1_index][node2];
            neighbour_data[id_index][coord_id_node1+1]++;
            neighbour_data[id_index][coord_id_node2+1]++;
            neighbour_data[node1_index][coord_node1_node2+1]++;
            neighbour_data[node1_index][coord_node1_node2+2]++;

            prefc[node1_index][id]++;

        }
    }

    // if(my_rank == 3)
    // {
    //     for(int i=0; i<my_num_nodes; i++)
    //     {
    //         for(int j=0; j<neighbour_data[i].size(); j++)
    //         {
    //             cout<<neighbour_data[i][j]<<endl;
    //         }
    //         cout<<"----------------------"<<endl;
    //     }
    // }

    int K = 3; // truss value

    bool all_edge_deleted;

    bool edges_remain;

    // for(int i=0; i<my_num_nodes; i++)
    // {
    //     for(int j=0; j<neighbour_data[i].size()/5; j++)
    //     {
    //         if(neighbour_data[i][5*j+1] <= 0)
    //         {
    //             neighbour_data[i][5*j+4] = 0;
    //             neighbour_data[i][5*j+3] = 2;
    //         }
    //     }
    // }

    do
    {
        do
        {
            bool edge_deleted = 0;
            vector<vector<int>> query_unroll_wedge(world_size);
            for(int i=0; i<my_num_nodes; i++)
            {
                int id = my_vertex_ids[i];
                int degree = neighbour_data[i].size()/5;
                if(degree <2)
                {
                    continue;
                }
                for(int j=0; j<degree; j++)
                {
                    int id_j = neighbour_data[i][5*j];
                    int deg_j = all_degrees[id_j];
                    if(neighbour_data[i][5*j+4] == 0)
                    {
                        continue;
                    }

                    for(int k=j+1; k<degree; k++)
                    {
                        if(neighbour_data[i][5*k+4] == 0)
                        {
                            continue;
                        }
                        int id_k = neighbour_data[i][5*k];
                        int deg_k = all_degrees[id_k];

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

                        if(neighbour_data[i][5*j+1]<K-2 || neighbour_data[i][5*k+1]<K-2)
                        {
                            query_unroll_wedge[owner].push_back(id);
                            query_unroll_wedge[owner].push_back(smaller);
                            query_unroll_wedge[owner].push_back(larger);
                        }

                    }
                }
            }
            int query_unroll_wedge_sizes[world_size];

            for(int i=0; i<world_size; i++)
            {
                query_unroll_wedge_sizes[i] = query_unroll_wedge[i].size();
            }

            int unroll_wedge_size_recv[world_size];

            MPI_Alltoall(query_unroll_wedge_sizes,1,MPI_INT, unroll_wedge_size_recv, 1, MPI_INT, MPI_COMM_WORLD);

            vector<vector<int>> unroll_recv(world_size);

            for(int i=0; i<world_size; i++)
            {
                if(i == my_rank)
                {
                    continue;
                }
                for(int j=0; j<unroll_wedge_size_recv[i]; j++)
                {
                    unroll_recv[i].push_back(0);
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
                        
                        MPI_Send(&(query_unroll_wedge[j][0]), query_unroll_wedge_sizes[j], MPI_INT, j, i, MPI_COMM_WORLD);
                    }
                }
                else
                {
                    MPI_Recv(&(unroll_recv[i][0]), unroll_wedge_size_recv[i], MPI_INT, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }

            }

            vector<vector<int>> wedge_response(world_size);

            for(int i=0; i<world_size; i++)
            {
                if(i == my_rank)
                {
                    continue;
                }
                for(int j=0; j<unroll_recv[i].size()/3; j++)
                {
                    int id = unroll_recv[i][3*j];
                    int node1 = unroll_recv[i][3*j+1];
                    int node2 = unroll_recv[i][3*j+2];

                    int node1_index = node1/world_size;

                    if(neighbour_coordinates[node1_index].find(node2) != neighbour_coordinates[node1_index].end())
                    {
                        int coord = neighbour_coordinates[node1_index][node2];
                        if(neighbour_data[node1_index][coord+4] == 1)
                        {
                            neighbour_data[node1_index][coord+1]--;
                            if(neighbour_data[node1_index][coord+1] <= 0)
                            {
                                edge_deleted = 1;
                                if(neighbour_data[node1_index][coord+1] == 0)
                                    neighbour_data[node1_index][coord+3] = K-1;
                                neighbour_data[node1_index][coord+4] = 0;
                            }
                            neighbour_data[node1_index][coord+2]--;

                            prefc[node1_index][id]--;

                            wedge_response[i].push_back(id);
                            wedge_response[i].push_back(node1);
                            wedge_response[i].push_back(node2);
                        }
                    }
                }
            }

            int wr_len[world_size];
            for(int i=0; i<world_size; i++)
            {
                wr_len[i] = wedge_response[i].size();
            }

            int wr_len_recv[world_size];

            MPI_Alltoall(wr_len,1,MPI_INT, wr_len_recv, 1, MPI_INT, MPI_COMM_WORLD);

            vector<vector<int>> wr_recv(world_size);
    
            for(int i=0; i<world_size; i++)
            {
                if(i == my_rank)
                {
                    continue;
                }
                for(int j=0; j<wr_len_recv[i]; j++)
                {
                    wr_recv[i].push_back(0);
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
                        MPI_Send(&(wedge_response[j][0]), wr_len[j], MPI_INT, j, i, MPI_COMM_WORLD);
                    }
                }
                else
                {
                    MPI_Recv(&(wr_recv[i][0]), wr_len_recv[i], MPI_INT, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }

            }

            for(int i=0; i<world_size; i++)
            {
                if(i == my_rank)
                {
                    continue;
                }
                for(int j=0; j<wr_recv[i].size()/3; j++)
                {
                    int id = wr_recv[i][3*j];
                    int node1 = wr_recv[i][3*j+1];
                    int node2 = wr_recv[i][3*j+2];

                    int node_index = id/world_size;

                    int coord1 = neighbour_coordinates[node_index][node1];
                    int coord2 = neighbour_coordinates[node_index][node2];
                    neighbour_data[node_index][coord1+1]--;
                    neighbour_data[node_index][coord2+1]--;
                    if(neighbour_data[node_index][coord1+1] <= 0)
                    {
                        edge_deleted = 1;
                        if(neighbour_data[node_index][coord1+1] == 0)
                            neighbour_data[node_index][coord1+3] = K-1;
                        neighbour_data[node_index][coord1+4] = 0;
                    }
                    if(neighbour_data[node_index][coord1+1] <= 0)
                    {
                        edge_deleted = 1;
                        if(neighbour_data[node_index][coord1+1] == 0)
                            neighbour_data[node_index][coord1+3] = K-1;
                        neighbour_data[node_index][coord1+4] = 0;
                    }
                }
            }

            for(int i=0; i<query_unroll_wedge[my_rank].size()/3; i++)
            {
                int id = query_unroll_wedge[my_rank][3*i];
                int node1 = query_unroll_wedge[my_rank][3*i+1];
                int node2 = query_unroll_wedge[my_rank][3*i+2];

                int id_index = id/world_size;

                int node1_index = node1/world_size;

                int coord_id_node1 = neighbour_coordinates[id_index][node1];
                int coord_id_node2 = neighbour_coordinates[id_index][node2];

                if(neighbour_coordinates[node1_index].find(node2) != neighbour_coordinates[node1_index].end())
                {
                    int coord_node1_node2 = neighbour_coordinates[node1_index][node2];

                    if(neighbour_data[node1_index][coord_node1_node2+4] == 1)
                    {
                        neighbour_data[id_index][coord_id_node1+1]--;
                        neighbour_data[id_index][coord_id_node2+1]--;
                        neighbour_data[node1_index][coord_node1_node2+1]--;
                        if(neighbour_data[id_index][coord_id_node1+1] <= 0)
                        {
                            edge_deleted = 1;
                            if(neighbour_data[id_index][coord_id_node1+1] == 0)
                                neighbour_data[id_index][coord_id_node1+3] = K-1;
                            neighbour_data[id_index][coord_id_node1+4] = 0;
                        }
                        if(neighbour_data[id_index][coord_id_node2+1] <= 0)
                        {
                            edge_deleted = 1;
                            if(neighbour_data[id_index][coord_id_node2+1] == 0)
                                neighbour_data[id_index][coord_id_node2+3] = K-1;
                            neighbour_data[id_index][coord_id_node2+4] = 0;
                        }
                        if(neighbour_data[node1_index][coord_node1_node2+1] <= 0)
                        {
                            edge_deleted = 1;
                            if(neighbour_data[node1_index][coord_node1_node2+1] == 0)
                                neighbour_data[node1_index][coord_node1_node2+3] = K-1;
                            neighbour_data[node1_index][coord_node1_node2+4] = 0;
                        }
                        neighbour_data[node1_index][coord_node1_node2+2]--;

                        prefc[node1_index][id]--;
                    }

                }
            }

            //At this point unroll wedge has been completed

            vector<vector<int>> qr_send(world_size);
            for(int i=0; i<my_num_nodes; i++)
            {
                int id = my_vertex_ids[i];
                int degree = neighbour_data[i].size()/5;

                for(int j=0; j<degree; j++)
                {
                    if(neighbour_data[i][5*j+4] == 0)
                    {
                        continue;
                    }
                    int id_j = neighbour_data[i][5*j];

                    if(neighbour_data[i][5*j+1]<K-2 && neighbour_data[i][5*j+2]>0)
                    {
                        for(int k=0; k<n; k++)
                        {
                            if(neighbour_data[i][5*k+4] == 1)
                            {
                                continue;
                            }
                            if(prefc[i][k]>0)
                            {
                                int owner = k%world_size;

                                qr_send[owner].push_back(k);
                                qr_send[owner].push_back(id);
                                qr_send[owner].push_back(id_j);
                            }
                        }
                    }
                }
            }

            int qr_send_sizes[world_size];

            for(int i=0; i<world_size; i++)
            {
                qr_send_sizes[i] = qr_send[i].size();
            }

            int qr_size_recv[world_size];

            MPI_Alltoall(qr_send_sizes,1,MPI_INT, qr_size_recv, 1, MPI_INT, MPI_COMM_WORLD);

            vector<vector<int>> qr_recv(world_size);

            for(int i=0; i<world_size; i++)
            {
                if(i == my_rank)
                {
                    continue;
                }
                for(int j=0; j<qr_size_recv[i]; j++)
                {
                    qr_recv[i].push_back(0);
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
                        MPI_Send(&(qr_send[j][0]), qr_send_sizes[j], MPI_INT, j, i, MPI_COMM_WORLD);
                    }
                }
                else
                {
                    MPI_Recv(&(qr_recv[i][0]), qr_size_recv[i], MPI_INT, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }

            }

            vector<vector<int>> qr_resp(world_size);

            for(int i=0; i<world_size; i++)
            {
                if(i == my_rank)
                {
                    continue;
                }
                for(int j=0; j<qr_recv[i].size()/3; j++)
                {
                    int id = qr_recv[i][3*j];
                    int node1 = qr_recv[i][3*j+1];
                    int node2 = qr_recv[i][3*j+2];

                    int id_index = id/world_size;

                    if(neighbour_coordinates[id_index].find(node2) != neighbour_coordinates[id_index].end())
                    {
                        int coord1 = neighbour_coordinates[id_index][node1];
                        int coord2 = neighbour_coordinates[id_index][node2];
                        if(neighbour_data[id_index][coord2+4] == 1)
                        {
                            neighbour_data[id_index][coord1+1]--;
                            neighbour_data[id_index][coord2+1]--;

                            if(neighbour_data[id_index][coord1+1] <= 0)
                            {
                                edge_deleted = 1;
                                if(neighbour_data[id_index][coord1+1] == 0)
                                    neighbour_data[id_index][coord1+3] = K-1;
                                neighbour_data[id_index][coord1+4] = 0;
                            }

                            if(neighbour_data[id_index][coord2+1] <= 0)
                            {
                                edge_deleted = 1;
                                if(neighbour_data[id_index][coord2+1] == 0)
                                    neighbour_data[id_index][coord2+3] = K-1;
                                neighbour_data[id_index][coord2+4] = 0;
                            }

                            qr_resp[i].push_back(id);
                            qr_resp[i].push_back(node1);
                            qr_resp[i].push_back(node2);
                        }
                    }
                }
            }

            int qr_resp_len[world_size];
            for(int i=0; i<world_size; i++)
            {
                qr_resp_len[i] = qr_resp[i].size();
            }
            
            int qr_resp_len_recv[world_size];

            MPI_Alltoall(qr_resp_len,1,MPI_INT, qr_resp_len_recv, 1, MPI_INT, MPI_COMM_WORLD);

            vector<vector<int>> qr_resp_recv(world_size);

            for(int i=0; i<world_size; i++)
            {
                if(i == my_rank)
                {
                    continue;
                }
                for(int j=0; j<qr_resp_len_recv[i]; j++)
                {
                    qr_resp_recv[i].push_back(0);
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
                        MPI_Send(&(qr_resp[j][0]), qr_resp_len[j], MPI_INT, j, i, MPI_COMM_WORLD);
                    }
                }
                else
                {
                    MPI_Recv(&(qr_resp_recv[i][0]), qr_resp_len_recv[i], MPI_INT, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }

            }

            for(int i=0; i<world_size; i++)
            {
                if(i == my_rank)
                {
                    continue;
                }
                for(int j=0; j<qr_resp_recv[i].size()/3; j++)
                {
                    int id = qr_resp_recv[i][3*j];
                    int node1 = qr_resp_recv[i][3*j+1];
                    int node2 = qr_resp_recv[i][3*j+2];

                    int node_index = node1/world_size;

                    int coord = neighbour_coordinates[node_index][node2];

                    neighbour_data[node_index][coord+1]--;
                    neighbour_data[node_index][coord+2]--;

                    if(neighbour_data[node_index][coord+1] <= 0)
                    {
                        edge_deleted = 1;
                        if(neighbour_data[node_index][coord+1] == 0)
                            neighbour_data[node_index][coord+3] = K-1;
                        neighbour_data[node_index][coord+4] = 0;
                    }

                }
            }

            for(int i=0; i<qr_send[my_rank].size()/3; i++)
            {
                int id = qr_send[my_rank][3*i];
                int node1 = qr_send[my_rank][3*i+1];
                int node2 = qr_send[my_rank][3*i+2];

                int id_index = id/world_size;

                int node1_index = node1/world_size;

                int coord_id_node1 = neighbour_coordinates[id_index][node1];
                int coord_node1_node2 = neighbour_coordinates[node1_index][node2];

                if(neighbour_coordinates[id_index].find(node2) != neighbour_coordinates[id_index].end())
                {
                    int coord_id_node2 = neighbour_coordinates[id_index][node2];

                    if(neighbour_data[id_index][coord_id_node2+4] == 1)
                    {

                        neighbour_data[id_index][coord_id_node1+1]--;
                        neighbour_data[id_index][coord_id_node2+1]--;
                        neighbour_data[node1_index][coord_node1_node2+1]--;
                        if(neighbour_data[id_index][coord_id_node1+1] <= 0)
                        {
                            edge_deleted = 1;
                            if(neighbour_data[id_index][coord_id_node1+1] == 0)
                                neighbour_data[id_index][coord_id_node1+3] = K-1;
                            neighbour_data[id_index][coord_id_node1+4] = 0;
                        }
                        if(neighbour_data[id_index][coord_id_node1+1] <= 0)
                        {
                            edge_deleted = 1;
                            if(neighbour_data[id_index][coord_id_node1+1] == 0)
                                neighbour_data[id_index][coord_id_node1+3] = K-1;
                            neighbour_data[id_index][coord_id_node1+4] = 0;
                        }
                        if(neighbour_data[node1_index][coord_node1_node2+1] <= 0)
                        {
                            edge_deleted = 1;
                            if(neighbour_data[node1_index][coord_node1_node2+1] == 0)
                                neighbour_data[node1_index][coord_node1_node2+3] = K-1;
                            neighbour_data[node1_index][coord_node1_node2+4] = 0;
                        }

                        neighbour_data[node1_index][coord_node1_node2+2]--;
                    }

                }
            }

            // At this point qr unroll has been completed
            MPI_Allreduce(&edge_deleted, &all_edge_deleted, 1, MPI_CXX_BOOL, MPI_LAND, MPI_COMM_WORLD);

        }
        while(all_edge_deleted);

        K++;

        bool my_remain = 0;
        for(int i=0; i<my_num_nodes; i++)
        {
            bool found;
            for(int j=0; j<neighbour_data[i].size()/5; j++)
            {
                if(neighbour_data[i][5*j+4] == 1 && neighbour_data[i][5*j+1>0])
                {
                    found = 1;
                    // cout<<my_rank<<" found node "<<i*world_size+my_rank<<" with neighbour "<<neighbour_data[i][5*j]<<" "<<neighbour_data[i][5*j+1]<<endl;
                    my_remain = 1;
                    break;
                }
            }
            if(found)
            {
                break;
            }
        }
        MPI_Allreduce(&my_remain, &edges_remain, 1, MPI_CXX_BOOL, MPI_LOR, MPI_COMM_WORLD);
        // cout<<edges_remain<<" edges_remain"<<endl;
    }
    while(edges_remain);

    if(my_rank ==4)
    {
        for(int i=0; i<my_num_nodes; i++)
        {
            for(int j=0; j<neighbour_data[i].size(); j++)
            {
                cout<<neighbour_data[i][j]<<endl;
            }
            cout<<"------------------"<<endl;
        }
    }

    MPI_Finalize();
    return(0);
}