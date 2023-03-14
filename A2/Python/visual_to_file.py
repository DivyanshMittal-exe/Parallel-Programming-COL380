import struct

def create_binary_file(n, edges, filename):
    # open the binary file in write-binary mode
    with open(filename, "wb") as file:
        # calculate the number of edges
        m = len(edges)
        
        # write n and m to the file
        file.write(struct.pack("ii", n, m))
        
        # create a dictionary to store the neighbors of each node
        neighbors = {i: [] for i in range(n)}
        
        # loop over the edges and add them to the neighbors dictionary
        for u, v in edges:
            neighbors[u].append(v)
            neighbors[v].append(u)
        
        # loop over the nodes and write their information to the file
        for i in range(n):
            # get the node ID and its neighbors
            node_id = i
            node_neighbors = neighbors[i]
            node_deg = len(node_neighbors)
            
            # write the node ID and degree to the file
            file.write(struct.pack("ii", node_id, node_deg))
            
            # write the IDs of the node's neighbors to the file
            for neighbor_id in node_neighbors:
                file.write(struct.pack("i", neighbor_id))

# example usage
n = 4
edges = [(0, 1), (0, 10), (0, 14), (0, 9), (1, 9), (1, 2), (1, 3), (1, 18), (2, 8), (3, 6), (4, 6), (4, 7), (5, 12), (6, 8), (7, 8), (7, 19), (8, 16), (8, 9), (8, 10), (8, 13), (9, 19), (9, 15), (10, 18), (10, 16), (10, 17), (12, 19), (14, 19), (15, 17)]
create_binary_file(n, edges, "graph.bin")
