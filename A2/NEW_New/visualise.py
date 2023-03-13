import struct
import networkx as nx

# open the binary file in read-binary mode
with open("0.gra", "rb") as file:
    # read n and m from the first 8 bytes of the file
    n, m = struct.unpack("ii", file.read(8))

    # create an empty graph
    G = nx.Graph()

    # loop over the nodes and their neighbors
    for i in range(n):
        # read the node ID and degree from the file
        node_id, deg = struct.unpack("ii", file.read(8))

        # loop over the node's neighbors and add edges to the graph
        for j in range(deg):
            neighbor_id = struct.unpack("i", file.read(4))[0]
            G.add_edge(node_id, neighbor_id)

    # output all the edges of the graph
    print(len(G.edges()))

    for u, v in G.edges():
        print(f"{u} {v}")
