import os
import sys
sys.path.append(os.getcwd())
import dendropy
import numpy as np
import networkx as nx
from ete3 import Tree

nuc_names = ['A', 'C', 'G', 'T']
alphabet_size = len(nuc_names)

# DNA/RNA alphabet
nuc_vec = {'A': [1., 0., 0., 0.], 'C': [0., 1., 0., 0.], 'G': [0., 0., 1., 0.], 'T': [0., 0., 0., 1.], 'U': [0., 0., 0., 1.],
           '-': [0., 0., 0., 0.], '.': [0., 0., 0., 0.], '?': [0., 0., 0., 0.], 'N': [0.25, 0.25, 0.25, 0.25],
            'X': [0.25, 0.25, 0.25, 0.25], 'V': [0.33, 0.33, 0.33, 0.], 'H': [0.33, 0.33, 0., 0.33], 'D': [0.33, 0., 0.33, 0.33],
            'B': [0., 0.33, 0.33, 0.33], 'M': [0.5, 0.5, 0., 0.], 'R': [0.5, 0., 0.5, 0.], 'W': [0.5, 0., 0., 0.5],
            'S': [0., 0.5, 0.5, 0.], 'Y': [0., 0.5, 0., 0.5], 'K': [0., 0., 0.5, 0.5] }


def convert_data_str_to_onehot(data):
    """
    Converts given data's alphabet into one-hot vectors.

    Args:
        data (np.ndarray): Data array with shape (n_taxa, n_sites).

    Returns:
        data_onehot (np.ndarray): Data array with shape (n_taxa, n_sites, 4).
    """
    n_leaves = len(data)
    n_sites = len(data[0])
    data_onehot = np.ones((n_leaves, n_sites, alphabet_size))
    for i in range(n_leaves):
        data_onehot[i] = np.array([nuc_vec[c] for c in data[i]])
    return data_onehot


def newick2nx(newick_str, n_leaves, scale=0.1):
    """ Converts a newick string to Networkx graph. Newick -> Dendropy Tree -> NetworkX graph.
        TODO Beware! If the number of nodes in the newick string is greater than 2*taxa-2, it causes problem.
        It might create a cycle!!"""
    # Create Dendropy Tree object.
    dendro_tree = dendropy.Tree.get_from_string(newick_str, "newick")

    # Add taxa to internal nodes and convert leaf taxa to integers.
    root_visited = False
    n_nodes = 2 * n_leaves - 2
    node_idx = n_leaves
    for node in dendro_tree.preorder_node_iter():
        # Root
        if not root_visited:
            node.taxon = 2 * n_leaves - 3
            root_visited = True
        else:
            # Internal node
            if node.is_internal():
                node.taxon = node_idx
                node_idx = node_idx + 1
            # Leaf
            else:
                try:
                    node.taxon = int(str(node.taxon)[1:-1]) # Dendropy's leaf taxa has the form: 'name'. We take substring.
                except:
                    node.taxon = int(str(node.taxon)[2:-1])  # Sometimes the node names have "V" as well.

    # Convert Dendropy Tree to Networkx graph.
    tree = nx.Graph()
    node_names = [n for n in range(n_nodes)]
    leaves = node_names[:n_leaves]

    # Add nodes to the graph
    for node in leaves:
        tree.add_node(node, type='leaf')
    for node in node_names[n_leaves:-1]:
        tree.add_node(node, type='internal')
    tree.add_node(2 * n_leaves - 3, type='root')

    # Add edges to the graph
    for node in dendro_tree.preorder_node_iter():
        if node.is_internal():
            children = []
            for child_node in node.child_nodes():
                if child_node.edge_length is not None:
                    t = child_node.edge_length
                    #print('ULA')
                else:
                    t = np.random.exponential(scale)
                    #print('NNE')

                tree.add_edge(node.taxon, child_node.taxon, t=t)
                tree.add_node(child_node.taxon, parent=node.taxon)
                children.append(child_node.taxon)
            tree.add_node(node.taxon, children=children)
    return tree


def nx2ete(graph, root):
    """
    Converts a given .nx graph to ete3 tree.

    Args:
        graph (networkx.Graph): Input graph to be converted.
        root (int): root node
    Returns:
        tree (ete3.Tree): the resulting Tree object
    """
    tree = Tree()
    # Setting up a root node for lvl-1 to attach to
    tree.add_child(name=root)
    # A copy in a list, because you may not want to edit the original graph
    edges = list(graph.edges)
    while len(edges) > 0:
        for parent, child in edges:
            t = graph.edges[(parent, child)]['t']
            # check if this edge's parent is in the tree
            added_tree_list = list(tree.traverse("preorder"))
            if len(added_tree_list) > 0:
                for leaf in added_tree_list:
                    if leaf.name == parent:
                        # if it is, add child and thus create an edge
                        leaf.add_child(name=child, dist=t)

                        edges.remove((parent, child))
                    elif leaf.name == child:
                        # if it is, add child and thus create an edge
                        leaf.add_child(name=parent, dist=t)
                        edges.remove((parent, child))
        # Now if there are edges still unplaced, try again.
    return tree


def nx_to_newick(tree, root):
    """
    Converts a given .nx graph to newick string.
    nx -> ete -> nw
    Args:
        tree (networkx.Graph): Input graph to be converted.
        root (int): root node
    Returns:
        newick_str (string): the resulting newick string
    """
    newick_str = nx2ete(tree, root).write()
    # ETE adds extra root char, remove root
    newick_str = newick_str[1:len(newick_str)-5]+";"
    return newick_str







