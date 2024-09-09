import os
import sys
import pickle

sys.path.append(os.getcwd())
import time
import numpy as np
import networkx as nx
from collections import defaultdict
from skbio import DistanceMatrix
from skbio.tree import nj
import copy
import argparse
from snj.utils import convert_data_str_to_onehot, newick2nx, nx_to_newick
import newick
from Bio import Align

print('pairwise version for unaligned sequences...')

def DFS(src, visited, subtree_size, subtree_leaves, parent_array, n, tree):
    """
    Perfoms depth-first search (DFS) on a given tree.

    Args:
        src(int): root node
        visited (np.ndarray): boolean array keeping track of nodes already visited
        subtree_size (np.ndarray): keeps the sizes of subtrees
        subtree_leaves (np.ndarray): keeps the leaves at subtrees
        parent_array (np.ndarray): keeps the parents of nodes
        n(int): the number of nodes
        tree (networkx.Graph): undirected graph representing the phylogenetic tree
    """
    visited[src] = True
    n[0] += 1
    subtree_size[src] = 1
    if src < n_leaves:
        subtree_leaves[src].append(src)
    for adj in tree.adj[src]:
        if not visited[adj] and not centroidMarked[adj]:
            DFS(adj, visited, subtree_size, subtree_leaves, parent_array, n, tree)
            subtree_size[src] += subtree_size[adj]
            parent_array[adj] = src
            subtree_leaves[src] += subtree_leaves[adj]


def getCentroid(src, visited, subtree_size, n, tree):
    """
    Finds a centroid for a given tree

    Args:
        src(int): root node
        visited (np.ndarray): boolean array keeping track of nodes already visited
        subtree_size (np.ndarray): keeps the sizes of subtrees
        n(int): the number of nodes
        tree (networkx.Graph): undirected graph representing the phylogenetic tree
    """
    is_centroid = True
    visited[src] = True
    heaviest_child = 0

    for adj in tree.adj[src]:
        if not visited[adj] and not centroidMarked[adj]:
            if subtree_size[adj] > n / 2:
                is_centroid = False

            if heaviest_child == 0 or subtree_size[adj] > subtree_size[heaviest_child]:
                heaviest_child = adj

    if is_centroid and n - subtree_size[src] <= n / 2:
        return src

    return getCentroid(heaviest_child, visited, subtree_size, n, tree)


# function to first do DFS (if necessary) and then find a centroid
def getCentroidTree(src, tree, subtree_size=None):
    """
    First does DFS (if necessary) and then finds a centroid

    Args:
        src(int): root
        tree (networkx.Graph): undirected graph representing the phylogenetic tree
        subtree_size (np.ndarray): keeps the sizes of subtrees

    Returns:
        centroid(int): centroid node for the given tree.
    """
    if subtree_size == None:
        visited = [False] * MAXN
        subtree_size = [0] * MAXN
        parent_array = [-1] * MAXN
        subtree_leaves = defaultdict(list)
        n = [0]

        DFS(src, visited, subtree_size, subtree_leaves, parent_array, n, tree)
    else:
        n = [subtree_size[src]]

    visited = [False] * MAXN
    centroid = getCentroid(src, visited, subtree_size, n[0], tree)
    centroidMarked[centroid] = True

    return centroid


def orient_pick(held_out, leaves_at_adj, num_of_ort):
    """
    Picks orienting leaves from a given subtree

    Args:
        held_out(int): the leaf to be placed
        leaves_at_adj (np.ndarray): leaves at a given subtree
        num_of_ort (int): the number of orienting leaves

    Returns:
        ort_leaves(list): orienting leaves.
    """
    num_leaf_at_adj = len(leaves_at_adj)
    # if less leaves than the number of orienting leaves, no need any heuristic
    if num_leaf_at_adj < num_of_ort:
        num_leaf_to_sample = num_of_ort - num_leaf_at_adj
        ort_leaves = list(leaves_at_adj) + list(np.random.choice(leaves_at_adj, num_leaf_to_sample))
    # if the number of leaves equals to the number of orienting leaves, no need any heuristic
    elif num_leaf_at_adj == num_of_ort:
        ort_leaves = list(leaves_at_adj.copy())
    # if more leaves than the number of orienting leaves, use heuristic
    else:
        # check if the number of leaves is greater than the number of sampling leaves
        if num_leaf_at_adj > upp_limit:
            leaves_at_adj_samp = np.random.choice(leaves_at_adj, upp_limit)
        else:
            leaves_at_adj_samp = leaves_at_adj.copy()

        # select the closest leaves as the orienting leaves
        dist_vector = np.zeros([len(leaves_at_adj_samp)])
        for or_i in range(len(leaves_at_adj_samp)):
            # align
            aligner = Align.PairwiseAligner(mode='global', match_score=2, mismatch_score=-1)
            aligner.open_gap_score = -5
            aligner.extend_gap_score = -1
            # aligner.target_end_gap_score = 0.0
            # aligner.query_end_gap_score = 0.0
            seq1 = "".join([str(i) for i in data[held_out, :]])
            seq2 = "".join([str(i) for i in data[leaves_at_adj_samp[or_i], :]])
            alignments = aligner.align(seq1, seq2)
            # choose best alignment
            alignment = sorted(alignments)[0]
            aligned_seqs = np.array(alignment, dtype=str)
            # compute the distance
            aligned_onehot = convert_data_str_to_onehot(aligned_seqs)
            dist_vector[or_i] = np.abs(aligned_onehot[0, :, :] - aligned_onehot[1, :, :]).sum()
        #data_leaves_onehot = data_onehot[leaves_at_adj_samp, :, :]
        #data_heldout_onehot = np.squeeze(data_onehot[held_out, :, :])
        #dist_vector = np.abs(data_leaves_onehot[:, :, :] - data_heldout_onehot[:, :]).sum(axis=(-1, -2))

        for i in range(len(leaves_at_adj_samp)):
            computed_distances[held_out][leaves_at_adj_samp[i]] = dist_vector[i]
        dist_vector = dist_vector / (2 * seq_length)
        dist_vector = (-3 / 4) * np.log(1 - dist_vector * 4 / 3)
        ort_leaves = leaves_at_adj_samp[np.argsort(dist_vector)[0: num_of_ort]]
        ort_leaves = list(ort_leaves)
    return ort_leaves


def decomposeTree(root, tree, held_out, current_degree=0, first_placement=False, subtree_leaves=None,
                  subtree_size=None, parent_array=None, root_parent=-1):
    """
    Quartet-based subtree selection, stops when the subtree reduces to a single edge.

    Args:
        root(int): root node of the given subtree
        tree (networkx.Graph): the given subtree
        held_out(int): the leaf to be placed
        current_degree(int): the number of recursive calls
        first_placement(boolean): if it is the first placement
        subtree_size (np.ndarray): keeps the sizes of subtrees
        subtree_leaves (np.ndarray): keeps the leaves at subtrees
        parent_array (np.ndarray): keeps the parents of nodes
        root_parent: the parent node of the root node of the subtree (if the subtree is the whole tree, then -1)

    Returns:
        selected_adj (int): one of the two nodes where the selected edge lies between
        next_selected_adj (int): the other node
        current_degree (int): the number of recursive calls
    """
    # if it is first recursion and centroid is not known, find the centroid
    if (
            root, current_degree, tuple(stop_node_dict[root]),
            root_parent) not in root_centroid_dict and current_degree == 0:
        cend_tree = getCentroidTree(root, tree)
    # if it is not first recursion and centroid is not known, find the centroid
    elif (
            root, current_degree, tuple(stop_node_dict[root]),
            root_parent) not in root_centroid_dict and current_degree != 0:
        cend_tree = getCentroidTree(root, tree, subtree_size)
        root_centroid_dict[(root, current_degree, tuple(stop_node_dict[root]), root_parent)] = cend_tree
    # if the centroid is known, check if it is still up
    elif (root, current_degree, tuple(stop_node_dict[root]), root_parent) in root_centroid_dict:
        cend_tree = root_centroid_dict[(root, current_degree, tuple(stop_node_dict[root]), root_parent)]
        # check if the centroid is still up, it may need to shift one node due to last placement
        if current_degree != 0:
            curr_adj_sizes = []
            curr_adj_list = []
            for adj in tree.adj[cend_tree]:
                curr_adj_sizes.append(subtree_size[adj])
                curr_adj_list.append(adj)

            sorted_adj_sizes = np.sort(np.array(curr_adj_sizes.copy()))
            sorted_adj_list = np.array(curr_adj_list.copy())
            sorted_adj_list = sorted_adj_list[np.argsort(np.array(curr_adj_sizes))]
            sorted_adj_sizes[-1] = subtree_size[root] - sorted_adj_sizes[0] - sorted_adj_sizes[1] - 1

            # fix the centroid if necessary
            if (sorted_adj_sizes > subtree_size[root] / 2).any():
                cend_tree = sorted_adj_list[np.argmax(sorted_adj_sizes)]
                root_centroid_dict[(root, current_degree, tuple(stop_node_dict[root]), root_parent)] = cend_tree

        centroidMarked[cend_tree] = True

    # if it is the first placement, save the DFS results for future use
    if first_placement and current_degree == 0:
        visited = [False] * MAXN
        subtree_size = [0] * MAXN
        n = [0]
        parent_array = [-1] * MAXN
        subtree_leaves = defaultdict(list)

        DFS(cend_tree, visited, subtree_size, subtree_leaves, parent_array, n, tree)
        root_centroid_dict[(root, current_degree, tuple(stop_node_dict[root]), root_parent)] = cend_tree
        global first_dfs_subtree_leaves
        global first_dfs_subtree_sizes
        global first_parent_array
        first_dfs_subtree_leaves = copy.deepcopy(subtree_leaves)
        first_dfs_subtree_sizes = copy.deepcopy(subtree_size)
        first_parent_array = copy.deepcopy(parent_array)

    # initialize orienting leaves and subtree sizes
    orienting_leaves = []
    subtree_leaf_size = []
    active_adj = []
    if current_degree == 0:
        subtree_all_leaves = subtree_leaves[cend_tree]
    else:
        subtree_all_leaves = subtree_leaves[root]

    # detect the subtrees for the centroid of interest, and get orienting leaves
    # do it first for the two subtrees that are children of the centroid
    for adj in tree.adj[cend_tree]:
        if adj != parent_array[cend_tree]:
            leaves_at_adj = np.array(subtree_leaves[adj])
            subtree_leaf_size.append(leaves_at_adj.shape[0])
            subtree_all_leaves = list([leaf for leaf in subtree_all_leaves if leaf not in leaves_at_adj])
            if leaves_at_adj.shape[0] == 0:  # if no leaves, then consider the partition of the whole tree by the 'cend_tree'
                leaves_at_adj = np.array(first_dfs_subtree_leaves[adj].copy())
            orienting_leaves.append(orient_pick(held_out, leaves_at_adj, num_of_ort))
            active_adj.append(adj)
        else:
            parent_adj = adj
    # now for the other subtree
    if len(subtree_all_leaves) > 0:  # if it has leaves
        subtree_leaf_size.append(len(subtree_all_leaves))
        orienting_leaves.append(orient_pick(held_out, np.array(subtree_all_leaves), num_of_ort))
        active_adj.append(parent_adj)
    elif current_degree != 0:  # else consider the partition of the whole tree
        subtree_leaf_size.append(len(subtree_all_leaves))
        grand_parent = parent_array[root]
        leaves_at_grand_parent = first_dfs_subtree_leaves[grand_parent].copy()
        leaves_at_sibling = list(
            [leaf for leaf in leaves_at_grand_parent if leaf not in first_dfs_subtree_leaves[root]])
        orienting_leaves.append(orient_pick(held_out, np.array(leaves_at_sibling), num_of_ort))
        active_adj.append(parent_adj)

    # take orienting leaves in array format
    orients = np.zeros([num_of_ort, 4], dtype=int)
    orients[:, :-1] = np.array(orienting_leaves, dtype=int).T
    orients[:, -1] = held_out

    # compute distance matrix for the quartet
    # vectorized form works slower...
    dist_matrix = np.zeros([num_of_ort, 4, 4])
    #data_orients = data_onehot[orients]
    for k in range(num_of_ort):
        for l in range(1, 4):
            for m in range(l):
                if int(orients[k, l]) in computed_distances[int(orients[k, m])]:
                    dist_matrix[k, l, m] = computed_distances[int(orients[k, m])][int(orients[k, l])]
                elif int(orients[k, m]) in computed_distances[int(orients[k, l])]:
                    dist_matrix[k, l, m] = computed_distances[int(orients[k, l])][int(orients[k, m])]
                else:
                    # align
                    aligner = Align.PairwiseAligner(mode='global', match_score=2, mismatch_score=-1)

                    aligner.open_gap_score = -5
                    aligner.extend_gap_score = -1
                    # aligner.target_end_gap_score = 0.0
                    # aligner.query_end_gap_score = 0.0
                    seq1 = "".join([str(i) for i in data[int(orients[k, l]), :]])
                    seq2 = "".join([str(i) for i in data[int(orients[k, m]), :]])

                    alignments = aligner.align(seq1, seq2)
                    # choose best alignment
                    alignment = sorted(alignments)[0]
                    aligned_seqs = np.array(alignment, dtype=str)
                    # compute the distance
                    aligned_onehot = convert_data_str_to_onehot(aligned_seqs)
                    dist_matrix[k, l, m] = np.abs(aligned_onehot[0, :, :] - aligned_onehot[1, :, :]).sum()

                    #dist_matrix[k, l, m] = np.abs(data_orients[k, l, :, :] - data_orients[k, m, :, :]).sum()

                    computed_distances[int(orients[k, m])][int(orients[k, l])] = dist_matrix[k, l, m]


    dist_matrix = dist_matrix / (2 * seq_length)
    dist_matrix = (-3 / 4) * np.log(1 - (dist_matrix * 4 / 3))
    dist_matrix = dist_matrix.mean(axis=0)
    dist_matrix += dist_matrix.T

    # build a quartet using NJ
    dm = DistanceMatrix(dist_matrix)
    NJ_tree = nj(dm)

    held_out_idx = 3
    selected_path = int(NJ_tree.find(str(held_out_idx)).siblings()[0].name)
    selected_adj = active_adj[selected_path]

    # if the selected subtree is not a child of the centroid, update the data structures accordingly
    if not centroidMarked[selected_adj]:
        if selected_adj == parent_array[cend_tree]:
            grand_parent = parent_array[cend_tree]
            while grand_parent != -1 and grand_parent != parent_array[root]:
                subtree_leaves[grand_parent] = list(
                    [leaf for leaf in subtree_leaves[grand_parent] if leaf not in subtree_leaves[cend_tree]])
                subtree_size[grand_parent] -= subtree_size[cend_tree]
                stop_node_dict[grand_parent].append(cend_tree)
                grand_parent = parent_array[grand_parent]

            subtree_leaves[cend_tree] = []
            subtree_size[cend_tree] = 0
            selected_adj = root

    if selected_adj != root:
        root_parent = cend_tree
    if current_degree == 0:
        root_parent = cend_tree

    current_degree += 1

    # if the root of the selected subtree is an internal node and not an centroid, continue recursion
    # otherwise, return the two nodes of the selected edge
    if selected_adj >= n_leaves and not centroidMarked[selected_adj]:
        return decomposeTree(selected_adj, tree, held_out, current_degree=current_degree, subtree_leaves=subtree_leaves,
                             subtree_size=subtree_size, parent_array=parent_array, root_parent=root_parent,
                                )
    else:
        if selected_adj < n_leaves:
            next_selected_adj = parent_array[selected_adj]
        else:
            next_selected_adj = cend_tree

        return selected_adj, next_selected_adj, current_degree


##### MAIN ######
def main():
    # get hyperparameters
    parser = argparse.ArgumentParser()
    parser.add_argument('-data', required=True, help=' ds1 | ds2 | ... | ds8 | virus ')
    parser.add_argument('-seed', required=True, help=' 5 | 10 | 20  ')
    parser.add_argument('-n_i', required=True, help=' 50, 150, 250, 350 ')
    parser.add_argument('-n_s', required=True, help=' 5, 10, 15, 20, 25 ')
    parser.add_argument('-n_o', required=True, help=' 1, 3, 5, 7, 9')
    args = parser.parse_args()

    # initialize the model
    global n_leaves, upp_limit, num_of_ort, data, data_onehot, computed_distances, first_dfs_subtree_sizes, \
        first_dfs_subtree_leaves, first_parent_array, root_centroid_dict, centroidMarked, stop_node_dict, MAXN, seq_length

    # get data
    data = np.load('data/' + args.data + '.npy', allow_pickle=False)
    print('Data obtained succesfully.')

    n_leaves = data.shape[0]
    initial_n_leaves = int(args.n_i)
    upp_limit = int(args.n_s)
    num_of_ort = int(args.n_o)
    num_of_heldout = n_leaves - initial_n_leaves
    num_of_seed = int(args.seed)
    time_array = np.zeros([num_of_seed])
    seq_length = data.shape[1]
    #data_onehot = convert_data_str_to_onehot(data)

    computed_distances = [{} for i in range(n_leaves)]

    # loop for different seeds
    for seed in range(num_of_seed):
        np.random.seed(seed)
        MAXN = 2 * n_leaves - 2
        # some data structures needs to be reset for each seed
        first_dfs_subtree_sizes = [0] * MAXN  # []#
        first_dfs_subtree_leaves = defaultdict(list)  # []#
        first_parent_array = [-1] * MAXN
        root_centroid_dict = {}

        # RANDOM OPTION - to select the leaves for the initial backbone tree construction
        perm_idx = np.random.permutation(n_leaves)
        initial_leaves = perm_idx[0:initial_n_leaves]
        held_out_leaves = perm_idx[initial_n_leaves:]

        our_tic = time.process_time()
        # construct the initial backbone tree:
        #data_onehot_initial = data_onehot[initial_leaves, :, :]
        data_initial = data[initial_leaves, :]
        dist_matrix = np.zeros([initial_n_leaves, initial_n_leaves])  # np.zeros([4, 4])

        for ki in range(initial_n_leaves):
            for kj in range(ki+1, initial_n_leaves):
                # align

                aligner = Align.PairwiseAligner(mode='global', match_score=2, mismatch_score=-1)
                aligner.open_gap_score = -5
                aligner.extend_gap_score = -1
                #aligner.target_end_gap_score = 0.0
                #aligner.query_end_gap_score = 0.0
                seq1 = "".join([str(i) for i in data_initial[ki, :]])
                seq2 = "".join([str(i) for i in data_initial[kj, :]])
                alg_tic = time.process_time()
                alignments = aligner.align(seq1, seq2)
                alg_toc = time.process_time()
                #print('alingment time: ',alg_toc-alg_tic)
                # choose best alignment
                alignment = sorted(alignments)[0]
                aligned_seqs = np.array(alignment, dtype=str)
                # compute the distance
                aligned_onehot = convert_data_str_to_onehot(aligned_seqs)
                dist_matrix[ki, kj] = np.abs(aligned_onehot[0, :, :] - aligned_onehot[1, :, :]).sum()


        dist_matrix += dist_matrix.T
        for i in range(1, initial_n_leaves):
            for j in range(i):
                computed_distances[initial_leaves[i]][initial_leaves[j]] = dist_matrix[i, j]

        dist_matrix = dist_matrix / (2 * seq_length)
        dist_matrix = (-3 / 4) * np.log(1 - dist_matrix * 4 / 3)

        dm = DistanceMatrix(dist_matrix)
        tree_newick = nj(dm, result_constructor=str)
        tree = newick2nx(tree_newick, initial_n_leaves)

        mapping = {}
        for k in range(0, initial_n_leaves):
            mapping[k] = initial_leaves[k]  # 2 * initial_n_leaves  - 2 -
        new_int_node = n_leaves
        for k in range(initial_n_leaves, 2 * initial_n_leaves - 2):
            mapping[k] = new_int_node
            new_int_node += 1
        tree = nx.relabel_nodes(tree, mapping)
        root = new_int_node - 1

        # insert the remaining leaves one by one
        first_insertion = True
        for held_out in held_out_leaves:
            centroidMarked = [False] * MAXN
            stop_node_dict = defaultdict(list)

            if first_insertion:  # we will save the DFS results from first insertion
                selected_adj, other_selected_adj, current_degree = decomposeTree(root, tree, held_out,
                                                                                 first_placement=True)
                first_insertion = False
            else:
                selected_adj, other_selected_adj, current_degree = decomposeTree(root, tree, held_out,
                                                                                 subtree_leaves=first_dfs_subtree_leaves.copy(),
                                                                                 subtree_size=first_dfs_subtree_sizes.copy(),
                                                                                 parent_array=first_parent_array.copy(),
                                                                                 )
            # now we have the two nodes of the selected edge, determine which is the parent of the other
            # that info is required for a correct update on data structures
            if other_selected_adj == first_parent_array[selected_adj]:
                upper_neighbour = other_selected_adj
                lower_neighbour = selected_adj
            else:
                upper_neighbour = selected_adj
                lower_neighbour = other_selected_adj

            # place the held-out leaf on the selected edge, and update data structures
            tree.remove_edge(upper_neighbour, lower_neighbour)
            first_centroid = root_centroid_dict[(root, 0, tuple([]), -1)]

            tree.add_node(new_int_node, parent=upper_neighbour)
            tree.add_edge(upper_neighbour, new_int_node, t=np.random.exponential(0.1))
            tree.add_node(lower_neighbour, parent=new_int_node)
            tree.add_edge(new_int_node, lower_neighbour, t=np.random.exponential(0.1))
            first_dfs_subtree_leaves[upper_neighbour].append(held_out)
            first_dfs_subtree_sizes[upper_neighbour] += 2
            grand_parent = first_parent_array[upper_neighbour]
            while grand_parent != -1:
                first_dfs_subtree_leaves[grand_parent].append(held_out)
                first_dfs_subtree_sizes[grand_parent] += 2
                grand_parent = first_parent_array[grand_parent]

            first_dfs_subtree_leaves[held_out].append(held_out)
            first_dfs_subtree_sizes[held_out] = 1

            first_dfs_subtree_leaves[new_int_node] = first_dfs_subtree_leaves[lower_neighbour].copy()
            first_dfs_subtree_leaves[new_int_node].append(held_out)
            first_dfs_subtree_sizes[new_int_node] = 2 + first_dfs_subtree_sizes[lower_neighbour]

            first_parent_array[lower_neighbour] = new_int_node
            first_parent_array[held_out] = new_int_node
            first_parent_array[new_int_node] = upper_neighbour

            tree.add_node(held_out, parent=new_int_node)
            tree.add_edge(new_int_node, held_out, t=np.random.exponential(0.1))

            new_int_node += 1
            # check the first centroid if it is affected by this insertion
            current_size = tree.number_of_nodes()
            adj_sizes = np.array([first_dfs_subtree_sizes[adj] for adj in tree.adj[first_centroid]])
            fc_adj_list = np.array([adj for adj in tree.adj[first_centroid]])
            if (adj_sizes > current_size / 2).any():  # if so, fix that
                # update the centroid and the data structures
                new_first_centroid = fc_adj_list[np.argmax(adj_sizes)]
                all_leaves = first_dfs_subtree_leaves[first_centroid].copy()

                first_dfs_subtree_sizes[first_centroid] -= first_dfs_subtree_sizes[new_first_centroid]
                first_dfs_subtree_leaves[first_centroid] = list(
                    [node for node in first_dfs_subtree_leaves[first_centroid] if
                     node not in first_dfs_subtree_leaves[new_first_centroid]])

                first_dfs_subtree_sizes[new_first_centroid] = current_size
                first_dfs_subtree_leaves[new_first_centroid] = all_leaves

                first_parent_array[first_centroid] = new_first_centroid
                first_parent_array[new_first_centroid] = -1
                # register the new centroid
                root_centroid_dict[(root, 0, tuple([]), -1)] = new_first_centroid

        our_toc = time.process_time()
        # pickle.dump(tree, open(
        #     'results/' + args.data + '_SNJ_tree_' + str(seed) + '.pickle', 'wb'))

        tree_nw_str = nx_to_newick(tree, np.max(tree.nodes))
        tree_nw = newick.loads(tree_nw_str)
        newick.write(tree_nw, 'results/' + args.data + '_SNJ_tree_' + str(seed) + '.nw')

        time_array[seed] = our_toc - our_tic

    print('*******************')
    print('DONE=>>>>> ' + args.data + ' - SNJ tree(s) ready!' )
    print('Elapsed time(s): ', time_array)
    print('*******************')


if __name__ == "__main__":
    main()


