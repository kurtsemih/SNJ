# SNJ
Sparse Neighbor Joining: rapid phylogenetic inference using a sparse distance matrix.

## Workflow
Please see pseudo code:

![snj_pseudo_code](https://github.com/kurtsemih/SNJ/assets/68644299/4ecdd511-a9c7-430e-bc8e-93f2e93786f8)


## Installation

For installation, simply download this repo and run the following commands. **Python >= 3.8** is recommended for running SNJ.

    cd SNJ/
    pip install -e .

## Usage

The input data must be in `.npy` format and placed in `./data`.

### Hyperparameters

`-seed` Number of seeds. 

`-n_i` Number of initial leaves.

`-n_s` Number of sampled leaves.

`-n_o` Number of orienting leaves.

 #### Example usage with test data:
 
 The following command will run the phylogenetic inference on test data.

`snj -data test -seed 1 -n_i 10 -n_s 5 n_o 3`

## Outputs
SNJ produces the following files in `./results`:

`inputdata_SNJ_tree_#.pickle` The inferred phylogenetic tree for `inputdata` and seed `#`. 

It is possible to obtain the tree in `.nx` format via:
 `pickle.load(open(inputdata_SNJ_tree_#.pickle, 'rb'))`


See [Wiki](https://github.com/kurtsemih/copyVAE/wiki) for more details.
