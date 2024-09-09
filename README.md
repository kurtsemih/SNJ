# SNJ
Sparse Neighbor Joining: rapid phylogenetic inference using a sparse distance matrix.


## Installation

For installation, simply download this repo and run the following commands. **Python >= 3.8** is recommended for running SNJ.

    cd SNJ/
    pip install -e .

## Usage

The input data must be in `.npy` format and placed in `./data`.

Use `-data` flag to input your data.

E.g., `-data inputdata` for `./data/inputdata.npy`

### Hyperparameters

`-seed` Number of seeds. 

`-n_i` Number of initial leaves. Recommended: `sqrt(n log n)` for `n` taxa

`-n_s` Number of sampled leaves. Recommended: `log n` for `n` taxa

`-n_o` Number of orienting leaves. Recommended: `3`

 #### Example usage with test data:
 
 The following command will run the phylogenetic inference on test data:

`snj -data test -seed 1 -n_i 10 -n_s 5 -n_o 3`

If the input data is unaligned, then use `snj_pw` command: 

`snj_pw -data test -seed 1 -n_i 10 -n_s 5 -n_o 3`

## Outputs
SNJ produces the following files in `./results`:

`inputdata_SNJ_tree_#.nw` The inferred phylogenetic tree for `inputdata` and seed `#`, in `.nw` format . 


See [Wiki](https://github.com/kurtsemih/copyVAE/wiki) for more details.
