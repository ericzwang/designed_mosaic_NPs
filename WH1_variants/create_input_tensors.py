import pandas as pd
import numpy as np
import tensorflow as tf
import pickle



def constants():
    """Return constants including amino acid order and single mutation effects."""
    aa_order = ['*', 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K',
                'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    single_muts = pd.read_csv('inputs/mutdata.csv')[['site_RBD', 'wildtype', 'mutant', 'expr_avg']]
    single_mut_tf = tf.constant(single_muts.pivot(
        index='site_RBD', columns='mutant', values='expr_avg'
    ).replace(np.nan, -5))
    return aa_order, single_mut_tf

def rename_mutation_columns_short_group(class_, group):
    """Return a dictionary of renamed short columns based on class and group."""
    names = ['Group1-A', 'Group1-B', 'Group1-C', 'Group2-A', 'Group2-B', 'Group2-C']
    if group not in [1, 2]:
        raise ValueError('Group must be 1 or 2')
    group_names = [name for name in names if str(group) in name]
    return {name: f"Cls{class_}-{name}" for name in group_names}

def get_product_tfs(list_of_muts, single_mut_tf):
    """Return Tensors for a list of mutations."""
    aa_order = constants()[0]
    RBD_length = 201
    num_aas = 21

    product_tfs = []
    for i in range(len(list_of_muts)):
        mutation_tf = np.zeros((RBD_length, num_aas))
        for mutation in list_of_muts[i]:
            wt = mutation[0]
            site = int(mutation[1:-1])
            mutant = mutation[-1]
            aa_index = aa_order.index(mutant)
            mutation_tf[site - 1, aa_index] = 1
        mutation_tf = tf.constant(mutation_tf)
        product_tf = tf.multiply(mutation_tf, single_mut_tf)
        product_tfs.append(product_tf)
    product_tfs = tf.stack(product_tfs)
    return product_tfs

def create_tensor(group, mosaic):
    """Create and save a tensor for a given group and mosaic data."""
    aa_order, single_mut_tf = constants()
    class1_columns = list(rename_mutation_columns_short_group(1, group).values())
    class2_columns = list(rename_mutation_columns_short_group(2, group).values())
    mut_columns = class1_columns + class2_columns
    mosaic_muts_start = mosaic[mut_columns].values.tolist()
    mosaic_muts = [
        [mut[0] + str(int(mut[1:-1]) - 330) + mut[-1] for mut in seq]
        for seq in mosaic_muts_start
    ]
    mosaic_seqs = get_product_tfs(mosaic_muts, single_mut_tf)
    pickle.dump(mosaic_seqs, open(f"pkl/mosaic_seqs_group{group}.pkl", "wb"), protocol=4)

def run():
    mosaic = pd.read_csv('mosaic-mutations/all.csv')
    create_tensor(1, mosaic)
    create_tensor(2, mosaic)

if __name__ == '__main__':
    run()

