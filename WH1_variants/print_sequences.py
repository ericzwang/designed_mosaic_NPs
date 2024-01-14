import pandas as pd
import numpy as np
import pickle



def rename_mutation_columns_short_group(class_, group):
    """Return a dictionary of renamed short columns based on class and group."""
    names = ['Group1-A', 'Group1-B', 'Group1-C', 'Group2-A', 'Group2-B', 'Group2-C']
    if group not in [1, 2]:
        raise ValueError('Group must be 1 or 2')
    group_names = [name for name in names if str(group) in name]
    return {name: f"Cls{class_}-{name}" for name in group_names}

def column_names_group(group):
    """Return a list of column names for a given group."""
    class1_names = list(rename_mutation_columns_short_group(1, group).values())
    class2_names = list(rename_mutation_columns_short_group(2, group).values())
    return class1_names + class2_names

def make_soluble_mosaic(threshold=-40.5916):
    """Filter soluble indices and create a soluble mosaic."""
    scores_indices = pd.read_csv('csv_files/solubleindices.csv', header=None)
    scores_indices_pivot = scores_indices.pivot(index=0, columns=1, values=2).iloc[1:]
    scores_indices_pivot_filter = scores_indices_pivot[scores_indices_pivot[1] < threshold]
    scores_indices_pivot_filter = scores_indices_pivot_filter[scores_indices_pivot_filter[2] < threshold]
    soluble_indices = scores_indices_pivot_filter.index.astype(int)

    mosaic = pd.read_csv('csv_files/all.csv')
    indices = np.loadtxt('csv_files/stableindices.csv')
    stable_mosaic = mosaic.iloc[indices]
    soluble_mosaic = stable_mosaic.iloc[soluble_indices].reset_index().drop(columns='index')
    soluble_mosaic.to_csv('csv_files/solublemosaic.csv', index=False)
    return soluble_mosaic

def get_entropy(group, sample):
    """Calculate entropy for a given group in the sample."""
    tot_entropy = 0
    for i, col in sample[column_names_group(group)].items():
        prob = np.divide(col.value_counts().values, col.value_counts().values.sum())
        entropy = np.sum(-1 * np.multiply(prob, np.log(prob)))
        tot_entropy += entropy
    return tot_entropy

def group_entropy(sample):
    """Calculate the average entropy for both groups in the sample."""
    return np.mean([get_entropy(1, sample), get_entropy(2, sample)])

def print_seqs(soluble_mosaic):
    """Print pairs of sequences with the highest group entropy."""
    group_entropies = []
    for seed in range(1000):
        np.random.seed(seed)
        sample = soluble_mosaic.iloc[0:20000].sample(5)
        group_entropies.append([sample, group_entropy(sample)])
    pickle.dump(group_entropies, open('pkl/group_entropies.pkl', 'wb'))

    dfs = [item[0] for item in group_entropies]
    best_muts = dfs[np.argmax([item[1] for item in group_entropies])]

    group1 = [list(row) for i, row in enumerate(best_muts[column_names_group(1)].to_numpy())]
    group2 = [list(row) for i, row in enumerate(best_muts[column_names_group(2)].to_numpy())]

    for i in range(len(group1)):
        print(f"Pair{i + 1}")
        print(", ".join(group1[i]))
        print(", ".join(group2[i]))

def run():
    soluble_mosaic = make_soluble_mosaic()
    print_seqs(soluble_mosaic)

if __name__ == '__main__':
    run()