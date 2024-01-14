import pandas as pd
import numpy as np
import itertools
from copy import deepcopy



def load_data(file='inputs/mutdata.csv'):
    """
    Load data from a CSV file and extract relevant columns.
    Returns mutations, residues, mutant amino acids, and amino acid mapping.
    """
    data = pd.read_csv(file)
    stable = data[data["expr_avg"] > -0.2]
    muts = stable["mutation"]
    resi = stable["site_SARS2"]
    mutaa = stable["mutant"]
    aa_map = {"G": 0, "A": 0, "V": 0, "C": 0, "P": 0, "L": 0, "I": 0, "M": 0, "W": 0, "F": 0,
              "K": 1, "R": 1, "H": 1,
              "D": 2, "E": 2,
              "S": 3, "T": 3, "Y": 3, "N": 3, "Q": 3,
              "X": 4, "J": 5, "B": 6, "Z": 7}
    return muts, resi, mutaa, aa_map

def get_escape_dataframe(ab_class, file='inputs/escape_data.csv'):
    """
    From a CSV file of DMS data, get a DataFrame of escape fractions for each mutation
    against antibodies of a certain class.
    """
    escape_data = pd.read_csv(file)
    escape_data = escape_data[escape_data['condition_subtype'].isin([f'class {ab_class}'])]
    escape_lst = []

    for condition in escape_data['condition'].unique():
        df_condition = escape_data[escape_data["condition"] == condition]
        df_condition = df_condition.sort_values(by=['mut_escape'], ascending=False)
        for ind, row in df_condition.iterrows():
            mutation = row["wildtype"] + str(row["site"]) + row["mutation"]
            escape = row["mut_escape"]
            site = row['site']
            escape_lst.append([condition, mutation, escape, site])
    df = pd.DataFrame(escape_lst)
    df_pivot = df.pivot(index=1, columns=0, values=2)
    df_pivot = df_pivot.div(df_pivot.sum(axis=0), axis=1)  # normalize each column
    return df_pivot, df

def generate_groups(lst, n):
    """Return all the ways to split a list into groups of size n."""
    if not lst:
        yield []
    else:
        for group in (((lst[0],) + xs) for xs in itertools.combinations(lst[1:], n-1)):
            for groups in generate_groups([x for x in lst if x not in group], n):
                yield [list(group)] + groups

def iterate_residues(all_comb, escape_fractions):
    """
    Go through list of residues, turn them into mutations, and determine escape fractions.
    """
    muts, resi, _, aa_map = load_data()
    all_muts = deepcopy(all_comb)
    all_fracs = np.zeros((len(all_muts), 3))

    for i, combination in enumerate(all_comb):
        fractions = []
        for j, sequence in enumerate(combination):
            for k, residue in enumerate(sequence):
                mutations = escape_fractions[escape_fractions.index.isin(
                    muts[resi == residue].values
                )].sort_values(ascending=False)
                for ind, escape_frac in enumerate(mutations):
                    mutation = mutations.index[ind]
                    wildtype_aa = mutation[0]
                    mutant_aa = mutation[-1]
                    if not (aa_map[wildtype_aa] in [1, 2] and aa_map[mutant_aa] == 0):
                        fractions.append(escape_frac)
                        all_muts[i][j][k] = mutation
                        break
        fractions = np.array(fractions)
        all_fracs[i] = np.array([fractions[0:2].mean(), fractions[3:5].mean(), fractions.mean()]) * 1000
    return all_muts, all_fracs

def make_sequences_df(all_muts, all_fracs):
    """Create dataframe from lists."""
    temp = pd.DataFrame(all_muts)
    df1 = pd.DataFrame(temp[0].to_list(), columns=['Group1-A', 'Group1-B', 'Group1-C'])
    df2 = pd.DataFrame(temp[1].to_list(), columns=['Group2-A', 'Group2-B', 'Group2-C'])
    df3 = pd.DataFrame(all_fracs)
    sequences_df = pd.concat([df1, df2, df3], axis=1).rename(
        columns={0: 'Group1 Frac', 1: 'Group2 Frac', 2: 'Average Frac'}
    )
    return sequences_df

def make_csv_files(ab_class):
    """
    Create class CSV files and write them out
    """
    df_pivot, df = get_escape_dataframe(ab_class)
    escape_fractions = df_pivot.mean(axis=1)
    escape_residues = df[[2, 3]].groupby([3]).sum().sort_values(
        by=[2], ascending=False
    ).iloc[0:20].index.to_numpy()
    iter_comb = list(itertools.combinations(escape_residues, 6))
    all_comb = sum([list(generate_groups(comb, 3)) for comb in iter_comb])

    all_muts, all_fracs = iterate_residues(all_comb, escape_fractions)

    sequences_df = make_sequences_df(all_muts, all_fracs)
    sequences_df.to_csv(f'csv_files/class{ab_class}.csv', index=False)

def run():
    make_csv_files(1)
    make_csv_files(2)

if __name__ == '__main__':
    run()
