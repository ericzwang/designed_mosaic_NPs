import tensorflow as tf
import tensorflow.keras as K
import numpy as np
import pandas as pd
import pickle
tf.get_logger().setLevel('ERROR')



def get_output_avg(seqs):
    """Get the average and standard deviation of model outputs for a given set of sequences."""
    outputs = []
    for fold in range(10):
        model = K.models.load_model(f'expression_predictor/model{fold}')
        output = np.reshape(model.predict(seqs), [-1])
        outputs.append(output)
    return np.vstack(outputs).mean(axis=0), np.vstack(outputs).std(axis=0)

def get_stable_indices(group, threshold=-0.2):
    """Get indices of stable sequences for a given group."""
    mosaic_seqs = pickle.load(open(f"mosaic_seqs_group{group}.pkl", "rb"))
    mosaic_outputs_avg, mosaic_outputs_std = get_output_avg(mosaic_seqs)
    indices = np.where(mosaic_outputs_avg > threshold)
    return indices

def rename_mutation_columns_short_group(class_, group):
    """Return a dictionary of renamed short columns based on class and group."""
    names = ['Group1-A', 'Group1-B', 'Group1-C', 'Group2-A', 'Group2-B', 'Group2-C']
    assert group in [1, 2]
    group_names = [name for name in names if str(group) in name]
    return {name: f"Cls{class_}-{name}" for name in group_names}

def column_names_group(group):
    """Return a list of column names for a given group."""
    class1_names = list(rename_mutation_columns_short_group(1, group).values())
    class2_names = list(rename_mutation_columns_short_group(2, group).values())
    return class1_names + class2_names

def run():
    mosaic = pd.read_csv('mosaic-mutations/all.csv')
    group1_ind = get_stable_indices(1)
    group2_ind = get_stable_indices(2)
    stable_ind = np.sort(np.intersect1d(group1_ind, group2_ind)).astype(int)
    stable_mosaic = mosaic.iloc[stable_ind]
    for group in [1, 2]:
        file = f'csv_files/stablegroup{group}.csv'
        stable_mosaic[column_names_group(group)].to_csv(file, index=False)

if __name__ == '__main__':
    run()