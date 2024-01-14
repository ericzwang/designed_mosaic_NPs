import pandas as pd
import sys

# Parameters
# There will be length^2 rows, but filtered down so that residues are not shared
length = 2000
filter_length = 10000

def rename_mutation_columns(class_):
    """Return a dictionary of renamed columns"""
    original_names = ['Group1-A', 'Group1-B', 'Group1-C', 'Group2-A',
                      'Group2-B', 'Group2-C', 'Group1 Frac',
                      'Group2 Frac', 'Average Frac']
    return {name: f"Cls{class_}-{name}" for name in original_names}

def rename_mutation_columns_short(class_):
    """Return a dictionary of renamed short columns"""
    original_names = ['Group1-A', 'Group1-B', 'Group1-C', 'Group2-A',
                      'Group2-B', 'Group2-C']
    return {name: f"Cls{class_}-{name}" for name in original_names}

def combine_classes(i, class1, class2):
    """Get a dataframe of both classes, filtered so that no residues are shared in both classes
    i is the row index"""
    mut_columns = rename_mutation_columns_short(1).tolist() + rename_mutation_columns_short(2).tolist()
    temp = pd.DataFrame(class1.iloc[i:i+1], columns=class1.columns)
    temp2 = pd.concat([temp]*(length)).reset_index().drop(columns='index')
    temp3 = pd.concat([temp2, class2], axis=1)
    temp4 = pd.concat([temp3[name].str[1:-1] for name in mut_columns], axis=1)
    indices = [i for i in range(length) if temp4.iloc[i].value_counts().sum() <= temp4.iloc[i].value_counts().shape[0]]
    return temp3.iloc[indices].reset_index()

def concat_frames(start_row, stop_row, class1, class2):
    """Concatenate dataframes using the indices of start and stop rows, then save to file"""
    frame_list = [combine_classes(i, class1, class2) for i in range(start_row, stop_row)]
    frames = pd.concat(frame_list, ignore_index=True).drop(columns='index')
    frames.to_csv(f'csv_files/{start_row}_{stop_row}.csv', index=False)

def read_class_csv(ab_class, length):
    """Read class CSV file and rename columns"""
    df = pd.read_csv(f'csv_files/class{ab_class}.csv').iloc[0:length]
    return df.rename(columns=rename_mutation_columns(ab_class))

def format_class_df(df, length, filter_length, ab_class):
    """Format class dataframe"""
    tmp = df.iloc[0:filter_length].sample(length)
    tmp2 = tmp.sort_values(by=f'Cls{ab_class}-Average Frac', ascending=False)
    return tmp2.reset_index().drop(columns='index')

def run():
    # Read and format class dataframes
    class1_original = read_class_csv(1, length)
    class2_original = read_class_csv(2, length)
    class1 = format_class_df(class1_original, length, filter_length, 1)
    class2 = format_class_df(class2_original, length, filter_length, 2)

    # Concatenate frames for specified row range
    for start_row in [int(sys.argv[1])]:
        stop_row = start_row + 50
        concat_frames(start_row, stop_row, class1, class2)

if __name__ == '__main__':
    run()

    