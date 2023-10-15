import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re

def parse_fasta_metadata(fasta_file):
    with open(fasta_file, 'r') as f:
        lines = f.readlines()

    metadata_lines = [line.strip() for line in lines if line.startswith('>')]

    data = []
    for line in metadata_lines:
        parts = line[1:].split()  # Remove the '>' and split by whitespace
        match = re.search(r':(\d+)-(\d+)', parts[2])
        assert(match)
        fr, to = int(match.group(1)), int(match.group(2))
        len = to - fr
        assert(len > 0)

        entry = {
            'query_name': parts[0],
            'from_ref_sim': fr,
            'to_ref_sim': to,
            #'sr': float(parts[4]),
            #'ir': float(parts[6]),
            #'dr': float(parts[8]),
            #'sd': int(parts[10])
            'read_len': len, 
        }
        data.append(entry)

    return pd.DataFrame(data)

def plot_all_columns(df):
    n = len(df.columns)
    # Calculate the number of rows and columns for the grid
    rows = int(np.ceil(np.sqrt(n)))
    cols = int(np.ceil(n / rows))
    
    fig, axes = plt.subplots(rows, cols, figsize=(15, 15))
    # If only one row, make sure axes is 2D array
    if rows == 1:
        axes = np.reshape(axes, (1, -1))
    
    for i, column in enumerate(df.columns):
        ax = axes[i // cols, i % cols]
        # If the column is numeric, plot a histogram
        if pd.api.types.is_numeric_dtype(df[column]):
            df[column].hist(ax=ax, bins=20)
            ax.set_title(f'Histogram of {column}')
            ax.set_xlabel(column)
            ax.set_ylabel('Frequency')
        # If the column is categorical and has fewer than 20 unique values, plot a bar chart
        elif df[column].nunique() < 20:
            df[column].value_counts().plot(kind='bar', ax=ax)
            ax.set_title(f'Bar plot of {column}')
            ax.set_ylabel('Count')
    
    # Remove any unused subplots
    for j in range(i+1, rows*cols):
        fig.delaxes(axes.flatten()[j])
    
    plt.tight_layout()
    plt.show()

# the intersection divided by the union of the reference intervals 
def get_jaccard_nucl_overlap(row1, row2):
    assert(row1['query_name'] == row2['query_name'])
    if (row1['ref_start'] > row2['ref_end']) or (row2['ref_start'] > row1['ref_end']):
        return 0
    intersection = (min(row1['ref_end'], row2['ref_end']) - max(row1['ref_start'], row2['ref_start']))
    union = (max(row1['ref_end'], row2['ref_end']) - min(row1['ref_start'], row2['ref_start']))
    if union == 0:
        return 0
    return intersection / union

# calculate the maximal jaccard of a row with all other rows with the same query_name in the dataframe
def get_max_jaccard(row, groundtruth_df):
    common = groundtruth_df[groundtruth_df['query_name'] == row['query_name']]
    if common.shape[0] == 0:
        return -1
    return common.apply(lambda row2: get_jaccard_nucl_overlap(row, row2), axis=1).max()

# add a column to the tested_df with the maximal intersection with a groundtruth row with the same query_name
def add_overlap_column(tested_df, groundtruth_df, intersection_column_name):
    tested_df[intersection_column_name] = tested_df.apply(lambda row: get_max_jaccard(row, groundtruth_df), axis=1)

# get the accuracy
def get_accuracy(tested_df, groundtruth_df):
    add_overlap_column(tested_df, groundtruth_df, 'overlap')
    accuracy = tested_df[tested_df['overlap'] > 0.1].shape[0] / tested_df.shape[0]
    return accuracy