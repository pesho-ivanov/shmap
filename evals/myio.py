import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from pathlib import Path
from IPython.display import display, HTML
from Bio import SeqIO
import os
#from astropy import units as u
import sys
from readpaf import parse_paf
from collections import Counter, defaultdict
from tqdm import tqdm
import pickle

#import sv

#from itables import init_notebook_mode
#init_notebook_mode(all_interactive=True)
pd.set_option('display.max_columns', None)


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

def to_latex(df, data, refname):
    latex = ""
    df.index = df.index.map(lambda x: f'\\{x}')
    df.columns = df.columns.str.replace(' ', '\\\\')
    df.columns = df.columns.str.replace('%', '\%')
    df.columns = df.columns.map(lambda x: '\makecell{' + x + '}')
    #df = df.astype(str).map(lambda x: x.rstrip('0').rstrip('.') if '.' in x else x)
    latex += df.to_latex(escape=False, label=f'tab:{refname}', caption=data, float_format = lambda x: '{:0.2f}'.format(x) if pd.notna(x) else '-')
    #latex += df.to_latex(float_format = lambda x: '{:0.2f}'.format(x) if pd.notna(x) else '-')
    latex += '\n'
    return latex

def plot_mapping_quality(df, ax, tool, reads):
    marker = { 'blend': 'o', 'minimap': 's', 'winnowmap': 'x', 'shmap': 'd', 'mm2-fast': 'v', 'mapquik': 'p' } 

    df_sorted = df.sort_values(by="mapping_quality", ascending=False).reset_index(drop=True)
    #display(df_sorted.head(10))

    total_reads = len(df)
    fraction_mapped = []
    error_rates = []
    wrongs = 0
    
    for i in range(total_reads):
        if df_sorted.iloc[i]["is_correct"] == False:
            wrongs += 1
        if i == total_reads-1 or df_sorted.iloc[i]["mapping_quality"] != df_sorted.iloc[i+1]["mapping_quality"]:
            error_rates.append(wrongs / (i+1))
            fraction_mapped.append((i+1) / reads)

    #print(reads, total_reads, wrongs)
    #print(error_rates, fraction_mapped)

    # Plot the results
    tool_marker = marker[ tool.split('-')[0] ]
    ax.plot(error_rates, fraction_mapped, marker=tool_marker, linestyle='-', label=tool)
    plt.xlim(1e-5, 1e-1)
    plt.ylim(0.9, 1)
    ax.set_xscale('log')  # Log scale for error rate
    ax.set_xlabel("Error rate of mapped reads")
    ax.set_ylabel("Fraction of mapped reads")
    ax.set_title("Fraction of Mapped Reads vs Error Rate")
    
    # Customize x-axis labels to show 10^-i without coefficients
    def log_format(x, pos):
        exponent = int(np.log10(x))
        return f'$10^{{{exponent}}}$'
    
    ax.xaxis.set_major_formatter(FuncFormatter(log_format))

    ax.grid(True, which="both", linestyle="--", linewidth=0.5)

    #def perc(a, b):
#    if b == 0:
#        return np.nan
#    return 100.0 * a / b

def fasta2df(fn):
    seqs = SeqIO.parse(fn, "fasta")
    df = pd.DataFrame((str(s.id), str(s.seq)) for s in seqs)
    df.columns = ["ID", "Sequence"]
    return df

#def is_overlapping(a, sv_row):
#    return a.GT_from <= sv_row['END'] and sv_row['POS'] <= a.GT_to 
    
min_overlap = 0.1

def get_overlap(a):
     if a.GT_ref != a.target_name:
         return False
     if a.GT_strand != a.strand:
         return False
     union_from = min(a.GT_from, a.target_start)
     union_to = max(a.GT_to, a.target_end)

     intersect_from = max(a.GT_from, a.target_start)
     intersect_to = min(a.GT_to, a.target_end)
     overlap = max(0.0, (intersect_to - intersect_from) / (union_to - union_from))
     return overlap

def read_paf(pref, reads, experiment, tool, ax):
    paf_file = pref.parent / (pref.name + '.paf')
    serialize_file = pref.parent / (pref.name + '.pkl')

    if serialize_file.exists():
        try:
            with open(serialize_file, 'rb') as f:
                res = pickle.load(f)
                #df_max = pickle.load(f)
                #plot_mapping_quality(df_max, ax, tool, reads)
                return res
        except Exception as e:
            print(f"Error reading pickle file {serialize_file}: {e}. Recomputing.")

    print(f"Computing {tool} for {experiment}...")

    #print(paf_file)
    no_GT = False
    if not paf_file.exists():
        raise Exception(f"File does not exist or is empty: {paf_file}")
    with open(paf_file) as handle:
        df = parse_paf(handle, dataframe=True)
        df['experiment'] = experiment
        df['tool'] = tool
        try:
            df[ ['read_name', 'GT_ref', 'GT_from', 'GT_to', 'GT_strand'] ] = df['query_name'].str.split('!', expand=True)
            df['GT_from'] = df['GT_from'].astype(int)
            df['GT_to'] = df['GT_to'].astype(int)
            df['overlap'] = df.apply(get_overlap, axis=1)
            df['is_correct'] = df['overlap'] >= min_overlap
            #df['is_correct_labels'] = df.apply(lambda x: is_correct_labels(x, orig_l, mutated_l), axis=1)
            #df['is_correct'] = df.apply(lambda x: is_correct_labels(x, orig_l, mutated_l), axis=1)
            #df['start_diff'] = df.target_start - df.GT_from  # TODO: different coordinate systems!
            #df['end_diff'] = df.target_end - df.GT_to  # TODO: different coordinate systems!
            #df['read_sv'] = 'none' # df.apply(lambda x: read_falls_on_what_sv(x, vcf_df), axis=1)
        except Exception as e:
            #display(e)
            df['read_name'] = df['query_name']
            df['GT_ref'] = np.NaN
            df['GT_from'] = np.NaN
            df['GT_to'] = np.NaN
            df['GT_strand'] = np.NaN
            df['overlap'] = np.NaN
            df['is_correct'] = True
            #df['is_correct_labels'] = True
            #df['start_diff'] = 0
            #df['end_diff'] = 0
            #df['read_sv'] = 'none'
            no_GT = True
    df = df.sort_values(['read_name', 'residue_matches'], ascending=[True, False], ignore_index=True)
    #display(df)

    df = df[['read_name', 'mapping_quality', 'is_correct', 'alignment_block_length']]  # speeding up and freeing memory

    paf = defaultdict(int)
    paf['Mapped Q60'] = 0
    paf['Q<60 or missed'] = np.nan
    paf['Wrong Q60'] = 0
    #mapped_reads = 0
    #wrong = 0

    #def process_group(group_first_index, group_last_index):
    #    nonlocal mapped_reads, wrong
    #    group = df.loc[group_first_index:group_last_index]
    #    mapped_reads += 1
    #    if (group.mapping_quality == 60).all():
    #        paf['Mapped Q60'] += 1
    #        if not group.is_correct.all():
    #            paf['Wrong Q60'] += 1
    #    else:
    #        wrong += 1

    #group_first_i, group_read_name = 0, df.loc[0, 'read_name']
    #for i, a in df.iterrows():
    #    if a.read_name != group_read_name:
    #        process_group(group_first_i, i-1)
    #        group_first_i, group_read_name = i, a.read_name
    #process_group(group_first_i, len(df)-1)

    #df['is_q60'] = (df['mapping_quality'] == 60)
    
    df_max = df.loc[df.groupby('read_name')['alignment_block_length'].idxmax()]
    df_max = df_max.reset_index(drop=True)

    #grouped = df.groupby('read_name').agg(
    #    all_q60=('is_q60', 'all'),       # True if every row in the group is Q60
    #    all_correct=('is_correct', 'all')# True if every row in the group is correct
    #)
    mapped_reads = len(df_max)  # total number of read_name groups
    mapped_q60 = (df_max['mapping_quality'] == 60).sum()
    mapped_ql60 = (df_max['mapping_quality'] < 60).sum()
    wrong_q60 = ((df_max['mapping_quality'] == 60) & (df_max['is_correct'] == False)).sum()
    #wrong = (df_max['is_correct'] == False).sum()

    paf['Mapped Q60'] = int(mapped_q60)
    paf['Wrong Q60']  = int(wrong_q60)

    missed = reads - mapped_reads
    Ql60_or_missed = mapped_ql60 + missed
    #print('reads=', reads, 'mapped_reads=', mapped_reads, 'mapped_q60=', mapped_q60, 'missed=', missed, 'Ql60_or_missed=', Ql60_or_missed)
    #paf['Q<60 or missed'] = wrong_or_missed
    paf['Q<60 or missed'] = "{:.1f}%".format(100 * Ql60_or_missed / reads)

    if no_GT:
        paf['Wrong Q60'] = 'n/a'

    #plot_mapping_quality(df_max, ax, tool, reads)

    res = pd.Series(paf, dtype='object')
    with open(serialize_file, 'wb') as f:
        pickle.dump(res, f)
        #pickle.dump(df_max, f)

    return res
    
index_time_col = 'Index [sec]'
map_time_col = 'Map [sec]'
memory_col = 'Memory [GB]'

def read_times(pref):
    times = {}
    with open(str(pref) + '.index.time') as f_index_time:
        index_time, index_mem = map(float, f_index_time.readline().split())
        times[index_time_col] = index_time #* u.second
        #times['index_mem'] = index_mem / 2**20
        with open(str(pref) + '.time') as f_time:
            total_time, total_mem = map(float, f_time.readline().split())
            #times['time total'] = total_time #* u.second
            times[map_time_col] = total_time - times[index_time_col]
            times[memory_col] = (total_mem / 2**20) #* u.GB
    return pd.Series(times, dtype='object')

def get_comparison_table(main_dir: Path, refname, experiment: Path, tools):
    empty_cell = -1
    alldf = pd.DataFrame()
    alldf.name = experiment
    ref = fasta2df(Path('refs') / (refname + '.fa'))
    reads = fasta2df(Path('reads') / (str(experiment) + '.fa'))

    rows = []
    fig, ax = None, None
    #fig, ax = plt.subplots()

    for tool in tqdm(tools, desc=f'Tools for {experiment}', leave=False):
        d = Path(main_dir) / tool / experiment / tool
        row = pd.Series({
            'tool': tool,
        })
        try:
            paf = read_paf(d, len(reads), experiment, tool, ax)
        except Exception as e:
            print(f"An error occurred while reading PAF {d}: {e}")
            continue
        row = pd.concat([row, paf])
        try:
            row = pd.concat([row, read_times(d).map('{:.1f}'.format)])   # .time, .index.time
        except Exception as e:
            print(f"An error occurred while reading times {d}: {e}")
            row[index_time_col] = empty_cell
            row[map_time_col] = empty_cell
            row[memory_col] = empty_cell
        rows.append(row)
        #except Exception as e:
        #    print(f"An error occurred while reading PAF {d}: {e}")
    #ax.legend()
    #plt.show()

    alldf = pd.DataFrame(rows)
    alldf = alldf.set_index('tool')
    alldf.index.name = None
    return alldf

def highlight_min(s):
    is_min = s == s.min()
    return ['font-weight: bold' if v else '' for v in is_min]

def highlight_max(s):
    is_max = s == s.max()
    return ['font-weight: bold' if v else '' for v in is_max]
