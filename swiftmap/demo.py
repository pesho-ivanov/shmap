import sketch
from Bio import SeqIO
import sys
import random
from tqdm import tqdm


def test_random():
    # Random sequence of A T G C
    #T = ''.join(random.choice('ATGC') for _ in range(100)).encode('utf-8')
    T = b'AGCAAGAAGTGGTTACTAGATCGCACTCGGAGGTAGCAACTAGACGGAAACGGGGCATGCCAAGGGTGTGTTTTTATTGAGTTTAGATTAGAGTCGGTAT'
    S = T[10:10+20]
    K = 5
    h_frac = 1

    T_sketch = sketch.Sketch(T, K, h_frac)
    T_index = sketch.Index()
    T_index.add(T_sketch)

    S_sketch = sketch.Sketch(S, K, h_frac)

    res = sketch.sweep_map(T_index, S_sketch, max_hashes=20, width=len(S)-K+1)
    print(f'width={len(S)-K+1}')
    print(res)

def test_real_data_one_read():
    ref = next(SeqIO.parse("../evals/refs/t2tChrY.fa", format="fasta"))
    #read = list(SeqIO.parse("../evals/reads/t2tChrY-a0.99-d1-l10000.oneread.fa", format="fasta"))[0]
    reads = list(SeqIO.parse("../evals/reads/t2tChrY-a0.99-d1-l10000.fa", format="fasta"))
    #reads_bytes = [bytes(r.seq) for r in reads]

    read = reads[17]

    T = bytes(ref.seq)
    P = bytes(read.seq)
    K = 22
    R = 0.1
    S = 300

    T_sketch = sketch.Sketch(T, K, R)
    T_index = sketch.Index()
    T_index.add(T_sketch)

    P_sketch = sketch.Sketch(P, K, R)
    res = sketch.sweep_map(T_index, P_sketch, max_hashes=S, width=len(P)-K+1)

    print(res, read.name)


def test_real_data():
    ref = next(SeqIO.parse("../evals/refs/t2tChrY.fa", format="fasta"))
    reads = list(SeqIO.parse("../evals/reads/t2tChrY-a0.99-d1-l10000.fa", format="fasta"))
    reads_bytes = [bytes(r.seq) for r in reads]

    T = bytes(ref.seq)
    K = 22
    R = 0.1
    S = 300
    M = 100

    T_sketch = sketch.Sketch(T, K, R)
    T_index = sketch.Index()
    T_index.add(T_sketch, max_matches=M)

    mapping = []
    overlaps = []
    for read, P in tqdm(zip(reads, reads_bytes)):
        res = sketch.sweep_map(T_index, sketch.Sketch(P, K, R), max_hashes=S, width=len(P)-K+1) 
        map_left, I = res
        map_right = map_left + len(P)

        name_parts = read.name.split('!')
        correct_left, correct_right = int(name_parts[2]), int(name_parts[3])
        overlap = max(0, min(correct_right, map_right) - max(correct_left, map_left))
        rel_overlap = overlap / (correct_right - correct_left)
        overlaps.append(rel_overlap)

        mapping.append(res)

    import numpy as np
    
    overlaps = np.array(overlaps)
    print(f'Average overlap: {np.mean(overlaps)}')
    print(f'Precision: {np.mean(overlaps > 0.1)}')




if __name__ == '__main__':
    #test_random()
    test_real_data()