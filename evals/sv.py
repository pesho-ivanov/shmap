import vcfpy
import pandas as pd
import itertools
from Bio import SeqIO
from Bio.Seq import Seq
from collections import OrderedDict

import warnings
from vcfpy.exceptions import FieldInfoNotFound

warnings.filterwarnings('ignore', '.*dup_num.*', category=FieldInfoNotFound)

def read_vcf(vcf_file_path):
    reader = vcfpy.Reader.from_path(vcf_file_path)
    records = []

    for record in reader:
        row = {
            'CHROM': record.CHROM,
            'POS': record.POS,
            'ID': record.ID,
            'REF': record.REF,
            'ALT': ','.join([str(alt.value) for alt in record.ALT]),
            'QUAL': record.QUAL,
            'FILTER': ','.join(record.FILTER) if record.FILTER else None,
            # FORMAT and SAMPLE fields could be added here similarly
        }
        for key, value in record.INFO.items():
            if isinstance(value, list):
                value = ','.join(map(str, value))
            row[key] = value
        
        records.append(row)

    return pd.DataFrame(records)

def read_fasta_file(fasta_file):
    sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences[record.id] = str(record.seq)
    return sequences

def are_equal(fa1: dict, fa2: dict):
    if fa1.keys() != fa2.keys():
        print('different chromosomes:')
        print('fa1: ' + ' '.join(fa1.keys()))
        print('fa2: ' + ' '.join(fa2.keys()))
        return False
    for ch in fa1.keys():
        if len(fa1[ch]) != len(fa2[ch]):
            print(f'chromosome {ch} has different lenghts:')
            print(len(fa1[ch]))
            print(len(fa2[ch]))
            return False
        if fa1[ch] != fa2[ch]:
            print(f'chromosome {ch} is different')
            return False
    return True

def sum_len(d: dict):
    return sum(len(value) for value in d.values())

def mutate(fa, l, vcf_df):
    assert(sum_len(fa) == sum_len(l))
    fa = OrderedDict(sorted(fa.items()))
    res_fa = {}
    res_l = {}
    new_labels = 1
    vcf_df = vcf_df.sort_values(by=['CHROM', 'POS'])
    display(vcf_df)

    iter_fa = iter(fa.items())
    fa_chrom, fa_seq = next(iter_fa)
    curr_pos = 0
    for index, sv in itertools.chain(vcf_df.iterrows(), [('!!fake_index!!', pd.Series({'CHROM': '!!fake_chrom!!'}))]):
        try:
            assert(sum_len(res_fa) == sum_len(res_l))
            display(sv)
            while sv.empty or fa_chrom != sv.CHROM:
                if fa_chrom not in res_fa:
                    res_fa[fa_chrom] = ''
                    res_l[fa_chrom] = []
                print(f'1. add {curr_pos}: of len {len(fa_seq[curr_pos:])}')
                res_fa[fa_chrom] += fa_seq[curr_pos:]
                res_l[fa_chrom] += l[fa_chrom][curr_pos:]
                fa_chrom, fa_seq = next(iter_fa)
                curr_pos = 0
            assert(fa_chrom == sv.CHROM)
            if sv.CHROM not in res_fa:
                res_fa[sv.CHROM] = ''
                assert(sv.CHROM not in res_l)
                res_l[sv.CHROM] = []
            print(f'2. add {curr_pos}:{sv.POS} of len {len(fa_seq[curr_pos:sv.POS])}')
            assert(curr_pos <= sv.POS)
            res_fa[sv.CHROM] += fa_seq[curr_pos:sv.POS]
            res_l[sv.CHROM] += l[sv.CHROM][curr_pos:sv.POS]
            curr_pos = sv.POS
            if sv.SVTYPE == 'DUP':
                segm = fa_seq[sv.POS:sv.END]
                print(f'3. dup at {sv.POS} {len(segm)}*{int(sv.dup_num)} of len {len(segm)*int(sv.dup_num)}')
                res_fa[sv.CHROM] += segm*int(sv.dup_num)
                res_l[sv.CHROM] += l[sv.CHROM][sv.POS:sv.END]*int(sv.dup_num)
                curr_pos = sv.POS
            elif sv.SVTYPE == 'INS':
                print(f'4. ins of len {len(sv.ALT)}')
                res_fa[sv.CHROM] += sv.ALT
                res_l[sv.CHROM] += range(-new_labels, -new_labels-len(sv.ALT), -1)
                new_labels += len(sv.ALT)
                curr_pos = curr_pos  # nothing changes
            elif sv.SVTYPE == 'DEL':
                if sv.REF != fa_seq[sv.POS:sv.END]:
                    print(sv.REF)
                    print(fa_seq[sv.POS:sv.END])
                    assert(False)
                print(f'5. delete {sv.POS}:{sv.END} of len {len(fa_seq[sv.POS:sv.END])}')
                curr_pos = sv.END
            elif sv.SVTYPE == 'INV':
                segm = Seq(fa_seq[sv.POS:sv.END])
                print(f'6. inv {sv.POS}:{sv.END} of len {len(segm)}')
                res_fa[sv.CHROM] += str(segm.reverse_complement())
                res_l[sv.CHROM] += l[sv.CHROM][sv.POS:sv.END][::-1]
                curr_pos = sv.END
            else:
                print(sv.SVTYPE, 'is not a supported SV type')
                assert(False)
                return Null
        except StopIteration:
            break
    assert(len(fa) == len(res_fa))
    print(sum_len(res_fa), sum_len(res_l))
    assert(sum_len(res_fa) == sum_len(res_l))
    return res_fa, res_l

def gen_unique_labels(d: dict):
    first_label = 1
    l = dict()
    for ch in d.keys():
        l[ch] = list(range(first_label, first_label + len(d[ch])))
        first_label += len(d[ch])
    return l