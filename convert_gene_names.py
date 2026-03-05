# After quantifying reads, with kallisto, the output of the abundance.tsv file looks like this:
# target_id	length	eff_length	est_counts	tpm
# ENST00000513300.5	1924	1746.98	102.328	11129.2
# Turn ENSMBL Transcript Identifiers into human readable gene names

# Import necessary packages
import pandas as pd

# Read in the RNA counts from Kallisto, stored in abundance.tsv
# The format of abundance.tsv should be:
# target_id	length	eff_length	est_counts	tpm
# ENST00000513300.5	1924	1746.98	102.328	11129.2
# ENST00000282507.7	2355	2177.98	1592.02	138884

counts = '/data/bulkRNA/IMG/kallisto/abundance.tsv'
counts_df = pd.read_table(counts, delimiter='\t', header=0)
# print(counts_df.head)

# print(counts_df.index)
# print(counts_df.columns)

# print(counts_df['target_id'])


# Read in the file containing the ENSMBL Transcript Identifiers and the gene names
# T2g = transcripts to genes
# ENSMUST00000132100.2    ENSMUSG00000086053.2    Gm15178 Gm15178-201     1       75368775        75373007        -
# ENSMUST00000185910.2    ENSMUSG00000100764.2    Gm29155 Gm29155-201     1       43782744        43783012        -
transcripts_to_genes = '/data/bulkRNA/IMG/kallisto/pre_built_index_patcher_lab_t2g.txt'
t2g_df = pd.read_table(transcripts_to_genes)
# print(t2g_df.head())

num_lines = 0

for transcript in counts_df['target_id']:
    if transcript.startswith('ENSMUST'):
        num_lines += 1

print(num_lines)