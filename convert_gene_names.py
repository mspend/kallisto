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
t2g_df = pd.read_table(transcripts_to_genes, header=None,
    names=[
        "target_id",
        "gene_id",
        "gene_name",
        "transcript_name",
        "chr",
        "start",
        "end",
        "strand",
    ],)
# print(t2g_df.head())

# Add a column to the counts_df on the far right containing the gene name from t2g_df
counts_merged= counts_df.merge(
    t2g_df[["target_id", "gene_name"]],
    on="target_id",
    how="left"
)

print(counts_merged)

# Quality control: check what percentage of the transcripts are missing a gene name
print(f'The percent of transcripts whose gene name is N/A is {counts_merged["gene_name"].isna().sum()/len(counts_df)*100}')

# Create a Series containing only the gene name and the TPM
gene_tpm = counts_merged.groupby("gene_name")["tpm"].sum().reset_index()
print(gene_tpm)

# Save as a CSV
# Will save to the present working directory
gene_tpm.to_csv('bulk_RNA_seq_counts_tpm.csv')
