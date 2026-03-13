# After quantifying reads, with kallisto, the output of the abundance.tsv file looks like this:
# target_id	length	eff_length	est_counts	tpm
# ENST00000513300.5	1924	1746.98	102.328	11129.2
# Turn ENSMBL Transcript Identifiers into human readable gene names

# Import necessary packages
import pandas as pd
from pathlib import Path

data = Path('/data/bulkRNA/IMG/kallisto')

# Read in the RNA counts from Kallisto, stored in abundance.tsv
counts = data / 'rep_1' / 'abundance.tsv'
counts_df = pd.read_table(counts, delimiter='\t', header=0)


# Read in the file containing the ENSMBL Transcript Identifiers and the gene names
# T2g = transcripts to genes
transcripts_to_genes = data / 'pre_built_index_patcher_lab_t2g.txt'
t2g_df = pd.read_table(transcripts_to_genes, header=None, names=["target_id", "gene_id", "gene_name", "transcript_name", "chr", "start", "end", "strand"])

# Add a column to the counts_df on the far right containing the gene name from t2g_df
counts_merged= counts_df.merge(t2g_df[["target_id", "gene_name"]], on="target_id", how="left")

# Quality control: check what percentage of the transcripts are missing a gene name
print(f'The percent of transcripts whose gene name is N/A is {counts_merged["gene_name"].isna().sum()/len(counts_df)*100}')

# Reorder the columns
reordered = counts_merged[["gene_name", "target_id", "tpm", "est_counts", "length", "eff_length"]]

# Create a Series containing only the gene name and the TPM
# For some genes, thare are multiple target_ids meaning multiple transcripts that map to the same gene. 
# These could be isoforms, splice variations, etc
# I summed the counts for all transcripts mapping to the same gene. 
summary = counts_merged.groupby("gene_name")["tpm"].sum().reset_index()

# Save as a CSV
# Will save to the present working directory
reordered.to_csv(data / 'rep_1' /'bulk_RNA_seq_counts.csv')

summary.to_csv(data / 'rep_1' /'bulk_RNA_seq_counts_summary.csv')
