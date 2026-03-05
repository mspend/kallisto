# After quantifying reads, with kallisto, the output of the abundance.tsv file looks like this:
# target_id	length	eff_length	est_counts	tpm
# ENST00000513300.5	1924	1746.98	102.328	11129.2
# Turn ENSMBL Transcript Identifiers into human readable gene names

input_path = Path('/data/bulkRNA/IMG/kallisto')

with open('abundance.tsv') as file:
    