### Sequence names are problematic, should be ASV1 or whole sequence? latter creates name duplication for some reason

barcode='ITS'
out_path="data/${barcode}"
# Alignment
module load mafft
mkdir -p "${out_path}/5_tree"
mafft --thread -72 "${out_path}/4_taxonomy/asv.fa" > "${out_path}/5_tree/asv_alignment.fa"

# Mask alignments 
trimal -in "${out_path}/5_tree/asv_alignment.fa" -out "${out_path}/5_tree/asv_alignment_masked.fa" -automated1


# Compute tree (raxml)
module load StdEnv/2023  gcc/12.3  openmpi/4.1.5 raxml-ng/1.2.0
raxml-ng --msa "${out_path}/5_tree/asv_alignment.fa" --model GTR+G --prefix "${out_path}/5_tree/tree" --threads 48 # Don't use more, oversubscription error

# Or iq-tree
ml iq-tree
iqtree2 -s "${out_path}/5_tree/asv_alignment_masked.fa" -m MFP -B 1000 --prefix "${out_path}/5_tree/tree_masked" -T AUTO
