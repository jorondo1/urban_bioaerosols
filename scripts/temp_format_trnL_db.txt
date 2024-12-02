trnl_ref=/home/def-ilafores/analysis/20240125_trnL_db/trnL_hits.lineage.filtered.fa

cat $trnl_ref | sed -E 's/^>(([^;]*;){5}[^;]*);.*/>\1;/' > /home/def-ilafores/analysis/urbanBio/data/trnL/trnL_hits_reference2024.fa

