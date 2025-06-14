trnl_ref=/fast2/def-ilafores/envBarcodeMiner/trnL_CGAAATCGGTAGACGCTACG/hits.lineage.noacc.fa

cat $trnl_ref | sed -E 's/^>(([^;]*;){5}[^;]*);.*/>\1;/' > /fast2/def-ilafores/envBarcodeMiner/trnL_CGAAATCGGTAGACGCTACG/hits.lineage.Genus.fa

