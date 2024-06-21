#1. CONVERT GFF INTO BED
sortBed -i GRCh37_latest_genomic.gff | gff2bed > GRCh37_latest_genomic.bed
sortBed -i homo_sapiens.GRCh37.Regulatory_Build.regulatory_features.20201218.gff | gff2bed > GRCh37_regulatory_features.bed

#2. RENAME CHROMOSOMES FROM REFSEQ TO NUMBERS FOR GENE GFF FILE
refseq=(NC_000001.10 NC_000002.11 NC_000003.11 NC_000004.11 NC_000005.9 NC_000006.11 NC_000007.13 NC_000008.10 NC_000009.11 NC_000010.10 NC_000011.9 NC_000012.11 NC_000013.10 NC_000014.8 NC_000015.9 NC_000016.9 NC_000017.10 NC_000018.9 NC_000019.9 NC_000020.10 NC_000021.8 NC_000022.10 )

numseq=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22)

for i in $(seq 0 21); do awk -v r="${refseq["$i"]}" -v c="${numseq["$i"]}" '{if($1==r)print c,$2,$3,$4,$5,$6,$7,$8,$9,$10}' OFS='\t' GRCh37_latest_genomic.bed; done | grep -v -e $'RefSeq\tregion' > GRCh37_latest_genomic_chromsRenamed.bed

#3. ADD 10KB FLANKS TO BED FILES
bedtools flank -i GRCh37_regulatory_features.bed -g GRCh37_chrom_lengths.txt -b 10000 | sortBed > GRCh37_regulatory_features_10kb_flanks.bed
bedtools flank -i GRCh37_latest_genomic_chromsRenamed.bed -g GRCh37_chrom_lengths.txt -b 10000 | sortBed > GRCh37_latest_genomic_chromsRenamed_10kb_flanks.bed

#4. MERGE REGIONS
bedtools merge -i GRCh37_regulatory_features_10kb_flanks.bed | sortBed -g GRCh37_chrom_lengths.txt > GRCh37_regulatory_features_10kb_flanks_merged_elements.bed
bedtools merge -i GRCh37_latest_genomic_chromsRenamed_10kb_flanks.bed | sortBed -g GRCh37_chrom_lengths.txt > GRCh37_latest_genomic_chromsRenamed_10kb_flanks_merged_elements.bed
  
#5. COMBINE BED FILES
cat GRCh37_latest_genomic_chromsRenamed_10kb_flanks_merged_elements.bed GRCh37_regulatory_features_10kb_flanks_merged_elements.bed | sortBed -g GRCh37_chrom_lengths.txt | bedtools merge > GRCh37_genomic_regulatory_10kb_flanks.bed

#6. MASK FUNCTIONAL REGIONS
bedtools subtract -a GRCh37_chrom_lengths.bed -b GRCh37_genomic_regulatory_10kb_flanks.bed | awk '{print $1,$2,$3,$3-$2}' OFS='\t' > GRCh37_functionalMasked.bed

#7. CONVERT PHASTCONS BIGWIG TO BED, INDEX, AND SPLIT BY CHROMOSOME, REMOVE PHASTCONS SCORES OF 0
bigWigToBedGraph hg19.100way.phastCons.bw hg19.100way.phastCons.bed
tabix hg19.100way.phastCons.bed.gz
for i in $(seq 1 22); do tabix hg19.100way.phastCons.bed.gz chr"$i" | bgzip > hg19.100way.phastCons_chr"$i".bed.gz; done
for i in $(seq 1 22); do tabix gff/hg19.100way.phastCons_chr"$i".bed.gz chr"$i" | awk '{if($4>0)print}' | bgzip > phastCons/hg19_chr"$i".bed.gz; done

#10. CREATE PHASTCONS BED FILES FOR EACH FUNCTIONAL REGION
for i in $(seq 0 ${#chroms[@]}); do tabix hg19_chr"${chroms[i]}".bed.gz chr"${chroms[i]}":"${starts[i]}"-"${ends[i]}" | awk -v start="${starts[i]}" '{print $1,$2-(start-2),$3-(start-2)}' OFS='\t' | awk '{if($2>0){print $1,$2,$3} else{print $1,1,$3}}' OFS='\t' | bgzip > functional_regions/chr"${chroms[i]}"_"${starts[i]}"_"${ends[i]}".bed.gz; done


#12. DOWNLOAD 1KG ACCESSIBILITY MASKS AND MASK USING PHASTCONS SCORES
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/accessible_genome_masks/20140520.strict_mask.autosomes.bed
for i in $(seq 1 22); do bedtools subtract -a <(awk -v chrom="$i" '{if($1=="chr"chrom)print}' accessibility_mask/20140520.strict_mask.autosomes.bed) -b phastCons/hg19_chr"$i".bed.gz -sorted; done | bgzip > accessibility_mask/20140520.strict_mask.autosomes_phastCons_masked.bed.gz

#13. ADD INFORMATION ABOUT HOW MANY SITES PASS STRICT MASK THRESHOLD (AND ARE NOT CONSTRAINED)
bedtools intersect -a gff/GRCh37_functionalMasked_500bp.bed -b accessibility_mask/20140520.strict_mask.autosomes_phastCons_masked.bed.gz -wo -sorted | bedtools sort | bedtools groupby -c 9 -o sum | bedtools intersect -a - -b gff/GRCh37_functionalMasked_500bp.bed -wo | awk '{print $1,$2,$3,$8,$4}' OFS='\t' > GRCh37_functionalMasked_500bp_1kG_accessibility.bed

#14. GET DECODE RECOMBINATION RATES FROM UCSC
wget -t 20 https://hgdownload.soe.ucsc.edu/gbdb/hg19/decode/SexAveraged.bw
bigWigToBed SexAveraged.bw SexAveraged.bedgraph

#15. INTERSECT NON-FUNCTIONAL REGIONS WITH REC RATES
bedtools intersect -a GRCh37_functionalMasked_500bp_1kG_accessibility.bed -b recombinationRates/SexAveraged.bedgraph -wo | awk '{print $1,$2,$3,$4,$5,$9,$10}' OFS='\t' > GRCh37_functionalMasked_500bp_1kG_accessibility_recRates.bed

#16. GET MUTATION RATES FROM MUTATION RATE BROWSER
wget https://download.molgeniscloud.org/downloads/gonl_public/mutation_rate_map/release2/local_mutation_rate.bias_corrected.SEXAVG.bed
wget -t 20 http://mutation.sph.umich.edu/hg19/*.bw
for i in *.bw; do bigWigToBed "$i" "$i".bed; done
for i in *.bed; do bgzip "$i"; done
for i in *.bed.gz; do tabix "$i"; done

#17 INTERSECT NON-FUNCTIONAL REGIONS WITH MUTATION RATES
bedtools intersect -a GRCh37_functionalMasked_500bp_1kG_accessibility.bed -b mutationRates/local_mutation_rate_mean.bed -wo | awk '{print $1,$2,$3,$4,$5,$9,$10}' OFS='\t' > GRCh37_functionalMa
sked_500bp_1kG_accessibility_mutRates.bed

#18. RUN human_demog_DFE_parse_nonFunc_regions.ipynb TO CALCULATE UNMASKED SEQUENCE LENGTHS AND ADD MUTATION AND RECOMBINATION RATES

#19. GET 20KB WINDOWS
bedtools makewindows -b GRCh37_chrom_lengths.bed -w 20000 | awk '{print "chr"$1,$2-1,$3-1}' OFS='\t' > GRCh37_chrom_lengths_20kb_wins.bed
files=(AT_CG AT_GC AT_TA GC_CG GC_TA GC_AT)
for f in "${files[@]}"; do for i in $(seq 1 22); do tabix
 "$f".bed.gz chr"$i" | bedtools intersect -a - -b ../gff/GRCh37_chrom_lengths_20kb_wins.bed -wo | awk '{print $5,$6,$7,$4}' OFS='\t' | sort -k1,1 -k2,2n | bedtools groupby -c 4 -o mean >> "$f"_20kb_wins.bed; done; done
cat *20kb_wins.bed | bedtools sort | bedtools groupby -c 4 -o mean > 20kb_wins_meanRates.bed


