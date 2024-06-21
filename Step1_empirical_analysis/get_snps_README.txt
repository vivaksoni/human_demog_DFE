#1. GET SAMPLES AND SAMPLE SIZES FOR RELEVANT POPULATIONS
bcftools query -l ALL.chr1.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz > samples.txt
join -j 1 <(sort -k1 samples.txt) <(sort -k1 igsr_samples.tsv) | awk -F ' ' '{print $1,$4,$6}' OFS='\t' > vcf_samples.txt

#2. DOWNLOAD GRCh37 VCF FILES
for i in $(seq 1 22); do wget -t 20 http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr"$i".phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz; done
for i in $(seq 1 22); do wget -t 20 http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr"$i".phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi; done

#3. GET 1KG STRICT MASK BEDFILE
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/accessible_genome_masks/20140520.strict_mask.autosomes.bed

#4. MASK PHASTCONS SITES FROM ACCESSIBILITY MASK FILE
for i in $(seq 1 22); do bedtools subtract -a <(awk -v chrom="$i" '{if($1=="chr"chrom)print}' accessibility_mask/20140520.strict_mask.autosomes.bed) -b phastCons/hg19_chr"$i".bed.gz -sorted; done | bgzip > accessibility_mask/20140520.strict_mask.autosomes_phastCons_masked.bed.gz

#5. GET VCF FILES FOR REGIONS OF INTEREST
chroms=($(awk '{print substr($1, 4)}' GRCh37_functionalMasked_15kb.bed))
starts=($(awk '{print $2}' GRCh37_functionalMasked_15kb.bed))
ends=($(awk '{print $3}' GRCh37_functionalMasked_15kb.bed))
regions=($(awk '{print $8}' GRCh37_functionalMasked_15kb.bed))

for i in $(seq 0 ${#chroms[@]}); do bcftools view -v snps -r "${chroms[i]}":${starts[i]}"-${ends[i]}" ALL.chr"${chroms[i]}".phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz | bgzip > regions_only/region"${regions[i]}".vcf.gz; done
for i in *.vcf.gz; do tabix "$i"; done

#6. REMOVE MASKED SITES (BOTH PHASTCONS AND ACCESSIBILITY MASKED)
for i in $(seq 0 ${#chroms[@]}); do cat <(bcftools view -h region"${regions[i]}".vcf.gz) <(bedtools intersect -a region"${regions[i]}".vcf.gz -b <(tabix ../../accessibility_mask/20140520.strict_mask.autosomes_phastCons_masked.bed.gz chr"${chroms[i]}" | awk -v chrom="${chroms[i]}" '{print chrom,$2,$3}' OFS='\t') -u -sorted) | bgzip > masked/region"${regions[i]}".vcf.gz; done
for i in *.vcf.gz; do tabix "$i"; done

#7. GET SNP COUNTS (REPEAT FOR MASKED VCFS)
for i in $(seq 0 ${#chroms[@]}); do tabix region"${regions[i]}".vcf.gz "${chroms[i]}" | wc -l | awk -v region="${regions[i]}" '{print region,$1}' OFS='\t' >> snp_counts.txt; done

#8. CREATE SAMPLES FILES FOR POPULATIONS OF INTEREST
awk '{if($2=="LWK")print $1}' vcf_samples.txt > vcf_samples_LWK.txt
awk '{if($3=="EUR")print $1}' vcf_samples.txt > vcf_samples_EUR.txt
awk '{if($3=="EAS")print $1}' vcf_samples.txt > vcf_samples_EAS.txt
awk '{if($3=="SAS")print $1}' vcf_samples.txt > vcf_samples_SAS.txt
cat vcf_samples_LWK.txt vcf_samples_EUR.txt > vcf_samples_LWK_EUR.txt
cat vcf_samples_LWK.txt vcf_samples_EAS.txt > vcf_samples_LWK_EAS.txt
cat vcf_samples_LWK.txt vcf_samples_SAS.txt > vcf_samples_LWK_SAS.txt
cat vcf_samples_EUR.txt vcf_samples_EAS.txt > vcf_samples_EUR_EAS.txt
cat vcf_samples_EUR.txt vcf_samples_SAS.txt > vcf_samples_EUR_SAS.txt
cat vcf_samples_EAS.txt vcf_samples_SAS.txt > vcf_samples_EAS_SAS.txt

#9. SUBSET VCFS FOR ONLY POPULATIONS OF INTEREST (FOR BOTH MASKED AND UNMASKED REGIONS)
pops=(4pops EAS EUR LWK SAS)
for p in "${pops[@]}"; do for i in $(seq 0 ${#regions[@]}); do bcftools view -S ../vcf_samples_"$p".txt -c 1 -i 'TYPE="SNP"' region"${regions[i]}".vcf.gz | bgzip > "$p"/region"${regions[i]}".vcf.gz; done; done
for p in "${pops[@]}"; do for i in $(seq 0 ${#regions[@]}); do tabix "$p"/region"${regions[i]}".vcf.gz; done; done
pops=(LWK_EUR LWK_EAS LWK_SAS EUR_EAS EUR_SAS EAS_SAS)
for p in "${pops[@]}"; do for i in $(seq 0 ${#regions[@]}); do bcftools view -S ../../vcf_samples_"$p".txt -c 1 -i 'TYPE="SNP"' region"${regions[i]}".vcf.gz | bgzip > "$p"/region"${regions[i]}".vcf.gz; done; done
for p in "${pops[@]}"; do for i in $(seq 0 ${#regions[@]}); do tabix "$p"/region"${regions[i]}".vcf.gz; done; done

#10. RUN get_empirical_Fst.py TO OBTAIN Fst STATS

#11. RUN get_empirical_stats.py TO OBTAIN FINAL STATS FILE




