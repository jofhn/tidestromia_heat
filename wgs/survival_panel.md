# Survival panel

WGS of 225 libraries.
SRA will be released upon publication.
Vendor: Novogene
Platform: Novaseq X 25B (PE150)
Run ID: 25363Srh

Sections 1-5 were run once on all samples. From section 6 the analysis splits into
two target sets, all populations together (all_pops) and DV/Perimeter (DVP). Each was
imputed separately against the same species-range phased reference. Association testing 
was run with DV/Perimeter only. Commands from section 6 onward were run once per target unless noted.

## 1. Read QC and trimming

Raw read files were verified against vendor MD5 checksums. Reads were concatenated
across lanes per sample. FastQC and MultiQC were run on the raw, concatenated, and
trimmed reads.

### Trimming

```bash
fastp \
  -i "${R1}" -I "${R2}" \
  -o "${OUT_R1}" -O "${OUT_R2}" \
  --detect_adapter_for_pe \
  --trim_poly_g --trim_poly_x \
  --thread "${THREADS}" \
  --json "${OUTDIR}/${BASE}.fastp.json" \
  --html "${OUTDIR}/${BASE}.fastp.html"
```

## 2. Alignment

Trimmed reads were aligned to the AmA10 chromosome-level reference assembly
(Prado et al. 2025) with bwa-mem2, read groups assigned during alignment, and
output piped directly to a coordinate-sorted BAM. Run per sample.
Alignment was assessed with `samtools flagstat` and `samtools stats`.

```bash
RG="@RG\tID:${SAMPLE}.${LANE}\tSM:${SAMPLE}\tLB:${SAMPLE}.lib1\tPL:ILLUMINA\tPU:${LANE}"

bwa-mem2 mem -t "${THREADS}" -R "${RG}" "${REF}" "${R1}" "${R2}" \
  | samtools sort -@ "${THREADS}" -m 2G -o "${BAM}" -
```

## 3. Duplicate marking

Duplicates were marked, not removed, with GATK MarkDuplicates.
Deduplicated BAMs were assessed with `samtools flagstat` and mosdepth v0.3.6.

```bash
gatk --java-options "-Xmx28g" MarkDuplicates \
  -I "${IN_BAM}" \
  -O "${OUT_BAM}" \
  -M "${METRICS}" \
  --TMP_DIR "${TMP_DIR}" \
  --VALIDATION_STRINGENCY SILENT \
  --CREATE_INDEX true
```

## 4. Variant calling

Chromosome intervals were split into 256 scattered interval lists for parallel
joint genotyping (see reference panel, section 4).

### Per-sample calling

HaplotypeCaller in GVCF mode, run per sample against the chromosome intervals.

```bash
gatk --java-options "-Xmx28g -Djava.io.tmpdir=${TMPDIR}" HaplotypeCaller \
  -R "${REF}" \
  -I "${BAM}" \
  -O "${GVCF}" \
  -ERC GVCF \
  --native-pair-hmm-threads "${THREADS}" \
  -L "${INTERVALS_CHR}"
```

### Consolidation

GenomicsDBImport with a sample map of all per-sample gVCFs, one workspace per
chromosome.

```bash
gatk --java-options "-Xmx110g -Djava.io.tmpdir=${TMPDIR}" GenomicsDBImport \
  --sample-name-map "${SAMPLE_MAP}" \
  --genomicsdb-workspace-path "${WORKSPACE}" \
  -L "${CHR}" \
  --reader-threads 4 \
  --batch-size 50
```

### Joint genotyping

Run across the 256 scattered intervals, one task per interval.

```bash
gatk --java-options "-Xmx30g -Djava.io.tmpdir=${TMPDIR}" GenotypeGVCFs \
  -R "${REF}" \
  -V "gendb://${WORKSPACE}" \
  -L "${INTERVAL_FILE}" \
  -O "${OUT_VCF}"
```

### Gathering

Scattered VCFs merged, then split into per-chromosome VCFs.

```bash
gatk GatherVcfs -I "${VCF_LIST}" -O "${MERGED_VCF}"

bcftools view -r "${CHR}" -Oz -o "${CHR_VCF}" "${MERGED_VCF}"
```

## 5. Hard filtering

Run on the gathered chromosome VCF: SNPs were extracted, labelled against GATK
Best Practices thresholds, then failing records removed. VariantRecalibrator and
ApplyVQSR were not used, as no validated variant set exists for this species.

```bash
gatk SelectVariants -R "${REF}" -V "${MERGED_VCF}" \
  --select-type-to-include SNP \
  -O "${SNPS}"
```

```bash
gatk VariantFiltration -R "${REF}" -V "${SNPS}" \
  -O "${SNPS_LABELLED}" \
  --filter-name "QD2"                --filter-expression "QD < 2.0" \
  --filter-name "FS60"               --filter-expression "FS > 60.0" \
  --filter-name "SOR3"               --filter-expression "SOR > 3.0" \
  --filter-name "MQ40"               --filter-expression "MQ < 40.0" \
  --filter-name "MQRankSum-12.5"     --filter-expression "MQRankSum < -12.5" \
  --filter-name "ReadPosRankSum-8"   --filter-expression "ReadPosRankSum < -8.0"
```

```bash
gatk SelectVariants -R "${REF}" -V "${SNPS_LABELLED}" \
  --exclude-filtered \
  -O "${SNPS_PASS}"
```

The PASS SNP set is the input to section 6.

## 6. Subset VCFs

The survival panel PASS SNP set was subset by sample ID list to the 
DV/Perimeter individuals. The all_pops target used the unsubsetted PASS SNP set. 

```bash
bcftools view --threads "${THREADS}" -S "${IDS}" -Ob -o "${OUT_BCF}" "${IN_VCF}"
bcftools index -f "${OUT_BCF}"
```

## 7. Intersect with reference

The survival panel and the reference panel were each reduced to the SNPs called
in both PASS SNP sets. Shared SNPs were identified with `bcftools isec` and the
resulting list used to subset each side. The reduced reference is the input to
phasing in section 8; the reduced survival panel is the imputation target in
section 9. Shared BCFs were converted to bgzipped VCF for input to BEAGLE.

```bash
bcftools isec --threads "${THREADS}" -p "${ISEC_DIR}" -n=2 -w1 \
  "${REF_BCF}" \
  "${TARGET_BCF}"

cp "${ISEC_DIR}/sites.txt" "${SHARED_SITES}"
```

```bash
bcftools view --threads "${THREADS}" -R "${SHARED_SITES}" -Ob -o "${REF_SHARED}" \
  "${REF_BCF}"
bcftools index -f "${REF_SHARED}"

bcftools view --threads "${THREADS}" -R "${SHARED_SITES}" -Ob -o "${TARGET_SHARED}" \
  "${TARGET_BCF}"
bcftools index -f "${TARGET_SHARED}"
```

## 8. Reference panel phasing

The reduced reference panel from section 7 was split per chromosome and 
phased with BEAGLE 5.4, run for Chr1-Chr11. No `ref=` argument is given, 
so BEAGLE phases the input rather than imputing against a panel. 
Phasing used a fixed effective population size of 10,000. The phased output
is the reference for both imputation runs in section 9.

```bash
bcftools view --threads "${THREADS}" -r "${CHR}" -Oz -o "${GT}" "${REF_SHARED}"
bcftools index -f "${GT}"
```

```bash
java -Xmx120g \
  -Djava.io.tmpdir="${TMPDIR}" \
  -jar "${EBROOTBEAGLE}/beagle.jar" \
  gt="${GT}" \
  out="${OUT_PREFIX}" \
  chrom="${CHR}" \
  nthreads="${THREADS}" \
  window=20.0 \
  overlap=2.0 \
  iterations=10 \
  ne=10000
```

## 9. Imputation

### Effective population size (Ne)

BEAGLE's Ne was set per target from nucleotide diversity, calculated in 10 kb windows
on the PASS SNP set of each target, then Ne = mean(pi) / (4 mu) with a mutation rate of
7e-9 per site per generation. This gave Ne = 82,000 for all_pops and 85,000 for DV/Perimeter.

```bash
vcftools --gzvcf "${VCF}" \
  --window-pi 10000 \
  --out "${PI_PREFIX}"
```

### Imputation

Both targets were imputed against the same phased species-range reference panel from
section 8, run per chromosome. All BEAGLE parameters are identical between the 
two runs except Ne. DR2 filtering is applied in section 10.

```bash
bcftools view --threads "${THREADS}" -r "${CHR}" -Oz -o "${GT}" "${TARGET_SHARED_VCF}"
bcftools index -f "${GT}"
```

```bash
java -Xmx28g \
  -Djava.io.tmpdir="${TMPDIR}" \
  -jar "${EBROOTBEAGLE}/beagle.jar" \
  ref="${REF_PHASED}" \
  gt="${GT}" \
  out="${OUT_PREFIX}" \
  chrom="${CHR}" \
  nthreads="${THREADS}" \
  window=40.0 \
  overlap=4.0 \
  iterations=20 \
  ne="${NE}"
```

## 10. Post-imputation processing

Imputed genotypes were filtered on BEAGLE's DR2 > 0.9, applied per chromosome. 
Filtered chromosomes were concatenated to a whole-genome VCF and restricted to biallelic
SNPs. Dosages were then converted to a PLINK 2 pgen and variants sorted. Run once per target.

```bash
bcftools filter --threads "${THREADS}" \
  -i "DR2>0.9" \
  -Oz -o "${FILT_VCF}" \
  "${IMPUTED_VCF}"
bcftools index -f "${FILT_VCF}"
```

```bash
bcftools concat --threads "${THREADS}" -Oz -o "${CONCAT_VCF}" "${CHR_VCFS[@]}"
bcftools index -f "${CONCAT_VCF}"

bcftools view --threads "${THREADS}" -m2 -M2 -v snps \
  -Oz -o "${FILTERED_VCF}" \
  "${CONCAT_VCF}"
bcftools index -f "${FILTERED_VCF}"
```

```bash
plink2 \
  --vcf "${FILTERED_VCF}" dosage=DS \
  --keep "${IIDS}" \
  --max-alleles 2 \
  --snps-only just-acgt \
  --maf 0.05 \
  --mind 0.20 \
  --set-all-var-ids '@:#:$r:$a' \
  --rm-dup force-first \
  --make-pgen \
  --threads "${THREADS}" \
  --out "${IMPUTED_PGEN}"

plink2 \
  --pfile "${IMPUTED_PGEN}" \
  --sort-vars \
  --make-pgen \
  --threads "${THREADS}" \
  --out "${IMPUTED_PGEN}_sorted"
```

## 11. Population structure

### Weir and Cockerham's estimator of Fst (1984)

Pairwise Fst was calculated between DV, Perimeter and Outside on the all_pops imputed SNP set. 

```bash
plink2 \
  --pfile "${PGEN_SORTED}" \
  --pheno "${POP_FILE}" \
  --pheno-name POP \
  --fst POP \
  --threads "${THREADS}" \
  --out "${OUT_PREFIX}"
```

### LD pruning

Each sorted pgen was LD-pruned in windows of 150 kb (all_pops) or 100 kb
(DV/Perimeter), stepping one variant at a time and dropping variants with r2 > 0.2
against a retained variant. The resulting `.prune.in` list is the input to PCA.

```bash
plink2 \
  --pfile "${PGEN_SORTED}" \
  --indep-pairwise "${WINDOW}" 1 0.2 \
  --threads "${THREADS}" \
  --out "${OUT_PREFIX}.prune"
```

### Principal component analysis

PCA was run on the LD-pruned SNP set of each sorted pgen.

```bash
plink2 \
  --pfile "${PGEN_SORTED}" \
  --extract "${PRUNE_IN}" \
  --pca 50 \
  --threads "${THREADS}" \
  --out "${OUT_PREFIX}"
```

## 12. Survival phenotype derivation

Relative survival was derived by Cox proportional-hazards regression on days to death, fitted in R (methods).

## 13. Association testing

Association testing was run with a linear model in PLINK 2, one chromosome per task,
on the sorted imputed DV/Perimeter pgen. Phenotype was relative survival from section 12 (time_to_four_rel). 
Phenotype and covariates were read from the same table. Covariates were experimental tray, coded
as dummy variables with one level dropped as the reference, plus principal components 1-25 
from the LD-pruned PCA for DV/Perimeter in section 11.

```bash
plink2 \
  --pfile "${PGEN_SORTED}" \
  --chr "${CHR}" \
  --pheno "${PHENO}" \
  --pheno-name time_to_four_rel \
  --covar "${PHENO}" \
  --covar-name "${COVARS[@]}" \
  --glm zs omit-ref pheno-ids hide-covar \
  --threads "${THREADS}" \
  --memory 32000 \
  --out "${OUT_PREFIX}_chr${CHR}"
```

## 14. Clumping and locus boundaries

SNPs were grouped into approximately independent loci by LD clumping in PLINK 2, 
with the sorted imputed pgen as the LD reference and the glm results from section 13 as input. 
The index threshold was Bonferroni-corrected at 1.16e-7, 
recruiting variants within 300 kb at r2 > 0.1 and p < 1e-2 from unphased genotypes.

```bash
plink2 \
  --pfile "${PGEN_SORTED}" \
  --clump "${GLM_RESULTS}" \
  --clump-snp-field ID \
  --clump-field P \
  --clump-unphased \
  --clump-p1 1.16e-7 \
  --clump-p2 1e-2 \
  --clump-r2 0.1 \
  --clump-kb 300 \
  --threads "${THREADS}" \
  --out "${OUT_PREFIX}"
```