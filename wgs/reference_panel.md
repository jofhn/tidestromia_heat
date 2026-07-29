# Reference panel

WGS of 48 individuals, of which 38 were used in the reference panel.
SRA  will be released upon publication.
Vendor: Novogene
Platform: NovaSeq X Plus (150 bp paired-end)
Run ID: SUS20241101029-JF20241029

Analysis steps from raw reads to a filtered SNP set and population substructure
analysis by PCA and F_ST.

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
Alignment rates were assessed with `samtools flagstat`.

```bash
RG="@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tLB:${SAMPLE}.lib1\tPL:ILLUMINA\tPU:lane1"

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

### Interval scattering

Chromosome intervals were split into 256 scattered interval lists for parallel
joint genotyping.

```bash
gatk SplitIntervals \
  -R "${REF}" \
  -L "${INTERVALS_CHR}" \
  --scatter-count 256 \
  -O "${SCATTER_DIR}"
```

### Per-sample calling

HaplotypeCaller in GVCF mode, run per sample against the chromosome intervals.

```bash
gatk --java-options "-Xmx28g" HaplotypeCaller \
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
gatk --java-options "-Xmx70g -Djava.io.tmpdir=${TMPDIR}" GenomicsDBImport \
  --sample-name-map "${SAMPLE_MAP}" \
  --genomicsdb-workspace-path "${WORKSPACE}" \
  -L "${CHR}" \
  --reader-threads 4 \
  --batch-size 50
```

### Joint genotyping

Run across the 256 scattered intervals, one task per interval.

```bash
gatk --java-options "-Xmx30g" GenotypeGVCFs \
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

The PASS SNP set is the input to sections 6-9 and to the survival panel
imputation.

## 6. PLINK conversion

Filter-passing SNPs were converted to PLINK 2 format, restricted to the DV,
Perimeter and Outside samples used in the reference panel (N = 38).

```bash
plink2 \
  --vcf "${SNPS_PASS}" \
  --keep "${ID_LIST}" \
  --max-alleles 2 \
  --snps-only just-acgt \
  --geno 0.05 \
  --mind 0.20 \
  --maf 0.05 \
  --set-all-var-ids '@:#:$r:$a' \
  --rm-dup force-first \
  --make-pgen \
  --threads "${THREADS}" \
  --out "${PFILE}"
``` 

Allele frequencies are not written automatically below 50 samples, so they were
generated separately and passed to the PCA with `--read-freq`.

```bash
plink2 \
  --pfile "${PFILE}" \
  --freq \
  --out "${PFILE}"
```

## 7. Differentiation

Pairwise Weir & Cockerham F_ST between collection regions, computed on the full
SNP set before any thinning.

```bash
plink2 \
  --pfile "${PFILE}" \
  --pheno "${POP_FILE}" \
  --fst pop method=wc report-variants \
  --out "${FST_OUT}"
```

## 8. Linkage disequilibrium decay

Unphased r2 between SNP pairs within a 1 Mb window, computed per chromosome on
the full SNP set. Run separately for Chr1-Chr11.

```bash
plink2 \
  --pfile "${PFILE}" \
  --chr "${CHR}" \
  --r2-unphased yes-really zs cols=chrom,pos,id \
  --ld-window-kb 1000 \
  --ld-window 999999 \
  --ld-window-r2 0 \
  --threads "${THREADS}" \
  --out "${LD_OUT}_${CHR}"
```

Pairs were binned by physical distance, with pair counts and summed r2 recorded
per bin, using the following edges (bp):

```
0 50 100 200 500 1000 2000 5000 10000 20000 50000 100000 150000 200000
300000 500000 750000 1000000 1500000 2000000 3000000 5000000
```

## 9. Population structure

### Distance thinning

SNPs were thinned by physical distance to one variant per 150 kb. LD-based
pruning was not used given the small sample size.

```bash
plink2 \
  --pfile "${PFILE}" \
  --bp-space 150000 \
  --write-snplist \
  --threads "${THREADS}" \
  --out "${THINNED}"
```

### Principal component analysis

Principal component analysis on the thinned set, using the allele frequencies
computed in section 6.

```bash
plink2 \
  --pfile "${PFILE}" \
  --extract "${THINNED}.snplist" \
  --read-freq "${FREQ}" \
  --pca 20 \
  --threads "${THREADS}" \
  --out "${PCA_OUT}"
```
