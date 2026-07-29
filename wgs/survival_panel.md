# Survival panel

WGS of 225 libraries.
SRA  will be released upon publication.
Vendor: Novogene
Platform: Novasex X 25B (PE150)
Run ID: 25363Srh

Analysis steps from raw reads to a filtered SNP set, imputation against the
reference panel, and association testing.

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

## 7. Intersect with reference

## 8. Reference panel phasing

## 9. Imputation

## 10. Post-imputation processing

## 11. Population structure

## 12. Survival phenotype derivation

## 13. Association testing

## 14. Clumping and locus boundaries
EOF