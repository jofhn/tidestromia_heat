# RNA-seq (heat/cold series)

RNA-seq of 33 libraries.
Raw reads and count matrices will be released at NCBI GEO upon publication.
Vendor: Novogene
Platform: NovaSeq X Plus (PE150)
Run ID: SUS20251122028-20250929hotcold

## 1. Read QC and trimming

FastQC and MultiQC were run on the raw and trimmed reads. Trimming used paired-end
adapter auto-detection, poly-G and poly-X trimming, and a 3' quality trim cutting from
the tail while the mean quality in the window fell below 20. Run per sample.

```bash
fastp \
  -i "${R1}" -I "${R2}" \
  -o "${OUT_R1}" -O "${OUT_R2}" \
  --detect_adapter_for_pe \
  --trim_poly_g --trim_poly_x \
  --cut_tail --cut_tail_mean_quality 20 \
  --thread "${THREADS}" \
  --json "${OUTDIR}/${SAMPLE}.fastp.json" \
  --html "${OUTDIR}/${SAMPLE}.fastp.html"
```

## 2. Genome index

The AmA10 PASA annotation was converted from GFF3 to GTF for STAR, then a STAR index
built from the chromosome-level assembly and that GTF. 

```bash
gffread "${GFF3}" -T -g "${FA}" -o "${GTF}"
```

```bash
STAR \
  --runThreadN "${THREADS}" \
  --runMode genomeGenerate \
  --genomeDir "${STAR_INDEX}" \
  --genomeFastaFiles "${FA}" \
  --sjdbGTFfile "${GTF}" \
  --sjdbOverhang 297
```

## 3. Alignment

Trimmed reads were aligned to the AmA10 chromosome-level reference assembly
(Prado et al. 2025) with STAR v2.7.11b, run per sample, writing a coordinate-sorted
BAM. Per-gene read counts were generated during alignment, giving a four-column
`ReadsPerGene.out.tab` per sample. BAMs were indexed with `samtools index` and
assessed with `samtools flagstat` and `samtools stats`.

```bash
STAR \
  --genomeDir "${STAR_INDEX}" \
  --readFilesIn "${R1}" "${R2}" \
  --readFilesCommand zcat \
  --runThreadN "${THREADS}" \
  --outFileNamePrefix "${SAMPLE_DIR}/${SAMPLE}." \
  --outSAMtype BAM SortedByCoordinate \
  --limitBAMsortRAM 16000000000 \
  --quantMode GeneCounts
```

## 4. Strandedness

Library strandedness was checked by summing the unstranded, forward and reverse count
columns of `ReadsPerGene.out.tab` across genes. The unstranded column was used in section 5.

```bash
awk 'BEGIN{u=0;f=0;r=0} $1!~/^N_/ {u+=$2; f+=$3; r+=$4} \
  END{printf("unstranded=%d forward=%d reverse=%d\n",u,f,r)}' "${TAB}"
```

## 5. Count matrix

Per-sample `ReadsPerGene.out.tab` files were assembled into a gene-by-sample matrix.
STAR's leading `N_*` summary rows were dropped and the unstranded count column taken.
Gene order was checked to be identical across samples before assembly.

Two libraries, JF11 and JF18, were found to be mislabelled during downstream QC (methods). 
The corrected assignments are JF11 = 55_hot rep3 and JF18 = 55_cold rep2. 
Column names in this matrix are the original library names, so the corrected labels are 
carried in the sample metadata rather than the counts file.

```python
files = sorted(glob.glob(f"{STARBASE}/JF*/JF*.ReadsPerGene.out.tab"))

for f, s in zip(files, samples):
    with open(f) as fh:
        for line in fh:
            gene, unstranded, fwd, rev = line.rstrip("\n").split("\t")
            if gene.startswith("N_"):
                continue
            g.append(gene)
            c.append(unstranded)
```