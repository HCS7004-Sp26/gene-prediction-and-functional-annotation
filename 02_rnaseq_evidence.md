# Module 2 — RNA-seq Evidence Preparation

## Downloading, Trimming, and Validating the SRP166999 Dataset

---

## Background

High-quality RNA-seq evidence is the single most effective way to improve
gene model accuracy in a eukaryotic genome annotation. Transcript boundaries
resolved from direct RNA sequencing are far more reliable than ab initio
predictions alone, particularly for:

* Defining precise exon–intron boundaries
* Distinguishing real genes from spurious ORFs predicted in intergenic regions
* Capturing alternatively spliced isoforms
* Providing expression evidence for lowly conserved, lineage-specific genes

The Wang et al. (2020) RNA-seq dataset (GEO accession **GSE121872**, SRA study
**SRP166999**, BioProject **PRJNA498715**) covers *P. fructicola* grown under
three biologically distinct conditions. This dataset has **no replicates** —
one RNA-seq library was generated per condition — which is common in annotation
support datasets but has important implications for downstream differential
expression analysis (discussed in the biological questions below).

> **Connection to Module 6 of the previous tutorial:** In that module you
> processed a single RNA-seq run through alignment and quantification.
> Here, you download and trim the complete 3-run dataset that will serve
> as transcript evidence for gene prediction.

---

## Learning Objectives

* Download all 3 RNA-seq runs from SRP166999 using SRA Toolkit
* Perform adapter trimming with Trimmomatic
* Validate read quality with FastQC
* Merge trimmed reads into a format ready for transcript assembly (Module 4b)

---

## The SRP166999 Dataset

| Condition | Code | GSM | SRR | Approx. size |
| --- | --- | --- | --- | --- |
| Potato dextrose (PD), 5 days | PD-5d | GSM3449211 | SRR8115198 | ~1.4 Gb |
| Potato dextrose (PD), 15 days | PD-15d | GSM3449212 | SRR8115199 | ~1.3 Gb |
| PD + PEG 6000 (osmotic stress), 15 days | PD-PEG | GSM3449213 | SRR8115200 | ~1.3 Gb |

> **Important — study design note:** This dataset has **no biological
> replicates**. Each condition is represented by a single RNA-seq library.
> This is acceptable for providing transcript evidence to a gene predictor
> (which only needs to detect expressed exon boundaries), but it means
> the data **cannot be used for differential expression analysis** —
> statistical tests of DE require replicates to estimate within-condition
> variance. Keep this distinction in mind when interpreting annotation
> quality in Module 6.

> **Verify before running:** Always confirm run accessions in the SRA Run
> Selector before submitting download jobs. Go to
> <https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA498715> to
> download the current accession list.

---

## Step 1: Create the Accession List File

```bash
cat > ${ANNOT}/01_rnaseq/SRP166999_accessions.txt << 'EOF'
SRR8115198
SRR8115199
SRR8115200
EOF

echo "Accession list created:"
cat ${ANNOT}/01_rnaseq/SRP166999_accessions.txt
```

---

## Step 2: Download All 3 Runs via SLURM Job Array

```bash
cat > ${ANNOT}/scripts/02a_download_sra_array.sh << 'EOF'
#!/bin/bash
#SBATCH --account=PAS3260
#SBATCH --job-name=sra_download
#SBATCH --output=logs/sra_download_%A_%a.out
#SBATCH --error=logs/sra_download_%A_%a.err
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=16G
#SBATCH --array=1-3

set -euo pipefail

user_name=Jonathan

CONTAINERS=/fs/scratch/PAS3260/${user_name}/Annotation/containers
ANNOT=/fs/scratch/PAS3260/${user_name}/Annotation
OUTDIR=${ANNOT}/01_rnaseq/raw
mkdir -p ${OUTDIR}/tmp

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${ANNOT}/01_rnaseq/SRP166999_accessions.txt)
echo "Downloading: ${SAMPLE}  (task ${SLURM_ARRAY_TASK_ID}/3)"

apptainer exec \
  --bind ${OUTDIR}:/data \
  ${CONTAINERS}/sratools_3.2.1.sif \
  prefetch ${SAMPLE} \
    --output-directory /data \
    --max-size 50G

apptainer exec \
  --bind ${OUTDIR}:/data \
  ${CONTAINERS}/sratools_3.2.1.sif \
  fasterq-dump /data/${SAMPLE}/${SAMPLE}.sra \
    --split-files \
    --outdir /data \
    --temp /data/tmp \
    --threads 4

gzip ${OUTDIR}/${SAMPLE}_1.fastq
gzip ${OUTDIR}/${SAMPLE}_2.fastq

echo "${SAMPLE} complete."
ls -lh ${OUTDIR}/${SAMPLE}_*.fastq.gz
EOF

cd ${ANNOT}
sbatch ${ANNOT}/scripts/02a_download_sra_array.sh
squeue -u ${USER}
```

Once all 3 jobs complete, verify all files:

```bash
echo "=== Downloaded FASTQ files ==="
ls -lh ${ANNOT}/01_rnaseq/raw/SRR8115*.fastq.gz | wc -l
echo "files (expected: 6 — R1 + R2 for each of 3 runs)"

echo ""
echo "=== Read counts per run ==="
for ACC in $(cat ${ANNOT}/01_rnaseq/SRP166999_accessions.txt); do
  R1="${ANNOT}/01_rnaseq/raw/${ACC}_1.fastq.gz"
  if [[ -f "${R1}" ]]; then
    N=$(zcat ${R1} | wc -l | awk '{print $1/4}')
    echo "${ACC}: ${N} read pairs"
  else
    echo "${ACC}: MISSING"
  fi
done
```
---

## Step 3: Adapter Trimming — SLURM Array

```bash
cat > ${ANNOT}/scripts/02b_trimmomatic_array.sh << 'EOF'
#!/bin/bash
#SBATCH --account=PAS3260
#SBATCH --job-name=trimmomatic
#SBATCH --output=logs/trimmomatic_%A_%a.out
#SBATCH --error=logs/trimmomatic_%A_%a.err
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=16G
#SBATCH --array=1-3

set -euo pipefail

user_name=Jonathan

CONTAINERS=/fs/scratch/PAS3260/${user_name}/Annotation/containers
ANNOT=/fs/scratch/PAS3260/${user_name}/Annotation

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${ANNOT}/01_rnaseq/SRP166999_accessions.txt)
OUTDIR=${ANNOT}/01_rnaseq/trimmed
mkdir -p ${OUTDIR}

echo "Trimming: ${SAMPLE} (task ${SLURM_ARRAY_TASK_ID}/3)"

ADAPTER=/opt/conda/pkgs/trimmomatic-0.40-hdfd78af_0/share/trimmomatic-0.40-0/adapters/TruSeq3-PE-2.fa

apptainer exec \
  --bind ${ANNOT}:/data \
  ${CONTAINERS}/trimmomatic_0.40.sif \
  trimmomatic PE \
    -threads 8 \
    /data/01_rnaseq/raw/${SAMPLE}_1.fastq.gz \
    /data/01_rnaseq/raw/${SAMPLE}_2.fastq.gz \
    /data/01_rnaseq/trimmed/${SAMPLE}_1_paired.fastq.gz \
    /data/01_rnaseq/trimmed/${SAMPLE}_1_unpaired.fastq.gz \
    /data/01_rnaseq/trimmed/${SAMPLE}_2_paired.fastq.gz \
    /data/01_rnaseq/trimmed/${SAMPLE}_2_unpaired.fastq.gz \
    ILLUMINACLIP:${ADAPTER}:2:30:10:2:True \
    LEADING:3 \
    TRAILING:3 \
    SLIDINGWINDOW:4:20 \
    MINLEN:36

echo "${SAMPLE} trimming complete."
ls -lh ${OUTDIR}/${SAMPLE}_*_paired.fastq.gz
EOF

cd ${ANNOT}
sbatch ${ANNOT}/scripts/02b_trimmomatic_array.sh
squeue -u ${USER}
```

---

## Step 4: Quality Assessment with FastQC

```bash
cat > ${ANNOT}/scripts/02c_fastqc.sh << 'EOF'
#!/bin/bash
#SBATCH --account=PAS3260
#SBATCH --job-name=fastqc
#SBATCH --output=logs/fastqc_%j.out
#SBATCH --error=logs/fastqc_%j.err
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=16G

set -euo pipefail

user_name=Jonathan

CONTAINERS=/fs/scratch/PAS3260/${user_name}/Annotation/containers
ANNOT=/fs/scratch/PAS3260/${user_name}/Annotation

mkdir -p ${ANNOT}/01_rnaseq/qc

apptainer exec \
  --bind ${ANNOT}:/data \
  ${CONTAINERS}/fastqc_0.12.1.sif \
  bash -c "fastqc \
    /data/01_rnaseq/trimmed/*_paired.fastq.gz \
    --outdir /data/01_rnaseq/qc \
    --threads 8"

echo "FastQC complete."
ls ${ANNOT}/01_rnaseq/qc/*.html | wc -l
echo "HTML reports generated (expected: 6 — R1 + R2 for each of 3 runs)"
EOF

cd ${ANNOT}
sbatch ${ANNOT}/scripts/02c_fastqc.sh
```

---

## Step 5: Merge Trimmed Reads for Transcript Assembly

Modules 4b and 5 use a single merged R1 and R2 file. Merging all three
conditions ensures the transcript assembly captures genes expressed under
any of the profiled conditions:

```bash
cat > ${ANNOT}/scripts/02d_merge_reads.sh << 'EOF'
#!/bin/bash
#SBATCH --account=PAS3260
#SBATCH --job-name=merge_reads
#SBATCH --output=logs/merge_reads_%j.out
#SBATCH --error=logs/merge_reads_%j.err
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --mem=16G

set -euo pipefail

user_name=Jonathan

ANNOT=/fs/scratch/PAS3260/${user_name}/Annotation
TRIMMED=${ANNOT}/01_rnaseq/trimmed

echo "Merging R1 reads from all 3 samples..."
cat ${TRIMMED}/SRR8115198_1_paired.fastq.gz \
    ${TRIMMED}/SRR8115199_1_paired.fastq.gz \
    ${TRIMMED}/SRR8115200_1_paired.fastq.gz \
    > ${ANNOT}/01_rnaseq/all_samples_R1.fastq.gz

echo "Merging R2 reads from all 3 samples..."
cat ${TRIMMED}/SRR8115198_2_paired.fastq.gz \
    ${TRIMMED}/SRR8115199_2_paired.fastq.gz \
    ${TRIMMED}/SRR8115200_2_paired.fastq.gz \
    > ${ANNOT}/01_rnaseq/all_samples_R2.fastq.gz

echo "=== Merge complete ==="
ls -lh ${ANNOT}/01_rnaseq/all_samples_R*.fastq.gz

echo ""
echo "Total read pairs in merged R1:"
zcat ${ANNOT}/01_rnaseq/all_samples_R1.fastq.gz | wc -l | awk '{print $1/4}'
EOF

cd ${ANNOT}
sbatch ${ANNOT}/scripts/02d_merge_reads.sh
```

---

## Biological Questions — Module 2

1. **No-replicate design:** This dataset has one library per condition.
   Why is this sufficient for genome annotation but insufficient for
   differential expression analysis? What statistical property is missing
   when there is only a single observation per condition?

2. **Three conditions vs. more:** The original module assumed 15 runs
   across 5 conditions, but the actual dataset has only 3 conditions with
   no replicates. For the purpose of transcript evidence in gene prediction,
   what is the advantage of having multiple growth conditions rather than
   a single condition, even without replicates? Which functional categories
   of genes are most likely to be missed if only PD-5d were used?

3. **Paired-end vs. single-end for annotation:** These RNA-seq libraries
   are paired-end. For the purpose of gene model building, what advantage
   do paired-end reads provide over single-end reads of the same length?
   Think specifically about intron spanning.

4. **Library strandedness:** The GEO record does not explicitly state
   whether these libraries are stranded or unstranded. What experiment
   would you run to determine this empirically from the alignment data?
   What is the consequence of incorrectly specifying strandedness in
   the HISAT2 and StringTie steps of Module 4b?

---

## Take-Home Messages

> **A no-replicate RNA-seq dataset is still valuable for annotation.**
> Gene predictors need transcript boundaries, not statistical estimates
> of expression differences. Even a single library per condition can
> reveal exon structure for most expressed genes.

> **Trimming is required before annotation, not just before alignment.**
> Adapter sequences in RNA-seq reads can create artifactual splice
> junctions that generate incorrect exon structures if passed directly
> to gene predictors.

> **Merging reads from multiple conditions maximizes gene model coverage.**
> In contrast to differential expression analysis (where samples must
> stay separate), annotation benefits from treating all expressed
> transcripts as a single pool of evidence.

---

*Previous: [01 Setup](01_setup.md) | Next: [Module 3 → Genome Cleaning](03_genome_clean.md)*