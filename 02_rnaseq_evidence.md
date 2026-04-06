# Module 2 — RNA-seq Evidence Preparation
## Downloading, Trimming, and Validating the SRP166999 Dataset

---

## Background

High-quality RNA-seq evidence is the single most effective way to improve
gene model accuracy in a eukaryotic genome annotation. Transcript boundaries
resolved from direct RNA sequencing are far more reliable than ab initio
predictions alone, particularly for:

- Defining precise exon–intron boundaries
- Distinguishing real genes from spurious ORFs predicted in intergenic regions
- Capturing alternatively spliced isoforms
- Providing expression evidence for lowly conserved, lineage-specific genes

The Wang et al. (2020) RNA-seq dataset (SRA study **SRP166999**) covers
*P. fructicola* grown under five biologically relevant conditions in
triplicate (15 runs total). Using all 15 conditions maximizes gene model
coverage: a gene that is silenced in liquid culture may be highly expressed
during apple fruit colonization. Providing transcripts from all conditions
ensures the annotation reflects the full expressed gene complement.

> **Connection to Module 6 of the previous tutorial:** In that module you
> processed a single run (SRR8119502) through alignment and quantification.
> Here, you download and trim the complete 15-run dataset that will serve
> as transcript evidence for gene prediction.

---

## Learning Objectives

- Download all 15 RNA-seq runs from SRP166999 using SRA Toolkit
- Perform adapter trimming with Trimmomatic
- Validate read quality with FastQC
- Merge paired reads into a format ready for Funannotate2 training

---

## The SRP166999 Dataset

| Condition | Code | SRR accessions |
|-----------|------|----------------|
| PDB liquid broth, 5 days | PDB-5d | SRR8119502, SRR8119503, SRR8119504 |
| PDB liquid broth, 15 days | PDB-15d | SRR8119505, SRR8119506, SRR8119507 |
| PDB + PEG 6000 (osmotic stress), 15 days | PDB-PEG | SRR8119508, SRR8119509, SRR8119510 |
| Potato dextrose agar (solid), 15 days | PDA-15d | SRR8119511, SRR8119512, SRR8119513 |
| Apple fruit inoculation, 15 days | Apple | SRR8119514, SRR8119515, SRR8119516 |

> **Verify before running:** Always confirm run accessions in the SRA Run
> Selector before submitting download jobs. Go to
> <https://www.ncbi.nlm.nih.gov/sra?term=SRP166999>, click
> "Send to → File → Accession List" to download the current list.

---

## Step 1: Create the Accession List File

```bash
cat > ${ANNOT}/01_rnaseq/SRP166999_accessions.txt << 'EOF'
SRR8119502
SRR8119503
SRR8119504
SRR8119505
SRR8119506
SRR8119507
SRR8119508
SRR8119509
SRR8119510
SRR8119511
SRR8119512
SRR8119513
SRR8119514
SRR8119515
SRR8119516
EOF

echo "Accession list created:"
cat ${ANNOT}/01_rnaseq/SRP166999_accessions.txt
```

---

## Step 2: Download All 15 Runs via SLURM Job Array

Downloading 15 SRA runs in series would be slow. A SLURM job array
runs all 15 downloads concurrently (up to 5 at a time):

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
#SBATCH --array=1-15%5

set -euo pipefail

CONTAINERS=/fs/scratch/PAS3260/${user_name}/Annotation/containers
ANNOT=/fs/scratch/PAS3260/${user_name}/Annotation
OUTDIR=${ANNOT}/01_rnaseq/raw
mkdir -p ${OUTDIR}/tmp

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${ANNOT}/01_rnaseq/SRP166999_accessions.txt)
echo "Downloading: ${SAMPLE}  (task ${SLURM_ARRAY_TASK_ID}/15)"

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

Once all 15 jobs complete, verify all files:

```bash
echo "=== Downloaded FASTQ files ==="
ls -lh ${ANNOT}/01_rnaseq/raw/SRR81195*.fastq.gz | wc -l
echo "files (expected: 30 — R1 + R2 for each of 15 runs)"

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
#SBATCH --array=1-15

set -euo pipefail

CONTAINERS=/fs/scratch/PAS3260/${user_name}/Annotation/containers
ANNOT=/fs/scratch/PAS3260/${user_name}/Annotation

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${ANNOT}/01_rnaseq/SRP166999_accessions.txt)
OUTDIR=${ANNOT}/01_rnaseq/trimmed
mkdir -p ${OUTDIR}

echo "Trimming: ${SAMPLE}"

ADAPTER_PATH=$(apptainer exec --contain ${CONTAINERS}/trimmomatic_0.40.sif \
  find / -name "TruSeq3-PE-2.fa" 2>/dev/null | head -1)

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
    ILLUMINACLIP:${ADAPTER_PATH}:2:30:10:2:True \
    LEADING:3 \
    TRAILING:3 \
    SLIDINGWINDOW:4:20 \
    MINLEN:36

echo "${SAMPLE} trimming complete."
EOF

cd ${ANNOT}
sbatch ${ANNOT}/scripts/02b_trimmomatic_array.sh
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

CONTAINERS=/fs/scratch/PAS3260/${user_name}/Annotation/containers
ANNOT=/fs/scratch/PAS3260/${user_name}/Annotation

mkdir -p ${ANNOT}/01_rnaseq/qc

apptainer exec \
  --bind ${ANNOT}:/data \
  ${CONTAINERS}/fastqc_0.12.1.sif \
  fastqc \
    /data/01_rnaseq/trimmed/*_paired.fastq.gz \
    --outdir /data/01_rnaseq/qc \
    --threads 8

echo "FastQC complete."
ls ${ANNOT}/01_rnaseq/qc/*.html | wc -l
echo "HTML reports generated (expected: 30)"
EOF

cd ${ANNOT}
sbatch ${ANNOT}/scripts/02c_fastqc.sh
```

---

## Step 5: Merge Paired Files for Funannotate2

Funannotate2's `train` command takes merged R1 and R2 files. Concatenating
all 15 pairs ensures all conditions contribute to ab initio predictor training:

```bash
cat > ${ANNOT}/scripts/02d_merge_reads.sh << 'EOF'
#!/bin/bash
#SBATCH --account=PAS3260
#SBATCH --job-name=merge_reads
#SBATCH --output=logs/merge_reads_%j.out
#SBATCH --error=logs/merge_reads_%j.err
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=32G

set -euo pipefail

ANNOT=/fs/scratch/PAS3260/${user_name}/Annotation
TRIMMED=${ANNOT}/01_rnaseq/trimmed

echo "Merging R1 reads from all 15 samples..."
cat ${TRIMMED}/SRR8119502_1_paired.fastq.gz \
    ${TRIMMED}/SRR8119503_1_paired.fastq.gz \
    ${TRIMMED}/SRR8119504_1_paired.fastq.gz \
    ${TRIMMED}/SRR8119505_1_paired.fastq.gz \
    ${TRIMMED}/SRR8119506_1_paired.fastq.gz \
    ${TRIMMED}/SRR8119507_1_paired.fastq.gz \
    ${TRIMMED}/SRR8119508_1_paired.fastq.gz \
    ${TRIMMED}/SRR8119509_1_paired.fastq.gz \
    ${TRIMMED}/SRR8119510_1_paired.fastq.gz \
    ${TRIMMED}/SRR8119511_1_paired.fastq.gz \
    ${TRIMMED}/SRR8119512_1_paired.fastq.gz \
    ${TRIMMED}/SRR8119513_1_paired.fastq.gz \
    ${TRIMMED}/SRR8119514_1_paired.fastq.gz \
    ${TRIMMED}/SRR8119515_1_paired.fastq.gz \
    ${TRIMMED}/SRR8119516_1_paired.fastq.gz \
    > ${ANNOT}/01_rnaseq/all_samples_R1.fastq.gz

echo "Merging R2 reads from all 15 samples..."
cat ${TRIMMED}/SRR8119502_2_paired.fastq.gz \
    ${TRIMMED}/SRR8119503_2_paired.fastq.gz \
    ${TRIMMED}/SRR8119504_2_paired.fastq.gz \
    ${TRIMMED}/SRR8119505_2_paired.fastq.gz \
    ${TRIMMED}/SRR8119506_2_paired.fastq.gz \
    ${TRIMMED}/SRR8119507_2_paired.fastq.gz \
    ${TRIMMED}/SRR8119508_2_paired.fastq.gz \
    ${TRIMMED}/SRR8119509_2_paired.fastq.gz \
    ${TRIMMED}/SRR8119510_2_paired.fastq.gz \
    ${TRIMMED}/SRR8119511_2_paired.fastq.gz \
    ${TRIMMED}/SRR8119512_2_paired.fastq.gz \
    ${TRIMMED}/SRR8119513_2_paired.fastq.gz \
    ${TRIMMED}/SRR8119514_2_paired.fastq.gz \
    ${TRIMMED}/SRR8119515_2_paired.fastq.gz \
    ${TRIMMED}/SRR8119516_2_paired.fastq.gz \
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

1. **Study design:** The RNA-seq data spans five growth conditions representing
   a progression from laboratory culture to host-associated growth. Why does
   including all five conditions, rather than just a single condition, improve
   genome annotation quality? Which functional categories of genes might be
   missed if you used only the PDB liquid culture samples?

2. **Read depth per sample:** Examine the read counts you printed in Step 2.
   Are the read counts approximately equal across all 15 samples? What would
   a large imbalance between samples imply for the merged transcript evidence?

3. **Paired-end vs. single-end for annotation:** These RNA-seq libraries are
   paired-end (125 bp). For the purpose of gene model building, what
   advantage do paired-end reads provide over single-end reads of the same
   length? Think about intron spanning.

4. **Evidence completeness:** The paper reports 8,072 protein-coding genes
   with 94.8% supported by ≥10 FPKM. If all 15 RNA-seq runs are used as
   evidence, what upper bound would you predict for the fraction of gene
   models with direct RNA-seq support in the final annotation?

---

## Take-Home Messages

> **RNA-seq evidence from multiple biological conditions is fundamentally
> different from single-condition data.** Genes that are condition-specific
> (e.g., expressed only on the apple surface) will only contribute transcript
> evidence if that condition is included.

> **Trimming is required before annotation, not just before alignment.**
> Adapter sequences in RNA-seq reads can create artifactual splice junctions
> that generate incorrect exon structures if passed directly to gene predictors.

> **Merging reads from multiple conditions is the correct approach for
> annotation.** In contrast to differential expression analysis (where
> samples must stay separate), annotation benefits from treating all
> expressed transcripts as a single pool of evidence.

---

*Previous: [01 Setup](01_setup.md) | Next: [Module 3 → Genome Cleaning](03_genome_clean.md)*
