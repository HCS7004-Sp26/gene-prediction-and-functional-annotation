# Module 4 — Training Ab Initio Gene Predictors
## Building Species-Specific Gene Models for Augustus, SNAP, GlimmerHMM, and GeneMark

---

## Background

Ab initio gene predictors work by applying statistical models of gene
structure — models that encode the length distributions, sequence
composition, and splice-site signals characteristic of genes in a
particular organism. Using a model trained on *Aspergillus* to annotate
*P. fructicola* would produce inaccurate gene structures because the two
organisms differ substantially in:

- **Intron length:** *P. fructicola* median intron = 50 bp vs. 80–150 bp
  in many other fungi
- **Intron density:** 1.36 introns per gene vs. 2–4 in typical Ascomycetes
- **Codon usage:** organism-specific biases affect the likelihood scores
  assigned by the predictor
- **Gene density:** compact intergenic regions change the prior probability
  of gene boundaries

`funannotate2 train` automates the complete predictor training workflow:
it uses **BUSCOlite** to identify a high-confidence set of conserved
single-copy orthologs, builds gene models from those orthologs, and uses
those models to train Augustus, SNAP, and GlimmerHMM. When RNA-seq data
are provided, transcript alignments augment the training set.

> **GeneMark-ES is included in this class installation.** The shared
> Funannotate2 container has GeneMark-ES pre-installed at
> `${SHARED_F2}/gmes_linux_64_4/`. This licensed tool provides a fourth
> ab initio prediction track, improving gene model consensus. It is
> bound into the container automatically using the standard bind pattern.

---

## Learning Objectives

- Understand why species-specific training improves annotation accuracy
- Run `funannotate2 train` with the soft-masked genome and RNA-seq reads
- Interpret the training output and verify that all four predictors trained
- Understand the relationship between training set quality and final
  annotation completeness

---

## The funannotate2 train Workflow

```
Soft-masked genome
        │
        ▼
BUSCOlite: find conserved orthologs
        │
        ▼
miniprot: align fungal proteins to genome
        │
        ├── High-confidence gene models → train Augustus
        ├── High-confidence gene models → train SNAP
        ├── High-confidence gene models → train GlimmerHMM
        └── High-confidence gene models → train GeneMark-ES
        │
        ▼  (RNA-seq provided)
HISAT2: align merged reads to genome
        │
        ▼
StringTie: assemble transcripts
        │
        ▼
Transcript models augment training set
        │
        ▼
params.json  ← stored training parameters for predict
```

---

## Step 1: Verify Inputs

```bash
# Confirm soft-masked genome exists
ls -lh ${ANNOT}/00_genome/Pf_assembly_softmasked.fasta

# Confirm merged RNA-seq reads exist
ls -lh ${ANNOT}/01_rnaseq/all_samples_R*.fastq.gz

# Quick count check
echo "R1 reads:"
zcat ${ANNOT}/01_rnaseq/all_samples_R1.fastq.gz | wc -l | awk '{print $1/4}'
echo "R2 reads:"
zcat ${ANNOT}/01_rnaseq/all_samples_R2.fastq.gz | wc -l | awk '{print $1/4}'

# Confirm augustus_config is writable in your directory
ls -la ${ANNOT}/augustus_config/
```

---

## Step 2: Run Transcript Assembly

Before training, assemble RNA-seq transcripts that will augment the
training set and later serve as evidence for gene prediction.

```bash
cat > ${ANNOT}/scripts/04a_transcript_assembly.sh << 'EOF'
#!/bin/bash
#SBATCH --account=PAS3260
#SBATCH --job-name=transcript_assembly
#SBATCH --output=logs/transcript_assembly_%j.out
#SBATCH --error=logs/transcript_assembly_%j.err
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=64G

set -euo pipefail

user_name=Jonathan

CONTAINERS=/fs/scratch/PAS3260/${user_name}/Annotation/containers
ANNOT=/fs/scratch/PAS3260/${user_name}/Annotation
GENOME=${ANNOT}/00_genome/Pf_assembly_softmasked.fasta
R1=${ANNOT}/01_rnaseq/all_samples_R1.fastq.gz
R2=${ANNOT}/01_rnaseq/all_samples_R2.fastq.gz
OUTDIR=${ANNOT}/01_rnaseq/transcripts
mkdir -p ${OUTDIR}

# --- Step 1: Build HISAT2 genome index ---
echo "=== Building HISAT2 index: $(date) ==="
apptainer exec \
  --bind ${ANNOT}:/data \
  ${CONTAINERS}/hisat2_2.2.2.sif \
  hisat2-build \
    -p 16 \
    /data/00_genome/Pf_assembly_softmasked.fasta \
    /data/01_rnaseq/transcripts/Pf_genome

# --- Step 2: Align merged reads ---
echo "=== HISAT2 alignment: $(date) ==="
apptainer exec \
  --bind ${ANNOT}:/data \
  ${CONTAINERS}/hisat2_2.2.2.sif \
  bash -c "hisat2 \
    -p 16 \
    --rna-strandness RF \
    -x /data/01_rnaseq/transcripts/Pf_genome \
    -1 /data/01_rnaseq/all_samples_R1.fastq.gz \
    -2 /data/01_rnaseq/all_samples_R2.fastq.gz \
    --dta \
    2>/data/01_rnaseq/transcripts/hisat2_summary.txt \
  | samtools sort -@ 16 -o /data/01_rnaseq/transcripts/Pf_rnaseq.bam \
  && samtools index /data/01_rnaseq/transcripts/Pf_rnaseq.bam"

echo "=== HISAT2 alignment rate ==="
grep "overall alignment rate" ${OUTDIR}/hisat2_summary.txt

# --- Step 3: Assemble transcripts with StringTie ---
echo "=== StringTie assembly: $(date) ==="
apptainer exec \
  --bind ${ANNOT}:/data \
  ${CONTAINERS}/stringtie_3.0.3.sif \
  stringtie \
    /data/01_rnaseq/transcripts/Pf_rnaseq.bam \
    --rf \
    -p 16 \
    -o /data/01_rnaseq/transcripts/Pf_transcripts.gtf

# --- Step 4: Extract transcript sequences with gffread ---
echo "=== gffread: extracting transcript FASTA: $(date) ==="
apptainer exec \
  --bind ${ANNOT}:/data \
  ${CONTAINERS}/gffread_0.12.7.sif \
  gffread \
    /data/01_rnaseq/transcripts/Pf_transcripts.gtf \
    -g /data/00_genome/Pf_assembly_softmasked.fasta \
    -w /data/01_rnaseq/transcripts/Pf_transcripts.fasta

echo "=== Transcript assembly complete: $(date) ==="
echo "Transcript count:"
grep -c "^>" ${OUTDIR}/Pf_transcripts.fasta
ls -lh ${OUTDIR}/Pf_transcripts.fasta
EOF

cd ${ANNOT}
sbatch ${ANNOT}/scripts/04a_transcript_assembly.sh
squeue -u ${USER}
```

---

## Step 3: Run funannotate2 train

This is the most computationally intensive step in the annotation pipeline.
Expect **2–4 hours** for a 19 Mb fungal genome with 15 RNA-seq samples
using 16 CPUs.

```bash
cat > ${ANNOT}/scripts/04b_f2_train.sh << 'EOF'
#!/bin/bash
#SBATCH --account=PAS3260
#SBATCH --job-name=f2_train
#SBATCH --output=logs/f2_train_%j.out
#SBATCH --error=logs/f2_train_%j.err
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=64G

set -euo pipefail

user_name=Jonathan

ANNOT=/fs/scratch/PAS3260/${user_name}/Annotation
SHARED_F2=/fs/scratch/PAS3260/Team_Project/Containers/Funannotate2
F2_CONTAINER=${SHARED_F2}/funannotate2.sif
F2_DB=${SHARED_F2}/databases
AUGUSTUS_CONFIG=${ANNOT}/augustus_config

apptainer exec \
  --bind ${ANNOT}:/data \
  --bind ${F2_DB}:/f2_db \
  --bind ${AUGUSTUS_CONFIG}:/opt/augustus/config \
  --bind ${SHARED_F2}/gmes_linux_64_4:/gmes_linux_64_4 \
  --bind ~/.gm_key:/root/.gm_key \
  --env AUGUSTUS_CONFIG_PATH=/opt/augustus/config \
  --env FUNANNOTATE2_DB=/f2_db \
  ${F2_CONTAINER} \
  funannotate2 train \
    -f /data/00_genome/Pf_assembly_softmasked.fasta \
    -s "Peltaster fructicola" \
    --strain LNHT1506 \
    -o /data/02_funannotate \
    --busco-lineage dothideomycetes \
    --genemark-mode guided \
    --cpus 16

echo "=== funannotate2 train complete: $(date) ==="
ls -lh ${ANNOT}/02_funannotate/
EOF

cd ${ANNOT}
sbatch ${ANNOT}/scripts/04b_f2_train.sh
```

**Key parameters:**

| Parameter | Value | Purpose |
|-----------|-------|---------|
| `-f` | soft-masked FASTA | Input genome |
| `-s` | `"Peltaster fructicola"` | Species name used in locus tags |
| `--strain` | `LNHT1506` | Strain identifier |
| `-o` | `02_funannotate` | Output directory (used by all Funannotate2 stages) |
| `--busco-lineage` | `dothideomycetes` | BUSCO lineage for training set construction |
| `--genemark-mode` | `guided` | GeneMark training mode: guided (full+hints) |

> **Strandedness note:** Most Illumina libraries prepared with the dUTP
> method (TruSeq Stranded kits) are **RF** stranded. If you are unsure,
> you can verify by running RSeQC's `infer_experiment.py` on a test
> alignment, or check the GEO record for SRP166999.

---

## Step 4: Monitor the Training Job

```bash
# Check job status
squeue -u ${USER}
```

**Expected stages in the log:**

```
[INFO] Running BUSCOlite to identify training loci...
[INFO] Found N complete BUSCO loci for training
[INFO] Training Augustus...
[INFO] Training SNAP...
[INFO] Training GlimmerHMM...
[INFO] Training GeneMark-ES...
[INFO] Saving parameters to params.json
[INFO] funannotate2 train complete
```

---

## Step 5: Inspect the Training Output

```bash
# Output directory structure
ls -la ${ANNOT}/02_funannotate/

# The key output — stored training parameters for predict
ls -lh ${ANNOT}/02_funannotate/train_results/*params.json

# Verify params.json is valid and contains expected fields
apptainer exec \
  --bind ${ANNOT}:/data \
  --bind ${F2_DB}:/f2_db \
  --bind ${AUGUSTUS_CONFIG}:/opt/augustus/config \
  --bind ${SHARED_F2}/gmes_linux_64_4:/gmes_linux_64_4 \
  --bind ~/.gm_key:/root/.gm_key \
  --env AUGUSTUS_CONFIG_PATH=/opt/augustus/config \
  --env FUNANNOTATE2_DB=/f2_db \
  ${F2_CONTAINER} \
  python3 -c "
import json
with open('${ANNOT}/02_funannotate/train_results/Peltaster_fructicola_LNHT1506.params.json') as f:
    params = json.load(f)
print('Keys:', list(params.keys()))
print('Species:', params.get('species', 'not found'))
print('Augustus species slug:', params.get('augustus_species', 'not found'))
print('Training loci used:', params.get('training_set', {}).get('n_loci', 'not found'))
"
```

> **Check that your trained Augustus species parameters were written**
> to your writable copy of augustus_config:
> ```bash
> ls ${ANNOT}/augustus_config/species/ | grep -i peltaster
> ```
> You should see a directory named something like `peltaster_fructicola`.
> If it is missing, the training step did not bind the correct
> `augustus_config` path.

---

## Step 6: Review the HISAT2 Alignment Summary

```bash
# Find and print the HISAT2 summary from training
grep -r "overall alignment rate" ${ANNOT}/logs/f2_train_*.out | tail -5
```

**Expected:** overall alignment rate >90% for a well-assembled genome
with high-quality RNA-seq data from the cognate organism.

---

## Biological Questions — Module 4

1. **Training set size:** For `dothideomycetes_odb12` (4,643 orthologs),
   what fraction would you expect to be found in *P. fructicola* based on
   your BUSCO results from the previous tutorial? Estimate the number of
   loci available for training.

2. **Stranded vs. unstranded RNA-seq:** The `--stranded RF` flag tells
   HISAT2 and StringTie the orientation of the cDNA library. What error
   could result if you used `--stranded RF` on an unstranded library? What
   would happen to the transcript assemblies and, consequently, to the
   training set?

3. **Why train instead of using a pre-trained model?** Funannotate2 includes
   pre-trained Augustus models for dozens of fungal species. For an
   annotation of *P. fructicola*, when would it be acceptable to use the
   pre-trained *Aspergillus nidulans* model rather than training a new one?
   What features of the target genome would make pre-trained models
   inadequate?

4. **GeneMark-ES vs. other ab initio predictors:** GeneMark-ES uses a
   different statistical framework (inhomogeneous Markov model) from
   Augustus (generalized hidden Markov model). Why does having four
   independent ab initio predictions improve the Evidence Modeler
   consensus, rather than just using the best single predictor?

---

## Take-Home Messages

> **Species-specific training is the single largest contributor to ab
> initio prediction accuracy.** For a compact genome with unusual intron
> statistics, using a heterologous training set produces systematically
> wrong intron boundary predictions that cascade through every downstream
> analysis.

> **The augustus_config directory must be writable during training.** This
> is why each student has their own copy. The trained species parameters
> written there are what `funannotate2 predict` will use in the next module.

> **The RNA-seq alignment rate during training is a diagnostic.** A low
> rate (< 70%) suggests a problem with the genome assembly, the RNA-seq
> data quality, or the strandedness setting — all of which should be
> investigated before proceeding to prediction.

---

*Previous: [03 Genome Cleaning](03_genome_clean.md) | Next: [Module 5 → Gene Prediction](05_predict.md)*