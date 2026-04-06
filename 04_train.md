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

# Confirm augustus_config is writable in your directory
ls -la ${ANNOT}/augustus_config/
```

---

## Step 2: Run funannotate2 train

This is the most computationally intensive step in the annotation pipeline.
Expect **2–4 hours** for a 19 Mb fungal genome with 15 RNA-seq samples
using 16 CPUs.

```bash
cat > ${ANNOT}/scripts/04a_f2_train.sh << 'EOF'
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

ANNOT=/fs/scratch/PAS3260/${user_name}/Annotation
SHARED_F2=/fs/scratch/PAS3260/Team_Project/Containers/Funannotate2
F2_CONTAINER=${SHARED_F2}/funannotate2.sif
F2_DB=${SHARED_F2}/databases
AUGUSTUS_CONFIG=${ANNOT}/augustus_config

apptainer exec \
  --bind ${ANNOT}:/data \
  --bind ${F2_DB}:/f2_db \
  --bind ${AUGUSTUS_CONFIG}:/opt/augustus_config \
  --bind ${SHARED_F2}/gmes_linux_64_4:/gmes_linux_64_4 \
  --bind ${SHARED_F2}/signalp-6-package:/signalp-6-package \
  --env AUGUSTUS_CONFIG_PATH=/opt/augustus_config \
  ${F2_CONTAINER} \
  funannotate2 train \
    -f /data/00_genome/Pf_assembly_softmasked.fasta \
    -s "Peltaster fructicola" \
    --strain LNHT1506 \
    -o /data/02_funannotate \
    --left  /data/01_rnaseq/all_samples_R1.fastq.gz \
    --right /data/01_rnaseq/all_samples_R2.fastq.gz \
    --stranded RF \
    --database /f2_db \
    --busco_db dothideomycetes_odb12 \
    --cpus 16

echo "=== funannotate2 train complete: $(date) ==="
EOF

cd ${ANNOT}
sbatch ${ANNOT}/scripts/04a_f2_train.sh
squeue -u ${USER}
```

**Key parameters:**

| Parameter | Value | Purpose |
|-----------|-------|---------|
| `-f` | soft-masked FASTA | Input genome |
| `-s` | `"Peltaster fructicola"` | Species name used in locus tags |
| `--strain` | `LNHT1506` | Strain identifier |
| `-o` | `02_funannotate` | Output directory (used by all Funannotate2 stages) |
| `--left / --right` | merged R1/R2 | Paired-end RNA-seq reads |
| `--stranded RF` | `RF` | dUTP library strand orientation (reverse-forward) |
| `--busco_db` | `dothideomycetes_odb12` | BUSCO lineage for training set construction |

> **Strandedness note:** Most Illumina libraries prepared with the dUTP
> method (TruSeq Stranded kits) are **RF** stranded. If you are unsure,
> you can verify by running RSeQC's `infer_experiment.py` on a test
> alignment, or check the GEO record for SRP166999.

---

## Step 3: Monitor the Training Job

```bash
# Check job status
squeue -u ${USER}

# Follow the log in real time (once the job starts)
tail -f ${ANNOT}/logs/f2_train_*.out
```

**Expected stages in the log:**

```
[INFO] Running BUSCOlite to identify training loci...
[INFO] Found N complete BUSCO loci for training
[INFO] Aligning RNA-seq reads with HISAT2...
[INFO] Assembling transcripts with StringTie...
[INFO] Training Augustus...
[INFO] Training SNAP...
[INFO] Training GlimmerHMM...
[INFO] Training GeneMark-ES...
[INFO] Saving parameters to params.json
[INFO] funannotate2 train complete
```

---

## Step 4: Inspect the Training Output

```bash
# Output directory structure
ls -la ${ANNOT}/02_funannotate/

# The key output — stored training parameters for predict
ls -lh ${ANNOT}/02_funannotate/params.json

# Verify params.json is valid and contains expected fields
apptainer exec \
  --bind ${ANNOT}:/data \
  --bind ${F2_DB}:/f2_db \
  --bind ${AUGUSTUS_CONFIG}:/opt/augustus_config \
  --bind ${SHARED_F2}/gmes_linux_64_4:/gmes_linux_64_4 \
  --bind ${SHARED_F2}/signalp-6-package:/signalp-6-package \
  --env AUGUSTUS_CONFIG_PATH=/opt/augustus_config \
  ${F2_CONTAINER} \
  python3 -c "
import json
with open('/data/02_funannotate/params.json') as f:
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

## Step 5: Review the HISAT2 Alignment Summary

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
