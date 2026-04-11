# Module 3 — Genome Cleaning
## Preparing the Assembly for Gene Prediction

---

## Background

Before gene prediction can begin, the genome assembly must be prepared in
three ways:

1. **Sorting and filtering:** Funannotate2 requires scaffolds to be sorted
   by length (longest first) and assigns locus tags based on position.
   Very short scaffolds below a minimum length threshold are removed because
   they rarely contain complete gene models and inflate false-positive
   prediction counts.

2. **Header standardization:** Gene prediction tools (Augustus, SNAP,
   GlimmerHMM) are sensitive to FASTA header formatting. Long, complex
   headers with special characters can cause failures. Funannotate2's
   `clean` command normalizes headers.

3. **Repeat masking:** Transposable elements (TEs) and other repetitive
   sequences must be identified and soft-masked (converted to lowercase)
   before gene prediction. Ab initio predictors trained on coding sequence
   will otherwise predict spurious genes within TE bodies. Soft masking
   (lowercase, not replacing with N) allows gene predictors to *span*
   repeats when building models but prevents them from being seeded within
   repeat sequences.

> **Key point for *P. fructicola*:** With only 0.34% total repeat content
> and no full-length transposons, repeat masking has minimal quantitative
> impact on this genome. However, the step is included in all annotation
> workflows because it is essential for any repeat-rich organism, and you
> will apply this workflow to your plant genome where repeat masking is
> critical.

---

## Learning Objectives

- Run `funannotate2 clean` to sort, filter, and standardize the assembly
- Understand the repeat masking step and its role in annotation quality
- Verify the cleaned assembly before training begins
- Compare pre- and post-cleaning statistics

---

## Step 1: Inspect the Raw Assembly

```bash
cd ${ANNOT}

# Number of sequences
grep -c "^>" ${ANNOT}/00_genome/Pf_assembly_raw.fasta

# Sequence headers
grep "^>" ${ANNOT}/00_genome/Pf_assembly_raw.fasta | head -10

# Sequence lengths (sorted)
awk '
  /^>/ { if (seq) print length(seq), name
         name = $0; seq = "" }
  !/^>/ { seq = seq $0 }
  END  { print length(seq), name }
' ${ANNOT}/00_genome/Pf_assembly_raw.fasta | sort -k1,1rn

# Total assembly size
awk '!/^>/{t+=length($0)} END{printf "Total: %d bp (%.2f Mb)\n", t, t/1e6}' \
  ${ANNOT}/00_genome/Pf_assembly_raw.fasta
```

---

## Step 2: Run funannotate2 clean

```bash
cat > ${ANNOT}/scripts/03a_clean_genome.sh << 'EOF'
#!/bin/bash
#SBATCH --account=PAS3260
#SBATCH --job-name=f2_clean
#SBATCH --output=logs/f2_clean_%j.out
#SBATCH --error=logs/f2_clean_%j.err
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --mem=8G

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
  --env AUGUSTUS_CONFIG_PATH=/opt/augustus/config \
  --env FUNANNOTATE2_DB=/f2_db \
  ${F2_CONTAINER} \
  funannotate2 clean \
    -f /data/00_genome/Pf_assembly_raw.fasta \
    -o /data/00_genome/Pf_assembly_cleaned.fasta \
    -m 1000

echo "=== funannotate2 clean complete ==="

echo ""
echo "Sequences after cleaning:"
grep -c "^>" ${ANNOT}/00_genome/Pf_assembly_cleaned.fasta

echo ""
echo "New sequence headers:"
grep "^>" ${ANNOT}/00_genome/Pf_assembly_cleaned.fasta

echo ""
echo "Total assembly size after cleaning:"
awk '!/^>/{t+=length($0)} END{printf "Total: %d bp (%.2f Mb)\n", t, t/1e6}' \
  ${ANNOT}/00_genome/Pf_assembly_cleaned.fasta
EOF

cd ${ANNOT}
sbatch ${ANNOT}/scripts/03a_clean_genome.sh
```

**`funannotate2 clean` key parameters:**

| Parameter | Effect |
|-----------|--------|
| `-f` | Input FASTA file |
| `-o` | Output FASTA file |
| `-m 1000` | Remove sequences shorter than 1,000 bp |

---

## Step 3: Verify the Cleaned Assembly

```bash
echo "=== Before cleaning ==="
grep -c "^>" ${ANNOT}/00_genome/Pf_assembly_raw.fasta
echo "sequences"

echo ""
echo "=== After cleaning ==="
grep -c "^>" ${ANNOT}/00_genome/Pf_assembly_cleaned.fasta
echo "sequences"

echo ""
echo "Cleaned sequence headers:"
grep "^>" ${ANNOT}/00_genome/Pf_assembly_cleaned.fasta

echo ""
echo "Sequence lengths (sorted):"
awk '
  /^>/ { if (seq) print length(seq), name
         name = $0; seq = "" }
  !/^>/ { seq = seq $0 }
  END  { print length(seq), name }
' ${ANNOT}/00_genome/Pf_assembly_cleaned.fasta | sort -k1,1rn
```

**Expected output:** 5 sequences, headers simplified to `scaffold_1`
through `scaffold_5` or similar, total size ~18.99 Mb.

---

## Step 4: Repeat Masking

Funannotate2 uses **RepeatMasker** with a species-specific repeat library
built by **RepeatModeler**. Both are available inside the Funannotate2
container. Running masking explicitly as its own step gives you a reusable
soft-masked FASTA and visibility into the repeat landscape.

```bash
cat > ${ANNOT}/scripts/03b_repeat_mask.sh << 'EOF'
#!/bin/bash
#SBATCH --account=PAS3260
#SBATCH --job-name=repeat_mask
#SBATCH --output=logs/repeat_mask_%j.out
#SBATCH --error=logs/repeat_mask_%j.err
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem=32G

set -euo pipefail

user_name=Jonathan

ANNOT=/fs/scratch/PAS3260/${user_name}/Annotation
SHARED_F2=/fs/scratch/PAS3260/Team_Project/Containers/Funannotate2
F2_CONTAINER=${SHARED_F2}/funannotate2.sif
F2_DB=${SHARED_F2}/databases
AUGUSTUS_CONFIG=${ANNOT}/augustus_config

mkdir -p ${ANNOT}/00_genome/repeat_masking

# Step 1: Build a de novo repeat library with RepeatModeler
apptainer exec \
  --bind ${ANNOT}:/data \
  --bind ${F2_DB}:/f2_db \
  --bind ${AUGUSTUS_CONFIG}:/opt/augustus/config \
  --env AUGUSTUS_CONFIG_PATH=/opt/augustus/config \
  --env FUNANNOTATE2_DB=/f2_db \
  ${F2_CONTAINER} \
  bash -c "
    cd /data/00_genome/repeat_masking && \
    BuildDatabase \
      -name Pf_db \
      /data/00_genome/Pf_assembly_cleaned.fasta && \
    RepeatModeler \
      -database Pf_db \
      -threads 12 \
      -LTRStruct
  "

echo "RepeatModeler complete."
ls -lh ${ANNOT}/00_genome/repeat_masking/

# Step 2: Run RepeatMasker — soft-mask with de novo library + Dfam fungi
apptainer exec \
  --bind ${ANNOT}:/data \
  --bind ${F2_DB}:/f2_db \
  --bind ${AUGUSTUS_CONFIG}:/opt/augustus/config \
  --env AUGUSTUS_CONFIG_PATH=/opt/augustus/config \
  --env FUNANNOTATE2_DB=/f2_db \
  ${F2_CONTAINER} \
  RepeatMasker \
    -pa 12 \
    -xsmall \
    -lib /data/00_genome/repeat_masking/Pf_db-families.fa \
    -dir /data/00_genome/repeat_masking \
    /data/00_genome/Pf_assembly_cleaned.fasta

# Copy the soft-masked output to a clean name
cp ${ANNOT}/00_genome/repeat_masking/Pf_assembly_cleaned.fasta.masked \
   ${ANNOT}/00_genome/Pf_assembly_softmasked.fasta

echo "=== Repeat masking complete ==="
cat ${ANNOT}/00_genome/repeat_masking/Pf_assembly_cleaned.fasta.tbl
EOF

cd ${ANNOT}
sbatch ${ANNOT}/scripts/03b_repeat_mask.sh
```

---

## Step 5: Inspect the Repeat Masking Output

```bash
# RepeatMasker summary table
cat ${ANNOT}/00_genome/repeat_masking/Pf_assembly_cleaned.fasta.tbl

# Calculate masked fraction
awk '!/^>/{
  for(i=1;i<=length($0);i++){
    c=substr($0,i,1); total++
    if(c~/[acgtn]/) masked++
  }
}
END {
  printf "Total bases:   %d\n", total
  printf "Masked bases:  %d\n", masked
  printf "Masked %%:     %.2f%%\n", (masked/total)*100
}' ${ANNOT}/00_genome/Pf_assembly_softmasked.fasta
```

**Expected result:** For *P. fructicola*, the masked fraction should be
close to the published 0.34% repeat content.

> **Compare to your plant genome:** A plant genome with 50–85% TEs would
> show a dramatically different repeat landscape: extensive LTR
> retrotransposons (Copia, Gypsy) and DNA transposons would dominate the
> table. In that context, masking is not a minor preprocessing step —
> it determines whether gene prediction is even meaningful.

---

## Biological Questions — Module 3

1. **Soft masking vs. hard masking:** Funannotate2 uses soft masking
   (lowercase letters) rather than hard masking (N characters). Explain
   why soft masking is preferable for gene prediction. Under what
   specific scenario would hard masking be appropriate?

2. **Minimum contig length:** The `clean` command removes sequences shorter
   than 1,000 bp. For a draft assembly of a similar-sized fungal genome
   with 2,000 scaffolds ranging from 200 bp to 5 Mb, how many gene models
   might be lost by this filter, and why would you retain the filter anyway?

3. **De novo vs. database repeat libraries:** We ran RepeatModeler to build
   a *de novo* library, then combined it with the Dfam fungal database in
   RepeatMasker. What is the advantage of each approach? For which element
   types is the de novo library essential?

4. **Repeat content and annotation quality:** The *P. fructicola* genome
   has 0.34% repeats. A related Dothideomycete like *Alternaria alternata*
   has ~25% repeats. Predict how annotation completeness (as measured by
   BUSCO proteins mode) would differ between these two organisms at the
   same sequencing depth and with the same annotation pipeline, and explain
   the mechanistic reasons.

---

## Take-Home Messages

> **Genome cleaning is not optional.** Funannotate2 relies on consistent,
> predictable FASTA headers and scaffold ordering. Skipping this step
> causes failures at multiple downstream points.

> **Soft masking is how gene prediction pipelines handle repeats.**
> The lowercase signal tells Augustus and SNAP: "do not start a gene
> model here, but you may extend through here."

> **The repeat landscape is a biological finding, not just a preprocessing
> step.** The RepeatMasker `.tbl` output tells you what classes of elements
> colonize the genome — or, in *P. fructicola*'s case, what has been lost.

---

*Previous: [02 RNA-seq Evidence](02_rnaseq_evidence.md) | Next: [Module 4 → Training](04_train.md)*