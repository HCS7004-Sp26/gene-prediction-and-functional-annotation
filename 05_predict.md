# Module 5 — Gene Prediction
## Generating Consensus Gene Models with Evidence Modeler

---

## Background

With the genome cleaned, repeats masked, and ab initio predictors trained,
the prediction step integrates all evidence streams — ab initio predictions,
RNA-seq transcript evidence, and protein homology — into a single weighted
consensus using **Evidence Modeler (EVM)**.

Evidence Modeler assigns different weights to each evidence type:

| Evidence | Weight rationale |
|----------|----------------|
| High-confidence ab initio (Augustus, GeneMark) | High — well-trained on this organism |
| RNA-seq transcript alignments | Very high — direct experimental evidence |
| Protein homology (miniprot alignments) | High — functionally validated sequence |
| Lower-confidence ab initio (SNAP, GlimmerHMM) | Moderate — useful as corroboration |

> **Critical concept:** `funannotate2 predict` requires the output
> directory from `funannotate2 train` (`-i` flag). It continues building
> within that same directory structure, adding prediction results alongside
> the training data. This is why the output directory is the same
> across all Funannotate2 commands.

---

## Learning Objectives

- Provide protein evidence sequences to `funannotate2 predict`
- Run the full Evidence Modeler consensus prediction
- Assess prediction quality with BUSCO (proteins mode)
- Inspect the output GFF3 and FASTA files
- Compare the Funannotate2 gene set to the published NCBI annotation

---

## Step 1: Verify the protein file exist

# Verify the protein evidence file is present in the shared database
```bash
ls -lh /fs/scratch/PAS3260/Team_Project/Containers/Funannotate2/databases/uniprot_sprot.fasta
grep -c "^>" /fs/scratch/PAS3260/Team_Project/Containers/Funannotate2/databases/uniprot_sprot.fasta
echo "UniProt/Swiss-Prot sequences available as protein evidence"
```

# Verify transcript evidence from Module 4b is present
```bash
ls -lh ${ANNOT}/01_rnaseq/transcripts/Pf_transcripts.fasta
grep -c "^>" ${ANNOT}/01_rnaseq/transcripts/Pf_transcripts.fasta
echo "StringTie transcript sequences available"
```
---

## Step 2: Run funannotate2 predict

```bash
cat > ${ANNOT}/scripts/05b_f2_predict.sh << 'EOF'
#!/bin/bash
#SBATCH --account=PAS3260
#SBATCH --job-name=f2_predict
#SBATCH --output=logs/f2_predict_%j.out
#SBATCH --error=logs/f2_predict_%j.err
#SBATCH --time=08:00:00
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
  --bind ${AUGUSTUS_CONFIG}:/opt/augustus_config \
  --bind ${SHARED_F2}/gmes_linux_64_4:/gmes_linux_64_4 \
  --env AUGUSTUS_CONFIG_PATH=/opt/augustus_config \
  --env FUNANNOTATE2_DB=/f2_db \
  ${F2_CONTAINER} \
  funannotate2 predict \
    -i /data/02_funannotate \
    -s "Peltaster fructicola" \
    --strain LNHT1506 \
    --bam /data/01_rnaseq/transcripts/Pf_rnaseq.bam \
    --bam-library RF \
    -ps /f2_db/uniprot_sprot.fasta \
    -ts /data/01_rnaseq/transcripts/Pf_transcripts.fasta \
    --min-intron 40 \
    --locus-tag PELF \
    --cpus 16

echo "=== funannotate2 predict complete: $(date) ==="
echo ""
echo "Output files:"
ls -lh ${ANNOT}/02_funannotate/predict_results/
EOF

cd ${ANNOT}
sbatch ${ANNOT}/scripts/05b_f2_predict.sh
squeue -u ${USER}
```

**Key parameters:**

| Parameter | Value | Purpose |
|---|---|---|
| `-i` | `02_funannotate` | Input/output directory from `train` |
| `-s` | `"Peltaster fructicola"` | Species name for locus tags and output naming |
| `--strain` | `LNHT1506` | Strain identifier |
| `--locus-tag` | `PELF` | Species-specific locus tag prefix; required unique for NCBI submission |
| `--bam` | `Pf_rnaseq.bam` | Coordinate-sorted RNA-seq BAM from Module 4b for splice hints and UTR refinement |
| `--bam-library` | `RF` | dUTP/reverse-stranded library orientation (provisional — verify with RSeQC if needed) |
| `-ps` | `uniprot_sprot.fasta` | UniProt/Swiss-Prot proteins for miniprot cross-species homology alignment |
| `-ts` | `Pf_transcripts.fasta` | StringTie transcripts from Module 4b for same-species exon boundary evidence |
| `--min-intron` | `40` | Excludes spurious very short intron calls inconsistent with *P. fructicola*'s 50 bp median intron |
| `--cpus` | `16` | Parallel threads for ab initio predictors and alignment steps |
---

## Step 3: Inspect Prediction Output

```bash
PREDICT_DIR=${ANNOT}/02_funannotate/predict_results

ls -lh ${PREDICT_DIR}/

echo "=== Key output files ==="
echo ""
echo "GFF3 annotation:"
ls -lh ${PREDICT_DIR}/*.gff3

echo ""
echo "Protein sequences:"
ls -lh ${PREDICT_DIR}/*.proteins.fa

echo ""
echo "Transcript sequences:"
ls -lh ${PREDICT_DIR}/*.transcripts.fa

echo ""
echo "GenBank format:"
ls -lh ${PREDICT_DIR}/*.gbk
```

---

## Step 4: Summary Statistics

```bash
PREDICT_DIR=${ANNOT}/02_funannotate/predict_results

echo "=== Gene model statistics ==="
grep -v "^#" ${PREDICT_DIR}/*.gff3 | awk '$3 == "gene"' | wc -l
echo "gene models"

echo ""
grep -v "^#" ${PREDICT_DIR}/*.gff3 | awk '$3 == "mRNA"' | wc -l
echo "mRNA (protein-coding) models"

echo ""
grep -c "^>" ${PREDICT_DIR}/*.proteins.fa
echo "protein sequences"

echo ""
echo "Introns per gene distribution:"
grep -v "^#" ${PREDICT_DIR}/*.gff3 | awk '$3 == "intron"' | \
  awk '{match($9, /Parent=([^;]+)/, arr); print arr[1]}' | \
  sort | uniq -c | awk '{print $1}' | \
  sort -n | uniq -c | \
  awk '{printf "  %s introns: %s genes\n", $2, $1}'

echo ""
echo "Gene length statistics:"
grep -v "^#" ${PREDICT_DIR}/*.gff3 | awk '$3 == "gene" {print $5 - $4 + 1}' | \
  sort -n | awk '
    NR==1{min=$1}
    {sum+=$1; n++}
    END{printf "  Count: %d  Min: %d bp  Mean: %.0f bp\n", n, min, sum/n}'
```

---

## Step 5: BUSCO Assessment of Predicted Proteins

```bash
cat > ${ANNOT}/scripts/05c_busco_predict.sh << 'EOF'
#!/bin/bash
#SBATCH --account=PAS3260
#SBATCH --job-name=busco_predict
#SBATCH --output=logs/busco_predict_%j.out
#SBATCH --error=logs/busco_predict_%j.err
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem=16G

set -euo pipefail

user_name=${USER}

CONTAINERS=/fs/scratch/PAS3260/${user_name}/Annotation/containers
ANNOT=/fs/scratch/PAS3260/${user_name}/Annotation
BUSCO_DL=${ANNOT}/busco_downloads

OUTDIR=${ANNOT}/02_funannotate/busco_predict
mkdir -p ${OUTDIR}
mkdir -p ${BUSCO_DL}
cd ${OUTDIR}

PROTEINS=$(ls ${ANNOT}/02_funannotate/predict_results/*.proteins.fa | head -1)
echo "Running BUSCO on: ${PROTEINS}"

apptainer exec \
  --bind ${ANNOT}:/data \
  --bind ${BUSCO_DL}:/busco_downloads \
  ${CONTAINERS}/busco_6.0.0.sif \
  busco \
    --in /data/02_funannotate/predict_results/$(basename ${PROTEINS}) \
    --out Pf_predict_proteins_busco \
    --mode proteins \
    --lineage_dataset dothideomycetes_odb12 \
    --download_path /busco_downloads \
    --cpu 12 \
    --force

echo "=== BUSCO of predicted proteins complete: $(date) ==="
cat ${OUTDIR}/Pf_predict_proteins_busco/short_summary*.txt
EOF

cd ${ANNOT}
sbatch ${ANNOT}/scripts/05c_busco_predict.sh
```

> **What to expect:** If training was successful and evidence was well-
> integrated, predicted protein BUSCO completeness should be **≥95%**
> for this well-studied genome.

---

## Step 6: Compare to the NCBI Reference Annotation

```bash
echo "=== Comparison: Funannotate2 vs. NCBI annotation ==="
echo ""
echo "Funannotate2 predicted gene count:"
grep -v "^#" ${ANNOT}/02_funannotate/predict_results/*.gff3 \
  | awk '$3 == "gene"' | wc -l

echo ""
echo "NCBI reference gene count:"
grep -v "^#" /fs/scratch/PAS3260/${user_name}/Peltaster/00_data/genome/Pf_annotation.gff3 \
  | awk '$3 == "gene"' | wc -l

echo ""
echo "Paper reports: 8,072 protein-coding genes"
```

---

## Biological Questions — Module 5

1. **Evidence weighting in EVM:** Funannotate2 assigns higher weight to
   RNA-seq transcript evidence than to ab initio predictions. Explain why
   this is biologically justified. Describe a specific scenario where
   over-weighting RNA-seq evidence could actually *reduce* annotation
   quality (hint: think about lowly expressed essential genes).

2. **Gene count discrepancy:** Your Funannotate2 prediction may differ from
   the NCBI-reported 8,072 genes by 200–500 models. List three reasons why
   two annotation pipelines applied to the same genome assembly can produce
   different gene counts, even when both use high-quality evidence.

3. **Interpreting BUSCO for predicted proteins:** Suppose your predicted
   protein BUSCO completeness is 94.2% (S=93.8%, D=0.4%, F=1.1%, M=4.7%).
   Compare this to the genome-mode BUSCO from the previous tutorial. What
   does the gap between genome completeness and proteins completeness tell
   you about your annotation?

4. **Short vs. long intron model:** *P. fructicola*'s 50 bp median intron
   means that the trained Augustus model has learned very different
   splice-site windows than models trained on most other fungi. If you
   applied your trained model to annotate a closely related but
   unsequenced Capnodiales species, would you expect better or worse
   performance than a generic *Aspergillus* model? What features would
   need to be similar for the transfer to work well?

5. **Locus tag design:** The `--strain LNHT1506` parameter sets the strain
   component of locus tags. NCBI requires globally unique locus tags for
   genome submissions. What information should a locus tag convey, and why
   does the inclusion of a strain identifier matter for comparative
   genomics databases?

---

## Take-Home Messages

> **Evidence Modeler is the intellectual core of gene prediction.** It is
> not a single predictor but a meta-predictor that translates your
> biological knowledge (which evidence is most reliable?) into mathematical
> weights that produce the best-supported consensus model.

> **BUSCO proteins completeness after prediction is the standard quality
> gate.** Do not proceed to functional annotation if BUSCO completeness is
> below ~85% — investigate whether training, repeat masking, or evidence
> inputs need adjustment.

> **The predicted GFF3 and GenBank files are primary research outputs.**
> Every downstream analysis depends on the accuracy of these gene models.

---

*Previous: [04 Training](04_train.md) | Next: [Module 6 → Functional Annotation](06_annotate.md)*
