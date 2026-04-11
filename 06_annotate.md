# Module 6 — Functional Annotation
## Assigning Biological Meaning to Gene Models

---

## Background

Structural annotation tells you *where* genes are. Functional annotation
tells you *what they do*. Funannotate2's `annotate` command searches the
predicted protein sequences against multiple curated databases to assign:

- **Protein domains:** Pfam-A (structural and functional domains)
- **Enzyme families:** dbCAN (CAZymes — carbohydrate-active enzymes)
- **Characterized proteins:** UniProtKB/Swiss-Prot (manually reviewed)
- **Proteases:** MEROPS (protease families and clans)
- **Secretome:** SignalP 6 (signal peptide prediction — available in this installation)
- **Conserved orthologs:** BUSCOlite (for completeness tracking)

Two additional tools extend the functional annotation further:

- **InterProScan:** the most comprehensive protein signature database,
  integrating 13 member databases (Pfam, TIGRFAM, PANTHER, SUPERFAMILY,
  CDD, etc.) and assigning Gene Ontology (GO) terms
- **EggNOG-mapper:** orthology-based functional annotation that transfers
  GO terms, KEGG pathways, COG categories, and gene names from characterized
  orthologs in the EggNOG5 database

Together, these three annotation layers provide complementary and largely
independent evidence for gene function.

> **SignalP 6 is included in this class installation** at
> `${SHARED_F2}/signalp-6-package/`. This tool predicts signal peptides
> and identifies secreted proteins — an important functional class for a
> fungal pathogen. It is bound into the container automatically using the
> standard bind pattern.

---

## Learning Objectives

- Run `funannotate2 annotate` to assign core functional annotations
- Run InterProScan for domain and GO term annotation
- Run EggNOG-mapper for orthology-based functional transfer
- Integrate InterProScan and EggNOG results into the Funannotate2 output
- Inspect the final annotated gene set and interpret functional categories
- Explore *P. fructicola*'s CAZyme repertoire as a biological case study

---

## Step 1: Run funannotate2 annotate

This command runs the core annotation databases (Pfam, dbCAN, UniProt,
MEROPS, SignalP, BUSCOlite) against the predicted proteins.

```bash
cat > ${ANNOT}/scripts/06a_f2_annotate.sh << 'EOF'
#!/bin/bash
#SBATCH --account=PAS3260
#SBATCH --job-name=f2_annotate
#SBATCH --output=logs/f2_annotate_%j.out
#SBATCH --error=logs/f2_annotate_%j.err
#SBATCH --time=04:00:00
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

apptainer exec \
  --bind ${ANNOT}:/data \
  --bind ${F2_DB}:/f2_db \
  --bind ${AUGUSTUS_CONFIG}:/opt/augustus_config \
  --bind ${SHARED_F2}/gmes_linux_64_4:/gmes_linux_64_4 \
  --bind ${SHARED_F2}/signalp-6-package:/signalp-6-package \
  --env AUGUSTUS_CONFIG_PATH=/opt/augustus_config \
  --env FUNANNOTATE2_DB=/f2_db \
  ${F2_CONTAINER} \
  funannotate2 annotate \
    -i /data/02_funannotate \
    -s "Peltaster fructicola" \
    -st LNHT1506 \
    --busco-lineage dothideomycetes_odb12 \
    --cpus 12

echo "=== funannotate2 annotate complete: $(date) ==="
ls -lh ${ANNOT}/02_funannotate/annotate_results/
EOF

cd ${ANNOT}
sbatch ${ANNOT}/scripts/06a_f2_annotate.sh
```

---

## Step 2: Run InterProScan

InterProScan searches predicted proteins against 13 integrated protein
signature databases and assigns Gene Ontology (GO) terms. Expect
**2–6 hours** for ~8,000 proteins.

```bash
cat > ${ANNOT}/scripts/06b_interproscan.sh << 'EOF'
#!/bin/bash
#SBATCH --account=PAS3260
#SBATCH --job-name=interproscan
#SBATCH --output=logs/interproscan_%j.out
#SBATCH --error=logs/interproscan_%j.err
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=64G

set -euo pipefail

user_name=Jonathan

CONTAINERS=/fs/scratch/PAS3260/${user_name}/Annotation/containers
ANNOT=/fs/scratch/PAS3260/${user_name}/Annotation

mkdir -p ${ANNOT}/03_iprscan

PROTEINS=$(ls ${ANNOT}/02_funannotate/predict_results/*.proteins.fa | head -1)
echo "Input proteins: ${PROTEINS}"

apptainer exec \
  --bind ${ANNOT}:/data \
  ${CONTAINERS}/interproscan_5.73.sif \
  interproscan.sh \
    -i /data/02_funannotate/predict_results/$(basename ${PROTEINS}) \
    -o /data/03_iprscan/Pf_interproscan.xml \
    -f XML \
    -appl Pfam,TIGRFAM,PANTHER,CDD,SUPERFAMILY,Gene3D \
    -goterms \
    -iprlookup \
    --cpu 16 \
    -dp

echo "=== InterProScan complete: $(date) ==="
ls -lh ${ANNOT}/03_iprscan/Pf_interproscan.xml
EOF

cd ${ANNOT}
sbatch ${ANNOT}/scripts/06b_interproscan.sh
```

**Key InterProScan flags:**

| Flag | Purpose |
|------|---------|
| `-appl` | Applications to run (comma-separated) |
| `-goterms` | Map matched signatures to GO terms |
| `-iprlookup` | Add InterPro entry descriptions |
| `-dp` | Disable pre-calculated match lookup (recommended for novel organisms) |
| `-f XML` | Output format required for Funannotate2 integration |

---

## Step 3: Run EggNOG-mapper

```bash
cat > ${ANNOT}/scripts/06c_eggnog.sh << 'EOF'
#!/bin/bash
#SBATCH --account=PAS3260
#SBATCH --job-name=eggnog_mapper
#SBATCH --output=logs/eggnog_mapper_%j.out
#SBATCH --error=logs/eggnog_mapper_%j.err
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=64G

set -euo pipefail

user_name=Jonathan

CONTAINERS=/fs/scratch/PAS3260/${user_name}/Annotation/containers
ANNOT=/fs/scratch/PAS3260/${user_name}/Annotation
EGGNOG_DB=/fs/scratch/PAS3260/Team_Project/Containers/Funannotate2/eggnog_db

mkdir -p ${ANNOT}/04_eggnog

PROTEINS=$(ls ${ANNOT}/02_funannotate/predict_results/*.proteins.fa | head -1)
echo "Input proteins: ${PROTEINS}"

apptainer exec \
  --bind ${ANNOT}:/data \
  --bind ${EGGNOG_DB}:/eggnog_db \
  ${CONTAINERS}/eggnog_mapper_2.1.13.sif \
  emapper.py \
    -i /data/02_funannotate/predict_results/$(basename ${PROTEINS}) \
    -o Pf_eggnog \
    --output_dir /data/04_eggnog \
    --data_dir /eggnog_db \
    --tax_scope Fungi \
    --go_evidence non_electronic \
    --target_orthologs all \
    --seed_ortholog_evalue 0.001 \
    --seed_ortholog_score 60 \
    --cpu 16 \
    --override

echo "=== EggNOG-mapper complete: $(date) ==="
ls -lh ${ANNOT}/04_eggnog/Pf_eggnog.*
EOF

cd ${ANNOT}
sbatch ${ANNOT}/scripts/06c_eggnog.sh
```

---

## Step 4: Integrate InterProScan and EggNOG into Funannotate2

After both external annotation tools finish, integrate their results
back using `funannotate2-addons` (`f2a`):

```bash
cat > ${ANNOT}/scripts/06d_integrate_annotations.sh << 'EOF'
#!/bin/bash
#SBATCH --account=PAS3260
#SBATCH --job-name=f2_integrate
#SBATCH --output=logs/f2_integrate_%j.out
#SBATCH --error=logs/f2_integrate_%j.err
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=16G

set -euo pipefail

user_name=${USER}

ANNOT=/fs/scratch/PAS3260/${user_name}/Annotation
SHARED_F2=/fs/scratch/PAS3260/Team_Project/Containers/Funannotate2
F2_CONTAINER=${SHARED_F2}/funannotate2.sif
F2_DB=${SHARED_F2}/databases
AUGUSTUS_CONFIG=${ANNOT}/augustus_config

# Parse InterProScan XML results
apptainer exec \
  --bind ${ANNOT}:/data \
  --bind ${F2_DB}:/f2_db \
  --env FUNANNOTATE2_DB=/f2_db \
  ${F2_CONTAINER} \
  f2a iprscan \
    -i /data/02_funannotate \
    --iprscan_xml /data/03_iprscan/Pf_interproscan.xml \
    --cpus 4

# Parse EggNOG results
apptainer exec \
  --bind ${ANNOT}:/data \
  --bind ${F2_DB}:/f2_db \
  --env FUNANNOTATE2_DB=/f2_db \
  ${F2_CONTAINER} \
  f2a emapper \
    -i /data/02_funannotate \
    --emapper_annotations /data/04_eggnog/Pf_eggnog.emapper.annotations \
    --cpus 4

# Re-run funannotate2 annotate incorporating all external annotations
apptainer exec \
  --bind ${ANNOT}:/data \
  --bind ${F2_DB}:/f2_db \
  --bind ${AUGUSTUS_CONFIG}:/opt/augustus_config \
  --bind ${SHARED_F2}/gmes_linux_64_4:/gmes_linux_64_4 \
  --bind ${SHARED_F2}/signalp-6-package:/signalp-6-package \
  --env AUGUSTUS_CONFIG_PATH=/opt/augustus_config \
  --env FUNANNOTATE2_DB=/f2_db \
  ${F2_CONTAINER} \
  funannotate2 annotate \
    -i /data/02_funannotate \
    --busco-lineage dothideomycetes_odb12 \
    --cpus 4

echo "=== Integration complete: $(date) ==="
ls -lh ${ANNOT}/02_funannotate/annotate_results/
EOF

cd ${ANNOT}
sbatch ${ANNOT}/scripts/06d_integrate_annotations.sh
```

---

## Step 5: Inspect the Final Annotated Output

```bash
ANNOTATE_DIR=${ANNOT}/02_funannotate/annotate_results

ls -lh ${ANNOTATE_DIR}/

echo "=== Output file types ==="
echo "*.gff3          — GFF3 with functional annotations in attributes"
echo "*.gbk           — GenBank format (NCBI submission or visualization)"
echo "*.proteins.fa   — Protein FASTA with product descriptions in headers"
echo "*.annotations.txt — Tab-delimited annotation table"
echo ""

echo "=== Annotation coverage ==="
TOTAL=$(grep -v "^#" ${ANNOTATE_DIR}/*.gff3 | awk '$3 == "gene"' | wc -l)
echo "Total gene models: ${TOTAL}"

echo ""
echo "Genes with GO terms:"
grep -v "^#" ${ANNOTATE_DIR}/*.gff3 | \
  awk '$3 == "mRNA" && $9 ~ /Ontology_term/' | wc -l

echo ""
echo "Genes with Pfam domains:"
grep -v "^#" ${ANNOTATE_DIR}/*.gff3 | \
  awk '$3 == "mRNA" && $9 ~ /Dbxref.*Pfam/' | wc -l

echo ""
echo "Genes with SignalP signal peptides (secreted proteins):"
grep -v "^#" ${ANNOTATE_DIR}/*.gff3 | \
  awk '$3 == "mRNA" && $9 ~ /SignalP/' | wc -l

echo ""
echo "Genes with named products (not hypothetical):"
grep -v "^#" ${ANNOTATE_DIR}/*.gff3 | \
  awk '$3 == "mRNA" && $9 ~ /product=/ && $9 !~ /hypothetical/' | wc -l
```

---

## Step 6: Explore the CAZyme Repertoire

One of the key findings of Wang et al. (2020) is that *P. fructicola* has
a reduced carbohydrate-active enzyme (CAZyme) repertoire compared to
plant-pathogenic relatives. Let us reproduce that analysis:

```bash
ANNOTATE_DIR=${ANNOT}/02_funannotate/annotate_results

echo "=== CAZyme annotation summary ==="
echo ""
grep -v "^#" ${ANNOTATE_DIR}/*.gff3 | \
  awk '$3 == "mRNA" && $9 ~ /CAZy/' | \
  grep -oP 'CAZy=[^;]+' | \
  sed 's/CAZy=//' | \
  tr ',' '\n' | \
  grep -v "^$" | \
  sort | uniq -c | sort -rn | \
  awk '{printf "%-8s %s\n", $1, $2}' | head -30

echo ""
echo "Summary by CAZyme class:"
grep -v "^#" ${ANNOTATE_DIR}/*.gff3 | \
  awk '$3 == "mRNA" && $9 ~ /CAZy/' | \
  grep -oP 'CAZy=[^;]+' | \
  sed 's/CAZy=//' | \
  tr ',' '\n' | \
  grep -oP '^[A-Z]+' | \
  sort | uniq -c | sort -rn | \
  awk '{printf "  %-6s: %s families\n", $2, $1}'
```

**CAZyme class codes:**

| Code | Enzyme class |
|------|-------------|
| GH | Glycoside hydrolases (cellulases, amylases, etc.) |
| GT | Glycosyltransferases |
| PL | Polysaccharide lyases (pectate lyases, etc.) |
| CE | Carbohydrate esterases |
| AA | Auxiliary activities (oxidoreductases) |
| CBM | Carbohydrate-binding modules |

> **Biological interpretation:** The paper reports that *P. fructicola*
> has only ~15 CAZymes total, compared to 250–400 in aggressive plant cell
> wall-degrading pathogens. Your annotation should confirm this finding.

---

## Step 7: Explore the Secretome

SignalP 6 (pre-installed in the shared container) predicts signal peptides,
identifying proteins likely to be secreted — a key functional class for
a surface-colonizing fungus:

```bash
echo "=== Secretome summary ==="

# Count predicted secreted proteins
grep -v "^#" ${ANNOTATE_DIR}/*.gff3 | \
  awk '$3 == "mRNA" && $9 ~ /SignalP/' | wc -l
echo "proteins with predicted signal peptides"

echo ""
# What functions do secreted proteins have?
echo "Top functional categories in predicted secretome:"
grep -v "^#" ${ANNOTATE_DIR}/*.gff3 | \
  awk '$3 == "mRNA" && $9 ~ /SignalP/ && $9 ~ /product=/' | \
  grep -oP 'product=[^;]+' | \
  sed 's/product=//' | \
  sort | uniq -c | sort -rn | \
  grep -v "hypothetical" | head -20
```

---

## Step 8: Generate an Annotation Summary Table

```bash
grep -v "^#" ${ANNOTATE_DIR}/*.gff3 | \
  awk '$3 == "mRNA"' | \
  awk '{
    match($9, /ID=([^;]+)/, id)
    match($9, /product=([^;]+)/, prod)
    match($9, /Ontology_term=([^;]+)/, go)
    match($9, /SignalP/, sp)
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
      id[1], $1, $4, $5,
      (prod[1] != "") ? prod[1] : "hypothetical protein",
      (go[1]  != "") ? go[1]  : "NA",
      (sp[0]  != "") ? "yes"  : "no"
  }' | \
  sort -k1,1 > ${ANNOTATE_DIR}/Pf_gene_annotation_table.tsv

echo "Annotation table created:"
echo "Columns: GeneID  Chromosome  Start  End  Product  GO_terms  SignalP"
head -5 ${ANNOTATE_DIR}/Pf_gene_annotation_table.tsv
echo ""
wc -l ${ANNOTATE_DIR}/Pf_gene_annotation_table.tsv
```

---

## Biological Questions — Module 6

1. **Annotation depth by evidence type:** You ran three independent
   annotation tools: Funannotate2 core databases, InterProScan, and
   EggNOG-mapper. For each of the following gene types, which tool would
   provide the most informative annotation and why?
   - a) A cutinase gene retained from a plant-pathogenic ancestor
   - b) A novel, *P. fructicola*-specific gene with no BLAST hits
   - c) A core metabolic enzyme (e.g., pyruvate kinase)
   - d) A CAZyme with a known substrate specificity

2. **Hypothetical proteins:** What percentage of gene models are annotated
   as "hypothetical protein" (no functional name assigned)? Is this higher
   or lower than you expected for a well-characterized fungal genome? What
   three approaches would you take to reduce the hypothetical protein fraction?

3. **CAZyme comparison:** The paper reports that *P. fructicola* has ~15
   total CAZymes vs. ~400 in an aggressive pathogen. Does your annotation
   recover a similar number? If your count differs substantially from the
   published value, identify which annotation step (prediction, Pfam, dbCAN)
   is most likely responsible for the discrepancy.

4. **Secretome and lifestyle:** Using your SignalP predictions, calculate
   what percentage of the predicted proteome is secreted. Compare this to
   the known biology of *P. fructicola* as a surface colonist that dissolves
   epicuticular wax but does not penetrate plant cells. What substrate
   categories would you expect to be enriched among the secreted proteins?

5. **NCBI submission readiness:** The GenBank output file (`.gbk`) from
   `funannotate2 annotate` is formatted for NCBI submission. What metadata
   would you need to add before submitting this annotation to the NCBI
   GenBank database? What is the role of the `--strain` parameter in locus
   tag uniqueness for public database submissions?

---

## Take-Home Messages

> **Functional annotation is an inference, not a measurement.** Every
> gene name and GO term in your annotation is a prediction based on
> sequence similarity to characterized proteins. The confidence varies
> from near-certain (ortholog of a well-studied enzyme with 95% identity)
> to speculative (remote Pfam domain match).

> **Multiple independent annotation databases provide complementary
> evidence.** A gene annotated consistently by Pfam, InterProScan, EggNOG,
> and UniProt is high-confidence. A gene annotated by only one database
> should be treated with more caution.

> **The CAZyme analysis and secretome are windows into ecological
> adaptation.** The dramatic reduction in plant cell wall-degrading enzymes,
> confirmed by your annotation, is the molecular signature of
> *P. fructicola*'s evolutionary transition from pathogen to ectophyte.

---

## Appendix: Final Output Files and Their Uses

| File | Format | Use for |
|------|--------|---------|
| `*.gff3` | GFF3 | Genome browser visualization, differential expression |
| `*.gbk` | GenBank | NCBI submission, Artemis visualization |
| `*.proteins.fa` | FASTA | Comparative genomics, OrthoFinder, BLAST |
| `*.transcripts.fa` | FASTA | DESeq2 reference, isoform analysis |
| `*.annotations.txt` | TSV | Pathway analysis, GO enrichment |
| `Pf_gene_annotation_table.tsv` | TSV | Custom analyses, filtering by function |

---

*Previous: [05 Gene Prediction](05_predict.md) | Back to: [00 Overview](00_overview.md)*
