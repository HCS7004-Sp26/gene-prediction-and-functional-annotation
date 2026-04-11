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
- **Secretome:** SignalP 6 (signal peptide prediction — baked into the container)
- **Conserved orthologs:** BUSCOlite (for completeness tracking)

Two additional tools extend the functional annotation further:

- **InterProScan:** the most comprehensive protein signature database,
  integrating 17 member databases (Pfam, NCBIfam, PANTHER, SUPERFAMILY,
  CDD, Gene3D, HAMAP, SMART, SFLD, and others) and assigning Gene Ontology
  (GO) terms. The shared installation at `${SHARED_IPS}` includes all
  databases pre-configured and ready to use.
- **EggNOG-mapper:** orthology-based functional annotation that transfers
  GO terms, KEGG pathways, COG categories, and gene names from characterized
  orthologs in the EggNOG5 database. The shared installation at
  `${SHARED_EGGNOG}` includes the container and the full reference database.

Together, these three annotation layers provide complementary and largely
independent evidence for gene function.

> **SignalP 6** model weights are baked into the `funannotate2.sif`
> container — no bind-mount is required. The tool predicts signal peptides
> and identifies secreted proteins, which is an important functional class
> for a fungal pathogen.

---

## Learning Objectives

- Run `funannotate2 annotate` to assign core functional annotations
- Run InterProScan for domain and GO term annotation against 17 databases
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
  --bind ${AUGUSTUS_CONFIG}:/opt/augustus/config \
  --env AUGUSTUS_CONFIG_PATH=/opt/augustus/config \
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

> **Note:** `funannotate2 annotate` does not call GeneMark, so the
> `gmes_linux_64_4` and `.gm_key` binds are not needed here.

---

## Step 2: Run InterProScan

InterProScan searches predicted proteins against 17 integrated protein
signature databases and assigns Gene Ontology (GO) terms. The container
and all databases are pre-installed in the shared class directory —
no setup is required. Expect **2–6 hours** for ~8,000 proteins.

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

# ── Paths ─────────────────────────────────────────────────────────────────────
user_name=Jonathan
ANNOT=/fs/scratch/PAS3260/${user_name}/Annotation
SHARED_IPS=/fs/scratch/PAS3260/Team_Project/Containers/InterProScan
SIF=${SHARED_IPS}/interproscan_5.77-108.0.sif
IPS_DATA=${SHARED_IPS}/data

# ── Output setup ──────────────────────────────────────────────────────────────
mkdir -p ${ANNOT}/03_iprscan

# Node-local temp — avoid overriding system $TMPDIR used by Java internally
JOB_TMP=/tmp/ips_${SLURM_JOB_ID}
mkdir -p ${JOB_TMP}

# ── Input ─────────────────────────────────────────────────────────────────────
PROTEINS=$(ls ${ANNOT}/02_funannotate/predict_results/*.proteins.fa | head -1)
echo "Input proteins: ${PROTEINS}"
echo "Container:      ${SIF}"
echo "IPS data:       ${IPS_DATA}"
echo "Job started:    $(date)"

# ── Applications — 17 databases confirmed present in shared data directory ────
# TMHMM excluded: the standard DTU binary is incompatible with IPS 5.77
# (IPS 5.77 requires a patched binary not distributed by DTU).
# Transmembrane topology is covered by Phobius, which is included.
APPS=AntiFam,CDD,FunFam,Gene3D,HAMAP,NCBIfam,PANTHER,Pfam,Phobius,PIRSF,PIRSR,PRINTS,ProSitePatterns,ProSiteProfiles,SFLD,SMART,SUPERFAMILY

# ── Run InterProScan ──────────────────────────────────────────────────────────
apptainer exec \
  --bind ${ANNOT}:/data \
  --bind ${IPS_DATA}:/opt/interproscan/data \
  --bind ${SHARED_IPS}/interproscan.properties:/opt/interproscan/interproscan.properties \
  ${SIF} \
  /opt/interproscan/interproscan.sh \
    -i /data/02_funannotate/predict_results/$(basename ${PROTEINS}) \
    -o /data/03_iprscan/Pf_interproscan.xml \
    -f XML \
    -appl ${APPS} \
    -goterms \
    -iprlookup \
    --cpu ${SLURM_NTASKS_PER_NODE} \
    --tempdir ${JOB_TMP} \
    -dp

echo "=== InterProScan complete: $(date) ==="
ls -lh ${ANNOT}/03_iprscan/Pf_interproscan.xml
EOF

cd ${ANNOT}
sbatch ${ANNOT}/scripts/06b_interproscan.sh
```

| Flag | Value in script | Purpose |
|------|----------------|---------|
| `-i` | `predict_results/*.proteins.fa` | Input protein FASTA file from Funannotate2 gene prediction |
| `-o` | `03_iprscan/Pf_interproscan.xml` | Output file path and base name |
| `-f` | `XML` | Output format; XML is required for Funannotate2 `annotate` integration |
| `-appl` | `AntiFam,CDD,FunFam,...` | Comma-separated list of member databases to search against |
| `-goterms` | *(flag only)* | Map matched signatures to Gene Ontology (GO) terms |
| `-iprlookup` | *(flag only)* | Add InterPro entry accessions and descriptions to output |
| `-dp` | *(flag only)* | Disable pre-calculated match lookup; forces local computation (recommended for non-model organisms) |
| `--cpu` | `${SLURM_NTASKS_PER_NODE}` | Number of parallel CPU cores to use; dynamically set from SLURM allocation |
| `--tempdir` | `/tmp/ips_${SLURM_JOB_ID}` | Directory for temporary files; node-local `/tmp` avoids shared filesystem I/O bottlenecks |

Databases used in this run:

| Database | Version | Type | What it detects |
|----------|---------|------|-----------------|
| **AntiFam** | 8.0 | HMM profiles | Spurious protein predictions; flags false-positive ORFs with homology to non-coding RNAs |
| **CDD** | 3.21 | RPS-BLAST profiles | NCBI conserved domain database; identifies functional domains and sites from multiple source databases (Pfam, SMART, COG, etc.) |
| **FunFam** | 4.3.0 | HMM profiles | Functional Families within CATH structural superfamilies; sub-classifies domains to predict specific molecular function |
| **Gene3D** | 4.3.0 | HMM profiles | Structural domains based on the CATH protein structure classification; assigns sequences to known 3D structural superfamilies |
| **HAMAP** | 2025_01 | HMM profiles | High-quality manually curated profiles for bacterial, archaeal, and plastid proteins; strong functional annotation in well-characterized families |
| **NCBIfam** | 18.0 | HMM profiles | NCBI-curated equivalents of TIGRFAMs; identifies protein families with precise functional or taxonomic specificity |
| **PANTHER** | 19.0 | HMM profiles | Large-scale family and subfamily classification linked to evolutionary history and GO annotations; broad phylogenetic coverage |
| **Pfam** | 38.1 | HMM profiles | The most widely used protein family database; covers domains, families, and repeats across all kingdoms of life |
| **Phobius** | 1.01 | Signal prediction | Combined transmembrane topology and signal peptide predictor; distinguishes genuine signal peptides from transmembrane helices |
| **PIRSF** | 3.10 | HMM profiles | Protein Information Resource SuperFamilies; whole-protein classification capturing evolutionary relationships at the family level |
| **PIRSR** | 2025_05 | HMM profiles | PIR Site Rules; detects functionally important residue-level features (active sites, binding sites) within PIRSF families |
| **PRINTS** | 42.0 | Fingerprints | Protein family fingerprints composed of conserved motif groups; useful for detecting distant homologs not captured by single-domain models |
| **ProSitePatterns** | 2025_01 | Regex patterns | Short, highly conserved sequence motifs linked to specific biochemical functions (e.g., active sites, post-translational modification sites) |
| **ProSiteProfiles** | 2025_01 | Scoring profiles | Generalized profiles for protein domains and families from PROSITE; more sensitive than patterns for detecting divergent members |
| **SFLD** | 4 | HMM profiles | Structure-Function Linkage Database; classifies enzymes by mechanistic and structural features, linking sequence to reaction chemistry |
| **SMART** | 9.0 | HMM profiles | Focuses on signaling, extracellular, and chromatin-associated domains; particularly valuable for eukaryotic domain architecture analysis |
| **SUPERFAMILY** | 1.75 | HMM profiles | Assigns sequences to SCOP structural superfamilies using hidden Markov models; provides evolutionary and structural context for domains |

---

## Step 3: Run EggNOG-mapper

EggNOG-mapper assigns functional annotations by finding the best ortholog
in the EggNOG5 database and transferring its GO terms, KEGG pathways, COG
categories, and gene name. Both the container and the reference database
are pre-installed in the shared class directory — no setup is required.

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

ANNOT=/fs/scratch/PAS3260/${user_name}/Annotation
SHARED_EGGNOG=/fs/scratch/PAS3260/Team_Project/Containers/eggNOG
EGGNOG_SIF=${SHARED_EGGNOG}/eggnog_mapper_2.1.13.sif
EGGNOG_DB=${SHARED_EGGNOG}/eggnog_db

mkdir -p ${ANNOT}/04_eggnog

PROTEINS=$(ls ${ANNOT}/02_funannotate/predict_results/*.proteins.fa | head -1)
echo "Input proteins: ${PROTEINS}"
echo "Container:      ${EGGNOG_SIF}"
echo "Database:       ${EGGNOG_DB}"
echo "Job started:    $(date)"

apptainer exec \
  --bind ${ANNOT}:/data \
  --bind ${EGGNOG_DB}:/eggnog_db \
  ${EGGNOG_SIF} \
  emapper.py \
    -i /data/02_funannotate/predict_results/$(basename ${PROTEINS}) \
    -o Pf_eggnog \
    --output_dir /data/04_eggnog \
    --data_dir /eggnog_db \
    --tax_scope Fungi \
    --go_evidence all \
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
back using `funannotate2-addons` (`f2a`). The `--parse` flag tells each
`f2a` command to skip re-running the tool and go straight to parsing the
pre-computed output file.

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

user_name=Jonathan
ANNOT=/fs/scratch/PAS3260/${user_name}/Annotation
SHARED_F2=/fs/scratch/PAS3260/Team_Project/Containers/Funannotate2
F2_CONTAINER=${SHARED_F2}/funannotate2.sif
F2_DB=${SHARED_F2}/databases
AUGUSTUS_CONFIG=${ANNOT}/augustus_config

echo "=== Step 1: Parse InterProScan XML results: $(date) ==="
apptainer exec \
  --bind ${ANNOT}:/data \
  --bind ${F2_DB}:/f2_db \
  --env FUNANNOTATE2_DB=/f2_db \
  ${F2_CONTAINER} \
  f2a iprscan \
    -i /data/02_funannotate \
    --parse /data/03_iprscan/Pf_interproscan.xml \
    --cpus 4

echo "=== Step 2: Parse EggNOG-mapper results: $(date) ==="
apptainer exec \
  --bind ${ANNOT}:/data \
  --bind ${F2_DB}:/f2_db \
  --env FUNANNOTATE2_DB=/f2_db \
  ${F2_CONTAINER} \
  f2a emapper \
    -i /data/02_funannotate \
    --parse /data/04_eggnog/Pf_eggnog.emapper.annotations \
    --cpus 4

echo "=== Step 3: Re-run funannotate2 annotate with all evidence: $(date) ==="
apptainer exec \
  --bind ${ANNOT}:/data \
  --bind ${F2_DB}:/f2_db \
  --bind ${AUGUSTUS_CONFIG}:/opt/augustus/config \
  --env AUGUSTUS_CONFIG_PATH=/opt/augustus/config \
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

SignalP 6 (baked into the shared container) predicts signal peptides,
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