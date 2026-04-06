# HCS 7004 — Genome Analytics
## Tutorial: Gene Prediction and Functional Annotation
### *Peltaster fructicola* LNHT1506 — From Assembly to Annotated Genome

---

## Overview

In the previous tutorial series you assembled the *Peltaster fructicola*
LNHT1506 genome, assessed it with BUSCO, and ran RNA-seq alignment and
quantification. You now have a polished, telomere-to-telomere (T2T) assembly
and a collection of short-read RNA-seq evidence. This tutorial takes you
through the next step: transforming a raw genome FASTA into a fully annotated
genome — with gene structures, protein sequences, and biological functions
assigned to every gene model.

We use **Funannotate2**, a comprehensive eukaryotic annotation pipeline that
integrates repeat masking, ab initio gene prediction, RNA-seq transcript
evidence, and several functional annotation databases into a single, coherent
workflow.

> **Working directory:** Each student works in their own directory:
> `/fs/scratch/PAS3260/${user_name}/Annotation`
> Replace `${user_name}` with your OSC username in every command and script.
> The first thing you do in Module 1 is set this variable.

---

## Tutorial Map

| Module | File | Topic | Est. time |
|--------|------|--------|----------|
| 0 | `00_overview.md` | Introduction & pipeline overview | — |
| 1 | `01_setup.md` | **Directory structure, containers, and shared resources** | 30 min |
| 2 | `02_rnaseq_evidence.md` | Downloading all RNA-seq evidence (SRP166999) | 30 min |
| 3 | `03_genome_clean.md` | Genome cleaning and repeat masking | 30 min |
| 4 | `04_train.md` | Training ab initio gene predictors | 60 min |
| 5 | `05_predict.md` | Gene prediction with Evidence Modeler | 90 min |
| 6 | `06_annotate.md` | Functional annotation and output inspection | 60 min |

> **Total active work:** ~5–6 hours + SLURM queue time.
> All bioinformatics tools run through **Apptainer containers**.
> The Funannotate2 container and its databases are pre-installed in a
> shared class directory. Students pull the remaining utility containers
> to their own space in Module 1.

---

## Learning Objectives

By the end of this tutorial you will be able to:

1. **Explain** the conceptual difference between structural annotation
   (gene prediction) and functional annotation
2. **Execute** each stage of the Funannotate2 pipeline on a real fungal genome
3. **Interpret** the role of RNA-seq, protein homology, and ab initio evidence
   in consensus gene model construction
4. **Critically evaluate** annotation quality using BUSCO and summary statistics
5. **Assign** functional annotations using Pfam, CAZyme, UniProt, EggNOG, and
   InterProScan databases
6. **Understand** why annotation decisions cascade downstream into every
   biological analysis that depends on gene models

---

## The Annotation Pipeline: Conceptual Overview

Gene annotation integrates three independent evidence streams into a
single consensus gene model set using **Evidence Modeler (EVM)**:

```
                    ┌─────────────────────────────────────────────────┐
                    │           Cleaned, repeat-masked genome          │
                    └──────────────────┬──────────────────────────────┘
                                       │
             ┌─────────────────────────┼─────────────────────────┐
             ▼                         ▼                         ▼
   Ab initio prediction       RNA-seq transcripts        Protein homology
   (Augustus, SNAP,           (all 15 conditions         (UniProt fungi,
    GlimmerHMM, GeneMark)      from SRP166999)            related species)
             └─────────────────────────┼─────────────────────────┘
                                       ▼
                            Evidence Modeler (EVM)
                          Weighted consensus gene models
                                       │
                                       ▼
                              BUSCO completeness check
                                       │
                                       ▼
                         ┌────────────────────────────┐
                         │   Functional annotation     │
                         │  Pfam · CAZymes · UniProt   │
                         │  MEROPS · InterProScan      │
                         │  EggNOG · GO terms · SignalP│
                         └────────────────────────────┘
```

Each evidence type compensates for the weaknesses of the others:

| Evidence type | Strength | Weakness |
|--------------|---------|---------|
| Ab initio (Augustus/SNAP/GlimmerHMM/GeneMark) | Finds novel genes with no homologs | Needs species-specific training |
| RNA-seq transcripts | Direct experimental evidence of expression | Misses unexpressed genes |
| Protein homology | Captures conserved genes accurately | Misses lineage-specific genes |

---

## Shared Resources for This Tutorial

The Funannotate2 container, annotation databases, GeneMark-ES, and SignalP 6
are pre-installed in a shared class directory and do **not** need to be
downloaded by individual students:

```
/fs/scratch/PAS3260/Team_Project/Containers/Funannotate2/
├── funannotate2.sif        ← Funannotate2 Apptainer container
├── databases/              ← Pfam, dbCAN, UniProt, MEROPS, gene names
├── augustus_config/        ← Augustus species models (read-only master copy)
├── gmes_linux_64_4/        ← GeneMark-ES (licensed, pre-installed)
└── signalp-6-package/      ← SignalP 6 (licensed, pre-installed)
```

Each student copies `augustus_config` to their own working directory
(required for the training step, which writes new species parameters).
All other shared resources are accessed read-only via bind mounts.

---

## Your Input Data

| Data | Source | Location |
|------|--------|----------|
| Genome assembly | GCA_001592805.2 | From previous tutorial |
| RNA-seq (15 runs) | SRP166999 (Wang et al. 2020) | Downloaded in Module 2 |

The RNA-seq dataset covers five growth conditions in triplicate:

| Condition | Replicates | SRR accessions |
|-----------|-----------|----------------|
| Potato dextrose broth, 5 days (PDB-5d) | 3 | SRR8119502–SRR8119504 |
| Potato dextrose broth, 15 days (PDB-15d) | 3 | SRR8119505–SRR8119507 |
| PDB + PEG 6000 osmotic stress, 15 days | 3 | SRR8119508–SRR8119510 |
| Potato dextrose agar (solid), 15 days | 3 | SRR8119511–SRR8119513 |
| Apple fruit inoculation, 15 days | 3 | SRR8119514–SRR8119516 |

> **Why use all 15 samples for annotation?** Gene models inferred from
> a single growth condition are biased toward genes expressed under that
> condition. Providing transcripts from all five conditions maximizes the
> number of gene models supported by direct RNA-seq evidence, reducing
> both false positives and false negatives in the final gene set.

---

## Pre-Tutorial Discussion Questions

Before starting Module 1, consider the following with your group:

1. What is the fundamental difference between a gene *prediction* and a gene
   *annotation*? At what point does a sequence feature become biologically
   interpretable?

2. *P. fructicola* has exceptionally compact genome architecture
   (median intron size 50 bp, 1.36 introns per gene). How would you expect
   these features to affect the performance of ab initio predictors trained
   on larger-intron organisms like *Aspergillus* or *Neurospora*?

3. You have RNA-seq evidence from five growth conditions. Would you expect
   the set of expressed genes to differ among conditions? What categories of
   genes might be condition-specific, and how would this affect the annotation
   if you used only a single condition?

---

*Proceed to:* **[Module 1 → Setup](01_setup.md)**
