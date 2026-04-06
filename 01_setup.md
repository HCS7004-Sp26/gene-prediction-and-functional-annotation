# Module 1 — Environment Setup
## Directory Structure, Containers, and Shared Resources

---

## Background

Genome annotation is computationally intensive and generates a large,
interconnected set of intermediate files. A well-organized directory structure
is not optional — Funannotate2 expects specific input locations and writes
output to a consistent folder hierarchy. Before any analysis, we establish
that structure, pull the utility containers each student needs, and configure
access to the shared Funannotate2 resources pre-installed for this class.

> **Relationship to the previous tutorial:** The genome assembly lives at
> `/fs/scratch/PAS3260/${user_name}/Peltaster/00_data/genome/`.
> This annotation tutorial uses a dedicated working directory so the two
> projects remain independent and reproducible.

---

## Learning Objectives

- Set the `user_name` variable that customizes every path in this tutorial
- Create the annotation project directory tree
- Pull utility containers (SRA Toolkit, Trimmomatic, FastQC, InterProScan,
  EggNOG-mapper) to your personal containers directory
- Copy the shared `augustus_config` to your working directory
- Verify all containers and shared resources are accessible
- Set persistent environment variables for the project

---

## Step 1: Set Your Username Variable

> **Do this first, every session.** All paths in this tutorial use
> `${user_name}` to point to your personal working directory.

```bash
# Replace "myusername" with your actual OSC username
export user_name=myusername

# Verify
echo "user_name = ${user_name}"
echo "Your working root = /fs/scratch/PAS3260/${user_name}/Annotation"
```

---

## Step 2: Define Path Variables

```bash
# Personal working directory and containers
export ANNOT=/fs/scratch/PAS3260/${user_name}/Annotation
export CONTAINERS=/fs/scratch/PAS3260/${user_name}/Annotation/containers

# Shared Funannotate2 resources (pre-installed, read-only)
export SHARED_F2=/fs/scratch/PAS3260/Team_Project/Containers/Funannotate2
export F2_CONTAINER=${SHARED_F2}/funannotate2.sif
export F2_DB=${SHARED_F2}/databases

# Your own writable augustus_config (copied from shared in Step 5)
export AUGUSTUS_CONFIG=${ANNOT}/augustus_config

# EggNOG-mapper database (downloaded per-student in Step 7)
export EGGNOG_DB=${ANNOT}/eggnog_db
```

Confirm the shared resources are accessible:

```bash
ls -lh ${SHARED_F2}/
# Expected output:
# funannotate2.sif
# databases/
# augustus_config/
# gmes_linux_64_4/
# signalp-6-package/
```

---

## Step 3: Create the Project Directory Structure

```bash
mkdir -p \
  ${ANNOT}/00_genome \
  ${ANNOT}/01_rnaseq/raw \
  ${ANNOT}/01_rnaseq/trimmed \
  ${ANNOT}/01_rnaseq/qc \
  ${ANNOT}/02_funannotate \
  ${ANNOT}/03_iprscan \
  ${ANNOT}/04_eggnog \
  ${ANNOT}/logs \
  ${ANNOT}/scripts \
  ${CONTAINERS}

tree -d ${ANNOT}
```

### Annotated directory tree

```
/fs/scratch/PAS3260/${user_name}/Annotation/
├── 00_genome/           ← cleaned, repeat-masked FASTA
├── 01_rnaseq/
│   ├── raw/             ← downloaded FASTQ files (SRP166999)
│   ├── trimmed/         ← adapter-trimmed paired FASTQ
│   └── qc/              ← FastQC and MultiQC reports
├── 02_funannotate/      ← Funannotate2 output (train, predict, annotate)
├── 03_iprscan/          ← InterProScan output and parsed results
├── 04_eggnog/           ← EggNOG-mapper output and parsed results
├── augustus_config/     ← your writable copy of the Augustus config
├── containers/          ← Apptainer .sif images (utility tools)
├── eggnog_db/           ← EggNOG-mapper reference database
├── logs/                ← SLURM stdout/stderr for every job
└── scripts/             ← all SLURM batch scripts
```

---

## Step 4: Copy the Genome Assembly

```bash
# Verify the source genome exists
ls -lh /fs/scratch/PAS3260/${user_name}/Peltaster/00_data/genome/Peltaster_fructicola_genome.fasta

# Copy to the annotation project directory
cp /fs/scratch/PAS3260/${user_name}/Peltaster/00_data/genome/Peltaster_fructicola_genome.fasta \
   ${ANNOT}/00_genome/Pf_assembly_raw.fasta

echo "Genome copied."
grep -c "^>" ${ANNOT}/00_genome/Pf_assembly_raw.fasta
```

---

## Step 5: Copy Your Personal augustus_config

Augustus writes trained species parameters into a subdirectory of
`$AUGUSTUS_CONFIG_PATH`. Because the shared copy is read-only, each student
needs their own writable copy. This is a one-time step:

```bash
cp -r ${SHARED_F2}/augustus_config ${ANNOT}/augustus_config

echo "augustus_config copied."
ls ${ANNOT}/augustus_config/
# Should contain: config/  extrinsic/  model/  profile/  species/
```

> **Why this matters:** During `funannotate2 train`, Augustus saves the
> trained *P. fructicola* parameters to
> `${AUGUSTUS_CONFIG}/species/peltaster_fructicola/`. Without a writable
> copy, training fails immediately.

---

## Step 6: Persist Environment Variables

Write all path variables to `~/.bashrc` so they survive log-out.
We write literal paths (not variable references) to avoid
expansion-order issues at login.

```bash
cat >> ~/.bashrc << EOF

# ---- HCS 7004 Annotation tutorial ----
export user_name=${user_name}
export ANNOT=/fs/scratch/PAS3260/${user_name}/Annotation
export CONTAINERS=/fs/scratch/PAS3260/${user_name}/Annotation/containers
export SHARED_F2=/fs/scratch/PAS3260/Team_Project/Containers/Funannotate2
export F2_CONTAINER=\${SHARED_F2}/funannotate2.sif
export F2_DB=\${SHARED_F2}/databases
export AUGUSTUS_CONFIG=/fs/scratch/PAS3260/${user_name}/Annotation/augustus_config
export EGGNOG_DB=/fs/scratch/PAS3260/${user_name}/Annotation/eggnog_db
EOF

source ~/.bashrc

# Verify
echo "user_name      = ${user_name}"
echo "ANNOT          = ${ANNOT}"
echo "F2_CONTAINER   = ${F2_CONTAINER}"
echo "F2_DB          = ${F2_DB}"
echo "AUGUSTUS_CONFIG= ${AUGUSTUS_CONFIG}"
```

---

## Step 7: Pull Utility Containers

The Funannotate2 container is pre-built and shared. You only need to pull
the five utility containers used for data download, read QC, and functional
annotation. Always `cd ${ANNOT}` before submitting so SLURM logs land
in `${ANNOT}/logs/`.

```bash
cat > ${ANNOT}/scripts/01_pull_containers.sh << 'EOF'
#!/bin/bash
#SBATCH --account=PAS3260
#SBATCH --job-name=pull_containers
#SBATCH --output=logs/pull_containers_%j.out
#SBATCH --error=logs/pull_containers_%j.err
#SBATCH --time=01:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=8G

set -euo pipefail

CONTAINERS=/fs/scratch/PAS3260/${user_name}/Annotation/containers

echo "=== Pulling utility containers ==="
echo "Started: $(date)"

# ---- sra-tools (RNA-seq download) ----
apptainer pull \
  ${CONTAINERS}/sratools_3.2.1.sif \
  oras://community.wave.seqera.io/library/sra-tools:3.2.1--846898724ee33c64

# ---- Trimmomatic (adapter trimming) ----
apptainer pull \
  ${CONTAINERS}/trimmomatic_0.40.sif \
  oras://community.wave.seqera.io/library/trimmomatic:0.40--7b5b7590373e6fc4

# ---- FastQC (read QC) ----
apptainer pull \
  ${CONTAINERS}/fastqc_0.12.1.sif \
  oras://community.wave.seqera.io/library/fastqc:0.12.1--104d26ddd9519960

# ---- InterProScan (protein domain annotation) ----
apptainer pull \
  ${CONTAINERS}/interproscan_5.72.sif \
  oras://community.wave.seqera.io/library/interproscan:5.72-103.0--hec16e2b_0

# ---- EggNOG-mapper (orthology-based functional annotation) ----
apptainer pull \
  ${CONTAINERS}/eggnog_mapper_2.1.12.sif \
  oras://community.wave.seqera.io/library/eggnog-mapper:2.1.12--pyhdfd78af_0

echo ""
echo "=== All utility containers pulled ==="
echo "Finished: $(date)"
ls -lh ${CONTAINERS}/
EOF

cd ${ANNOT}
sbatch ${ANNOT}/scripts/01_pull_containers.sh
squeue -u ${USER}
```

> **Note:** Container tags change as new builds are released. If any
> `oras://` URI returns a 404, visit <https://seqera.io/containers/> and
> search for the tool name to find the current tag.

---

## Step 8: Verify All Containers

Once the pull job completes, verify each container:

```bash
echo "=== Shared Funannotate2 container ==="
echo -n "funannotate2:   "
apptainer exec \
  --bind ${F2_DB}:/f2_db \
  --bind ${AUGUSTUS_CONFIG}:/opt/augustus_config \
  --env AUGUSTUS_CONFIG_PATH=/opt/augustus_config \
  ${F2_CONTAINER} \
  funannotate2 --version 2>&1 | head -1

echo ""
echo "=== Utility containers ==="
echo -n "sra-tools:      "
apptainer exec ${CONTAINERS}/sratools_3.2.1.sif \
  fastq-dump --version 2>&1 | grep "fastq-dump" | head -1

echo -n "trimmomatic:    "
apptainer exec ${CONTAINERS}/trimmomatic_0.40.sif \
  trimmomatic -version 2>&1 | head -1

echo -n "fastqc:         "
apptainer exec ${CONTAINERS}/fastqc_0.12.1.sif \
  fastqc --version 2>&1 | head -1

echo -n "interproscan:   "
apptainer exec ${CONTAINERS}/interproscan_5.72.sif \
  interproscan.sh --version 2>&1 | grep "InterProScan" | head -1

echo -n "eggnog-mapper:  "
apptainer exec ${CONTAINERS}/eggnog_mapper_2.1.12.sif \
  emapper.py --version 2>&1 | head -1

echo ""
echo "Verification complete."
```

---

## Step 9: Set Up the EggNOG-mapper Database

The EggNOG-mapper reference database (~9 GB) is not shared because each
student needs read-write access during mapping. Download it once to your
own space:

```bash
mkdir -p ${EGGNOG_DB}

cat > ${ANNOT}/scripts/01b_setup_eggnog.sh << 'EOF'
#!/bin/bash
#SBATCH --account=PAS3260
#SBATCH --job-name=eggnog_db
#SBATCH --output=logs/eggnog_db_%j.out
#SBATCH --error=logs/eggnog_db_%j.err
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --mem=16G

set -euo pipefail

CONTAINERS=/fs/scratch/PAS3260/${user_name}/Annotation/containers
EGGNOG_DB=/fs/scratch/PAS3260/${user_name}/Annotation/eggnog_db

apptainer exec \
  --bind ${EGGNOG_DB}:/eggnog_db \
  ${CONTAINERS}/eggnog_mapper_2.1.12.sif \
  download_eggnog_data.py \
    --data_dir /eggnog_db \
    -y

echo "EggNOG database download complete."
ls -lh ${EGGNOG_DB}/
EOF

cd ${ANNOT}
sbatch ${ANNOT}/scripts/01b_setup_eggnog.sh
```

---

## How Funannotate2 Is Run in This Tutorial

Every Funannotate2 command uses the following bind-mount pattern so the
container can access the shared databases, your writable `augustus_config`,
and the licensed tools (GeneMark-ES, SignalP 6):

```bash
apptainer exec \
  --bind ${ANNOT}:/data \
  --bind ${F2_DB}:/f2_db \
  --bind ${AUGUSTUS_CONFIG}:/opt/augustus_config \
  --bind ${SHARED_F2}/gmes_linux_64_4:/gmes_linux_64_4 \
  --bind ${SHARED_F2}/signalp-6-package:/signalp-6-package \
  --env AUGUSTUS_CONFIG_PATH=/opt/augustus_config \
  ${F2_CONTAINER} \
  funannotate2 <command> ...
```

| Bind mount | Purpose |
|-----------|---------|
| `${ANNOT}:/data` | Your project data, accessible as `/data` inside container |
| `${F2_DB}:/f2_db` | Shared annotation databases (Pfam, dbCAN, UniProt, MEROPS) |
| `${AUGUSTUS_CONFIG}:/opt/augustus_config` | Your writable Augustus config |
| `gmes_linux_64_4:/gmes_linux_64_4` | GeneMark-ES binary (requires license) |
| `signalp-6-package:/signalp-6-package` | SignalP 6 for secretome prediction |

> This pattern is repeated verbatim in every Funannotate2 script throughout
> the tutorial. All input/output paths inside container commands start
> with `/data/`, `/f2_db/`, etc.

---

## SLURM Header Template for This Tutorial

Every script uses the following header. Always `cd ${ANNOT}` before
calling `sbatch` so SLURM logs land in `${ANNOT}/logs/`.

```bash
#!/bin/bash
#SBATCH --account=PAS3260
#SBATCH --job-name=descriptive_name
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --time=HH:MM:SS
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=N
#SBATCH --mem=XG

# ---- Path variables (hardcoded — do not use $user_name inside scripts) ----
ANNOT=/fs/scratch/PAS3260/${user_name}/Annotation
CONTAINERS=/fs/scratch/PAS3260/${user_name}/Annotation/containers
SHARED_F2=/fs/scratch/PAS3260/Team_Project/Containers/Funannotate2
F2_CONTAINER=${SHARED_F2}/funannotate2.sif
F2_DB=${SHARED_F2}/databases
AUGUSTUS_CONFIG=${ANNOT}/augustus_config
EGGNOG_DB=${ANNOT}/eggnog_db

set -euo pipefail
```

---

## Checkpoint: Before You Proceed

- [ ] `echo ${user_name}` prints your OSC username (not empty)
- [ ] `echo ${ANNOT}` prints `/fs/scratch/PAS3260/<your_username>/Annotation`
- [ ] `ls ${ANNOT}` shows the full directory tree including `augustus_config/`
- [ ] `ls ${SHARED_F2}/` lists all five shared resources (`.sif`, `databases/`, etc.)
- [ ] `ls ${ANNOT}/augustus_config/` shows the Augustus configuration files
- [ ] All 5 utility container `.sif` files exist in `${CONTAINERS}`
- [ ] Funannotate2 container prints a version string (Step 8)
- [ ] All utility containers print version strings
- [ ] `${ANNOT}/00_genome/Pf_assembly_raw.fasta` exists and has 5 sequences

---

*Previous: [00 Overview](00_overview.md) | Next: [Module 2 → RNA-seq Evidence](02_rnaseq_evidence.md)*
