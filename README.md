# GATK Joint-Genotyping Workflow

A portable [gwf](https://gwf.app/) workflow for short-read variant calling on a
SLURM cluster. Starting from paired-end FASTQ files and a reference genome, it
produces a filtered, multi-sample VCF using the GATK joint-genotyping
(GVCF → GenomicsDB → GenotypeGVCFs) best-practices flow.

To adapt it to a new dataset you only edit two text files — a YAML config and a
sample-sheet TSV — and (once) recreate a handful of conda environments. No code
changes are required.

---

## 1. Overview

The pipeline runs in five stages. Every step is one or more SLURM jobs managed by
gwf; dependencies between jobs are tracked with `.DONE` sentinel files written to
the `logs/` directory.

```
                      ┌─────────────────────────────┐
  reference.fasta ──► │ index_ref  (bwa-mem2 index)  │
                      │ gatk_dict  (.fai + .dict)    │
                      └──────────────┬──────────────┘
                                     │
   per sample (from sample sheet):   ▼
   FASTQ r1/r2 ─► align ─► addRG ─► fixmate ─► sort ─► markdup ─► MQ-filter ─► index
                                                                          │  (<sample>_final.bam)
                                                                          ▼
                         per contig:  HaplotypeCaller (-ERC GVCF)  ──►  merge (bcftools concat)
                                                                          │  (<sample>.g.vcf.gz)
                                     ┌────────────────────────────────────┘
                                     ▼   all samples
   per contig:  GenomicsDBImport ─► GenotypeGVCFs (-all-sites) ─► VariantFiltration
                                                                          │  (<contig>_GATK_filtered.vcf.gz)
                                                                          ▼
                         concat all contigs (bcftools)  ──►  joint/<joint_name>.vcf.gz
```

**Stage summary**

1. **Reference prep** — build the `bwa-mem2` index and create the `.fai` / `.dict`
   companions GATK requires.
2. **Per-sample alignment** — `bwa-mem2 mem` → add read groups → `fixmate` →
   coordinate `sort` → `markdup` (duplicates removed) → MAPQ/proper-pair filter →
   index. Intermediate BAMs are deleted as the chain progresses, leaving one
   `<sample>_final.bam`.
3. **Per-sample GVCFs** — `HaplotypeCaller` is run once per contig (massively
   parallel) and the per-contig GVCFs are concatenated into one
   `<sample>.g.vcf.gz`.
4. **Joint genotyping** — per contig, all sample GVCFs are imported into a
   GenomicsDB workspace, jointly genotyped (`-all-sites`, i.e. invariant sites are
   emitted too), and hard-filtered with the GATK germline expression.
5. **Final merge** — all per-contig VCFs are concatenated into the final indexed
   `joint/<joint_name>.vcf.gz`.

---

## 2. Requirements

- A **SLURM** cluster (jobs are submitted with `--account`).
- **conda / mamba** (e.g. miniforge). The workflow activates named conda envs
  inside each job.
- The **gwf** workflow manager (provided as one of the env files).
- **Scratch space** at `/scratch/$SLURM_JOBID/` on compute nodes — HaplotypeCaller
  and the GenomicsDB/GenotypeGVCFs steps write temporary data there. Adjust the
  `--tmp-dir` paths in `workflow_sources/workflow_templates.py` if your cluster
  uses a different scratch location.
- Enough memory for joint genotyping: the `joint_<contig>` jobs request **200 GB**
  (60 GB JVM heap). Tune in the templates if needed.

> The default resource requests (cores / memory / walltime) were sized for
> mammalian-genome, ~50-sample WGS. Smaller datasets will run fine but
> over-request; very large datasets may need more. See **Customizing**.

---

## 3. Repository layout

```
GATK_genotype_calling/
├── README.md
├── workflow.py                     # gwf entry point: loads every configurations/*.config.yaml
├── configurations/
│   └── config.template.yaml        # copy -> <name>.config.yaml and edit
├── sample_sheet.template.tsv       # copy -> sample_sheet.tsv and edit
├── envs/                           # conda environment.yml files (see Install)
│   ├── bwa2.yml
│   ├── samtools117.yml
│   ├── gatk4.yml
│   ├── bcftools.yml
│   └── gwf_new.yml
└── workflow_sources/
    ├── workflow_sources.py         # builds the target graph from the config + sample sheet
    ├── workflow_templates.py       # one function per job (the actual shell commands)
    └── sample_map_writer.py        # writes the GenomicsDBImport sample-name map

# created at run time (git-ignored):
├── steps/   index/  align/<sample>/  gvcf/<sample>/  joint/
└── logs/    *.DONE sentinels + sample_map.tsv
```

---

## 4. Install the conda environments

The job scripts activate environments **by name**, so the names must match
exactly: `bwa2`, `samtools117`, `gatk4`, `bcftools`. Create them once:

```bash
conda env create -f envs/bwa2.yml
conda env create -f envs/samtools117.yml
conda env create -f envs/gatk4.yml
conda env create -f envs/bcftools.yml
conda env create -f envs/gwf_new.yml   # the gwf manager (name can differ; see below)
```

(`mamba env create -f ...` works too and is faster.)

The pinned versions (bwa-mem2 2.2.1, samtools 1.19/1.17, GATK 4.5.0.0,
bcftools 1.16, gwf 2.1.1) are the ones the workflow was validated with. You may
relax the pins, but keep the **env names** unchanged for the four tool
environments. The gwf env name is only used in the run commands below, so you can
call it whatever you like.

---

## 5. Configure your run

### 5a. Reference genome

You need a reference FASTA and its `.fai` index **before** running anything (the
contig list is read from the `.fai`):

```bash
conda activate samtools117
samtools faidx /path/to/genome.fasta      # creates genome.fasta.fai
```

The `.dict` is created automatically by the `gatk_dict` target, so you do not need
to make it yourself.

### 5b. Sample sheet

Copy the template and fill in one row per sample (tab-separated):

```bash
cp sample_sheet.template.tsv sample_sheet.tsv
```

```
sample_id	r1	r2
sampleA	/data/sampleA_r1.fq.gz	/data/sampleA_r2.fq.gz
sampleB	/data/sampleB_r1.fastq.gz	/data/sampleB_r2.fastq.gz
```

- Columns are **tab-separated**: `sample_id`, `r1`, `r2` (absolute FASTQ paths).
- The header row is optional (skipped automatically); `#` comments and blank lines
  are ignored.
- `sample_id` becomes the read-group sample name (`RGSM`) and is used in every
  output path, so keep it filename-safe (`[A-Za-z0-9_-]`).
- Paired-end reads only; R1 and R2 must both be present.

### 5c. Config YAML

Copy the template to a name matching `*.config.yaml` (only those are loaded) and
edit every value with **absolute** paths:

```bash
cp configurations/config.template.yaml configurations/myproject.config.yaml
```

| Key | Meaning |
|-----|---------|
| `TAG` | Free-text label (bookkeeping only). |
| `account` | SLURM account jobs are charged to. |
| `project_name` | Human-readable project name. |
| `working_directory_path` | Absolute path to this checkout (where `workflow.py` is). |
| `output_directory_path` | Where results go; sub-dirs `index/ align/ gvcf/ joint/` are created here. |
| `log_directory_path` | Where `.DONE` sentinels and `sample_map.tsv` are written. |
| `scripts_path` | Absolute path to `workflow_sources/` (kept for project-layout parity). |
| `reference_genome` | Reference FASTA. |
| `reference_fai` | Its `.fai` index (must already exist — used to list contigs). |
| `sample_sheet` | Path to the TSV from step 5b. |
| `index_basename` | Prefix for the bwa-mem2 index files. |
| `joint_name` | Basename of the final VCF → `joint/<joint_name>.vcf.gz`. |
| `batch_size` | `GenomicsDBImport --batch-size` (lower it if that step runs out of memory). |

You can drop **multiple** `*.config.yaml` files into `configurations/` (e.g. one
per species/dataset) and they will all be registered in a single `gwf` graph.

---

## 6. Run

From the repo root:

```bash
conda activate gwf_new

gwf status            # show the full target graph and what would run
gwf run               # submit everything to SLURM (respecting dependencies)
gwf status -f         # follow progress; or: gwf logs <target>
```

`gwf run` submits only the jobs whose outputs are missing or out of date, so it is
safe to re-run after a partial failure — completed targets (those with a `.DONE`
sentinel) are skipped.

**Expected outputs** (under `output_directory_path`):

| Path | Contents |
|------|----------|
| `index/<index_basename>.*` | bwa-mem2 index. |
| `align/<sample>/<sample>_final.bam(.bai)` | Final filtered, dedup'd alignment + `.markdup.stat`. |
| `gvcf/<sample>/<sample>.g.vcf.gz(.tbi)` | Per-sample GVCF (all contigs). |
| `joint/<contig>_GATK_filtered.vcf.gz` | Per-contig joint-genotyped + filtered VCF. |
| `joint/<joint_name>.vcf.gz(.tbi)` | **Final multi-sample VCF.** |

---

## 7. Pipeline details

Each row is a gwf target template in `workflow_sources/workflow_templates.py`.
Targets fan out per sample (`<sample>_*`) and per contig.

| Target | Tool(s) | conda env | Cores / Mem / Time | Key parameters |
|--------|---------|-----------|--------------------|----------------|
| `index_ref` | `bwa-mem2 index` | bwa2 | 12 / 128g / 12h | index prefix = `index_basename` |
| `gatk_dict` | `samtools faidx`, `gatk CreateSequenceDictionary` | samtools117 → gatk4 | 2 / 8g / 2h | creates `.fai`/`.dict` if missing |
| `<s>_align` | `bwa-mem2 mem` \| `samtools view` | bwa2 | 20 / 64g / 24h | `-t 20`; SAM→BAM |
| `<s>_addRG` | `gatk AddOrReplaceReadGroups` | gatk4 | 4 / 16g / 12h | `RGSM=<s>`, `RGPL=illumina`, `RGLB=lib0`, `RGPU=unknown` |
| `<s>_fixmate` | `samtools fixmate -rm` | samtools117 | 16 / 32g / 12h | |
| `<s>_sort` | `samtools sort` | samtools117 | 16 / 32g / 12h | coordinate sort |
| `<s>_markdup` | `samtools markdup -r` | samtools117 | 16 / 32g / 12h | duplicates **removed**; stats → `.markdup.stat` |
| `<s>_mq` | `samtools view` | samtools117 | 16 / 32g / 12h | `-bq 60 -f 0x2 -F 0x4` (MAPQ ≥ 60, proper pairs, mapped) |
| `<s>_bai` | `samtools index` | samtools117 | 8 / 8g / 4h | |
| `<s>_<contig>_gvcf` | `gatk HaplotypeCaller` | gatk4 | 4 / 24g / 36h | `-ERC GVCF`, `-L <contig>`, `-Xmx20g`, scratch tmp |
| `<s>_gvcf_merge` | `bcftools concat`, `gatk IndexFeatureFile` | bcftools → gatk4 | 6 / 24g / 12h | merges all per-contig GVCFs |
| `joint_<contig>` | `GenomicsDBImport`, `GenotypeGVCFs`, `VariantFiltration` | gatk4 | 6 / 200g / 24h | `-all-sites`; hard filter (below) |
| `<joint_name>_concat` | `bcftools concat`, `gatk IndexFeatureFile` | bcftools → gatk4 | 6 / 24g / 12h | final VCF |

**Hard-filter expression** (FILTER tag `gatk_germline`, applied in `joint_<contig>`):

```
QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0
```

Variants are flagged, not dropped — failing records keep `gatk_germline` in the
FILTER column so you can subset later (e.g. `bcftools view -f PASS`).

---

## 8. Customizing

All of the following live in `workflow_sources/workflow_templates.py`:

- **Resources** — edit the `options = {'cores', 'memory', 'walltime', 'account'}`
  dict in the relevant template.
- **Filter thresholds** — edit the `--filter-expression` in `gatk_consolidate`.
- **Read-group fields** — edit `picard_addRG` (e.g. set a real `RGPL`/`RGLB` or
  per-lane `RGPU`).
- **MAPQ / flag filter** — edit the `samtools view` flags in `samtools_mq_filter`.
- **Invariant sites** — remove `-all-sites true` in `gatk_consolidate` to emit
  variant sites only.
- **Scratch path** — change the `/scratch/$SLURM_JOBID/` `--tmp-dir` paths if your
  cluster differs.
- **batch_size** — set in the config (`GenomicsDBImport --batch-size`).

The four tool conda-env names (`bwa2`, `samtools117`, `gatk4`, `bcftools`) are
hard-coded in the templates. If you must rename an env, update the matching
`conda activate <name>` lines as well.

---

## 9. Troubleshooting

- **`gwf status` errors about a missing `.fai`** — run `samtools faidx` on the
  reference first (step 5a); the contig list is read from it at graph-build time.
- **A target failed** — inspect logs with `gwf logs <target>` (add `-e` for
  stderr), fix the cause, then `gwf run` again. Only failed/incomplete targets are
  resubmitted.
- **Force a rerun** — delete the corresponding `logs/<target>.DONE` sentinel (and
  any stale output), then `gwf run`.
- **GenomicsDBImport OOM** — lower `batch_size`, or raise the `joint_<contig>`
  memory in the template.
- **Out of scratch space** — point `--tmp-dir` at a larger scratch filesystem.
- **Config not picked up** — it must match `configurations/*.config.yaml`
  (the shipped `config.template.yaml` is intentionally ignored).

---

## 10. Tools used

[bwa-mem2](https://github.com/bwa-mem2/bwa-mem2) ·
[samtools](https://www.htslib.org/) ·
[GATK4](https://gatk.broadinstitute.org/) ·
[bcftools](https://samtools.github.io/bcftools/) ·
[gwf](https://gwf.app/)

Please cite the respective tools when publishing results produced with this
workflow.
