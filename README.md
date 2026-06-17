# PD Metagenomics — Gut Microbiome Dysbiosis in Parkinson's Disease.

![Pipeline](https://img.shields.io/badge/pipeline-Kraken2%20%7C%20Bracken%20%7C%20Krona-blue)
![Dataset](https://img.shields.io/badge/dataset-PRJNA834801%20Wallen%202022-green)
![Status](https://img.shields.io/badge/status-complete-brightgreen)
![Platform](https://img.shields.io/badge/platform-Linux%20%2F%20WSL-orange)

## Overview

Shotgun metagenomic analysis of Parkinson's Disease gut microbiome using PRJNA834801 (Wallen et al. 2022).
Identifies depleted and enriched microbial taxa in PD vs healthy controls, computes alpha/beta diversity, and generates interactive Krona visualisations.

**Dataset:** PRJNA834801 — Wallen et al. 2022  
**Scope:** 2 PD + 3 Healthy Control samples  
**Sequencing:** Shotgun metagenomics, Illumina NovaSeq 6000, Paired-end

---

## Biological Question

> *Is the gut microbiome significantly different in Parkinson's Disease patients vs healthy controls? Which taxa are depleted or enriched?*

---

## Pipeline Summary
---

## Key Findings

| Metric | PD | HC |
|--------|----|----|
| Shannon Index | 2.84 | 3.26 |
| Bray-Curtis (within group) | 0.825 | 0.537 |

**Top depleted in PD:** *Blautia massiliensis*, *Sutterella wadsworthensis*, *Bifidobacterium pseudocatenulatum*  
**Top enriched in PD:** *Ligilactobacillus salivarius*, *Klebsiella quasipneumoniae*, *Citrobacter freundii*

---

## Tools

| Tool | Purpose |
|------|---------|
| Kraken2 2.1.7 | Taxonomic classification |
| Bracken | Abundance re-estimation |
| Krona | Interactive visualisation |
| FastQC 0.12.1 | Quality control |
| Trimmomatic | Adapter trimming |
| Python 3.11 | Diversity analysis |

---

## Repository Structure
---

## How to Reproduce

```bash
git clone https://github.com/Ganapathirajan/pd-metagenomics
cd pd-metagenomics
conda env create -f environment.yml
conda activate metagenomics_env
bash data/download_samples.sh
bash workflow/run_pipeline.sh SRR19064316
python workflow/diversity_analysis.py
python workflow/beta_diversity.py
```

---

## Samples

| SRR | Group | Size |
|-----|-------|------|
| SRR19064316 | PD | 8.34 GB |
| SRR19064317 | PD | 10.54 GB |
| SRR19064453 | HC | 6.8 GB |
| SRR19064454 | HC | 6.1 GB |
| SRR19064809 | HC | — |

---

## Limitations

- Pilot: 5 samples only (full cohort = 724 samples)
- p = 0.20 due to small n — directionally consistent with Wallen 2022
- Raw FASTQ + Kraken2 DB excluded (too large) — use download script

