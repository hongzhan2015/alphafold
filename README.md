# FHV protein A modeling on CHTC

This workspace is a deployment and experiment starter for AlphaFold 3 (AF3)
on UW-Madison CHTC, tailored to Flock House virus (FHV) protein A and PDB
8FM9.

## Important modeling limitation

AF3 does **not** accept an EM density map or arbitrary spatial restraints as an
inference input. Its supported inputs include sequences, MSAs, single-chain
templates, ligands, covalent bonds, and random seeds. The mature-crown map must
therefore be used after AF3, during density fitting and real-space refinement.

The 2023 paper also distinguishes two assemblies:

- Proto-crown: one C12 ring (12 x 998 aa = 11,976 protein tokens).
- Mature crown: two stacked C12 rings, hence 24 protein A molecules in two
  conformations, plus an asymmetric central density.

A direct full C12 prediction is far above the 5,120-token size that the AF3
authors verified on one A100/H100 80 GB GPU. A mature C24 prediction is larger
again. Start with monomer, adjacent dimer, and four-chain local-wedge tests,
then assemble with C12 symmetry and fit/refine against the map.

## Recommended scientific workflow

1. Validate AF3 on the full-length monomer and compare it with 8FM9.
2. Run dimer and tetramer local-interface experiments with several seeds.
3. Treat residues 1-452 (the N-terminal floor/leg segment), 379-452 (the
   flexible Pol linker), 453-896 (Pol), and 897-998 (predicted disordered tail)
   explicitly when interpreting confidence.
4. Use 8FM9 as the experimentally anchored proto-crown model. Do not expect AF3
   to rediscover its C12 geometry from a 12-chain run.
5. For the mature crown, fit basal Pol/floor and apical Pol/leg domains into the
   cryo-ET map, impose C12 non-crystallographic symmetry where justified, and
   refine conservatively against the map. Keep the asymmetric central density
   separate.
6. When half-maps exist, cross-validate against the unused half-map. For the
   current sharpened map, report the missing half-map limitation and keep the
   interpretation at domain/secondary-structure level.

Suitable density-fitting tools include ChimeraX for initial placement and
Phenix real-space refinement or ISOLDE for restrained refinement. AF3 confidence
is not evidence of agreement with the EM map.

## CHTC storage and hardware prerequisites

The current official AF3 installation uses about 630 GB for expanded databases
(about 252 GB downloaded) and recommends up to 1 TB total disk, at least 64 GB
RAM for long data-pipeline targets, Linux, and NVIDIA compute capability 8.0 or
newer. CHTC `/home` is too small; arrange a `/staging` allocation of at least
700 GB, preferably about 1 TB, with CHTC facilitators before downloading data.

The supplied inference submit file requests an 80 GB, capability >= 8.0 GPU.
At CHTC this can match A100 80 GB, H100, or H200 nodes. It intentionally excludes
P100, RTX 2080 Ti, A100 40 GB, and L40/L40S nodes.

## Installation outline

On a machine that can build Docker images (or in the CHTC-supported container
build workflow):

```bash
git clone https://github.com/google-deepmind/alphafold3.git
cd alphafold3
docker build -t alphafold3 -f docker/Dockerfile .
```

Convert the image to `alphafold3.sif` following CHTC's Docker-to-Apptainer
guide, then place it at:

```text
/staging/u/NETID/alphafold3/alphafold3.sif
```

Download the official compressed model parameter file (subject to Google's
terms) into `models/`; AF3 reads `af3.bin.zst` directly. From the cloned AF3
source, download the databases into `databases/`:

```bash
mkdir -p /staging/u/NETID/alphafold3/models
wget -P /staging/u/NETID/alphafold3/models \
  https://storage.googleapis.com/alphafold3/af3.bin.zst
./fetch_databases.sh /staging/u/NETID/alphafold3/databases
```

The intended staging layout is:

```text
/staging/u/NETID/alphafold3/
  alphafold3.sif
  databases/
  models/
  runs/
  jax-cache/
```

Edit `AF3_ROOT` in both files under `chtc/` to replace `REPLACE_NETID`.

## Generate and submit experiments

For the newly supplied G4WWB0 FASTA and staged C24 strategy, follow
[`FHV_C24_WORKFLOW.md`](FHV_C24_WORKFLOW.md). Generate its recommended inputs:

```bash
python3 scripts/make_fhv_c24_inputs.py
```

The older 8FM9-construct JSONs for 1, 2, and 4 full-length copies can be
regenerated separately:

```bash
python3 scripts/make_af3_inputs.py
```

Copy this project to your CHTC access point, then run the CPU data pipeline for
one input:

```bash
condor_submit chtc/af3_data.sub input_json=inputs/fhv_a_monomer.json job_name=fhv_a_monomer
```

After the data job writes `fhv_a_monomer_data.json`, submit inference:

```bash
condor_submit chtc/af3_inference.sub job_name=fhv_a_monomer
```

Repeat for `fhv_a_dimer` and, only after those succeed, `fhv_a_tetramer`.
Use `condor_tail`, `.out`, `.err`, and `.log` to monitor and diagnose jobs.

## What would improve validation

The supplied sharpened map is sufficient for provisional rigid-domain fitting.
Unsharpened half-maps, if they can ever be recovered, would permit overfitting
checks. Also preserve the nominal/local resolution, applied symmetry,
sharpening/filtering history, and map origin/axis convention.
