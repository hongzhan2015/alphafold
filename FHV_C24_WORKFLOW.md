# FHV mature-crown C24 modeling workflow

## What AlphaFold 3 can and cannot do

The mature crown may be represented biologically as 24 protein A molecules,
but a direct full-length C24 AlphaFold 3 input contains 23,952 protein tokens.
That is far beyond the official 5,120-token A100/H100 80-GB benchmark and is
not a sensible first CHTC calculation. AlphaFold 3 also has no input field for
an MRC map or arbitrary density/distance restraints.

Use AlphaFold 3 to sample monomers and local interfaces. Use the experimental
8FM9 C12 proto-crown as the basal structural anchor, assemble a second C12 ring
from well-supported local/domain hypotheses, and then fit/refine the resulting
C24 model against the mature-crown map.

The supplied sharpened 8.5-A map is useful for domain placement. Because no
unsharpened half-maps are available, it cannot support independent half-map
cross-validation or a strongly restrained atomic refinement. Treat the result
as a provisional, map-consistent model rather than a new high-resolution
atomic structure.

## Sequence choice

`reference/FHV_proteinA_G4WWB0.fasta` is an exact workspace copy of the FASTA
you supplied. It is 998 aa long. It differs at 57 positions from
`reference/fhv_protein_a.fasta`, the sequence previously used for the 8FM9
construct analysis (94.3% identity). Confirm which isolate/construct produced the
mature map before final refinement. The AF3 input names below explicitly use
`g4wwb0` to prevent mixing them up.

## Generate the staged inputs

From the project root:

```bash
python3 scripts/make_fhv_c24_inputs.py
column -ts $'\t' inputs/fhv_g4wwb0/manifest.tsv
```

The default set is:

| Input | Tokens | Purpose |
|---|---:|---|
| full monomer | 998 | baseline/domain confidence |
| full dimer | 1,996 | adjacent-interface alternatives |
| full tetramer | 3,992 | local crown neighborhood |
| residues 55-396, C12 | 4,104 | complete floor-core hypothesis |
| residues 1-452, C4 | 1,808 | local floor/leg alternatives |
| residues 453-896, C4 | 1,776 | local polymerase-ring alternatives |
| residues 379-896, C4 | 2,072 | linker-to-Pol context |

The generator can write two slightly oversized C12 domain jobs for an H200 or
unified-memory experiment:

```bash
python3 scripts/make_fhv_c24_inputs.py --include-risky-c12
```

It can also emit a full C24 JSON only when explicitly requested:

```bash
python3 scripts/make_fhv_c24_inputs.py --include-impractical-c24
```

This switch documents the input; it is not a recommendation to submit it.

## Submit on CHTC

Start with the monomer. Run the CPU data pipeline first:

```bash
condor_submit chtc/af3_data.sub \
  input_json=inputs/fhv_g4wwb0/fhv_g4wwb0_full_monomer.json \
  job_name=fhv_g4wwb0_full_monomer
```

When that job succeeds, run inference:

```bash
condor_submit chtc/af3_inference.sub job_name=fhv_g4wwb0_full_monomer
```

Run the dimer next. Inspect chain-interface PAE, ipTM, clashes, and whether the
interface repeats consistently across seeds. Run the tetramer and domain C4/C12
jobs only after the smaller cases work. A high AF3 ranking score alone does not
show that an oligomer agrees with the mature EM density.

## Assemble against the mature map

1. Rigid-body fit the experimentally determined 8FM9 C12 ring into the basal
   crown density.
2. Compare the C4 predictions across seeds; retain interfaces that recur and
   fit the density without clashes.
3. Propagate one accepted asymmetric unit/interface around the map with C12
   symmetry rather than asking AF3 to infer all 24 chains simultaneously.
4. Fit the N-terminal floor/leg and polymerase portions separately where the
   mature map indicates different conformations. Keep residues 397-452 and
   897-998 flexible/low-confidence unless density supports them.
5. Refine conservatively with C12 NCS, secondary-structure, Ramachandran, and
   reference-model restraints. At 8.5 A, emphasize rigid domains and helices;
   do not optimize side-chain rotamers into sharpened noise.
6. Report map-model FSC and geometry, but state explicitly that half-map
   cross-validation was unavailable.

Do not apply C24 symmetry to the map: the architecture is two stacked C12
rings in distinct conformations, and the central density is asymmetric.
