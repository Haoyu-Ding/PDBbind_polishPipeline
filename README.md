# PDBbind Protein-Ligand Pipeline

This workspace contains a staged pipeline for building a clean protein-ligand 3D complex dataset from the PDBbind v2020R1 protein-ligand index and structure package.

## Goals

- Split ligands into `sugar` and `nonsugar`.
- Keep sugar ligands with minimal filtering because the expected count is small.
- Split nonsugar ligands into `rigid` and `flexible`.
- Reduce similarity for nonsugar entries in a reproducible way.
- Export a final `complex.pdb` for each retained entry.
- Keep only complexes with exactly one protein chain and one target ligand item.
- Remove water, metal ions, salts, cofactors, and any non-target hetero species from the exported complex.

## Dataset Inputs

- Structure archive: `/Users/teihiroshisakai/Documents/Codex/2026-06-17/pdbbind_v2020r1_official/P-L.tar`
- Index archive: `/Users/teihiroshisakai/Documents/Codex/2026-06-17/pdbbind_v2020r1_official/index.tar`
- Main index file: `index/INDEX_general_PL.2020R1.lst`

Although the item set is the PDBbind v2020R1 protein-ligand collection, the bundled structure files are reprocessed PDBbind v2024 structures according to the package README.

## Current Rule Set

The default dataset profile is `strict_v1`.

Hard exclusions:

- `covalent complex`
- `incomplete ligand` or `incomplete ligand structure`
- peptide-like ligands annotated as `(<n>-mer)`
- `isomer`
- entries whose ligand files cannot be parsed
- entries that fail the one-protein-chain or one-ligand-item structural constraints

Retained but tracked:

- `redundant to XXXX` annotations
- `NMR` structures
- lower-quality structures, unless later removed by representative selection

Nonsugar processing:

- Flexibility metric: RDKit strict rotatable bond count
- `rigid` if rotatable bonds `< 7`
- `flexible` if rotatable bonds `>= 7`
- Ligand deduplication within each flexibility bucket by Morgan fingerprint Tanimoto similarity
- Protein deduplication starts with exact sequence identity only in v1

Sugar processing:

- Sugar entries do not undergo flexibility grouping or similarity reduction in `strict_v1`
- Ambiguous sugar assignments are kept in a review bucket instead of being forced into a class

## Planned Execution Stages

1. Extract or stream-read the index archive and normalize the protein-ligand index.
2. Build a master manifest that links index metadata with structure file paths.
3. Classify ligands into `sugar`, `nonsugar`, or `ambiguous_sugar`.
4. Audit protein and ligand structures for chain count, ligand item count, hetero content, parseability, and sequence information.
5. Score nonsugar ligands for flexibility.
6. Deduplicate nonsugar entries by ligand similarity and then by exact protein sequence.
7. Export cleaned `complex.pdb` files.
8. Validate outputs and produce loss / summary reports.

## Directory Layout

```text
pdbbind_pl_pipeline/
  config/
  docs/
  scripts/
  src/pdbbind_pl/
  data/
    raw/
    extracted/
    interim/
    final/
    reports/
  reports/
  logs/
```

## Notes

- The pipeline is designed to be stageable and restartable.
- Intermediate outputs should be written as machine-readable tables such as CSV or Parquet.
- Filtering decisions should be recorded explicitly so that future relaxed profiles can be derived without recomputing all earlier stages.
- Project environment: `/Users/teihiroshisakai/Documents/Codex/2026-06-17/pdbbind-index-users-teihiroshisakai-documents-codex/pdbbind_pl_pipeline/.venv`
- Preferred runner: `/Users/teihiroshisakai/Documents/Codex/2026-06-17/pdbbind-index-users-teihiroshisakai-documents-codex/pdbbind_pl_pipeline/.venv/bin/python`
