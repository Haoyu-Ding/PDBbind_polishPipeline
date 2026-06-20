# Master Manifest Schema

The master manifest is the central table for all downstream filtering, scoring, deduplication, and export steps. Each row corresponds to one PDBbind protein-ligand entry from `INDEX_general_PL.2020R1.lst`.

## Identity and Paths

| Field | Type | Description |
| --- | --- | --- |
| `pdb_id` | string | Four-character PDB identifier in lowercase. |
| `release_year` | int | Release year from the index line. |
| `structure_archive_member_dir` | string | Directory prefix inside `P-L.tar` for the entry. |
| `protein_pdb_member` | string | Path to `<pdb_id>_protein.pdb` inside the structure archive. |
| `ligand_sdf_member` | string | Path to `<pdb_id>_ligand.sdf` inside the structure archive. |
| `ligand_mol2_member` | string | Path to `<pdb_id>_ligand.mol2` inside the structure archive. |
| `source_index_member` | string | Path to the source index file inside `index.tar`. |

## Raw Index Metadata

| Field | Type | Description |
| --- | --- | --- |
| `resolution_raw` | string | Raw resolution token from the index, such as `2.10` or `NMR`. |
| `resolution_value` | float nullable | Numeric resolution if available. |
| `structure_method_class` | string | Normalized method class, such as `xray` or `nmr`. |
| `binding_data_raw` | string | Raw affinity token. |
| `reference_raw` | string | Raw reference token after `//`. |
| `ligand_name_raw` | string | Original ligand annotation from the index comment. |
| `index_comment_raw` | string | Full comment section from the line after `//`. |

## Parsed Annotation Flags

| Field | Type | Description |
| --- | --- | --- |
| `flag_covalent_complex` | bool | True if the index comment contains `covalent complex`. |
| `flag_incomplete_ligand` | bool | True if the index comment indicates incomplete ligand information. |
| `flag_isomer_annotation` | bool | True if the index comment contains `isomer`. |
| `flag_redundant_annotation` | bool | True if the index comment contains `redundant to`. |
| `flag_peptide_like_mer_annotation` | bool | True if the ligand annotation matches `(<n>-mer)`. |
| `peptide_like_mer_size` | int nullable | Extracted `<n>` value for `(<n>-mer)` entries. |

## Ligand Parsing and Classification

| Field | Type | Description |
| --- | --- | --- |
| `ligand_parse_status` | string | `ok`, `failed_sdf`, `failed_mol2`, or `failed_all`. |
| `ligand_formula` | string nullable | Normalized ligand formula if parsed. |
| `ligand_heavy_atom_count` | int nullable | Heavy atom count. |
| `ligand_mol_wt` | float nullable | Molecular weight. |
| `ligand_rotatable_bonds` | int nullable | Strict rotatable bond count. |
| `ligand_component_count` | int nullable | Number of disconnected components in the ligand graph. |
| `ligand_component_rule` | string nullable | Classification note for single vs multiple components. |
| `ligand_class` | string | `sugar`, `nonsugar`, or `ambiguous_sugar`. |
| `ligand_class_reason` | string | Short machine-readable explanation of the classification. |
| `ligand_primary_code` | string nullable | Primary residue or ligand code inferred from the ligand file. |

## Protein Structure Audit

| Field | Type | Description |
| --- | --- | --- |
| `protein_parse_status` | string | `ok` or parse failure code. |
| `protein_chain_count` | int nullable | Number of protein chains based on ATOM records only. |
| `protein_chain_ids` | string nullable | Delimited ordered list of chain IDs. |
| `protein_selected_chain_id` | string nullable | Chain selected for export if valid. |
| `protein_sequence` | string nullable | Sequence of the selected protein chain. |
| `protein_sequence_sha1` | string nullable | Stable hash used for exact sequence deduplication. |
| `protein_atom_count` | int nullable | Number of ATOM records retained before export cleaning. |
| `protein_water_count` | int nullable | Count of `HOH` residues in the source protein file. |
| `protein_metal_count` | int nullable | Count of metal ions in the source protein file. |
| `protein_other_hetero_count` | int nullable | Count of other non-water HETATM groups in the source protein file. |

## Filtering and Dataset Membership

| Field | Type | Description |
| --- | --- | --- |
| `hard_filter_pass` | bool | True if the entry passes all hard filters for the active profile. |
| `hard_filter_fail_reasons` | string | Delimited list of fail reason codes. |
| `nonsugar_flexibility_bucket` | string nullable | `rigid`, `flexible`, or null for non-applicable entries. |
| `ligand_dedup_cluster_id` | string nullable | Cluster ID assigned during ligand deduplication. |
| `ligand_dedup_is_representative` | bool nullable | Whether the entry is the representative of its ligand cluster. |
| `protein_dedup_cluster_id` | string nullable | Cluster ID assigned during protein deduplication. |
| `protein_dedup_is_representative` | bool nullable | Whether the entry is the representative of its protein cluster. |
| `final_dataset_bucket` | string nullable | Final bucket such as `sugar`, `nonsugar_rigid`, or `nonsugar_flexible`. |
| `final_export_status` | string nullable | `pending`, `exported`, or failure code. |

## Export Artifacts

| Field | Type | Description |
| --- | --- | --- |
| `final_complex_pdb_path` | string nullable | Absolute path to the exported `complex.pdb`. |
| `final_protein_chain_pdb_path` | string nullable | Optional path to the cleaned protein chain file. |
| `final_ligand_structure_path` | string nullable | Optional path to the normalized ligand structure file. |
| `validation_status` | string nullable | Result of post-export validation. |
| `validation_notes` | string nullable | Free-text or code notes from validation. |

## Design Notes

- The manifest should preserve enough information to derive both `strict` and future `relaxed` dataset profiles.
- Every filtering or deduplication step should append explicit reason codes rather than mutating or overwriting upstream evidence.
- If Parquet is available it should be the primary format; CSV snapshots can be exported for inspection.
