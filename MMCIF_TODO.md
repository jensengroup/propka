# mmCIF migration TODO

Goal: add native mmCIF input support without silently inheriting PDB fixed-width
limits. The risky part is not parsing `_atom_site`; it is preserving identity
through the existing PROPKA pipeline.

## Phase 0: protect current behavior

- [ ] Run the existing regression suite before structural changes.
- [ ] Keep `.pdb` input behavior byte-for-byte compatible unless a bug fix is
      explicitly covered by tests.
- [ ] Add focused tests before each identity refactor so PDB regressions are
      visible immediately.

## Phase 1: remove internal PDB identity limits

- [x] Replace single-character chain sorting.
      - Current hotspot: `ConformationContainer.sort_atoms_key()` uses
        `ord(atom.chain_id)`, which breaks for mmCIF multi-character chains.
      - Target: use tuple sorting, e.g. `(chain_id, res_num, icode, atom name)`.
- [x] Add structured atom/residue identity helpers on `Atom`.
      - Target residue key: `(chain_id, res_num, icode)`.
      - Target atom key: `(chain_id, res_num, icode, name)`.
      - Keep `residue_label` as display text only, not as the source of truth.
- [x] Include insertion code in all residue identity comparisons.
      - Hotspots:
        - `ConformationContainer.top_up_from_atoms()`
        - `ConformationContainer.find_group()`
        - `MolecularContainer.top_up_conformations()`
        - `energy.radial_volume_desolvation()`
        - `iterative.Iterative.__eq__()`
- [x] Update generated hydrogen identity.
      - Hotspots:
        - `hydrogens.make_new_h()`
        - `protonate.Protonate` residue-label updates
      - Generated atoms must inherit `icode` and structured keys from the parent.
- [ ] Audit labels and output formatting.
      - `Group.label` and `Atom.__str__()` use PDB-width formatting.
      - Output may remain PDB-like, but internal identity must not depend on
        truncated display labels.

## Phase 2: add native mmCIF atom-site parsing

- [x] Add `.cif` and `.mmcif` dispatch in `read_molecule_file()`.
- [x] Use `gemmi` for native mmCIF parsing.
      - Do not implement a fragile hand-written `_atom_site` tokenizer.
      - Current development environment has been checked with `gemmi 0.7.5`.
- [x] Map mmCIF fields explicitly:
      - record type: `group_PDB`
      - atom serial: `id`
      - element: `type_symbol`
      - atom name: prefer `auth_atom_id`, fallback `label_atom_id`
      - residue name: prefer `auth_comp_id`, fallback `label_comp_id`
      - chain: prefer `auth_asym_id`, fallback `label_asym_id`
      - residue number: `auth_seq_id`
      - insertion code: `pdbx_PDB_ins_code`
      - altloc: `label_alt_id`
      - model: `pdbx_PDB_model_num`
      - coordinates: `Cartn_x`, `Cartn_y`, `Cartn_z`
      - occupancy/B-factor: `occupancy`, `B_iso_or_equiv`
- [x] Treat mmCIF null markers `.` and `?` consistently.
      - Chain: choose a stable placeholder only if both chain ids are null.
      - Insertion code: normalize to `" "`.
      - Altloc: normalize to `"A"` for no-altloc only after preserving model
        separation.
- [x] Do not infer elements from atom names for mmCIF.
      - Use `_atom_site.type_symbol` directly.

## Phase 3: terminal and conformation semantics

- [ ] Rework N-terminus/C-terminus detection for mmCIF.
      - PDB path currently depends on `TER`, first `ATOM`, and `OXT`.
      - mmCIF may not provide equivalent TER records in `_atom_site`.
- [ ] Preserve multi-model handling.
      - Current conformation names are model + one-character altloc, e.g. `1A`.
      - mmCIF parsing must avoid collisions if altloc/model identifiers are
        unusual.
- [ ] Decide how to handle altloc groups.
      - Current code creates separate conformations from altloc.
      - Tests should cover `.`/`?`, `A/B`, and non-standard numeric altloc.

## Phase 4: number-limit tests

- [x] Multi-character chain test:
      - `auth_asym_id` values like `AA` and `AB` must remain separate.
- [x] Insertion-code test:
      - residues `(A, 10, A)` and `(A, 10, B)` must not merge.
- [x] Large residue number test:
      - `auth_seq_id > 9999` must parse and sort without truncation.
- [x] Large atom id test:
      - `_atom_site.id > 99999` must parse without hybrid36 conversion.
- [ ] Negative or non-contiguous residue number test if mmCIF input includes it.
- [ ] Chain filter and `--titrate_only` tests using multi-character chain IDs.

## Phase 4.5: cross-format verification

- [x] Convert a real BioLip `receptor_fixed.cif` to PDB with `gemmi` and compare
      PROPKA results.
      - Checked 4gc2 fixed receptor: CIF path and gemmi-converted PDB path both
        produced 111 summary groups with no pKa differences above 0.01.
      - Checked additional fixed receptors: 101m, 102m, 10gs, 5d46, 4yls, and
        4yml had identical summary groups and no pKa differences above 0.01.
      - Checked 1i94: summary group sets matched exactly; two GLU pKa values
        were swapped between GLU 85 A and GLU 245 A, consistent with
        order-sensitive coupled/equivalent group assignment rather than missing
        atoms or label loss.
- [x] Stress-test huge mmCIF input.
      - Checked 5fci fixed receptor: `_atom_site` has 365,380 atoms and max atom
        id 365,380. Native mmCIF parsing preserved the max atom id and loaded
        178,568 non-hydrogen atoms in about 5.6 s with default `keep_protons=False`.
      - Full all-chain PROPKA calculation did not finish within a 5-minute
        guard, which appears to be an algorithmic scale limit rather than an
        mmCIF parsing failure.
      - Checked 5fci chain A end-to-end: CIF and gemmi-converted PDB both
        produced 8,380 atoms, 2,168 groups, 284 summary groups, and no pKa
        differences above 0.01.
      - The gemmi-converted 5fci PDB uses hybrid36 atom serials past 99,999
        (e.g. A0000), which the existing PDB parser can decode.
- [x] Fetch an official RCSB PDB by PDB ID for smoke testing.
      - Checked 4GC2 from RCSB: PDB path runs successfully.
      - It is not expected to be group-count identical to the BioLip fixed
        receptor because the official PDB and fixed receptor differ in
        preprocessing/numbering/conformations.

## Phase 5: documentation and user-facing behavior

- [ ] Update README and command docs to mention `.cif` / `.mmcif`.
- [ ] Clarify whether output remains `.pka` only or whether converted PDB output
      helpers support long chain IDs.
- [ ] Add changelog entry describing native mmCIF support and remaining PDB-like
      output limitations.

## Review follow-up backlog

- [ ] Fix PDB N-terminus detection for insertion-code residues.
      - Current PDB parser still tracks the N-terminal residue with only
        `line[22:26]`; it does not include `line[26]`.
      - Risk: residues such as `A:1` and `A:1A` can both mark their `N` atom as
        `N+`, changing protonation, group type assignment, and pKa results.
      - Hotspot: `propka/input.py::get_atom_lines_from_pdb()`.
- [ ] Add an end-to-end mmCIF process-flow regression.
      - Current tests cover `read_mmcif()` field parsing, but not the full
        `read_molecule_file()` -> `setup_molecule()` -> top-up -> bonding and
        protonation -> group extraction -> pKa flow.
      - Target: one small PDB/mmCIF equivalent fixture with matching chains,
        groups, key labels, and selected pKa values.
- [ ] Add auth/label fallback edge-case tests for mmCIF.
      - Cover `auth_asym_id='?'`, `auth_seq_id='?'`, label-only chains,
        non-polymer rows with null `label_seq_id`, and chain filtering against
        the chosen chain namespace.
      - Current implementation is auth-first and label-fallback; that is a
        reasonable default, but the user-visible chain filter must be
        documented and tested.
- [ ] Decide how to represent mmCIF terminal semantics without PDB `TER`.
      - Current mmCIF parser marks only the first `ATOM` residue per
        `(model, altloc, chain)` as `N+` and uses `OXT`/`O''` for `C-`.
      - Risk: discontinuous polymer segments sharing one chain id may miss
        additional N-termini that PDB `TER` would have exposed.
- [ ] Explain or isolate the `3SGB.dat` golden-result drift.
      - Expected cause is insertion-code-aware residue identity, but this should
        be documented or covered by a focused test so unrelated PDB regressions
        are not hidden by broad golden updates.
- [ ] Audit PDB output helpers after adding internal mmCIF identity.
      - `Atom.make_pdb_line()` currently drops insertion code and still uses
        fixed PDB-style atom serial, chain, and residue-number formatting.
      - This may be acceptable if output remains legacy PDB-like, but it must be
        documented as an output limitation separate from native mmCIF input.

## Upstream existing issue: `ConformationContainer.get_chain()`

- [ ] Deep-dive and either fix or deprecate `get_chain()`.
      - Location: `propka/conformation_container.py::get_chain()`.
      - Current implementation returns
        `[atom for atom in self.atoms if atom.chain_id != chain]`.
      - The docstring says "Get atoms associated with a specific chain", so the
        predicate appears reversed; expected behavior is likely
        `atom.chain_id == chain`.
      - Confirmed with `git show HEAD~2:propka/conformation_container.py` that
        this predates the two mmCIF commits and is not newly introduced.
      - Confirmed with `rg` that no current code path calls `get_chain()`;
        output code filters chains directly with `g.atom.chain_id == chain`.
      - Risk level today: low runtime risk because it is unused, but high API
        surprise if downstream or future code calls it by name.
      - Before changing, add a focused unit test with two chains to lock the
        intended behavior.
