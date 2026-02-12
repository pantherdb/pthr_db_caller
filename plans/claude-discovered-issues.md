# Plan: Fix 9 Confirmed Issues in pthr_db_caller

## Context

A thorough codebase walkthrough identified 9 confirmed issues ranging from critical init bugs to DRY violations. All issues are verified with exact line numbers. This plan addresses all 9 in priority order, grouped by file to minimize churn.

---

## Step 1: Fix `__int__` typo in RefProtEntry and RefProtFastaEntry

**File:** `pthr_db_caller/models/refprot_file.py`

**Problem:** `RefProtEntry` defines `__int__` (line 71) instead of `__init__`. `RefProtFastaEntry.__init__` (line 112) calls `RefProtEntry.__int__(self, ...)` — parent init is never properly called.

**Fix:**
- Line 71: Rename `def __int__` → `def __init__`
- Line 112: Change `RefProtEntry.__int__(self, ...)` → `super().__init__(uniprot_id, taxon_id)`

---

## Step 2: Fix double-self in GeneDatEntry.__init__

**File:** `pthr_db_caller/models/panther.py`

**Problem:** Line 179: `super().__init__(self, long_id, description, synonym, mod_id)` passes `self` as first arg to `object.__init__`, which is incorrect.

**Fix:**
- Change to `super().__init__()` (DatEntry has no `__init__`, so just call object's)

---

## Step 3: Deduplicate `extract_clade_name()`

**Files:** `pthr_db_caller/panther_tree_graph.py` (lines 29-53), `pthr_db_caller/taxon_validate.py` (lines 24-43)

**Problem:** Near-identical function in two files. The `panther_tree_graph.py` version has an extra `Opisthokonts → Opisthokonta` mapping.

**Fix:**
- Keep the more complete version in `panther_tree_graph.py` (it has the Opisthokonts fix)
- In `taxon_validate.py`: remove the local definition and import from `panther_tree_graph`:
  ```python
  from pthr_db_caller.panther_tree_graph import extract_clade_name
  ```

---

## Step 4: Convert global mutable state to instance state

### 4a: `THE_REST` in taxon_validate.py

**File:** `pthr_db_caller/taxon_validate.py`

**Problem:** Module-level `THE_REST = []` (line 21) accumulates species across calls.

**Fix:**
- Remove module-level `THE_REST`
- Make `get_all_species_from_tree()` return its results instead of mutating global state
- Pass the list as a parameter or return value through `get_all_species_from_tree()` → `append_species_to_table()`
- Update callers in `bin/align_taxon_term_table_species.py`

### 4b: `PC_NAMES` in protein_class.py

**File:** `pthr_db_caller/models/protein_class.py`

**Problem:** `PthrToPc.PC_NAMES = {}` (line 89) is a class-level dict that leaks across instances.

**Fix:**
- Move `PC_NAMES` to be an instance variable on the collection/caller rather than a class variable
- Or scope it to a module-level dict that gets explicitly cleared — the simplest change is to pass it as a parameter to `handle_pc_annots()` and store on the parsing context

---

## Step 5: Add index to GeneDatFile for O(1) lookups

**File:** `pthr_db_caller/models/panther.py`

**Problem:** `GeneDatFile.find_long_id()` (line 197) iterates all entries. Other similar classes use index dicts.

**Fix:**
- Add `self.by_long_id: dict = {}` to `GeneDatFile`
- Override `add_entry()` to populate the index (keyed by `str(entry.long_id)`)
- Rewrite `find_long_id()` to use the dict

---

## Step 6: Fix `OrganismDatFile.parse_organism_dat()` return type

**File:** `pthr_db_caller/models/panther.py`

**Problem:** `parse_organism_dat()` (line 251) returns a raw dict instead of following the `DatFile.parse()` pattern.

**Fix:**
- Keep the existing `parse_organism_dat()` static method as-is (it's used by callers that expect a dict)
- Add a proper `parse()` classmethod that returns an `OrganismDatFile` instance, with the dict accessible as an attribute
- This avoids breaking existing callers while bringing the API in line

---

## Step 7: Add validation to PthrSequence

**File:** `pthr_db_caller/models/panther.py`

**Problem:** `PthrSequence.__init__()` (line 8) does no format validation — crashes with unhelpful `ValueError` on malformed input.

**Fix:**
- Wrap the `split("|")` in a try/except that raises `ValueError` with a descriptive message including the malformed input
- Same for the `split("=")` on the uniprot part

---

## Step 8: Improve error handling in PaintIbaXmlParser.parse_xml()

**File:** `pthr_db_caller/models/paint.py`

**Problem:** `AssertionError` is caught and only printed (line 499), silently returning an empty collection.

**Fix:**
- Replace `print()` with `logging.warning()` so it integrates with the logging framework
- Add `import logging` at top of file if not present
- Keep the behavior of returning empty collection (callers depend on this) but make it visible via logging

---

## Step 9: Add `slim_maker.py` to setup.py scripts

**File:** `setup.py`

**Fix:**
- Add `"bin/slim_maker.py"` to the `scripts=[]` list

---

## Step 10: Make config paths resolve relative to package

**Files:** `pthr_db_caller/db_caller.py` (line 19), `pthr_db_caller/config.py` (line 4)

**Problem:** Default paths `config/config.yaml` and `config/build.yaml` assume CWD is project root.

**Fix:**
- Use `os.path.dirname(__file__)` to resolve paths relative to the package location:
  ```python
  _PACKAGE_DIR = os.path.dirname(os.path.abspath(__file__))
  _DEFAULT_CONFIG = os.path.join(_PACKAGE_DIR, "..", "config", "config.yaml")
  ```
- This way the defaults work regardless of CWD, while still allowing explicit override via the `config_path` parameter

---

## Verification

```bash
# Run full test suite after all changes
python -m unittest test -v

# Specifically verify the init fixes don't break anything
python -m unittest test.TestRefProtFastaFile -v
python -m unittest test.TestRefProtPantherIdMapping -v
python -m unittest test.TestPantherTreeGraph -v
python -m unittest test.TestXmlToGaf -v
python -m unittest test.TestOrthoXml -v

# Verify imports work after dedup
python -c "from pthr_db_caller.taxon_validate import extract_clade_name"
python -c "from pthr_db_caller.panther_tree_graph import extract_clade_name"

# Verify PthrSequence validation
python -c "from pthr_db_caller.models.panther import PthrSequence; PthrSequence('bad')" 2>&1 | grep -q ValueError

# Verify config path resolution
cd /tmp && python -c "from pthr_db_caller.db_caller import DBCallerConfig"
```

## Files Modified (summary)

| File | Steps |
|------|-------|
| `pthr_db_caller/models/refprot_file.py` | 1 |
| `pthr_db_caller/models/panther.py` | 2, 5, 6, 7 |
| `pthr_db_caller/panther_tree_graph.py` | 3 (keep function here) |
| `pthr_db_caller/taxon_validate.py` | 3, 4a |
| `pthr_db_caller/models/protein_class.py` | 4b |
| `pthr_db_caller/models/paint.py` | 8 |
| `setup.py` | 9 |
| `pthr_db_caller/db_caller.py` | 10 |
| `pthr_db_caller/config.py` | 10 |
| `bin/align_taxon_term_table_species.py` | 4a (caller update) |
