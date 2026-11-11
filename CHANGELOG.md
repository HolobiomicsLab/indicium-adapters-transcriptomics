# Changelog

## Unreleased

- 2026-11-11: CI and development installs now use the published dependencies from PyPI.

## v0.14.0 — Full indicium-adapters #1 alignment: SSSOMMapping + workflow_coupling (2026-06-04)

### SKOS predicate routing + fallback-kind confidence

Pre-release refinement of the SSSOMMapping emission: every resolver now signals its sub-kind via an internal `_resolver_kind` sentinel inside the returned xrefs (popped by the dispatcher before xrefs is stored on the claim, so downstream consumers never see it). The sentinel routes the SSSOM fragment to the correct confidence floor AND the correct SKOS predicate:

* **Strict resolutions** (regex-anchored accession, gene-symbol lookup, anatomy lookup, canonical concept lookup) emit `skos:exactMatch` because the raw text is interchangeable with the CURIE.
* **LLM-tolerant fallbacks** emit `skos:closeMatch` because we either modified the input (trim-tail) or extracted the id from a longer wrapping string (embedded-pathway substring scan).

Without this distinction, every fragment got `skos:closeMatch` plus the `accession` / `concept` confidence floor regardless of which path actually fired — downstream SPARQL queries could not tell a heuristic substring hit from a regex-anchored accession match. The fix lets consumers filter on `predicate = skos:exactMatch` or `confidence >= 0.90` to drop the soft matches.

| Path | Predicate | Confidence |
|---|---|---|
| accession (regex hit) | `skos:exactMatch` | 0.99 |
| gene / anatomy / concept (lookup hit) | `skos:exactMatch` | 0.95 / 0.95 / 0.90 |
| embedded_pathway (substring fallback) | `skos:closeMatch` | 0.85 |
| trimmed_concept (trim-tail fallback) | `skos:closeMatch` | 0.80 |

5 additional pinning tests in `test_bulk.py` cover the routing behaviour: embedded-pathway emits closeMatch + 0.85; strict accession emits exactMatch + 0.99; trimmed concept emits closeMatch + 0.80; strict concept emits exactMatch + 0.90; the `_resolver_kind` sentinel never leaks into the claim's `xrefs`.

### Initial v0.14 work


Closes the two remaining "Leverage current indicium features" items of [indicium-adapters #1](https://github.com/HolobiomicsLab/indicium-adapters/issues/1) now that the upstream support shipped in `indicium` v1.11.0 (SSSOMMapping) and v1.10.0 (Evidence.workflow_coupling) on 2026-06-03. Additive: every claim shape from v0.13 still works; the new `sssom_mappings` field is emitted alongside the existing `ontology_terms` / `ontology_term_justifications` / `xrefs`.

### SSSOMMapping emission (indicium >= 1.11.0)

Every resolver hit in `enrich_claim()` now produces a first-class SSSOMMapping fragment in `claim["sssom_mappings"]` (a dict keyed by mapping id). Each fragment is shaped like `indicium.generated.model.SSSOMMapping`'s pydantic-constructor kwargs so Perspicacite can fold the entries straight into `IndiciumDocument.sssom_mappings` without round-tripping.

Fragment fields populated by the helper `_enrichment.build_sssom_mapping`:

| Field | Value | Notes |
|---|---|---|
| `id` | `"{subject_id}--{object_id}"` | Matches the convention from indicium's v1.11.0 test suite |
| `subject_id` | `"asb.lex:{slug(raw_text)}"` | Synthetic namespace for free-text LLM mentions; stable across runs for dedup |
| `subject_label` | the raw text | Carries the un-normalised mention |
| `predicate_id` | `"skos:closeMatch"` | Single predicate for all our resolver outputs |
| `object_id` | the resolved CURIE | Mirrors `ontology_terms[slot]` |
| `object_label` | optional | Pulled from the resolver's `xrefs.label` when available |
| `mapping_justification` | the SEMAPV string | Same value as `ontology_term_justifications[slot]` |
| `confidence` | per-resolver floor | accession=0.99, gene/anatomy=0.95, concept=0.90, embedded_pathway=0.85, trimmed_concept=0.80 |

Confidence floors are deliberately conservative — they let downstream consumers filter heuristic-fallback hits via `confidence >= 0.90` without losing the strict matches.

### `default_workflow_coupling` advertisement (indicium >= 1.10.0)

Both adapters now expose `default_workflow_coupling: str = "loose"` as a class attribute. Bulk RNA-seq pipelines (DESeq2 / edgeR / Salmon / STAR) and single-cell pipelines (CellRanger / STARsolo / Seurat / Scanpy) are file-based and individually reproducible per Suter et al. 2025 (DOI 10.1016/j.future.2025.107974), so the spec recommends `"loose"`. Perspicacite reads this attribute when constructing `Evidence` nodes for transcriptomics claims; runtime overrides remain possible.

### indicium dependency bump

`indicium >= 1.11.0` (was `indicium >= 1.4`). The bump unlocks the `SSSOMMapping` class and the `WorkflowCoupling` enum used by the new tests.

### Tests

12 new pinning tests across `test_bulk.py` and `test_singlecell.py`:

* `test_enrich_claim_emits_sssom_mapping_for_gene_subject` — fragment shape on the gene-subject path
* `test_enrich_claim_sssom_fragments_round_trip_through_indicium` — every fragment constructs through `SSSOMMapping(**fragment)`
* `test_enrich_claim_emits_sssom_for_accession_object` — accession-object emits with the 0.99 floor
* `test_enrich_claim_does_not_emit_sssom_when_no_resolution` — no emission on unresolved claims
* `test_default_workflow_coupling_is_loose` (bulk + sc) — pin the class attribute
* `test_workflow_coupling_value_is_in_indicium_enum` (bulk + sc) — round-trip through `Evidence(workflow_coupling=...)`
* `test_singlecell_enrich_emits_sssom_for_anatomy_object` — multi-fragment emission on a localized-in-tissue claim
* `test_singlecell_sssom_fragments_round_trip_through_indicium`

### Real-LLM smoke verification

`examples/perspicacite_real_llm_smoke.py` updated to also print the SSSOM fragments. Verified end-to-end with claude-sonnet-4-5: 5 fragments emitted across the 4 v0.13 smoke passages (signed-DE, marker, pathway), each with valid id / predicate / justification / confidence, all round-tripping through the pydantic constructor.


## v0.13.0 — indicium-adapters #1 spec alignment (2026-06-04)

Additive alignment with the [indicium-adapters #1](https://github.com/HolobiomicsLab/indicium-adapters/issues/1) contract. No breaking changes: every qualifier and CURIE prefix shipped before v0.13 still works exactly as it did. The deliberate divergence from the spec (two siblings `bulk-transcriptomics` + `single-cell-transcriptomics` instead of one `transcriptomics` adapter) is now documented in the README; rationale is the disjoint-qualifier structural invariant enforced by `tests/test_composition.py`.

### Qualifiers added to `BulkTranscriptomicsAdapter.qualifiers`

- `upregulated_in` / `downregulated_in` — signed specialisations of `differentially_expressed_in`. Extractors should prefer the signed form when log2FC sign is known and fall back to the generic form otherwise.
- `marker_for` — spec-canonical naming alongside the legacy `marker_of`. Both stay supported; new claims should prefer `marker_for`.
- `enriched_in_pathway` — narrows `enriched_for` for claims whose object resolves into REACTOME / KEGG / WP. The generic `enriched_for` is kept for GO / phenotype / disease enrichment.

The four new qualifiers fall through to the default resolver priority (`accession > gene > concept`) since they are feature-level claims; the study-level priority map is unchanged.

### CURIE prefixes added to `TRANSCRIPTOMICS_PREFIXES`

| Prefix | Expansion |
|---|---|
| `REACTOME` | `https://identifiers.org/reactome:` |
| `KEGG` | `https://identifiers.org/kegg.pathway:` |
| `WP` | `https://identifiers.org/wikipathways:` |

### Pathway-DB identifier recognition in `accession_to_curie`

The accession resolver now recognises pathway-database identifiers in the same way it recognises repository accessions, so `enrich_claim` can lift a bare `R-HSA-9012999` / `hsa04310` / `WP554` to its prefixed CURIE without the LLM extractor having to pre-prefix:

- Reactome: `R-{3-letter species}-{4-9 digits}` (e.g. `R-HSA-9012999`)
- KEGG pathway: `{2-4 letter org prefix}{5 digits}` (e.g. `hsa04310`, `map04310`, `ko04310`)
- WikiPathways: `WP{1-7 digits}` with optional `_r{digits}` revision suffix (e.g. `WP554`, `WP554_r118500`)

### `extraction_context()` updates

The prompt fragment now documents the direction-of-regulation rule (signed forms when log2FC sign is known, generic form otherwise), the pathway-vs-general enrichment rule (`enriched_in_pathway` for REACTOME/KEGG/WP, `enriched_for` for GO/phenotype/disease) and the `marker_for` vs `marker_of` rule (spec-canonical form preferred for new claims, legacy form kept for backward compatibility). It also notes that bare pathway IDs are auto-lifted to CURIEs.

### Tests

8 new pinning tests:

- `tests/test_bulk.py::test_qualifiers_include_spec_aligned_signed_and_pathway_terms`
- `tests/test_bulk.py::test_enrich_signed_de_routes_through_default_priority`
- `tests/test_bulk.py::test_enrich_marker_for_resolves_gene_subject`
- `tests/test_bulk.py::test_pathway_accession_to_curie` (12 parametrised cases)
- `tests/test_bulk.py::test_non_pathway_strings_do_not_misroute_to_pathway_prefix` (8 edge cases)
- `tests/test_bulk.py::test_enriched_in_pathway_resolves_reactome_object_end_to_end`
- `tests/test_bulk.py::test_enriched_in_pathway_resolves_kegg_object_end_to_end`
- `tests/test_bulk.py::test_enriched_in_pathway_resolves_wikipathways_object_end_to_end`
- `tests/test_namespaces.py::test_pathway_database_prefixes_present`


## v0.12.0 — Multi-species anatomy, broader vocabulary, layered concept index (2026-06-03)

Anatomy resolution becomes species-aware. The concept index grows from ~400 to ~547 verified CURIEs across 21 ontology sources. The single-cell adapter's anatomy lookup gains a `species` claim hint that routes terms to organism-specific ontologies when available.

### Multi-species `AnatomyIndex`

The anatomy lookup, previously UBERON-only, now merges UBERON (35 species-agnostic terms) with four organism-paired ontologies:

| Ontology | Species | Entries |
|---|---|---|
| UBERON | species-agnostic | 35 |
| FBbt | *Drosophila melanogaster* | 34 |
| WBbt | *Caenorhabditis elegans* | 21 |
| ZFA | zebrafish | 34 |
| EMAPA | mouse developmental | 28 |

Lookup is species-aware: a `species` hint on the claim routes to the organism-specific ontology when available, with UBERON winning by default. Cnidarian-specific terms (*physa*, *mesoglea*, *gastrodermis*) are documented as triaged out because UBERON has no entry for them.

### Vocabulary expansion (CL / GO / ZP / UBERON)

- **CL (Cell Ontology)** — 65 cell-type entries spanning cnidocyte, T cell / B cell subtypes, alveolar macrophage, Kupffer cell, pancreatic A/B/D, podocyte, keratinocyte, melanocyte, neuronal subtypes, stem cells.
- **GO biological process / molecular function** — ~56 concrete terms (Wnt / MAPK / BMP / Hippo signalling, apoptotic process, cell cycle, DNA repair, gene expression, oxidative phosphorylation, immune response, regeneration, ATP binding, kinase activity).
- **ZP zebrafish phenotype** — swim bladder / neural tube / notochord / fin / heart / circulation morphology phenotypes in the canonical `<entity> <quality>, abnormal` form.

### Concept-layer table

Concept lookup is now layered with a documented priority ladder (paper-modelling vocab > OLS4 corpus > bio.tools > single-ontology indexes). The README now lists every layer + its source ontologies in a table.

### `xrefs[slot]` shape for anatomy hits

Anatomy resolvers now return `{"curie", "label", "ontology", "description"}` matching the concept-resolver shape, plus an additional `"species"` key when the hit came from an organism-specific ontology via a species hint.


## v0.11.0 — Real-corpus smoke test on 16 papers, 45 claims (2026-05-28)

Expands ``examples/test_on_papers.py`` from the 2-paper / 10-claim
minimum-viable smoke test (v0.10.0) to a **diverse 16-paper /
45-claim corpus** covering the six supported organisms across bulk,
single-cell, spatial and long-read modalities — plus four
methodology reviews for vocab-density coverage. 4 new PDFs added to
``articles/`` (Brunet Avalos 2019 Drosophila, Packer 2019 C. elegans,
Wagner 2018 zebrafish, Tirosh 2016 melanoma) extend the prior 6-PDF
mini-corpus from v0.3.1.

### Corpus composition

**12 primary research papers** (10 with PDFs in ``articles/``, 2
based on prior corpus extraction — Warner / Cheung):

| Paper | Organism | Modality | Source |
|---|---|---|---|
| Warner et al. 2025 | *Nematostella* | bulk regeneration | prior |
| Cheung et al. 2024 | *Nematostella* | spatial Tomo-seq | prior |
| Sebé-Pedrós et al. 2018 (PMID 29856953) | *Nematostella* | scRNA atlas | PDF |
| Layden et al. 2016 | *Nematostella* | bulk microarray | PDF |
| Tabula Muris 2018 (PMID 30283141) | mouse | scRNA atlas | PDF |
| Tabula Sapiens 2022 (PMID 35549404) | human | scRNA atlas | PDF |
| Brunet Avalos et al. 2019 (PMID 31746739) | *Drosophila* | scRNA brain | PDF (new) |
| Packer et al. 2019 (PMID 31488706) | *C. elegans* | scRNA embryogenesis | PDF (new) |
| Wagner et al. 2018 (PMID 29700227) | zebrafish | scRNA embryo | PDF (new) |
| Tirosh et al. 2016 | melanoma | scRNA cancer | PDF (new) |
| Stickels et al. 2021 (Slide-seqV2) | mouse | spatial method | PDF |
| Pardo-Palacios et al. 2024 (LRGASP) | multi | long-read benchmark | PDF |

**4 methodology reviews** (PubMed efetch abstracts):

| Review | PMID |
|---|---|
| Conesa et al. 2016 — *Survey of best practices for RNA-seq* | 26813401 |
| Hwang et al. 2018 — *scRNA-seq technologies and pipelines* | 30089861 |
| Marx 2021 — *Method of the Year: spatially resolved transcriptomics* | 33408395 |
| Stark et al. 2019 — *RNA sequencing: the teenage years* | 31341269 |

### Methodology

Claims hand-curated with traceable sources:

- **PDF papers**: accessions, instruments and methods extracted
  via regex / abstract reading from the actual PDF text. The most
  rigorous tier. Notable verified data this release: Brunet Avalos
  GSE134722 + Illumina NovaSeq 6000, Packer GSE126954 + NextSeq,
  Wagner GSE112294, Tirosh GSE72056 + GSE77940 + Smart-seq2,
  Sebé-Pedrós GSE95723 + MARS-seq, Layden E-MTAB-4831,
  Tabula Sapiens GSE201333 (which had been a memory-based guess —
  now verified).
- **Memory-based** (Warner / Cheung): established from prior
  v0.3.1 extraction work; claims are verifiable against the
  published papers.
- **Review abstracts**: pulled via PubMed efetch; claims drawn
  from explicit mentions of tools / protocols in the abstract
  text.

Method explicitly: **not** full-text LLM extraction. The point of
v0.11.0 is reproducibility — every claim has a traceable source.

### Results

- **44 / 45 claims enriched (~98%)**.
- **18 distinct ontologies fired**: EFO, GEO, HGNC, biotools, NCIT,
  WBbt, CHEBI, CL, UBERON, MGI, ZFA, FlyBase, WormBase, ZFIN, DOID,
  ArrayExpress, FBbt, EMAPA. Two new this release vs v0.10.0
  (the minimum-viable smoke test):
  - **ArrayExpress** — via Layden 2016 E-MTAB-4831 accession.
  - **DOID** — via Tirosh 2016 *melanoma* clinical claim.
- **SHACL validation conforms** on the synthesised RDF data graph
  (116 triples).

### Single unresolved term — honest signal

- ``scVI`` — known absent from bio.tools, documented at v0.6.x.
  Tabula Sapiens uses scVI; LLM extraction will surface this on
  any single-cell paper using probabilistic modelling. Add a
  ``scvi-tools`` direct-ID attempt in a v0.11.x bio.tools harvest
  refresh if benchmark frequency justifies.

### Added

- ``examples/test_on_papers.py`` rewritten with the 16-paper corpus
  + per-paper + per-ontology coverage matrix + traceable per-claim
  sourcing notes. Distinct from ``examples/enrich_warner2018.py``
  (the 25-claim hand-built resolution-coverage showcase).

### Notes

- No data or production code changed in v0.11.0 — concept index
  stays at 19 sources / ~447 CURIEs from v0.9.0. Only the example
  script + 4 new PDFs in ``articles/`` (gitignored).
- The script's coverage matrix is the v0.11.0 deliverable: a
  one-page empirical report on what the adapter catches on a
  realistic 16-paper transcriptomics corpus, reproducible by
  anyone with the package installed.

## v0.10.0 — Integration smoke tests + Protocol-conformance suite (2026-05-28)

Quality + v1.0-prep release inspired by an audit of the sibling
``indicium-adapters-epigenomics`` repo. Brings in the patterns the
sibling adapter had at v0.1 that we hadn't yet adopted, without
deferring v1.0 itself.

### Added — Tests

- ``tests/test_protocol_conformance.py`` (15 tests) — centralises
  the cross-cutting Protocol invariants previously scattered in
  ``test_bulk.py`` / ``test_singlecell.py``, **plus** depth checks
  the sibling adapter has at v0.1:
  - Every adapter exported by the package is discoverable via the
    ``indicium.adapters`` entry-point group (catches missing entry-
    point registration).
  - ``domain_id`` is kebab-case (lowercase + hyphens).
  - All ``qualifiers`` are snake_case — extraction prompts and SPARQL
    templates depend on the convention.
  - Every ``ontology_prefixes`` map entry dereferences to a valid
    http(s) URI ending with the expected separator (``:`` / ``/`` /
    ``_`` / ``#``).
  - ``shacl_shapes()`` returns an ``rdflib.Graph`` with at least one
    ``sh:NodeShape``.
- Parametrised over ``ALL_ADAPTER_CLASSES`` — new sibling adapters
  added to that tuple automatically inherit the full suite.

### Added — Examples

- ``examples/quickstart.py`` — 30-line "Hello World" runnable demo:
  discover adapters → enrich one synthetic claim per adapter →
  validate the composed SHACL shapes graph parses cleanly under
  pyshacl. Distinct from the Warner-flavoured fuller demo. Use this
  as the bring-up smoke test in a new environment.
- ``examples/test_on_papers.py`` — pipeline integration smoke test
  on hand-extracted claims from real papers (Warner et al. 2025
  and Tabula Muris 2018 in v0.10.0). 10 claims across both
  adapters; each goes through ``enrich_claim()`` → synthesised
  RDF data graph → pyshacl validation against the composed shapes
  graph. Proves the full claim-dict → CURIE → RDF → SHACL round-
  trip is plumbed. **v0.11.0 will expand this to a SciLEx-driven
  diverse 10-paper corpus across all 6 supported species.**

### Notes

- Sibling-audit verdict (vs ``indicium-adapters-epigenomics``):
  we ship more (CHANGELOG, ROADMAP, OLS4 audit, 4 parsers, 6
  species, ``scripts/check.ps1``) than they do at v0.1.0, but had
  the conformance / quickstart / real-papers gaps that this
  release closes. v1.0.0 will add: CITATION.cff, README badges,
  ontology-provenance table, final tag.
- No data added in v0.10.0 — concept index stays at 19 sources /
  ~447 CURIEs from v0.9.0.

### Tests

- 126 / 126 pass, 92.35% coverage. The 15 new parametrised
  conformance tests run in <1s.

## v0.9.0 — Broader ontology coverage for benchmark consistency (2026-05-28)

Following Dina's call to "do it properly" before the upcoming ASB
benchmark phase, v0.9.0 extends the v0.8.0 audit-driven set with
**Tier-1 universal** (vocabularies that fire on almost any transcripto
paper) and **Tier-2 organism-paired** (per-species anatomy / phenotype
that complete the v0.6.x gene-index coverage). Tier-3 niche stays
deferred.

### Added — Tier-1 universal indexes

- **MS** (Mass Spectrometry Ontology, PSI-MS) — **21 entries**
  covering ionisation methods (electrospray, MALDI, APCI),
  mass analysers (quadrupole, time-of-flight, orbitrap, FT-ICR),
  combined methods (GC-MS), acquisition modes (DDA / DIA / SRM),
  sample-prep terms (trypsin, in-gel digestion). Important for the
  proteomics overlap with the metabolomics adapter.
- **UO** (Units of measurement Ontology) — **35 entries**
  covering mass (g/kg/mg/µg/ng/pg), volume (l/ml/µl/nl),
  concentration (molar / millimolar / micromolar / nanomolar / picomolar),
  time, temperature, counts, length. Essential for normalising
  quantitative claims.
- **MOD** (Protein Modifications, PSI-MOD) — **19 entries** after
  candidate-list refinement. PSI-MOD uses ``<verb-form> residue`` /
  ``<verb-form> L-<aa>`` patterns rather than the ``-ation`` PTM-
  class names common in literature; candidates aligned with that
  convention catch ``phosphorylated residue``, ``acetylated residue``,
  ``methylated residue``, ``ubiquitinated residue``, ``glycosylated
  residue``, ``palmitoylated residue``, ``myristoylated residue``,
  ``farnesylated residue``, etc.

### Added — Tier-2 organism-paired anatomy / phenotype

Each ontology completes the per-species v0.6.x gene index. Single-
species papers get end-to-end resolution (organism gene + organism
anatomy + organism phenotype) without falling back to UBERON for
non-mammalian anatomy.

- **FBbt** (Drosophila anatomy) — **34 entries**. Pairs with v0.6.2
  FlyBase gene index. Covers imaginal discs (wing / eye / leg /
  antennal), adult organ systems (ovary / testis / brain / ventral
  nerve cord / gut), life stages, muscle, fat body.
- **WBbt** (C. elegans anatomy) — **21 entries**. Pairs with v0.6.2
  WormBase gene index. Covers reproductive system (germ line, gonad,
  spermatheca), digestive (pharynx, intestine), nervous (neuron
  types, nerve ring), muscle (body wall, pharyngeal).
- **ZFA** (Zebrafish anatomy) — **34 entries**. Pairs with v0.6.2
  ZFIN gene index. Covers early-development structures (somite,
  notochord, neural tube), adult organs (heart, brain, retina,
  liver, swim bladder), fins (pectoral, pelvic, dorsal, caudal),
  pharyngeal arches. Note: life-stage candidates resolve via the
  imported **ZFS** zebrafish-stage namespace.
- **EMAPA** (Mouse Developmental Anatomy) — **28 entries**. Pairs
  with v0.6.0 MGI gene index. Covers embryo / fetal / neonate /
  adult life stages, somite / neural-tube / neural-crest, limb-bud
  development, brain regions, organ primordia, placenta / yolk-sac.
- **XPO** (Xenopus Phenotype) — **9 entries**. Regeneration /
  morphology phenotypes for the canonical regeneration model.

**Total new in v0.9.0**: **201 CURIEs across 8 ontologies**.
Concept index grows from ~246 rows (v0.8.0) to **~447 rows** loaded
from 19 source TSVs.

### Deferred to v0.9.x (poor OLS4 fit)

Four ontologies in the v0.8.0 audit shortlist had yields too low to
ship at v0.9.0 quality bar. Each is annotated in
``build_extra_ontology_indexes.py`` with the per-ontology rationale
and will need a dedicated harvester or a different source:

- **EC** (Enzyme Commission, ~13k terms in OLS4) — re-exports GO
  molecular-function terms without human-readable labels. 3 / 30
  yield. Defer until we wire ExPASy enzyme.dat directly.
- **RXNO** (Reaction Ontology) — mostly CHEBI / BFO imports without
  named reaction labels. 2 / 16. Defer.
- **HCAO** (Human Cell Atlas Ontology) — re-publishes CARO / CL /
  CLO / GO terms; redundant with our existing v0.5.0 CL index.
  0 / 17. Defer.
- **ZP** (Zebrafish Phenotype) — labels follow a peculiar
  ``<entity> abnormal, abnormal`` pattern that doesn't match the
  natural-language ``abnormal <entity>`` form. 1 / 16. Refactor
  the harvester to flip word order before re-attempting.

### Changed — runtime + harvester

- ``_load_one_concept_tsv`` is now **tolerant of missing files**:
  returns ``[]`` rather than raising ``FileNotFoundError`` when a
  bundled TSV is absent. Still strict on malformed headers
  (raises). This lets the package import cleanly during maintainer
  build steps that re-generate a subset of TSVs sequentially. Fix
  for a real bug that surfaced during the v0.9.0 harvest sequence.
- ``build_extra_ontology_indexes.REGISTRY`` extended with 8 new
  entries (3 Tier-1 + 5 Tier-2). The 4 deferred ontologies are
  commented out in-place with their rationale, so a future
  harvester author can see what was tried.

### Tests

- 8 new pinned tests, one per new ontology. Includes a regression
  guard for the somite collision between EMAPA and ZFA (the loader's
  priority ladder puts mouse EMAPA first, so a bare ``somite``
  resolves via EMAPA; ZFA-specific tests use ZFA-unique terms).
- Size invariant updated to ≥447 rows across 19 source TSVs.

### Demo

- 4 new claims in ``examples/enrich_warner2018.py`` showing the
  cross-axis multi-ontology coverage (EMAPA limb bud, FBbt wing
  disc, MS orbitrap, PSI-MOD phosphorylated residue). 25/25 claims
  enrich.

## v0.8.0 — Top-7 OLS4 ontology indexes for benchmark exhaustiveness (2026-05-28)

### Added — OLS4 audit + 7 high-priority ontology indexes

Following an audit of the OLS4 catalogue (280 ontologies, see
``articles/_ols4_audit_report.md``), seven high-priority ontologies
were harvested into dedicated bundled indexes. The selection comes
from a curated short-list of ontologies most likely to fire on a
diverse transcriptomics corpus in the upcoming ASB benchmark phase.

Harvest mechanism: same per-ontology cascade as v0.5.0 CL — curated
candidate names → OLS4 ``exact label/synonym`` then ``prefix`` →
re-verify by IRI lookup. Hand-curation is limited to the candidate
name list per ontology; the matching CURIE comes from OLS4.

- **CHEBI** — chemicals / drugs / ligands. **36 entries** including
  ``U0126`` → CHEBI:90693 (the MEK inhibitor used in Warner 2025,
  now resolvable end-to-end), ``cisplatin``, ``paclitaxel``,
  ``estradiol``, ``dopamine``, ``ATP``, ``glucose``,
  ``lipopolysaccharide``, etc.
- **HP** — Human Phenotype Ontology. **30 entries** for clinical
  transcripto vocabulary: ``intellectual disability``,
  ``microcephaly``, ``seizure``, ``cardiomyopathy``,
  ``hepatomegaly``, ``anemia``, etc.
- **MP** — Mammalian Phenotype Ontology. **19 entries** for mouse
  studies: ``abnormal body weight``, ``abnormal heart morphology``,
  ``decreased fertility``, ``abnormal T cell physiology``, etc.
- **DOID** — Human Disease Ontology. **29 entries** for clinical /
  cancer transcriptomics: ``cancer``, ``breast cancer``,
  ``lung cancer``, ``Alzheimer's disease``, ``COVID-19``,
  ``inflammatory bowel disease``, etc.
- **PR** — Protein Ontology. **18 entries** for protein-level
  claims: ``alpha-synuclein``, ``tubulin``, ``myosin``,
  ``ubiquitin``, etc. Note: PR is heavily nested and generic
  candidate names (``transcription factor``, ``actin``) sometimes
  land on specific sub-class entries via synonym matching —
  the CURIEs are valid but not always the most canonical.
- **OBA** — Ontology of Biological Attributes. **7 entries** for
  cross-species traits: ``body size``, ``heart rate``,
  ``blood pressure``, ``immune response``, ``circadian rhythm``.
  Note: OBA imports terms from PATO and GO; some rows surface
  PATO / GO CURIEs while the ``ontology`` column records OBA as
  the harvesting source.
- **ENVO** — Environment Ontology. **15 entries** for holobiont /
  sample-context vocabulary: ``marine biome``, ``coral reef``,
  ``intertidal zone``, ``sediment``, ``salinity``, etc.

**Total**: 154 new CURIEs across the 7 ontologies. Concept index
grows from ~92 rows (v0.5.0) to **~246 rows** (v0.8.0).

### Added — runtime + audit artifact

- Loader extended in ``_enrichment._load_concept_index_rows`` to
  read all 7 new TSVs in priority order after the existing four
  (static / OLS4 corpus / bio.tools / CL). The whole concept index
  is now 11 sources merged.
- ``data/build_extra_ontology_indexes.py`` — single REGISTRY-driven
  harvester that handles all 7 ontologies. CLI flag ``--only <id>``
  lets you re-harvest one ontology in isolation when its candidate
  list grows.
- ``articles/_ols4_ontology_audit.py`` + ``_ols4_audit_report.md``
  — the OLS4 catalogue audit committed as a reference document.
  ~28 curated picks plus a heuristic bucket-dump of the full
  280-ontology catalogue. Future ontology additions should consult
  the audit's queue.
- ``articles/`` directory excluded from ruff (dev workbench, not
  package).

### Tests

- 7 new tests in ``test_bulk.py`` pinning representative resolutions
  per ontology — CHEBI ``U0126``, HP ``intellectual disability``,
  MP ``abnormal body weight``, DOID ``breast cancer``, PR
  ``alpha-synuclein``, OBA ``heart rate``, ENVO ``marine biome``.
- ``test_concept_index_carries_all_11_concept_sources`` updates the
  size invariant to ≥246 rows.
- ``test_enrich_claim_is_safe_noop_for_unknown_concept_terms`` had
  to be updated: pre-v0.8.0 it used ``U0126`` as a "supposedly
  unresolvable" example, but v0.8.0 now resolves it via ChEBI —
  the test pivots to two genuinely bogus terms.

### Demo

- Three new claims in ``examples/enrich_warner2018.py`` showing
  end-to-end resolution through ChEBI / HP / DOID. 21/21 claims
  enrich.

### Notes on scope

The OLS4 audit identified 28 high-priority candidate ontologies.
v0.8.0 ships the top-7 absolute priority subset (CHEBI / HP / MP /
DOID / PR / OBA / ENVO). The remaining 21 candidates are queued
in ROADMAP — they get added incrementally as benchmark misses
surface them, **not** preemptively. Adding an ontology after the
fact is cheap (1-2h per ontology); adding dead weight before
benchmark signal is wasteful. See ``project_asb_benchmark_exhaustiveness``
in the memory bank for the principle.

### Known anomalies

- PR generic candidates (``transcription factor`` → PR:000000032
  ``transcription factor Sp``; ``actin`` → ``actin (yeast)``) —
  documented in the v0.8.0 test docstring.
- CHEBI ``hydrocortisone`` resolves to CHEBI:135746
  ``hydrocortisone aceponate`` (the ester) because the OLS4 search
  ranks the synonym match before the bare ``hydrocortisone`` term.
- OBA ``fertility`` resolves to PATO:0000274 because OBA imports
  PATO terms.

These are all documented in the per-ontology unresolved logs
(``<onto>_unresolved.txt``) and the v0.8.x roadmap.

## v0.7.0 — Repository metadata parsers (2026-05-28)

### Added — `parsers/` subpackage

Four parsers, one per common repository metadata format, sharing a
uniform contract: take the raw input (text / JSON / directory) and
return a ``list[dict]`` of claim dicts directly consumable by
``BulkTranscriptomicsAdapter.enrich_claim``.

- **`parse_geo_soft(text)`** — GEO SOFT format (the GEO's native text
  serialisation). Tokenises ``^TYPE = id`` records and
  ``!field = value`` lines, then walks ``^SERIES`` + ``^SAMPLE``
  records. Emits ``consolidates_dataset`` claims for each
  ``!Series_sample_id``, ``sequenced_on`` claims (deduplicated)
  per distinct ``!Sample_instrument_model``, and
  ``uses_database PubMed:<id>`` from ``!Series_pubmed_id``.
- **`parse_ena_xml(xml_text)`** — ENA / SRA study & experiment XML.
  Walks all ``<EXPERIMENT>`` elements, links each to its
  ``<STUDY_REF>``, and reads the vendor-named child of ``<PLATFORM>``
  to extract ``<INSTRUMENT_MODEL>``. Emits ``consolidates_dataset``
  + ``sequenced_on`` claims.
- **`parse_biostudies_json(text_or_obj)`** — BioStudies JSON
  documents. Recursively flattens the ``section.attributes`` tree
  into ``{name → value}``, recognises ``Organism`` / ``Pubmed`` /
  ``Instrument`` / ``Platform``, emits study-level claims under the
  document's ``accno`` accession.
- **`parse_isa_tab(directory)`** — ISA-Tab investigation/study/assay
  bundle. Reads ``i_*.txt`` for ``Study Identifier`` and
  ``Study PubMed ID``, walks all ``s_*.txt`` for sample
  ``Characteristics[organism]``, walks all ``a_*.txt`` for assay
  ``Parameter Value[instrument model]``. Emits study-level claims
  under the study identifier.

### Added — species hint propagation

All four parsers detect a recognised model-organism name in the
relevant field (``!Sample_organism_ch1`` for GEO,
``Characteristics[organism]`` for ISA-Tab, ``Organism`` attribute for
BioStudies — ENA's organism is in the SAMPLE descriptor which v0.7.0
doesn't yet parse) and attach a ``species`` hint to every emitted
claim. The hint then drives the v0.6.0 species-aware gene resolver
downstream, so gene mentions in a Tabula Muris-style claim emit
correctly to MGI even when the LLM normalised casing.

Recognised organisms: human / mouse / Drosophila / C. elegans /
zebrafish / Nematostella.

### Tests

- 12 new tests in ``tests/test_parsers.py`` covering: per-parser
  field extraction (organism, pubmed, instrument), the species-hint
  detection path, end-to-end enrichment through
  ``BulkTranscriptomicsAdapter`` for GEO SOFT (asserts
  ``GSE168938 → GEO:GSE168938`` + ``Illumina NovaSeq 6000 →
  EFO:0008637`` chain works), JSON-string vs dict input for
  BioStudies, missing investigation fallback for ISA-Tab.

### Notes

- v0.7.0 is the **pipeline-entry** for repository metadata —
  upstream of the adapter, downstream of the registry's native
  format. Field coverage is minimum-viable: every parser extracts
  the 3-4 fields most commonly populated; deepening to e.g.
  ``Sample_characteristics_ch1`` sub-fields or ISA-Tab's
  ``Comment[...]`` columns is left for v0.7.x as real use-cases
  surface them.
- The parsers are **stateless and pure** — no network calls, no
  files written. Each function takes its input and returns its
  list; the caller chooses how to feed the claims onward (most
  often into ``adapter.enrich_claim`` one at a time).

## v0.6.3 — UniProt cross-refs for Drosophila / C. elegans / Zebrafish (2026-05-28)

### Added

- `build_gene_index.py` now fetches UniProt's **per-organism
  idmapping** files (gzip TSVs, ~5-12 MB each) for the v0.6.2
  species and joins them on the organism-specific id to populate
  the ``uniprot`` column. URLs:
  - DROME_7227_idmapping.dat.gz (Drosophila)
  - CAEEL_6239_idmapping.dat.gz (C. elegans)
  - DANRE_7955_idmapping.dat.gz (zebrafish)
- New helper ``_load_uniprot_by_organism(path, organism_type)``
  parses the three-column TSV and returns a
  ``{organism-id → UniProt accession}`` dict, filtered to the
  ``FlyBase`` / ``WormBase`` / ``ZFIN`` TYPE column.
- Per-cache CLI options (``--drosophila-uniprot-cache``,
  ``--celegans-uniprot-cache``, ``--zebrafish-uniprot-cache``)
  plus ``--no-other-species-uniprot`` to skip.
- Yields after the v0.6.3 rebuild:
  - Drosophila: 13,941 rows enriched with UniProt (~55.6%).
  - C. elegans: 19,974 rows (~42.6%).
  - Zebrafish:  20,983 rows (~38.5%).
  - The non-100% coverage reflects UniProt's catalog itself —
    not every annotated NCBI / FlyBase / WormBase / ZFIN gene
    has a UniProt entry (especially non-coding RNAs, pseudogenes,
    predicted ORFs).

### Changed

- ``_load_ncbi_gene_info`` now accepts ``Ensembl_Metazoa`` as
  well as ``Ensembl`` in the ``dbXrefs`` column, in anticipation
  of NCBI surfacing Ensembl Genomes cross-refs for invertebrate
  organisms. As of 2026-05-28 NCBI's gene_info files for
  Drosophila and C. elegans do **not** include any Ensembl /
  Ensembl_Metazoa entries (only FlyBase / WormBase + AllianceGenome
  + miRBase appear in dbXrefs), so the parser change is defensive
  and yields zero new Ensembl cross-refs on this run. Backfilling
  invertebrate Ensembl ids requires a separate source — Ensembl
  Genomes FTP files — scheduled for v0.6.4.

### Tests

- Two existing v0.6.2 pins extended: ``dpp`` now also asserts
  ``uniprot=P07713``; ``unc-22`` asserts ``uniprot=Q23551``.

### Notes

- Mouse UniProt coverage from v0.6.1 (~31%) could similarly be
  improved using UniProt's MOUSE_10090_idmapping.dat.gz — same
  pattern as the new species, swapping ``MGI`` as the TYPE column.
  Not done in v0.6.3 because the MRK_SwissProt_TrEMBL.rpt source
  already covers it; the UniProt mapping file would mostly add
  TrEMBL accessions and isoforms with diminishing returns.

## v0.6.2 — Drosophila / C. elegans / Zebrafish gene indexes (2026-05-28)

### Added — three model organisms

- `build_gene_index.py` now fetches NCBI's per-organism
  ``gene_info`` files from Entrez (gzip-compressed TSV, ~1-2 MB
  each) and joins each into the bundled `gene_index.tsv`:
  - **Drosophila melanogaster** — 25,079 rows from
    ``Drosophila_melanogaster.gene_info.gz``. 99.9% carry a
    FlyBase accession (``FBgn…``).
  - **Caenorhabditis elegans** — 46,927 rows from
    ``Caenorhabditis_elegans.gene_info.gz``. 99.9% carry a
    WormBase accession (``WBGene…``).
  - **Danio rerio (zebrafish)** — 54,563 rows from
    ``Danio_rerio.gene_info.gz``. ~50% carry a ZFIN accession,
    ~40% carry an Ensembl gene id (``ENSDARG…``).
- ``--no-other-species`` CLI flag skips the v0.6.2 step (regresses
  to v0.6.1 6-source build). Per-organism cache paths configurable.
- Bundled `gene_index.tsv` grows from 128k (v0.6.1) to **254,922 rows**.

### Added — schema + runtime extension

- Three new columns in the index schema: ``flybase``, ``wormbase``,
  ``zfin``. Empty for rows of other species; populated for the
  matching species rows.
- ``_KNOWN_SPECIES`` extended to ``("human", "mouse", "zebrafish",
  "drosophila", "celegans", "nematostella")``. The order encodes
  the load-priority fallback when neither hint nor heuristic
  disambiguates a symbol collision.
- New direct-accession lookup paths (``_by_flybase``, ``_by_wormbase``,
  ``_by_zfin``). All three accept both the bare id (``FBgn0000490``)
  and the prefixed form (``FlyBase:FBgn0000490``).
- ``_GENE_PRINCIPAL_CURIE_PREFIX`` priority extended to:
  HGNC > MGI > FlyBase > WormBase > ZFIN > Ensembl > NCBIGene.
- ``_GENE_XREF_COLUMNS`` extended likewise so xrefs surface the
  organism-specific id.
- New ``TRANSCRIPTOMICS_PREFIXES`` entries: ``FlyBase``,
  ``WormBase``, ``ZFIN`` (mapped to identifiers.org URLs).

### Changed — casing heuristic

- New C. elegans rule: 3-4 lowercase letters + hyphen +
  alphanumeric (``unc-22``, ``daf-16``, ``let-7``) → ``celegans``.
  Most distinctive pattern across organisms.
- **Lowercase-only token no longer maps to a species.** Pre-v0.6.2 a
  lowercase ``mmp9`` fell through to human via the cross-listed
  Nematostella row; v0.6.2 keeps that behaviour by treating
  all-lowercase as ambiguous (could be a casual human/mouse typo
  or a canonical zebrafish/drosophila symbol). Returns ``None`` so
  the resolver falls through to the load-order rule (human first).
  Zebrafish and Drosophila claims that aren't covered by the
  single-candidate shortcut must pass ``species="zebrafish"`` or
  ``species="drosophila"`` to disambiguate.

### Tests

- 6 new tests in ``test_bulk.py``:
  - drosophila-unique symbol (``dpp``) → FlyBase via single-candidate
  - ``Egfr`` with vs without ``species="drosophila"`` hint
  - C. elegans heuristic on ``unc-22`` → WormBase
  - zebrafish-unique symbol (``tncb``) → ZFIN
  - zebrafish via ``species="zebrafish"`` hint (``nanog`` collision)
  - direct-accession resolution for FlyBase / WormBase / ZFIN
    accessions in claim slots.

### Updated — demo

- 3 new claims in ``examples/enrich_warner2018.py`` covering dpp,
  unc-22, nanog+species. 18/18 claims enrich.

### Notes

- ``Ensembl`` is extracted from NCBI gene_info ``dbXrefs`` for
  zebrafish (~40% of rows) but the field is absent for Drosophila
  and C. elegans NCBI files (those use Ensembl Metazoa with a
  different dbXrefs prefix). Backfilling Ensembl for the two
  invertebrates is on a v0.6.3 follow-up.
- UniProt is not consistently surfaced in NCBI gene_info for any
  of the three organisms; ``uniprot`` is empty in v0.6.2 for all
  three. The UniProt id mapping files give a clean cross-ref but
  parsing them is another v0.6.3+ task.

## v0.6.1 — Mouse cross-references from MGI (Ensembl / Entrez / UniProt) (2026-05-28)

### Added

- `data/build_gene_index.py` now fetches three additional MGI
  cross-reference reports and joins them on the MGI accession:
  - `MRK_ENSEMBL.rpt` (16 MB) — mouse Ensembl gene id
    (``ENSMUSG…``). Yield: 63,787 mouse rows enriched (~85%).
  - `MGI_EntrezGene.rpt` (88 MB) — mouse NCBI Gene (Entrez) id.
    Filtered to ``Status="O" AND Marker Type="Gene"``. Yield: 45,496
    (~60%).
  - `MRK_SwissProt_TrEMBL.rpt` (2 MB) — mouse UniProt accession
    (first accession listed, may be either SwissProt-reviewed or
    TrEMBL). Yield: 23,534 (~31%) — UniProt is dominated by
    SwissProt coverage which is smaller than the full proteome.
- All three files are downloaded to the same cache directory as
  HGNC and MRK_List2 (``~/.cache/indicium-adapters-transcriptomics/``).
  ``--mgi-ensembl-cache`` / ``--mgi-entrez-cache`` /
  ``--mgi-uniprot-cache`` CLI options accept custom paths;
  ``--no-mgi-xrefs`` skips the v0.6.1 enrichment entirely
  (regression to the v0.6.0 minimal mouse rows).
- ``build_gene_index._load_mgi_xrefs`` reads the three header-less,
  positional MGI reports and merges them into a single
  ``{MGI accession → {ensembl, ncbigene, uniprot}}`` dict consumed
  by step 3 of ``build()``.

### Tests

- 1 new test pinning ``Mapk1``'s full xref map (MGI + Ensembl +
  Entrez + UniProt) — catches regressions in any of the three MGI
  report parsers.

### Notes

- Mouse rows now look semantically equivalent to human rows in
  terms of xref coverage. A typical mouse claim like
  ``Mapk1 differentially_expressed_in regenerating_tissue`` produces:
  ``ontology_terms[subject] = MGI:1346858`` +
  ``xrefs[subject] = {mgi=MGI:1346858, ensembl=ENSMUSG00000063358,
  ncbigene=26413, uniprot=Q3UF82}``.
- ``v0.6.0`` mouse rows that lacked these xrefs would now silently
  re-enrich on adapter import — no API change.

## v0.6.0 — Mouse (MGI) gene index + species-aware lookup (2026-05-28)

### Added — MGI mouse gene catalogue

- `data/build_gene_index.py` extended to fetch the **MGI marker list**
  (MRK_List2.rpt, ~78 MB). Filters to ``Marker Type == "Gene"`` to
  exclude BAC/YAC ends and syntenic markers. **75,250 mouse rows**
  added to the bundled index — every protein-coding (and most
  non-coding) mouse marker known to MGI.
- Mouse rows carry `species=mouse` and an `mgi=MGI:NNN` accession;
  ensembl / ncbigene / uniprot cross-refs are left empty in v0.6.0
  (those live in separate MGI report files — integrating them is on
  the v0.6.1 roadmap).
- Bundled `gene_index.tsv` grows from ~53k to **128,353 rows**:
  14,150 Nematostella-anchored + 38,953 pure-HGNC human + 75,250 MGI
  mouse.

### Added — species-aware lookup (B + C disambiguation)

- New schema columns ``mgi`` and ``species`` (``human`` / ``mouse`` /
  ``nematostella``). The species tag drives the v0.6.0 lookup rules.
- New optional `species` field on the claim dict (``"human"`` /
  ``"mouse"`` / ``"nematostella"``). The dispatcher forwards it to
  every resolver; the gene resolver uses it to disambiguate
  symbol-level collisions (``MAPK1`` ↔ ``Mapk1`` both lower-case
  to ``mapk1``).
- **Casing heuristic** for claims with no species hint:
  ``MAPK1`` (all uppercase) → human (HGNC); ``Mapk1`` (Title case)
  → mouse (MGI). This is the universal naming convention in the
  field and works without any LLM-side change.
- **Fallback priority** when neither hint nor heuristic disambiguates:
  human > mouse > nematostella (load order, encoded in
  ``_KNOWN_SPECIES``). Most transcriptomics claims are human-centric.
- **Cross-listed rows.** A Nematostella-anchored row that carries a
  human HGNC ortholog under ``hgnc`` is treated as "serving human"
  too. So ``MMP9`` (all-uppercase → heuristic says human) resolves
  to ``HGNC:7176`` via the Nematostella row's ortholog field, not
  to the mouse ``Mmp9`` row. Rule lives in
  ``GeneIndex._pick_for_species``.
- MGI accessions (``MGI:1346858``) in claim slots resolve directly
  via the new ``_by_mgi`` index — bypassing the symbol-collision
  logic entirely.

### Added — runtime + resolver-signature uniformity

- Principal CURIE priority extended: HGNC > **MGI** > Ensembl >
  NCBIGene. Mouse rows with empty ``hgnc`` naturally surface
  ``MGI:NNN`` as principal; Nematostella rows with a human ortholog
  still surface HGNC first.
- ``xrefs[slot]`` now also surfaces the ``mgi`` column for any row
  that carries it.
- All three resolvers (accession / gene / concept / anatomy) now
  share a uniform signature ``(name, *, species=None)`` so the
  dispatcher can forward the hint without branching. Non-gene
  resolvers explicitly ignore the hint (anatomy and concept entries
  are not species-specific in the current bundles).

### Added — tests

- 5 new tests in ``test_bulk.py``:
  - mouse via Title-case heuristic (``Mapk1`` → MGI:1346858)
  - human via UPPER-case heuristic (``MAPK1`` → HGNC:6871)
  - species hint overrides heuristic (``MAPK1``+``species=mouse``
    → MGI:1346858)
  - MGI accession resolves directly (``MGI:1346858`` → MGI:1346858)
  - mouse-only convenience pin (``Foxp3`` → MGI:1891436)

### Updated — demo

- ``examples/enrich_warner2018.py`` gains 3 v0.6.0 claims: Foxp3 →
  MGI + Treg → CL (cross-mechanism), Mapk1 → MGI via heuristic,
  MAPK1+species=mouse → MGI override. 15/15 claims enrich.

### Notes

- v0.6.1 follow-up — wire MGI Ensembl / MGI EntrezGene /
  MGI SwissProt cross-ref files into ``build_gene_index.py`` so
  mouse rows ship full xref maps (currently only ``mgi`` is
  populated).
- v0.6.2 follow-up — zebrafish (ZFIN), C. elegans (WormBase),
  Drosophila (FlyBase) following the same multi-species pattern.

## v0.5.0 — Cell Ontology enrichment + OLS4 cascade refinements (2026-05-28)

### Added — Cell Ontology (CL) harvest

- `data/build_cl_concept_index.py` maintainer script that resolves a
  curated list of 37 cell-type candidate names against the **Cell
  Ontology** via OLS4 search with the same cascade match modes used
  for the corpus harvest (exact label → exact synonym → prefix).
  Substring mode is omitted to avoid false positives — cell-type
  labels in CL are already minimal.
- Candidate list drawn from the three cell-atlas papers in the
  6-paper mini-corpus: **Sebé-Pedrós 2018** (*Nematostella vectensis*
  whole-organism cell atlas), **Tabula Muris 2018** (mouse), and
  **Tabula Sapiens 2022** (human). Plus canonical mammalian immune /
  nervous / epithelial cell types that recur in any transcriptomics
  study.
- Bundled `data/cl_concept_index.tsv` — **35 verified entries**
  (~95% yield), including: `cnidocyte` (CL:0011029, the cnidarian-
  specific cell type — CL covers it even though UBERON has no
  cnidarian anatomy), `T cell` / `B cell` / `macrophage` / `neutrophil`
  / `monocyte` (mammalian immune), `astrocyte` / `oligodendrocyte`
  (CNS), `fibroblast` / `endothelial cell` / `hepatocyte` (epithelial /
  connective), `enterocyte` / `goblet cell` / `Paneth cell` (GI), etc.
- `data/cl_concept_unresolved.txt` — 2 unresolved candidates
  (`gland cell`, `nematocyte` — the latter deduped against `cnidocyte`
  since both map to CL:0011029).
- Loader now reads four TSVs: static + OLS4 corpus + bio.tools + CL.
  Same priority ladder (first match wins) with CL placed last.

### Changed — OLS4 cascade matcher fixes

- **Label-first cascade pass in exact mode.** Previously the exact
  matcher returned the first OLS4 result whose label OR synonym
  matched the corpus term. OLS4's relevance ranking sometimes puts a
  synonym-match doc ahead of the label-match doc — e.g. searching
  ``macrophage`` against CL returns CL:0000394 plasmatocyte (which
  has macrophage as a synonym) before CL:0000235 macrophage (whose
  label *is* macrophage). The v0.5.0 cascade does label-exact across
  all docs first, then synonym-exact second; the canonical label-
  match wins.
- **Fix OLS4 search-endpoint synonym field name.** Search docs expose
  synonyms under ``exact_synonyms``, not ``synonym`` (the latter is
  used by the term-by-IRI endpoint). Pre-v0.5.0 the synonym pass was
  silently dead because the code looked at the wrong field; resolutions
  that should have used synonym matching fell through to prefix /
  substring matching instead. The fix is defensive — both
  ``exact_synonyms`` and ``synonym`` are now consulted.
- Net effect on the v0.3.1 corpus index (re-built with the fixes):
  one new resolution (``10x Genomics`` → ``EFO:0008995`` 10x technology)
  plus three canonicalisations:
  - ``pacbio`` from ``EFO:0920006`` PacBio Revio → ``NCIT:C146819``
    **SMRT Sequencing** (the technology, not a specific model).
  - ``oxford nanopore`` from ``OBI:0002750`` ONT MinION →
    ``NCIT:C146818`` **Nanopore Sequencing**.
  - ``bulk rna-seq`` from ``OBI:0003090`` bulk RNA-seq assay →
    ``NCIT:C219552`` **Bulk RNA Sequencing**.
  All three new mappings are semantically more canonical — the
  corpus term is the technology, not a specific instrument / assay
  realisation. Pinned tests updated accordingly. Corpus yield rises
  from 16 → 17 / 38 (~45%).

### Added — tests

- `test_cl_cell_types_resolve_via_concept_index` pins 12 CL
  resolutions including the macrophage-not-plasmatocyte assertion
  (proves the label-first refinement is in effect end-to-end).
- v0.4.0 test `test_singlecell_resolves_gene_subject_for_marker_of_cell_type`
  updated: pre-v0.5.0 it asserted ``cnidocyte`` stayed unenriched
  (CL not wired); v0.5.0 asserts the object slot resolves to
  ``CL:0011029``. The test now exercises both slots of the
  ``marker_of_cell_type`` qualifier — gene subject + CL object.
- `test_corpus_concept_terms_resolve_via_concept_index` updated for
  the v0.5.0 canonicalised mappings (pacbio / oxford nanopore /
  bulk rna-seq) plus the new ``10x genomics`` resolution. Pins the
  match-mode tag in the source column.
- Size invariant now ≥92 rows (was ≥56): 31 static + 17 OLS4-corpus
  + 9 bio.tools + 35 CL.

### Notes

- Single-cell now has **two-slot resolution** for all five of its
  cell-resolution qualifiers — both the gene subject (v0.4.0 wiring)
  and the cell-type object (v0.5.0 wiring) get CURIEs. The
  `marker_of_cell_type` qualifier is now fully expressive.
- The 24 corpus terms unresolved at v0.3.1 → 21 at v0.5.0. The four
  changes: +1 new (`10x Genomics`), the three canonicalisations don't
  change the unresolved count (they were already resolved, just to
  a different CURIE).
- Demo `examples/enrich_warner2018.py` updated to reflect the new
  v0.5.0 behaviour on claim 11 (cnidocyte resolves to CL:0011029
  instead of staying bare).

## v0.4.0 — qualifier-aware resolution priority (2026-05-27)

### Changed

- `BulkTranscriptomicsAdapter.enrich_claim()` is now **qualifier-aware**.
  Previously the resolver chain was hardcoded
  ``accession`` → ``gene`` → ``concept`` for all qualifiers; this always
  resolved short ambiguous names like ``STAR`` to the HGNC gene
  (HGNC:11359, steroidogenic acute regulatory protein) even when the
  claim's qualifier made it obvious the term was the aligner.
- Two priority orders are now defined and picked per claim:
  - **Default** (``accession`` → ``gene`` → ``concept``) — used for the
    claim-level feature qualifiers (``differentially_expressed_in``,
    ``co_expressed_with``, ``enriched_for``, ``expressed_in``,
    ``marker_of``, ``regulated_by``) and for ``consolidates_dataset``.
    Gene mentions dominate these claims, so gene-lookup wins for
    collision-prone names.
  - **Study-level** (``concept`` → ``accession`` → ``gene``) — used for
    the six study-metadata qualifiers ``sequenced_on``,
    ``aligned_with``, ``quantified_with``, ``normalized_with``,
    ``tested_with``, ``uses_database``. Object slot of these is a tool
    / platform / database name, so the concept index (including
    bio.tools) wins; the gene resolver only fires as a fallback.
- Internally, the three previously-inlined resolvers are now
  first-class functions (``_resolve_accession``, ``_resolve_gene``,
  ``_resolve_concept``) dispatched through a ``_RESOLVERS`` map. Each
  returns ``(curie, xrefs_dict) | None`` so the resolver-chain is now
  a simple loop. Same shape, half the code.

### Added

- 4 regression tests pinning the v0.4.0 contract:
  - ``aligned_with: <study> STAR`` → ``biotools:star`` (was HGNC:11359
    pre-v0.4.0).
  - ``expressed_in: STAR <something>`` → ``HGNC:11359`` (gene-priority
    preserved for feature qualifiers).
  - Missing or unknown qualifier → default priority (safe fallback).
  - Study-level priority still falls through to gene-lookup when the
    concept index has no hit, so a ``aligned_with`` claim whose object
    is a gene name (sloppy LLM output) still gets a CURIE rather than
    a dropped slot.

### Changed — SingleCellTranscriptomicsAdapter (symmetric refactor)

- ``SingleCellTranscriptomicsAdapter`` now uses the same
  resolver-dispatcher pattern and **gains a gene resolver**. Previously
  single-cell had only anatomy + concept lookup wired, so a slot like
  ``marker_of_cell_type: WNT4 cnidocyte`` left ``WNT4`` unenriched
  even though it's a perfectly valid HGNC symbol. Now it resolves to
  ``HGNC:12783``.
- Single-cell uses one default priority order (``gene → anatomy →
  concept``) for all its qualifiers — they're all claim-level
  (cell-resolution feature assertions), so no study-level override is
  needed. Anatomy still wins for tissue terms (gene-lookup of
  ``mesentery`` returns ``None``, the cascade falls through to
  anatomy).
- 2 new tests pin the contract: gene resolution for the cell-type
  qualifiers, and the gene-doesn't-override-anatomy fallback.

### Added — shared resolver helpers

- The three resolvers (gene / anatomy / concept) moved into
  ``_enrichment.py`` as module-level callables (``resolve_gene``,
  ``resolve_anatomy``, ``resolve_concept``), keeping the row → CURIE
  rendering logic next to the corresponding Index classes. Bulk and
  single-cell both import and dispatch through them — one source of
  truth for the rendering rules.

### Added — local CI parity helper

- ``scripts/check.ps1`` runs the same checks the GitHub Actions CI
  runs (ruff + pytest with coverage). Use it before every push;
  pytest alone misses lint failures. Direct response to the v0.3.1
  CI red caused by a ruff E741 that didn't show up locally.

### Notes

- The 13 corpus terms unresolved at v0.3.1 (Cell Ranger, scVI, FLAIR,
  TALON, IsoQuant, SQANTI, StringTie2, Bambu, FLAMES, LRGASP, 10x
  Genomics, Iso-Seq, ONT direct RNA) are unaffected by v0.4.0 — they
  were absent from both OLS4 and bio.tools, not blocked by priority
  ordering.
- A monkeypatch test that previously patched ``bulk.CONCEPT_INDEX`` /
  ``singlecell.CONCEPT_INDEX`` directly was updated to patch
  ``_enrichment.CONCEPT_INDEX`` instead — since v0.4.0, both adapters
  read the index through the shared ``resolve_concept`` helper and
  the singleton's home is ``_enrichment``.

## v0.3.1 — corpus-driven concept harvest: OLS4 (cascade modes) + bio.tools (2026-05-27)

### Added — OLS4 corpus harvester with cascade match modes

- `data/build_corpus_concept_candidates.py` maintainer script that
  ingests the per-paper vocabulary JSONs produced by
  `articles/_extract_vocab.py` over the 6-paper mini-corpus and uses
  the **OLS4 search-by-label** endpoint to propose CURIEs, then
  re-verifies each by IRI lookup. Hand-curation of CURIEs is avoided;
  the only curation is the axis → ontology preference list and the
  axis → match-mode declaration.
- Per-axis **cascade of match modes** tried in order; the first hit
  wins. Strictness beats preference: an exact match in any preferred
  ontology beats a relaxed match in the first preferred ontology.
  - `exact` — label or synonym equals the corpus term verbatim
    (case-insensitive). Default for all software axes and for high-
    precision lookups.
  - `prefix` — corpus term is a whole-word **prefix** of the OLS4
    label (lets descriptive suffixes through: "Wilcoxon" → "Wilcoxon
    Signed Rank Test"; "bulk RNA-seq" → "bulk RNA-seq assay";
    "Visium" → "Visium Spatial Gene Expression"). Used by
    `protocol` and `stat_test` axes.
  - `substring` — corpus term is a whole-word **substring** of the
    OLS4 label (handles instrument families whose canonical label
    prepends a vendor name: "NovaSeq" → "Illumina NovaSeq 6000").
    Used by the `platform` axis only; OLS4 ranking picks the
    canonical model.
- Bundled `data/corpus_concept_index.tsv` — **16 verified entries**:
  - exact-match (6): `Illumina` (NCIT:C215820), `MARS-seq`
    (EFO:0008796), `Smart-seq2` (EFO:0008931), `Drop-seq`
    (EFO:0008722), `Slide-seq` (EFO:0009920), `t-test` (NCIT:C53231,
    via synonym).
  - prefix-match (6): `PacBio` (EFO:0920006 PacBio Revio),
    `Oxford Nanopore` (OBI:0002750 Oxford Nanopore MinION),
    `bulk RNA-seq` (OBI:0003090), `Visium` (EFO:0010961),
    `Wilcoxon` (NCIT:C53246), `Benjamini-Hochberg` (NCIT:C61596).
  - substring-match (4): `NextSeq` (EFO:0008566 Illumina NextSeq 550),
    `NovaSeq` (EFO:0008637 Illumina NovaSeq 6000), `HiSeq`
    (EFO:0004203 Illumina HiSeq 2000), `Sequel` (EFO:0008630
    PacBio Sequel system).
- The match mode is tagged in the TSV's `source` column
  (`corpus extraction (<axis>, <mode>-match)`) so the provenance of
  each row is auditable.

### Added — bio.tools registry harvester (with name normalisation)

- `data/build_biotools_index.py` maintainer script picking up where
  OLS4 left off: for software-axis terms still unresolved, queries
  the **bio.tools registry**, the canonical home for bioinformatics
  tools. Two-step resolution per term:
  1. Direct API lookup `GET https://bio.tools/api/<id>/?format=json`
     with a normalised ID (lowercased, spaces → hyphens).
  2. If 404, search endpoint with paginated results, filtered
     client-side for **exact name match (case-insensitive)** — and
     additionally for a **normalised-name match** that drops
     whitespace / hyphens / underscores so corpus `Bowtie2` matches
     the canonical bio.tools name `Bowtie 2`. Without this
     normalisation step the OLS4 harvester's substring trick has no
     bio.tools equivalent.
- Bundled `data/biotools_concept_index.tsv` — **9 verified entries**:
  `STAR` (biotools:star), `Bowtie` (biotools:bowtie), **`Bowtie2`**
  (biotools:bowtie2, via name normalisation), `Minimap2`
  (biotools:minimap2), `RSEM` (biotools:rsem), `StringTie`
  (biotools:stringtie), `DESeq2` (biotools:deseq2), `Seurat`
  (biotools:seurat), `Picard` (biotools:picard_tools).
- New `biotools:` prefix in `TRANSCRIPTOMICS_PREFIXES` (resolves to
  `https://bio.tools/<id>`).
- `data/biotools_unresolved.txt` — final log of the 13 corpus terms
  not findable in OLS4 nor bio.tools. Almost all are long-read tools
  (FLAIR, TALON, IsoQuant, SQANTI, StringTie2, Bambu, FLAMES, LRGASP)
  that genuinely do not exist in either registry — they live in
  their authors' GitHub / Zenodo only.

### Added — runtime + tests

- **Triple-TSV loader** in `_enrichment._load_concept_index_rows()`:
  loads `concept_index.tsv` → `corpus_concept_index.tsv` →
  `biotools_concept_index.tsv` in that order. A hypothetical
  collision resolves to the first occurrence: paper-modelling-
  grounded (static) wins over OLS4 method-mention wins over
  bio.tools-registry.
- 6 regression tests covering: pinned resolution for each of the 16
  OLS4-corpus terms (with the match-mode tag asserted in the source
  column), pinned resolution for each of the 9 bio.tools terms, two
  end-to-end `enrich_claim()` paths (one OLS4-resolved platform via
  `Illumina`, one bio.tools-resolved aligner via `minimap2` — `STAR`
  is deliberately avoided because it also matches HGNC:11359
  steroidogenic acute regulatory protein, which gene-lookup correctly
  resolves first), and the triple-TSV size lower bound (≥56 rows).

### Notes

- Total v0.3.1 yield: **25/38 corpus terms resolved** (~66%) —
  16 via OLS4 (6 exact + 6 prefix + 4 substring) + 9 via bio.tools.
- The 13 still-unresolved terms fall into one clean bucket: niche
  long-read tools that are genuinely outside both OLS4 (no ontology
  covers them) and bio.tools (not registered). Resolution would
  require either upstream registry contributions or per-tool
  GitHub/Zenodo CURIE schemes — neither is on the v0.3.x roadmap.
- A surprise found while writing the tests: some short bio.tools
  software names (`STAR`, `Picard`) **also** exist as HGNC gene
  symbols. The current resolver runs gene-lookup before
  concept-lookup, so a bare `STAR` in a claim resolves to HGNC:11359.
  This is correct default behaviour (gene mentions are more frequent
  in transcriptomics claims), but a qualifier-aware resolver
  (`aligned_with` → prefer bio.tools; `expressed_in` → prefer HGNC)
  would be a useful v0.4.0 refinement.

## v0.3.0 — TranscriptoKG concept-index enrichment (NCIT / INO / OBI / STATO / GO) (2026-05-27)

### Added

- New ontology axes wired into `TRANSCRIPTOMICS_PREFIXES`: **INO**
  (regulation processes), **NCIT** (analysis concepts), **SIO**, **RO** and
  **IAO** — drawn from the TranscriptoKG vocabulary work so adapters and
  the upstream KG stay aligned on the same CURIE namespace set.
- `data/build_concept_index.py` maintainer script that verifies a curated
  list of 31 CURIEs against EBI OLS4 and writes a TSV with the canonical
  label and definition. The candidate set is drawn from the TranscriptoKG
  mapping & vocab (Boukhajou, 2026): 2 INO regulation terms, 20 NCIT
  analysis concepts (Fold Change, log2, Timepoint, Multiple Comparison
  Adjustment, RNA Sequencing Data File, ...), 3 GO roots, 5 OBI
  investigation / role / protocol / DE-transformation / RNA-sequencing
  terms, and 1 STATO statistical-variable anchor. The build retries on
  OLS4 5xx and throttles 0.3s between requests.
- Bundled `data/concept_index.tsv` — **31 verified entries**, every one
  re-checked against OLS4 at build time (no hand-curation in the shipped
  row). Two OBI entries (`OBI:0000249`, `OBI:0000650`) had been tagged
  incorrectly in the upstream TranscriptoKG vocab — OLS4 verification
  caught the drift and the source notes were corrected to reflect the
  canonical labels ("technical replicate role" and "differential
  expression analysis data transformation"); the terms themselves stay
  relevant to the transcriptomics domain so the rows are kept.
- `_enrichment.ConceptIndex` real lookup class loading the bundled TSV
  at import time. Accepts a term (case-insensitive) and returns the full
  row dict (CURIE / label / description / ontology) or `None`. Strict
  header check at load time: a malformed TSV raises before silent
  mis-mapping can occur.
- Concept-lookup branch wired into both adapters' `enrich_claim()` — runs
  after gene / accession / anatomy lookup, so domain-specific resolution
  always takes precedence. Adds the principal CURIE under
  `ontology_terms[slot]` and `{curie, label, ontology, description}`
  under `xrefs[slot]`.
- 7 regression tests covering: index loadability, `None` for unknown /
  empty / non-str input, safe-no-op for unknown terms, end-to-end CURIE +
  xref attachment, single-cell-side wiring via monkeypatched index.

### Changed

- Test docstrings that referenced "while the TSV is empty (pending OLS4
  verification)" updated to reflect the populated v0.3.0 state.
- `BulkTranscriptomicsAdapter` and `SingleCellTranscriptomicsAdapter` now
  expose 17 / 12 effective lookup axes respectively (gene + accession +
  anatomy + concept).

### Notes

- A second concept-harvesting pass (corpus-driven: tools / platforms /
  protocols extracted from the 6-paper mini-corpus — STAR, Seurat,
  MARS-seq, Illumina, 10x, etc.) is scheduled for v0.3.1 / v0.4.0 once a
  `build_corpus_concept_candidates.py` script using OLS4
  search-by-label is written. Hand-curation of CURIEs is explicitly
  avoided.
- The earlier ROADMAP note about contributing cnidarian-specific anatomy
  upstream to UBERON is unchanged — still on the v0.4.0 wish-list.

## v0.2.0 — UBERON anatomy enrichment for SingleCellTranscriptomicsAdapter (2026-05-26)

### Added

- `data/build_anatomy_index.py` maintainer script that resolves a curated
  list of anatomy terms (drawn from the Cheung et al. 2024 Tomo-seq paper
  on *Nematostella vectensis*) against the EBI OLS4 UBERON index. Only
  entries whose returned label exactly matches the query are written —
  hyponym promotion (e.g. *"nasal tentacle"* for query *"tentacle"*)
  cannot slip through.
- Bundled `data/anatomy_index.tsv` (4 verified UBERON anchors for v0.2.0):
  `mesentery` (UBERON:0002095), `pharynx` (UBERON:0006562), `body wall`
  (UBERON:0000309), `organism subdivision` (UBERON:0000475). The bundle
  is deliberately small — UBERON has no entry for cnidarian-specific
  concepts (physa, mesoglea, gastrodermis, polyp body form, oral / aboral
  pole, body column), so the build script documents them in
  `TRIAGED_OUT_KNOWN_CNIDARIAN` rather than emitting wrong CURIEs.
- `_enrichment.AnatomyIndex` real lookup class loading the bundled TSV at
  import time. Accepts an anatomical term (case-insensitive) and also
  matches against the UBERON `label` field; returns the full row dict
  or `None`.
- `SingleCellTranscriptomicsAdapter.enrich_claim()` now resolves anatomy
  terms in `subject` / `object` to their UBERON CURIE under
  `ontology_terms[slot]` and attaches the label + description under
  `xrefs[slot]` so consumers can render without an OLS4 round-trip.
- 7 regression tests covering: case-insensitive anatomy lookup, UBERON
  CURIE attachment, xref shape, both subject / object slots resolved,
  cnidarian-specific terms left untouched (no invented CURIEs).

### Notes

- Cell-type enrichment via the Cell Ontology (CL) for the cell-type
  qualifiers (`marker_of_cell_type`, `enriched_in_cluster`,
  `co_expressed_in_cell_type`) is scheduled for v0.2.1 once an
  scRNA-seq reference paper is selected.
- The four cnidarian-specific terms triaged out in v0.2.0 are candidates
  for a v0.3.0 effort to contribute them upstream to UBERON, or to allow
  alternative ontology prefixes (e.g. CARO, GO for ECM) in the anatomy
  index.

## v0.1.2 — comprehensive HGNC + Nematostella gene index (2026-05-26)

### Added

- `data/build_gene_index.py` maintainer script that regenerates the bundled
  TSV from two authoritative sources: the **HGNC complete set** (downloaded
  with a fallback mirror) and the **G2V2 Nematostella cross-reference table**
  (Nv2g / Nv2t / JGI / NvERTx / NVE alignment with UniProt-derived
  orthology, GO and KEGG annotations). The script accepts both sources as
  local paths; the HGNC TSV is cached under `~/.cache/`.
- Two-pass join strategy:
  1. Nematostella-anchored rows from G2V2 — HGNC ortholog attached via
     UniProt accession when it matches HGNC, otherwise via a
     case-insensitive `Preferred_name → HGNC symbol` fallback.
  2. Pure-HGNC rows for every human gene that was *not* already paired with
     a Nematostella entry, so plain human-gene mentions (e.g. an MMP in a
     cnidarian-vs-mammal comparison) still resolve to their HGNC CURIE.
- Build-time-only dependencies (`pandas`, `openpyxl`, `requests`) declared
  under `[project.optional-dependencies] build`; the runtime keeps its
  minimal closure (indicium, indicium-adapters, rdflib, pyshacl).

### Changed

- Bundled `data/gene_index.tsv` rebuilt: **~53,000 entries** (up from 22
  hand-curated). ~7,100 Nematostella-anchored rows with HGNC ortholog
  identified, ~7,000 Nematostella-only, ~39,000 pure-HGNC human.
- Behaviour-based test refactor: gene-index tests now assert lookup
  semantics (case-insensitive symbol resolution, principal-CURIE order,
  xref attachment, both lookup paths exercised) rather than specific
  hand-curated entries. Two new regression tests cover the v0.1.2-specific
  capabilities (Nemve1 → HGNC via symbol fallback; pure-HGNC coverage for
  human-only mentions).

## v0.1.1 — bundled gene-name index (2026-05-26)

### Added

- Bundled `data/gene_index.tsv` — 22 hand-curated entries drawn from the two
  benchmark reference papers (Warner et al. 2025 supplementary DEG tables;
  Cheung et al. 2024 main text). Covers core regeneration signalling (WNT4,
  FGF8, BMP1, SPRY2, MEOX1, HMX2), MEK/ERK pathway (MAPK1, MAPK3, MAP2K1),
  cnidarian matrix and ECM (MMP9/14/17, TIMP2, COL4A1, COL24A1, VCAN, MFAP4),
  β-catenin (CTNNB1), Brachyury (TBXT) and a handful of transcription
  factors and adhesion proteins. Each row carries HGNC, Ensembl, NCBI Gene,
  Nemve1 / JGI, NvERTx and UniProt cross-references plus a description.
- `_enrichment.GeneIndex` real lookup class loading the bundled TSV at
  import time. Accepts gene symbols (case-insensitive), Nemve1 / JGI ids,
  and UniProt accessions; returns the full row dict or `None`.
- `BulkTranscriptomicsAdapter.enrich_claim()` now resolves gene names in
  `subject` / `object` to a principal CURIE (HGNC > Ensembl > NCBIGene)
  and attaches the full xref map under `xrefs[slot]`. Repository accession
  enrichment (v0.1.0) is unchanged and still takes priority over the gene
  lookup.
- Regression tests covering all three lookup paths (symbol case-insensitive,
  Nemve1 / JGI, UniProt), the principal-CURIE order, and the xref
  attachment. Coverage: 95.17%.

## v0.1.0 — initial scaffold (2026-05-26)

### Added

- **`BulkTranscriptomicsAdapter`** (api_version=1) with 13 transcriptomics
  qualifiers (claim-level: `differentially_expressed_in`, `co_expressed_with`,
  `enriched_for`, `expressed_in`, `marker_of`, `regulated_by`; study-level:
  `sequenced_on`, `aligned_with`, `quantified_with`, `normalized_with`,
  `tested_with`, `uses_database`, `consolidates_dataset`).
- **`SingleCellTranscriptomicsAdapter`** (api_version=1) — structural stub
  with 5 cell-resolution qualifiers (`marker_of_cell_type`,
  `enriched_in_cluster`, `co_expressed_in_cell_type`, `localized_in_tissue`,
  `localized_in_spatial_domain`). `enrich_claim()` is a no-op pass-through in
  v0.1.0; real CL + UBERON lookup arrives in v0.2.0.
- `TRANSCRIPTOMICS_PREFIXES` IRI map (EFO, OBI, STATO, NCBITaxon, CL, UBERON,
  GO, SO, MONDO, PATO, EDAM plus HGNC, MGI, Ensembl, NCBI Gene, RefSeq,
  UniProt and the 8 major transcriptomics repository accession prefixes),
  shared by both adapters.
- `BulkTranscriptomicsAdapter.enrich_claim()` resolves repository accessions
  to CURIEs (GEO, SRA, ENA, ArrayExpress, ENCODE, BioProject, dbGaP,
  BioStudies) via the `accession_to_curie` helper. Gene-symbol enrichment is
  wired but no-op (real bundled index lands in v0.1.1).
- `shacl_shapes()` on both adapters returns a minimal SHACL graph validating
  CURIE format on ASB claim graphs.
- Entry-point registration under `indicium.adapters` for both
  `bulk-transcriptomics` and `single-cell-transcriptomics`.
- Pytest suite covering Protocol conformance, attribute presence, enrichment
  behaviour, accession-pattern table and the cross-adapter disjoint-qualifier
  invariant (`tests/test_composition.py`).
