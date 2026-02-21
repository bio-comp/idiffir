PyML to scikit-learn Migration Design
=====================================

Summary
-------
This document defines the migration contract for SpliceGrapher classifiers currently bound to
PyML ``.svm`` assets. The migration target is scikit-learn with serialized ``.joblib`` model
artifacts and explicit metadata to preserve feature parity.

Scope for this phase:

1. lock the model API contract used by callers
2. define feature mapping parity requirements
3. define artifact format and metadata policy
4. define compatibility/deprecation sequencing for legacy ``.svm`` assets

Target Model API Contract
-------------------------
Model loading and inference paths should converge on a single runtime surface:

1. ``load_model(path, metadata_path) -> model_handle``
2. ``predict_scores(model_handle, features) -> ndarray[float]``
3. ``predict_labels(model_handle, features, threshold) -> ndarray[int]``

Contract requirements:

1. binary labels are ``0/1`` and threshold defaults stay compatible with legacy classifier config
2. both site-classification (acceptor/donor) and gap-classification tasks share the same
   metadata-driven load path
3. feature extraction remains outside the model object and must be versioned in metadata

Feature Extraction Parity Mapping
---------------------------------
Parity is defined by named feature mappings, not by positional assumptions.

Each migrated model metadata file must include:

1. legacy feature identifier (from existing config or code path)
2. new scikit-learn feature name emitted by extractor
3. optional notes when one legacy feature expands into multiple normalized features

For site classifiers:

1. preserve k-mer window semantics from legacy config (``mink``, ``maxk``, mismatch profile)
2. preserve sequence window boundaries (intron/exon offsets)

For gap classifiers:

1. preserve legacy feature names from ``GapClassifier.py`` constants
2. preserve feature ordering through explicit metadata to avoid hidden regressions

Serialized Model Format Policy
------------------------------
Canonical persisted format is:

1. model: ``*.joblib``
2. metadata: ``*.json`` matching ``ModelMetadata`` schema in
   ``iDiffIR.SpliceGrapher.predict.model_schema``

Metadata requirements:

1. ``schema_version`` (currently ``1``)
2. ``classifier_task`` (site_acceptor, site_donor, gap)
3. ``artifact_format`` (joblib)
4. ``artifact_path``
5. ``sklearn_version``
6. ``feature_extractor``
7. ``feature_mappings`` list
8. optional legacy pointers (``legacy_config_path``, ``legacy_svm_path``)

Compatibility and Deprecation Plan
----------------------------------
Phase A (current):

1. keep legacy ``.svm`` assets present
2. add metadata schema and migration design artifacts
3. implement migration tooling in a follow-up issue

Phase B (dual-read):

1. default to ``.joblib`` when metadata exists
2. allow explicit fallback to legacy ``.svm`` path for unmigrated species bundles
3. emit one warning per process when fallback path is used

Phase C (deprecation enforcement):

1. add release note announcing legacy fallback removal timeline
2. fail closed for missing metadata in classifier bundles designated as migrated
3. remove runtime loading from PyML in a dedicated cleanup issue

Validation Criteria for Migration PRs
-------------------------------------
Each migration PR should satisfy:

1. task-level parity tests for feature extraction outputs
2. classifier loading tests using metadata schema round-trip
3. deterministic inference checks for fixed synthetic inputs
4. no user-facing CLI flag changes
