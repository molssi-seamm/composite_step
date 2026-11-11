Phase 1 notes
=============

Done
----

- ``composite_step.recipes``: the YAML format and its checks; the ``combine``
  language, parsed with :mod:`ast` and evaluated over a whitelist (never
  ``eval``) with ``cbs_power``, ``hf_extrapolate_exp`` (G4 and Karton-Martin
  form), ``schwenke``, ``scale``, ``hf_part`` and ``corr_part``; the library of
  shipped, user and job recipes, in which a local recipe may not shadow a shipped
  name and version; the capability check; fingerprints (the descriptive fields
  left out); and the thermal terms, H0 = E0 + s ZPE and H = E0 + s ZPE +
  [H - E - ZPE] at the frequency level.
- The recipes ``thermo-quick`` (r2SCAN-3c) and ``thermo-default`` (r2SCAN-3c //
  revDSD-PBEP86-D4(2021)/def2-TZVPP), ZPE x 0.9899 pinned with its source. Tests
  check that the factor is still seamm_thermochemistry's, and that the recipes
  reproduce the thermochemistry step's arithmetic bit for bit. The names may
  change (science suggests naming by content: ``r2SCAN-3c-HO``,
  ``revDSD-r2SCAN-3c-HO``); Paul decides.
- ``CompositeStep.get_model_chemistry_options``: ``SEAMM:CMP@<recipe>`` for the
  newest version of each recipe, with its geometry level, components and ZPE
  factor in the option. Whether the components' codes are installed is not
  checked there (asking model-chemistry discovery from inside it would recurse);
  the step checks before it runs, with a matcher that lets the basis be any, as
  ``match_model_chemistry`` does.

The task route (next)
---------------------

seamm_exec's Evaluator takes one task per structure from one provider, and a
composite is several single points per structure, each on its component's own
code. For phase 1 the components run as one Evaluator per component, over all the
structures at once, through the component's own model chemistry; the composite
then combines the energies per structure. This needs no change to seamm_exec, and
it keeps every component an ordinary calculation with its ordinary settings
fingerprint, so the fragment store finds and stores components, never composites
(design §6), and a recipe's component of ``thermo-default`` is exactly the
default preset's energy level (science's pass criterion (a)).

A later change can let the Evaluator take several tasks per structure from one
provider (``get_tasks`` / ``analyze_tasks``), so that every consumer of a model
chemistry gets composites without knowing about them, and all components run in
one task set. That is a seamm_exec change for its own PR.

Pilot (science's criteria, in order)
------------------------------------

(a) The recipe's component settings fingerprints equal those the preset path
    produces, so the same store rows are hit.
(b) H0 and H(298.15 K) reproduce the step's to round-off on the 20 demo molecules.
(c) Three to five molecules computed fresh through the task route, end to end.

Only today's presets (from 2026-11-08: VERYTIGHTSCF, DEFGRID3, NoCOSX rules) are
a valid reference.
