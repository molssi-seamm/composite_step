=======
History
=======
2026.11.12.3 -- Quasi-restricted orbitals for open-shell species
    * Recipes can say which reference their correlated components use for
      open-shell species, ``open_shell_reference: QRO`` (quasi-restricted orbitals)
      or ``UHF``, and a component can override it. Both halves of a difference
      must use the same reference.
    * The three tier-A recipes use QRO. Closed-shell species are not affected: their
      calculations, fingerprints and stored results are unchanged.
    * Needs orca-step 2026.11.12.2 or later for QRO.

2026.11.12.2 -- Core-valence and relativistic corrections, and the tier-A recipes
    * Recipe components can be given their frozen core (named sets of frozen
      orbitals per element), a scalar-relativistic Hamiltonian with its nucleus
      model, and a decontracted basis, so recipes can add core-valence and
      scalar-relativistic corrections. A difference of two components must change
      one thing only, or the recipe is refused.
    * Three shipped recipes give accurate reference values after W4: CCSD(T)-F12
      (H-Ar), CCSD(T) extrapolated from triple and quadruple zeta (H-Ar, Ga-Br),
      and the same with a decontracted relativistic correction for species with
      calcium. Their core-valence and relativistic corrections agree with W4's to
      0.13 and 0.01 kJ/mol on test reactions.
    * Components of one model chemistry with different settings are kept apart, in
      the job and the fragment store.
    * Needs orca-step 2026.11.12.1 or later for these settings.

2026.11.12.1 -- Recipes that use the parts of an energy
    * Recipes can use the parts of a component's energy: ``part(c, "name")`` for
      the HF energy, the CABS singles, the correlation energies and the triples that
      the program reports, and parts derived from them (the CCSD-F12 correlation, the
      scaled triples (T*), ...). ``hf_part(c)`` is the HF energy and ``corr_part(c)``
      the correlation alone, so a recipe can, for example, extrapolate the
      correlation energy and not the HF energy.
    * A recipe that splits an F12 calculation into HF and correlation must account
      for its CABS singles, and the HF parts of components with different references
      cannot be combined unless the recipe allows it; such recipes are refused when
      they are read.
    * A stored energy without the parts a recipe needs is computed again; parts are
      never derived.

2026.11.12 -- Recipes that combine several calculations run in every step
    * A composite recipe that combines several calculations now runs wherever a
      model chemistry is used: the Energy step, the Thermochemistry step's CBH
      approach, and the others that run model chemistries through seamm-exec's
      Evaluator. Each component runs as one batch of its own program over all the
      structures, the cheapest first, with that program's own memory and cores.
    * A structure that fails one component is not sent to the later ones, and its
      result names the component that failed.
    * Every result carries the composite record: the recipe's fingerprint and every
      component's settings fingerprint, in the same form whether the recipe has one
      component or several, and whether a component was computed or found in a
      store.
    * A user's or a job's own recipe (in a ``recipes`` directory of SEAMM's data
      path) is found wherever its model chemistry runs.
    * Needs seamm-exec 2026.11.12 or later for recipes with several components.

2026.11.11 -- First release: composite methods as model chemistries
    * Composite methods are described by recipes, short YAML files that say how the
      energy is made from ordinary calculations: the geometry and frequency level,
      single points at that geometry, and a formula that combines them, with
      basis-set extrapolations, scaling and empirical terms.
    * Every recipe is a model chemistry, ``SEAMM:CMP@<recipe>``, for the steps that
      take one. Its components are ordinary calculations with their own settings
      fingerprints, so the fragment store reuses them.
    * The recipes ``r2SCAN-3c-HO`` and ``revDSD-r2SCAN-3c-HO`` are shipped; they are
      the Thermochemistry step's ``quick`` and ``default`` presets.
    * The Composite step itself does nothing yet and is not offered in the editor.
