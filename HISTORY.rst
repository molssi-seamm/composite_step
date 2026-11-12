=======
History
=======
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
