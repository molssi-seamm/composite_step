=======
History
=======

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
