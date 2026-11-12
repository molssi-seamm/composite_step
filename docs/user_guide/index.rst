.. _user-guide:

**********
User Guide
**********

A composite method makes an energy from several ordinary calculations: a geometry
and frequencies at one level, single points at others at that geometry, and a
formula that combines them (extrapolations to the basis-set limit, additive
corrections, empirical terms). The Composite plug-in describes each method as a
*recipe*, a short YAML file, and offers every recipe as a model chemistry,
``SEAMM:CMP@<recipe>``.

Using a composite method
========================

Composite methods are used through a **Model Chemistry** step, like any other
model chemistry, by the steps that take one. For example::

    SEAMM:CMP@revDSD-r2SCAN-3c-HO

The **Thermochemistry** step's CBH presets are recipes: ``quick`` is
``r2SCAN-3c-HO`` and ``default`` is ``revDSD-r2SCAN-3c-HO``, and the step records
the recipe, its version and its fingerprint with the results.

The components of a recipe are ordinary calculations, each with its own model
chemistry and settings fingerprint. A component already in the fragment store is
reused, and one shared by two recipes is computed once; the store keeps
components, never composites.

The Composite step itself, in a flowchart, does nothing yet.

The shipped recipes
===================

``r2SCAN-3c-HO``
    r2SCAN-3c geometry, frequencies and energy (ORCA), with the zero-point energy
    scaled by 0.9899 and scaled harmonic thermal terms ("HO").

``revDSD-r2SCAN-3c-HO``
    r2SCAN-3c geometry and frequencies, as above, and the
    revDSD-PBEP86-D4(2021)/def2-TZVPP energy at that geometry.

Recipes are named by what they contain, not by the role they play in another step,
so a result labelled with a recipe stays true when presets change.

Writing a recipe
================

A recipe is a YAML file::

    recipe: revDSD-r2SCAN-3c-HO
    version: 1
    description: >-
      r2SCAN-3c geometry and frequencies (ZPE scaled), and the
      revDSD-PBEP86-D4(2021)/def2-TZVPP energy at that geometry.
    reference: Santra, Sylvetsky, Martin, J. Phys. Chem. A 123, 5129 (2019)

    geometry:
      level: ORCA:CMP@R2SCAN-3C
    frequencies:
      scale_zpe: 0.9899
      scale_source: where the factor comes from

    components:
      revdsd:
        level: ORCA:DFT@REVDSD-PBEP86-D4_2021/def2-TZVPP

    combine: |
      E0 = revdsd

``recipe``, ``version``, ``geometry`` and ``combine`` are required. The name may not
contain ``:``, ``@``, ``/`` or ``|``, which separate the parts of a model
chemistry.

``geometry``
    The level of the geometry and frequencies, as a model chemistry.

``frequencies``
    ``scale_zpe``, the factor for the zero-point energy and the harmonic thermal
    terms, and ``scale_source``, where it comes from.

``components``
    The single points at the recipe's geometry, by name. Each gives its ``level``,
    a complete model chemistry. The code of every component must be installed;
    any basis set may be given for a method the code offers.

``combine``
    Assignments, one per line, in a small language: numbers, the names of the
    components, ``geometry`` (the geometry level's own energy), the values
    assigned on earlier lines, the constants under ``empirical``, ``+ - * /``,
    parentheses and these functions:

    * ``cbs_power(e_x, e_y, x, y, alpha=3)`` -- two-point inverse-power
      extrapolation from the cardinal numbers ``x`` and ``y``;
    * ``hf_extrapolate_exp(e_x, e_y, n, alpha)`` -- the HF energy as
      E(n) = E_CBS + B (n + 1) exp(-alpha sqrt(n)), from cardinal numbers ``n``
      and ``n + 1`` (alpha 1.63 in G4; 9 in Karton and Martin's form);
    * ``schwenke(e_x, e_y, factor)`` -- Schwenke's two-point extrapolation;
    * ``scale(e, factor)``;
    * ``hf_part(c)`` and ``corr_part(c)`` -- the HF and correlation parts of a
      correlated component.

    ``E0`` must be assigned: it is the composite energy. Every value must scale
    like an energy, so a combination gives the same result in any unit; this is
    checked when the recipe is read. Empirical constants are given in hartree.

``description``, ``reference``, ``doi``, ``uncertainty`` and ``validation``
describe a recipe without changing its energy. A recipe's fingerprint covers every
other field, so a result records exactly which recipe made it.

The thermal terms follow from the geometry level: H(0 K) = E0 + s ZPE, and
H(T) = E0 + s ZPE + [H(T) - E - ZPE], with s the ZPE scale factor.

Where recipes come from
=======================

The shipped recipes come with the plug-in. Recipes of your own go in a
``recipes`` directory of SEAMM's data path (the job's data, or yours). A recipe of
your own may not reuse the name and version of a shipped one.

What is not there yet
=====================

Each component needs a complete ``level``; choosing the code from a method and a
basis set is not implemented yet. ``hf_part`` and ``corr_part`` need the HF part
of a component, which the codes do not return yet. A recipe whose energy comes from one
calculation runs anywhere a model chemistry is used; a recipe that combines
several runs through the steps that use ``composite_step.evaluate.evaluate_recipe``,
such as Thermochemistry. Gradients of composite methods are not available.

Indices and tables
==================

* :ref:`genindex`
* :ref:`search`
