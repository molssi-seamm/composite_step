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

Three recipes give accurate reference values, after W4 (Karton et al., J. Chem.
Phys. 125, 144108 (2006)), on r2SCAN-3c geometries with scaled harmonic thermal
terms. Each adds to a valence CCSD(T) energy the core-valence correlation (all of
the outer core correlated, minus a frozen valence core, with cc-pwCVTZ) and the
scalar relativistic correction (DKH2 with a point nucleus, minus the
non-relativistic calculation):

``CCSD(T)-F12+CV+DKH2-r2SCAN-3c-HO``
    The valence energy from CCSD(T)-F12D/RI with cc-pVTZ-F12 (ORCA's F12, not
    Molpro's F12b), and the relativistic correction with cc-pVTZ-DK. For H-Ar:
    cc-pVTZ-F12 exists for no heavier element.

``CCSD(T)-CBS(TQ)+CV+DKH2-r2SCAN-3c-HO``
    The HF energy with cc-pVQZ and the CCSD(T) correlation extrapolated from
    cc-pVTZ and cc-pVQZ, and the relativistic correction with cc-pVTZ-DK. For H-Ar
    and Ga-Br, with all-electron basis sets.

``CCSD(T)-CBS(TQ)+CV+DKH2(decontracted)-r2SCAN-3c-HO``
    As the previous one, but the relativistic correction uses decontracted
    cc-pVTZ in both halves, for species with Ca, which has no DK basis. It agrees
    with the cc-pVTZ-DK correction within 0.003 kJ/mol on two W4 reactions and
    within 0.06 kJ/mol on GeH4 -> GeH2 + H2. X2C with x2c-TZVPall was not used:
    on CaH2 its correction to the atomization energy differed from decontracted
    DKH2 by 0.82 kJ/mol.

The core-valence and relativistic corrections agree with W4's to 0.13 and 0.01
kJ/mol on the reactions CO + 3 H2 -> CH4 + H2O and CO2 + 4 H2 -> CH4 + 2 H2O.

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
      correlated component, and ``part(c, "name")`` -- any named part (see
      `Parts of an energy`_).

    ``E0`` must be assigned: it is the composite energy. Every value must scale
    like an energy, so a combination gives the same result in any unit; this is
    checked when the recipe is read. Empirical constants are given in hartree.

``description``, ``reference``, ``doi``, ``uncertainty`` and ``validation``
describe a recipe without changing its energy. A recipe's fingerprint covers every
other field, so a result records exactly which recipe made it.

The thermal terms follow from the geometry level: H(0 K) = E0 + s ZPE, and
H(T) = E0 + s ZPE + [H(T) - E - ZPE], with s the ZPE scale factor.

Parts of an energy
==================

A correlated calculation reports the parts of its energy with it, so a recipe can
treat them separately: extrapolate the correlation energy but not the HF energy,
or scale the triples. ORCA reports them for MP2, RI-MP2, CCSD(T), DLPNO-CCSD(T)
and their F12 forms:

=====================  ============================================================
Run                    Parts
=====================  ============================================================
MP2                    ``HF``, ``MP2 correlation``
RI-MP2                 ``HF``, ``RI-MP2 correlation``
CCSD(T)                ``HF``, ``CCSD correlation``, ``(T)``; ``MP2 correlation``
                       (closed-shell RHF only)
CCSD(T)-F12D/RI        ``HF``, ``CABS singles``, ``CCSD correlation``,
                       ``F12 correction``, ``(T)``; ``MP2 correlation`` (RHF only)
DLPNO-CCSD(T)          ``HF``, ``DLPNO-MP2 correlation``, ``DLPNO-CCSD
                       correlation``, ``DLPNO-(T0)``
DLPNO-CCSD(T)-F12D     ``HF``, ``CABS singles``, ``DLPNO-CCSD correlation``,
                       ``F12 correction``, ``DLPNO-(T0)``, ``DLPNO-MP2
                       correlation``
=====================  ============================================================

Some parts are computed from others when asked for: ``CCSD-F12 correlation``
(CCSD correlation + F12 correction), ``HF+CABS``, ``(T*)`` (the triples scaled by f,
the ratio of the CCSD correlation with and without F12) and their DLPNO forms.

* ``HF`` is the energy of the reference determinant the correlation was built on:
  the SCF energy, or for quasi-restricted orbitals (QRO) the QRO determinant's
  energy, not the UHF one. For F12 it does not include the CABS singles.
* ``(T)`` is never scaled. f is ORCA's own value for canonical F12; for DLPNO-F12,
  where ORCA prints none, it is computed with the same definition.
* The DLPNO correlation is everything ORCA's total includes. ORCA's printed
  "uncorrected" DLPNO-F12 correlation leaves some of it out, so the CCSD
  correlation is the total correlation less the F12 correction and the triples.
* The parts of a calculation always sum to its energy, within 2e-9 E_h. If they
  do not, no parts are reported.

``hf_part(c)`` is ``HF``, and ``corr_part(c)`` is the correlation alone: the energy
less HF and the CABS singles, which belong to neither. A recipe that splits an F12
component with ``hf_part`` and ``corr_part`` must use its ``part(c, "CABS
singles")``, or say ``omit: ["CABS singles"]`` for that component; otherwise it is
refused when it is read. HF parts of components with different references (RHF,
UHF, QRO) cannot be combined unless the recipe says ``allow_mixed_reference:
true``.

The parts belong to their calculation. They are stored with it in the fragment
store, and a recipe reuses them only from the same calculation (the same settings
fingerprint). A stored energy without parts, needed by a recipe that uses parts, is
computed again; parts are never derived.

Frozen cores and relativistic corrections
=========================================

A component can be given:

``frozen_core``
    The name of a set in ``frozen_core_definitions`` (the number of frozen
    ORBITALS for each element), ``none`` (all electrons correlated), or
    ``default`` (the program's choice; not allowed in shipped recipes). A
    structure with an element the set does not cover is refused. ORCA counts core
    electrons, and its defaults are not the usual conventions for Li, Be, Na, Mg,
    K, Ca or Ga-Kr, so recipes give the counts explicitly. After each calculation
    the number of frozen orbitals is checked against the count asked for.

``relativistic``, ``nucleus``
    A scalar-relativistic Hamiltonian, ``DKH2`` or ``X2C``, with the nucleus model
    stated explicitly, ``point`` or ``finite``. The shipped recipes use a point
    nucleus: our explicit choice; W4 does not state one.

``decontract``
    ``true`` decontracts the orbital basis.

The shipped recipes' sets are ``valence`` (the cc-pVnZ conventions: 1s frozen for
Li-Ne, 1s-2p for Na-Ar, 1s-3p for K and Ca, 1s-3d for Ga-Kr) and ``inner`` (the
cc-pwCVnZ conventions: all electrons for Li-Ne, 1s frozen for Na-Ar, 1s-2p for K
and Ca and for Ga-Kr), each with its source in ``frozen_core_sources``.

A difference of two components, ``a - b``, must change one thing only: the
method, the basis set, the frozen core, or the relativistic Hamiltonian (with its
basis), so that it is the correction it is meant to be. Anything else is refused
when the recipe is read.

Open-shell species
==================

``open_shell_reference`` in a recipe sets the reference its correlated components
use for open-shell species (radicals, triplets, open-shell atoms): ``QRO``,
quasi-restricted orbitals built from the UHF solution, or ``UHF``. A component can
override it with its own ``reference``; both halves of a difference must use the
same one. The geometry level keeps its program's default. ROHF is not offered:
ORCA's canonical coupled cluster refuses it.

The three tier-A recipes use QRO. It changes nothing for closed-shell species: their
calculations, fingerprints and stored results are the same with or without it. For
an open-shell species the calculation's settings fingerprint tells QRO and UHF
apart, and after each run the reference the program used is checked against the
one asked for.

Where recipes come from
=======================

The shipped recipes come with the plug-in. Recipes of your own go in a
``recipes`` directory of SEAMM's data path (the job's data, or yours). A recipe of
your own may not reuse the name and version of a shipped one.

What is not there yet
=====================

A recipe runs anywhere a model chemistry is used: the Energy step, the
Thermochemistry step's CBH approach, and the others. A recipe that combines several
calculations runs each component as one batch of its own program over all the
structures, the cheapest components first. A structure that fails one component is
not sent to the later ones, and its result names the component that failed. Each
component's memory and cores are its own program's; note that a program can use
more memory than it is told: ORCA's canonical CCSD(T) with a quadruple-zeta basis
used about 2.7 times its ``%maxcore``, so size such components by the memory they
actually use. Each
result carries the composite record: the recipe's fingerprint together with every
component's settings fingerprint.

Not there yet:

* Each component needs a complete ``level``; choosing the code from a method and a
  basis set is not implemented yet.
* Only ORCA reports the parts of its energies so far.
* Empirical terms that depend on the structure (higher-level corrections,
  spin-orbit tables) are not implemented; a recipe that needs them is refused.
* Composite methods give energies only: no gradients, so they cannot be used to
  optimize structures or by steps that need forces.

Indices and tables
==================

* :ref:`genindex`
* :ref:`search`
