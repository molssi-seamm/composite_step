***************
Getting Started
***************

Installation
============
The Composite step is probably already installed in your SEAMM
environment, but if not or if you wish to check, follow the directions for the `SEAMM
Installer`_. The graphical installer is the easiest to use. In the SEAMM conda
environment, simply type:: 

  seamm-installer

or use the shortcut if you installed one. Switch to the second tab, `Components`, and
check for `composite-step`. If it is not installed, or
can be updated, check the box next to it and click `Install selected` or `Update
selected` as appropriate.

The non-graphical installer is also straightforward::

  seamm-installer install --update composite-step

will ensure both that it is installed and up-to-date.

.. _SEAMM Installer: https://molssi-seamm.github.io/installation/index.html

An example
==========
To use a composite method, give it to a **Model Chemistry** step as
``SEAMM:CMP@<recipe>``. For example, the enthalpy of formation of a molecule from
the Thermochemistry step's CBH approach with r2SCAN-3c geometries and
revDSD-PBEP86-D4(2021)/def2-TZVPP energies::

    Model Chemistry: SEAMM:CMP@revDSD-r2SCAN-3c-HO

(the Thermochemistry step's ``default`` preset does the same). The recipes that are
installed are listed in the Model Chemistry step's choices.

That should be enough to get started. For more detail about the functionality in this plug-in, see the :ref:`User Guide <user-guide>`.
