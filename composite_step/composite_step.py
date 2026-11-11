# -*- coding: utf-8 -*-

import logging
from pathlib import Path

import composite_step

logger = logging.getLogger(__name__)


class CompositeStep(object):
    """Helper class needed for the stevedore integration.

    This must provide a `description()` method that returns a dict containing a
    description of this node, and `create_node()` and `create_tk_node()` methods
    for creating the graphical and non-graphical nodes.

    The dictionary for the description is the class variable just below these
    comments. The felds are as follows:

        my_description : {str, str}
            A human-readable description of this step. It can be
            several lines long, and needs to be clear to non-expert users.
            It contains the following keys: description, group, name.

        my_description["description"] : tuple
            A description of the Composite step. It must be
            clear to non-experts.

        my_description["group"] : str
            Which group in the menus to put this step. If the group does
            not exist it will be created. Common groups are "Building",
            "Control", "Custom", "Data", and "Simulations".

        my_description["name"] : str
            The name of this step, to be displayed in the menus.
    """

    my_description = {
        "description": "An interface for Composite",
        "group": "Simulations",
        "name": "Composite",
    }

    def __init__(self, flowchart=None, gui=None):
        """Initialize this helper class, which is used by
        the application via stevedore to get information about
        and create node objects for the flowchart
        """
        pass

    @classmethod
    def get_model_chemistry_options(
        cls, periodic_only=False, mdi_only=False, data_path=None
    ):
        """The composite methods SEAMM provides: ``SEAMM:CMP@<recipe>`` for the
        newest version of each recipe found.

        The recipes are the shipped set plus those in a ``recipes/`` directory of
        each ``data_path`` entry (a job's data, the user's). Composites are
        molecular and have no MDI engine or gradients yet, so a periodic-only or
        MDI-only request gets nothing.

        Whether the codes for a recipe's components are installed is not checked
        here: that needs the other steps' model chemistries, and asking for them
        from inside model-chemistry discovery would recurse. The step checks it
        before it runs anything (:func:`composite_step.recipes.check_capabilities`).
        """
        from .recipes import RecipeLibrary

        if periodic_only or mdi_only:
            return {}
        extra = []
        for directory in data_path or ():
            extra.append(Path(directory) / "recipes")
        library = RecipeLibrary(extra_directories=extra)
        for problem in library.problems:
            logger.warning("Composite recipes: %s", problem)
        options = {}
        for name in library.names():
            recipe = library.get(name)
            scale, source = recipe.zpe_scale
            options[recipe.model_chemistry] = {
                "model_chemistry": recipe.model_chemistry,
                "type": "CMP",
                "method": recipe.name,
                "basis": None,
                "description": recipe.description or recipe.name,
                "periodic": False,
                "mdi_capable": False,
                "analytic_hessian": False,
                "parallel_tasks": True,
                "prefers_batch": True,
                "recipe": recipe.name,
                "recipe_version": recipe.version,
                "recipe_source": recipe.source,
                "geometry_level": recipe.geometry.get("level"),
                "components": {k: v.get("level") for k, v in recipe.components.items()},
                "zpe_scale": scale,
                "zpe_scale_source": source,
            }
        return options

    @classmethod
    def can_run_task(cls, configuration, model_chemistry, **kwargs):
        """Whether one task per structure can compute this composite: only when
        E0 depends on one energy (the component's own task), for a molecule."""
        from .evaluate import delegate
        from .recipes import RecipeError

        if int(getattr(configuration, "periodicity", 0) or 0) != 0:
            return False
        try:
            _, _, mc, provider = delegate(model_chemistry)
        except RecipeError:
            return False
        check = getattr(provider, "can_run_task", None)
        return True if check is None else check(configuration, mc, **kwargs)

    @classmethod
    def get_task(cls, configuration, model_chemistry, **kwargs):
        """The task of a one-energy recipe: its component's task. A recipe that
        combines several calculations runs through
        :func:`composite_step.evaluate.evaluate_recipe` instead."""
        from .evaluate import delegated_get_task

        return delegated_get_task(configuration, model_chemistry, **kwargs)

    @classmethod
    def analyze_task(cls, result, model_chemistry, configuration, **kwargs):
        """The component's results with the energy replaced by the recipe's E0
        (kJ/mol), and the composite described under ``composite``."""
        from .evaluate import delegated_analyze_task

        return delegated_analyze_task(result, model_chemistry, configuration, **kwargs)

    @classmethod
    def expected_fingerprints(cls, configuration, model_chemistry, **kwargs):
        """The component's own expected fingerprints, for a one-energy recipe."""
        from .evaluate import delegate
        from .recipes import RecipeError

        try:
            _, _, mc, provider = delegate(model_chemistry)
        except RecipeError:
            return None
        predict = getattr(provider, "expected_fingerprints", None)
        return None if predict is None else predict(configuration, mc, **kwargs)

    def create_node(self, flowchart=None, **kwargs):
        """Create and return the new node object.

        Parameters
        ----------
        flowchart: seamm.Node
            A non-graphical SEAMM node

        **kwargs : keyword arguments
            Various keyword arguments such as title, namespace or
            extension representing the title displayed in the flowchart,
            the namespace for the plugins of a subflowchart and
            the extension, respectively.

        Returns
        -------
        Composite
        """

        return composite_step.Composite(flowchart=flowchart, **kwargs)

    def create_tk_node(self, canvas=None, **kwargs):
        """Create and return the graphical Tk node object.

        Parameters
        ----------
        canvas : tk.Canvas
            The Tk Canvas widget

        **kwargs : keyword arguments
            Various keyword arguments such as tk_flowchart, node, x, y, w, h
            representing a graphical flowchart object, a non-graphical node for
            a step, and dimensions of the graphical node.

        Returns
        -------
        TkComposite
        """

        return composite_step.TkComposite(canvas=canvas, **kwargs)

    def description(self):
        """Return a description of what this step does.

        Returns
        -------
        description : dict(str, str)
        """
        return CompositeStep.my_description
