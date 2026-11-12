# -*- coding: utf-8 -*-

"""Non-graphical part of the Composite step in a SEAMM flowchart"""

import logging
import importlib.resources
import pprint  # noqa: F401

import composite_step
import molsystem
import seamm
from seamm_util import ureg, Q_  # noqa: F401
import seamm_util.printing as printing
from seamm_util.printing import FormattedText as __

# In addition to the normal logger, two logger-like printing facilities are
# defined: "job" and "printer". "job" send output to the main job.out file for
# the job, and should be used very sparingly, typically to echo what this step
# will do in the initial summary of the job.
#
# "printer" sends output to the file "step.out" in this steps working
# directory, and is used for all normal output from this step.

logger = logging.getLogger(__name__)
job = printing.getPrinter()
printer = printing.getPrinter("Composite")

# Add this module's properties to the standard properties
path = importlib.resources.files("composite_step") / "data"
csv_file = path / "properties.csv"
if path.exists():
    molsystem.add_properties_from_file(csv_file)


class Composite(seamm.Node):
    """
    The non-graphical part of a Composite step in a flowchart.

    Attributes
    ----------
    parser : configargparse.ArgParser
        The parser object.

    options : tuple
        It contains a two item tuple containing the populated namespace and the
        list of remaining argument strings.

    subflowchart : seamm.Flowchart
        A SEAMM Flowchart object that represents a subflowchart, if needed.

    parameters : CompositeParameters
        The control parameters for Composite.

    See Also
    --------
    TkComposite,
    Composite, CompositeParameters
    """

    def __init__(
        self, flowchart=None, title="Composite", extension=None, logger=logger
    ):
        """A step for Composite in a SEAMM flowchart.

        You may wish to change the title above, which is the string displayed
        in the box representing the step in the flowchart.

        Parameters
        ----------
        flowchart: seamm.Flowchart
            The non-graphical flowchart that contains this step.

        title: str
            The name displayed in the flowchart.
        extension: None
            Not yet implemented
        logger : Logger = logger
            The logger to use and pass to parent classes

        Returns
        -------
        None
        """
        logger.debug(f"Creating Composite {self}")

        super().__init__(
            flowchart=flowchart,
            title="Composite",
            extension=extension,
            module=__name__,
            logger=logger,
        )  # yapf: disable

        self._metadata = composite_step.metadata
        self.parameters = composite_step.CompositeParameters()

    @property
    def version(self):
        """The semantic version of this module."""
        return composite_step.__version__

    @property
    def git_revision(self):
        """The git version of this module."""
        return composite_step.__git_revision__

    def description_text(self, P=None):
        """Create the text description of what this step will do.
        The dictionary of control values is passed in as P so that
        the code can test values, etc.

        Parameters
        ----------
        P: dict
            An optional dictionary of the current values of the control
            parameters.
        Returns
        -------
        str
            A description of the current step.
        """
        if not P:
            P = self.parameters.values_to_dict()

        text = (
            "Composite methods are used through a Model Chemistry step, as "
            "SEAMM:CMP@<recipe> (e.g. SEAMM:CMP@revDSD-r2SCAN-3c-HO), by the steps "
            "that take a model chemistry, such as Thermochemistry. This step "
            "itself does nothing yet."
        )

        return self.header + "\n" + __(text, **P, indent=4 * " ").__str__()

    def run(self):
        """Run a Composite step: print what it does, which is nothing yet.

        Returns
        -------
        seamm.Node
            The next node object in the flowchart.
        """
        next_node = super().run(printer)
        printer.important(__(self.description_text({}), indent=self.indent))
        return next_node

    def analyze(self, indent="", **kwargs):
        """Nothing to analyze yet."""
