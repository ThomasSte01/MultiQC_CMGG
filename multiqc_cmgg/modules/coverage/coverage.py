import logging
from multiqc import config
from multiqc.base_module import BaseMultiqcModule

# Initialise the main MultiQC logger
log = logging.getLogger("multiqc")

log.info("started module")
class MultiqcModule(BaseMultiqcModule):
    def __init__(self):

         # Initialise the parent module Class object
        super(MultiqcModule, self).__init__(
            name="Coverage",
            anchor="Coverage",
            href="?",
            info="Pannel selection ",
        )