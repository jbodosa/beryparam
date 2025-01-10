# Imports
import os
from pathlib import Path
import subprocess
import shutil
import logging

from openmm.app import *

from parse_param import *

####################
####################
# System (class) #
####################
####################

logger = logging.getLogger(__name__)

class System:
    def __init__(self):
        """
            System class which has all the attributes.
        """
