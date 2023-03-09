import sys
import os
from tqdm import tqdm
import numpy as np
import pandas as pd
import copy
import shutil


# Read paths
path_to_parameters = sys.argv[1]
path_to_log = sys.argv[2]
path_to_cavity = sys.argv[3]
path_to_rbcavity = sys.argv[4]


# create rbcavity function
def create_rbcavity(path_to_parameters, path_to_log, path_to_cavity, path_to_rbcavity):
    """
    Create cavity before docking with rDock.
    
    Args:
        
        path_to_parameters (str): Path to structure parameters file
        path_to_log (str): Path to log
        path_to_cavity (str): Path to cavity (GRD extension)
        path_to_rbcavity (str): Path to rbcavity (default utils/rDock/rbcavity)
    
    
    """

    # Run rbcavity
    command = path_to_rbcavity + ' -was -d -r ' + path_to_parameters + ' > ' + path_to_log
    out = os.popen(command).read()
    sys.stderr.write(out)




# os.environ["RBT_ROOT"] = "utils/rDock"
# os.environ["RBT_HOME"] = "utils/rDock"
# os.environ["LD_LIBRARY_PATH"] = "$LD_LIBRARY_PATH:utils/rDock/lib"
# os.environ["PATH"] = "$PATH:utils/rDock/bin"

# sys.stderr.write("\n\n" + str(os.getenv('RBT_ROOT')) + "\n\n")
# sys.stderr.write("\n\n" + str(os.getenv('RBT_HOME')) + "\n\n")
# sys.stderr.write("\n\n" + str(os.getenv('LD_LIBRARY_PATH')) + "\n\n")
# sys.stderr.write("\n\n" + str(os.getenv('PATH')) + "\n\n")
# sys.stderr.flush()



create_rbcavity(path_to_parameters, path_to_log, path_to_cavity, path_to_rbcavity)