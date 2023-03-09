import sys
import os
from tqdm import tqdm
import numpy as np
import pandas as pd
import copy
import shutil
from subprocess import Popen, PIPE


def prepare_pdb(infile, outfile, logfile, moebatch = "/aloy/home/acomajuncosa/programs/MOE/moe2020/bin/moebatch"):
    """
    Prepare a PDB file for rDock docking. Requires MOE license & bin files.
    
    Args:
    
        infile (str): Path to input structure file (PDB format)
        outfile (str): Path to output structure file (MOL2 format)
        logfile (str): Path to LOG file. Check the logs for further details on protein preparation.
        moebatch (str): Path to your moebatch file. 
        
        
        

    """
    
    # Paths to moe functions
    moefunct = "utils/moefunctions.svl"
    
    # Prepare structures 
    command = moebatch + " -load " + moefunct + " -exec \"proteinprep['" + infile + "','" + outfile + "','" + logfile + "']\""
    
    process=Popen([command], stdout=PIPE, stderr=PIPE, shell=True)
    stdout, stderr = process.communicate()
    
    if "error" not in str(stderr).lower():
        sys.stderr.write("Preparation completed!")
        
    else:
        sys.stderr.write(str(stdout) + "\n\n")
        sys.stderr.write(str(stderr) + "\n\n")
        
        
        
        
def create_parameter_file(outfile, path_to_st, path_to_ctr, radius=str(12.0)):
    """
    Create structure parameter file.
    
    Args:
    
        outfile (str): Path to output structure parameter file (rDock-PRM format)
        path_to_st (str): Path to receptor (MOL2 format).
        path_to_ctr (str): Path to pocket centroid (SD format). 
        radius (float/str): Pocket radius (default 12.0).
        
    
    """
    
    text = """RBT_PARAMETER_FILE_V1.00
TITLE Change this if you want

RECEPTOR_FILE """ + path_to_st + """
RECEPTOR_FLEX 3.0

#################################################################
## CAVITY DEFINITION: REFERENCE LIGAND METHOD
#################################################################
SECTION MAPPER
    SITE_MAPPER RbtLigandSiteMapper
    REF_MOL """ + path_to_ctr + """
    RADIUS """ + radius + """
    SMALL_SPHERE 1.0
    MIN_VOLUME 100
    MAX_CAVITIES 1
END_SECTION


SECTION CAVITY
    SCORING_FUNCTION RbtCavityGridSF
    WEIGHT 1.0
END_SECTION

#################################################
# WARNING!!! WE ARE DOING RIGID LIGAND DOCKING!!
#################################################
SECTION LIGAND
    TRANS_MODE FREE
    ROT_MODE FREE
    DIHEDRAL_MODE """ + "FIXED" + """
END_SECTION"""
    
    with open(outfile, "w") as f:
        f.write(text)
