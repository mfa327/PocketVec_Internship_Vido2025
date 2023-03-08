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
