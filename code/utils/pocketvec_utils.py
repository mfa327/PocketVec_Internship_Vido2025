from subprocess import Popen, PIPE
import scipy.stats as ss
from tqdm import tqdm
import pandas as pd
import numpy as np
import pickle
import shutil
import copy
import pybel
import sys
import os



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
    moefunct = "code/utils/moefunctions.svl"
    
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
        
        
        
        
def environmental_variables(path_to_rDock):
    """
    Set environmental variables to run rDock.
    
    Args:
       
        path_to_rDock (str): Path to rDock installation directory. You may use the precompiled files we
        provide in code/utils/rDock_compiled or your local installed version of rDock.
    
    
    """
    
    if path_to_rDock[-1] == '/':
        path_to_rDock = path_to_rDock[:-1]
    
    os.environ["RBT_ROOT"] = path_to_rDock
    os.environ["LD_LIBRARY_PATH"] = path_to_rDock + "/lib"
    os.environ["PATH"] = os.environ["PATH"] + ":" + path_to_rDock + "/bin"
    
    
    

    
def create_cavity(outpath, path_to_log, path_to_cavity, path_to_rbcavity='rbcavity'):
    """
    Create cavity before docking with rDock. 3 files are generated:
        - log file
        - .grd file (docking cavity)
        - .as file
    
    Args:
        
        path_to_parameters (str): Path to structure parameters file
        path_to_log (str): Path to log
        path_to_cavity (str): Path to cavity (GRD extension)
        path_to_rbcavity (str): Path to rbcavity (if environment variables are set properly, 'rbcavity' should be enough)
    
    
    """
    
    # Path to parameters file
    path_to_parameters = os.path.join(outpath, 'st_parameters.prm')
    
    # Run rbcavity
    command = path_to_rbcavity + ' -was -d -r ' + path_to_parameters + " | tee " + path_to_log
    out = os.popen(command).read()
    sys.stderr.write(out)
    
    # Rename grid file
    os.rename(os.path.join(outpath, "st_parameters_cav1.grd"), os.path.join(outpath, "cavity.grd"))
    
    
    
    
    
    
def run_rDock(outpath, path_to_lib='../data/libs/TOP_128_rDock_LLM.sdf', nruns=25, seed=42):
    
    """
    Run rDock. Results are placed in the 'results.sd' file located in 'outpath'.
    
    Args:
        
        outpath (str): Path where all parameter files (2 x .prm) are located. All generated files will also be located in this path.
        path_to_lib (str): Path to the library of molecules to dock. By default, 128 LLMs are used.
        nruns (int/str): Number of docking runs for each molecule. By default, 25 runs are performed.
        seed (int/str): Seed to make docking calculations reproducible.
    
    
    """

    # Define variables
    results = os.path.join(outpath, 'results')  # SD file with results
    path_to_st_parameters = os.path.join(outpath, 'st_parameters.prm') # Parameters (st)
    path_to_rdock_parameters = os.path.join(outpath, 'dock.prm') # Parameters (dock)
    rdock_logs = os.path.join(outpath, 'rDock_log.log')
    nruns = str(nruns)
    seed = str(seed)

    command = 'rbdock -i ' + path_to_lib + ' -o ' + results + ' -r ' + path_to_st_parameters + ' -p ' + path_to_rdock_parameters + ' -n ' + str(nruns) + ' -s ' + seed + ' | tee ' + rdock_logs
    os.system(command)
    
    
    
    
def create_file_scores(path_to_scores, path_to_results):
    
    """
    Given a file with docking results (results.sd), create a file summarizing the minimum score for each docked molecule.
    
    Args:
        
        path_to_scores (str): Summary file with minimum docking scores.
        path_to_results (str): File with docking results (sd format)
    
    """
    
    
    with open(path_to_scores, "w") as out_score:
        
        # Read file
        file = pybel.readfile("sd", path_to_results)
        
        # Min affinity is set to an arbitrary large number
        aff = 1000000
        affinities = {}
        
        # Iterate over molecules and store the minimum score for each one
        for molecule in file:
            if molecule.title not in affinities:
                affinities[molecule.title] = aff
            value = float(molecule.data['SCORE.INTER'])
            if value < affinities[molecule.title]:
                affinities[molecule.title] = value

        # Print results
        for molecule in sorted(affinities):
            out_score.write(molecule + "\t" + str(affinities[molecule]) + "\n")
    

def rank_fp(fp_raw):
    
    """
    Given a fp with raw docking scores, convert it to ranks and deal with outlier molecules appropriately
    
    Args:
    
        fp_raw (list/np array): fp with raw docking scores
    
    """
    
    # Rank
    fp = ss.rankdata(fp_raw, method='max')
    # Increase ranking for positive scores
    for c in range(len(fp)):
        if fp_raw[c] > 0 and fp_raw[c] < 50:
            fp[c] = len(fp) + 1
        elif fp_raw[c] > 50 and fp_raw[c] < 100:
            fp[c] = len(fp) + 2
        elif fp_raw[c] > 100:
            fp[c] = len(fp) + 3
    return fp
    

def read_scores(path_to_scores):
    
    """
    Function to read tab-separated scores (molecule \t score)
    
    Args:
    
        path_to_scores (str): Path to file.
    
    """
    
    scores = {}
    with open(path_to_scores, "r") as f:
        for l in f:
            scores[l.split("\t")[0]] = float(l.split("\t")[1])
    return scores




def raw_fp(dict_scores, file_order):
    
    """
    Given a dict with the minimum docking score for each molecule and file with the order of molecules, create a fp of raw docking scores
    
    Args:
    
        dict_scores (dict): molecule-docking score
        file_order (str): path to a pickle file having a list with the ordered molecules
    
    """
    
    molecules = pickle.load(open(file_order, "rb"))
    return np.array([dict_scores[i] for i in molecules])
