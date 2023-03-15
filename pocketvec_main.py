import sys
sys.path.insert(1, './code/utils')
from pocketvec_utils import create_parameter_file
from pocketvec_utils import environmental_variables
from pocketvec_utils import create_cavity
from pocketvec_utils import run_rDock
from pocketvec_utils import create_file_scores
from pocketvec_utils import read_rDock_scores
from pocketvec_utils import raw_fp
from pocketvec_utils import rank_fp
import argparse
import shutil
import pickle
import os



def run_pocketvec(receptor, pocket_centroid, outpath, lib, n_runs, radius, seed):


    print("###########################################")
    print("#######    SETTING PARAMETERS    ##########")
    print("###########################################\n")


    # 1. Create st parameters file for rDock docking
    outfile = os.path.join(outpath, "st_parameters.prm")
    create_parameter_file(outfile, receptor, pocket_centroid, radius=radius)


    # 2. Set environmental variables
    path_to_rDock = "code/utils/rDock_compiled/"
    environmental_variables(path_to_rDock)


    print("###########################################")
    print("#####    CREATING POCKET CAVITY    ########")
    print("###########################################\n")


    # 3. Run rbcavity to define the pocket
    path_to_parameters = os.path.join(outpath, 'st_parameters.prm')
    path_to_log = os.path.join(outpath, 'cavity_log.log')
    path_to_cavity = os.path.join(outpath, 'cavity.grd')
    path_to_rbcavity = 'rbcavity'
    create_cavity(outpath, path_to_log, path_to_cavity, path_to_rbcavity)



    # 4. Copy docking parameters file
    _ = shutil.copyfile("code/utils/rDock_compiled/dock.prm", os.path.join(outpath, "dock.prm"))


    print("\n\n")
    print("###########################################")
    print("####   RUNNING DOCKING CALCULATIONS   #####")
    print("###########################################\n")


    # 5. Run rDock
    run_rDock(outpath, lib, nruns=n_runs, seed=seed)



    print("\n\n")
    print("###########################################")
    print("####  CREATING POCKETVEC DESCRIPTOR   #####")
    print("###########################################\n")


    # 6. Create file with docking scores
    path_to_scores = os.path.join(outpath, "scores.tsv")
    path_to_results = os.path.join(outpath, "results.sd")
    create_file_scores(path_to_scores, path_to_results)


    # 8. Create PocketVec descriptor

    path_to_sorted_lib = 'data/libs/order/' + lib.split("/")[-1].split(".sdf")[0] + ".pkl"

    # Read scores
    scores = read_rDock_scores(path_to_scores)

    # Create fingerprint with raw docking scores
    raw = raw_fp(scores, path_to_sorted_lib)

    # Rank scores fp
    rank = rank_fp(raw)

    # Dump PocketVec descriptor in pickle format
    pickle.dump(rank, open(os.path.join(outpath, "PocketVec_fp.pkl"), "wb"))



    print("\n\n")
    print("###########################################")
    print("#######    PIPELINE COMPLETED !!  #########")
    print("###########################################\n")






def get_parser():
    parser = argparse.ArgumentParser(description='Runs PocketVec')

    # Required
    parser.add_argument('-r','--receptor', required=True, type=str, nargs='?', help='Path to protein receptor (MOL2 format)')
    parser.add_argument('-pc','--pocket_centroid', required=True, type=str, nargs='?', help='Path to pocket centroid (3D coordinate, SD format)')
    parser.add_argument('-o','--outpath', required=True, type=str, nargs='?', help='Path to output')


    # Optional
    parser.add_argument('-lib','--lib', required=False, type=str, default='data/libs/TOP_128_rDock_LLM.sdf', nargs='?', help='Path to compound library to dock (default data/libs/TOP_128_rDock_LLM.sdf)')
    parser.add_argument('-nr','--n_runs', required=False, type=str, default='25', nargs='?', help='Number of rDock docking runs (default 25)')
    parser.add_argument('-radius','--radius', required=False, type=str, default='12', nargs='?', help='Cavity radius (default 12 A)')
    parser.add_argument('-s','--seed', required=False, type=str, default='42', nargs='?', help='Seed (default 42)')


    return parser



if __name__ == "__main__":
    args = vars(get_parser().parse_args())
    print("\n\n\nINPUT PARAMETERS:\n")
    print("\n".join([i + " -- " + str(args[i]) for i in args]))
    print("\n")
    run_pocketvec(**args)

