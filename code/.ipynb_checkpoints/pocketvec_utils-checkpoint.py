import sys
import os
from tqdm import tqdm
import numpy as np
import pandas as pd
import copy
import shutil


def traverse_mpaths_from_mpath(metapath):
    """Return all the possible metapath starting for each node in a given metapath"""
    if type(metapath) == str:
        mp = metapath.split('-')
    else:
        mp = metapath
    metapaths = []
    for i in range(0,len(mp)-1,2):
        q = mp[i:] + mp[1:i+1]
        metapaths.append(q)
    return metapaths

def extend_sources_through_mpaths(mpath,sources):
    medges  = mpath2medges(mpath)
    srcs_iv = sources[::-1]
    v = sources + sources[::-1]
    while len(v) < len(medges):
        v+= sources[::-1]
    return v[:len(medges)]
