# plotter/io.py

import os
import uproot
from .config import base_path, files_SM, files_EFT

tree_name = "PolarizationTree"

def load_tree(file_path):
    return uproot.open(file_path)[tree_name]

def load_all_trees(files_dict):
    """Returns a dict of {polarization: tree}"""
    trees = {}
    for pol, fname in files_dict.items():
        full_path = os.path.join(base_path, fname)
        trees[pol] = load_tree(full_path)
    return trees

def load_SM_and_EFT_trees():
    return load_all_trees(files_SM), load_all_trees(files_EFT)

