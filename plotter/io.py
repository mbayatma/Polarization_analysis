# plotter/io.py

import os
import uproot
from .config import ( 
             base_path, base_path_SM, base_path_BSM, 
             files_SM, files_EFT,
             files_SM_spin, files_BSM_spin,
            
               )

#tree_name = "PolarizationTree"

#adding the class in order to take the new madspin files to account

class WrappedTree:
      def __init__ (self, tree, num_entries):
          self.tree = tree
          self.num_entries = num_entries
          
      def __getitem__(self, key):
          return self.tree[key]
          
      def keys(self):
          return list(self.tree.keys())


def load_tree(file_path, tree_name):    
    return uproot.open(file_path)[tree_name]
    
    

def load_all_trees(files_dict, base_path, tree_name):
    """Returns a dict of {polarization: tree}"""
    trees = {}
    for pol, fname in files_dict.items():
        full_path = os.path.join(base_path, fname)
        trees[pol] = load_tree(full_path, tree_name)
    return trees

def load_SM_and_EFT_trees(source ="old"):
    """
    Load SM and EFT trees based on source type:
    - "old": Pythia-showered (PolarizationTree)
    - "spin": MadSpin/noMadSpin (mytree)
    """
    if source == "old":
       sm_dict = files_SM
       eft_dict = files_EFT
       base_sm = base_path
       base_eft = base_path
       tree_name = "PolarizationTree"
       
    elif source == "spin":
         sm_dict = files_SM_spin
         eft_dict = files_BSM_spin
         base_sm = base_path_SM
         base_eft = base_path_BSM
         tree_name= "mytree"
         
    else:
         raise ValueError(f"Unknown source type: {source}") 
    
    trees_SM = load_all_trees(sm_dict, base_sm, tree_name)
    trees_EFT = load_all_trees(eft_dict, base_eft, tree_name)
    return trees_SM, trees_EFT 
            
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

