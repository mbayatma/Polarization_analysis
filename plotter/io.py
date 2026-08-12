import os
import uproot


def load_SM_and_EFT_trees(source="spin"):
    base_dir = "/eos/user/m/mbayatma/polarization_files"

    sm_dir = os.path.join(base_dir, "SMRootFiles")
    eft_dir = os.path.join(base_dir, "BSMRootFiles")

    # Change this if your tree has another name
    tree_name = "mytree"

    pol_map = {
        "$W_L Z_L$": "WL_ZL",
        "$W_T Z_T$": "WT_ZT",
        "$W_T Z_L$": "WT_ZL",
        "$W_L Z_T$": "WL_ZT",
    }

    operators = [
        "cHWB",
        "cHDD",
        "cHl3",
        "cHW",
        "cWtil",
        "cHbox",
        "cHq1",
        "cHq3",
        "cHu",
        "cHd",
        "cll1",
    ]

    trees_SM = {}

    for pol_label, pol_tag in pol_map.items():
        sm_file = os.path.join(
            sm_dir,
            f"WZ_3L_{pol_tag}_SM.root"
        )

        if not os.path.exists(sm_file):
            print(f"[WARN] SM file not found for {pol_label}: {sm_file}")
            continue

        trees_SM[pol_label] = uproot.open(sm_file)[tree_name]

    trees_EFT = {}

    for op in operators:
        trees_EFT[op] = {}

        for pol_label, pol_tag in pol_map.items():
            eft_file = os.path.join(
                eft_dir,
                f"WZ_3L_{pol_tag}_{op}.root"
            )

            if not os.path.exists(eft_file):
                print(f"[WARN] EFT file not found for {op}, {pol_label}: {eft_file}")
                continue

            trees_EFT[op][pol_label] = uproot.open(eft_file)[tree_name]

        if len(trees_EFT[op]) == 0:
            del trees_EFT[op]

    print("Loaded SM samples:")
    for pol in trees_SM:
        print("  ", pol)

    print("Loaded EFT samples:")
    for op, trees in trees_EFT.items():
        print(f"  {op}: {list(trees.keys())}")

    return trees_SM, trees_EFT
    
    
    
    
    
    
    
    
    
    

