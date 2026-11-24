import os

import uproot

from plotter.config import (
    base_path_EFT,
    base_path_EFT_REGEN,
    base_path_SM,
    eft_operators,
    pol_to_short,
)

TREE_NAME = "mytree"

# The LHE-level production wrote cW as "NPcW"; the showered samples use "cW".
# Config and plot labels use the physics name, so the filename alias is kept here.
LHE_OPERATOR_ALIASES = {"cW": "NPcW"}


def resolve_eft_file(operator, pol_tag):
    """Locate the EFT ROOT file for one operator/polarization.

    The regenerated directory is searched first. 23 of the original 44 EFT
    files were byte-identical duplicates of each other (a production-script
    bug: the param_card wasn't rewritten between runs, so different Wilson
    coefficients produced the same events, `eventWeight` included). Those 23
    were regenerated; the rest of the originals were verified independent and
    are still used from the original directory. Preferring the regenerated
    copy here means callers get the fixed file wherever one exists, without
    having to know which operators were affected.
    """
    names = [f"WZ_3L_{pol_tag}_{operator}.root"]
    alias = LHE_OPERATOR_ALIASES.get(operator)
    if alias:
        names.append(f"WZ_3L_{pol_tag}_{alias}.root")

    for directory in (base_path_EFT_REGEN, base_path_EFT):
        for name in names:
            candidate = os.path.join(directory, name)
            if os.path.exists(candidate):
                return candidate, directory
    return None, None


def load_SM_and_EFT_trees(source="spin", operator=None):
    """Open the LHE-level SM and EFT trees.

    `operator` selects a single EFT operator; when omitted, every operator in
    `eft_operators` is loaded and the EFT result is nested as
    `trees_EFT[operator][polarization]` (which is what export_ml.py expects).
    Passing an operator returns the flat `trees_EFT[polarization]` form that
    main.py's per-operator loop expects.
    """
    operators = [operator] if operator is not None else list(eft_operators)

    trees_SM = {}
    for pol_label, pol_tag in pol_to_short.items():
        sm_file = os.path.join(base_path_SM, f"WZ_3L_{pol_tag}_SM.root")
        if not os.path.exists(sm_file):
            print(f"[WARN] SM file not found for {pol_label}: {sm_file}")
            continue
        trees_SM[pol_label] = uproot.open(sm_file)[TREE_NAME]

    trees_EFT = {}
    regen_used = []
    for op in operators:
        per_pol = {}
        for pol_label, pol_tag in pol_to_short.items():
            eft_file, directory = resolve_eft_file(op, pol_tag)
            if eft_file is None:
                print(f"[WARN] EFT file not found for {op}, {pol_label}")
                continue
            per_pol[pol_label] = uproot.open(eft_file)[TREE_NAME]
            if directory == base_path_EFT_REGEN:
                regen_used.append(f"{pol_tag}_{op}")
        if per_pol:
            trees_EFT[op] = per_pol

    print(f"Loaded SM samples: {sorted(trees_SM)}")
    for op, trees in trees_EFT.items():
        print(f"Loaded EFT {op}: {sorted(trees)}")
    if regen_used:
        print(f"Using regenerated EFT files ({len(regen_used)}): {sorted(regen_used)}")

    # A single requested operator collapses to the flat form.
    if operator is not None:
        return trees_SM, trees_EFT.get(operator, {})
    return trees_SM, trees_EFT
