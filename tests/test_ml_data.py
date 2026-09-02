import math

import numpy as np
import pandas as pd
import pytest

from ml.data import (
    add_labels,
    clean_features,
    get_feature_columns,
    get_wz_binary_labels,
    stratify_key,
)


@pytest.mark.parametrize(
    "polarization,expected",
    [
        (r"$W_L Z_L$", (1, 1)),
        (r"$W_L Z_T$", (1, 0)),
        (r"$W_T Z_L$", (0, 1)),
        (r"$W_T Z_T$", (0, 0)),
    ],
)
def test_get_wz_binary_labels_pure_states(polarization, expected):
    assert get_wz_binary_labels(polarization) == expected


def test_get_wz_binary_labels_madspin_excluded():
    w_is_l, z_is_l = get_wz_binary_labels(r"$W_T Z_T$ (Madspin)")
    assert math.isnan(w_is_l)
    assert math.isnan(z_is_l)


def test_add_labels_drops_unlabelable_rows_and_sets_is_eft():
    df = pd.DataFrame(
        {
            "polarization": [r"$W_L Z_L$", r"$W_T Z_T$ (Madspin)", r"$W_T Z_L$"],
            "theory": ["SM", "SM", "EFT"],
        }
    )
    out = add_labels(df)

    assert len(out) == 2
    assert out["W_is_L"].tolist() == [1, 0]
    assert out["Z_is_L"].tolist() == [1, 1]
    assert out["is_EFT"].tolist() == [0, 1]


def test_get_feature_columns_excludes_metadata():
    df = pd.DataFrame(
        {
            "M_mumu": [1.0],
            "MET": [2.0],
            "polarization": ["x"],
            "theory": ["SM"],
            "operator": ["SM"],
            "W_is_L": [1],
            "Z_is_L": [1],
            "is_EFT": [0],
        }
    )
    assert get_feature_columns(df) == ["M_mumu", "MET"]


def test_clean_features_fills_nan_with_train_medians():
    X_train = pd.DataFrame({"a": [1.0, 2.0, 3.0], "b": [10.0, np.nan, 30.0]})
    X_train_clean, medians = clean_features(X_train)

    assert X_train_clean["b"].isna().sum() == 0
    assert X_train_clean["b"].iloc[1] == pytest.approx(20.0)  # median of [10, 30]

    # Applying the *training* medians to a different (e.g. test) split must
    # not leak information from that split's own values: inf/nan there fall
    # back to the median computed on the training split, not this split's own.
    X_other = pd.DataFrame({"a": [np.inf, -np.inf], "b": [np.nan, 5.0]})
    X_other_clean, _ = clean_features(X_other, feature_medians=medians)
    assert X_other_clean["a"].iloc[0] == pytest.approx(2.0)  # train median of a
    assert X_other_clean["a"].iloc[1] == pytest.approx(2.0)
    assert X_other_clean["b"].iloc[0] == pytest.approx(20.0)  # train median of b


def test_clean_features_clips_extreme_finite_values():
    X_train = pd.DataFrame({"a": [1.0, 2.0, 1e9]})
    X_train_clean, _ = clean_features(X_train)
    assert X_train_clean["a"].iloc[2] == pytest.approx(1e6)


def test_stratify_key_combines_theory_and_polarization_labels():
    df = pd.DataFrame(
        {
            "theory": ["SM", "EFT"],
            "W_is_L": [1, 0],
            "Z_is_L": [0, 1],
        }
    )
    assert stratify_key(df).tolist() == ["SM_1_0", "EFT_0_1"]
