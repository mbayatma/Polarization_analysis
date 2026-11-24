"""DNN architectures for polarization / EFT-vs-SM studies."""

import tensorflow as tf
from tensorflow.keras import layers, models


def build_two_output_dnn(input_dim, learning_rate=1e-3):
    """Build the SM-only two-output polarization classifier.

    A single Dense(2, sigmoid) output: P(W boson longitudinal), P(Z boson
    longitudinal). No theory (SM-vs-EFT) supervision at all -- this is
    intended to be trained only on SM events, then applied unchanged
    (inference only) to EFT events to see whether the SM-trained scores
    drift for a given EFT operator.
    """
    model = models.Sequential([
        layers.Input(shape=(input_dim,)),

        layers.Dense(256, activation="relu"),
        layers.BatchNormalization(),
        layers.Dropout(0.3),

        layers.Dense(128, activation="relu"),
        layers.BatchNormalization(),
        layers.Dropout(0.3),

        layers.Dense(64, activation="relu"),
        layers.BatchNormalization(),
        layers.Dropout(0.3),

        layers.Dense(32, activation="relu"),

        # output 0 = P(W_L), output 1 = P(Z_L)
        layers.Dense(2, activation="sigmoid"),
    ])

    model.compile(
        optimizer=tf.keras.optimizers.Adam(learning_rate=learning_rate),
        loss="binary_crossentropy",
        metrics=[
            tf.keras.metrics.BinaryAccuracy(name="binary_accuracy"),
            tf.keras.metrics.AUC(name="auc"),
        ],
    )

    return model


def build_binary_dnn(input_dim, learning_rate=1e-3):
    """Build the single-output W_L Z_L vs. rest classifier.

    This is the simplest of the three models: one sigmoid asking "is this
    event fully longitudinal". It collapses the four polarization states into
    one bit, so unlike the two-output model it cannot say *which* boson is
    longitudinal -- but it is the natural baseline, and the quantity it
    targets (the fully-longitudinal fraction) is the one most directly tied
    to electroweak symmetry breaking.

    Kept architecturally identical to build_two_output_dnn apart from the
    output layer, so the comparison between them reflects the choice of
    target rather than a difference in capacity.
    """
    model = models.Sequential([
        layers.Input(shape=(input_dim,)),

        layers.Dense(256, activation="relu"),
        layers.BatchNormalization(),
        layers.Dropout(0.3),

        layers.Dense(128, activation="relu"),
        layers.BatchNormalization(),
        layers.Dropout(0.3),

        layers.Dense(64, activation="relu"),
        layers.BatchNormalization(),
        layers.Dropout(0.3),

        layers.Dense(32, activation="relu"),

        layers.Dense(1, activation="sigmoid"),
    ])

    model.compile(
        optimizer=tf.keras.optimizers.Adam(learning_rate=learning_rate),
        loss="binary_crossentropy",
        metrics=[
            tf.keras.metrics.BinaryAccuracy(name="binary_accuracy"),
            tf.keras.metrics.AUC(name="auc"),
        ],
    )

    return model


def build_multitask_dnn(input_dim, learning_rate=1e-3):
    """Build the shared-trunk, two-head polarization/theory classifier.

    Head 1 (pol_output, 2 sigmoids): P(W boson longitudinal), P(Z boson longitudinal).
    Head 2 (theory_output, 1 sigmoid): P(event is EFT, vs. SM).

    The trunk is trained jointly on both objectives, so it has to learn
    kinematic structure that is useful both for identifying polarization
    and for separating EFT from SM kinematics.
    """
    inputs = layers.Input(shape=(input_dim,), name="features")

    x = layers.Dense(256, activation="relu")(inputs)
    x = layers.BatchNormalization()(x)
    x = layers.Dropout(0.3)(x)

    x = layers.Dense(128, activation="relu")(x)
    x = layers.BatchNormalization()(x)
    x = layers.Dropout(0.3)(x)

    x = layers.Dense(64, activation="relu")(x)
    x = layers.BatchNormalization()(x)
    x = layers.Dropout(0.3)(x)

    shared = layers.Dense(32, activation="relu", name="shared_representation")(x)

    pol_head = layers.Dense(32, activation="relu")(shared)
    pol_output = layers.Dense(2, activation="sigmoid", name="pol_output")(pol_head)

    theory_head = layers.Dense(32, activation="relu")(shared)
    theory_output = layers.Dense(1, activation="sigmoid", name="theory_output")(theory_head)

    model = models.Model(
        inputs=inputs,
        outputs={"pol_output": pol_output, "theory_output": theory_output},
    )

    model.compile(
        optimizer=tf.keras.optimizers.Adam(learning_rate=learning_rate),
        loss={
            "pol_output": "binary_crossentropy",
            "theory_output": "binary_crossentropy",
        },
        loss_weights={"pol_output": 1.0, "theory_output": 1.0},
        metrics={
            "pol_output": [
                tf.keras.metrics.BinaryAccuracy(name="binary_accuracy"),
                tf.keras.metrics.AUC(name="auc"),
            ],
            "theory_output": [
                tf.keras.metrics.BinaryAccuracy(name="binary_accuracy"),
                tf.keras.metrics.AUC(name="auc"),
            ],
        },
    )

    return model
