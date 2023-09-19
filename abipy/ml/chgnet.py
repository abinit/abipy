"""
"""
from __future__ import annotations

from abipy.ml.tools import get_structures_labels_from_files
from chgnet.model.model import CHGNet
from chgnet.data.dataset import StructureData, get_train_val_test_loader
from chgnet.trainer import Trainer


class ChgnetSystem:
    """
    See also https://github.com/CederGroupHub/chgnet/blob/main/examples/fine_tuning.ipynb
    """

    def __init__(self, filepaths, verbose):
        self.filepaths = filepaths
        self.verbose = verbose

    def train(self):
        structures, labels = get_structures_labels_from_files(self.filepaths)

        magmoms = None
        #magmoms = labels["magmoms"]

        dataset = StructureData(
            structures=structures,
            energies=labels["energies"],
            forces=labels["forces"],
            stresses=labels["stresses"],  # can be None
            magmoms=magmoms,              # can be None
        )

        train_loader, val_loader, test_loader = get_train_val_test_loader(
            dataset, batch_size=8, train_ratio=0.9, val_ratio=0.05
        )

        model = CHGNet()

        # Define Trainer
        trainer = Trainer(
            model=model,
            targets="efs" if magmoms is None else "efsm",
            optimizer="Adam",
            scheduler="CosLR",
            criterion="MSE",
            epochs=5,
            learning_rate=1e-2,
            use_device="cpu",
            print_freq=6,
        )


        trainer.train(train_loader, val_loader, test_loader)
