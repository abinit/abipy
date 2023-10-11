"""
"""
from __future__ import annotations

from abipy.ml.tools import get_structures_labels_from_files
from monty.string import list_strings
from chgnet.model.model import CHGNet
from chgnet.data.dataset import StructureData, get_train_val_test_loader
from chgnet.trainer import Trainer


class ChgnetSystem:
    """
    See also https://github.com/CederGroupHub/chgnet/blob/main/examples/fine_tuning.ipynb
    """
    def __init__(self, filepaths, workdir, verbose):
        """
        Args:
            filepaths: List of files with ab-initio results.
            workdir: Working directory.
            verbose: Verbosity level.
        """
        self.filepaths = list_strings(filepaths)
        self.workdir = workdir
        self.verbose = verbose

    def train(self):
        structures, labels = get_structures_labels_from_files(self.filepaths)
        magmoms = labels.pop("magmoms", None)

        dataset = StructureData(
            structures=structures,
            energies=labels["energies"],
            forces=labels["forces"],
            stresses=labels["stresses"],  # can be None
            magmoms=magmoms,              # can be None
        )

        train_loader, val_loader, test_loader = get_train_val_test_loader(
            dataset, batch_size=8, train_ratio=0.8, val_ratio=0.1,
        )

        model = CHGNet()

        # Define Trainer
        trainer = Trainer(
            model=model,
            targets="efs" if magmoms is None else "efsm",
            optimizer="Adam",
            scheduler="CosLR",
            criterion="MSE",
            #epochs=5,
            epochs=25,
            learning_rate=1e-2,
            use_device="cpu",
            print_freq=6,
        )

        trainer.train(train_loader, val_loader, test_loader,
                      save_dir=None)
