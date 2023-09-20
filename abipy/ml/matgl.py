"""
"""
from __future__ import annotations

#import os
#import shutil
#import warnings
import numpy as np
import pytorch_lightning as pl

from dgl.data.utils import split_dataset
from mp_api.client import MPRester
from pytorch_lightning.loggers import CSVLogger

import matgl
from matgl.ext.pymatgen import Structure2Graph, get_element_list
from matgl.graph.data import M3GNetDataset, MGLDataLoader, collate_fn_efs
from matgl.models import M3GNet
from matgl.utils.training import PotentialLightningModule
from monty.string import list_strings
from abipy.ml.tools import get_structures_labels_from_files


class MatglSystem:
    """
    See also: https://github.com/materialsvirtuallab/matgl/blob/main/examples/Training%20a%20M3GNet%20Potential%20with%20PyTorch%20Lightning.ipynb
    """

    def __init__(self, filepaths, workdir, verbose):
        self.filepaths = list_strings(filepaths)
        self.workdir = workdir
        self.verbose = verbose

    def train(self):
        structures, labels = get_structures_labels_from_files(self.filepaths)
        magmoms = labels.pop("magmoms", None)
        #print(labels)

        element_types = get_element_list(structures)
        converter = Structure2Graph(element_types=element_types, cutoff=5.0)

        dataset = M3GNetDataset(
            threebody_cutoff=4.0,
            structures=structures,
            converter=converter,
            labels=labels,
        )

        train_data, val_data, test_data = split_dataset(
            dataset,
            frac_list=[0.8, 0.1, 0.1],
            shuffle=True,
            random_state=42,
        )

        train_loader, val_loader, test_loader = MGLDataLoader(
            train_data=train_data,
            val_data=val_data,
            test_data=test_data,
            collate_fn=collate_fn_efs,
            batch_size=2,
            num_workers=1,
        )

        model = M3GNet(
            element_types=element_types,
            is_intensive=False,
        )

        lit_module = PotentialLightningModule(model=model)

        # If you wish to disable GPU or MPS (M1 mac) training, use the accelerator="cpu" kwarg.
        logger = CSVLogger("logs", name="M3GNet_training")
        # Inference mode = False is required for calculating forces, stress in test mode and prediction mode
        trainer = pl.Trainer(max_epochs=1, accelerator="cpu", logger=logger, inference_mode=False)
        trainer.fit(model=lit_module, train_dataloaders=train_loader, val_dataloaders=val_loader)

        # save trained model
        #model_export_path = "./trained_model/"
        #model.save(model_export_path)
        # load trained model
        #model = matgl.load_model(path = model_export_path)
        #
        # This code just performs cleanup for this notebook.
        #for fn in ("dgl_graph.bin", "dgl_line_graph.bin", "state_attr.pt", "labels.json"):
        #    try:
        #        os.remove(fn)
        #    except FileNotFoundError:
        #        pass

        #shutil.rmtree("logs")
        #shutil.rmtree("trained_model")
        #shutil.rmtree("finetuned_model")
