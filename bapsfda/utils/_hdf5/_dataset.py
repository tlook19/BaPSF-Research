import h5py
import numpy as np


def create_experiment_set(file_name, experiment_name, **kwargs):
    with h5py.File(f"{file_name}.hdf5", "w") as f:
        f.create_group("{experiment_name}")
        f["experiment"].attrs["name"] = experiment_name
        for key, value in kwargs.items():
            f["experiment"].attrs[key] = value


def create_dataset(file_name, experiment_name, dataset_name, data, **kwargs):
    with h5py.File(f"{file_name}.hdf5", "a") as f:
        f[f"{experiment_name}"].create_dataset(f"{dataset_name}", data=data)
        for key, value in kwargs.items():
            f[f"{experiment_name}/{dataset_name}"].attrs[key] = value
