# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

import antipickle
import torch


class TorchAntipickleAdapter(antipickle.AbstractAdapter):
    typestring = "torch"

    def __init__(self):
        self.cpu_device = torch.device("cpu")

    def check_type(self, obj):
        return type(obj) is torch.Tensor  # ignore inherited classes

    def to_dict(self, obj):
        assert obj.device == self.cpu_device, "serializing only cpu tensors"
        return {"data": antipickle.wrap(obj.numpy())}  # use numpy serialization

    def from_dict(self, d):
        return torch.from_numpy(d["data"])
