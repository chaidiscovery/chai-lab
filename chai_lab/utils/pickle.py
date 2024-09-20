# Copyright (c) 2024 Chai Discovery, Inc.
# This source code is licensed under the Chai Discovery Community License
# Agreement (LICENSE.md) found in the root directory of this source tree.

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
