# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

from dataclasses import dataclass

import torch
from torch import Tensor

from chai_lab.utils.typing import Float, typecheck


@dataclass(frozen=True)
class InferenceNoiseSchedule:
    s_max: float = 160.0
    s_min: float = 4e-4
    p: float = 7.0
    sigma_data: float = 16.0

    @typecheck
    def get_schedule(
        self,
        device,
        num_timesteps: int = 200,
    ) -> Float[Tensor, "{num_timesteps}"]:
        times = torch.linspace(0, 1, 2 * num_timesteps + 1, device=device)[1::2]
        return self.get_noise_for_times(times)

    @typecheck
    def get_noise_for_times(
        self, times: Float[Tensor, "n_samples"]
    ) -> Float[Tensor, "n_samples"]:
        if times.min() < 0 or times.max() > 1:
            raise ValueError("times must be in [0, 1]")

        sigmas = self.sigma_data * _power_interpolation(
            times, val_0=self.s_max, val_1=self.s_min, p=self.p
        )
        return sigmas


@typecheck
def _power_interpolation(
    t: Float[Tensor, "n_samples"], val_0: float, val_1: float, p: float
) -> Float[Tensor, "n_samples"]:
    # val0 at t=0, and val1 at t=1
    assert t.min() >= 0 and t.max() <= 1, f"0 <= t <= 1, but {t=}"
    return (t * val_1 ** (1 / p) + (1 - t) * val_0 ** (1 / p)) ** p
