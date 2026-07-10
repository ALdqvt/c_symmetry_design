# MotifConfig dataclass + sweep generator

from dataclasses import dataclass
from itertools import product
import numpy as np
from scipy.spatial.transform import Rotation


@dataclass(frozen=True)
class MotifConfig:
    r: float  # Center of Mass distance from origin
    theta: float  # Longitude angle of u (COM -> Lig vector)
    phi: float  # Latitude angle of u
    zeta: float  # Rotation around u


def motif_extent_worst_case(chain_atoms, n_samples=72):
    """Max y-z extent of the motif over all rotations about x (u),
    i.e. the safe radius to use regardless of zeta."""
    coords = np.array([a.get_coord() for a in chain_atoms])
    max_extent = 0.0
    for angle_deg in np.linspace(0, 360, n_samples, endpoint=False):
        rot = Rotation.from_euler('x', angle_deg, degrees=True)
        rotated = rot.apply(coords)
        yz_extent = np.max(np.linalg.norm(rotated[:, 1:], axis=1))
        max_extent = max(max_extent, yz_extent)
    return max_extent


def generate_configs(r_min, r_max, r_step, grid_step_deg):
    r_values = np.arange(r_min, r_max + r_step, r_step)
    phis = np.arange(30, 150 + grid_step_deg, grid_step_deg)

    theta_max = 60
    theta_min = -theta_max
    for r in r_values:
        for phi in phis:
            n_theta_equator = (theta_max - theta_min) // grid_step_deg
            n_theta = max(1, round(n_theta_equator * np.sin(np.deg2rad(phi))))
            thetas = np.linspace(theta_min, theta_max, n_theta, endpoint=True)

            for theta in thetas:
                zetas = np.arange(0, 360, grid_step_deg)
                for zeta in zetas:
                    yield MotifConfig(r=float(r), theta=float(theta), phi=float(phi), zeta=float(zeta))
