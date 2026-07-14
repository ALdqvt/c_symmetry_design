import numpy as np
import pytest
from motif_config import MotifConfig, generate_configs


def test_yields_motif_config_instances():
    configs = list(generate_configs(r_min=10, r_max=10, r_step=1, grid_step_deg=90))
    assert len(configs) > 0
    assert all(isinstance(c, MotifConfig) for c in configs)


def test_fields_are_floats():
    configs = list(generate_configs(r_min=10, r_max=10, r_step=1, grid_step_deg=90))
    c = configs[0]
    assert isinstance(c.r, float)
    assert isinstance(c.theta, float)
    assert isinstance(c.phi, float)
    assert isinstance(c.zeta, float)


def test_r_values_within_bounds():
    configs = list(generate_configs(r_min=5, r_max=10, r_step=5, grid_step_deg=90))
    r_values = {c.r for c in configs}
    assert r_values == {5.0, 10.0}


def test_phi_values_within_bounds():
    configs = list(generate_configs(r_min=10, r_max=10, r_step=1, grid_step_deg=90))
    for c in configs:
        assert 0 <= c.phi <= 180


def test_theta_values_within_bounds():
    configs = list(generate_configs(r_min=10, r_max=10, r_step=1, grid_step_deg=90))
    for c in configs:
        assert 0 <= c.theta < 360


def test_zeta_values_within_bounds():
    configs = list(generate_configs(r_min=10, r_max=10, r_step=1, grid_step_deg=90))
    for c in configs:
        assert 0 <= c.zeta < 360


def test_poles_have_fewer_theta_steps_than_equator():
    # at phi=0 (pole) there should be exactly 1 theta value;
    # at phi=90 (equator) there should be the full n_theta_equator count
    configs = list(generate_configs(r_min=10, r_max=10, r_step=1, grid_step_deg=90))
    theta_counts_by_phi = {}
    for c in configs:
        theta_counts_by_phi.setdefault(c.phi, set()).add(c.theta)

    n_theta_equator = 360 // 90  # = 4
    assert len(theta_counts_by_phi[0.0]) == 1        # pole: only 1 theta
    assert len(theta_counts_by_phi[90.0]) == n_theta_equator  # equator: full count


def test_generator_is_lazy():
    # generate_configs should not eagerly compute everything -
    # calling it should return a generator, not a list
    import types
    gen = generate_configs(r_min=10, r_max=10, r_step=1, grid_step_deg=10)
    assert isinstance(gen, types.GeneratorType)


def test_configs_are_frozen():
    configs = list(generate_configs(r_min=10, r_max=10, r_step=1, grid_step_deg=90))
    with pytest.raises(Exception):
        configs[0].r = 999.0