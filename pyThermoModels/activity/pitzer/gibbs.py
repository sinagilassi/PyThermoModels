"""Canonical molality-basis Pitzer excess Gibbs energy and derivatives.

The kernel uses coefficients in the convention documented in docs/pitzer-v1.md.
Its inputs are resolved before differentiation; ionic strength remains variable.
"""
from itertools import combinations
from math import isfinite, sqrt
from typing import Any, Callable, cast

import autograd.numpy as anp
import numpy as np
from autograd import grad

from .components import normalize_molalities, normalize_pitzer_components
from .interactions import (make_binary_interaction_key, make_theta_key,
    make_psi_cc_a_key, make_psi_c_aa_key, require_parameter)
from .parameters import PitzerMultiComponentBinaryParameters, PitzerThetaParameter, PitzerPsiParameter


def gibbs_gradient(kernel: Callable[[Any], Any], molalities: np.ndarray) -> np.ndarray:
    """Differentiate a scalar kernel using autograd's curried public API.

    Autograd decorates ``grad`` at runtime; its undecorated source signature
    is not the signature exposed by the decorator.
    """
    make_gradient = cast(Callable[[Callable[[Any], Any]], Callable[[Any], Any]], grad)
    return np.asarray(make_gradient(kernel)(molalities), dtype=float)


def _g(x):
    # Branching on a scalar is supported by autograd. The series avoids the
    # subtraction of nearly equal numbers in the closed expression.
    if abs(x) < 1e-3:
        return 1 - 2*x/3 + x*x/4 - x**3/15 + x**4/72 - x**5/420
    return 2 * (-anp.expm1(-x) - x*anp.exp(-x)) / x**2


def _j(x):
    # Pitzer (1975), eq. 47; see Wolery et al. (2010), Appendix A.
    return x / (4 + 4.581*x**(-0.7237)*anp.exp(-0.012*x**0.528))


def _lookup(mapping, key, family, strict, typ=None):
    value = require_parameter(mapping, key, family, strict)
    if value is not None and typ is not None and not isinstance(value, typ):
        raise TypeError(f"{family} parameter for {key!r} must be {typ.__name__}")
    return value


def prepare_terms(components, values, params, strict=True):
    """Resolve and validate interaction keys without discarding zero species."""
    cations = [i for i, c in enumerate(components) if c.is_cation()]
    anions = [i for i, c in enumerate(components) if c.is_anion()]
    neutrals = [i for i, c in enumerate(components) if c.is_neutral()]
    binary, pairs, triples, neutral_pairs, zetas = [], [], [], [], []
    for i in cations:
        for j in anions:
            p = _lookup(params.binary, make_binary_interaction_key(components[i], components[j]), "binary", strict, PitzerMultiComponentBinaryParameters)
            if p is not None:
                binary.append((i, j, p))
    for group, opposite, cation_pair in ((cations, anions, True), (anions, cations, False)):
        for i, j in combinations(group, 2):
            p = _lookup(params.theta, make_theta_key(components[i], components[j]), "theta", strict, PitzerThetaParameter)
            if p is not None:
                zi, zj = components[i].get_net_charge(), components[j].get_net_charge()
                if zi != zj and params.j_strategy is None:
                    raise ValueError("Unequal-charge mixing requires an explicit source-matched j_strategy")
                pairs.append((i, j, p.theta))
            for k in opposite:
                key = make_psi_cc_a_key(components[i], components[j], components[k]) if cation_pair else make_psi_c_aa_key(components[k], components[i], components[j])
                p = _lookup(params.psi_cc_a if cation_pair else params.psi_c_aa, key, "psi_cc_a" if cation_pair else "psi_c_aa", strict, PitzerPsiParameter)
                if p is not None:
                    triples.append((i, j, k, p.psi))
    formulas = [c.get_key("Formula") for c in components]
    for i in neutrals:
        for j in cations + anions:
            p = _lookup(params.lambda_ni, (formulas[i], formulas[j]), "lambda_ni", strict)
            if p is not None: neutral_pairs.append((i, j, float(p), 2))
        for j in neutrals:
            if j < i: continue
            key = tuple(sorted((formulas[i], formulas[j])))
            p = _lookup(params.lambda_nn, key, "lambda_nn", strict)
            if p is not None: neutral_pairs.append((i, j, float(p), 1 if i == j else 2))
        for j in cations:
            for k in anions:
                p = _lookup(params.zeta_g, (formulas[i], formulas[j], formulas[k]), "zeta_g", strict)
                if p is not None: zetas.append((i, j, k, float(p)))
    if any(not isfinite(p) for *_, p in zetas) or any(not isfinite(p) for _, _, p, _ in neutral_pairs):
        raise ValueError("Neutral interaction coefficients must be finite")
    return binary, pairs, triples, neutral_pairs, zetas


def excess_gibbs_molal(m, charges, params, terms):
    """Return Gex/(RT kg solvent), in mol/kg, for independent molalities."""
    I = anp.dot(m, charges * charges) / 2
    Z = anp.dot(m, anp.abs(charges))
    # Pure neutral systems have I identically zero and no electrostatic term.
    if I > 0:
        root = anp.sqrt(I)
        result = -4 * params.A_phi * I / params.b * anp.log1p(params.b * root)
    else:
        root = 0.0
        result = anp.sum(m) * 0.0
    binary, pairs, triples, neutral_pairs, zetas = terms
    for i, j, p in binary:
        B = p.beta0 + p.beta1 * _g(p.alpha1 * root)
        if p.beta2: B += p.beta2 * _g(p.alpha2 * root)
        C = p.c_phi / (2 * sqrt(abs(charges[i] * charges[j])))
        result += m[i]*m[j]*(2*B + Z*C)
    for i, j, theta in pairs:
        phi = theta
        if I > 0 and charges[i] != charges[j] and params.j_strategy == "pitzer_1975":
            factor = 6*params.A_phi*root
            ji, jj = charges[i], charges[j]
            phi += ji*jj/(4*I)*(_j(factor*ji*jj) - 0.5*_j(factor*ji*ji) - 0.5*_j(factor*jj*jj))
        result += 2*m[i]*m[j]*phi
    for i, j, k, p in triples: result += m[i]*m[j]*m[k]*p
    for i, j, p, factor in neutral_pairs: result += factor*m[i]*m[j]*p
    for i, j, k, p in zetas: result += m[i]*m[j]*m[k]*p
    return result


def evaluate_gibbs(components, molalities, params, strict=True, water_molar_mass=0.01801528):
    normalize_pitzer_components(components)
    values = normalize_molalities(components, molalities)
    m = np.array([values[c.get_key("Formula-State")] for c in components], dtype=float)
    charges = np.array([c.get_net_charge() for c in components], dtype=float)
    charge = float(np.dot(m, charges))
    if abs(charge) > 1e-12 + 1e-10*np.dot(m, np.abs(charges)):
        raise ValueError(f"Pitzer ionic state is not electroneutral: net molal charge = {charge:.6e}")
    if not isfinite(water_molar_mass) or water_molar_mass <= 0:
        raise ValueError("water_molar_mass must be finite and positive")
    terms = prepare_terms(components, values, params, strict)
    total = float(sum(m))
    if total == 0:
        G, derivatives, phi = 0.0, np.zeros(len(m)), 1.0
    else:
        kernel = lambda x: excess_gibbs_molal(x, charges, params, terms)
        G = float(kernel(m))
        derivatives = gibbs_gradient(kernel, m)
        phi = 1 + (float(np.dot(m, derivatives)) - G)/total
    ln_aw = -water_molar_mass*total*phi
    return G, derivatives, phi, ln_aw
