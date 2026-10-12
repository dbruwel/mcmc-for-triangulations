"""New Script"""

import random

import numpy as np
from regina import Triangulation3


def get_pi(f3, gamma):

    return np.exp(-gamma * f3**2)


def get_aut(iso):
    t = Triangulation3.fromIsoSig(iso)  # type: ignore
    return len(t.findAllIsomorphisms(t))


def get_f_vector(iso):
    t = Triangulation3.fromIsoSig(iso)  # type: ignore
    f = t.fVector()
    return f


def perform_move(iso, a, i):
    if a == 1:
        target = Triangulation3.fromIsoSig(iso)  # type: ignore
        if target.pachner(target.triangle(i), True, False):
            target.pachner(target.triangle(i), False, True)
            tiso = target.isoSig_RidgeDegrees()

            return tiso
    if a == 2:
        target = Triangulation3.fromIsoSig(iso)  # type: ignore
        if target.pachner(target.edge(i), True, False):
            target.pachner(target.edge(i), False, True)
            tiso = target.isoSig_RidgeDegrees()

            return tiso


def get_sites(iso, f):
    up_sites = []
    down_sites = []
    target = Triangulation3.fromIsoSig(iso)  # type: ignore

    for t in range(f[2]):
        if target.pachner(target.triangle(t), True, False):
            up_sites.append(t)

    for e in range(f[1]):
        if target.pachner(target.edge(e), True, False):
            down_sites.append(e)

    return {"up": up_sites, "down": down_sites}


def choosemove(iso, gamma, sites=None, f=None, aut=None):
    # Prepare the triangulations details if not provided.
    if f is None:
        f = get_f_vector(iso)
    if sites is None:
        sites = get_sites(iso, f)
    if aut is None:
        aut = get_aut(iso)

    N_up = len(sites["up"])
    N_down = len(sites["down"])
    weight_up = N_up * get_pi(f[3] + 1, gamma)
    weight_down = N_down * get_pi(f[3] - 1, gamma)

    x = random.random()
    if x <= weight_up / (weight_up + weight_down):
        if N_up == 0:
            raise ValueError("No up moves available.")
        a = 1
        i = random.choice(sites["up"])
        iso_prop = perform_move(iso, a, i)
    else:
        if N_down == 0:
            raise ValueError("No down moves available.")
        a = 2
        i = random.choice(sites["down"])
        iso_prop = perform_move(iso, a, i)

    if iso_prop is None:
        raise ValueError("Move could not be performed.")

    f_prop = get_f_vector(iso_prop)
    aut_prop = get_aut(iso_prop)
    sites_prop = get_sites(iso_prop, f_prop)

    N_up_prop = len(sites_prop["up"])
    N_down_prop = len(sites_prop["down"])
    weight_up_prop = N_up_prop * get_pi(f_prop[3] + 1, gamma)
    weight_down_prop = N_down_prop * get_pi(f_prop[3] - 1, gamma)

    alpha = min(
        (weight_up + weight_down)
        / (weight_up_prop + weight_down_prop)
        * (aut_prop / aut),
        1,
    )

    u = random.random()
    if u <= alpha:
        return iso_prop, f_prop, sites_prop, aut_prop, alpha
    else:
        return iso, f, sites, aut, alpha
