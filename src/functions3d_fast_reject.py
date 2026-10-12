"""Metropolis-Hastings variant of `choosemove` from `functions3d.py`.

Unlike the original (and `functions3d_fast_accept_all.py`), the target weight is
set explicitly by `get_pi`, currently exp(-gamma * f_3^2). The direction of the
move is proposed in proportion to (number of valid sites) * (target weight after
the move), a site is then chosen uniformly, and the move is accepted or rejected
with the Metropolis-Hastings ratio, including the automorphism group correction.

Valid sites for the current state are returned by `choosemove` and can be passed
back in on the next step, so each step enumerates sites only for the proposed
triangulation.
"""

import random

import numpy as np
from regina import Triangulation3


def get_pi(f3, gamma):
    """Return the (unnormalised) target weight for a triangulation with `f3` tetrahedra.

    Parameters
    ----------
    f3 : int
        Number of tetrahedra.
    gamma : float
        Strength of the penalty on the number of tetrahedra.

    Returns
    -------
    float
        exp(-gamma * f3**2).
    """

    return np.exp(-gamma * f3**2)


def get_aut(iso):
    """Return the size of the automorphism group of a triangulation.

    Parameters
    ----------
    iso : str
        Isomorphism signature of the triangulation.

    Returns
    -------
    int
        The number of combinatorial automorphisms.
    """
    t = Triangulation3.fromIsoSig(iso)  # type: ignore
    return len(t.findAllIsomorphisms(t))


def get_f_vector(iso):
    """Return the f-vector of a triangulation.

    Parameters
    ----------
    iso : str
        Isomorphism signature of the triangulation.

    Returns
    -------
    list[int]
        The f-vector (f_0, f_1, f_2, f_3).
    """
    t = Triangulation3.fromIsoSig(iso)  # type: ignore
    f = t.fVector()
    return f


def perform_move(iso, a, i):
    """Perform a Pachner move at a given site.

    Parameters
    ----------
    iso : str
        Isomorphism signature of the initial triangulation.
    a : int
        The move type. 1 is a 2-3 move (up, on triangle `i`), 2 is a 3-2 move (down, on edge `i`).
    i : int
        The index of the triangle (a=1) or edge (a=2) to move on.

    Returns
    -------
    str | None
        Isomorphism signature (`isoSig_RidgeDegrees`) of the new triangulation, or None if the move is not valid at that site.
    """
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
    """List all sites at which an up or down move can be performed.

    Parameters
    ----------
    iso : str
        Isomorphism signature of the triangulation.
    f : list[int]
        The f-vector of the triangulation.

    Returns
    -------
    dict[str, list[int]]
        `"up"` maps to the indices of triangles where a 2-3 move is valid, `"down"` to the indices of edges where a 3-2 move is valid.
    """
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
    """Perform one Metropolis-Hastings step of the MCMC chain.

    With W(T) = N_up * pi(f_3 + 1) + N_down * pi(f_3 - 1), an up move is
    proposed with probability N_up * pi(f_3 + 1) / W(T) (else a down move), at a
    uniformly random valid site. The proposal T' is accepted with probability

        alpha = min(1, W(T) / W(T') * |Aut(T')| / |Aut(T)|).

    Passing `sites`, `f` and `aut` from the previous step avoids recomputing
    them.

    Parameters
    ----------
    iso : str
        Isomorphism signature (`isoSig_RidgeDegrees`) of the current triangulation.
    gamma : float
        Parameter passed to `get_pi`.
    sites : dict[str, list[int]] | None, optional
        Output of `get_sites` for the current triangulation. Computed if None.
    f : list[int] | None, optional
        f-vector of the current triangulation. Computed if None.
    aut : int | None, optional
        Automorphism group size of the current triangulation. Computed if None.

    Returns
    -------
    tuple[str, list[int], dict[str, list[int]], int, float]
        The isomorphism signature, f-vector, sites and automorphism group size of the next state, and the acceptance probability `alpha`.

    Raises
    ------
    ValueError
        If the chosen direction has no valid sites, or the move could not be performed.
    """
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
