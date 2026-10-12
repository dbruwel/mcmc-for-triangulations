"""New Script"""

import math
import random

from regina.engine import Triangulation3


def get_aut(tri: Triangulation3) -> int:
    """Gets the size of the automorphism group for a given triangulation.

    Args:
        tri (Triangulation3): The triangulation

    Returns:
        int: The size of the automorphism group
    """
    return len(tri.findAllIsomorphisms(tri))


def get_f_vector(tri: Triangulation3) -> list[int]:
    """Compute the fVector for a given triangulation

    Args:
        tri (Triangulation3): The triangulation

    Returns:
        list[int]: The fVector
    """
    f = tri.fVector()
    return f


def perform_move(tri: Triangulation3, a: int, i: int) -> Triangulation3:
    """Performs a given Pachner move on a specified site of a triangulation.

    Args:
        tri (Triangulation3): The initial triangulation
        a (int): The type of move. 1 is a (3-2) move, 2 is a (2-3) move
        i (int): The id of the face / edge to perform the move around

    Returns:
        Triangulation3: The modified triangulation
    """
    target = Triangulation3(tri)

    if a == 1:
        target.pachner(target.triangle(i), False, True)
    if a == 2:
        target.pachner(target.edge(i), False, True)

    return target


def get_sites(tri: Triangulation3, f: list[int], a: int) -> list[int]:
    """Enumerates all sites where a given move type can be performed.

    Args:
        tri (Triangulation3): Triangulation3
        f (list[int]): fVector of the triangulation
        a (int): The move type, see `perform_move` for details

    Returns:
        list[int]: A list of the site ids where the given move type an be performed
    """
    sites = []
    target = Triangulation3(tri)

    if a == 1:
        for t in range(f[2]):
            if target.pachner(target.triangle(t), True, False):
                sites.append(t)

    if a == 2:
        for e in range(f[1]):
            if target.pachner(target.edge(e), True, False):
                sites.append(e)

    return sites


### Choose move from proposal
def choosemove(
    iso: str,
    gamma: float,
    f: list[int] | None = None,
    aut: int | None = None,
) -> tuple[str, list[int], int]:
    """Performes a step in the MCMC sequence

    Args:
        iso (str): The isomorhpism signature of the triangulation
        gamma (float): The up vs down acceptance rate
        f (list[int] | None, optional): The fVector of the triangulation associated to iso. Defaults to None.
        aut (int | None, optional): The size of the automorhpism group of the triangulation associated to iso. Defaults to None.

    Raises:
        RuntimeError: If the function passed without exiting

    Returns:
        tuple[str, list[int], int]: The new isomorphism signature, the new fVector, and the new automorphism group
    """
    tri = Triangulation3.fromIsoSig(iso)

    if aut is None:
        aut = get_aut(tri)
    if f is None:
        f = get_f_vector(tri)

    x = random.random()
    if x < math.exp((-1) * gamma * f[3]):
        a = 1
    else:
        a = 2

    sites = get_sites(tri, f, a)
    num_sites = len(sites)
    i = random.random()
    u = random.random()

    # go up (2-3)
    if a == 1:
        if i > float(num_sites) / float(f[2]):
            return tri.isoSig_RidgeDegrees(), f, aut
        else:
            tri_prop = perform_move(tri, a, random.choice(sites))
            aut_prop = get_aut(tri_prop)
            accept_prob = min(1, aut_prop / aut)
            if u > accept_prob:
                return tri.isoSig_RidgeDegrees(), f, aut
            return (
                tri_prop.isoSig_RidgeDegrees(),
                [
                    f[0],
                    f[1] + 1,
                    f[2] + 2,
                    f[3] + 1,
                ],
                aut_prop,
            )

    # go down
    elif a == 2:
        if f[3] <= 2:
            return tri.isoSig_RidgeDegrees(), f, aut
        if i > float(num_sites) / float(f[2] - 2):
            return tri.isoSig_RidgeDegrees(), f, aut
        else:
            tri_prop = perform_move(tri, a, random.choice(sites))
            aut_prop = get_aut(tri_prop)
            accept_prob = min(1, aut_prop / aut)
            if u > accept_prob:
                return tri.isoSig_RidgeDegrees(), f, aut
            return (
                tri_prop.isoSig_RidgeDegrees(),
                [
                    f[0],
                    f[1] - 1,
                    f[2] - 2,
                    f[3] - 1,
                ],
                aut_prop,
            )

    raise RuntimeError("Should not reach this point in choosemove")
