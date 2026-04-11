#!/usr/bin/env python
"""

"""
import json
from math import cos, sin

import numpy as np

from ..utils import logger


def two(acoor, paras):
    i = np.array(acoor[0])
    j = np.array(acoor[1])
    a = paras[0]
    v = (1 - a) * i + a * j
    # logger.info(v)
    return v


def two_fd(acoor, paras):
    i = np.array(acoor[0])
    j = np.array(acoor[1])
    fixed_dist = paras[0]
    a = fixed_dist
    r_ij_length = np.sqrt(np.sum((i - j) ** 2, axis=0))
    unit_diff = (i - j) / r_ij_length
    v = i + a * unit_diff
    # logger.info(v)
    return v


def three(acoor, paras):
    i = np.array(acoor[0])
    j = np.array(acoor[1])
    k = np.array(acoor[2])
    a = paras[0]
    b = paras[1]
    v = (1 - a - b) * i + a * j + b * k
    logger.info(v)
    return v


def three_fd(acoor, paras):
    a = paras[0]
    b = paras[1]
    i = np.array(acoor[0])
    j = np.array(acoor[1])
    k = np.array(acoor[2])
    r_ij = i - j
    r_jk = i - k
    diff = (1 - a) * r_ij + a * r_jk
    length_of_diff = np.sqrt(np.sum(diff**2, axis=0))
    unit_diff = diff / length_of_diff
    v = i + b * unit_diff
    logger.info(v)
    return v


def three_fad(acoor, paras):
    fixed_dist = paras[0]
    fixed_angle = paras[1]
    cos_theti = cos(np.radians(fixed_angle))
    sin_theti = sin(np.radians(fixed_angle))
    i = np.array(acoor[0])
    j = np.array(acoor[1])
    k = np.array(acoor[2])
    r_ij = i - j
    r_jk = j - k

    quotient_of_dots = np.dot(r_ij, r_jk) / np.dot(r_ij, r_ij)
    r_ppd = r_jk - quotient_of_dots * r_ij

    ij_length = np.sqrt(np.sum(r_ij**2, axis=0))
    ppd_length = np.sqrt(np.sum(r_ppd**2, axis=0))

    cos_dist_term = fixed_dist * cos_theti * r_ij / ij_length
    sin_dist_term = fixed_dist * sin_theti * r_ppd / ppd_length

    v = i + cos_dist_term + sin_dist_term
    # logger.info(v)
    return v


def three_out(acoor, paras):
    i = np.array(acoor[0])
    j = np.array(acoor[1])
    k = np.array(acoor[2])
    a = paras[0]
    b = paras[1]
    c = paras[2]
    r_ij = i - j
    r_ik = i - k
    v = i + a * r_ij + b * r_ik + c * np.cross(r_ij, r_ik)
    logger.info(v)
    return v


def four_fd(acoor, paras):
    i = np.array(acoor[0])
    j = np.array(acoor[1])
    k = np.array(acoor[2])
    l = np.array(acoor[3])  # noqa
    a = paras[0]
    b = paras[1]
    c = paras[2]
    r_ij = i - j
    r_ik = i - k
    r_il = i - l
    r_ja = a * r_ik - r_ij
    r_jb = b * r_il - r_ij
    r_m = np.cross(r_ja, r_jb)
    r_m_length = np.sqrt(np.sum(r_m**2, axis=0))
    r_m_unit = r_m / r_m_length
    v = i + c * r_m_unit
    logger.info(v)
    return v


def double_three_out(acoor, paras):
    v1 = three_out(acoor, paras)
    paras[-1] = paras[-1] * -1
    v2 = three_out(acoor, paras)
    return [v1, v2]


__func_id = {
    "style1": two,
    "style2": two_fd,
    "style3": three,
    "style4": three_fd,
    "style5": three_fad,
    "style6": three_out,
    "style7": four_fd,
    "style8": double_three_out,
}


def vs_generate_coor(style, acoor, paras):
    return __func_id[style](acoor, paras)


def vs_generate(mole, vsfile):
    #vs_setting = json.loads(open(vsfile).read())
    vs_setting = {"cl_1":{"paras":[0.5],"style":"style2","binc":0.0403}}
    for aa in mole.Atoms:
        if aa.atom_type_name in vs_setting:
            aa.vs_setting = vs_setting[aa.atom_type_name]
        if hasattr(aa, "atom_type_name_m2"):
            if aa.atom_type_name_m2 in vs_setting:
                aa.vs_setting_m2 = vs_setting[aa.atom_type_name_m2]
    mole.create_vs()
    #mole.set_gromacs_atom_info("LIG")

    return mole
