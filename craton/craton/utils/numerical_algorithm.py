import math

import numpy as np
import scipy.optimize as opt
from scipy.stats import linregress

# 'fmin', 'fmin_powell', 'fmin_bfgs', 'fmin_ncg', 'fmin_cg', 'fminbound', 'fmin_l_bfgs_b',
# 'fmin_tnc', 'fmin_cobyla', 'fmin_slsqp',


def rmse_calculate(x0, x1):
    assert len(x0) == len(x1)
    diff = np.abs(np.array(x0) - np.array(x1))
    return np.sqrt(np.average(np.square(diff))), np.average(diff)

def rmse_r_calculate(x0, x1):
    assert len(x0) == len(x1)
    diff = np.abs((np.array(x0) - np.array(x1))/np.array(x0))
    return np.sqrt(np.average(np.square(diff))), np.average(diff)


def linear_fitting(X, Y):
    assert len(X) == len(Y)
    if len(X) <= 2:
        # Manual defined trivial case solution per old implementation.
        return (1.0, 1.0, 1.0)
    a, b, r = linregress(X, Y)[:3]
    return (a, b, r * r)


def block_average(
    arr,
    interval=None,
    blocks=None,
):
    """
    块平均的统计方法
    输入： arr       ->    list  用于处理的数据
          interval  ->    None or int   间隔的距离
          blocks    ->    None or int   块的数目
    """
    if interval is None:
        if blocks is None:
            interval = int(len(arr) / 10)
        else:
            interval = int(len(arr) / blocks)
    tmp = []
    n = int(len(arr) / interval)
    for ii in range(n):
        tmp.append(sum(arr[ii * interval : (ii + 1) * interval]) / interval)
    if n * interval < len(arr):
        tmp.append(sum(arr[n * interval :]) / (len(arr) - n * interval))
    return tmp


def histogram_analysis(
    arr,
    interval=None,
    center_type="zero",
    uplimit=None,
    downlimit=None,
    rs=None,
    ls=None,
    move_half=True,
):
    """
    Note(wenyu): Incomplete method.

    柱状图数据的分析方法
    输入： arr         ->    list  用于处理的数据
          interval    ->    None or float   间隔的距离
          center_type ->    str
                                "zero":以0为中心；
                                "mean": 以平均值为中心；
                                “median”：最大值和最小值的中间;
                                "max":从最大值开始向下划分区间;
                                "min":从最小值开始向上划分区间;
          uplimit     ->  None or float  数个上限
          downlimit   ->  None or float  数值下限
          rs          ->   None or int   右边半边的区间数目
          ls          ->   None or int   左边半边的区间数目
          move_half   ->  True or False  区别边界是否移动一半的间隔值
    """
    if interval is None:
        if ls is None and rs is None:
            interval = (max(arr) - min(arr)) / 10.0
        elif ls is None and rs is not None:
            interval = (max(arr) - min(arr)) / rs
        elif ls is not None and rs is None:
            interval = (max(arr) - min(arr)) / ls
        else:
            interval = (max(arr) - min(arr)) / (ls + rs)

    if center_type == "zero":
        center = 0.0
    elif center_type == "mean":
        center = int((sum(arr) / len(arr)) / interval)
    elif center_type == "median":
        center = (max(arr) - min(arr)) / 2.0
    elif center_type == "max":
        center = max(arr)
    elif center_type == "min":
        center = min(arr)

    if uplimit is None:
        uplimit = max(arr)
    if downlimit is None:
        downlimit = min(arr)
    if rs is None:
        if move_half:
            rs = math.ceil((uplimit - (center + (interval) * 0.5)) / interval)
        else:
            rs = math.ceil((uplimit - center) / interval)
    if ls is None:
        if move_half:
            ls = math.ceil(((center - (interval) * 0.5) - downlimit) / interval)
        else:
            ls = math.ceil((center - downlimit) / interval)
