import math

import scipy as sp
import numpy as np
import uncertainties as uc
from uncertainties import ufloat as uf
import uncertainties.unumpy as unp
from uncertainties.umath import *
from molmass import Formula
import statsmodels.api as sm
import matplotlib.pyplot as plt

from solution import Solution

# 1


ret_t = np.loadtxt('Data/retention_times.csv', delimiter=',', skiprows=1)

chromo_mean = np.mean(ret_t)
chromo_std_dev = np.std(ret_t, ddof=1)
chromo_std_err = chromo_std_dev / np.sqrt(ret_t.size)

ci_99 = sp.stats.norm.interval(0.99, loc=chromo_mean, scale=chromo_std_dev)
ci_95 = sp.stats.norm.interval(0.95, loc=chromo_mean, scale=chromo_std_dev)

# 2

cal_data = np.loadtxt('Data/calorimeter_results.csv', delimiter=',', skiprows=1)
cal_calibration = uf(2422, 3)


def energy(m, t_i, t_f): return m * (cal_calibration * (t_f - t_i))


cal_energy = np.array([energy(*data) for data in cal_data])
cal_mean = np.mean(np.array([uc.nominal_value(data) for data in cal_energy]))
cal_std = np.sqrt(sum([dev ** 2 for dev in unp.std_devs(cal_energy)]))
cal_std_err = cal_std / np.sqrt(len(cal_energy))

cal_ci_99 = sp.stats.norm.interval(0.99, loc=cal_mean, scale=cal_std_err)
cal_ci_95 = sp.stats.norm.interval(0.95, loc=cal_mean, scale=cal_std_err)

# 3

na_oh = Solution(solute_id='NaOH', solute_mass=4.6556, volume=500)
v_i = uf(43.74, 0.05)
v_f = uf(2.38, 0.05)


def unknown_mol(molarity, v_f, v_i): return molarity * (v_f - v_i)


benzoic_acid_mol = unknown_mol(na_oh.concentration.mol, v_f, v_i)

# 4

tube_vol = uf(623.2, 0.5)
p_i = uf(202.3, 0.5)
p_f = uf(123.6, 0.5)


def bulb_volume(tube_vol, p_i, p_f): return tube_vol * (p_i / p_f) - tube_vol


bulb_vol = bulb_volume(tube_vol, p_i, p_f)

# 5

arr = np.loadtxt('Data/arrhenius.csv', delimiter=',', skiprows=1)

arr_e = np.empty((0, 2))
for row in arr:
    ls = []
    for i in range(0, 3, 2):
        ls.append(uf(row[i], row[i + 1]))
    arr_e = np.vstack([arr_e, ls])

x_w_err = (1 / (arr_e[:, 0]))
y_w_err = arr_e[:,1]
for i in range(len(y_w_err)):
    y_w_err[i] = log(y_w_err[i]) # np.log won't accept non-float

x_mean = x_w_err.mean()
y_mean = y_w_err.mean()

m_nom = 0
for i in range(len(x_w_err)):
    m_nom = m_nom + (x_w_err[i] - x_w_err.mean()) * (y_w_err[i] - y_mean)

m_denom = 0
for i in range(len(x_w_err)):
    m_denom = m_denom + (x_w_err[i] - x_mean) ** 2

m = m_nom / m_denom
c = y_mean - m * x_mean

plt.plot(uc.nominal_value(x_w_err), uc.std_dev(y_w_err))
