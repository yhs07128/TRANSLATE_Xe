import getopt
import sys
import numpy as np
import scipy.interpolate as interpolate
import math
from matplotlib import pyplot as plt

save_files = False
plot_data = False
plot_probabilities = False

m_e = 9.1e-31

def liquid_sigma_e(eV, gas_sigma_e):
    if eV < 2.85:
        return 5.5e-16
    return gas_sigma_e(eV)

def liquid_sigma_p(eV, S, gas_sigma_p):
    if eV < 3.95:
        return S(eV) * 5.5e-16
    return gas_sigma_p(eV)

def E_to_v(eV):
    return math.sqrt((2 * eV * 1.60218e-19) / (m_e))

if __name__ == '__main__':
    opts, args = getopt.getopt(sys.argv[1:], "sdp")

    for o, a in opts:
        if o == '-s':
            save_files = True
        if o == '-d':
            plot_data = True
        if o == '-p':
            plot_probabilities = True

    sigma_e_eVs, sigma_e = np.loadtxt('energy_xsec.txt', unpack=True)
    sigma_e *= 1e4
    sigma_e_func = interpolate.interp1d(sigma_e_eVs, sigma_e, kind='slinear')

    sigma_p_eVs, sigma_p = np.loadtxt('momentum_xsec.txt', unpack=True)
    sigma_p *= 1e4
    sigma_p_func = interpolate.interp1d(sigma_p_eVs, sigma_p, kind='slinear')

    sigma_i_eVs = np.array([0, 15, 17, 18.5, 20, 21, 22.5, 25, 27.5, 30, 32.5, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 90, 100, 110, 120, 140, 160, 180, 200, 225, 250, 275, 300, 350, 400, 500, 600, 700, 800, 900, 1000])
    sigma_i = np.array([0, 0, 0.159, 0.419, 0.604, 0.769, 1.00, 1.25, 1.58, 1.75, 2.07, 2.21, 2.41, 2.49, 2.53, 2.53, 2.51, 2.51, 2.52, 2.51, 2.54, 2.55, 2.51, 2.48, 2.42, 2.33, 2.25, 2.17, 2.09, 2.01, 1.91, 1.80, 1.73, 1.58, 1.47, 1.27, 1.13, 1.01, 0.914, 0.850, 0.783])
    sigma_i *= 1e-16
    sigma_i_func = interpolate.interp1d(sigma_i_eVs, sigma_i, kind='slinear')

    S_eVs = np.array([0, 1, 2, 3, 4, 5, 6, 7])
    S = np.array([0.048, 0.055, 0.125, 0.4, 1.115, 1.29, 1.16, 1.01])
    S_func = interpolate.interp1d(S_eVs, S, kind='quadratic', fill_value=(0.048, 1), bounds_error=False)

    x = np.logspace(-3, np.log10(290), 3000)

    if plot_data:
        plt.plot(x, sigma_e_func(x), label=r'GAr $\sigma_E$')
        plt.plot(x, sigma_p_func(x), label=r'GAr $\sigma_p$')
        plt.plot(x, sigma_i_func(x), label=r'GAr $\sigma_I$')
        plt.plot(x, np.full(len(x), 5.5e-16), label=r'LAr $\sigma_E$')
        plt.plot(x, S_func(x) * 5.5e-16, label=r'LAr $\sigma_p$')

        plt.ylabel(r'$\sigma$ (cm$^2$)')

    if plot_probabilities:
        y_vals1 = []
        for i in x:
            y_vals1.append(sigma_i_func(i) * E_to_v(i) * 1e2)
        plt.plot(x, y_vals1, label=r'$v\sigma_I(v)$')

        y_vals2 = []
        for i in x:
            y_vals2.append(liquid_sigma_e(i, sigma_e_func) * E_to_v(i) * 1e2)
        plt.plot(x, y_vals2, label=r'$v\sigma_E(v)$')

        y_vals3 = []
        for i in x:
            y_vals3.append(liquid_sigma_p(i, S_func, sigma_p_func) * E_to_v(i) * 1e2)
        plt.plot(x, y_vals3, label=r'$v\sigma_p(v)$')

        y_tots = []
        for i in range(len(x)):
            y_tots.append((y_vals1[i] + y_vals2[i] + y_vals3[i]))
        plt.plot(x, y_tots, label=r'$v\sigma_{total}(v)$')

        plt.ylabel(r'$v\sigma(v)$ (vcm$^2$)')

    if save_files:
        np.savetxt('sigma_i.txt', y_vals1, newline=',')
        np.savetxt('sigma_e.txt', y_vals2, newline=',')
        np.savetxt('sigma_p.txt', y_vals3, newline=',')
        np.savetxt('k_max.txt', y_tots, newline=',')

    plt.yscale('log')
    plt.xscale('log')

    plt.xlabel('eV')

    plt.legend()
    plt.grid()
    plt.show()
