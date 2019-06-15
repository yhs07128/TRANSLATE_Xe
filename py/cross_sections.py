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
    return np.sqrt((2 * eV * 1.60218e-19) / (m_e)) # m / s

if __name__ == '__main__':
    opts, args = getopt.getopt(sys.argv[1:], "sdp")

    for o, a in opts:
        if o == '-s':
            save_files = True
        if o == '-d':
            plot_data = True
        if o == '-p':
            plot_probabilities = True

    sigma_e_eVs, sigma_e = np.loadtxt('energy_xsec.txt', unpack=True)       # m^2
    sigma_e *= 1e4                                                          # changing to cm^2
    sigma_e_func = interpolate.interp1d(sigma_e_eVs, sigma_e, kind='slinear')

    sigma_p_eVs, sigma_p = np.loadtxt('momentum_xsec.txt', unpack=True)     # m^2
    sigma_p *= 1e4                                                          # changing to cm^2
    sigma_p_func = interpolate.interp1d(sigma_p_eVs, sigma_p, kind='slinear')

    sigma_i_eVs = np.array([0, 15.75, 17, 18.5, 20, 21, 22.5, 25, 27.5, 30, 32.5, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 90, 100, 110, 120, 140, 160, 180, 200, 225, 250, 275, 300, 350, 400, 500, 600, 700, 800, 900, 1000])                                                      # eV
    sigma_i = np.array([0, 0, 0.159, 0.419, 0.604, 0.769, 1.00, 1.25, 1.58, 1.75, 2.07, 2.21, 2.41, 2.49, 2.53, 2.53, 2.51, 2.51, 2.52, 2.51, 2.54, 2.55, 2.51, 2.48, 2.42, 2.33, 2.25, 2.17, 2.09, 2.01, 1.91, 1.80, 1.73, 1.58, 1.47, 1.27, 1.13, 1.01, 0.914, 0.850, 0.783])     # 1e-16 cm^2
    sigma_i *= 1e-16                                                                                                                                                                                                                                                                # changing to cm^2
    sigma_i_func = interpolate.interp1d(sigma_i_eVs, sigma_i, kind='slinear')

    S_eVs = np.array([0, 1, 2, 3, 4, 5, 6, 7])                              # eV
    S = np.array([0.048, 0.055, 0.125, 0.4, 1.115, 1.29, 1.16, 1.01])
    S_func = interpolate.interp1d(S_eVs, S, kind='quadratic', fill_value=(0.048, 1), bounds_error=False)

    x = np.logspace(-3, np.log10(290), 3000)

    if plot_data:
        plt.plot(x, sigma_e_func(x), label=r'GAr $\sigma_E$')
        plt.plot(x, sigma_p_func(x), label=r'GAr $\sigma_p$')
        plt.plot(x, sigma_i_func(x), label=r'GAr $\sigma_I$')
        plt.plot(x, np.full(len(x), 5.5e-16), label=r'LAr $\sigma_E$')      # Low velocity energy cross section
        plt.plot(x, S_func(x) * 5.5e-16, label=r'LAr $\sigma_p$')           # Low velocity momentum cross section

        plt.ylabel(r'$\sigma$ (cm$^2$)')

    sigma_i_final = sigma_i_func(x)
    sigma_e_final = np.array([])
    sigma_p_final = np.array([])
    sigma_sum = np.array([])

    for i in x:
        sigma_e_final = np.append(sigma_e_final, liquid_sigma_e(i, sigma_e_func)) # cm^3 / s

    for i in x:
        sigma_p_final = np.append(sigma_p_final, liquid_sigma_p(i, S_func, sigma_p_func)) # cm^3 / s

    sigma_sum = sigma_i_final + sigma_e_final + sigma_p_final # cm^3 / s

    if plot_probabilities:
        plt.plot(x, sigma_i_final * E_to_v(x) * 1e2, label=r'$v\sigma_I(v)$')
        plt.plot(x, sigma_e_final * E_to_v(x) * 1e2, label=r'$v\sigma_E(v)$')
        plt.plot(x, sigma_p_final * E_to_v(x) * 1e2, label=r'$v\sigma_p(v)$')
        plt.plot(x, sigma_sum * E_to_v(x) * 1e2 * 3, label=r'$v\sigma_{total}(v)$')
        plt.ylabel(r'$v\sigma(v)$ (vcm$^2$)')

    if save_files:
        eVs = np.linspace(0, 300, 3000)

        sigma_e_l = []
        sigma_p_l = []

        for i in eVs:
            sigma_e_l.append(liquid_sigma_e(i, sigma_e_func))

        for i in eVs:
            sigma_p_l.append(liquid_sigma_p(i, S_func, sigma_p_func))

        np.savetxt('sigma_e_gas.txt', sigma_e_func(eVs), newline=',')
        np.savetxt('sigma_p_gas.txt', sigma_p_func(eVs), newline=',')
        np.savetxt('K_max_gas.txt', (sigma_e_func(eVs) + sigma_p_func(eVs) + sigma_i_func(eVs)) * E_to_v(eVs) * 1e2, newline=',')
        np.savetxt('sigma_i.txt', sigma_i_func(eVs), newline=',')
        np.savetxt('sigma_e.txt', sigma_e_l, newline=',')
        np.savetxt('sigma_p.txt', sigma_p_l, newline=',')
        np.savetxt('k_max.txt', sigma_sum * E_to_v(eVs) * 1e2 * 3, newline=',')

    if plot_data or plot_probabilities:
        plt.yscale('log')
        plt.xscale('log')

        plt.xlabel('eV')

        plt.legend()
        plt.grid()
        plt.show()