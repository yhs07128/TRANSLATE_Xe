import sys
import getopt
from math import sqrt
import matplotlib.pyplot as plt
import csv
import numpy as np
import glob
import re

T = 87
n = 2.11e22

plot_mobility = False

def mobility(E):
    a0 = 551.6
    a1 = 7953.7
    a2 = 4440.43
    a3 = 4.29
    a4 = 43.63
    a5 = 0.2053
    T0 = 89

    num = a0 + a1*E + a2*E**(3/2) + a3*E**(5/2)
    den = 1 + (a1/a0)*E + a4*E**2 + a5*E**3

    return (num / den) * (T / T0)**(-3/2)

class Graph:
    def __init__(self, file_list):
        self.n = 0
        self.tail_percent = 0.05
        self.name = re.search(r'(?:- )(\d*)(?:V)', file_list[0]).group(1)
        
        self.t = []
        self.x = []
        self.y = []
        self.z = []
        self.ke = []
        self.drift = []
        self.ionized = []

        for file in file_list:
            if file.find('1D') is not -1:
                _t, _x, _ke, _drift, _ionized = np.loadtxt(file, delimiter=',', unpack=True)
                self.t.append(_t)
                self.x.append(_x)
                self.ke.append(_ke)
                self.drift.append(_drift)
                self.ionized.append(_ionized)
                self.n += 1
            elif file.find('3D') is not -1:
                _t, _x, _y, _z, _ke, _drift, _ionized = np.loadtxt(file, delimiter=',', unpack=True)
                self.t.append(_t)
                self.x.append(_x)
                self.y.append(_y)
                self.z.append(_z)
                self.ke.append(_ke)
                self.drift.append(_drift)
                self.ionized.append(_ionized)
                self.n += 1

    def drift_mean(self):
        mean = 0
        for i in range(self.n):
            tail_steps = int(len(self.t[i]) * self.tail_percent)
            mean += sum(self.drift[i][-tail_steps:]) / tail_steps
        mean /= self.n
        return mean

    def drift_std_dev(self):
        mean = self.drift_mean()
        sigma = 0
        for i in range(self.n):
            tail_steps = int(len(self.t[i]) * self.tail_percent)
            sigma += ((sum(self.drift[i][-tail_steps:]) / tail_steps) - mean)**2
        sigma /= (self.n - 1)
        return sqrt(sigma)
    
    def plot_drift_velocity(self):
        plt.figure()
        for i in range(self.n):
            plt.plot(self.t[i], self.drift[i])

        mean = self.drift_mean()
        plt.plot(self.t[0], np.full(len(self.t[0]), mean), 'b--', label=r'$\bar{v}$')
        plt.plot(self.t[0], np.full(len(self.t[0]), mobility(int(self.name) * 1e-3) * int(self.name) * 1e-2), 'r--', label='Data from [1]')

        plt.xlabel('t (ns)')
        plt.ylabel('v (m/s)')
        plt.title(self.name + 'V $-$ Drift Velocity ($\mu = $' + str(int(mean)) + r', $\sigma = $' + str(int(self.drift_std_dev())) + ')')
        plt.figtext(0.08, 0.01, '[1] Y. Li, et al. "Measurement of Longitudinal Electron Diffusion in Liquid Argon", NIMA 816, 160 (2016).', fontsize='small')

        plt.grid()
        plt.legend()
    
    def plot_energy(self):
        plt.figure()
        plt.hist(self.ke, rwidth=0.9)

        plt.xlabel('kinetic energy (eV)')
        plt.ylabel('instances')
        plt.title(self.name + 'V $-$ Energy')

        plt.grid()
        plt.yscale('log')

if __name__ == '__main__':
    file_list = sorted(glob.glob('sim_data/*.txt'))
    opts, args = getopt.getopt(sys.argv[1:], "mg:v")
    file_name = None

    for o, a in opts:
        if o == '-m':
            plot_mobility = True
        if o == '-g':
            file_name = a

    if file_name is not None:
        new_list = []
        for file in file_list:
            if (file.find(file_name) is not -1):
                new_list.append(file)

        graph = Graph(new_list)

        graph.plot_drift_velocity()
        graph.plot_energy()

        plt.show()
    else:
        files_1d = {}
        files_3d = {}

        for file in file_list:
            key = int(re.search(r'(?:- )(\d*)(?:V)', file).group(1))
            if (file.find('1D') is not -1):
                if (key in files_1d):
                    files_1d[key].append(file)
                else:
                    files_1d[key] = [file]
            elif (file.find('3D') is not -1):
                if (key in files_3d):
                    files_3d[key].append(file)
                else:
                    files_3d[key] = [file]

        volts_1d = []
        means_1d = []
        std_devs_1d = []

        for key, value in sorted(files_1d.items()):
            graph = Graph(value)
            volts_1d.append(key)
            means_1d.append(graph.drift_mean())
            std_devs_1d.append(graph.drift_std_dev())

        volts_3d = []
        means_3d = []
        std_devs_3d = []

        for key, value in sorted(files_3d.items()):
            graph = Graph(value)
            volts_3d.append(key)
            means_3d.append(graph.drift_mean())
            std_devs_3d.append(graph.drift_std_dev())

        if plot_mobility is True:
            means_1d = [means_1d[i] * (n / float(volts_1d[i])) * 1e2 for i in range(len(means_1d))]
            std_devs_1d = [std_devs_1d[i] * (n / float(volts_1d[i])) * 1e2 for i in range(len(means_1d))]
            volts_1d = [float(i) / n for i in volts_1d]

            means_3d = [means_3d[i] * (n / float(volts_3d[i])) * 1e2 for i in range(len(volts_3d))]
            std_devs_3d = [std_devs_3d[i] * (n / float(volts_3d[i])) * 1e2 for i in range(len(volts_3d))]
            volts_3d = [float(i) / n for i in volts_3d]

            x = np.logspace(np.log10(volts_1d[0]), np.log10(volts_1d[-1]), 15)
            y = mobility(x * n * 1e-3) * n

            a = [1.4e-21, 2.8e-21, 6.2e-21, 1.2e-20, 2.5e-20, 5.1e-20, 1e-19, 2.2e-19, 4.6e-19, 1e-18]
            b = [9e24, 10e24, 10e24, 9e24, 6e24, 4.4e24, 3e24, 2e24, 1e24, 7e23]

            plt.plot(x, y, label='Data from [1]')
            plt.plot(a, b, label='Data from [2]')

            plt.xscale('log')
            plt.yscale('log')

            plt.title('e$^-$ mobility')
            plt.figtext(0.08, 0.03, '[1] Y. Li, et al. "Measurement of Longitudinal Electron Diffusion in Liquid Argon", NIMA 816, 160 (2016).', fontsize='small')
            plt.figtext(0.08, 0.01, '[2] S.S.-S. Huang, G.R. Freeman, Phys Rev. A 24 (1981) 714.', fontsize='small')

            plt.xlim(0.1e-20, 100e-20)
            plt.ylim(0.1e24, 100e24)

            plt.xlabel('F/n $($Vcm$^2)$')
            plt.ylabel('$\mu$n $($cmVs$)^{-1}$')
        else:
            x = np.logspace(1, np.log10(2e4), 15)
            y = mobility(x * 1e-3) * x * 1e-2

            plt.plot(x, y, label='Data from [1]')

            plt.xscale('log')
            plt.yscale('log')

            plt.title('Drift Velocity')
            plt.figtext(0.08, 0.01, '[1] Y. Li, et al. "Measurement of Longitudinal Electron Diffusion in Liquid Argon", NIMA 816, 160 (2016).', fontsize='small')
            
            plt.xlabel('V/cm')
            plt.ylabel('m/s')

        plt.errorbar(volts_1d, means_1d, yerr=std_devs_1d, capsize=3, label='1D')
        plt.errorbar(volts_3d, means_3d, yerr=std_devs_3d, capsize=3, label='3D')

        plt.grid()
        plt.legend()
        plt.show()