import sys
import getopt
import scipy.interpolate as interpolate
from math import sqrt, log
import matplotlib.pyplot as plt
import csv
import numpy as np
import glob
import re

class Graph:
    def __init__(self, file_list):
        self.n = 0
        self.tail_percent = 0.05
        
        self.t = []
        self.x = []
        self.y = []
        self.z = []
        self.ke = []
        self.drift = []
        self.ionized = []

        for file in file_list:
            _t, _x, _y, _z, _ke, _drift, _ionized = np.loadtxt(file, delimiter=',', unpack=True)
            self.t.append(_t)
            self.x.append(_x)
            self.y.append(_y)
            self.z.append(_z)
            self.ke.append(_ke)
            self.drift.append(_drift)
            self.ionized.append(_ionized)
            self.n += 1

    def M(self, recursive=False):
        M_over_d = 0
        for i in range(self.n):
            if not recursive:
                M_over_d += (2**(self.ionized[i][-1] * 10 / self.x[i][-1])) # If using estimated number of ionization electrons
            else:
                M_over_d += (self.ionized[i][-1] + 1)**(10 / self.x[i][-1]) # If using actual number of ionization electrons
        M_over_d /= self.n
        return M_over_d

    def alpha(self):
        alph = 0
        for i in range(self.n):
            alph += self.ionized[i][-1] * log(2) / (self.x[i][-1] * 1e-4)
        alph /= self.n
        return alph

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
    
    def drift_velocity_plot(self):
        return self.t, self.drift
    
    def energy_plot(self):
        return self.ke