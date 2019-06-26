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

        self.t_end_avg = 0
        
        self.t = []
        self.x = []
        self.y = []
        self.z = []
        self.ke = []
        self.drift = []
        self.ionized = []

        for file in file_list:
            _t, _x, _y, _z, _ke, _drift, _ionized = np.loadtxt(file, delimiter=',', unpack=True, ndmin=1)
            self.t.append(_t)
            self.x.append(_x)
            self.y.append(_y)
            self.z.append(_z)
            self.ke.append(_ke)
            self.drift.append(_drift)
            self.ionized.append(_ionized)
            self.n += 1
            self.t_end_avg += _t[-1]
        
        self.t_end_avg /= self.n

    def M(self, recursive=False):
        M_over_d = 0
        d = 10

        N = 0
        avg_dist = 0

        for i in range(self.n):
            N += self.ionized[i][-1]
            avg_dist += (self.x[i][-1]**2 + self.y[i][-1]**2 + self.z[i][-1]**2)**0.5

        N /= self.n
        avg_dist /= self.n

        if not recursive:
            M_over_d = 2**(N * (d / avg_dist))
        else:
            M_over_d = (N + 1)**(d / avg_dist)
        
        return M_over_d

    def alpha(self, recursive=False):
        alph = 0

        N = 0
        avg_dist = 0

        for i in range(self.n):
            N += self.ionized[i][-1]
            avg_dist += (self.x[i][-1]**2 + self.y[i][-1]**2 + self.z[i][-1]**2)**0.5

        N /= self.n
        avg_dist /= self.n

        if not recursive:
            alph = N * log(2) / (avg_dist * 1e-4)
        else:
            alph += log(N + 1) / (avg_dist * 1e-4)

        return alph

    def avg_x(self, k):
        x = []
        for i in range(self.n):
            x.append(self.x[i][k])
        return sum(x) / len(x)

    def variance_x(self, k):
        mean = self.avg_x(k)
        variance = []
        for i in range(self.n):
            variance.append((self.x[i][k] - mean)**2)
        return sum(variance) / (len(variance) - 1)

    def avg_y(self, k):
        y = []
        for i in range(self.n):
            y.append(self.y[i][k])
        return sum(y) / len(y)

    def variance_y(self, k):
        mean = self.avg_y(k)
        variance = []
        for i in range(self.n):
            variance.append((self.y[i][k] - mean)**2)
        return sum(variance) / (len(variance) - 1)

    def avg_z(self, k):
        z = []
        for i in range(self.n):
            z.append(self.z[i][k])
        return sum(z) / len(z)

    def variance_z(self, k):
        mean = self.avg_z(k)
        variance = []
        for i in range(self.n):
            variance.append((self.z[i][k] - mean)**2)
        return sum(variance) / (len(variance) - 1)

    def avg_t(self, k):
        t = []
        for i in range(self.n):
            t.append(self.t[i][k])
        return sum(t) / len(t)

    def eV_mean(self):
        mean = 0
        for i in range(self.n):
            mean += sum(self.ke[i]) / len(self.ke[i])
        mean /= self.n
        return mean

    def eV_std_dev(self):
        mean = self.eV_mean()
        sigma = 0
        for i in range(self.n):
            sigma += (sum(self.ke[i]) / len(self.ke[i]) - mean)**2
        sigma /= (self.n - 1)
        return sqrt(sigma)

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