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

        self.x_bounds = [np.inf, -np.inf]
        self.y_bounds = [np.inf, -np.inf]
        self.z_bounds = [np.inf, -np.inf]

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

            min_x = np.amin(_x)
            max_x = np.amax(_x)

            if min_x < self.x_bounds[0]:
                self.x_bounds[0] = min_x
            if max_x > self.x_bounds[1]:
                self.x_bounds[1] = max_x

            min_y = np.amin(_y)
            max_y = np.amax(_y)

            if min_y < self.y_bounds[0]:
                self.y_bounds[0] = min_y
            if max_y > self.y_bounds[1]:
                self.y_bounds[1] = max_y

            min_z = np.amin(_z)
            max_z = np.amax(_z)

            if min_z < self.z_bounds[0]:
                self.z_bounds[0] = min_z
            if max_z > self.z_bounds[1]:
                self.z_bounds[1] = max_z

            self.t_end_avg += _t[-1]
        
        self.t_end_avg /= self.n

    def area_xy(self):
        x_width = self.x_bounds[1] - self.x_bounds[0] * 1e-6
        y_width = self.y_bounds[1] - self.y_bounds[0] * 1e-6
        return x_width * y_width

    def area_xz(self):
        x_width = self.x_bounds[1] - self.x_bounds[0] * 1e-6
        z_width = self.z_bounds[1] - self.z_bounds[0] * 1e-6
        return x_width * z_width

    def area_yz(self):
        y_width = self.y_bounds[1] - self.y_bounds[0] * 1e-6
        z_width = self.z_bounds[1] - self.z_bounds[0] * 1e-6
        return y_width * z_width

    def diff_xy(self):
        return self.area_xy() / self.t_end_avg

    def diff_xz(self):
        return self.area_xz() / self.t_end_avg

    def diff_yz(self):
        return self.area_yz() / self.t_end_avg

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
            M_over_d = N**(d / avg_dist)
        
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
            alph += log(N) / (avg_dist * 1e-4)

        return alph

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