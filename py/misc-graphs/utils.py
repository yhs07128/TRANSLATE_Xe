import sys
import getopt
import scipy.interpolate as interpolate
from math import sqrt, log
import matplotlib.pyplot as plt
import csv
import numpy as np
import glob
import re
from scipy import signal

def smooth_sig(xn, order, cutoff):
    b, a = signal.butter(order, cutoff)
    zi = signal.lfilter_zi(b, a)
    z, _ = signal.lfilter(b, a, xn, zi=zi*xn[0])
    z2, _ = signal.lfilter(b, a, z, zi=zi*z[0])
    y = signal.filtfilt(b, a, xn)
    return y

class Graph:
    def __init__(self, file_list, smooth=False, order=3, cutoff=0.01):
        self.n = 0
        self.tail_percent = 0.05

        #self.t_end_avg = 0
        
        self.t = []
        self.x = []
        self.y = []
        self.z = []
        self.ke = []
        self.drift = []
        self.ionized = []

        for file in file_list:
            _t, _x, _y, _z, _ke, _drift, _ionized = np.loadtxt(file, delimiter=',', unpack=True, ndmin=1)
            if smooth is True:
                self.t.append(smooth_sig(_t, order, cutoff))
                self.x.append(smooth_sig(_x, order, cutoff))
                self.y.append(smooth_sig(_y, order, cutoff))
                self.z.append(smooth_sig(_z, order, cutoff))
                self.ke.append(smooth_sig(_ke, order, cutoff))
                self.drift.append(smooth_sig(_drift, order, cutoff))
                self.ionized.append(smooth_sig(_ionized, order, cutoff))
            else:
                self.t.append(_t)
                self.x.append(_x)
                self.y.append(_y)
                self.z.append(_z)
                self.ke.append(_ke)
                self.drift.append(_drift)
                self.ionized.append(_ionized)
            self.n += 1
            #self.t_end_avg += _t[-1]

        #self.t_end_avg /= self.n

    def M(self, recursive=False):
        '''
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
        '''

        d = 10

        est_ending_elecs_log = np.log(2) * np.array([k[-1] for k in self.ionized])
        mu = np.mean(est_ending_elecs_log)
        sigma2 = np.var(est_ending_elecs_log)
        avg_ending_dist = sum([(self.x[i][-1]**2 + self.y[i][-1]**2 + self.z[i][-1]**2)**0.5 for i in range(self.n)]) / self.n
        mult_factor = np.exp((mu + 0.5 * sigma2) * (d / avg_ending_dist))

        return mult_factor


    def alpha(self, recursive=False):
        '''
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
        '''

        avg_ending_dist = sum([(self.x[i][-1]**2 + self.y[i][-1]**2 + self.z[i][-1]**2)**0.5 for i in range(self.n)]) / self.n
        est_ending_elecs_log = np.log(2) * np.array([k[-1] for k in self.ionized])
        mu = np.mean(est_ending_elecs_log)
        sigma2 = np.var(est_ending_elecs_log)
        return (mu + 0.5 * sigma2) / (avg_ending_dist * 1e-4)


    def avg_x(self, k):
        return sum([self.x[i][k] for i in range(self.n)]) / self.n

    def variance_x(self, k):
        mean = self.avg_x(k)
        mean2 = sum([self.x[i][k]**2 for i in range(self.n)]) / self.n
        return mean2 - mean**2

    def avg_y(self, k):
        return sum([self.y[i][k] for i in range(self.n)]) / self.n

    def variance_y(self, k):
        mean = self.avg_y(k)
        mean2 = sum([self.y[i][k]**2 for i in range(self.n)]) / self.n
        return mean2 - mean**2

    def avg_z(self, k):
        return sum([self.z[i][k] for i in range(self.n)]) / self.n

    def variance_z(self, k):
        mean = self.avg_z(k)
        mean2 = sum([self.z[i][k]**2 for i in range(self.n)]) / self.n
        return mean2 - mean**2

    def avg_t(self, k):
        t = []
        for i in range(self.n):
            t.append(self.t[i][k])
        return sum(t) / len(t)

    '''
    def eV_mean(self, records):
        mean = sum([sum(k[-records:]) for k in self.ke]) / (self.n * records)
        return mean

    def eV_std_dev(self, records):
        mean = self.eV_mean(records)
        mean2 = sum([sum([i**2 for i in k[-records:]]) for k in self.ke]) / (self.n * records)
        return sqrt(mean2 - mean**2)
    '''

    def eV_mean(self):
        mean = sum([sum(k) / len(k) for k in self.ke]) / self.n
        return mean

    def eV_std_dev(self):
        mean = self.eV_mean()
        mean2 = sum([sum([i**2 for i in k]) / len(k) for k in self.ke]) / self.n
        return sqrt((self.n / (self.n - 1)) * mean2 - mean**2)

    def drift_mean(self):
        mean = sum([sum(k[-int(len(k) * self.tail_percent):]) / int(len(k) * self.tail_percent) for k in self.drift]) / self.n
        return mean

    def drift_std_dev(self):
        mean = self.drift_mean()
        mean2 = sum([sum([i**2 for i in k[-int(len(k) * self.tail_percent):]]) / int(len(k) * self.tail_percent) for k in self.drift]) / self.n
        return sqrt((self.n / (self.n - 1)) * mean2 - mean**2)
    
    def drift_velocity_plot(self):
        return self.t, self.drift
    
    def energy_plot(self):
        return self.ke