import sys
import getopt
import csv
import numpy as np
import glob
import re
from math import sqrt, log

import scipy.interpolate as interpolate
import matplotlib.pyplot as plt
from scipy import signal

def smooth_sig(xn, order, cutoff):
    b, a = signal.butter(order, cutoff)
    zi = signal.lfilter_zi(b, a)
    z, _ = signal.lfilter(b, a, xn, zi=zi*xn[0])
    z2, _ = signal.lfilter(b, a, z, zi=zi*z[0])
    y = signal.filtfilt(b, a, xn)
    return y

class Graph:
    def __init__(self, file_list, smooth=False, order=3, cutoff=0.01, nvar=10):
        self.n = 0
        self.tail_percent = 0.05
        
        self.t = []
        self.x = []
        self.y = []
        self.z = []
        self.ke = []
        self.drift = []
        self.angle = []
        self.ionized = []
        self.distance = [] # distance since the previous step
        self.interactions = []
        self.nvar = nvar

        for file in file_list:
            #print ('file name : ',file)
            #_angle = 0.
            if (self.nvar == 10):
                _t,  _x, _y, _z, _ke, _drift, _angle, _distance, _interactions, _ionized = np.loadtxt(file, delimiter=',', unpack=True, ndmin=2)
            #_t, _x, _y, _z, _ke, _drift, _angle, _ionized, _distance = np.loadtxt(file, delimiter=',', unpack=True, ndmin=2)
            if (self.nvar == 7):
                _t, _x, _y, _z, _ke, _drift,  _ionized = np.loadtxt(file, delimiter=',', unpack=True, ndmin=2)

            if smooth is True:
                self.t.append(smooth_sig(_t, order, cutoff))
                self.x.append(smooth_sig(_x, order, cutoff))
                self.y.append(smooth_sig(_y, order, cutoff))
                self.z.append(smooth_sig(_z, order, cutoff))
                self.ke.append(smooth_sig(_ke, order, cutoff))
                self.drift.append(smooth_sig(_drift, order, cutoff))
                self.ionized.append(smooth_sig(_ionized, order, cutoff))
                if (self.nvar == 10):
                    self.angle.append(smooth_sig(_angle, order, cutoff))
                    self.distance.append(smooth_sig(_distance, order, cutoff))
                    self.interactions.append(smooth_sig(_interactions, order, cutoff))

            else:
                self.t.append(_t)
                self.x.append(_x)
                self.y.append(_y)
                self.z.append(_z)
                self.ke.append(_ke)
                self.drift.append(_drift)
                self.ionized.append(_ionized)

                if (self.nvar == 10):
                    self.angle.append(_angle)
                    self.distance.append(_distance)
                    self.interactions.append(_interactions)

                self.n += 1

    def M(self, recursive=False):
        avg_ending_dist = sum([(self.x[i][-1]**2 + self.y[i][-1]**2 + self.z[i][-1]**2)**0.5 for i in range(self.n)]) / self.n

        if (recursive):
            mult_factor = sum([self.ionized[k][-1] + 1 for k in range(self.n)]) / self.n
            return mult_factor

        d = 10
        est_ending_elecs_log = np.log(2) * np.array([k[-1] for k in self.ionized])
        mu = np.mean(est_ending_elecs_log)
        variance = np.var(est_ending_elecs_log)
        mult_factor = np.exp((mu + 0.5 * variance) * (d / avg_ending_dist))

        return mult_factor

    def alpha(self, recursive=False):
        avg_ending_dist = sum([(self.x[i][-1]**2 + self.y[i][-1]**2 + self.z[i][-1]**2)**0.5 for i in range(self.n)]) / self.n

        if (recursive):
            return np.log(self.M(recursive)) / avg_ending_distance

        est_ending_elecs_log = np.log(2) * np.array([k[-1] for k in self.ionized])
        mu = np.mean(est_ending_elecs_log)
        variance = np.var(est_ending_elecs_log)

        return (mu + 0.5 * variance) / (avg_ending_dist * 1e-4)

    def drift_mean(self):
        mean = sum([sum(k[-int(len(k) * self.tail_percent):]) / int(len(k) * self.tail_percent) for k in self.drift]) / self.n
        return mean

    def drift_std_dev(self):
        mean = self.drift_mean()
        mean2 = sum([sum([i**2 for i in k[-int(len(k) * self.tail_percent):]]) / int(len(k) * self.tail_percent) for k in self.drift]) / self.n
        return sqrt((self.n / (self.n - 1)) * mean2 - mean**2)
