# comparison of spectral properties of schemes
# Chiel van Heerwaarden, 2011

import numpy
from pylab import *


# first derivative collocated
def modk(alpha, a, b, c, kx):
    kxmod = (a * sin(kx) + b / 2.0 * sin(2.0 * kx) + c / 3.0 * sin(3.0 * kx)) / (
        1.0 + 2.0 * alpha * cos(kx)
    )
    return kxmod


# first derivative staggered
def modks(alpha, a, b, c, kx):
    kxmod = (
        a * sin(1.0 / 2.0 * kx)
        + b / 3.0 * sin(3.0 / 2.0 * kx)
        + c / 5.0 * sin(5.0 / 2.0 * kx)
    ) / (1.0 / 2.0 + alpha * cos(kx))
    return kxmod


# transfer function interpolation
def itf(alpha, a, b, c, kx):
    tf = (
        a * cos(1.0 / 2.0 * kx) + b * cos(3.0 / 2.0 * kx) + c * cos(5.0 / 2.0 * kx)
    ) / (1.0 + 2.0 * alpha * cos(kx))
    return tf


kx = arange(0, pi + 0.0001, pi / 100.0)

dx2nd = modk(0.0, 1.0, 0.0, 0.0, kx)
dx2nds = modks(0.0, 1.0, 0.0, 0.0, kx)
dx4th = modk(0.0, 4.0 / 3.0, -1.0 / 3.0, 0.0, kx)
dx4ths = modks(0.0, 9.0 / 8.0, -1.0 / 8.0, 0.0, kx)
dx6thc = modk(1.0 / 3.0, 14.0 / 9.0, 1.0 / 9.0, 0.0, kx)

int2nd = itf(0.0, 1.0, 0.0, 0.0, kx)
int4th = itf(0.0, 9.0 / 8.0, -1.0 / 8.0, 0.0, kx)
int6th = itf(0.0, 75.0 / 64.0, -25.0 / 128.0, 3.0 / 128.0, kx)
int6thwrf = itf(0.0, 74.0 / 60.0, -16.0 / 60.0, 2.0 / 60.0, kx)

close("all")
figure(1)
plot(kx, kx, "k-")
plot(kx, dx2nd, label="2nd")
plot(kx, dx2nds, label="2nds")
plot(kx, dx4th, label="4th")
plot(kx, dx4ths, label="4ths")
plot(kx, dx6thc, label="6thc")
xlim(0, pi)
ylim(0, pi)
xlabel("k dx")
ylabel("k* dx")
legend(loc=0, frameon=False)

figure(2)
plot(kx, kx / kx, "k-")
plot(kx, dx2nd / kx, label="2nd")
plot(kx, dx2nds / kx, label="2nds")
plot(kx, dx4th / kx, label="4th")
plot(kx, dx4ths / kx, label="4ths")
plot(kx, dx6thc / kx, label="6thc")
xlim(0, pi)
ylim(0.0, 1.1)
xlabel("k dx")
ylabel("k*/k")
legend(loc=3, frameon=False)

figure(3)
plot(kx, kx / kx, "k-")
plot(kx, int2nd, label="2nds")
plot(kx, int4th, label="4ths")
plot(kx, int6th, label="6ths")
# plot(kx, int6thwrf, label='6thswrf')
xlim(0, pi)
ylim(0, 1.1)
xlabel("k dx")
ylabel("T")
legend(loc=3, frameon=False)

ming2i2 = numpy.min(array([int2nd, dx2nds / kx]), 0)
ming2i4 = numpy.min(array([int4th, dx2nds / kx]), 0)
ming4i4 = numpy.min(array([int4th, dx4ths / kx]), 0)
ming4i6 = numpy.min(array([int6th, dx4ths / kx]), 0)

figure(4)
plot(kx, kx / kx, "k-")
plot(kx, int2nd, label="int2nds")
plot(kx, int4th, label="int4ths")
plot(kx, int6th, label="int6ths")
plot(kx, dx2nds / kx, label="diff2nds")
plot(kx, dx4ths / kx, label="diff4ths")
plot(kx, dx6thc / kx, label="diff6thc")
xlabel("k dx")
ylabel("accuracy")
xlim(0, pi)
ylim(0, 1.1)
legend(loc=3, frameon=False)

figure(5)
plot(2.0 * pi / kx, kx / kx, "k-")
plot(2.0 * pi / kx, int2nd, label="int2nds")
plot(2.0 * pi / kx, int4th, label="int4ths")
plot(2.0 * pi / kx, int6th, label="int6ths")
plot(2.0 * pi / kx, dx2nds / kx, label="diff2nds")
plot(2.0 * pi / kx, dx4ths / kx, label="diff4ths")
plot(2.0 * pi / kx, dx6thc / kx, label="diff6thc")
xlabel("grid points per wave")
ylabel("accuracy")
xlim(2, 25)
ylim(0.95, 1.0)
grid()
legend(loc=0, frameon=False)
