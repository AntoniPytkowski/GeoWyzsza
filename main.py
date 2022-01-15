from math import *
import numpy as np

# grs 80
e2 = 0.00669437999013
a = 6378137
b = a * (1 - e2) ** 0.5
f = 1 / 298.257222101

# krasowski
akr = 6378245
bkr = 6356863
e2kr = (akr ** 2 - bkr ** 2) / (akr ** 2)


def dms(rad):
    dg = rad * 180 / pi
    d = int(dg)
    dg -= d
    dg *= 60
    m = int(dg)
    dg -= m
    dg *= 60
    s = round(dg, 5)

    def f(i):
        if i >= 10:
            return str(i)
        elif 0 < i < 10:
            return '0' + str(i)
        else:
            return '00'

    d = f(d)
    m = f(m)
    s = f(s)

    return d + '°' + m + "'" + s + "\""


def to_xyz(fi, lam, h):
    N = a / (1 - e2 * sin(fi) ** 2) ** 0.5
    x = (N + h) * cos(fi) * cos(lam)
    y = (N + h) * cos(fi) * sin(lam)
    z = ((N * (1 - e2) + h) * sin(fi))
    return round(x, 3), round(y, 3), round(z, 3)


def Hirvonen(x, y, z):
    r = ((x ** 2) + (y ** 2)) ** 0.5
    fi1 = atan(z / r * (1 - e2kr) ** -1)
    N = akr / (1 - e2kr * sin(fi1) ** 2) ** 0.5
    h = r / cos(fi1) - N
    fi2 = atan(z / r * (1 - e2kr * N / (N + h)) ** -1)

    while abs(fi2 - fi1) > 2.42406840554768e-10:
        N = akr / (1 - e2kr * sin(fi1) ** 2) ** 0.5
        h = r / cos(fi1) - N
        fi1 = fi2
        fi2 = atan(z / r * (1 - e2kr * N / (N + h)) ** -1)

    lam = atan(y / x)
    N = akr / (1 - e2kr * sin(fi2) ** 2) ** 0.5
    h = r / cos(fi2) - N

    return fi2, lam, round(h, 3)


def BursyWolf(xp, yp, zp):
    K = 0.8407728e-6
    al = -1.786877784465417e-06
    be = -2.5612706773016787e-07
    g = 4.089597325631379e-06

    x0 = -33.4297
    y0 = 146.5746
    z0 = 76.2865

    m1 = np.array([xp, yp, zp])
    m2 = np.array(
        [[K, g, -be],
         [-g, K, al],
         [be, -al, K]
         ])
    m3 = np.array([x0, y0, z0])

    ret = m1 + m2 @ m1 + m3
    for i in range(0, len(ret)):
        ret[i] = round(ret[i], 3)
    return ret


if __name__ == '__main__':
    A = {'fi': 50.25, 'lam': 20.75}
    B = {'fi': 50.00, 'lam': 20.75}
    C = {'fi': 50.25, 'lam': 21.25}
    D = {'fi': 50.00, 'lam': 21.25}

    pktSredni = {'fi': 50.125, 'lam': 21.00}
    pktSrodkowy = {'fi': 50.12527044824481, 'lam': 21.000651088258433}
    points = [A, B, C, D, pktSredni, pktSrodkowy]
    for P in points:
        print('punkt:', P)
        print("elipsoida grs 80:")
        x, y, z = to_xyz(radians(P['fi']), radians(P['lam']), 100)
        print(x, y, z)
        print("elipsoida krasowskiego:")
        xk, yk, zk = BursyWolf(x, y, z)
        print(xk, yk, zk)
        print('fi, lambda, h:')
        fp, lp, hp = Hirvonen(xk, yk, zk)
        print(dms(fp), dms(lp), hp)
