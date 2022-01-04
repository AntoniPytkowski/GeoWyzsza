# todo: dostajemy dane fi la
# todo: najpierw przeliczamy na G-K
# todo: mamy przeliczyc na 2k i 92 z G-K
# todo: i w druga strone xd

from math import *

# globals
a = 6378137
e2 = 0.00669437999013
b = sqrt((1 - e2) * a**2)
ep2 = (a**2 - b**2)/b**2


def liczM(fi):
    return a * (1 - e2) / (sqrt(1 - e2 * (sin(fi)) ** 2)) ** 3
def liczN(fi):
    return a / sqrt(1 - e2 * (sin(fi)) ** 2)

def naRad(st, m=0.0, s=0.0):
    r = st + m / 60 + s / 3600
    return r * pi / 180


# dane
phi = naRad(51, 42, 3.70702)
lam = naRad(18, 10, 31.66445)
print(f'phi: {round(phi*180/pi, 12)}\nlam: {round(lam*180/pi, 12)}')

class Sigma:
    def __init__(self, fi=None):
        self.A0 = 1 - e2/4 - 3*e2**2/64 - 5*e2**3/256
        self.A2 = 3*(e2 + e2**2/4 + 15*e2**3/128)/8
        self.A4 = 15*(e2**2 + 3*e2**3/4)/256
        self.A6 = 35*e2**3/3072
        if fi is not None:
            self.sigma = a*(self.A0*fi - self.A2*sin(2*fi) + self.A4*sin(4*fi) - self.A6*sin(6*fi))

def GaussKruger(fi, la, L0, acc=3, show=False):
    L0 *= pi/180
    o = Sigma(fi).sigma
    t = tan(fi)
    n2 = ep2 * cos(fi)**2
    l = la - L0
    Xgk = o + 0.5*l**2*liczN(fi)*sin(fi)*cos(fi)*(1 + l**2/12*cos(fi)**2*(5 - t**2 + 9*n2 + 4*n2**2) + l**4/360*cos(fi)**4*(61 - 58*t**2 + t**4 + 270*n2 - 330*n2*t**2))
    Ygk = l*liczN(fi)*cos(fi)*(1 + l**2/6*cos(fi)**2*(1 - t**2 + n2) + l**4/120*cos(fi)**4*(5 - 18*t**2 + t**4 + 14*n2 - 58*n2*t**2))
    Xgk = round(Xgk, acc)
    Ygk = round(Ygk, acc)
    if show:
        print(f'Xgk: {Xgk}\nYgk: {Ygk}')
    return Xgk, Ygk


def FLto92(fi, la, acc=3, show=False):
    Xgk, Ygk = GaussKruger(fi, la, 19, show=show)
    x92 = 0.9993*Xgk - 5300000
    y92 = 0.9993*Ygk + 500000
    x92 = round(x92, acc)
    y92 = round(y92, acc)
    if show:
        print(f'x92: {x92}\ny92: {y92}')
    return x92, y92  # obie wartosci sa za duze o tyle samo: 0.000999...


def strefa(la):
    la *= 180/pi
    la += 0.5
    la = int(la)
    return la/3

def FLto2k(fi, la, acc=3, show=False):
    nrS = strefa(lam)
    Xgk, Ygk = GaussKruger(fi, la, nrS*3, show=show)
    x2k = 0.999923*Xgk
    y2k = 0.999923*Ygk + 500000 + nrS*1000000
    x2k = round(x2k, acc)
    y2k = round(y2k, acc)
    if show:
        print(f'x2k: {x2k}\ny2k: {y2k}')
    return x2k, y2k


def GKtoFL(x, y, L0, acc=12, show=False):
    mianow = a*Sigma().A0
    fi1 = x/mianow
    o = Sigma(fi1).sigma

    eps = naRad(0, 0, 0.000001)
    while True:
        old = fi1
        fi1 += (x - o)/mianow
        o = Sigma(fi1).sigma

        if abs(old-fi1) < eps:
            break

    N = liczN(fi1)
    M = liczM(fi1)
    t = tan(fi1)
    n2 = ep2 * cos(fi1) ** 2

    fi = fi1 - y**2*t/(2*M*N)*(1 - y**2/(12*N**2)*(5+3*t**2+n2-9*n2*t**2-4*n2**2) + y**4/(360*N**4)*(61+90*t**2+45*t**4))
    la = L0*pi/180 + y/(N*cos(fi1))*(1 - y**2/(6*N**2)*(1+2*t**2+n2) + y**4/(120*N**4)*(5+28*t**2+24*t**4+6*n2+8*n2*t**2))
    fi *= 180/pi
    la *= 180/pi

    fi = round(fi, acc)
    la = round(la, acc)
    if show:
        print(f'fi: {fi}\nla: {la}')
    return fi, la


def u92toFL(x, y, acc=12, show=False):
    Xgk = (x+5300000)/0.9993
    Ygk = (y-500000)/0.9993
    return GKtoFL(Xgk, Ygk, 19, acc=acc, show=show)


# TODO
def u2ktoFL(x, y, acc=12, show=False):
    Xgk = x / 0.999923
    L0 = int((y-500000)/1000000)*3
    Ygk = (y - 500000) % 1000000 / 0.999923
    return GKtoFL(Xgk, Ygk, L0, acc=acc, show=show)


xx, yy = FLto2k(phi, lam, show=True)
ff, ll = u2ktoFL(xx, yy, show=True)

# φ=51°42'3.70702"
# λ=18°10'31.66445"
#
# współrzędne przeliczone:
#
# 2000:
# x=5729652.129 y=6512129.578
#
# 1992:
# x=426389.392 y=443036.288
