from math import *

# globals
a = 6378137
e2 = 0.00669437999013
b = sqrt((1 - e2) * a**2)
ep2 = (a**2 - b**2)/b**2

def dms(rad):
    dg = rad*180/pi
    d = int(dg); dg -= d; dg *= 60
    m = int(dg); dg -= m; dg *= 60
    s = round(dg, 5)
    def f(i):
        if i >= 10:
            return str(i)
        elif 0 < i < 10:
            return '0'+str(i)
        else:
            return '00'
    d = f(d); m = f(m); s = f(s)
    return d+'°'+m+"'"+s+"\""


def liczM(fi):
    return a * (1 - e2) / (sqrt(1 - e2 * (sin(fi)) ** 2)) ** 3
def liczN(fi):
    return a / sqrt(1 - e2 * (sin(fi)) ** 2)

def pole(f1, l1, f2, l2):
    def c(f):
        e1 = sqrt(e2)
        sf = sin(f)
        return sf/(1-e2*sf**2) + (1/(2*e1))*log((1+e1*sf)/(1-e1*sf))

    return round(abs(((b**2)*(l2-l1)/2) * (c(f2) - c(f1))), 3)

def naRad(st, m=0.0, s=0.0):
    r = st + m / 60 + s / 3600
    return r * pi / 180


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
        print(f'Xgk: {Xgk}, Ygk: {Ygk}')
    return Xgk, Ygk


def FLto92(fi, la, acc=3, show=False):
    Xgk, Ygk = GaussKruger(fi, la, 19, show=show)
    x92 = 0.9993*Xgk - 5300000
    y92 = 0.9993*Ygk + 500000
    x92 = round(x92, acc)
    y92 = round(y92, acc)
    if show:
        print(f'x92: {x92}, y92: {y92}')
    return x92, y92


def strefa(la):
    la *= 180/pi
    la += 0.5
    la = int(la)
    return la/3

def FLto2k(fi, la, acc=3, show=False):
    nrS = strefa(la)
    Xgk, Ygk = GaussKruger(fi, la, nrS*3)
    x2k = 0.999923*Xgk
    y2k = 0.999923*Ygk + 500000 + nrS*1000000
    x2k = round(x2k, acc)
    y2k = round(y2k, acc)
    if show:
        print(f'x2k: {x2k}, y2k: {y2k}')
    return x2k, y2k


def GKtoFL(x, y, L0, show=False):
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
    fi = dms(fi)
    la = dms(la)

    if show:
        print(f'fi: {fi}, la: {la}')
    return fi, la


def u92toFL(x, y, show=False):
    Xgk = (x+5300000)/0.9993
    Ygk = (y-500000)/0.9993
    return GKtoFL(Xgk, Ygk, 19, show=show)


def u2ktoFL(x, y, show=False):
    Xgk = x / 0.999923
    L0 = int((y-500000)/1000000)*3
    Ygk = (y - 500000) % 1000000 / 0.999923
    return GKtoFL(Xgk, Ygk, L0, show=show)


if __name__ == '__main__':
    nr = 9  # z poprzedniego zadania
    fiA = naRad(50, 15 + nr * 15)
    laA = naRad(20, 45)
    fiC = naRad(50, 15 + nr * 15)
    laC = naRad(21, 15)
    fiB = naRad(50, 0 + nr * 15)
    laB = naRad(20, 45)
    fiD = naRad(50, 0 + nr * 15)
    laD = naRad(21, 15)

    fiSrednie = naRad(52, 22, 30)
    laSrednie = naRad(21)
    fiSrodkowe = naRad(52, 22, 30.91097)
    laSrodkowe = naRad(21, 0, 02.48014)
    
    class Punkt:
        def __init__(self, fi, la):
            self.fi, self.la = fi, la
            self.xgk, self.ygk = 0.0, 0.0
            self.x92, self.y92 = 0.0, 0.0
            self.x2k, self.y2k = 0.0, 0.0
            

    pkts = [Punkt(fiA, laA), Punkt(fiB, laB), Punkt(fiC, laC), Punkt(fiD, laD), 
            Punkt(fiSrednie, laSrednie), Punkt(fiSrodkowe, laSrodkowe)]

    for p in pkts:
        print(f"fi: {dms(p.fi)}, la: {dms(p.la)}")
        p.x92, p.y92 = FLto92(p.fi, p.la, show=True)
        p.x2k, p.y2k = FLto2k(p.fi, p.la, show=True)
    print('')  # \n

    print('poleFL:', 947260271.645)
    pkts[0].xgk, pkts[0].ygk = GaussKruger(fiA, laA, 19)
    pkts[3].xgk, pkts[3].ygk = GaussKruger(fiD, laD, 19)

    def simplePole(x1, y1, x2, y2):
        return round(abs((x2-x1)*(y2-y1)), 3)

    print('poleGK:', simplePole(pkts[0].xgk, pkts[0].ygk, pkts[3].xgk, pkts[3].ygk))
    print('pole92:', simplePole(pkts[0].x92, pkts[0].y92, pkts[3].x92, pkts[3].y92))
    print('pole2k:', simplePole(pkts[0].x2k, pkts[0].y2k, pkts[3].x2k, pkts[3].y2k))
    print('')  # \n

    def mkappa(m0, L0):
        B, L = pkts[-1].fi, pkts[-1].la
        t2 = tan(B)**2
        l = L - L0
        n2 = ep2 * cos(B)**2
        m = 1 + l**2/2 * cos(B)**2 * (1+n2) + l**4/24 * cos(B)**4 * (5-4*t2)
        m *= m0
        K = (1 - m) * 1000
        return round(m, 3), round(K, 3)

    mgk, Kgk = mkappa(1, 19)
    m92, K92 = mkappa(0.9993, 19)
    m2k, K2k = mkappa(0.999923, 3*strefa(pkts[-1].la))
    print(f'mgk: {mgk}, Kgk: {Kgk}')
    print(f'm92: {m92}, K92: {K92}')
    print(f'm2k: {m2k}, K2k: {K2k}')
    print('')  # \n
    print(f'm^2 gk: {mgk**2}, K^2 gk: {(1-mgk**2)*10000}')
    print(f'm^2 92: {m92**2}, K^2 92: {(1-m92**2)*10000}')
    print(f'm^2 2k: {m2k**2}, K^2 2k: {(1-m2k**2)*10000}')
    
