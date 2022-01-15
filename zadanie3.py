from numpy import *

a = 6378137
e2 = 0.00669437999013
b = a * sqrt((1 - e2))

def naRad(st, m=0.0, s=0.0):
    return (st + m/60 + s/3600) * pi/180

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

def kivioj(fi, la, s, az):
    n = int(s/1000)
    ds = s / n
    for i in range(0, n):
        M = liczM(fi)
        N = liczN(fi)
        Dfi = ds * cos(az) / M
        fiM = fi + 0.5*Dfi
        M = liczM(fiM)
        N = liczN(fiM)
        azM = az + ds*sin(az)*tan(fiM)/N

        Dfi = ds*cos(azM)/M
        fi += Dfi
        Dla = ds*sin(azM)/(N*cos(fiM))
        la += Dla
        Daz = ds*sin(azM)*tan(fiM)/N
        az += Daz

    return fi, la, az+pi


def vincent(fiA, laA, fiB, laB):
    f = 1 - b / a
    dLa = laB - laA
    UA = arctan((1 - f) * tan(fiA))
    UB = arctan((1 - f) * tan(fiB))
    L = dLa

    while True:
        sinO = sqrt((cos(UB) * sin(L)) ** 2 + (cos(UA)*sin(UB) - sin(UA) * cos(UB) * cos(L)) ** 2)
        cosO = sin(UA) * sin(UB) + cos(UA) * cos(UB) * cos(L)
        O = arctan((sinO / cosO))
        sina = (cos(UA) * cos(UB) * sin(L)) / sinO
        cos2a = 1 - sina ** 2
        cos2Om = cosO - 2 * sin(UA) * sin(UB) / cos2a
        C = f * cos2a * (4 + f * (4 - 3 * cos2a)) / 16
        oldL = L
        L = dLa + (1 - C) * f * sina * (O + C * sinO * (cos2Om + C * cosO * (2 * cos2Om ** 2 - 1)))
        if abs(L - oldL) < naRad(0, 0, 0.000001):
            break

    u2 = (a**2 - b**2)*cos2a/b**2
    A = 1 + u2*(4096 + u2*(-768 + u2*(320 - 175*u2))) / 16384
    B = u2*(256 + u2*(-128 + u2*(74 - 47*u2))) / 1024
    dO = B*sinO*(cos2Om + (1/4)*B*(cosO*(-1 + 2*cos2Om**2) - (1/6)*B*cos2Om*(-3 + 4*sinO**2)*(-3 + 4*cos2Om**2)))
    s = b*A*(O - dO)

    def azymut(dy, dx):
        ret = arctan(dy/dx)
        if dx < 0:
            ret += pi
        elif dy < 0:
            ret += 2*pi
        return ret

    Az1 = azymut(cos(UB)*sin(L), (cos(UA)*sin(UB) - sin(UA)*cos(UB)*cos(L)))
    Az2 = azymut(cos(UA)*sin(L), (-sin(UA)*cos(UB) + cos(UA)*sin(UB)*cos(L))) + pi
    return round(s, 3), round(Az1, 6), round(Az2, 6)


if __name__ == '__main__':
    def dms(rad):
        dg = rad*180/pi
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
                return '0'+str(i)
            else:
                return '00'

        d = f(d)
        m = f(m)
        s = f(s)

        return d+'°'+m+"'"+s+"\""

    from random import randint
    nr = randint(0, 15)
    print('nr', nr)
    fiA = naRad(50, 15 + nr * 15)
    laA = naRad(20, 45)
    fiC = naRad(50, 15 + nr * 15)
    laC = naRad(21, 15)
    fiB = naRad(50, 0 + nr * 15)
    laB = naRad(20, 45)
    fiD = naRad(50, 0 + nr * 15)
    laD = naRad(21, 15)

    sAD, AzAD, AzDA = vincent(fiA, laA, fiD, laD)
    Fk, Lk, Ak = kivioj(fiA, laA, sAD/2, AzAD)

    print(f'pkt średni   fi: {dms((fiA + fiD) / 2)}, lam: {dms((laA+laD)/2)}')
    print(f'pkt środkowy fi: {dms(Fk)}, lam: {dms(Lk)}')

    d, A1, A2 = vincent((fiA+fiD)/2, (laA+laD)/2, Fk, Lk)
    print(f'odległość: {d}\nazymuty: {dms(A1)}, {dms(A2)}')
    print('pole:', pole(fiA, laA, fiD, laD))
