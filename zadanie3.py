from numpy import *

a = 6378137
e2 = 0.00669437999013


# AD powinno byc 45295.37417179245

def naRad(st=0, m=0, s=0):
    r = st + m / 60 + s / 3600
    return r * pi / 180


def naSt(rad):
    return rad / pi * 180


###########
# VERSION 1
###########

def liczM(fi):
    return a * (1 - e2) / (sqrt(1 - e2 * (sin(fi)) ** 2)) ** 3


def liczN(fi):
    return a / sqrt(1 - e2 * (sin(fi)) ** 2)


def Kivioj(fi, la, s, az):
    ds = s / int(s)

    for i in range(0, int(s)):
        # N = liczN(fi)
        # M = liczM(fi)
        # nie wiem czy powinno być tak jak jest czy jak wyżej

        dB = ds * cos(az) / liczM(fi)
        fi += dB / 2
        az += ds * sin(az) * tan(fi) / liczN(fi)
        dB = ds * cos(az) / liczM(fi)
        fi += dB
        dL = ds * sin(az) / (liczN(fi) * cos(fi))
        la += dL
        dA = ds * sin(az) * tan(fi) / liczN(fi)
        az += dA

        if i + 5 > int(s):
            print(fi, la)

    dA = ds * sin(az) * tan(fi) / liczN(fi)
    az += dA
    az -= pi
    if az < 0:
        az += 2 * pi

    return fi, la, az


def Vincent(fiA, laA, fiB, laB):
    b = a * sqrt(1 - e2)
    f = 1 - b / a
    dLa = laB - laA
    UA = arctan((1 - f) * tan(fiA))
    UB = arctan((1 - f) * tan(fiB))
    L = dLa

    while True:
        sinO = sqrt((cos(UB) * sin(L)) ** 2 + (sin(UA) * cos(UB) * cos(L)) ** 2)
        cosO = sin(UA) * sin(UB) + cos(UA) * cos(UB) * cos(L)
        O = arctan((sinO / cosO))
        sina = (cos(UA) * cos(UB) * sin(L)) / sinO
        cos2a = 1 - sina ** 2
        cos2Om = cosO - 2 * sin(UA) * sin(UB) / cos2a
        C = f * cos2a * (4 + f * (4 - 3 * cos2a)) / 16
        oldL = L
        L = dLa + (1 - C) * f * sina * (O + C * sinO * (cos2Om + C * cosO * (2 * cos2Om ** 2 - 1)))
        if abs(L - oldL) < 4.84814e-12:
            break

    u2 = (a**2 - b**2)*cos2a/b**2
    A = 1 + u2*(4096 + u2*(-768 + u2*(320 - 175*u2))) / 16384
    B = u2*(256 + u2*(-128 + u2*(74 - 47*u2))) / 1024
    dO = B*sinO*(cos2Om + (1/4)*B*(cosO*(-1 + 2*cos2Om**2) - (1/6)*B*cos2Om*(-3 + 4*sinO**2)*(-3 + 4*cos2Om**2)))
    s = b*A*(O - dO)
    Az1 = arctan(cos(UB)*sin(L) / (cos(UA)*sin(UB) - cos(UA)*sin(UB)*cos(L)))
    Az2 = arctan((-1)*cos(UA)*sin(L) / (cos(UA)*sin(UB) - cos(UA)*sin(UB)*cos(L))) + pi
    return [s, Az1, Az2]


nr = 0
fiA = naRad(50, 15 + nr * 15)
laA = naRad(20, 45)
fiC = naRad(50, 15 + nr * 15)
laC = naRad(21, 15)
fiB = naRad(50, 0 + nr * 15)
laB = naRad(20, 45)
fiD = naRad(50, 0 + nr * 15)
laD = naRad(21, 15)

# Kivioj(fiA, laA, 100, naRad(45))
# print(fiA, laA)

# print(Vincent(fiA, laA, fiD, laD))

print(Vincent(naRad(10), naRad(30), naRad(20), naRad(40)))
# 1434648.68930 0.640047 3.882303
