import math

import matplotlib.pyplot as plt

d = 0.530 - 2 * 0.008
ro = 890
T0 = 63 + 273
v1 = 75 * 10 ** (-6)
T1 = 4 + 273
v2 = 35 * 10 ** (-6)
T2 = 20 + 273
delta = 0.2 * 10 ** (-3)
py = 20 * 10 ** 3
Kt = 3.5
Cv = 2000
Tn = 10 + 273

k = 1/(T2 - T1) * math.log(v1/v2)

v0 = v1 * math.exp((-k) * (T0 - T1))

x = [0, 20000, 40000, 60000, 80000, 100000, 120000]
z = [100, 325, 280, 275, 225, 195, 50]

p0 = 4.7 * 10 ** 6
pk = 0.3 * 10 ** 6

H0 = p0 / (ro * 9.81) + z[0]
Hk = pk / (ro * 9.81) + z[6]


def func_f(q: float):

    speed = 4 * q / (math.pi * d ** 2)
    Re = speed * d / v0
    eps = delta / d

    if Re < 10 / eps:
        lamb0 = 0.3164 / Re ** (1 / 4)
    elif 10 / eps < Re < 500 / eps:
        lamb0 = 0.11 * (68 / Re + eps) ** (1 / 4)
    else:
        lamb0 = 0.11 * (eps) ** (1 / 4)

    ik = lamb0 * speed ** 2 / (d * 2 * 9.81)
    H = []
    H.insert(0, pk / (ro * 9.81) + z[6])

    for i in range(len(x) - 1):
        H.insert(0, H[0] + ik * (x[6 - i] - x[5 - i]))
        if H[0] < z[5 - i] + (py / (ro * 9.81)):
            H[0] = (py / (ro * 9.81)) + z[5 - i]

    T = []
    T.append(T0)

    for i in range(len(x) - 1):
        vx = v1 * math.exp((-k) * (T0 - T[i]))
        Re = speed * d / vx
        eps = delta / d
        if Re < 10 / eps:
            lambx = 0.3164 / Re ** (1 / 4)
        elif 10 / eps < Re < 500 / eps:
            lambx = (0.11 * (68 / Re + eps) ** (1 / 4))
        else:
            lambx = (0.11 * (eps) ** (1 / 4))
        i0x = (lambx / d * speed ** 2 / (2 * 9.82))
        T.append(T[i] + (- math.pi * Kt * d / (ro * q * Cv) * (T[i] - Tn) + 9.81 * i0x / Cv) * (x[i + 1] - x[i]))

    Hi = []
    xn = []
    Hi.insert(0, H[6])
    xn.insert(0, x[6])

    for i in range(len(x) - 1):
        Ts = Tn + (T[len(x) - 1] - T[len(x) - 2])/(math.log((T[len(x) - 1]-Tn)/(T[len(x) - 2]-Tn)))
        vs = v1 * math.exp((-k) * (T0 - Ts))
        Re = speed * d / vs
        eps = delta / d
        if Re < 10 / eps:
            lambs = 0.3164 / Re ** (1 / 4)
        elif 10 / eps < Re < 500 / eps:
            lambs = 0.11 * (68 / Re + eps) ** (1 / 4)
        else:
            lambs = 0.11 * (eps) ** (1 / 4)
        i0s = lambs * speed ** 2 / (d * 2 * 9.81)
        Hi.insert(0, Hi[0] + i0s * (x[6 - i] - x[5 - i]))
        xn.insert(0, x[5 - i])
        if Hi[0] < z[5 - i] + (py / (ro * 9.81)):
            Hi[0] = (py / (ro * 9.81)) + z[5 - i]

    return Hi[0] - H0


Qn = 0.1
Qk = 2.5
Qs = 0
FQs = 10

while abs(FQs) >= 10 ** -10:
    FQn = func_f(Qn)
    FQk = func_f(Qk)
    if FQn * FQk >= 0:
        print( f'Нет решения на заданном отрезке')
        sys.exit()
    Qs = (Qn + Qk) / 2
    FQs = func_f(Qs)
    if FQs * FQn < 0:
        Qk = Qs
    else:
        Qn = Qs

print(f'Расход перекачки равен {round(Qs, 3)} м3/с или {round(Qs*3600, 3)} м3/ч')
speeds = 4 * Qs / (math.pi * d ** 2)
Tnew = []
Tnew.append(T0)

for i in range(len(x) - 1):
    vx = v1 * math.exp((-k) * (T0 - Tnew[i]))
    Re = speeds * d / vx
    eps = delta / d
    if Re < 10 / eps:
        lambx = 0.3164 / Re ** (1 / 4)
    elif 10 / eps < Re < 500 / eps:
        lambx = (0.11 * (68 / Re + eps) ** (1 / 4))
    else:
        lambx = (0.11 * (eps) ** (1 / 4))
    i0s = (lambx / d * speeds ** 2 / (2 * 9.82))
    Tnew.append(Tnew[i] + (- math.pi * Kt * d / (ro * Qs * Cv) * (Tnew[i] - Tn) + 9.81 * i0s / Cv) * (x[i + 1] - x[i]))

Hnew = []
xnew = []
znew = []
Hnew.insert(0, pk / (ro * 9.81) + z[6])
xnew.insert(0, x[6])
znew.insert(0, z[6])
nsam = 0


def sam(i_s, fi: float):
    f_lev = (fi - math.sin(fi)) ** (5 / 3) / (fi ** (2 / 3))
    f_prav = Qs / (4.137 * d ** (8 / 3) * math.sqrt(abs(i_s)))
    return f_lev - f_prav


def s_sam (i_s: float):
    fn = 1 ** -10
    fk = 2 * math.pi
    Fs = 1
    fs = 0
    Fn = sam(i_s, fn)
    Fk = sam(i_s, fk)

    while abs(Fs) >= 10 ** -10:
        if Fn * Fk > 0:
            print(f'Нет решения на заданном отрезке')
            sys.exit()
        fs = (fn + fk) / 2
        Fs = sam(i_s, fs)
        if Fs * Fn < 0:
            fk = fs
        else:
            fn = fs

    samot = fs

    return samot


for i in range(len(x)-1):
    Ts = Tn + (Tnew[len(x) - 1] - Tnew[len(x) - 2]) / (math.log((Tnew[len(x) - 1] - Tn) / (Tnew[len(x) - 2] - Tn)))
    vs = v1 * math.exp((-k) * (T0 - Ts))
    Re = speeds * d / vs
    eps = delta / d
    if Re < 10 / eps:
        lambs = 0.3164 / Re ** (1 / 4)
    elif 10 / eps < Re < 500 / eps:
        lambs = 0.11 * (68 / Re + eps) ** (1 / 4)
    else:
        lambs = 0.11 * (eps) ** (1 / 4)
    i0s = lambs * speeds ** 2 / (d * 2 * 9.81)
    Hnew.insert(0, Hnew[0] + i0s * (x[6 - i] - x[5 - i]))
    xnew.insert(0, x[5 - i])
    znew.insert(0, z[5 - i])
    if Hnew[0] < z[5 - i] + (py / (ro * 9.81)):
        Hnew[0] = (py / (ro * 9.81)) + z[5 - i]
        i_s = (z[5 - i] - z[6 -i ]) / (x[6 - i] - x[5 - i])
        X_s = (Hnew[0] - Hnew[1] + x[5 - i] * i_s - i0s * x[6 - i]) / (i_s - i0s)
        L_s = X_s - x[5 - i]
        znew.insert(1, znew[0] - i0s * L_s)
        Hnew.insert(1, Hnew[1] + i0s * (x[6 - i] - X_s))
        xnew.insert(1, X_s)
        S_sam = s_sam(float(i_s))
        stepen = (S_sam - math.sin(S_sam)) / (2 * math.pi)
        nsam = nsam + 1
        print(f'Самотечный участок на координатах {x[5 - i] / 1000} - {round(X_s / 1000, 2)} км протяженностью '
              f'{round(L_s / 1000, 2)} км и степенью заполнения {round(stepen * 100, 3)} %')

P = []

for k in range(len(Hnew)):
    P.append((Hnew[k] - znew[k]) * ro * 9.81)

Tnnew = [Tn, Tn, Tn, Tn, Tn, Tn, Tn]

plt.plot(x, z, '-o', c = 'b', mec='k', mfc='b', mew=1, ms=4, label='Профиль')
plt.plot(xnew, Hnew, '-.o', c = 'r', mec='k', mfc='r', mew=1, ms=4, label='Гидроуклон')
plt.title('График гидравлического уклона')
plt.xlabel('Координата, м')
plt.ylabel('Напор, м')
plt.legend(fontsize=10)
plt.show()

plt.subplot(1, 2, 1)
plt.plot(xnew, P, '-o', c = 'r', mec='k', mfc='r', mew=1, ms=4, label='Давление')
plt.title('График распределения давления')
plt.xlabel('Координата, м')
plt.ylabel('Давление, Па')
plt.legend(fontsize=8)

plt.subplot(1, 2, 2)
plt.plot(x, Tnew, '-o', c = 'r', mec='k', mfc='r', mew=1, ms=4, label='Температура')
plt.plot(x, Tnnew, '-.o', c = 'b', mec='k', mfc='r', mew=1, ms=4, label='Температура окружающая')
plt.title('График распределения температуры')
plt.xlabel('Координата, м')
plt.ylabel('Температура, К')
plt.legend(fontsize=8)

plt.show()

