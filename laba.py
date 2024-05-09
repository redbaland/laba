import math
import matplotlib.pyplot as plt
import numpy as np

#инициализация всех входных данных
PROFILE=[150.2, 153.3, 171.6, 96, 125.1,155,158,134.5,134.5,112.7,83.1,102.2,151.4,109.5,148.4,83.2,180.1,180.1,191.8,177, 77.1,161.3, 203.3,118.3,171.5,203.1]
DLINA=[0,20,40,60,80,100,120,140,140,160,180,200,220,240,260,280,300,300,320,340,360,380,400,420,440,460]
DLINA[:] = [x * 1000 for x in DLINA]
NPS2 = 7
NPS3 = 16
n = 26
Pmax = 6.3*10**(6)
Do = 630/1000 # внешний 0.36371343685456636
b = 8/1000 #стенка
Di = Do-2*b
E = 0.2/1000/Di # шероховатость
plotnost = 850
viscosity = 15*10**(-6)
Py = 10*(10**3)
Hk = 35
Podpor = 40
Hmin = 40
Hmax = Pmax/(9.81*plotnost)
Head_max = [x+Hmax for x in PROFILE]
Hy = Py/(9.81*plotnost)
kavitacia = [0 for i in range(n)]
for i in range(len(PROFILE)-1):
    kavitacia[i] = Hy + PROFILE[i]
Head_x = [0 for i in range(n)]
Head_x[-1] = Hk+PROFILE[-1]
Head = [0 for i in range(n)]
Head[-1] = Hk+PROFILE[-1]
NUMBER_SAM = 0
FINAL_SAMOTEK = [[0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]]

def Pumps(Q_FROM_PUMP, HEAD_FROM_PUMP): #МНК для насосов
    a = np.vstack([Q_FROM_PUMP**2, np.ones_like(Q_FROM_PUMP)]).T
    c = np.linalg.lstsq(a, HEAD_FROM_PUMP, rcond=None)[0]
    c[0] = c[0]*(-1)
    return c

Q_FROM_PUMP2 = np.array([1000, 1250, 1500, 1750, 2000]) # НПС-2 - НМ-4
HEAD_FROM_PUMP2 = np.array([260, 253, 245, 235, 222])
NASOS4 = Pumps(Q_FROM_PUMP2, HEAD_FROM_PUMP2) # коэф аппроксимации

Q_FROM_PUMP1 = np.array([1000, 1250, 1500, 1750, 2000]) # НПС-3 И ГНПС - НМ-2
HEAD_FROM_PUMP1 = np.array([243, 237, 233, 225, 218])
NASOS2 = Pumps(Q_FROM_PUMP1, HEAD_FROM_PUMP1) # коэф аппроксимации


def PUMPS_HEAD(Q, coeff_aprocs): # функция расчета дифф напора насосов
    diff_head = coeff_aprocs[1]-coeff_aprocs[0]*(Q*3600)**2
    return diff_head

def Soprot(Q): # считаю Lumbda
    Re = (Q * 4) / (math.pi*viscosity*Di)
    if Re <= 2320:
        gidroSoprot = 64 / Re
    if 2320 <= Re <= 10 * E**(-1):
        gidroSoprot = 0.3164 / (Re)**0.25
    if 10 * E**(-1) < Re < 500 * E**(-1):
        gidroSoprot = 0.11 *(E + (68 / Re))**0.25
    if Re >= 500 * E**(-1):
        gidroSoprot = 0.11 * E**0.25
    return gidroSoprot


def Head_loses (distance,Q):
    lymbda = Soprot(Q)
    lose_h = lymbda*distance*8*Q**2/(9.81*Di**5*math.pi**2)
    return lose_h


def delenie_popolam(foo, a, b, tol=0, max_iter=500):
    if foo(a) * foo(b) >= 0:
        print("Метод бисекции не применим")
        return None
    iter_count = 0
    while (b - a) / 2 > tol and iter_count < max_iter:
        c = (a + b) / 2
        if foo(c) == 0:
            return c
        elif foo(c) * foo(a) < 0:
            b = c
        else:
            a = c
        iter_count += 1
    return (a + b) / 2


def FI(F):
#нахождение fi для расчета самотека
    alpha = abs(math.atan(sam_angle))
    a1 = (F - math.sin(F)) ** (5 / 3) / F ** (2 / 3)
    a2 = 0.2419 * Q_last / (Di ** (8 / 3) * math.sqrt(math.sin(alpha)))
    otkl = a1 - a2
    return otkl


def Paschet(Q):
    NO_ITER=0
    #перезапись массивов двойников для каждого захода в foo()
    Head = [0 for i in range(n)]
    Head[-1] = Hk + PROFILE[-1]
    L = [0 for i in range(n)]
    for i in range(0, len(DLINA)):
        L[i] = DLINA[i]
    Z = [0 for i in range(n)]
    for i in range(0, len(DLINA)):
        Z[i] = PROFILE[i]
    Kav = [0 for i in range(n)]
    for i in range(0, len(DLINA)):
        Kav[i] = kavitacia[i]
    MAX = [0 for i in range(n)]
    for i in range(0, len(DLINA)):
        MAX[i] = Head_max[i]


    for i in range(n-2, -1, -1):
        Head[i] = Head[i + 1] + Head_loses(L[i + 1] - L[i], Q)

        if i==NPS2: # рассчет входа в станцию NPS2
            Head[i] = Head[i+1]-2*PUMPS_HEAD(Q, NASOS4)

        if i==NPS3: # рассчет входа в станцию
            Head[i] = Head[i+1]-2*PUMPS_HEAD(Q, NASOS2)

        # Проверка на отрицательное давление
        if Head[i]<0:
            NO_ITER=1

        if Head[i] < kavitacia[i] and NO_ITER==0: #ЕСЛИ НАШЛИ САМОТЕК ПРОВОДИМ ДОМ ПОСТРОЕНИЕ ТОЧКИ
            if (i == 0 or i == NPS2 or i == NPS3):
                print("Нарушено условие подпора")
            else:
                x_sam = L[i + 1] + (L[i + 1] - L[i]) * (Head[i + 1] - Z[i + 1] - Hy) / (
                            (Z[i + 1] - Z[i]) + (L[i + 1] - L[i]) * Head_loses(1, Q))

                h_sam = (Z[i + 1] - Z[i]) / (L[i + 1] - L[i]) * (x_sam - L[i]) + Z[i]

                L.insert(i + 1, x_sam)

                Z.insert(i + 1, h_sam)

                Head.insert(i + 1, Hy + Z[i + 1])

                MAX.insert(i + 1, Hmax + Z[i + 1])

                Kav.insert(i + 1, Hy + Z[i + 1])

                Head[i] = Hy + Z[i]

    a1=Head[0]

    a2=Podpor + 2*PUMPS_HEAD(Q, NASOS2)+ PROFILE[0]

    otkl=a2-a1

    print("несходимость метода деление отрезка пополам =", otkl)

    return otkl


# Поиск реального расхода с учетом самотечности
Q_last = delenie_popolam(Paschet, 0.2, 1.2)

for i in range(n - 2, -1, -1):
    Head[i] = Head[i + 1] + Head_loses(DLINA[i + 1] - DLINA[i], Q_last)

    if i == NPS2:
        Head[i] = Head[i + 1] - 2 * PUMPS_HEAD(Q_last, NASOS4)

    if i == NPS3:
        Head[i] = Head[i + 1] - 2 * PUMPS_HEAD(Q_last, NASOS2)

    # Проверка на отрицательное давление
    if Head[i] < 0:
        NO_ITER = 1

    if Head[i] < kavitacia[i]:
       NUMBER_SAM+=1
       x_sam = DLINA[i+1] + (DLINA[i + 1] - DLINA[i]) * (Head[i + 1] - PROFILE[i + 1] - Hy) / ((PROFILE[i + 1] - PROFILE[i])+(DLINA[i + 1] - DLINA[i])*Head_loses(1,Q_last))

       h_sam = (PROFILE[i + 1] - PROFILE[i]) / (DLINA[i + 1] - DLINA[i]) * (x_sam - DLINA[i]) + PROFILE[i]

       DLINA.insert(i + 1, x_sam)

       PROFILE.insert(i + 1, h_sam)

       Head.insert(i + 1, Hy + PROFILE[i + 1])

       Head_max.insert(i + 1, Hmax + PROFILE[i + 1])

       kavitacia.insert(i + 1, Hy + PROFILE[i + 1])

       Head[i] = Hy + PROFILE[i]

       sam_angle = (PROFILE[i] - h_sam)/(DLINA[i] - x_sam)

       F = delenie_popolam(FI,0.1, 2*math.pi)

       PART = (F - math.sin(F))/(2*math.pi)*100 # cтепень заполнения

       FINAL_SAMOTEK[NUMBER_SAM-1] = [DLINA[i],x_sam,PROFILE[i],h_sam, PART]



print("Расход [м3/c] =", Q_last)

if NUMBER_SAM > 0:
    print("Количество самотечных участков = ", NUMBER_SAM)

    print("(координата конца, координата начала, высотка конца, высотка начала, степень заполнения): 1 = ", FINAL_SAMOTEK[0][:])

    print("(координата конца, координата начала, высотка конца, высотка начала, степень заполнения): 2 = ", FINAL_SAMOTEK[1][:])

    print("(координата конца, координата начала, высотка конца, высотка начала, степень заполнения): 3 = ", FINAL_SAMOTEK[2][:])

    print("(координата конца, координата начала, высотка конца, высотка начала, степень заполнения): 4 = ", FINAL_SAMOTEK[3][:])

    print("(координата конца, координата начала, высотка конца, высотка начала, степень заполнения): 5 = ", FINAL_SAMOTEK[4][:])


if NUMBER_SAM == 0:
    print("Самотечные участки отсутствуют")

# График ГУ
fig = plt.figure()
ax = fig.add_subplot()
ax.set_xlabel('Расстояние, км')
ax.set_ylabel('Напор, м')

plt.plot(DLINA, PROFILE, label = "Высотные отметки Z(X)")
plt.plot(DLINA, Head_max, label = "Несущая способность Hmax(X)")
plt.plot(DLINA, Head, label = "Напор H(X)")

plt.grid()
plt.legend()

# график зависимости давления от координаты
Pmax_mass = [Pmax/1000000 for i in range(len(Head))]
Pressure = [0 for i in range(len(Head))]
for i in range(len(Head)-1, -1, -1):
    Pressure[i] = (Head[i]-PROFILE[i])*(plotnost*9.81)/1000000
fig = plt.figure()
ax = fig.add_subplot()
ax.set_xlabel('Расстояние, км')
ax.set_ylabel('Давление, МПа')

plt.plot(DLINA, Pressure, label = "Давление P(X)")
plt.plot(DLINA, Pmax_mass, label = "Несущая способность Pmax(X)")

plt.grid()
plt.legend()
plt.show()
