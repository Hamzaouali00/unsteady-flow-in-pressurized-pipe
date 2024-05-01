import pandas as pd
import numpy as np
#initialisation
H0 = 90  # initial charge in the tank
Q0 = 0.15 # initial flow rate
f = .05
K = .5 # local loss coefficient
a = 1250 # velocity
L = 3750 #length of the channel
D = .4 # diameter of the channel
A = (3.14*D**2)/4
tf = 5
q = 4
#---etape1---
dx = L/3
dt = dx/a
n = 3
#---etape2---
x = np.linspace(0,L,4)
Hbar = H0-(1+K+(f*x)/D)*((Q0**2)/(2*9.81*A**2))
Qbar = [Q0,Q0,Q0,Q0]
#---etape3---
Hp = Hbar
Qp = Qbar
R = f/(2*D*A)
Ca = 9.81*A/a
h_list = []
q_list = []
tt = 17
for t in range(tt):
    row_h = []
    row_q = []
    for i in range(4):
        if i == 3:
            Qa = Qp[i - 1]
            Ha = Hp[i - 1]
            if t <= 5:
                tau = 1-t/5
            else:
                tau = 0
            Cv = ((tau*Q0)**2)/(Ca*Hbar[i])
            Cp = Qa + 9.81*(A/a)*Ha-R*dt*Qa*abs(Qa)
            qp = .5*(-Cv+(Cv**2+q*Cp*Cv)**.5)
            hp = (Cp-qp)/(Ca)
            row_q.append(qp)
            row_h.append(hp)
        if i == 0:
            Qb = Qp[i + 1]
            Hb = Hp[i + 1]
            Cn = Qb - 9.81 * (A / a) * Hb - R * dt * Qb*abs(Qb)
            if Qp[i] > 0:
                k1 = Ca * (1 + K) / (2 * 9.81 * A ** 2)
            else:
                k1 = Ca*(1-K)/(2*9.81*A**2)
            qp = (-1+(1+4*k1*(Cn+Ca*H0))**.5)/(2*k1)
            if Qp[i] > 0:
                hp = H0 - (1 + K) * (qp ** 2) / (2 * 9.81 * A ** 2)
            else:
                hp = H0 - (1 - K) * (qp ** 2) / (2 * 9.81 * A ** 2)
            row_q.append(qp)
            row_h.append(hp)
        elif 0 < i < 3:
            Qa = Qp[i - 1]
            Ha = Hp[i - 1]
            Qb = Qp[i + 1]
            Hb = Hp[i + 1]
            Cn = Qb - 9.81 * (A / a) * Hb - R * dt * Qb * abs(Qb)
            Cp = Qa + 9.81 * (A / a) * Ha - R * dt * Qa * abs(Qa)
            qp = .5*(Cp+Cn)
            hp = (Cp-Cn)/(2*Ca)
            row_q.append(qp)
            row_h.append(hp)
    Qp = row_q
    q_list.append(Qp)
    Hp = row_h
    h_list.append(Hp)
for rr in range(tt):
    print("----------\nt: ", rr)
    df = pd.DataFrame({"X":x, "Hp": h_list[rr], "Qp":q_list[rr]})
    print(df)

