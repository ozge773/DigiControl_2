clear all 
close all
clc
s = tf('s')
Ts = 0.2
G = 1 / (s^2 + 0.4*s + 0.04)
G_ZOH = 1 / (1+ s*Ts/2)

Kc = 1 %for every |Kc|
C_SS = Kc /s

L1 = G * C_SS
L11 = L1 * G_ZOH 
%zeta = 0.56
Tp = mag2db(1.07)
Sp = mag2db(1.39)
M_T_LF = mag2db(0.008/ 0.1)
figure
nichols(L11, 'b')
hold on
T_grid(Tp)
T_grid(M_T_LF)
S_grid(Sp)
%axis([-360 -90 -150 150])
wc_des = 1  %1.95 / 3.55
w_norm =4 % one zero is +70 phase
wz = wc_des / w_norm 
C_Z = (1 + s/wz)^2
L2 = L11* C_Z
nichols(L2, 'r')
K = 10 ^(-27/20)
L3 = L2 * K
nichols(L3, 'g')

%one colsure pole
wp = 4
C_P = 1 / (1 + s/wp)
L4 = L3 *C_P
nichols(L4, 'm')

C0 = C_SS * C_Z * C_P * K
Cd = c2d(C0, Ts, 'Tustin')

%%
rho = 1
delta_a = 0
delta_y = 0 
delta_t = 0 

wt = 0
out = sim("my2sim.slx")

stepinfo(out.y.data, out.y.time, 1, 0, 'RiseTimeLimits', [0 1], 'SettlingTimeStreshold', 0.05)
figure
plot(out.e.time, out.e.data, LineWidth=0.5)
%%
rho = 0
delta_a = 0
delta_y = 0 
delta_t = 0.1

wt = 6
out = sim("my2sim.slx")
plot(out.y.time, out.y.data, LineWidth=0.5)
