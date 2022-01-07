#%% Hic Calculations
#import control
import numpy as np
import math as mt
from matplotlib import pyplot as plt
import control

from IPython import get_ipython
get_ipython().magic('reset -sf')
#%% Parameters
m = 5 # kg
k = 500.0 * 10**3 #N/m
v0 = 9 #m/s
wn = mt.sqrt(k/m) #rad/s
fn = wn / (2*np.pi) # Hz
t_del_max = 0.25/fn # s
steps = 500
time_step = t_del_max / (steps/2)
g = 9.81 #m/s^2

#%% Analytical Solution
t = np.linspace(0, 2*t_del_max, steps) #s
delta = v0/wn * np.sin(wn*t) #m
v_a = v0 * np.cos(wn*t) #m/s
a = -k * delta / m #m/s^2
ag_a = -a/g
F = k * delta #N 

#%% Euler Method
t = np.linspace(0, 2*t_del_max, steps) #s

#initial time step
v_e = [v0]
delta = [0]
F = [0]

#second time step
delta.append(1/2 * v_e[0] * time_step + delta[0])
F.append(k * delta[1])
v_e.append(v_e[0] - 1/(2 * m) * (F[0] + F[1]) * time_step)

#rest of time steps
for i in range(2, steps):
    delta.append(1/2 * (v_e[i-1] + v_e[i-2]) * time_step + delta[i-1])
    F.append(k * delta[i])
    v_e.append(v_e[i-1] - 1/(2 * m) * (F[i] + F[i-1]) * time_step)
F = np.asarray(F)
a = F/m
ag_e = a/g 

#%% Numerical vs Analytical Plots 
#Decleration vs Time 
plt.figure(0) 
plt.plot(t, ag_a, t, ag_e)
plt.title("Accel, Analytical vs Euler")
plt.xlabel("Time (s)")
plt.ylabel("Acceleration") 
plt.legend(('ag_a', 'ag_e'))
 
#Velocity vs Crush depth
plt.figure(1)
plt.plot(delta, v_a, delta, v_e)   
plt.title("Veloctiy, Analytical vs Euler")
plt.xlabel("Crush depth")
plt.ylabel("Velocity") 
plt.legend(('v_a', 'v_e'))

#Force vs Time
plt.figure(2)
plt.plot(t, F)   
plt.title("Force vs Time")
plt.xlabel("Time (s)")
plt.ylabel("Force") 
impulse = np.trapz(F,t)
#plt.legend(('v_a', 'v_e'))

#Force vs Delta
plt.figure(3)
plt.plot(delta, F)   
plt.title("Force vs Crush Depth")
plt.xlabel("Crush Depth")
plt.ylabel("Force") 
work = np.trapz(F,delta)
#%% HIC Calculations
a_avg = np.zeros([steps, steps]).tolist()
hic = np.zeros([steps, steps]).tolist()
for i in range(steps-1):
    for j in range(i+1,steps):
        a_avg[i][j] = (a_avg[i][j-1]*(t[j-1]-t[i]) + 0.5*(ag_e[j]+ag_e[j-1])*time_step) / (t[j]-t[i])
        hic[i][j] = a_avg[i][j]**2.5 * (t[j] - t[i]) #add up the average values between the two timesteps       

#a_avg = np.transpose(a_avg)   

#max HIC value
max_hic_map = map(max, hic)    
max_hic_list = list(max_hic_map) 
max_hic = max(max_hic_list)
max_hic_column = max_hic_list.index(max_hic)

a_avg = np.asarray(a_avg).T
hic = np.asarray(hic).T

#%% HIC vs ag Plots
fig, host = plt.subplots()
par1 = host.twinx()

p1, = host.plot(t, ag_e, "b-", label="ag")
p2, = par1.plot(t, hic[:, max_hic_column], "r-", label="HIC")

host.set_xlabel("Time (s)")
host.set_ylabel("ag")
par1.set_ylabel("HIC")

lines = [p1, p2]

host.legend(lines, [l.get_label() for l in lines])
plt.title('HIC vs ag')
plt.show()