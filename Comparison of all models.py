
# coding: utf-8

#%%
import numpy as np
from astropy import units as u
from astropy.table import Table
#from astropy.table.column import Column
import scipy.special as sps
import matplotlib.pyplot as plt
#%%

H0 = 95*u.pc 
R0_adams = 2.9e3*u.pc
R_sn = 0.03*u.yr**-1
R_sun = 8700 * u.pc
z = 24*u.pc

#%%

alpha = 1.09
beta =  3.87
R0_green = R_sun/beta
#%%
#comparision of all models within this cell


a = R_sn/(4*np.pi*sps.gamma(alpha+2)*R0_green**2*H0)
b = (R_sun/R0_green)**alpha
c = np.exp(-R_sun/R0_adams)
d = np.exp(-np.abs(z)/H0)


Hg, Rg = 100*u.pc, 15000*u.pc
q_uniform = R_sn/(2*np.pi*Hg*Rg**2)
q_uni = q_uniform
q_adams = (R_sn)/(4*np.pi*H0*R0_adams**2)*c*d
q_ad = q_adams
c = np.exp(-R_sun/R0_green)
q_green = a*b*c*d
q_gr = q_green
q_s = [q_uniform, q_adams, q_green]    

#returns the rate in units of SN/Gyr
#enter d in units of pc 
def rate(q, d):
    return q*(4/3)*np.pi*(d)**3

#adams dist 175.131 pc
#enter rate in units of yr^-1
sn_rate = 2.6e-6*u.yr**-1
def mean_dist(q, rate):
    return ((3*rate)/(4*np.pi*q))**(1/3)
 
#adams: 2.066e9 yrs
#enter d in units of pc
def mean_time(q, d):
    return 1/(q*(4/3)*np.pi*(d)**3)

#gives the rates of SN/yr for our favorite distances as a list depending on the input q
#unit in SN/years, unit takes Myr, Gyr 
distances = [10, 50, 100, 150]*u.pc 
def SN_per_yr(q, unit = u.yr):
    rates = []
    for dist in distances:
        r = rate(q,dist)
        #"%.4g" % 
        rates.append(r.value)
    rates = rates*u.yr**-1
    if unit == u.Myr:
        rates = rates.to(u.Myr**-1)
    if unit == u.Gyr:
        rates = rates.to(u.Gyr**-1)
    return rates

def all_SN_information():
    rates=[]
    avg_distance=[]
    avg_time = []
    for q in q_s:
        d = mean_dist(q, sn_rate)
        t = mean_time(q, 10*u.pc)
        avg_distance.append(d.value)
        avg_time.append(t.value)
    for d in distances:
        r = rate(q, 10*u.pc)
        rates.append(r.value)
    return (rates, avg_distance, avg_time)    

green = SN_per_yr(q_green, u.Gyr)
adams = SN_per_yr(q_adams, u.Gyr)
uniform = SN_per_yr(q_uniform, u.Gyr)

all_rates = Table([distances, uniform, adams, green], names=('Distances', 'Uniform', 'Adams', 'Green'), meta={'name': 'SN rates'})

all_mean_dist = Table

#%%
#Killer GRBs
#short bursts from neutron stars merging, long from supernovae
d = 2000*u.pc
#factor of 10**-5 comes from rate of grb/rate of sn find that number 
q_GRB = q_adams*10**-5
rate_GRB = rate(q_GRB, 2000*u.pc)
rate_GRB = rate(q_GRB, 200*u.pc)
'''
#%%

radius = np.linspace(0.1, 15000, num = 100)*u.pc
yvals = np.linspace(0,10**-9, num = 100)
sun = []
for value in yvals:
    sun.append(R_sun.value)
    
a = R_sn/(4*np.pi*sps.gamma(alpha+2)*R0_green**2*H0)
b = (R_sun/R0_green)**alpha

d = np.exp(-np.abs(z)/H0)


Hg, Rg = 100*u.pc, 15000*u.pc


q_adams = []
q_uniform = []
q_green = []
for value in radius:
    q_uniform.append((R_sn/(2*np.pi*Hg*Rg**2)).value)
    q_adams.append(((R_sn)/(4*np.pi*H0*R0_adams**2)*(np.exp(-value/R0_adams))*d).value)
    q_green.append((a*b*(np.exp(-value/R0_green))*d).value)
plt.figure
plt.plot(radius, q_uniform, color = 'cornflowerblue', label = 'Uniform Model')
plt.plot(radius, q_green, color = 'seagreen', label = 'Green Model')
plt.plot(radius, q_adams, color = 'darkviolet', label = 'Adams Model')
plt.xlabel('Distance from Milky Way Center (R)[pc]')
plt.ylabel('Local Supernova Rate Density (q)')
plt.xscale('log')
plt.yscale('log')
plt.plot(sun, radius,'r--', label = 'Earth + Sun location')
plt.title("Comparison of All Supernova Models", fontsize = 15)
plt.legend()
plt.ylim(0, 10**-10.9)
plt.text(10**3.9, 10**-10.6, 'Preliminary')
'''
#%%

gamma_uniform = []
gamma_adams = []
gamma_green = []
for (a, b, c, d) in zip(radius, q_adams, q_uniform, q_green):
    gamma_adams.append(rate(b, a).value)
    gamma_uniform.append(rate(c, a).value)
    gamma_green.append(rate(d, a).value)


log_uniform = []
log_adams = []
log_green = []
log_radius = []
for (a, b, c, d) in zip(radius, gamma_adams, gamma_uniform, gamma_green):
    log_radius.append(np.abs(np.log(a.value)))
    log_adams.append(np.abs(np.log(b)))
    log_uniform.append(np.abs(np.log(c)))
    log_green.append(np.abs(np.log(d)))


plt.figure()
plt.plot(radius, gamma_uniform, color = 'cornflowerblue', label = 'Uniform Model')
plt.plot(radius, gamma_green, color = 'seagreen', label = 'Green Model')
plt.plot(radius, gamma_adams, color = 'darkviolet', label = 'Adams Model')
plt.xlabel('log(Radius [pc])')
plt.ylabel('log(Local Supernova Rate (Γ))')
plt.xscale('log')
plt.yscale('log')
#plt.plot(sun, fluence,'r--', label = 'Earth + Sun location')
plt.title("Comparison of All Supernova Models", fontsize = 15)
plt.legend()
#plt.ylim(0, 10**-10.9)
#plt.text(10**3.9, 10**-10.6, 'Preliminary')
#%%

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax2 = ax1.twiny()

J_kill = 1e5*(u.J/u.m**2)
def E(r_kill):
    4*np.pi*r_kill**2*J_kill

SN_kill = 10*u.pc
SGRB_kill = 200*u.pc
LGRB_kill = 2000*u.pc

r = []
rate = []
f = np.linspace(0, 1e10*(1.05e33), num = 100)
for value in f:
    r.append(np.sqrt(E_sn/4*np.pi*value))
    

SN = []
sgrb = []
lgrb = []
#calculate r from fluence array and use below
for value in r: 
    #f_sn.append(F(value, E_sn).value)
    #f_sgrb.append(F(value, E_sgrb).value)
    #f_lgrb.append(F(value, E_lgrb).value)
    SN.append(rate(q_adams, SN_kill).value)
    sgrb.append(Adams_rate(value, SGRB_kill).value)
    lgrb.append(Adams_rate(value, LGRB_kill).value)

#pylab.plot(4.84e-10,6.08E35, color = 'cornflowerblue', label = 'Supernovae')
ax1.plot(f, SN, color = 'cornflowerblue', label = 'Supernovae')
ax1.plot(f, sgrb, color = 'orange', label = 'SGRB')
ax1.plot(f, lgrb, color = 'red', label = 'LGRB')
ax1.axvline (x=9.523809523809e37,color='green', linestyle='dashed', label = 'Kill Fluence')

#plt.axhline(y=4.8403911e-10, color='cornflowerblue', linestyle='-')

#ax1.set_xlim(10e34,10e50)
#plt.xlim(10,10e50)

#ax1.set_ylim(10e-30,10)
ax1.set_xscale('log')
ax2.set_xscale('log')
ax1.set_xlabel ('Fluence (J/pc^2)')
ax1.set_ylabel ('Rate (yr)')
ax1.set_title("Fluence versus Rate (Adams Model)", y = 1.15, fontsize = 15)

ax1.legend()

new_ticks = []
for x in ax1.get_xticks():
    new_ticks.append(x*1e-33)
x = []
for value in new_ticks:
    x.append('%1.e'% value)
    
ax2.set_xticks(ax1.get_xticks())
ax2.set_xbound(ax1.get_xbound())
ax2.set_xticklabels(x)
ax2.set_xlabel("(J/m^2)")
plt.yscale('log')
pylab.show()





#%%
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax2 = ax1.twiny()
E_sn = 1.14e40 *u.J
#=10^47 erg
#same exercise as Supernovae for gamma ray bursts energy
    
#E_sn = E(10*u.pc).to(u.J)

r = []
rate = []
f = np.linspace(0, 1e10*(1.05e33), num = 100)
for value in f:
    r.append(np.sqrt(E_sn/4*np.pi*value))
    rate.append((4/3)*np.pi*q_adams*(E/(4*np.pi*value)))
    
SN = []
sgrb = []
lgrb = []
#calculate r from fluence array and use below
for value in r: 
    #f_sn.append(F(value, E_sn).value)
    #f_sgrb.append(F(value, E_sgrb).value)
    #f_lgrb.append(F(value, E_lgrb).value)
    SN.append((rate).value)
    sgrb.append((rate(value, R_sgrb).value)
    lgrb.append((rate(value, R_lgrb).value)

#pylab.plot(4.84e-10,6.08E35, color = 'cornflowerblue', label = 'Supernovae')
ax1.plot(f, SN, color = 'cornflowerblue', label = 'Supernovae')
ax1.plot(f, sgrb, color = 'orange', label = 'SGRB')
ax1.plot(f, lgrb, color = 'red', label = 'LGRB')
ax1.axvline (x=9.523809523809e37,color='green', linestyle='dashed', label = 'Kill Fluence')

#plt.axhline(y=4.8403911e-10, color='cornflowerblue', linestyle='-')

#ax1.set_xlim(10e34,10e50)
#plt.xlim(10,10e50)

#ax1.set_ylim(10e-30,10)
ax1.set_xscale('log')
ax2.set_xscale('log')
ax1.set_xlabel ('Fluence (J/pc^2)')
ax1.set_ylabel ('Rate (yr)')
ax1.set_title("Fluence versus Rate (Adams Model)", y = 1.15, fontsize = 15)

ax1.legend()

new_ticks = []
for x in ax1.get_xticks():
    new_ticks.append(x*1e-33)
x = []
for value in new_ticks:
    x.append('%1.e'% value)
    
ax2.set_xticks(ax1.get_xticks())
ax2.set_xbound(ax1.get_xbound())
ax2.set_xticklabels(x)
ax2.set_xlabel("(J/m^2)")
plt.yscale('log')
pylab.show()

#%%
E_sn = 5.46e32 *u.J

fluence = []
for value in radius:
    fluence.append(E_sn/(4*np.pi*value**2))

gamma_uniform = []
gamma_adams = []
gamma_green = []
for (a, b, c, d) in zip(radius, q_adams, q_uniform, q_green):
    gamma_adams.append(rate(b, a).value)
    gamma_uniform.append(rate(c, a).value)
    gamma_green.append(rate(d, a).value)


plt.figure()
plt.plot(log_radius, log_uniform, color = 'cornflowerblue', label = 'Uniform Model')
plt.plot(log_radius, log_green, color = 'seagreen', label = 'Green Model')
plt.plot(log_radius, log_adams, color = 'darkviolet', label = 'Adams Model')
plt.xlabel('log(Radius [pc])')
plt.ylabel('log(Local Supernova Rate (Γ))')
plt.xscale('log')
plt.yscale('log')
#plt.plot(sun, fluence,'r--', label = 'Earth + Sun location')
plt.title("Comparison of All Supernova Models", fontsize = 15)
plt.legend()
#plt.ylim(0, 10**-10.9)
#plt.text(10**3.9, 10**-10.6, 'Preliminary')





'''
#%%
#Create a plot of fluence vs rate like the Melott & Thomas paper
#comment on dust effects

#E_sn= 10*u.MeV 
#Extreme UV energy: E_sn= 35*u.eV
#Find the opacity of E UV 
#given the r_kill & f_kill find the E_sn

E_sn = (10**37)*u.J
E_sgrb = (10e43*u.J)
E_lgrb = (5e45*u.J)

fluence = []
for value in radius:
    fluence.append(E_sn/(4*np.pi*value**2))
    
    
def gamma_fluence(q, f):
    return q*(4/3)*np.pi*(E_sn/(4*np.pi*f))**(3/2)

gamma_uniform = []
gamma_adams = []
gamma_green = []
for (a, b, c, d) in zip(fluence, q_adams, q_uniform, q_green):
    gamma_adams.append(gamma_fluence(b, a).value)
    gamma_uniform.append(gamma_fluence(c, a).value)
    gamma_green.append(gamma_fluence(d, a).value)
    
log_uniform = []
log_adams = []
log_green = []
log_fluence = []
for (a, b, c, d) in zip(fluence, gamma_adams, gamma_uniform, gamma_green):
    fluence.append(np.abs(np.log(a)))
    gamma_adams.append(np.abs(np.log(b)))
    gamma_uniform.append(np.abs(np.log(c)))
    gamma_green.append(np.abs(np.log(d)))


plt.figure()
plt.plot(fluence, gamma_uniform, color = 'cornflowerblue', label = 'Uniform Model')
plt.plot(fluence, gamma_green, color = 'seagreen', label = 'Green Model')
plt.plot(fluence, gamma_adams, color = 'darkviolet', label = 'Adams Model')
plt.xlabel('log(Fluence) ()[J/pc^2]')
plt.ylabel('log(Local Supernova Rate (Γ))')
plt.xscale('log')
plt.yscale('log')
#plt.plot(sun, fluence,'r--', label = 'Earth + Sun location')
plt.title("Comparison of All Supernova Models", fontsize = 15)
plt.legend()
#plt.ylim(0, 10**-10.9)
#plt.text(10**3.9, 10**-10.6, 'Preliminary')
    



'''
#%%
plt.figure()
explosions = ['ccSN', 'Short Hard GRB', 'Long Soft GRB']
years_ago = np.linspace(0, (500e6)*10**-9)
'''
Ordovician extinction: ~445 million years ago
Devonian extinction: ~375 million years ago
Permian extinction: ~250 million years ago
Triassic extinction: ~200 million years ago
Cretaceous extinction: ~65 million years ago
'''
ind = np.arange(4)
zeroes = []
ordovician = []
denovian = []
permian = []
triassic = []
cretaceous = []
for value in ind:
    zeroes.append(0)
    ordovician.append((445e6)*10**-9)
    denovian.append((375e6)*10**-9)
    permian.append((250e6)*10**-9)
    triassic.append((200e6)*10**-9)
    cretaceous.append((65e6)*10**-9)

plt.plot(ind, ordovician, 'b--', label = 'Ordovician Extinction')
plt.plot(ind, denovian, 'm--', label = 'Denovian Extinction')
plt.plot(ind, permian, 'c--', label = 'Permian Extinction')
plt.plot(ind, triassic, 'g--', label = 'Triassic Extinction')
plt.plot(ind, cretaceous, 'r--', label = 'Cretaceous Extinction')
'''
plt.text(0.25,(460e6)*10**-9, 'Ordovician Extinction')
plt.text(0.25,(390e6)*10**-9, 'Denovian Extinction')
plt.text(0.25,(265e6)*10**-9, 'Permian Extinction')
plt.text(0.25,(215e6)*10**-9, 'Triassic Extinction')
plt.text(0.25,(80e6)*10**-9, 'Cretaceous Extinction')
'''
plt.legend()
plt.xlim(0,3)
plt.yscale('log')

ind = np.arange(4)
plt.ylabel('Years Ago (Gyr)')
ind = [0.5, 1.5, 2.5]
plt.xticks(ind, (explosions))

t_ccSN = mean_time(q_ad, 150*u.pc)*10**-9
t_softGRB = mean_time(q_GRB, 1000*u.pc)*10**-9
t_hardGRB = mean_time(q_GRB, 500*u.pc)*10**-9


plt.text(0.25, t_ccSN.value+.0002, 'd = 150 pc')
plt.text(0.25, 1.25e-4, 'd' + 'ₖᵢₗₗ'+ ' = 10 pc')

plt.text(1.25, t_hardGRB.value-0.75, 'd = 500 pc')
plt.text(1.25, 1.25e-4, 'd' + 'ₖᵢₗₗ'+ ' = 200 pc')

plt.text(2.25, t_softGRB.value-0.1, 'd = 1000 pc')
plt.text(2.25, 1.25e-4, 'd' + 'ₖᵢₗₗ'+ ' = 2000 pc')

plt.text(2.65, 5.6, 'Preliminary')
plt.xlabel('Types of Explosions', labelpad = 20)
plt.title('Extinction Level Explosions vs Mass Extinctions', fontsize = 15)

plt.plot(0.5, t_ccSN, marker = 'o', markersize = 5, color='darkred')
plt.plot(1.5, t_hardGRB, marker = 'o', markersize = 5, color='darkred')
plt.plot(2.5, t_softGRB, marker = 'o', markersize = 5, color='darkred')