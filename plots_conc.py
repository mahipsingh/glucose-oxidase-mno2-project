from numpy import *
from pylab import *

name = 'NaCl_1'
data = loadtxt(name + '.txt', skiprows=2)
Concentration = data[:,0]
Absorbance = data[:,1]
figure(1)
plot(Concentration, Absorbance, label = name, markersize=2)
plot(Concentration, Absorbance, 'bo', markersize=6)
ylabel('Absorbance')
xlabel('Concentration (ppm)')
grid(which='major',linewidth=0.2,ls='-',alpha=0.8)
legend()
savefig(name + '.pdf')

name2 = 'KCl_1'
data2 = loadtxt(name2 + '.txt', skiprows=2)
Concentration2 = data2[:,0]
Conductivity2 = data2[:,1]
figure(2)
plot(Concentration2, Conductivity2, label = name2, markersize=2)
plot(Concentration2, Conductivity2, 'bo', markersize=6)
ylabel('Absorbance')
xlabel('Concentration (ppm)')
grid(which='major',linewidth=0.2,ls='-',alpha=0.8)
legend()
savefig(name2 + '.pdf')

figure(3)
plot(Concentration, Absorbance, 'b', label = name, markersize=2)
plot(Concentration, Absorbance, 'bo', markersize=6)
plot(Concentration2, Conductivity2, 'r', label = name2, markersize=2)
plot(Concentration2, Conductivity2, 'ro', markersize=6)
ylabel('Absorbance')
xlabel('Concentration (ppm)')
grid(which='major',linewidth=0.2,ls='-',alpha=0.8)
legend()
savefig('combined' + '.pdf')

show()
