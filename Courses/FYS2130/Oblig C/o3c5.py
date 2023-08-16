from scitools.std import *
from o3c4 import *

version = '03.02.2012'

# Constants
k = 10.0   
m = 0.1
b = 2./3
T = 5.0
dt = 0.01
D = b/abs(max(v))
omega = sqrt(k/m)
gamma = omega/3

# Numerical initialization 
N = T/dt
x_notlinear = zeros(N)
v_notlinear = zeros(N)
t = zeros(N)

# Initial conditions
x_notlinear[0] = 1.0
v_notlinear[0] = 0.1
t[0] = 0.0

# Differential equation
def diffEq(xNow,vNow,tNow):
	aNow = -(k/m)*xNow-(b/m)*vNow-(D/m)*abs(vNow)*vNow 
	return aNow
diffligning = 'a = -(k/m)*x-(b/m)*v-(D/m)|v|*v'

# Runge-Kutta 4th order method
def rk4(xStart,vStart,tStart):
	a1 = diffEq(xStart,vStart,tStart)
	v1 = vStart
	xHalf1 = xStart + v1 * dt/2.0
	vHalf1 = vStart + a1 * dt/2.0

	a2 = diffEq(xHalf1,vHalf1,tStart+dt/2.0)
	v2 = vHalf1
	xHalf2 = xStart + v2 * dt/2.0
	vHalf2 = vStart + a2 * dt/2.0

	a3 = diffEq(xHalf2,vHalf2,tStart+dt/2.0)
	v3 = vHalf2
	xEnd = xStart + v3 * dt
	vEnd = vStart + a3 * dt

	a4 = diffEq(xEnd,vEnd,tStart + dt)
	v4 = vEnd
	aMiddle = 1.0/6.0 * (a1 + 2*a2 + 2*a3 + a4)
	vMiddle = 1.0/6.0 * (v1 + 2*v2 + 2*v3 + v4)

	xEnd = xStart + vMiddle * dt
	vEnd = vStart + aMiddle * dt
	
	return xEnd, vEnd

# Implementing our differential equation into the Runge-Kutta function
for i in range(len(t)-1):
	x_notlinear[i+1],v_notlinear[i+1] = rk4(x_notlinear[i],v_notlinear[i],t[i])
	t[i+1] = t[i] + dt

# Plot
figure(1)
plot(t, x, t, x_notlinear, '--',
	xlabel='Time [s]', ylabel='Position [m]',
	legend=('Linear','Not linear'),
	title=('Position vs time; D=%.2f Ns/m' % D),
	hardcopy=('o3c5x%.2f.png') % D)
figure(2)
plot(x, v, x_notlinear, v_notlinear, '--',
	xlabel='Position [m]', ylabel='Velocity [m/s]',
	title=('Velocity vs time; D=%.2f Ns/m' % D),
	legend=('Linear','Not linear')
	hardcopy=('o3c5v%.2f.png') % D)
raw_input('Press Enter: ')

# Documentation
outfile = open('oblig3c5_1.txt', 'w')
outfile.write('Versjon: %s \n' % version)
outfile.write('Fjarstivhet: k=%.2f N/m \n' % k)
outfile.write('Masse: m=%.2f kg \n' % m)
outfile.write('Friksjonsparameter: b=%.2f Ns/m \n' % b)
outfile.write('Friksjonsparameter: D=%.2f  Ns^2/m^2 \n' % D)
outfile.write('Tid: T=%.2f s \n' % T)
outfile.write('Tidssteg: dt=%.2f s \n' % k)
outfile.write('Omega: %.2f s^-1 \n' %  omega)
outfile.write('Gamma: %.2f s^-1 \n' %  gamma)
outfile.write('Maksimal hastighet: v_max=%.2f m/s \n' % max(v))
outfile.write('Diffligning: %s \n' %  diffligning)
outfile.write('Beregnede data: t(i)[s]   x(i)[m]   v(i)[m/s] \n')
for i in range(len(t)-1):
	outfile.write('                %.2f       %.2f       %.2f \n' % (t[i],x_notlinear[i],v_notlinear[i]))
outfile.close()
