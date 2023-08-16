from scitools.std import *

version = '03.02.2012'

# Constants
k = 10.0   
m = 0.1
b = 2.5
T = 5.0
dt = 0.01
omega = sqrt(k/m)
gamma = 0.5*(b/m)

# Numerical initialization 
N = T/dt
x = zeros(N)
v = zeros(N)
t = zeros(N)

# Initial conditions
x[0] = 1.0
v[0] = 0.1
t[0] = 0.0

# Differential equation
def diffEq(xNow,vNow,tNow):
	aNow = -(k/m)*xNow-(b/m)*vNow 
	return aNow
diffligning = 'a = -(k/m)*x-(b/m)*v'

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
	x[i+1],v[i+1] = rk4(x[i],v[i],t[i])
	t[i+1] = t[i] + dt

# Analytical solutions
if gamma>omega:
	A1 = x[0]/2 + (v[0]+gamma*x[0])/(2*sqrt(gamma**2-omega**2))
	A2 = x[0]/2 - (v[0]+gamma*x[0])/(2*sqrt(gamma**2-omega**2))
	x_exact = A1*exp((-gamma+sqrt(gamma**2-omega**2))*t)+A2*exp((-gamma-sqrt(gamma**2-omega**2))*t)
elif gamma==omega:
	A = x[0]
	B = v[0] + gamma
	x_exact = A*exp(-gamma*t)+B*t*exp(-gamma*t)
else:
	phi = arctan(-(v[0]+gamma)/sqrt(omega**2-gamma**2))
	A = x[0]/cos(phi)
	x_exact = exp(-gamma*t)*A*cos(sqrt(omega**2-gamma**2)*t+phi)

# Plot
figure(1)
plot(t, x, t[0:-1:10], x_exact[0:-1:10], 'x',
	xlabel='Time [s]', ylabel='Position [m]',
	legend=('Numerical','Analytical'),
	title=('Position vs time; b=%.2f Ns/m' % b),
	hardcopy=('o3c4x%.2f.png') % b)
figure(2)
plot(t, x_exact-x,
	xlabel='Time [s]', ylabel='Position [m]',
	legend=('Differance: x-x_exact'),
	title=('Position vs time; b=%.2f Ns/m' % b),
	hardcopy=('o3c4diff%.2f.png') % b)

# Documentation
outfile = open('oblig3c4_1.txt', 'w')
outfile.write('Versjon: %s \n' % version)
outfile.write('Fjarstivhet: k=%.2f N/m \n' % k)
outfile.write('Masse: m=%.2f kg \n' % m)
outfile.write('Friksjonsparameter: b=%.2f Ns/m \n' % b)
outfile.write('Tid: T=%.2f s \n' % T)
outfile.write('Tidssteg: dt=%.2f s \n' % k)
outfile.write('Omega: %.2f s^-1 \n' %  omega)
outfile.write('Gamma: %.2f s^-1 \n' %  gamma)
outfile.write('Maksimal hastighet: v_max=%.2f m/s \n' % max(v))
outfile.write('Diffligning: %s \n' %  diffligning)
outfile.write('Beregnede data: t(i)[s]   x(i)[m]   v(i)[m/s] \n')
for i in range(len(t)-1):
	outfile.write('                %.2f       %.2f       %.2f \n' % (t[i],x[i],v[i]))
outfile.close()

