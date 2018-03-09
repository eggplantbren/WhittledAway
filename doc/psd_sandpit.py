from numba import jit
from pylab import *

w0 = 0.6
Q = 2.0 #5.0
eta = sqrt(1 - 1/(4*Q**2))
A = 2.1 #3.5
L = 30.0

t = linspace(-300, 300, 10001)
C = A**2*exp(-w0*abs(t)/(2*Q))*(cos(eta*w0*abs(t)) + sin(eta*w0*abs(t))/(2*eta*Q))
C = A**2*exp(-abs(t)/L)

@jit
def fourier_transform(f, w):
    F = zeros(len(w))
    for i in range(len(w)):
        F[i] = real(trapz(f*exp(-1j*w[i]*t), x=t))
        print(i+1)
    return F

figure(1)
plot(t, C)
xlabel("$t$")
ylabel("Covariance")

w = linspace(-10, 10, 10001)
S0 = A**2/(w0*Q)
#S = 2*S0*w0**4/((w**2 - w0**2)**2 + w0**2*w**2/Q**2)
S = 2.0*A*A*L/(L*L*w*w + 1.0);

F = fourier_transform(C, w)
print("Ratio = {r}".format(r=(F/S)[len(F)//2]))

figure(2)
plot(w, F)
plot(w, S)
xlabel("$\\omega$")
ylabel("PSD")
show()

