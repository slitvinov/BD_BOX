
C0(t) = exp(-t/tau_G0)
C(t)=1.0 / cos(phi) * cos(omega*t + phi) * exp(-t/tau_G)

fileX="<awk '{print NR-1, $1}' ~/Dropbox/dna_pwell.vel_grad2e-6fr-corr-X"
fileY="<awk '{print NR-1, $1}' ~/Dropbox/dna_pwell.vel_grad2e-6fr-corr-Y"
tau_G=40
phi = 1.0
omega = 0.025

maxt=300
fit [0:maxt] C(x) fileX via tau_G, omega, phi
fit [0:maxt] C0(x) fileX via tau_G0

plot [:2.0*maxt][:1.0] \
     fileX w p, \
     fileY w p, \
     C(x), C0(x)
