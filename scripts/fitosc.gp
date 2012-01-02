
C0(t) = exp(-t/tau_G)
C(t)=1.0 / cos(phi) * cos(omega*t + phi) * exp(-t/tau_G)

file="<awk '{print NR-1, $1}' examples/DNA_PWELL/dna_pwell.fr-corr-2"
tau_G=40
phi = 1.0
omega = 0.025
fit [0:200] C(x) file via tau_G, omega, phi

plot [:200][:1.0] file w l, C(x), C0(x)