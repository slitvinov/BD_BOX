
C0(t) = exp(-t/tau_G)
C(t)=1.0 / cos(phi) * cos(omega*t + phi) * exp(-t/tau_G)

file="<awk '{print NR-1, $1}' examples/DNA_PWELL/dna_pwell.fr-corr"
tua_G=40
phi = 10
omega = 0.025
plot [0:350][-0.4:1.0] file w l, C(x)