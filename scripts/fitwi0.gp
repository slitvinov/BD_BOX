C0(t) = exp(-t/tau)

fileX="<awk '{print NR-1, $1}' ~/Dropbox/dna_pwell.vel_grad0fr-corr-X"
fileY="<awk '{print NR-1, $1}' ~/Dropbox/dna_pwell.vel_grad0fr-corr-Y"

fit [0:maxt] C0(x) fileX via tau

maxt=300
plot [:2.0*maxt][:1.0] \
     fileX w p, \
     fileY w p, \
     C0(x)