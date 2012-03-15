#set size 0.5, 0.5
set xlabel "k_{theta}"
set ylabel "corr. of bond orientation"
set log
set key bottom
plot \
     [0.001:10][0.5:200] \
     "s-vs-th.dat" w lp ps 4 t "", \
     49.0*x t "k^1", \
     54.0*x**1.1 t "k^{1.1}", \
     59.0 t ""

# k_theta for persistance length 1.0     
print 59.0/49.0
call "../../scripts/saver.gp" "s-vs-th"