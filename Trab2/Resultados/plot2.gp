#!/usr/bin/gnuplot -c
set grid
set style data point
set style function line
set style line 1 lc 3 pt 7 ps 0.3
set boxwidth 1
set xtics
set xrange ["0":]
set xlabel  "Tamanho (N)"

#
# ALTERNATIVA 2: Tabela com 3 colunas 
#
set ylabel  "Tempo (s)"
set title   "MFLOP/s"
set terminal qt 1 title "MFLOP/s"
plot '< sort -nk1 flops_triang.csv' using 1:2 title "com otimização" with linespoints, \
     '' using 1:3 title "sem otimização" with linespoints

# Gerando figura PNG
set terminal png
set output "flops_triang.png"
plot '< sort -nk1 flops_triang.csv' using 1:2 title "com otimização" with linespoints, \
     '' using 1:3 title "sem otimização" with linespoints

replot
unset output

