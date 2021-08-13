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
# ALTERNATIVA 1: Tabelas em arquivos separados (2 colunas)
#
set ylabel  "Tempo (s)"
set title   "MFLOP/s"
set terminal qt 0 title "MFLOP/s"
plot '< sort -nk1 flops_ajc.csv' using 1:2 title "Tempo x Tamanho" with linespoints

# Gerando figura PNG
set terminal png
set output "flops_ajc.png"
plot '< sort -nk1 flops_ajc.csv' using 1:2 title "Tempo x Tamanho" with linespoints


replot
unset output
