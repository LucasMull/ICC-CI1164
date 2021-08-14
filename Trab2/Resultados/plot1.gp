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
set ylabel  "Memory Bandwidth[MBytes/s]"
set title   "Banda de Memória"
set terminal qt 0 title "Banda de Memória"
plot '< sort -nk1 mem_ajc.csv' using 1:2 title "MBytes/s x Tamanho" with linespoints

# Gerando figura PNG
set terminal png
set output "mem_ajc.png"
plot '< sort -nk1 mem_ajc.csv' using 1:2 title "MBytes/s x Tamanho" with linespoints


replot
unset output
