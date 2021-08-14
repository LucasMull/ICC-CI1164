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
set ylabel  "Miss Ratio"
set title   "Cache Miss Ratio L2"
set terminal qt 0 title "Cache Miss Ratio L2"
plot '< sort -nk1 cmiss_ajc.csv' using 1:2 title "Miss Ratio x Tamanho" with linespoints

# Gerando figura PNG
set terminal png
set output "cmiss_ajc.png"
plot '< sort -nk1 cmiss_ajc.csv' using 1:2 title "Miss Ratio x Tamanho" with linespoints


replot
unset output
