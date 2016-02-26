set term postscript eps enhanced color
set output "box.eps"
splot 'box.dat' using 1:2:3 title 'Grafico'