#!/bin/bash
PATH=.:/bin:/usr/bin


mkconf -n 1024 test.conf
md2 -T 2.0 -d 0.8 -t 1000 -a 1024 -b 0 -z 4.0 -r 1.5 -R 1.5 \
    -o test.data -c test.conf -C final.conf -H 0.0005 -Q 0.5 \
    -n 10 -l 10 -u
pltfile=$(mktemp /tmp/$0.XXXXXX)
cat > $pltfile <<EOF
set terminal pdfcairo
set output 'test.pdf'
plot 'test.data' using 1:2 title 'Kinetic energy' with lines, 'test.data' using 1:3 title 'Potential energy' with lines, 'test.data' using 1:6 title 'Total energy' with lines
plot 'test.data' using 1:5 title 'eta' with lines
EOF
gnuplot $pltfile
rm -f $pltfile
