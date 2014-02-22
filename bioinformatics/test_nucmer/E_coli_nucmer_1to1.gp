set terminal png tiny size 800,800
set output "wga_output/E_coli_nucmer_1to1.png"
set size 1,1
set grid
unset key
set border 15
set tics scale 0
set xlabel "gi|15829254|ref|NC_002695.1|"
set ylabel "gi|26245917|ref|NC_004431.1|"
set format "%.0f"
set mouse format "%.0f"
set mouse mouseformat "[%.0f, %.0f]"
set mouse clipboardformat "[%.0f, %.0f]"
set xrange [1:5498450]
set yrange [1:5231428]
set style line 1  lt 1 lw 3 pt 6 ps 1
set style line 2  lt 3 lw 3 pt 6 ps 1
set style line 3  lt 2 lw 3 pt 6 ps 1
plot \
 "wga_output/E_coli_nucmer_1to1.fplot" title "FWD" w lp ls 1, \
 "wga_output/E_coli_nucmer_1to1.rplot" title "REV" w lp ls 2
