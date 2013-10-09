set term pngcairo enh color size 1280, 960
rs(n) = sprintf("d_rs.dat%d",n)
dtrue(n) = sprintf("d_true.dat%d",n)
ddots(n) = sprintf("d_dots.dat%d",n)
do for [i=1:6] {
  	outfile = sprintf('~/w-reports/tr/oc1/drag-ns%d.png',i)
  	set output outfile
	splot rs(i) with lines t "Cokriging", dtrue(i) with lines t "Analytical", ddots(i) with points pt 7 t "Snapshots"
}


