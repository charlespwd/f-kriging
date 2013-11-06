set term pngcairo enh color size 1280, 960
if (! exists("funcname")) funcname='img'
if (! exists("imgpath")) imgpath='.'
if (! exists("nmax")) nmax = 0
rs(n) = sprintf("d_rs.dat%03d",n)
dtrue(n) = sprintf("d_true.dat000")
ddots(n) = sprintf("d_dots.dat%03d",n)
dstar(n) = sprintf("d_min.dat%03d",n)
do for [i=0:nmax] {
  	outfile = sprintf('%s/%s-ns%d.png',imgpath,funcname,i)
  	set output outfile
	if (exists("optimize")) {
		splot rs(i) with lines t "Cokriging", dtrue(i) with lines t "Analytical", \
		ddots(i) with points pt 7 t "Snapshots", dstar(i) with points pt 3 ps 3 t "ls-min"
	} else {
		splot rs(i) with lines t "Cokriging", dtrue(i) with lines t "Analytical", ddots(i) with points pt 7 t "Snapshots"
	}
}


