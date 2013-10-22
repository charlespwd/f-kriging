set term pngcairo enh color size 1280, 960
if (! exists("funcname")) funcname='img'
if (! exists("imgpath")) imgpath='.'
if (! exists("nmax")) nmax = 0
rs(n) = sprintf("d_rs.dat%03d",n)
dtrue(n) = sprintf("d_true.dat000")
ddots(n) = sprintf("d_dots.dat%03d",n)
do for [i=0:nmax] {
  	outfile = sprintf('%s/%s-ns%d.png',imgpath,funcname,i)
	if(exists("order")) { 
		outfile = sprintf('%s/%s-o%d-ns%d.png',imgpath,funcname,order,i)
	}
  	set output outfile
	splot rs(i) with lines t "Cokriging", dtrue(i) with lines t "Analytical", ddots(i) with points pt 7 t "Snapshots"
}


