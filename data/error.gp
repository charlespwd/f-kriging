# call with funcname=
set term pngcairo enh color size 1280,800
if (! exists("imgpath")) imgpath='.'
if (! exists("order")) {
	emse(n) = sprintf('e%smse.dat',n)
	esens(n) = sprintf('e%ssensitivity.dat',n)
} else {
	emsea(n) = sprintf('e%smse%darthur.dat',n,order)
	esensa(n) = sprintf('e%ssensitivity%darthur.dat',n,order)
	emsem(n) = sprintf('e%smse%dmartin.dat',n,order)
	esensm(n) = sprintf('e%ssensitivity%dmartin.dat',n,order)
	emsel(n) = sprintf('e%smse%dlophaven.dat',n,order)
	esensl(n) = sprintf('e%ssensitivity%dlophaven.dat',n,order)
}
set log y
set xlabel 'Number of Snapshots'
set ylabel 'L1 Norm of Error'
if (!exists("funcname")) {
	outfile = "error-.png"
	set output outfile
	plot "e.dat" with lp 
} else {
	outfile = sprintf('%s/error-%s.png',imgpath,funcname)
	if(exists("order")) {outfile = sprintf('%s/error-%s-%d.png',imgpath,funcname,order)}
	set output outfile
	plot emsea(funcname) with lp lw 2 ps 2 t 'Arthur-MSE', \
		emsel(funcname) with lp lw 2 ps 2 t 'Lophaven-MSE', \
		emsem(funcname) with lp lw 2 ps 2 t 'Martin-MSE', \
		esensa(funcname) with lp lw 2 ps 2 t 'Arthur-MSE+S', \
		esensl(funcname) with lp lw 2 ps 1 lt 7 t 'Lophaven-MSE+S', \
		esensm(funcname) with lp lw 2 ps 2  t 'Martin-MSE+S'
}


