rs(n) = sprintf("d_rs.dat%03d",n)
dtrue(n) = sprintf("d_true.dat")
ddots(n) = sprintf("d_dots.dat%03d",n)
splot rs(i) with lines t "Cokriging", dtrue(i) with lines t "Analytical", ddots(i) with points pt 7 t "Snapshots"
