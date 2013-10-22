default: 
	cd src; make

sequential:
	cd src; make sensitivity

clear :
	cd src; make clear
	cd data; rm -f *.dat*

