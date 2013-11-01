default: 
	cd src; make

sequential:
	cd src; make sensitivity

modular: 
	cd src; make modular

clear :
	cd src; make clear
	cd data; rm -f *.dat*
	rm *.dat

