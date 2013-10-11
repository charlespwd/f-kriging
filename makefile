default: 
	cd src; make

sensitivity:
	cd src; make sensitivity

clear :
	cd src; make clear
	cd data; rm -f *.dat*

