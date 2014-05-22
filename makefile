backopa: background.o chemical.o example.o
	ifort -w -O3 background.o chemical.o example.o -o backopa

clean:
	find . -maxdepth 2 -name "*.o" -delete ; find . -maxdepth 1 -name "*.mod" -delete ;
	find . -maxdepth 2 -name "*~" -delete ; find . -maxdepth 2 -name "backopa" -delete
        
chemical.o: chemical.f90
	ifort -c -O3 chemical.f90
	
background.o: background.f90
	ifort -w -c -O3 background.f90

example.o: background.o chemical.o example.f90
	ifort -w -c -O3 example.f90
