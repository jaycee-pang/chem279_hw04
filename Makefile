CPP = g++
CPPFLAGS = -std=c++14 -I/opt/homebrew/opt/armadillo/include
ARMA_INCL = -L/opt/homebrew/opt/armadillo/lib -larmadillo -llapack -lblas
OBJS = molecule.o AO.o
LIB = mylib.a 
EXECS = try_h2
molecule.o: molecule.cpp molecule.h AO.h
	$(CPP) $(CPPFLAGS) -c molecule.cpp
AO.o: AO.cpp AO.h molecule.h
	$(CPP) $(CPPFLAGS) -c AO.cpp

$(LIB): $(OBJS)
	ar rcs $(LIB) $(OBJS) 
try_h2: $(LIB) try_h2.cpp $(OBJS)
	$(CPP) $(CPPFLAGS) -o try_h2 try_h2.cpp $(OBJS) $(LIB) $(ARMA_INCL)
all: $(LIB) $(EXECS)
clean:
	rm -f $(OBJS)
	rm -f $(EXECS)
	rm $(LIB)
