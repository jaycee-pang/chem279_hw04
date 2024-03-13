CPP = g++
CPPFLAGS = -std=c++17 -I/opt/homebrew/opt/armadillo/include
ARMA_INCL = -L/opt/homebrew/opt/armadillo/lib -larmadillo -llapack -lblas
# OBJS = molecule.o AO.o
LIB = mylib.a 
# EXECS = try_h2
OBJS = AO_Zhe.o 
EXECS = try_h2

all: $(EXECS)
AO_Zhe.o: AO_Zhe.cpp AO_Zhe.h
	$(CPP) $(CPPFLAGS) -c AO_Zhe.cpp 
# molecule.o: molecule.cpp molecule.h AO.h
# 	$(CPP) $(CPPFLAGS) -c molecule.cpp 

# AO.o: AO.cpp AO.h molecule.h
# 	$(CPP) $(CPPFLAGS) -c AO.cpp

$(LIB): $(OBJS)
	ar rcs $(LIB) $(OBJS) 

# try_h2: try_h2.cpp $(OBJS) $(LIB)
# 	$(CPP) $(CPPFLAGS) -o try_h2 try_h2.cpp $(OBJS) $(LIB) $(ARMA_INCL)

clean:
	rm -f $(OBJS) $(EXECS) $(LIB)