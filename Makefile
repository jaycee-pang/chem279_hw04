CPP = g++
CPPFLAGS = -std=c++17 -I/opt/homebrew/opt/armadillo/include
ARMA_INCL = -L/opt/homebrew/opt/armadillo/lib -larmadillo -llapack -lblas
OBJS = molecule.o AO.o cndo.o
LIB = mylib.a 

# OBJS = AO_Zhe.o EH_Zhe.o util.o
EXECS = main explore

all: $(EXECS)
# util.o: util.cpp util.h
# 	$(CPP) $(CPPFLAGS) -c util.cpp 
# AO_Zhe.o: AO_Zhe.cpp AO_Zhe.h
# 	$(CPP) $(CPPFLAGS) -c AO_Zhe.cpp 
# EH_Zhe.o: EH_Zhe.cpp EH_Zhe.h
# 	$(CPP) $(CPPFLAGS) -c EH_Zhe.cpp 

molecule.o: molecule.cpp molecule.h AO.h
	$(CPP) $(CPPFLAGS) -c molecule.cpp 

AO.o: AO.cpp AO.h
	$(CPP) $(CPPFLAGS) -c AO.cpp 

cndo.o: cndo.cpp cndo.h
	$(CPP) $(CPPFLAGS) -c cndo.cpp 

$(LIB): $(OBJS)
	ar rcs $(LIB) $(OBJS) 

main: main.cpp $(OBJS) $(LIB)
	$(CPP) $(CPPFLAGS) -o main main.cpp $(OBJS) $(LIB) $(ARMA_INCL)
explore: explore.cpp $(OBJS) $(LIB)
	$(CPP) $(CPPFLAGS) -o explore explore.cpp $(OBJS) $(LIB) $(ARMA_INCL)

clean:
	rm -f $(OBJS) $(EXECS) $(LIB)