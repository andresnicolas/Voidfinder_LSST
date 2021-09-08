#__________________________________________________________
#                                                  Compiler
                                                  
CC = g++

#__________________________________________________________
#                                                    voro++

PIXLIB = -L${HEALPIX}/lib -lhealpix_cxx 
PIXINC = -I${HEALPIX}/include/healpix_cxx 

#__________________________________________________________
#                                                    OpenMP

OMP = -fopenmp

#__________________________________________________________
#

EXEC = main.x

OBJS = tools.o io.o main.o

INCL = Makefile proto.h global.h

CFLAGS += $(OMP) $(PIXINC) 
   	
LIBS = $(PIXLIB) -lm 

$(EXEC): $(OBJS) 
	$(CC) $(OMP) $(OBJS) $(LIBS) -o $(EXEC)  

$(OBJS): $(INCL) 

clean:
	rm -f $(OBJS) $(EXEC) *~ 
