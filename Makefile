TARGET = simctm
OBJS_DIR = build

PARDISO = -L/home/lunzhy/PhD_Study/numerical/pardiso -lpardiso500-GNU481-X86-64
LAPACK_BLAS = -L/usr/lib64/ -llapack -lblas
OMP = -L/usr/local/lib64 -lgomp
OTHERS_PARDISO = -lgfortran -fopenmp -lpthread -lm


.PHONY: default all clean

default: $(TARGET)
all: default

OBJECTS = $(patsubst %.cpp, $(OBJS_DIR)/%.o, $(wildcard *.cpp))
HEADERS = $(wildcard *.h)

$(OBJS_DIR)/%.o: %.cpp $(HEADERS)
	g++ -std=c++11 -c $< -o $@

.PRECIOUS: $(TARGET) $(OBJECTS)
$(TARGET): $(OBJECTS)
	g++ $(OBJECTS) $(PARDISO) $(LAPACK_BLAS) $(GOMP) $(OTHERS_PARDISO) -o $@

clean:
	-rm -f $(OBJS_DIR)/*.o
	-rm -f $(TARGET)


#INCLUDE = -I /usr/local
#g++ -std=c++11 -c *.cpp  $(INCLUDE)
#g++ -o simctm *.o
