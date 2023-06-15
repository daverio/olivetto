PYTHON_MODULE := ../oliveto.so

SOURCE_PATH := /home/daverio/Documents/olivetto/cpp_module/src


CXX := g++
CXXFLAGS_COMPILE := -std=c++17 -O3 -fPIC
INCLUDES := -I$(SOURCE_PATH) -I/usr/include/python3.11/

CXXFLAGS_LINKAGE := -L /lib64 -shared -std=c++17 -O3
LIBS := -lpython3 -lboost_python3 -lboost_numpy3 -lgsl -lgslcblas -lm

OBJECTS := $(SOURCE_PATH)/oliveto.o \
		   $(SOURCE_PATH)/point.o \
		   $(SOURCE_PATH)/polygone_tests.o \
		   $(SOURCE_PATH)/wrapper.o


$(PYTHON_MODULE): $(OBJECTS)
	$(CXX) $(CXXFLAGS_LINKAGE) $(OBJECTS) -o $@ $(LIBS)

%.o: %.cpp
	$(CXX) $(INCLUDES) $(CXXFLAGS_COMPILE)  -c $< -o $@

clean:
	rm $(OBJECTS) $(PYTHON_MODULE)
