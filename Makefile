CC = g++ -fdiagnostics-color=always
FLAGS = -std=c++17 -g -O3 -Wall -I.
LIBS = -pthread

BUILD = build
OBJ = build/obj

GRAPH = graph
MATCH = matching
UTILS = utils
BUILD_TOOLS = build/tools

all : dir $(BUILD)/csm

dir: $(OBJ)

$(OBJ) :
	mkdir -p $(OBJ)

#################### start ####################

$(BUILD)/csm: $(OBJ)/main.o \
		$(OBJ)/matching.o \
		$(OBJ)/sj_tree.o $(OBJ)/graphflow.o \
		$(OBJ)/turboflux.o $(OBJ)/symbi.o \
		$(OBJ)/iedyn.o \
		$(OBJ)/graph.o $(OBJ)/induced_graph.o \
		$(OBJ)/globals.o
	$(CC) $(FLAGS) $(OBJ)/main.o \
		$(OBJ)/matching.o \
		$(OBJ)/sj_tree.o $(OBJ)/graphflow.o \
		$(OBJ)/turboflux.o $(OBJ)/symbi.o \
		$(OBJ)/iedyn.o \
		$(OBJ)/graph.o $(OBJ)/induced_graph.o \
		$(OBJ)/globals.o \
		-o $(BUILD)/csm $(LIBS)

$(OBJ)/main.o: $(MATCH)/main.cpp \
		$(UTILS)/CLI11.hpp \
		$(UTILS)/globals.h $(UTILS)/types.h \
		$(GRAPH)/graph.h \
		$(MATCH)/sj_tree.h $(MATCH)/graphflow.h \
		$(MATCH)/turboflux.h $(MATCH)/symbi.h \
		$(MATCH)/iedyn.h
	$(CC) -c $(FLAGS) $(MATCH)/main.cpp -o $(OBJ)/main.o

#################### matching ####################

$(OBJ)/iedyn.o: $(MATCH)/iedyn.cpp \
		$(UTILS)/types.h $(UTILS)/globals.h \
		$(GRAPH)/graph.h \
		$(MATCH)/matching.h \
		$(MATCH)/iedyn.h
	$(CC) -c $(FLAGS) $(MATCH)/iedyn.cpp \
		-o $(OBJ)/iedyn.o

$(OBJ)/symbi.o: $(MATCH)/symbi.cpp \
		$(UTILS)/types.h $(UTILS)/globals.h \
		$(GRAPH)/graph.h \
		$(MATCH)/matching.h \
		$(MATCH)/symbi.h
	$(CC) -c $(FLAGS) $(MATCH)/symbi.cpp \
		-o $(OBJ)/symbi.o

$(OBJ)/turboflux.o: $(MATCH)/turboflux.cpp \
		$(UTILS)/types.h $(UTILS)/globals.h \
		$(GRAPH)/graph.h \
		$(MATCH)/matching.h \
		$(MATCH)/turboflux.h
	$(CC) -c $(FLAGS) $(MATCH)/turboflux.cpp \
		-o $(OBJ)/turboflux.o

$(OBJ)/graphflow.o: $(MATCH)/graphflow.cpp \
		$(UTILS)/types.h $(UTILS)/utils.h \
		$(UTILS)/globals.h \
		$(GRAPH)/graph.h \
		$(MATCH)/matching.h \
		$(MATCH)/graphflow.h
	$(CC) -c $(FLAGS) $(MATCH)/graphflow.cpp \
		-o $(OBJ)/graphflow.o

$(OBJ)/sj_tree.o: $(MATCH)/sj_tree.cpp \
		$(UTILS)/types.h $(UTILS)/globals.h \
		$(GRAPH)/graph.cpp $(GRAPH)/induced_graph.h \
		$(MATCH)/matching.h \
		$(MATCH)/sj_tree.h
	$(CC) -c $(FLAGS) $(MATCH)/sj_tree.cpp \
		-o $(OBJ)/sj_tree.o

$(OBJ)/matching.o: $(MATCH)/matching.cpp \
		$(UTILS)/types.h \
		$(GRAPH)/graph.h \
		$(MATCH)/matching.h
	$(CC) -c $(FLAGS) $(MATCH)/matching.cpp \
		-o $(OBJ)/matching.o

#################### graph ####################

$(OBJ)/graph.o: $(GRAPH)/graph.cpp \
		$(UTILS)/types.h $(UTILS)/utils.h \
		$(GRAPH)/graph.h
	$(CC) -c $(FLAGS) $(GRAPH)/graph.cpp \
		-o $(OBJ)/graph.o

$(OBJ)/induced_graph.o: $(GRAPH)/induced_graph.cpp \
		$(UTILS)/types.h \
		$(GRAPH)/induced_graph.h $(GRAPH)/graph.h
	$(CC) -c $(FLAGS) $(GRAPH)/induced_graph.cpp \
		-o $(OBJ)/induced_graph.o

#################### utils ####################

$(OBJ)/globals.o: $(UTILS)/globals.cpp $(UTILS)/globals.h
	$(CC) -c $(FLAGS) $(UTILS)/globals.cpp \
		-o $(OBJ)/globals.o

#################### end ####################

.PHONY: clean

clean: 
	rm -r ${BUILD}
