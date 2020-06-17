# compile cpp code with HEP libraries
# christof sauer 2019

# compiler
CXX = g++
CXXFLAGS = -std=c++0x -Wall -g

SRCS =  main.cpp PatternEngine.cpp SuperPixelHit.cpp
OBJS = ${SRCS:.cpp=.o}

# HEP flags
ROOTFLAGS      := $(shell root-config --cflags)
ROOTLDFLAGS    := $(shell root-config --ldflags)
ROOTLIBS       := $(shell root-config --glibs)

# define variables
HEPLIBS += $(ROOTLIBS)

# Add HEP flags to common stuff
CXXFLAGS += $(ROOTFLAGS)

LIBS     += $(HEPLIBS)
LDFLAGS   = $(ROOTLDFLAGS)

MAIN = out

all: ${MAIN}
	@echo	"Start"

${MAIN}: ${OBJS}
	${CXX} ${CXXFLAGS} ${LDFLAGS} ${OBJS} -o ${MAIN} ${LIBS}

.cpp.o:
	${CXX} ${CXXFLAGS} -c $< -o $@

clean:
	${RM} ${PROGS} ${OBJS} *.o *~.

