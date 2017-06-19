#ROOTSYS	     := /opt/products/root/5.20.00/
ROOTCONFIG   := $(ROOTSYS)/bin/root-config
ROOTCFLAGS   := $(shell $(ROOTCONFIG) --cflags)
ROOTLDFLAGS  := $(shell $(ROOTCONFIG) --ldflags)
ROOTLIBS     := $(shell $(ROOTCONFIG) --libs)
#ROOTLIBS  = $(shell root-config --libs) -lTreePlayer -lfortranstat

ROOTGLIBS    := $(shell $(ROOTCONFIG) --glibs)
ROOTCINT=$(ROOTSYS)/bin/rootcint

INCLUDEDIR=-I.

# CXXFLAGS=-isystem $(shell $(ROOTCONFIG) --incdir) -I$(ROOTSYS)/htmldoc $(INCLUDEDIR) -O2 -g $(ROOTCFLAGS)
CXXFLAGS= -I$(ROOTSYS)/htmldoc -fPIC $(INCLUDEDIR) -O2 -g $(ROOTCFLAGS)
LDFLAGS=$(ROOTLDFLAGS) -L.

#EXTRA_LIBDIRS = -L/user/helebran/3.1.15/lqlimit/h1method/weights/ -lweights
# DEFINES = -Df2cFortran -D__NOBOS -DDEBUG=0 -DINLINE=inline


LIB=stslimit
LB=lib$(LIB).so
DICTLB=lib$(LIB)Dict.so

# OBJC=
# OBJC += TSysLimit.o
# OBJC += TSysLimitChannel.o
# OBJC += TSysLimitResult.o
# OBJC += TSysLimitResultPBock.o
# OBJC += TSysLimitResultLogL.o
# OBJC += TSysLimitScan.o
#OBJC += weight.o
#OBJC += lqme.o

# CINTOBJC=
# CINTOBJC += TSysLimitDict.o
# CINTOBJC += TSysLimitChannelDict.o
# CINTOBJC += TSysLimitResultDict.o
# CINTOBJC += TSysLimitResultPBockDict.o
# CINTOBJC += TSysLimitResultLogLDict.o
# CINTOBJC += TSysLimitScanDict.o

#include ./H1Basic.mk
ROOTLIBS     := $(shell $(ROOTCONFIG) --libs)

#bin:   testlimit compute_syst_impact print_event_table
bin:   testlimit 
clean:
	rm -f $(DICTLB) $(LB) *.o *Dict.C *Dict.h

depend:
	makedepend $(INCLUDEDIR) *.C

%Dict.C: %.C %.h
	rm -f $@ $(@:.C=.h)
	$(ROOTCINT) $@ -c $<

.C.o:
	$(CXX) $(CXXFLAGS) -shared $< -o $@

$(LB): $(OBJC)
	ld --export-dynamic -shared -rpath . $(OBJC) -o $(LB)

$(DICTLB): $(CINTOBJC)
	ld --export-dynamic -shared -rpath . $(CINTOBJC) -o $(DICTLB)

docu: docu.C $(LB) $(DICTLB) Makefile
	$(CXX) $(CXXFLAGS) $< -o  $@ $(LDFLAGS) -l$(LIB)Dict -l$(LIB) \
	$(ROOTLIBS) -lHtml
	./docu

# testlimit: testlimit.C $(LB) $(DICTLB) Makefile
# 	$(CXX) $(CXXFLAGS) $< -o  $@ $(LDFLAGS) -l$(LIB)Dict -l$(LIB) \
# 	$(ROOTLIBS) -lrt

testlimit: testlimit.C Makefile
	$(CXX) $(CXXFLAGS) $< -o  $@ $(LDFLAGS) $(ROOTLIBS) -lrt


#  -lpdflib804 -llhapdf        
# DO NOT DELETE

#LQlimit.o: ./TSysLimitChannel.h ./TSysLimit.h ./TSysLimitResult.h
#LQlimit.o: ./TSysLimitScan.h LQlimit.h weight.h
#LQlimit2.o: ./TSysLimitChannel.h ./TSysLimit.h ./TSysLimitResult.h
#LQlimit2.o: ./TSysLimitScan.h LQlimit2.h weight2.h
# TSysLimit.o: ./TSysLimit.h ./TSysLimitChannel.h ./TSysLimitResultPBock.h
# TSysLimit.o: ./TSysLimitResult.h ./TSysLimitResultLogL.h
# TSysLimitChannel.o: ./TSysLimitChannel.h
# TSysLimitResult.o: ./TSysLimitResult.h
# TSysLimitResultLogL.o: ./TSysLimitResultLogL.h ./TSysLimitResult.h
# TSysLimitResultPBock.o: ./TSysLimitResultPBock.h ./TSysLimitResult.h
# TSysLimitScan.o: ./TSysLimit.h ./TSysLimitScan.h ./TSysLimitResult.h
# testlimit.o: ./TSysLimitChannel.h ./TSysLimit.h ./TSysLimitResult.h
testlimit.o: 
