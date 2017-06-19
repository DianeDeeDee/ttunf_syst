
ROOTCFLAGS    = $(shell root-config --cflags)

FJVER=3.0.3
FJCONTRIBVER=1.011
MYCF=""

ifneq ($(findstring m32,$(ROOTCFLAGS)),)
        MYCF= CXXFLAGS=-m32 CFLAGS=-m32
endif

FASTJETDIR = $(shell fastjet-install/bin/fastjet-config --prefix)
FASTJETLIB = $(shell fastjet-install/bin/fastjet-config --libs --plugins)

ROOTLIB    = $(shell root-config --libs)

SDDIR=..
SDLIB=-L$(SDDIR)/SD-tosvn/.libs 
#-lDeconstruction

CXXFLAGS   = -g -I.
CXXFLAGS  += -Wno-long-long -fPIC
CXXFLAGS  += $(shell root-config --cflags)
CXXFLAGS  += -I$(FASTJETDIR)/include

LDFLAGS    =
LDFLAGS   += $(ROOTLIB) -lCintex
LDFLAGS   += $(FASTJETLIB)

ROOTCOREDIR   = packages/RootCore

PACKAGES      = GoodRunsLists egammaAnalysisUtils CalibrationDataInterface MuonMomentumCorrections MissingETUtility ApplyJetCalibration ApplyJetResolutionSmearing JetEffiProvider JetResolution JetUncertainties MuonEfficiencyCorrections MuonMomentumCorrections PileupReweighting TopElectronSFUtils TopJetUtils TopMuonSFUtils egammaEvent TopDataPreparation WjetsCorrections ElectronEfficiencyCorrection JVFUncertaintyTool ttResoSingleLepton

#PACKAGES      = GoodRunsLists egammaAnalysisUtils CalibrationDataInterface MuonMomentumCorrections MissingETUtility ApplyJetCalibration ApplyJetResolutionSmearing JetEffiProvider JetResolution JetUncertainties MuonEfficiencyCorrections MuonMomentumCorrections PileupReweighting TopElectronSFUtils TopJetUtils TopMuonSFUtils egammaEvent TopDataPreparation WjetsCorrections ttResoSingleLepton ElectronEfficiencyCorrection JVFUncertaintyTool

PACKAGES     += TrigMuonEfficiency
CXXFLAGS     += -I$(ROOTCOREBIN)/include $(shell $(ROOTCOREDIR)/scripts/get_cxxflags.sh $(PACKAGES))
LDFLAGS      += $(shell $(ROOTCOREDIR)/scripts/get_ldflags.sh $(PACKAGES))

OBJS       = Correction.o Event.o LargeJet.o Muon.o
OBJS      += ParseUtils.o
OBJS      += SkimReader.o
OBJS      += Particle.o
OBJS      += Electron.o Jet.o MObject.o Reader.o RawReader.o
OBJS      += EventCutter.o EventCutterReco.o EventCutterPart.o EventCutterRecoDJ.o EventCutterPartDJ.o  EventCutterRecoBB.o EventCutterPartBB.o
OBJS      += MiniTree.o
OBJS      += Tools.o
OBJS      += AllCorrections.o
OBJS      += WeightCalc.o AllWeightCalcs.o
OBJS_READ += $(OBJS)
#OBJS_READ += RawReader.cxx SkimReader.cxx
#OBJS_READ += Event.o LargeJet.o Muon.o
#OBJS_READ += ParseUtils.o Particle.o Electron.o
#OBJS_READ += Jet.o MObject.o MiniTree.o Tools.o
OBJS_READ += EventCount.o
OBJS_READ += Plot.o PlotSemilep.o HistogramService.o
OBJS_READ += read.o 
OBJS_READ += FakeQCDMMWeight.o FakeQCDMMResolved.o FakeQCDMMBoosted.o

OBJS_RUNSD += $(OBJS)
OBJS_RUNSD += EventCount.o
OBJS_RUNSD += Plot.o PlotSemilep.o HistogramService.o
#PlotTTJet.o
OBJS_RUNSD += runsd.o
#OBJS_RUNSD += $(SDDIR)/SD/.libs/libDeconstruction.a
OBJS_RUNSD += FakeQCDMMWeight.o FakeQCDMMResolved.o FakeQCDMMBoosted.o
OBJS_RUNSD += fastjet-install/lib/libNsubjettiness.a

OBJS_PRE  += $(OBJS)
OBJS_PRE  += preselect.o

%.o: %.cxx
	g++ -c $(CXXFLAGS) -o $@ $<


all: fastjet preselect read runsd 

SkimReader.o: SkimReader.cxx SkimReader.h
	g++ -c $(CXXFLAGS) -o $@ SkimReader.cxx


fastjet: fastjet-install/lib/libfastjet.so.0.0.0 fjcontrib

fastjet-install/lib/libfastjet.so.0.0.0 : fastjet3.tar.gz
	tar -zxvf fastjet3.tar.gz
	cd fastjet-$(FJVER); ./configure --prefix `pwd`/../fastjet-install $(MYCF)
	cd fastjet-$(FJVER) ; make; make install

fjcontrib: fastjet-install/lib/libNsubjettiness.a

fastjet-install/lib/libNsubjettiness.a: fjcontrib-1.011.tar.gz
	rm -rf fjcontrib-$(FJCONTRIBVER)
	tar -zxvf fjcontrib-1.011.tar.gz
	cd fjcontrib-$(FJCONTRIBVER) ; ./configure --fastjet-config=`pwd`/../fastjet-install/bin/fastjet-config
	cd fjcontrib-$(FJCONTRIBVER) ; make ; make install


preselect: $(OBJS_PRE) fastjet-install/lib/libfastjet.so.0.0.0 fjcontrib 
	g++ $(CXXFLAGS) -o preselect $(OBJS_PRE) $(LDFLAGS)

read: $(OBJS_READ) fastjet-install/lib/libfastjet.so.0.0.0 fjcontrib
	g++ $(CXXFLAGS) -o read $(OBJS_READ) $(LDFLAGS)

runsd: $(OBJS_RUNSD) fastjet-install/lib/libfastjet.so.0.0.0 fjcontrib
	g++ $(CXXFLAGS) -o runsd $(OBJS_RUNSD) $(LDFLAGS) $(SDLIB)




clean:
	rm -rf *.o preselect read runsd


