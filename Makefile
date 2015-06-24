
all: Rivet2To2Jets.so
Rivet2To2Jets.so: 
	rivet-buildplugin RivetMC_2To2Jets.so MC_2To2Jets.cc
generate_data:
	Herwig++ read LHC-2to2.in
	Herwig++ run LHC-2to2.run -N1000 -d1
check_data:
ifeq ($(wildcard LHC-2to2.hepmc),)
	$(MAKE) generate_data
endif
run: check_data
	RIVET_ANALYSIS_PATH=$(PWD) rivet -a MC_2To2Jets -l Rivet.Analysis=DEBUG -H LHC-2to2.yoda LHC-2to2.hepmc
	rivet-mkhtml LHC-2to2.yoda