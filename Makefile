
all: Rivet2To2Jets.so
Rivet2To2Jets.so: 
	rivet-buildplugin RivetMC_2To2Jets.so MC_2To2Jets.cc
run:
	RIVET_ANALYSIS_PATH=$(PWD) Herwig++ read LHC-2to2.in
	RIVET_ANALYSIS_PATH=$(PWD) Herwig++ run LHC-2to2.run -N1000 -d1
	rivet-mkhtml LHC-2to2.yoda