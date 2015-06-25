all: RivetMC_2To2Jets.so

RivetMC_2To2Jets.so: MC_2To2Jets.cc
	rivet-buildplugin RivetMC_2To2Jets.so MC_2To2Jets.cc -Wextra -g

LHC-2to2.hepmc: LHC-2to2.in
	Herwig++ read LHC-2to2.in
	Herwig++ run LHC-2to2.run -N100000 -d1

run: RivetMC_2To2Jets.so LHC-2to2.hepmc
	RIVET_ANALYSIS_PATH=$(PWD) rivet -a MC_2To2Jets -l Rivet.Analysis=DEBUG -H LHC-2to2.yoda LHC-2to2.hepmc
	rivet-mkhtml LHC-2to2.yoda

debug: RivetMC_2To2Jets.so
	lldb -- python $(shell which rivet) --analysis-path=$(PWD) -a MC_2To2Jets -l Rivet.Analysis=DEBUG LHC-2to2.hepmc