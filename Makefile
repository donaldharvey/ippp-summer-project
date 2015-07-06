all: RivetMC_2To2Jets.so RivetMC_VBF.so

RivetMC_DM_VBF.so: MC_DM_VBF.cc
	rivet-buildplugin RivetMC_DM_VBF.so MC_DM_VBF.cc -Wextra -g

RivetMC_2To2Jets.so: MC_2To2Jets.cc
	rivet-buildplugin RivetMC_2To2Jets.so MC_2To2Jets.cc -Wextra -g

RivetMC_VBF.so: MC_VBF.cc
	rivet-buildplugin RivetMC_VBF.so MC_VBF.cc -Wextra -g

LHC-2to2.hepmc: LHC-2to2.in
	Herwig++ read LHC-2to2.in
	Herwig++ run LHC-2to2.run -N10000 -d1

LHC-VBF.hepmc: LHC-VBF.in
	Herwig++ read LHC-VBF.in
	Herwig++ run LHC-VBF.run -N10000 -d1	

run: RivetMC_2To2Jets.so LHC-2to2.hepmc
	RIVET_ANALYSIS_PATH=$(PWD) rivet -a MC_2To2Jets -l Rivet.Analysis=DEBUG -H LHC-2to2.yoda LHC-2to2.hepmc
	rivet-mkhtml LHC-2to2.yoda

runvbf: RivetMC_VBF.so LHC-VBF.hepmc
	RIVET_ANALYSIS_PATH=$(PWD) rivet -a MC_VBF -H LHC-VBF.yoda LHC-VBF.hepmc
	rivet-mkhtml LHC-VBF.yoda

rundmvbf: RivetMC_DM_VBF.so LHC-VBF.hepmc
	RIVET_ANALYSIS_PATH=$(PWD) rivet -a MC_DM_VBF -H LHC-DM-VBF.yoda LHC-VBF.hepmc
	rivet-mkhtml LHC-DM-VBF.yoda

debug: RivetMC_2To2Jets.so
	lldb -- python $(shell which rivet) --analysis-path=$(PWD) -a MC_2To2Jets -l Rivet.Analysis=DEBUG LHC-2to2.hepmc

debugdmvbf: RivetMC_DM_VBF.so
	lldb -- python $(shell which rivet) --analysis-path=$(PWD) -a MC_DM_VBF -l Rivet.Analysis=DEBUG LHC-VBF.hepmc	