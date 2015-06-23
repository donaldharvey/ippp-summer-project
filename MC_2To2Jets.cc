// -*- C++ -*-
#include "Rivet/Analyses/MC_JetAnalysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {
  /// @brief MC validation analysis for jet events
  class MC_2To2Jets : public Analysis {
  public:

    MC_2To2Jets()
      : Analysis("MC_2To2Jets")
    {    }


  public:

    void init() {
      FinalState fs;
      FastJets jetpro(fs, FastJets::ANTIKT, 0.4);
      addProjection(jetpro, "Jets");

      _histDY_12 = bookHisto1D("DY_12", 100, 0, 100);
      _histDPhi_12 = bookHisto1D("DPhi_12", 100, 0, TWOPI);

    }


    void analyze(const Event& event) {
      const double weight = event.weight();
      const Jets jets = applyProjection<FastJets>(event, "Jets").jetsByPt(0*GeV);
      const double dy_12 = jets[0].pseudorapidity() - jets[1].pseudorapidity();
      const double dPhi_12 = jets[0].phi() - jets[1].phi();
      _histDY_12->fill(dy_12, weight);
      _histDPhi_12->fill(dPhi_12, weight);

    }


    void finalize() {
      normalize(_histDY_12);
      normalize(_histDPhi_12);
    }

  private:
    Histo1DPtr _histDY_12;
    Histo1DPtr _histDPhi_12;

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_2To2Jets);

}
