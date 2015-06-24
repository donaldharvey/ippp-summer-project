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
    {  
      setNeedsCrossSection(true);
    }


  public:

    void init() {
      FinalState fs;
      FastJets jetpro(fs, FastJets::ANTIKT, 0.4);
      addProjection(jetpro, "Jets");

      for (unsigned long i=0;i<3;++i) {
        const double pTmax = 1.0/(double(i)+2.0) * sqrtS()/GeV/2.0;
        const unsigned long nbins_pT = 100/(i+1);
        const string pTname = "Jet_" + to_str(i+1) + "_p_T";
        MSG_DEBUG("pTmax = " << pTmax);
        MSG_DEBUG("nBins_pT = " << nbins_pT);
        MSG_DEBUG("pTname = " << pTname);
        
        _histJetP_T.push_back(bookHisto1D(pTname, logspace(nbins_pT, 10.0, pTmax)));
        for (unsigned long j=i+1;j<3;++j) {
          const std::pair<int, int> ij = std::make_pair(i, j);
          _histInterJetDeltaEta.insert(make_pair(ij, bookHisto1D(
            "Jet_" + to_str(i+1) + to_str(j+1) + "_DeltaEta", 100, 0, 10
          )));
          _histInterJetDeltaPhi.insert(make_pair(ij, bookHisto1D(
            "Jet_" + to_str(i+1) + to_str(j+1) + "_DeltaPhi", 100, 0, M_PI
          )));
          _histInterJetDeltaR.insert(make_pair(ij, bookHisto1D(
            "Jet_" + to_str(i+1) + to_str(j+1) + "_DeltaR", 100, 0, 15
          )));
        }
      }
    }


    void analyze(const Event& event) {
      const double weight = event.weight();
      const Jets jets = applyProjection<FastJets>(event, "Jets").jetsByPt(0*GeV);
      for (unsigned long i=0;i<3;i++) {
        if (jets.size() < i+1) continue;
        _histJetP_T[i]->fill(jets[i].pt());
        for (unsigned long j=i+1;j<3;++j) {
          const std::pair<int, int> ij = std::make_pair(i, j);
          _histInterJetDeltaEta[ij]->fill(jets[i].pseudorapidity() - jets[j].pseudorapidity(), weight);
          _histInterJetDeltaPhi[ij]->fill(deltaPhi(jets[i].momentum(), jets[j].momentum()), weight);
          _histInterJetDeltaR[ij]->fill(deltaR(jets[i].momentum(), jets[j].momentum()), weight);
        }
      }
    }


    void finalize() {
      for (unsigned long i=0;i<3;i++) {
        normalize(_histJetP_T[i]);
      }
      // Scale the d{eta,phi,R} histograms
      typedef map<pair<int, int>, Histo1DPtr> HistMap;
      foreach (HistMap::value_type& it, _histInterJetDeltaEta) scale(it.second, crossSection()/sumOfWeights());
      foreach (HistMap::value_type& it, _histInterJetDeltaPhi) scale(it.second, crossSection()/sumOfWeights());
      foreach (HistMap::value_type& it, _histInterJetDeltaR) scale(it.second, crossSection()/sumOfWeights());
    }

  private:
    std::map<std::pair<int, int>, Histo1DPtr> _histInterJetDeltaEta;
    std::map<std::pair<int, int>, Histo1DPtr> _histInterJetDeltaPhi;
    std::map<std::pair<int, int>, Histo1DPtr> _histInterJetDeltaR;
    std::vector<Histo1DPtr> _histJetP_T;
  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_2To2Jets);

}
