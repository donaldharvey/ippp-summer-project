// -*- C++ -*-
#include "Rivet/Analyses/MC_JetAnalysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"

// Cuts:

// pT1 > 100 GeV
// pT3 > 30 GeV
#define MIN_PT_CUT 30.0*GeV
#define PT1_CUT 100*GeV

// |eta_1|, |eta_2| <= 2.5
#define ETA_CUT 2.5

// M_12 > 220 GeV
#define M_CUT 220*GeV

// 0.5 < R_23 < 1.5
#define R_23_CUT_MIN 0.5
#define R_23_CUT_MAX 1.5

namespace Rivet {
  double calculateDeltaPhi(const Jet& jet1, const Jet& jet2) {
    return jet1.pseudojet().delta_phi_to(jet2.pseudojet());
  }
  double calculateBeta(const Jets& jets) {
    const double deltaEta = jets[2].pseudorapidity() - jets[1].pseudorapidity();
    const double dPhi = calculateDeltaPhi(jets[2].momentum(), jets[1].momentum());
    if (deltaEta == 0) {
      return HALFPI;
    }
    if (dPhi == 0) {
      return deltaEta > 0 ? 0 : PI;
    }
    return atan2(std::abs(dPhi), deltaEta);
  }
  double calculateDeltaR(const Jet& jet1, const Jet& jet2) {
    return jet1.pseudojet().delta_R(jet2.pseudojet());
  }
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
      cutEventCount = 0;
      VetoedFinalState fs;
      fs.addVetoID(PID::HIGGS);
      FastJets jetpro(fs, FastJets::ANTIKT, 0.5);
      addProjection(jetpro, "Jets");
      _histNumJets = bookHisto1D("NumJets", 5, 3, 5);
      _histDeltaPhi23Central = bookHisto1D("DeltaPhi23Central", 50, -2, 2);
      _histDeltaPhi23Forward = bookHisto1D("DeltaPhi23Forward", 50, -2, 2);
      _histDeltaEta23Central = bookHisto1D("DeltaEta23Central", 50, -2, 2);
      _histDeltaEta23Forward = bookHisto1D("DeltaEta23Forward", 50, -2, 2);
      _histBetaCentral = bookHisto1D("BetaCentral", 50, 0, PI);
      _histBetaForward = bookHisto1D("BetaForward", 50, 0, PI);

      for (unsigned long i=0;i<3;++i) {
        const double pTmax = 1.0/(double(i)+2.0) * sqrtS()/GeV/2.0;
        const unsigned long nbins_pT = 100/(i+1);
        const string pTname = "Jet_" + to_str(i+1) + "_p_T";
        
        _histJetP_T.push_back(bookHisto1D(pTname, 50, 0, 700));
        for (unsigned long j=i+1;j<3;++j) {
          const std::pair<int, int> ij = std::make_pair(i, j);
          _histInterJetDeltaEta.insert(make_pair(ij, bookHisto1D(
            "Jet_" + to_str(i+1) + to_str(j+1) + "_DeltaEta", 50, -5, 5
          )));
          _histInterJetDeltaPhi.insert(make_pair(ij, bookHisto1D(
            "Jet_" + to_str(i+1) + to_str(j+1) + "_DeltaPhi", 50, -PI, PI
          )));
          _histInterJetDeltaR.insert(make_pair(ij, bookHisto1D(
            "Jet_" + to_str(i+1) + to_str(j+1) + "_DeltaR", 100, 0, 6
          )));
        }
      }
    }


    void analyze(const Event& event) {
      const double weight = event.weight();
      const Jets jets = applyProjection<FastJets>(event, "Jets").jetsByPt(MIN_PT_CUT);
      // apply the cuts
      if(jets.size() < 3) {
        MSG_DEBUG("Vetoing: not enough jets. (Got " << jets.size() << ")");
        vetoEvent;
      }
      if(jets[0].pt() <= PT1_CUT) {
        MSG_DEBUG("Vetoing: first jet has insufficient pt.");
        vetoEvent;
      }
      if(std::abs(jets[0].pseudorapidity()) > ETA_CUT or std::abs(jets[1].pseudorapidity()) > ETA_CUT) {
        MSG_DEBUG("Vetoing: pseudorapidities too large.");
        vetoEvent;
      }
      if((jets[0].pseudojet() + jets[1].pseudojet()).m() <= M_CUT) {
        MSG_DEBUG("Vetoing: insufficient mass.");
        vetoEvent;
      }
      double deltaR23 = calculateDeltaR(jets[1], jets[2]);
      if(deltaR23 >= R_23_CUT_MAX or deltaR23 <= R_23_CUT_MIN) {
        MSG_DEBUG("Vetoing: deltaR23 not in correct range.");
        vetoEvent;
      }
      cutEventCount += 1;
      _histNumJets->fill(jets.size(), weight);
      for (unsigned long i=0;i<3;i++) {
        if (jets.size() < i+1) continue;
        _histJetP_T[i]->fill(jets[i].pt());
        for (unsigned long j=i+1;j<3;++j) {
          const std::pair<int, int> ij = std::make_pair(i, j);
          double deltaEta = sign(jets[i].pseudorapidity()) * (jets[j].pseudorapidity() - jets[i].pseudorapidity());

          _histInterJetDeltaEta[ij]->fill(deltaEta, weight);
          if (i == 1 and j == 2) {
            // fill the 23 histograms, separating on eta2
            double eta2 = jets[i].pseudorapidity();
            MSG_DEBUG("Got deltaEta of " << deltaEta);
            if (std::abs(eta2) > 0.8) {
              _histDeltaEta23Forward->fill(deltaEta, weight);
              _histBetaForward->fill(calculateBeta(jets), weight);
              _histDeltaPhi23Forward->fill(calculateDeltaPhi(jets[j], jets[i]), weight);
            }
            else {
              _histDeltaEta23Central->fill(deltaEta, weight);
              _histBetaCentral->fill(calculateBeta(jets), weight);
              _histDeltaPhi23Central->fill(calculateDeltaPhi(jets[j], jets[i]), weight);
            }
          }
          _histInterJetDeltaPhi[ij]->fill(calculateDeltaPhi(jets[j], jets[i]), weight);
          _histInterJetDeltaR[ij]->fill(calculateDeltaR(jets[j], jets[i]), weight);
        }
      }
    }


    void finalize() {
      MSG_DEBUG("Got " << cutEventCount << " unvetoed events.");
      for (unsigned long i=0;i<3;i++) {
        normalize(_histJetP_T[i]);
      }
      normalize(_histNumJets);
      // Scale the d{eta,phi,R} histograms
      typedef map<pair<int, int>, Histo1DPtr> HistMap;
      foreach (HistMap::value_type& it, _histInterJetDeltaEta) normalize(it.second);
      foreach (HistMap::value_type& it, _histInterJetDeltaPhi) normalize(it.second);
      foreach (HistMap::value_type& it, _histInterJetDeltaR) normalize(it.second);
      normalize(_histDeltaEta23Central);
      normalize(_histDeltaEta23Forward);
      normalize(_histDeltaPhi23Central);
      normalize(_histDeltaPhi23Forward);
    }

  private:
    std::map<std::pair<int, int>, Histo1DPtr> _histInterJetDeltaEta;
    std::map<std::pair<int, int>, Histo1DPtr> _histInterJetDeltaPhi;
    std::map<std::pair<int, int>, Histo1DPtr> _histInterJetDeltaR;
    std::vector<Histo1DPtr> _histJetP_T;
    Histo1DPtr _histNumJets;
    Histo1DPtr _histDeltaEta23Central;
    Histo1DPtr _histDeltaEta23Forward;
    Histo1DPtr _histDeltaPhi23Central;
    Histo1DPtr _histDeltaPhi23Forward;
    Histo1DPtr _histBetaCentral;
    Histo1DPtr _histBetaForward;
    unsigned long cutEventCount;
  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_2To2Jets);

}
