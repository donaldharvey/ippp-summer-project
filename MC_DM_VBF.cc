// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
/// @todo Include more projections as required, e.g. ChargedFinalState, FastJets, ZFinder...

#define MAIN_JET_PT_CUT 30.0*GeV
#define VBF_JET_1_MIN_PT 75.0*GeV
#define VBF_JET_2_MIN_PT 50.0*GeV
#define ELECTRON_MAX_PT 10.0*GeV
#define MUON_MAX_PT 5.0*GeV
#define TAUON_MAX_PT 20.0*GeV
#define MAX_VBF_DELTAPHI 2.5
#define MIN_VBF_DELTAETA 4.8
#define MIN_VBF_INVT_MASS 1*TeV
// add jet cut - max eta 4.8ish!
namespace Rivet {
  double calculateDeltaPhi(const Jet& jet1, const Jet& jet2) {
    return std::abs(jet1.pseudojet().delta_phi_to(jet2.pseudojet()));
  }
  class MC_DM_VBF : public Analysis {
  public:

    /// Constructor
    MC_DM_VBF()
      : Analysis("MC_DM_VBF")

    {    }


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      

      /// @todo Initialise and register projections here
      VetoedFinalState fs;
      fs.addVetoId(PID::HIGGS);
      FastJets jetpro(fs, FastJets::ANTIKT, 0.5);
      addProjection(jetpro, "Jets");
      /// @todo Book histograms here, e.g.:
      // _h_XXXX = bookProfile1D(1, 1, 1);
      _h_eta1 = bookHisto1D("Eta1", 50, -5, 5);
      _h_eta2 = bookHisto1D("Eta2", 50, -5, 5);
      _h_eta3 = bookHisto1D("Eta3", 50, -5, 5);
      _h_phi1 = bookHisto1D("Phi1", 50, 0, TWOPI);
      _h_phi2 = bookHisto1D("Phi2", 50, 0, TWOPI);
      _h_phi3 = bookHisto1D("Phi3", 50, 0, TWOPI);
      _h_pt1 = bookHisto1D("pT_1", 80, 0, 200);
      _h_pt2 = bookHisto1D("pT_2", 80, 0, 200);
      _h_pt3 = bookHisto1D("pT_3", 80, 0, 200);
      _h_vbf_mass = bookHisto1D("VBFDijetInvtMass", 50, 0, 3000);
      _h_deltaphi_vbf = bookHisto1D("DeltaPhi_VBF", 50, 0, PI);
      _h_deltaeta_vbf = bookHisto1D("DeltaEta_VBF", 73, 4.7, 7);
      _count = 0;
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();
      const Jets jets = applyProjection<FastJets>(event, "Jets").jetsByPt(MAIN_JET_PT_CUT);
      if (jets.size() < 2) {
        MSG_INFO("Too few jets (" << jets.size() << ")");
        vetoEvent;
      }
      // if(jets[0].pt() <= VBF_JET_1_MIN_PT or jets[1].pt() <= VBF_JET_2_MIN_PT) {
      //   MSG_INFO("Leading jet pts too small");
      //   vetoEvent;
      // }
      // check out leading jet properties.
      if (calculateDeltaPhi(jets[0], jets[1]) > MAX_VBF_DELTAPHI) {
        MSG_INFO("Leading jet delta phi too large");
        vetoEvent;
      }
      double vbfDeltaEta = std::abs(jets[0].pseudorapidity() - jets[1].pseudorapidity());
      if (vbfDeltaEta <= MIN_VBF_DELTAETA) {
        MSG_INFO("Leading jet delta eta too small");
        vetoEvent;
      }
      if (signbit(jets[0].pseudorapidity()) * signbit(jets[1].pseudorapidity()) != 0) {
        MSG_INFO("Not in oppposite hemispheres:" << signbit(jets[0].pseudorapidity()) * signbit(jets[1].pseudorapidity()));
        vetoEvent;
      }
      // if ((jets[0].pseudojet() + jets[1].pseudojet()).m() < MIN_VBF_INVT_MASS) {
      //   MSG_INFO("Dijet mass too low");
      //   vetoEvent;
      // }
      // jet masses!
      // for(size_t i=0; i<jets.size(); ++i) {
      
      // }
      // TODO: check for reconstructed leptons.
      _h_eta1->fill(jets[0].pseudorapidity(), weight);
      _h_eta2->fill(jets[1].pseudorapidity(), weight);
      if (jets.size() > 2) {
        _h_eta3->fill(jets[2].pseudorapidity(), weight);
        _h_phi3->fill(jets[2].phi(), weight);
        _h_pt3->fill(jets[2].pt(), weight);
      }
      _h_phi1->fill(jets[0].phi(), weight);
      _h_phi2->fill(jets[1].phi(), weight);
      _h_pt1->fill(jets[0].pt(), weight);
      _h_pt2->fill(jets[1].pt(), weight);
      _h_deltaeta_vbf->fill(vbfDeltaEta, weight);
      _h_deltaphi_vbf->fill(calculateDeltaPhi(jets[1], jets[0]), weight);
      _h_vbf_mass->fill((jets[0].pseudojet() + jets[1].pseudojet()).m());
      _count += 1;
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      /// @todo Normalise, scale and otherwise manipulate histograms here

      // scale(_h_YYYY, crossSection()/sumOfWeights()); // norm to cross section
      // normalize(_h_YYYY); // normalize to unity
      normalize(_h_eta1);
      normalize(_h_eta2);
      normalize(_h_eta3);
      normalize(_h_phi1);
      normalize(_h_phi2);
      normalize(_h_phi3);
      normalize(_h_deltaeta_vbf);
      normalize(_h_deltaphi_vbf);
      normalize(_h_vbf_mass);
      MSG_INFO(_count << " unvetoed events.");
    }

    //@}


  private:

    // Data members like post-cuts event weight counters go here


    /// @name Histograms
    //@{
    Histo1DPtr _h_eta1;
    Histo1DPtr _h_eta2;
    Histo1DPtr _h_eta3;
    Histo1DPtr _h_phi1;
    Histo1DPtr _h_phi2;
    Histo1DPtr _h_phi3;
    Histo1DPtr _h_pt1;
    Histo1DPtr _h_pt2;
    Histo1DPtr _h_pt3;
    Histo1DPtr _h_deltaeta_vbf;
    Histo1DPtr _h_deltaphi_vbf;
    Histo1DPtr _h_vbf_mass;
    unsigned long _count;
    //@}


  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_DM_VBF);


}
