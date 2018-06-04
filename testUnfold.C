#define testUnfold_cxx
#include "testUnfold.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>


#include "TUnfoldDensity.h"


using namespace std;



void testUnfold::Loop()
{
//   In a ROOT session, you can do:
//      root> .L testUnfold.C
//      root> testUnfold t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch


 TFile* output = new TFile( "test.root", "recreate");


// switch on histogram errors
TH1::SetDefaultSumw2();

// show fit result
gStyle->SetOptFit(1111);

Int_t const nDet=8;
Int_t const nGen=4;
Double_t const xminDet=10.0;
Double_t const xmaxDet=150.0;
Double_t const xminGen=10.0;
Double_t const xmaxGen=150.0;


//===========================================
// Define Histograms
//
 TH1D *histNgenJets = new TH1D("NgenJets",";njets",4, 0, 4);
 TH1D *histGenRecoEff=new TH1D("GenRecoEff",";ptl1(gen)",nGen,xminGen,xmaxGen);
 TH1D *histMgenMC=new TH1D("MgenMC",";ptl1(gen)",nGen,xminGen,xmaxGen);
 TH1D *histMgenData=new TH1D("MgenData",";ptl1(gen)",nGen,xminGen,xmaxGen);
 TH1D *histMdetData=new TH1D("MdetData",";ptl1(gen)",nDet,xminDet,xmaxDet);
 TH1D *histMdetMC=new TH1D("MdetMC",";ptl1(det)",nDet,xminDet,xmaxDet);
 TH2D *histMdetGenMC=new TH2D("MdetgenMC",";ptl1(det);ptl1(gen)",
			      nDet,xminDet,xmaxDet,nGen,xminGen,xmaxGen);



//=========================================================================
// divide by bin withd to get density distributions

 TH1D *histDensityGenData=new TH1D("DensityGenData",";ptl1(gen)",
                                    nGen,xminGen,xmaxGen);
 TH1D *histDensityGenMC=new TH1D("DensityGenMC",";ptl1(gen)",
				 nGen,xminGen,xmaxGen);

 

 if (fChain == 0) return;

 Long64_t nentries = fChain->GetEntries();


 Long64_t nbytes = 0, nb = 0;
 
 for (Long64_t jentry=0; jentry<nentries;jentry++) {

   Long64_t ientry = LoadTree(jentry);
  
   if (ientry < 0) break;
   
   nb = fChain->GetEntry(jentry);   nbytes += nb;
   
   // if (Cut(ientry) < 0) continue;
  


// Define weights 

 bool _ismc = 1; 

 Float_t luminosity = 39.0;

 bool trigger = 1; 
 //PassTrigger();

 // Float_t metFilter = (_ismc) ? METFilter_MC : METFilter_DATA;
 Float_t metFilter =  METFilter_MC;

 Float_t genWeight = 1.0; 
 if (GEN_weight_SM) genWeight = GEN_weight_SM / abs(GEN_weight_SM);

 float sf_idiso = 1.0;
 sf_idiso = LepSF2l__ele_cut_WP_Tight80X__mu_cut_Tight80x;

 float sf_reco = 1;//std_vector_lepton_recoW->at(0) * std_vector_lepton_recoW->at(1);

 Double_t totalWGen = baseW * luminosity * genWeight; //Removing PU from gen level

 Double_t totalWReco = sf_idiso * sf_reco * genWeight * metFilter * effTrigW * veto_EMTFBug * puW * baseW * luminosity;// * effW * triggW ;//puW *  mybaseW * luminosity;//puW * effW * triggW * mybaseW * luminosity;

 Double_t totalW = puW * baseW * luminosity * genWeight;//efficiencyW;


   // The GEN selection begins here
   //--------------------------------------------------------------------------
   
    /// ---> 1) Need deressed leptons to define the same fiducial region
    /// ---> 2) Count how many GEN leptons we have in each bin, applying the fidual region cuts
    /// ---> 3) Apply also, OF, jetbin and opposite-charged cuts.
   
 
   float genPtl1 = std_vector_leptonGen_pt->at(0);
   float genPtl2 = std_vector_dressedLeptonGen_pt->at(1);

   float genEta1 = std_vector_dressedLeptonGen_eta->at(0);
   float genEta2 = std_vector_dressedLeptonGen_eta->at(1);
  
   float genPID1 = std_vector_dressedLeptonGen_pid->at(0); 
   float genPID2 = std_vector_dressedLeptonGen_pid->at(1); 



 // DO JET VETO 

 int nGenJets = 0;

 for (unsigned int j= 0; j < std_vector_jetGen_pt->size(); j++ ) {

   if ( std_vector_jetGen_pt->at(j) < 30 ) continue;
   //if ( fabs(std_vector_jetGen_eta->at(j)) > 4.7 ) continue;
  
   nGenJets++;
 }

 histNgenJets->Fill(nGenJets);

 bool isGenEvent = 0;

 if ( genPtl1 > 10  && 
      genPtl2 > 10  && 
      //fabs(lepGenM1)== 24 &&  
      //fabs(lepGenM2)== 24 &&
      fabs(genPID1) != 15 &&
      fabs(genPID2) != 15 &&
      fabs(genPID1) != fabs(genPID2) &&  
      (genPID1*genPID2) < 0 && 
      ((fabs(genPID1) == 13 && fabs(genEta1) < 2.4) || 
       (fabs(genPID1) == 11 && fabs(genEta1) < 2.5)) &&  
      ((fabs(genPID2) == 13 && fabs(genEta2) < 2.4) || 
       (fabs(genPID2) == 11 && fabs(genEta2) < 2.5)) ) {

   isGenEvent = 1; 
}

 if (!isGenEvent) continue;


 if ( genPtl1 >= 150) histMgenMC->Fill(149, totalWGen);

 // generated MC distribution (for comparison only)
 histMgenMC->Fill(genPtl1, totalWGen);



 /// ---> 3) Going to apply the selection analysis cuts on RECO objects

 float recoPtl1 = std_vector_lepton_pt->at(0);
 float recoPtl2 = std_vector_lepton_pt->at(1);
 float recoPtl3 = std_vector_lepton_pt->at(2);

 bool pass_os = (std_vector_lepton_flavour->at(0) * std_vector_lepton_flavour->at(1) < 0);

 bool pass_bveto = true;

 for (unsigned  int h=0; h < std_vector_jet_pt->size(); h++)
   {
     pass_bveto &= (std_vector_jet_pt->at(h) < 20 || std_vector_jet_cmvav2->at(h) < -0.5884);
   }

 bool isRecoEvent = 0;


 if(  recoPtl1 > 10 && recoPtl2 > 10  && pass_os &&
      recoPtl3 < 10                              && 
      metPfType1 > 20                            &&
      mll > 12                                   &&
      mth > 55                                   &&
      ptll > 30                                  &&     
      pass_bveto                                 ) { 
	 
   isRecoEvent = 1; 	   
 }


 if (!isRecoEvent) continue;



 // reconstructed MC distribution (for comparison only)
 histMdetMC->Fill(recoPtl1, totalWReco);

 histMdetData->Fill(recoPtl1,totalWReco);

 histMdetGenMC->Fill(recoPtl1, genPtl1, totalWReco);

 histMgenData->Fill(genPtl1, totalWReco);


}  //End of loop 


 histGenRecoEff->Divide(histMgenMC,histMgenData);

// divide by bin withd to get density distributions

for(Int_t i=1;i<=nGen;i++) {

  histDensityGenData->SetBinContent(i,histMdetData->GetBinContent(i)/
                                       histMdetData->GetBinWidth(i));
  histDensityGenMC->SetBinContent(i,histMgenMC->GetBinContent(i)/
				  histMgenMC->GetBinWidth(i));
 } 

//=========================================================================
// set up the unfolding
// define migration matrix


 TUnfoldDensity unfold(histMdetGenMC,TUnfold::kHistMapOutputVert);


  // define input and bias scame
  // do not use the bias, because MC peak may be at the wrong place
  // watch out for error codes returned by the SetInput method
  // errors larger or equal 10000 are fatal:
  // the data points specified as input are not sufficient to constrain the
  // unfolding process
 ///  if(
 unfold.SetInput(histMdetData);//(>=10000) {
  //  std::cout<<"Unfolding result may be wrong\n";
  // }


//========================================================================
  // the unfolding is done here
  //
  // scan L curve and find best point
  Int_t nScan=30;
  // use automatic L-curve scan: start with taumin=taumax=0.0
  Double_t tauMin=0.0;
  Double_t tauMax=0.0;
  Int_t iBest;
  TSpline *logTauX,*logTauY;
  TGraph *lCurve;

 // if required, report Info messages (for debugging the L-curve scan)
#ifdef VERBOSE_LCURVE_SCAN
  Int_t oldinfo=gErrorIgnoreLevel;
  gErrorIgnoreLevel=kInfo;
#endif


 // this method scans the parameter tau and finds the kink in the L curve
  // finally, the unfolding is done for the best choice of tau
  iBest=unfold.ScanLcurve(nScan,tauMin,tauMax,&lCurve,&logTauX,&logTauY);


 // if required, switch to previous log-level
#ifdef VERBOSE_LCURVE_SCAN
  gErrorIgnoreLevel=oldinfo;
#endif


  //==========================================================================
  // define a correlated systematic error
  // for example, assume there is a 10% correlated error for all reconstructed
  // masses larger than 7

  //==========================================================================
  // print some results
  //
  std::cout<<"tau="<<unfold.GetTau()<<"\n";
  std::cout<<"chi**2="<<unfold.GetChi2A()<<"+"<<unfold.GetChi2L()
           <<" / "<<unfold.GetNdf()<<"\n";
  std::cout<<"chi**2(sys)="<<unfold.GetChi2Sys()<<"\n";

  //==========================================================================
  // create graphs with one point to visualize the best choice of tau
  //
  Double_t t[1],x[1],y[1];
  logTauX->GetKnot(iBest,t[0],x[0]);
  logTauY->GetKnot(iBest,t[0],y[0]);
  TGraph *bestLcurve=new TGraph(1,x,y);
  TGraph *bestLogTauLogChi2=new TGraph(1,t,x);

  //==========================================================================
  // retreive results into histograms

  // get unfolded distribution
  TH1 *histMunfold=unfold.GetOutput("Unfolded");

  TH1 *histMunfoldEff = NULL;


  histMunfold->Multiply(histGenRecoEff);

  //  histMunfoldEff->Multiply(histMunfold, histGenRecoEff);

  // get unfolding result, folded back
  TH1 *histMdetFold=unfold.GetFoldedOutput("FoldedBack");

  // get error matrix (input distribution [stat] errors only)
  // TH2D *histEmatData=unfold.GetEmatrix("EmatData");

  // get total error matrix:
  //   migration matrix uncorrelated and correlated systematic errors
  //   added in quadrature to the data statistical errors
  TH2 *histEmatTotal=unfold.GetEmatrixTotal("EmatTotal");

  // create data histogram with the total errors
  TH1D *histTotalError=
     new TH1D("TotalError",";mass(gen)",nGen,xminGen,xmaxGen);
  for(Int_t bin=1;bin<=nGen;bin++) {
    histTotalError->SetBinContent(bin,histMunfold->GetBinContent(bin));
    histTotalError->SetBinError
       (bin,sqrt(histEmatTotal->GetBinContent(bin,bin)));
  }


// get global correlation coefficients
  // for this calculation one has to specify whether the
  // underflow/overflow bins are included or not
  // default: include all bins
  // here: exclude underflow and overflow bins

  TH2 *gHistInvEMatrix;

  TH1 *histRhoi=unfold.GetRhoItotal("rho_I",
                                    0, // use default title
                                    0, // all distributions
                                    "*[UO]", // discard underflow and overflow bins on all axes
                                    kTRUE, // use original binning
                                    &gHistInvEMatrix // store inverse of error matrix
                                    );


//=====================================================================
  // plot some histograms
  TCanvas outputC;
  outputC.Divide(3,2);

  // Show the matrix which connects input and output
  // There are overflow bins at the bottom, not shown in the plot
  // These contain the background shape.
  // The overflow bins to the left and right contain
  // events which are not reconstructed. These are necessary for proper MC
  // normalisation
  outputC.cd(1);
  histMdetGenMC->Draw("BOX");

  // draw generator-level distribution:
  //   data (red) [for real data this is not available]
  //   MC input (black) [with completely wrong peak position and shape]
  //   unfolded data (blue)
  outputC.cd(2);
  histTotalError->SetLineColor(kBlue);
  histTotalError->Draw("E");
  histMunfold->SetLineColor(kGreen);
  histMunfold->Draw("SAME E1");
  //  histDensityGenData->SetLineColor(kRed);
  //histDensityGenData->Draw("SAME");
  //histDensityGenMC->Draw("SAME HIST");

  // show detector level distributions
  //    data (red)
  //    MC (black) [with completely wrong peak position and shape]
  //    unfolded data (blue)
  outputC.cd(3);
  histMdetFold->SetLineColor(kBlue);
  histMdetFold->Draw();
  histMdetMC->Draw("SAME HIST");

  TH1 *histInput=unfold.GetInput("Minput",";mass(det)");

  histInput->SetLineColor(kRed);
  histInput->Draw("SAME");

  // show correlation coefficients
  outputC.cd(4);
  histRhoi->Draw();

  // show tau as a function of chi**2
  outputC.cd(5);
  logTauX->Draw();
  bestLogTauLogChi2->SetMarkerColor(kRed);
  bestLogTauLogChi2->Draw("*");

  // show the L curve
  outputC.cd(6);
  lCurve->Draw("AL");
  bestLcurve->SetMarkerColor(kRed);
  bestLcurve->Draw("*");

  gROOT->GetListOfCanvases()->Draw();


  // Save the histograms
  //----------------------------------------------------------------------------
  output->cd();
  output->Write("", TObject::kOverwrite);
  output->Close();

}


//------------------------------------------------------------------------------
// PassTrigger
//
// https://github.com/latinos/PlotsConfigurations/blob/master/Configurations/ControlRegions/WW/Full2016/samples.py#L50-L56
//------------------------------------------------------------------------------
bool testUnfold::PassTrigger()
{

  bool _ismc = 1; 

  //  if (_ismc) 

  return true;

  /*
  if      (_sample.Contains("MuonEG"))         return ( trig_EleMu);
  else if (_sample.Contains("DoubleMuon"))     return (!trig_EleMu &&  trig_DbleMu);
  else if (_sample.Contains("SingleMuon"))     return (!trig_EleMu && !trig_DbleMu &&  trig_SnglMu);
  else if (_sample.Contains("DoubleEG"))       return (!trig_EleMu && !trig_DbleMu && !trig_SnglMu &&  trig_DbleEle);
  else if (_sample.Contains("SingleElectron")) return (!trig_EleMu && !trig_DbleMu && !trig_SnglMu && !trig_DbleEle && trig_SnglEle);
  else                                         return true;
  */

}
