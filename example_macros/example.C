#include <iostream>
#include <vector>
#include <algorithm>
#include "TCanvas.h"
#include "TError.h"
#include "TPad.h"
#include "TString.h"
#include "TRandom1.h"
#include "TH1F.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TString.h"
#include "TCut.h"
#include "TNtuple.h"
#include "TLine.h"
#include "/afs/cern.ch/user/d/dgulhan/workDir/missingPt/ntupler/trackTree.C"

void example(){
 TH1D::SetDefaultSumw2();
 
 //input file
 //TString directory="/mnt/hadoop/cms/store/user/yenjie/HiForest_v27/";
 // TString infname="Dijet100_HydjetDrum_v27_mergedV1";
 TString algo="akVs3Calo";
 TString directory="root://eoscms//eos/cms/store/group/phys_heavyions/dgulhan";
 TString infname="HiForest_Pythia_Hydjet_Jet80_Track8_Jet19_STARTHI53_LV1_merged_forest_0";
 
 trackTree * ftrk = new trackTree(Form("%s/%s.root",directory.Data(),infname.Data()));
 HiTree * fhi = new HiTree(Form("%s/%s.root",directory.Data(),infname.Data()));
 t * fjet = new t(Form("%s/%s.root",directory.Data(),infname.Data()),algo);
 
 //pt bins for track efficiency correction

 int npt=14; 
 double ptmin[]={0.4,0.4,0.4,0.4,0.4, 1, 1, 1,  1,  1, 3, 3,  3,  8};
 double ptmax[]={  1,  1,  1,  1,  1, 3, 3, 3,  3,  3, 8, 8,  8,300};
 
 int cent_min[]={  0, 10, 20, 30, 50, 0,10,20, 30, 50, 0,10, 20,  0};
 int cent_max[]={ 10, 20, 30, 50,100,10,20,30, 50,100,10,20,100,100};
 
 //getting histograms for track efficiency correction 
 TFile *f_eff[npt];
 TProfile *p_eff_cent[npt]; 
 TProfile2D *p_eff_accept[npt]; 
 TProfile *p_eff_pt[npt]; 
 TProfile *p_eff_rmin[npt];
 TFile *f_fake[npt];
 TProfile *p_fake_cent[npt]; 
 TProfile2D *p_fake_accept[npt]; 
 TProfile *p_fake_pt[npt]; 
 TProfile *p_fake_rmin[npt];  
 for(int ipt=0; ipt<npt;ipt++){
   f_eff[ipt]= new TFile(Form("../%s/eff/eff_pt%d_%d_cent%d_%d_%s_dogenjet0.root",algo.Data(),(int)ptmin[ipt],(int)ptmax[ipt],cent_min[ipt],cent_max[ipt],algo.Data()));
   p_eff_cent[ipt]=(TProfile*)f_eff[ipt]->Get("p_eff_cent");
   p_eff_pt[ipt]=(TProfile*)f_eff[ipt]->Get("p_eff_pt");
   p_eff_accept[ipt]=(TProfile2D*)f_eff[ipt]->Get("p_eff_acceptance");
   p_eff_rmin[ipt]=(TProfile*)f_eff[ipt]->Get("p_eff_rmin");
   f_fake[ipt]= new TFile(Form("../%s/fake/fake_pt%d_%d_cent%d_%d_%s_dogenjet0.root",algo.Data(),(int)ptmin[ipt],(int)ptmax[ipt],cent_min[ipt],cent_max[ipt],algo.Data()));
   p_fake_cent[ipt]=(TProfile*)f_fake[ipt]->Get("p_fake_cent");
   p_fake_pt[ipt]=(TProfile*)f_fake[ipt]->Get("p_fake_pt");
   p_fake_accept[ipt]=(TProfile2D*)f_fake[ipt]->Get("p_fake_acceptance");
   p_fake_rmin[ipt]=(TProfile*)f_fake[ipt]->Get("p_fake_rmin");
 }

 
 //output file and tree
 TFile *outf= new TFile(Form("track_ntuple_%s_%s_testforfake.root",infname.Data(),algo.Data()),"recreate");

 std::string trackVars="pt:eta:phi:rmin:trackselect:cent:eff:fake";
 TNtuple *nt_track = new TNtuple("nt_track","",trackVars.data());

 //loop over events
 int nentries = ftrk->GetEntriesFast();
 // for(int jentry=0;jentry<nentries;jentry++){
 for(int jentry=0;jentry<100;jentry++){
  if((jentry%1000)==0) std::cout<<jentry<<"/"<<nentries<<std::endl;

  ftrk->GetEntry(jentry);
  fhi->GetEntry(jentry);
  fjet->GetEntry(jentry);

  float cent=fhi->hiBin;
  //loop over tracks
  
  for(int itrk=0;itrk<ftrk->nTrk;itrk++){
   float trackselect=(ftrk->highPurity[itrk] && (ftrk->trkDxy1[itrk]/ftrk->trkDxyError1[itrk])<3.0 && (ftrk->trkDz1[itrk]/ftrk->trkDzError1[itrk])<3 && (ftrk->trkPtError[itrk]/ftrk->trkPt[itrk])<0.1);
   float eta=ftrk->trkEta[itrk];

   if(fabs(eta)>2.4) continue; //acceptance of the tracker   
   float pt=ftrk->trkPt[itrk];

   if(pt<0.5) continue; 
   float phi=ftrk->trkPhi[itrk];   
   
   float rmin=99;
   
   for(int ijet=0;ijet<fjet->ngen;ijet++){
    if(fabs(fjet->geneta[ijet])>2 || fjet->genpt[ijet]<30) continue;
    float r_reco=sqrt(pow(eta-fjet->geneta[ijet],2)+pow(acos(cos(phi-fjet->genphi[ijet])),2));
    if(r_reco<rmin)rmin=r_reco;
   }
   
   //get efficiency correction for the track   
   float fake_pt,fake_cent,fake_accept,fake_rmin;
   fake_pt=fake_cent=fake_accept=fake_rmin=0;
   
   float eff_pt,eff_cent,eff_accept,eff_rmin;
   eff_pt=eff_cent=eff_accept=eff_rmin=0;
   
   //given the pt,centrality,eta,phi and rmin of the track find the factorized efficiencies
   for(int ipt=0;ipt<npt;ipt++){
    if(pt>=ptmin[ipt] && pt<ptmax[ipt] && (0.5*cent)>=cent_min[ipt] && (0.5*cent)<cent_max[ipt]){
      eff_pt=p_eff_pt[ipt]->GetBinContent(p_eff_pt[ipt]->FindBin(pt));
      eff_cent=p_eff_cent[ipt]->GetBinContent(p_eff_cent[ipt]->FindBin(cent));
      eff_accept=p_eff_accept[ipt]->GetBinContent(p_eff_accept[ipt]->GetXaxis()->FindBin(phi),p_eff_accept[ipt]->GetYaxis()->FindBin(eta));
      if(rmin<5)eff_rmin=p_eff_rmin[ipt]->GetBinContent(p_eff_rmin[ipt]->FindBin(rmin));
      
      fake_pt=p_fake_pt[ipt]->GetBinContent(p_fake_pt[ipt]->FindBin(pt));
      fake_cent=p_fake_cent[ipt]->GetBinContent(p_fake_cent[ipt]->FindBin(cent));
      fake_accept=p_fake_accept[ipt]->GetBinContent(p_fake_accept[ipt]->GetXaxis()->FindBin(phi),p_fake_accept[ipt]->GetYaxis()->FindBin(eta));
      if(rmin<5)fake_rmin=p_fake_rmin[ipt]->GetBinContent(p_fake_rmin[ipt]->FindBin(rmin));      
     }     
   } 
   
   //multiply the factorized corrections to get the overall efficiency
   float eff=1;
   eff=eff_accept*eff_cent*eff_pt*eff_rmin;
   float fake=0;
   if(pt<100)fake=fake_accept+fake_cent+fake_pt+fake_rmin;   //the if statements are temporary next corrections won't have these
   if(eff==0){ //the if statements are temporary next corrections won't have these
    cout<<"zero efficiency"<<" eta="<<eta<<" pt="<<pt<<" phi="<<phi<<" cent="<<cent<<endl;
	  if(pt>100)eff=0.8;
	  else eff=1;
   }
   //fill in the output tree
   float entry[]={pt,eta,phi,rmin,trackselect,cent,eff,fake};
   nt_track->Fill(entry);
  }
   
 }
  cout<<"writing to file"<<endl;
  nt_track->Write();
  outf->Close();
}
