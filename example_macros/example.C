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
#include "../../../missingPt/ntupler/trackTree.C"
 
void example(){
 TH1D::SetDefaultSumw2();
 
 TString algo="akVs3Calo";
 
 TString directory="root://eoscms//eos/cms//store/group/phys_heavyions/dgulhan/HIMC/Track8_Jet24_FixedJEC/";
 TString infname="HydjetDrum_Pyquen_Dijet80_FOREST_Track8_Jet24_FixedPtHatJES_v0";
 trackTree * ftrk = new trackTree(Form("%s/%s.root",directory.Data(),infname.Data()));
 HiTree * fhi = new HiTree(Form("%s/%s.root",directory.Data(),infname.Data()));
 t * fjet = new t(Form("%s/%s.root",directory.Data(),infname.Data()),algo.Data());
 
 //pt bins for track efficiency correction
  int npt=29; 
 double ptmin[]={0.5,0.5,0.5,0.5,0.5,0.55,0.55,0.55,0.55,0.55,0.65,0.65,0.65,0.65,0.65,0.8,0.8,0.8,0.8,0.8, 1, 1, 1,  1,  1, 3, 3,  3,  8};
 double ptmax[]={ 0.55,0.55,0.55,0.55,0.55,0.65,0.65,0.65,0.65,0.65,0.8,0.8,0.8,0.8,0.8, 1,  1,  1,  1,  1, 3, 3, 3,  3,  3, 8, 8,  8,300};
 
 int cent_min[]={  0, 20, 40, 60,100,0, 20, 40, 60,100,0, 20, 40, 60,100,0, 20, 40, 60,100, 0,20,40, 60,100, 0,20, 40,  0};
 int cent_max[]={ 20, 40, 60,100,200,20, 40, 60,100,200,20, 40, 60,100,200,20, 40, 60,100,200,20,40,60,100,200,20,40,200,200};
 
   cout<<0<<endl;

 //getting histograms for track efficiency correction 
  
 TFile *f_eff[npt];
 TProfile *p_eff_cent[npt]; 
 TProfile2D *p_eff_accept[npt]; 
 TProfile *p_eff_pt[npt]; 
 TProfile *p_eff_rmin[npt]; 
 for(int ipt=0; ipt<npt;ipt++){
   f_eff[ipt]= new TFile(Form("/afs/cern.ch/user/d/dgulhan/workDir/JetTrack/TrackCorrectionTables/akVs3Calo_fineptbin/eff/eff_pt%d_%d_cent%d_%d.root",(int)(ptmin[ipt]*100),(int)(ptmax[ipt]*100),(int)(0.5*cent_min[ipt]),(int)(0.5*cent_max[ipt])));
   p_eff_cent[ipt]=(TProfile*)f_eff[ipt]->Get("p_eff_cent");
   p_eff_pt[ipt]=(TProfile*)f_eff[ipt]->Get("p_eff_pt");
   p_eff_accept[ipt]=(TProfile2D*)f_eff[ipt]->Get("p_eff_acceptance");
   p_eff_rmin[ipt]=(TProfile*)f_eff[ipt]->Get("p_eff_rmin");
 }

 TFile *f_fake[npt];
 TProfile *p_fake_cent[npt]; 
 TProfile2D *p_fake_accept[npt]; 
 TProfile *p_fake_pt[npt]; 
 TProfile *p_fake_rmin[npt]; 
 for(int ipt=0; ipt<npt;ipt++){
   f_fake[ipt]= new TFile(Form("/afs/cern.ch/user/d/dgulhan/workDir/JetTrack/TrackCorrectionTables/akVs3Calo_fineptbin/fake/fake_pt%d_%d_cent%d_%d.root",(int)(ptmin[ipt]*100),(int)(ptmax[ipt]*100),(int)(0.5*cent_min[ipt]),(int)(0.5*cent_max[ipt]),algo.Data()));
   p_fake_cent[ipt]=(TProfile*)f_fake[ipt]->Get("p_fake_cent");
   p_fake_pt[ipt]=(TProfile*)f_fake[ipt]->Get("p_fake_pt");
   p_fake_accept[ipt]=(TProfile2D*)f_fake[ipt]->Get("p_fake_acceptance");
   p_fake_rmin[ipt]=(TProfile*)f_fake[ipt]->Get("p_fake_rmin");
 }
 
 //output file and tree
 TFile *outf= new TFile(Form("MB.root",infname.Data(),algo.Data()),"recreate");
 std::string particleVars="pt:matchedpt:eta:phi:rmin:trackselect:cent:eff";

 TNtuple *nt_particle = new TNtuple("nt_particle","",particleVars.data());
 std::string trackVars="pt:eta:phi:rmin:trackselect:trackstatus:cent:eff:trkfake:fake";

 TNtuple *nt_track = new TNtuple("nt_track","",trackVars.data());

 //loop over events
 int nentries = ftrk->GetEntriesFast();
 // for(int jentry=0;jentry<nentries;jentry++){
 for(int jentry=0;jentry<1000;jentry++){
  if((jentry%1000)==0) std::cout<<jentry<<"/"<<nentries<<std::endl;
  ftrk->GetEntry(jentry);
  fhi->GetEntry(jentry);
  fjet->GetEntry(jentry);

  float cent=fhi->hiBin;
  //loop over tracks
  for(int itrk=0;itrk<ftrk->nParticle;itrk++){

   float trackselect=(ftrk->mtrkQual[itrk] && fabs(ftrk->mtrkDxy1[itrk]/ftrk->mtrkDxyError1[itrk])<3.0 && fabs(ftrk->mtrkDz1[itrk]/ftrk->mtrkDzError1[itrk])<3 && (ftrk->mtrkPtError[itrk]/ftrk->mtrkPt[itrk])<0.1);
   float eta=ftrk->pEta[itrk];
   if(fabs(eta)>2.4) continue; 
   float pt=ftrk->pPt[itrk];
   float mpt=ftrk->mtrkPt[itrk];
   float phi=ftrk->pPhi[itrk];
   float rmin=99;
 

   for(int ijet=0;ijet<fjet->nref;ijet++){
    if(fabs(fjet->jteta[ijet])>2 || fjet->jtpt[ijet]<30) continue;
    float r_reco=sqrt(pow(eta-fjet->jteta[ijet],2)+pow(acos(cos(phi-fjet->jtphi[ijet])),2));
    if(r_reco<rmin)rmin=r_reco;
   }
   
   //get efficiency correction for the track
   float eff_accept=1;
   float eff_pt=1;
   float eff_cent=1;
   float eff_rmin=1;

   for(int ipt=0;ipt<npt;ipt++){
    if(pt>=ptmin[ipt] && pt<ptmax[ipt] && cent>=cent_min[ipt] && cent<cent_max[ipt]){
      eff_pt=p_eff_pt[ipt]->GetBinContent(p_eff_pt[ipt]->FindBin(pt));
      eff_cent=p_eff_cent[ipt]->GetBinContent(p_eff_cent[ipt]->FindBin(cent));
      eff_accept=p_eff_accept[ipt]->GetBinContent(p_eff_accept[ipt]->GetXaxis()->FindBin(phi),p_eff_accept[ipt]->GetYaxis()->FindBin(eta));
      if(rmin<5)eff_rmin=p_eff_rmin[ipt]->GetBinContent(p_eff_rmin[ipt]->FindBin(rmin));
     }     
   } 
   
   float eff=1;
   eff=eff_accept*eff_cent*eff_pt*eff_rmin;
   if(eff==0){
    cout<<"zero efficiency"<<" eta="<<eta<<" pt="<<pt<<" phi="<<phi<<" cent="<<cent<<" rmin="<<rmin<<" eff_cent="<<eff_cent<<"eff_pt"<<eff_pt<<" eff_accept"<<eff_accept<<" eff_rmin"<<eff_rmin<<endl;
	  if(pt>100)eff=0.8;
	  else eff=1;
   }
   // if(eff_cent!=0) cout<<"cent= "<<cent<<endl;
   //fill in the output tree

   float entry[]={pt,mpt,eta,phi,rmin,trackselect,cent,eff};

   nt_particle->Fill(entry);

  }
  
  for(int itrk=0;itrk<ftrk->nTrk;itrk++){

   float trackselect=(ftrk->highPurity[itrk] && fabs(ftrk->trkDxy1[itrk]/ftrk->trkDxyError1[itrk])<3.0 && fabs(ftrk->trkDz1[itrk]/ftrk->trkDzError1[itrk])<3 && (ftrk->trkPtError[itrk]/ftrk->trkPt[itrk])<0.1);
   float eta=ftrk->trkEta[itrk];

   if(fabs(eta)>2.4) continue;   
   float pt=ftrk->trkPt[itrk];

   if(pt<0.5) continue; 
   float phi=ftrk->trkPhi[itrk];
   float trkfake=ftrk->trkFake[itrk];
   float trackstatus=ftrk->trkStatus[itrk];
   float rmin=99;

   for(int ijet=0;ijet<fjet->nref;ijet++){
    if(fabs(fjet->jteta[ijet])>2 || fjet->jtpt[ijet]<30) continue;
    float r_reco=sqrt(pow(eta-fjet->jteta[ijet],2)+pow(acos(cos(phi-fjet->jtphi[ijet])),2));
    if(r_reco<rmin)rmin=r_reco;
   }

   //get efficiency correction for the track
    float eff_accept=1;
    float eff_pt=1;
    float eff_cent=1;
    float eff_rmin=1;
   
    float fake_pt,fake_cent,fake_accept,fake_rmin;
    fake_pt=fake_cent=fake_accept=fake_rmin=0;

   for(int ipt=0;ipt<npt;ipt++){
    if(pt>=ptmin[ipt] && pt<ptmax[ipt] && cent>=cent_min[ipt] && cent<cent_max[ipt]){
      eff_pt=p_eff_pt[ipt]->GetBinContent(p_eff_pt[ipt]->FindBin(pt));
      eff_cent=p_eff_cent[ipt]->GetBinContent(p_eff_cent[ipt]->FindBin(cent));
      eff_accept=p_eff_accept[ipt]->GetBinContent(p_eff_accept[ipt]->GetXaxis()->FindBin(phi),p_eff_accept[ipt]->GetYaxis()->FindBin(eta));
      if(rmin<5)eff_rmin=p_eff_rmin[ipt]->GetBinContent(p_eff_rmin[ipt]->FindBin(rmin));
     }     
   } 
   
   for(int ipt=0;ipt<npt;ipt++){
    if(pt>=ptmin[ipt] && pt<ptmax[ipt] && cent>=cent_min[ipt] && cent<cent_max[ipt]){
      fake_pt=p_fake_pt[ipt]->GetBinContent(p_fake_pt[ipt]->FindBin(pt));
      fake_cent=p_fake_cent[ipt]->GetBinContent(p_fake_cent[ipt]->FindBin(cent));
      fake_accept=p_fake_accept[ipt]->GetBinContent(p_fake_accept[ipt]->GetXaxis()->FindBin(phi),p_fake_accept[ipt]->GetYaxis()->FindBin(eta));
      if(rmin<5) fake_rmin=p_fake_rmin[ipt]->GetBinContent(p_fake_rmin[ipt]->FindBin(rmin));
     }     
   }
  
   float eff=1;
   eff=eff_accept*eff_cent*eff_pt*eff_rmin;
   if(eff==0){
    // cout<<"zero efficiency"<<" eta="<<eta<<" pt="<<pt<<" phi="<<phi<<" cent="<<cent<<endl;
    if(pt>100)eff=0.8;
	  else eff=1;
   }
   float fake=0;
   fake=fake_accept+fake_cent+fake_pt+fake_rmin;
   
   // float eff=1;
   // eff=eff_accept*eff_cent*eff_pt*eff_rmin;
   // float fake=0;
   // if(pt<100)fake=fake_accept+fake_cent+fake_pt+fake_rmin;
   // if(eff==0){
	  // if(pt>100)eff=0.8;
	  // else eff=1;
   // }

   //fill in the output tree
   float entry[]={pt,eta,phi,rmin,trackselect,trackstatus,cent,eff,trkfake,fake};
   nt_track->Fill(entry);

  }
   
 }
 
  nt_track->Write();
  nt_particle->Write();
  outf->Close();
}
