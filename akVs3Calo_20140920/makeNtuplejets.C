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
#include "/data/dgulhan/missingPt/2014_07_03/ntupler/trackTree.C"

void makeNtuplejets(){
 TH1D::SetDefaultSumw2();

//all statistics 
 float pthatWeight[5] = {4.29284e-01,2.99974e-02,3.38946e-4,1.06172e-4,2.79631e-5};
//100k events
// float pthatWeight[5] = {0.429284,0.0299974,0.000949812,0.000232709,7.61038e-05};
 //float vertexShift = 0.501501;
//100k events high pthat

 // float pthatWeight[5] = {0,0,0.000654317,0.000156607,5.07966e-05};
 float vertexShift = 0.406408;

 const int nevents = 50000;

 TString directory="root://eoscms//eos/cms/store/group/phys_heavyions/dgulhan/HIMC/Track8_Jet24_FixedJEC/";
 const char* infname[5];
 infname[0] = "/HydjetDrum_Pyquen_Dijet30_FOREST_Track8_Jet24_FixedPtHatJES_v0";
 infname[1] = "/HydjetDrum_Pyquen_Dijet50_FOREST_Track8_Jet24_FixedPtHatJES_v0";
 infname[2] = "/HydjetDrum_Pyquen_Dijet80_FOREST_Track8_Jet24_FixedPtHatJES_v0";
 infname[3] = "/HydjetDrum_Pyquen_Dijet100_FOREST_Track8_Jet24_FixedPtHatJES_v0";
 infname[4] = "/HydjetDrum_Pyquen_Dijet120_FOREST_Track8_Jet24_FixedPtHatJES_v0";

 trackTree * ftrk[5];
 HiTree * fhi[5];
 HltTree * fhlt[5];
 t * fjet[5];
 TFile * evtSelFile[5];
 TTree * evtSel[5];
 int pcoll[5];
 
 
 const int n_rmin_bins=36;
 double rmin_bins[n_rmin_bins+1] = {0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.,2.2,2.4,2.6,2.8,3,100};
 TH1D * hrmin_gen = new TH1D("hrmin_gen","",n_rmin_bins,rmin_bins);
 TH1D * hrmin_reco = new TH1D("hrmin_reco","",n_rmin_bins,rmin_bins);
 
 for(int ifile=0; ifile<5; ifile++){
   ftrk[ifile] = new trackTree(Form("%s/%s.root",directory.Data(),infname[ifile]));
   fhi[ifile] = new HiTree(Form("%s/%s.root",directory.Data(),infname[ifile]));
   fjet[ifile] = new t(Form("%s/%s.root",directory.Data(),infname[ifile]),"akVs3Calo");
   fhlt[ifile] = new HltTree(Form("%s/%s.root",directory.Data(),infname[ifile]));
   // evtSelFile[ifile] = new TFile(Form("%s/%s.root",directory.Data(),infname[ifile]),"read");
   // evtSel[ifile] = (TTree*) evtSelFile[ifile]->Get("skimanalysis/HltTree");
   // evtSel[ifile]->SetBranchAddress("pcollisionEventSelection", &pcoll[ifile]);
  }

 TFile * centWeightsFile = new TFile("centrality_weights.root","read");
 TH1F * centWeights = new TH1F("centWeight","centWeight",100,0,200);
 centWeights = (TH1F*)centWeightsFile->Get("centrality_weight");
 
 //pt bins for track efficiency correction Yen-Jie
 int npt_fake=29; 
 double ptmin_fake[]={0.5 ,0.5 ,0.5 ,0.5 ,0.5 ,0.55 ,0.55 ,0.55 ,0.55 ,0.55 ,0.65,0.65,0.65,0.65,0.65,0.8,0.8,0.8,0.8,0.8,1 ,1 ,1 ,1 ,1  ,3 ,3 ,3  ,8};
 double ptmax_fake[]={0.55,0.55,0.55,0.55,0.55,0.65 ,0.65 ,0.65 ,0.65 ,0.65 ,0.8 ,0.8 ,0.8 ,0.8 ,0.8 ,1  ,1  ,1  ,1  ,1  ,3 ,3 ,3 ,3 ,3  ,8 ,8 ,8  ,300};
 
 int cent_min_fake[]={0   ,20  ,40  ,60  ,100  ,0    ,20   ,40   ,60   ,100   ,0   ,20  ,40  ,60  ,100  ,0  ,20 ,40 ,60 ,100 ,0 ,20,40,60,100 ,0 ,20,40 ,0};
 int cent_max_fake[]={20  ,40  ,60  ,100  ,200 ,20   ,40   ,60   ,100   ,200  ,20  ,40  ,60  ,100  ,200 ,20 ,40 ,60 ,100 ,200,20,40,60,100,200,20,40,200,200};
 
 int npt_eff=29; 
 double ptmin_eff[]={0.5 ,0.5 ,0.5 ,0.5 ,0.5 ,0.55 ,0.55 ,0.55 ,0.55 ,0.55 ,0.65,0.65,0.65,0.65,0.65,0.8,0.8,0.8,0.8,0.8,1 ,1 ,1 ,1 ,1  ,3 ,3 ,3  ,8};
 double ptmax_eff[]={0.55,0.55,0.55,0.55,0.55,0.65 ,0.65 ,0.65 ,0.65 ,0.65 ,0.8 ,0.8 ,0.8 ,0.8 ,0.8 ,1  ,1  ,1  ,1  ,1  ,3 ,3 ,3 ,3 ,3  ,8 ,8 ,8  ,300};
 
 int cent_min[]={0   ,20  ,40  ,60  ,100  ,0    ,20   ,40   ,60   ,100   ,0   ,20  ,40  ,60  ,100  ,0  ,20 ,40 ,60 ,100 ,0 ,20,40,60,100 ,0 ,20,40 ,0};
 int cent_max[]={20  ,40  ,60  ,100  ,200 ,20   ,40   ,60   ,100   ,200  ,20  ,40  ,60  ,100  ,200 ,20 ,40 ,60 ,100 ,200,20,40,60,100,200,20,40,200,200};
  cout<<0<<endl;

 //getting histograms for track efficiency correction 
 TFile *f_eff[npt_eff];
 TProfile *p_eff_cent[npt_eff]; 
 TProfile2D *p_eff_accept[npt_eff]; 
 TProfile *p_eff_pt[npt_eff]; 
 TProfile *p_eff_rmin[npt_eff]; 
 for(int ipt=0; ipt<npt_eff;ipt++){
   // f_eff[ipt]= new TFile(Form("/afs/cern.ch/user/d/dgulhan/workDir/trackCorrection/akVs3Calo_04302014/eff_pt%d_%d_cent%d_%d.root",(int)(ptmin_eff[ipt]*100),(int)(ptmax_eff[ipt]*100),(int)(0.5*cent_min[ipt]),(int)(0.5*cent_max[ipt])));
   f_eff[ipt]= new TFile(Form("/afs/cern.ch/user/d/dgulhan/workDir/JetTrack/TrackCorrectionTables/akVs3Calo_fineptbin/eff/eff_pt%d_%d_cent%d_%d.root",(int)(ptmin_eff[ipt]*100),(int)(ptmax_eff[ipt]*100),(int)(0.5*cent_min[ipt]),(int)(0.5*cent_max[ipt])));
   p_eff_cent[ipt]=(TProfile*)f_eff[ipt]->Get("p_eff_cent");
   p_eff_pt[ipt]=(TProfile*)f_eff[ipt]->Get("p_eff_pt");
   p_eff_accept[ipt]=(TProfile2D*)f_eff[ipt]->Get("p_eff_acceptance");
   p_eff_rmin[ipt]=(TProfile*)f_eff[ipt]->Get("p_eff_rmin");
 }

 TFile *f_fake[npt_fake];
 TProfile *p_fake_cent[npt_fake]; 
 TProfile2D *p_fake_accept[npt_fake]; 
 TProfile *p_fake_pt[npt_fake]; 
 TProfile *p_fake_rmin[npt_fake]; 
 for(int ipt=0; ipt<npt_fake;ipt++){
   // f_fake[ipt]= new TFile(Form("/afs/cern.ch/user/d/dgulhan/workDir/trackCorrection/akVs3Calo_04302014/fake_pt%d_%d_cent%d_%d.root",(int)(ptmin_fake[ipt]*100),(int)(ptmax_fake[ipt]*100),(int)(0.5*cent_min[ipt]),(int)(0.5*cent_max[ipt])));
   f_fake[ipt]= new TFile(Form("/afs/cern.ch/user/d/dgulhan/workDir/JetTrack/TrackCorrectionTables/akVs3Calo_fineptbin/fake/fake_pt%d_%d_cent%d_%d.root",(int)(ptmin_fake[ipt]*100),(int)(ptmax_fake[ipt]*100),(int)(0.5*cent_min[ipt]),(int)(0.5*cent_max[ipt])));
   p_fake_cent[ipt]=(TProfile*)f_fake[ipt]->Get("p_fake_cent");
   p_fake_pt[ipt]=(TProfile*)f_fake[ipt]->Get("p_fake_pt");
   p_fake_accept[ipt]=(TProfile2D*)f_fake[ipt]->Get("p_fake_acceptance");
   p_fake_rmin[ipt]=(TProfile*)f_fake[ipt]->Get("p_fake_rmin");
 }

 //output file and tree
 TFile *outf= new TFile("git_nojetcut_track_ntuple_pthatCombo_50k_jets_fixedeta_8_09_2014.root","recreate");
 
 std::string particleVars="pt:matchedpt:eta:phi:rmin:trackselect:cent:eff:cent_weight:pthat_weight:weight:pt1:pt2:dphi:asym:eta1:eta2:jteta_rmin:isleadrmin:issubleadrmin";
 TNtuple *nt_particle = new TNtuple("nt_particle","",particleVars.data());
 
 std::string trackVars="pt:eta:phi:rmin:trackselect:trackstatus:cent:eff:trkfake:fake:cent_weight:pthat_weight:weight:pt1:pt2:dphi:asym:eta1:eta2:jteta_rmin:isleadrmin:issubleadrmin";
 TNtuple *nt_track = new TNtuple("nt_track","",trackVars.data());


 //loop over events

 for(int ifile=0; ifile<5; ifile++){
 std::cout<<ifile<<std::endl;
 int nentries = ftrk[ifile]->GetEntriesFast();
 if(nevents<nentries) nentries = nevents; 
for(int jentry=0;jentry<nentries;jentry++){
 //for(int jentry=0;jentry<5000;jentry++){
  if((jentry%1000)==0) std::cout<<jentry<<"/"<<nentries<<std::endl;
  ftrk[ifile]->GetEntry(jentry);
  fhi[ifile]->GetEntry(jentry);
  fjet[ifile]->GetEntry(jentry);
  fhlt[ifile]->GetEntry(jentry);
  // evtSel[ifile]->GetEntry(jentry);

  float cent=fhi[ifile]->hiBin;
  float vz = fhi[ifile]->vz;
 
  if(fabs(vz-vertexShift)>15 || !(fhlt[ifile]->pcollisionEventSelection)) continue;

  float weight = 0;
  float pthat_weight = 0;
  float cent_weight = 0;

  if(fjet[ifile]->pthat <50)       pthat_weight = pthatWeight[0];
  else if(fjet[ifile]->pthat <80)  pthat_weight = pthatWeight[1];
  else if(fjet[ifile]->pthat <100) pthat_weight = pthatWeight[2];
  else if(fjet[ifile]->pthat <120) pthat_weight = pthatWeight[3];
  else                             pthat_weight = pthatWeight[4];

  cent_weight = centWeights->GetBinContent(centWeights->FindBin(cent));
  weight = pthat_weight*cent_weight;


  float pt1=-99;
  float phi1=-99;
  float eta1=-99;
  float refpt1=-99;
  float refeta1=-99;
  float refphi1=-99;
  float matchedpt1=-99;
  float matchedR1=-99;
  float trackMax1=-99;
  float pt2=-99;
  float phi2=-99;
  float eta2=-99;
  float refpt2=-99;
  float refphi2=-99;
  float refeta2=-99;
  float matchedpt2=-99;
  float matchedR2=-99;
  float trackMax2=-99;
  float pt3=-99;
  float phi3=-99;
  float eta3=-99;
  float refpt3=-99;
  float refeta3=-99;
  float refphi3=-99;
  float matchedpt3=-99;
  float matchedR3=-99;
  float trackMax3=-99;
  float dphi=-99;
  float ptratio=-99;
  float asym = -1;

std::vector<std::pair<double, std::pair<double,std::pair<double, std::pair<double,std::pair<double,std::pair<double,std::pair<double,std::pair<double,double> > > > > > > > > jets;
  int njet=0;
  for(int ijet=0;ijet<fjet[ifile]->nref;ijet++){

   if(fabs(fjet[ifile]->jteta[ijet])>2) continue;
   jets.push_back(std::make_pair(fjet[ifile]->jtpt[ijet],std::make_pair(fjet[ifile]->jteta[ijet], std::make_pair(fjet[ifile]->jtphi[ijet], std::make_pair(fjet[ifile]->refpt[ijet],std::make_pair(fjet[ifile]->refeta[ijet],std::make_pair(fjet[ifile]->refphi[ijet],std::make_pair(fjet[ifile]->matchedPt[ijet],std::make_pair(fjet[ifile]->matchedR[ijet],fjet[ifile]->trackMax[ijet])))))))));
   njet++;
  }

  std::sort(jets.begin(),jets.end());
  if(njet>0){
   pt1=       jets[njet-1].first;
   eta1=      jets[njet-1].second.first;
   phi1=      jets[njet-1].second.second.first;
   refpt1=    jets[njet-1].second.second.second.first;
   refeta1=   jets[njet-1].second.second.second.second.first;
   refphi1=   jets[njet-1].second.second.second.second.second.first;
   matchedpt1=jets[njet-1].second.second.second.second.second.second.first;
   matchedR1= jets[njet-1].second.second.second.second.second.second.second.first;
   trackMax1= jets[njet-1].second.second.second.second.second.second.second.second;
if(njet>1){
    pt2=jets[njet-2].first;
    eta2=jets[njet-2].second.first;
    phi2=jets[njet-2].second.second.first;
    refpt2=jets[njet-2].second.second.second.first;
    refeta2=jets[njet-2].second.second.second.second.first;
    refphi2=jets[njet-2].second.second.second.second.second.first;
    matchedpt2=jets[njet-2].second.second.second.second.second.second.first;
    matchedR2=jets[njet-2].second.second.second.second.second.second.second.first;
    trackMax2=jets[njet-2].second.second.second.second.second.second.second.second;
    dphi=acos(cos(phi1-phi2));
    ptratio=pt2/pt1;
    asym = (pt1-pt2)/(pt1+pt2);
    if(njet>2){
     pt3=jets[njet-3].first;
     eta3=jets[njet-3].second.first;
     phi3=jets[njet-3].second.second.first;
     refpt3=jets[njet-3].second.second.second.first;
     refeta3=jets[njet-3].second.second.second.second.first;
     refphi3=jets[njet-3].second.second.second.second.second.first;
     matchedpt3=jets[njet-3].second.second.second.second.second.second.first;
     matchedR3=jets[njet-3].second.second.second.second.second.second.second.first;
     trackMax3=jets[njet-3].second.second.second.second.second.second.second.second;
    }
   }
  }

  // if(pt1<120 || pt2<50 ||  fabs(eta1)>0.5 || fabs(eta2)>0.5 || dphi<(5*TMath::Pi()/6)) continue; 
  //loop over tracks
  for(int itrk=0;itrk<ftrk[ifile]->nParticle;itrk++){

   float trackselect=(ftrk[ifile]->mtrkQual[itrk] && fabs(ftrk[ifile]->mtrkDxy1[itrk]/ftrk[ifile]->mtrkDxyError1[itrk])<3.0 && fabs(ftrk[ifile]->mtrkDz1[itrk]/ftrk[ifile]->mtrkDzError1[itrk])<3 && (ftrk[ifile]->mtrkPtError[itrk]/ftrk[ifile]->mtrkPt[itrk])<0.1);
   float eta=ftrk[ifile]->pEta[itrk];

   if(fabs(eta)>2.4) continue; //acceptance of the tracker
   float pt=ftrk[ifile]->pPt[itrk];
   float mpt=ftrk[ifile]->mtrkPt[itrk];
   float phi=ftrk[ifile]->pPhi[itrk];
   float rmin=199;
   float isleadrmin=0;
   float issubleadrmin=0;
   float jteta_rmin=-999;
 
   for(int ijet=0;ijet<fjet[ifile]->nref;ijet++){
     if(fabs(fjet[ifile]->jteta[ijet])>2 || fjet[ifile]->jtpt[ijet]<50) continue;
     float r_reco=sqrt(pow(eta-fjet[ifile]->jteta[ijet],2)+pow(acos(cos(phi-fjet[ifile]->jtphi[ijet])),2));
     if(r_reco<rmin){
	  rmin=r_reco;
	  jteta_rmin=fjet[ifile]->jteta[ifile];
	  if(pt1==fjet[ifile]->jtpt[ifile])isleadrmin=1;
	  if(pt2==fjet[ifile]->jtpt[ifile])issubleadrmin=1;
	 }
    }

   
   //get efficiency correction for the track
   float eff_accept=1;
   float eff_pt=1;
   float eff_cent=1;
   float eff_rmin=1;
   for(int ipt=0;ipt<npt_eff;ipt++){
    if(pt>=ptmin_eff[ipt] && pt<ptmax_eff[ipt] && cent>=cent_min[ipt] && cent<cent_max[ipt]){
      eff_pt=p_eff_pt[ipt]->GetBinContent(p_eff_pt[ipt]->FindBin(pt));
      eff_cent=p_eff_cent[ipt]->GetBinContent(p_eff_cent[ipt]->FindBin(cent));
      eff_accept=p_eff_accept[ipt]->GetBinContent(p_eff_accept[ipt]->GetXaxis()->FindBin(phi),p_eff_accept[ipt]->GetYaxis()->FindBin(eta));
      if(rmin<=100) eff_rmin=p_eff_rmin[ipt]->GetBinContent(p_eff_rmin[ipt]->FindBin(rmin));
     }     
   } 

   float eff=eff_accept*eff_cent*eff_pt*eff_rmin;
   
   //fill in the output tree
  
   float entry[]={pt,mpt,eta,phi,rmin,trackselect,cent,eff,cent_weight,pthat_weight,weight,pt1,pt2,dphi,asym,eta1,eta2,jteta_rmin,isleadrmin,issubleadrmin};

   nt_particle->Fill(entry);
   hrmin_gen->Fill(rmin);
  }
  
  for(int itrk=0;itrk<ftrk[ifile]->nTrk;itrk++){

   float trackselect=(ftrk[ifile]->highPurity[itrk] && fabs(ftrk[ifile]->trkDxy1[itrk]/ftrk[ifile]->trkDxyError1[itrk])<3.0 && fabs(ftrk[ifile]->trkDz1[itrk]/ftrk[ifile]->trkDzError1[itrk])<3 && (ftrk[ifile]->trkPtError[itrk]/ftrk[ifile]->trkPt[itrk])<0.1);
   float eta=ftrk[ifile]->trkEta[itrk];

   if(fabs(eta)>2.4) continue; //acceptance of the tracker   
   
   float pt=ftrk[ifile]->trkPt[itrk];
   float phi=ftrk[ifile]->trkPhi[itrk];
   float trkfake=ftrk[ifile]->trkFake[itrk];
   float trackstatus=ftrk[ifile]->trkStatus[itrk];
   float rmin=199;
   
   float jteta_rmin,isleadrmin,issubleadrmin;
   jteta_rmin=-999;
   isleadrmin=issubleadrmin=0;

   //find rmin; 
     for(int ijet=0;ijet<fjet[ifile]->nref;ijet++){
     if(fabs(fjet[ifile]->jteta[ijet])>2 || fjet[ifile]->jtpt[ijet]<50) continue;
     float r_reco=sqrt(pow(eta-fjet[ifile]->jteta[ijet],2)+pow(acos(cos(phi-fjet[ifile]->jtphi[ijet])),2));
     if(r_reco<rmin){
	  rmin=r_reco;
	  jteta_rmin=fjet[ifile]->jteta[ifile];
	  if(pt1==fjet[ifile]->jtpt[ifile])isleadrmin=1;
	  if(pt2==fjet[ifile]->jtpt[ifile])issubleadrmin=1;
	 }
    }

   //get efficiency and fake rate correction for the track Yen-Jie
   float eff_accept=1;
   float eff_pt=1;
   float eff_cent=1;
   float eff_rmin=1;

   float fake_pt,fake_cent,fake_accept,fake_rmin;
   fake_pt=fake_cent=fake_accept=fake_rmin=0;

   for(int ipt=0;ipt<npt_eff;ipt++){
    if(pt>=ptmin_eff[ipt] && pt<ptmax_eff[ipt] && cent>=cent_min[ipt] && cent<cent_max[ipt]){
      eff_pt=p_eff_pt[ipt]->GetBinContent(p_eff_pt[ipt]->FindBin(pt));
      eff_cent=p_eff_cent[ipt]->GetBinContent(p_eff_cent[ipt]->FindBin(cent));
      eff_accept=p_eff_accept[ipt]->GetBinContent(p_eff_accept[ipt]->GetXaxis()->FindBin(phi),p_eff_accept[ipt]->GetYaxis()->FindBin(eta));
      if(rmin<=100) eff_rmin=p_eff_rmin[ipt]->GetBinContent(p_eff_rmin[ipt]->FindBin(rmin)); 
     }     
   } 
   
   for(int ipt=0;ipt<npt_fake;ipt++){
    if(pt>=ptmin_fake[ipt] && pt<ptmax_fake[ipt] && cent>=cent_min_fake[ipt] && cent<cent_max_fake[ipt]){
      fake_pt=p_fake_pt[ipt]->GetBinContent(p_fake_pt[ipt]->FindBin(pt));
      fake_cent=p_fake_cent[ipt]->GetBinContent(p_fake_cent[ipt]->FindBin(cent));
      fake_accept=p_fake_accept[ipt]->GetBinContent(p_fake_accept[ipt]->GetXaxis()->FindBin(phi),p_fake_accept[ipt]->GetYaxis()->FindBin(eta));
      if(rmin<=100) fake_rmin=p_fake_rmin[ipt]->GetBinContent(p_fake_rmin[ipt]->FindBin(rmin));
     }     
   }

   float eff=1;
   eff=eff_accept*eff_cent*eff_pt*eff_rmin;
   float fake=0;
   if(pt<100)fake=fake_accept+fake_cent+fake_pt+fake_rmin;
   if(eff==0){
          //cout<<"zero efficiency"<<" eta="<<eta<<" pt="<<pt<<" phi="<<phi<<" cent="<<cent<<endl;
	  if(pt>100)eff=0.8;
	  else eff=1;
   }

   //if(fake<0) fake=0;

   //fill in the output tree
   float entry[]={pt,eta,phi,rmin,trackselect,trackstatus,cent,eff,trkfake,fake,cent_weight,pthat_weight,weight,pt1,pt2,dphi,asym,eta1,eta2,jteta_rmin,isleadrmin,issubleadrmin};
   nt_track->Fill(entry); 
   hrmin_reco->Fill(rmin);
  }
 }
}
 
  //nt_track->Write(); 
 // nt_particle->Write();
    outf->Write();
    outf->Close();
    
    TFile *outf2=new TFile("weight_rmin.root","recreate");
    hrmin_reco->Write();
    hrmin_gen->Write();
    TH1D *hclone=(TH1D*)hrmin_gen->Clone("hclone");
    hclone->Divide(hrmin_reco);
    hclone->Write();
    outf2->Close();
}
