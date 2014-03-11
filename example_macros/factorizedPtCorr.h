

if(pt < .5){
  if(eff)
    return 1;
  else
    return 0;
 }

if(pt > 100){
  if(eff)
    return .8;
  else
    return 0;
 }

corrFactor = corrFactor*(centProf_p->GetBinContent(hiBin+1));
corrFactor = corrFactor*(etaPhiProf_p->GetBinContent(etaPhiProf_p->FindBin(phi, eta)));
corrFactor = corrFactor*(ptProf_p->GetBinContent(ptProf_p->FindBin(pt)));

if(rmin < rMinCut)
  corrFactor = corrFactor*(rminProf_p->GetBinContent(rminProf_p->FindBin(rmin)));

return corrFactor;
}



/*
  Above implemented in code as follows:

      Int_t hiBinDiv[5] = {20, 40, 60, 100, 200};
          Int_t hiSet[15] = {0, 5, 10, 1, 6, 11, 2, 7, 12, 3, 8, 12, 4, 9, 12};

	      for(Int_t hiBinIter = 0; hiBinIter < 5; hiBinIter++){
	            if(hiBin_ < hiBinDiv[hiBinIter]){
		            for(Int_t trkEntry = 0; trkEntry < nTrk_; trkEntry++){
			              Int_t ptPos = getPtBin(trkPt_[trkEntry], hiSet[hiBinIter*3], hiSet[hiBinIter*3 + 1], hiSet[hiBinIter*3 + 2]);

				                trkPtCorrPF_[trkEntry] = trkPt_[trkEntry]/factorizedPtCorr(hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], trkRMinPF_[trkEntry], PFcent_p[ptPos], PFphiEta_p[ptPos], PFpt_p[ptPos], PFdelR_p[ptPos], 3);
						          trkPtFactPF_[trkEntry] = factorizedPtCorr(hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], trkRMinPF_[trkEntry], PFcent_p[ptPos], PFphiEta_p[ptPos], PFpt_p[ptPos], PFdelR_p[ptPos], 3);

							            trkPtCorrCalo_[trkEntry] = trkPt_[trkEntry]/factorizedPtCorr(hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], trkRMinCalo_[trkEntry], Calocent_p[ptPos], CalophiEta_p[ptPos], Calopt_p[ptPos], CalodelR_p[ptPos], 5);
								              trkPtFactCalo_[trkEntry] = factorizedPtCorr(hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], trkRMinCalo_[trkEntry], Calocent_p[ptPos], CalophiEta_p[ptPos], Calopt_p[ptPos], CalodelR_p[ptPos], 5);
									              }
										              break;
											            }
												        }

*/
