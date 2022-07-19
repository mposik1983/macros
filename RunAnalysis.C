#define Analyze_cxx
#include <TVector3.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <fstream>


void RunAnalysis( double pTmin, double pTmax, double etamin, double etamax)
{
//Files
TString Tag        ="MiniTPC_800";
TString inRootFile =Form("./%s/RootFiles/%s_G4EICDetector",Tag.Data(),Tag.Data());
TString outHists =  Form("./%s/RootFiles/Hists/Histogram_%s",Tag.Data(),Tag.Data());
TString outHists2 =  Form("./%s/RootFiles/Hists/Histogram_nHitLayer_%s",Tag.Data(),Tag.Data());
TString outDat_p =  Form("./%s/dat/%s_p",Tag.Data(),Tag.Data());
TString outDat_pT = Form("./%s/dat/%s_pT",Tag.Data(),Tag.Data());
TString outDat_dca2d = Form("./%s/dat/%s_dca2d",Tag.Data(),Tag.Data());
//extract information from the filename
//Alternatively, take the information as input, and construct the filename from it.
TString infile = Form("%s_%.2f_%.2f_%.2f_%.2f_g4tracking_eval.root", inRootFile.Data(),pTmin, pTmax, etamin, etamax);

//check if file exists
if(gSystem->AccessPathName(infile))
{
        std::cout << "file " << infile << " does not exist" << std::endl;
	return;
} else 
{
        std::cout << "Found File: " << infile << std::endl;
}

//get the file
TFile myFile(infile.Data());

   
   Int_t           event;
   Int_t           gtrackID;
   Int_t           gflavor;
   Float_t         gpx;
   Float_t         gpy;
   Float_t         gpz;
   Float_t         gvx;
   Float_t         gvy;
   Float_t         gvz;
   Float_t         gvt;
   Int_t           trackID;
   Int_t           charge;
   Int_t           nhits;
   Float_t         px;
   Float_t         py;
   Float_t         pz;
   Float_t         pcax;
   Float_t         pcay;
   Float_t         pcaz;
   Float_t         dca2d;

   Float_t         hpDIRC_px;
   Float_t         hpDIRC_py;
   Float_t         hpDIRC_x;
   Float_t         hpDIRC_y;
   Float_t         hpDIRC_z;
   Float_t         hpDIRC_pz;
   Float_t         hpDIRC_proj_px;
   Float_t         hpDIRC_proj_py;
   Float_t         hpDIRC_proj_pz;

   Int_t           nHit_G4HIT_SVTX_0;
   Int_t           nHit_G4HIT_SVTX_1;
   Int_t           nHit_G4HIT_SVTX_2;
   Int_t           nHit_G4HIT_BARR_0;
   Int_t           nHit_G4HIT_BARR_1;
   Int_t           nHit_G4HIT_RWELL_0;
   Int_t           nHit_G4HIT_RWELL_1;
   Int_t           nHit_G4HIT_CTTL_0;
   Int_t           nHit_G4HIT_miniTPC_GAS;
//get the tree
TTree *tree;
myFile.GetObject("tracks", tree);
   
//get number of entries
   Int_t nentries = (Int_t)tree->GetEntries();
   
//set the branch addresses
   tree->SetBranchAddress("event", &event);
   tree->SetBranchAddress("gtrackID", &gtrackID);
   tree->SetBranchAddress("gflavor", &gflavor);
   tree->SetBranchAddress("gpx", &gpx);
   tree->SetBranchAddress("gpy", &gpy);
   tree->SetBranchAddress("gpz", &gpz);
   tree->SetBranchAddress("gvx", &gvx);
   tree->SetBranchAddress("gvy", &gvy);
   tree->SetBranchAddress("gvz", &gvz);
   tree->SetBranchAddress("gvt", &gvt);
   tree->SetBranchAddress("trackID", &trackID);
   tree->SetBranchAddress("charge", &charge);
   tree->SetBranchAddress("nhits", &nhits);
   tree->SetBranchAddress("px", &px);
   tree->SetBranchAddress("py", &py);
   tree->SetBranchAddress("pz", &pz);
   tree->SetBranchAddress("pcax", &pcax);
   tree->SetBranchAddress("pcay", &pcay);
   tree->SetBranchAddress("pcaz", &pcaz);
   tree->SetBranchAddress("dca2d", &dca2d);

   tree->SetBranchAddress("hpDIRC_px", &hpDIRC_px);
   tree->SetBranchAddress("hpDIRC_py", &hpDIRC_py);
   tree->SetBranchAddress("hpDIRC_pz", &hpDIRC_pz);
   tree->SetBranchAddress("hpDIRC_x", &hpDIRC_x);
   tree->SetBranchAddress("hpDIRC_y", &hpDIRC_y);
   tree->SetBranchAddress("hpDIRC_z", &hpDIRC_z);
   tree->SetBranchAddress("hpDIRC_proj_px", &hpDIRC_proj_px);
   tree->SetBranchAddress("hpDIRC_proj_py", &hpDIRC_proj_py);
   tree->SetBranchAddress("hpDIRC_proj_pz", &hpDIRC_proj_pz);

   tree->SetBranchAddress("nHit_G4HIT_SVTX_0", &nHit_G4HIT_SVTX_0);
   tree->SetBranchAddress("nHit_G4HIT_SVTX_1", &nHit_G4HIT_SVTX_1);
   tree->SetBranchAddress("nHit_G4HIT_SVTX_2", &nHit_G4HIT_SVTX_2);
   tree->SetBranchAddress("nHit_G4HIT_BARR_0", &nHit_G4HIT_BARR_0);
   tree->SetBranchAddress("nHit_G4HIT_BARR_1", &nHit_G4HIT_BARR_1);
   tree->SetBranchAddress("nHit_G4HIT_RWELL_0", &nHit_G4HIT_RWELL_0);
   tree->SetBranchAddress("nHit_G4HIT_RWELL_1", &nHit_G4HIT_RWELL_1);
   tree->SetBranchAddress("nHit_G4HIT_CTTL_0", &nHit_G4HIT_CTTL_0);
   tree->SetBranchAddress("nHit_G4HIT_miniTPC_GAS", &nHit_G4HIT_miniTPC_GAS);

//Run the Analysis
	

//Defining Histograms and variables
	TH1D *h_dca2d = new TH1D("h_dca2d", "dca2d", 1000, -0.03, 0.03);
	TH1D *h_momRes = new TH1D("h_momRes", "(reco_p - truth_p)/truth_p", 1000, -1, 1);
	TH1D *h_ptRes = new TH1D("h_ptRes", "(reco_pt - truth_pt)/truth_pt", 1000, -1, 1);
	TH1D *h_pTruth = new TH1D("h_pTruth", "Truth Momentum", 1000, 0, 30);
	TH1D *h_etaTruth = new TH1D("h_etaTruth", "Truth Eta", 1000, -4.5, 4.5);
	TH1D *h_ptTruth = new TH1D("h_ptTruth", "Truth pT", 1000, 0, 30);
	TH1D *h_thetaRes = new TH1D("h_thetaRes", "(reco_theta - truth_theta)", 1000, -1, 1);
	TH1D *h_phiRes = new TH1D("h_phiRes", "(reco_phi - truth_phi)", 1000, -1, 1);
	TH2D *h_phi_theta_true = new TH2D("h_phi_theta_true", "theta_true,phi_true", 1000, -TMath::Pi(), TMath::Pi(),100,-2.0*TMath::Pi(),2.0*TMath::Pi());
	TH2D *h_phi_theta_reco = new TH2D("h_phi_theta_reco", "theta_reco,phi_reco", 1000, -TMath::Pi(), TMath::Pi(),100,-2.0*TMath::Pi(),2.0*TMath::Pi());
	
	TH1D *h_DIRC_thetaRes = new TH1D("h_DIRC_thetaRes", "DIRC Theta Resolution", 1000, -0.02, 0.02);
	TH1D *h_DIRC_phiRes   = new TH1D("h_DIRC_phiRes", "DIRC Phi Resolution", 1000, -0.02, 0.02);

        TH1D *h_hitlayer = new TH1D(Form("h_hitlayer_%.2f_%.2f",pTmin,pTmax),"Number of Layers Hit",100,0,25);
        //TH2D *h_hits = new TH2D("h_hits","Number of Hit Layers",11,0,10,4,-1,1);//hits,p,eta
         
        Int_t nHitLayers = 0;        	
for (int i = 0; i < nentries; i++)
{
	tree->GetEntry(i);
        nHitLayers = 0;
        //cout << "Entry: " << i << endl;
        //cout << "event: " << event << " nhit: " << nHit_G4HIT_BARR_0 << endl;
        /*
        if(nHit_G4HIT_SVTX_0 > 0) ++nHitLayers;
        if(nHit_G4HIT_SVTX_1 > 0) ++nHitLayers;
        if(nHit_G4HIT_SVTX_2 > 0) ++nHitLayers;
        if(nHit_G4HIT_BARR_0 > 0) ++nHitLayers;
        if(nHit_G4HIT_BARR_1 > 0) ++nHitLayers;
        if(nHit_G4HIT_RWELL_0 > 0) ++nHitLayers;
        if(nHit_G4HIT_RWELL_1 > 0) ++nHitLayers;
        if(nHit_G4HIT_CTTL_0 > 0) ++nHitLayers;
        h_hitlayer->Fill(nHitLayers);
        cout << "nHitLayers: " << nHitLayers << endl;
        */
	//fill in analysis here
	TVector3 truthP(gpx, gpy, gpz);
	TVector3 recoP(px, py, pz);
        Double_t eta_true = truthP.Eta();
        Double_t theta_true = truthP.Theta();
        Double_t phi_true = truthP.Phi();
        Double_t theta_reco = recoP.Theta();
        Double_t phi_reco = recoP.Phi();
 
        TVector3 DIRC_truthP(hpDIRC_px,hpDIRC_py,hpDIRC_pz);
        TVector3 DIRC_recoP(hpDIRC_proj_px,hpDIRC_proj_py,hpDIRC_proj_pz);
        TVector3 DIRC_z(0,0,hpDIRC_z);
        //Double_t DIRC_dth = TMath::ATan((DIRC_recoP - DIRC_truthP).Dot( (DIRC_truthP).Cross(DIRC_z) )/(DIRC_truthP.Mag()*( (DIRC_truthP.Cross(DIRC_z)).Mag()) ));
        //Double_t DIRC_dphi = TMath::ATan((DIRC_recoP - DIRC_truthP).Dot( ( (DIRC_truthP).Cross(DIRC_z) ).Cross(DIRC_truthP) )/(DIRC_truthP.Mag()*(DIRC_truthP.Cross(DIRC_z).Cross(DIRC_truthP)).Mag()));
        //Double_t DIRC_dth = DIRC_recoP.Theta() - DIRC_truthP.Theta();
        //Double_t DIRC_dphi = DIRC_recoP.Phi() - DIRC_truthP.Phi();
  
	if (trackID != -9999)
	{
                h_dca2d->Fill(dca2d);
		h_momRes->Fill( (recoP.Mag() - truthP.Mag())/truthP.Mag());  	
		h_ptRes->Fill( (recoP.Pt() - truthP.Pt())/truthP.Pt());  	
		h_pTruth->Fill( (truthP.Mag()));  	
                h_etaTruth->Fill(eta_true);
		h_ptTruth->Fill( (truthP.Pt()));  	
		h_thetaRes->Fill( (theta_reco - theta_true));  	
		h_phiRes->Fill( (phi_reco - phi_true));  	
		h_phi_theta_true->Fill( theta_true,phi_true);  	
		h_phi_theta_reco->Fill(theta_reco,phi_reco);  			

               // h_DIRC_phiRes->Fill(DIRC_dphi);
               // h_DIRC_thetaRes->Fill(DIRC_dth);

                if(nHit_G4HIT_SVTX_0 > 0) ++nHitLayers;
                if(nHit_G4HIT_SVTX_1 > 0) ++nHitLayers;
                if(nHit_G4HIT_SVTX_2 > 0) ++nHitLayers;
                if(nHit_G4HIT_BARR_0 > 0) ++nHitLayers;
                if(nHit_G4HIT_BARR_1 > 0) ++nHitLayers;
                if(nHit_G4HIT_RWELL_0 > 0) ++nHitLayers;
                if(nHit_G4HIT_RWELL_1 > 0) ++nHitLayers;
                if(nHit_G4HIT_CTTL_0 > 0) ++nHitLayers;
                if(nHit_G4HIT_miniTPC_GAS > 0) ++nHitLayers;
                h_hitlayer->Fill(nHitLayers);
                //cout << "TPC Hits: " << nHit_G4HIT_miniTPC_GAS << endl;
	}	

}
        //h_hits->Fill(h_etaTruth->GetMean(),h_pTruth->GetMean(),h_hitlayer->GetMean());
	TFile outfile2(Form("%s_%.2f_%.2f.root",outHists2.Data(),etamin,etamax), "UPDATE");
        h_hitlayer->Write();
        //h_hits->Write();
  
        outfile2.Close();   

	TFile outfile(Form("%s_%.2f_%.2f.root",outHists.Data(),etamin,etamax), "UPDATE");
	fstream outfile_p_text;
         outfile_p_text.open(Form("%s_%.2f_%.2f.dat",outDat_p.Data(),etamin,etamax), ofstream::out | ofstream::app);
	fstream outfile_pT_text;
         outfile_pT_text.open(Form("%s_%.2f_%.2f.dat",outDat_pT.Data(),etamin,etamax), ofstream::out | ofstream::app);
	fstream outfile_dca2d_text;
         outfile_dca2d_text.open(Form("%s_%.2f_%.2f.dat",outDat_dca2d.Data(),etamin,etamax), ofstream::out | ofstream::app);

	//Here I should insert lines Getting the 2D histograms I want to update
	
        //rough fit to get mean

	//add line to fit 1D histogram of delta p to get sigmap
	TF1 *gausfitmom = new TF1("gausfitmom", "gaus",-0.02,0.02);
	//TF1 *gausfitpt = new TF1("gausfitpt", "gaus",roughfitpt->GetParameter(1)-2.0*roughfitpt->GetParameter(2), roughfitpt->GetParameter(1)+2.0*roughfitpt->GetParameter(2));
	TF1 *gausfitpt = new TF1("gausfitpt", "gaus",-0.02,0.02);
	h_momRes->Fit(gausfitmom, "RQ");
	h_ptRes->Fit(gausfitpt, "RQ");

	TF1 *gausfit_dca2d = new TF1("gausfit_dca2d", "gaus",-0.01,0.01);
	h_dca2d->Fit(gausfit_dca2d, "RQ");

        double sigma_dca2d = gausfit_dca2d->GetParameter(2);
        double sigma_dca2d_error = gausfit_dca2d->GetParError(2);

	double sigmaP = gausfitmom->GetParameter(2);
	double sigmaPt = gausfitpt->GetParameter(2);
	double sigmaP_error = gausfitmom->GetParError(2);
	double sigmaPt_error = gausfitpt->GetParError(2);

        double pAvg = h_pTruth->GetMean();			
        double pRMS = h_pTruth->GetRMS();

        double pT_Avg = h_ptTruth->GetMean();
        double pT_RMS = h_ptTruth->GetRMS();

        cout << Form("sigmaP: %.5f sigmaP error: %.5f\n",sigmaP,sigmaP_error);
        outfile_p_text << Form("%.5f  %.5f %.5f %.5f \n",pAvg, sigmaP, pRMS, sigmaP_error);

        cout << Form("sigmaPt: %.5f sigmaPt error: %.5f\n",sigmaPt,sigmaPt_error);
        outfile_pT_text << Form("%.5f  %.5f %.5f %.5f \n",pT_Avg, sigmaPt, 0.0, sigmaP_error);

        cout << Form("sigma_dca2d: %.5f sigma_dca2d error: %.5f\n",sigma_dca2d,sigma_dca2d_error);
        outfile_dca2d_text << Form("%.5f  %.5f %.5f %.5f \n",pT_Avg, 10000.0*sigma_dca2d, 0.0, sigma_dca2d_error); //um

	h_momRes->SetNameTitle(Form("h_momRes_P_%0.1f_pT_%0.1lf_Eta_%0.1lf_%0.1lf", pAvg, pTmin,  etamin, etamax), Form("(reco_p - truth_p)/truth_p for tracks with pt = %0.1f GeV/c and %0.1lf < eta < %0.1lf -Smear", pTmin, etamin, etamax) );
	
	h_ptRes->SetNameTitle(Form("h_ptRes_P_%0.1f_pT_%0.1lf_Eta_%0.1lf_%0.1lf",pT_Avg, pTmin, etamin, etamax), Form("(reco_pt - truth_pt)/truth_pt for tracks with pt = %0.1f GeV/c and %0.1lf < eta < %0.1lf - Smear", pTmin,etamin, etamax) );
	h_dca2d->SetNameTitle(Form("h_dca2d_pT_%0.1f_pT_%0.1f_Eta_%0.1f_%0.1f",pTmin, pTmax, etamin, etamax), Form("dca2d for tracks with %.1f < pt < %0.1f GeV/c and %0.1f < eta < %0.1f", pTmin,pTmax,etamin, etamax) );

         

	h_dca2d->Write();
	h_momRes->Write();
	h_ptRes->Write();
        h_pTruth->Write();
        h_thetaRes->Write();
        h_phiRes->Write();
        h_phi_theta_true->Write();
        h_phi_theta_reco->Write();
        
        //h_DIRC_phiRes->Write();
        //h_DIRC_thetaRes->Write();

	outfile.Close();
	//infile.Close();
	outfile_p_text.close();
	outfile_pT_text.close();
	outfile_dca2d_text.close();

//Close the input file
myFile.Close();

}
