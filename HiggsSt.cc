#include "Pythia8/Pythia.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TMath.h"

using namespace Pythia8;

int main(){ 
Pythia pythia;
pythia.readString("Beams:idA =  11");  
pythia.readString("Beams:idB = -11");
pythia.readString("Beams:eCM = 1000.0");  
pythia.readString(" HiggsSM:ffbar2HZ = on"); 
//pythia.readString("Stat:showProcessLevel = off");
// Higgs-Strahlung process
  
pythia.init();
  
TFile * out_file = new TFile("output_HiggsSt_1000.root","RECREATE");


//pT distribution histogram
TH1F * hist1d_pT_mu_pos = new TH1F("hist1d_pT_mu+", "Z^{0} #rightarrow #mu^{+} + #mu^{-}", 100, 0., 500.);
TH1F * hist1d_pT_mu_neg = new TH1F("hist1d_pT_mu-", "Z^{0} #rightarrow #mu^{+} + #mu^{-}", 100, 0., 500.);
TH1F * hist1d_pT_e_pos = new TH1F("hist1d_pT_e+", "Z^{0} #rightarrow e^{+} + e^{-}", 100, 0., 500.);
TH1F * hist1d_pT_e_neg = new TH1F("hist1d_pT_e-", "Z^{0} #rightarrow e^{+} + e^{-}", 100, 0., 500.);
TH1F * hist1d_pT_b = new TH1F("hist1d_pT_b", "H #rightarrow b + #bar b", 100, 0., 500.);
TH1F * hist1d_pT_bbar = new TH1F("hist1d_pT_bbar", " H #rightarrow b + #bar b", 100, 0., 500.);
TH1F * hist1d_pT_Z = new TH1F("hist1d_pT_Z", "e^{+} + e^{-} #rightarrow Z^{0} + H", 100, 0., 500.);
TH1F * hist1d_pT_H = new TH1F("hist1d_pT_H", "e^{+} + e^{-} #rightarrow Z^{0} + H", 100, 0., 500.);



//pseudo-rapidity distribution histogram
TH1F * hist1d_eta_mu_pos = new TH1F("hist1d_eta_mu+", "Z^{0} #rightarrow #mu^{+} + #mu^{-}", 100, -5., 5.);
TH1F * hist1d_eta_mu_neg = new TH1F("hist1d_eta_mu-", "Z^{0} #rightarrow #mu^{+} + #mu^{-}", 100, -5., 5.);
TH1F * hist1d_eta_e_pos = new TH1F("hist1d_eta_e+", "Z^{0} #rightarrow e^{+} + e^{-}", 100, -5., 5.);
TH1F * hist1d_eta_e_neg = new TH1F("hist1d_eta_e-", "Z^{0} #rightarrow e^{+} + e^{-}", 100, -5., 5.);
TH1F * hist1d_eta_b = new TH1F("hist1d_eta_b", "H #rightarrow b + #bar b", 100, -5., 5.);
TH1F * hist1d_eta_bbar = new TH1F("hist1d_eta_bbar", " H #rightarrow b + #bar b", 100, -5., 5.);
TH1F * hist1d_eta_Z = new TH1F("hist1d_eta_Z", "e^{+} + e^{-} #rightarrow Z^{0} + H", 100, -5., 5.);
TH1F * hist1d_eta_H = new TH1F("hist1d_eta_H", "e^{+} + e^{-} #rightarrow Z^{0} + H", 100, -5., 5.);



//rapidity distribution histogram
TH1F * hist1d_y_mu_pos = new TH1F("hist1d_y_mu+", "Z^{0} #rightarrow #mu^{+} + #mu^{-}", 100, -5., 5.);
TH1F * hist1d_y_mu_neg = new TH1F("hist1d_y_mu-", "Z^{0} #rightarrow #mu^{+} + #mu^{-}", 100, -5., 5.);
TH1F * hist1d_y_e_pos = new TH1F("hist1d_y_e+", "Z^{0} #rightarrow e^{+} + e^{-}", 100, -5., 5.);
TH1F * hist1d_y_e_neg = new TH1F("hist1d_y_e-", "Z^{0} #rightarrow e^{+} + e^{-}", 100, -5., 5.);
TH1F * hist1d_y_b = new TH1F("hist1d_y_b", "H #rightarrow b + #bar b", 100, -5., 5.);
TH1F * hist1d_y_bbar = new TH1F("hist1d_y_bbar", " H #rightarrow b + #bar b", 100, -5., 5.);
TH1F * hist1d_y_Z = new TH1F("hist1d_y_Z", "e^{+} + e^{-} #rightarrow Z^{0} + H", 100, -5., 5.);
TH1F * hist1d_y_H = new TH1F("hist1d_y_H", "e^{+} + e^{-} #rightarrow Z^{0} + H", 100, -5., 5.);



//theta distribution histogram
TH1F * hist1d_theta_mu_pos = new TH1F("hist1d_theta_mu+", "Z^{0} #rightarrow #mu^{+} + #mu^{-}", 100, -5., 5.);
TH1F * hist1d_theta_mu_neg = new TH1F("hist1d_theta_mu-", "Z^{0} #rightarrow #mu^{+} + #mu^{-}", 100, -5., 5.);
TH1F * hist1d_theta_e_pos = new TH1F("hist1d_theta_e+", "Z^{0} #rightarrow e^{+} + e^{-}", 100, -5., 5.);
TH1F * hist1d_theta_e_neg = new TH1F("hist1d_theta_e-", "Z^{0} #rightarrow e^{+} + e^{-}", 100, -5., 5.);
TH1F * hist1d_theta_b = new TH1F("hist1d_theta_b", "H #rightarrow b + #bar b", 100, -5., 5.);
TH1F * hist1d_theta_bbar = new TH1F("hist1d_theta_bbar", " H #rightarrow b + #bar b", 100, -5., 5.);
TH1F * hist1d_theta_Z = new TH1F("hist1d_theta_Z", "e^{+} + e^{-} #rightarrow Z^{0} + H", 100, -5., 5.);
TH1F * hist1d_theta_H = new TH1F("hist1d_theta_H", "e^{+} + e^{-} #rightarrow Z^{0} + H", 100, -5., 5.);



//phi distribution histogram
TH1F * hist1d_phi_mu_pos = new TH1F("hist1d_phi_mu+", "Z^{0} #rightarrow #mu^{+} + #mu^{-}", 100, -5., 5.);
TH1F * hist1d_phi_mu_neg = new TH1F("hist1d_phi_mu-", "Z^{0} #rightarrow #mu^{+} + #mu^{-}", 100, -5., 5.);
TH1F * hist1d_phi_e_pos = new TH1F("hist1d_phi_e+", "Z^{0} #rightarrow e^{+} + e^{-}", 100, -5., 5.);
TH1F * hist1d_phi_e_neg = new TH1F("hist1d_phi_e-", "Z^{0} #rightarrow e^{+} + e^{-}", 100, -5., 5.);
TH1F * hist1d_phi_b = new TH1F("hist1d_phi_b", "H #rightarrow b + #bar b", 100, -5., 5.);
TH1F * hist1d_phi_bbar = new TH1F("hist1d_phi_bbar", " H #rightarrow b + #bar b", 100, -5., 5.);
TH1F * hist1d_phi_Z = new TH1F("hist1d_phi_Z", "e^{+} + e^{-} #rightarrow Z^{0} + H", 100, -5., 5.);
TH1F * hist1d_phi_H = new TH1F("hist1d_phi_H", "e^{+} + e^{-} #rightarrow Z^{0} + H", 100, -5., 5.);



//Angle between mu+ and mu- , e+ and e- in Z decay  
TH1F * hist1d_angle_mumu = new TH1F("hist1d_Angle_mumu", "Z^{0} #rightarrow #mu^{+} + #mu^{-}", 100, 0., TMath::Pi());
TH1F * hist1d_angle_ee = new TH1F("hist1d_Angle_ee", "Z^{0} #rightarrow e^{+} + e^{-}", 100, 0., TMath::Pi());
TH1F * hist1d_angle_bbbar = new TH1F("hist1d_Angle_bbbar", "H #rightarrow b + #bar b ", 100, 0., TMath::Pi());
  
  
//Mass invarient in  Z decay 
TH1F * hist1d_mass_mumu = new TH1F("hist1d_mass_Zmumu", "Z^{0} #rightarrow #mu^{+} + #mu^{-}", 100, 50., 150.);
TH1F * hist1d_mass_ee = new TH1F("hist1d_mass_Zee", "Z^{0} #rightarrow e^{+} + e^{-}", 100, 50., 150.);
TH1F * hist1d_mass_bbbar = new TH1F("hist1d_mass_bbbar", "H #rightarrow b + #bar b ", 100, 50., 150.);
  
  
//Angle between Z and H in e- + e+ -> H + Z process 
TH1F * hist1d_angle_ZH = new TH1F("hist1d_Angle_ZH", "e^{+} + e^{-} #rightarrow Z^{0} + H ", 300, 2.8 , 3.4);
TH1F * hist1d_pT_ZH = new TH1F("hist1d_pT_ZH", "e^{+} + e^{-} #rightarrow Z^{0} + H ", 300, -100, 100);


for (int iEvent = 0; iEvent < 200000; ++iEvent) // 200000 event generated
{
	
	if (!pythia.next()) continue;
	
    for (int i = 0; i < pythia.event.size(); ++i)
    
    {
    	
	if(pythia.event[i].id() == 23) 
    {
    	
		int dau1 = pythia.event[i].daughter1();
		int dau2 = pythia.event[i].daughter2();
    	
		if (pythia.event[dau1].id() == 13 and pythia.event[dau2].id() == -13 ) 
		// Sort out the process : Z0-> mu-(13) + mu+(-13)
    	{	
    		int j = i;
    		while(j>0)
    		{	
    			int mot1 = pythia.event[j].mother1();

    			if (pythia.event[mot1].id() != 23)
    			//actual parent 
    			
    			{	
					
    				//cout << pythia.event[j].name()<<"->" <<pythia.event[dau1].name() << " + " << pythia.event[dau2].name() << endl;		
    				TLorentzVector pmu1; // for mu-
    				double pt_mu1 = pythia.event[dau1].pT();
    				double eta_mu1 = pythia.event[dau1].eta();
    				double phi_mu1 = pythia.event[dau1].phi();
    				double m_mu1 = pythia.event[dau1].mCalc();
    				double y_mu1 = pythia.event[dau1].y();
    				double theta_mu1 = pythia.event[dau1].theta();  
    				pmu1.SetPtEtaPhiM(pt_mu1, eta_mu1, phi_mu1, m_mu1);
    				
    				TLorentzVector pmu2;// for mu+
    				double pt_mu2 = pythia.event[dau2].pT();
    				double eta_mu2 = pythia.event[dau2].eta();
    				double phi_mu2 = pythia.event[dau2].phi();
    				double m_mu2 = pythia.event[dau2].mCalc();
    				double y_mu2 = pythia.event[dau2].y(); 
    				double theta_mu2 = pythia.event[dau2].theta();
    				pmu2.SetPtEtaPhiM(pt_mu2, eta_mu2, phi_mu2, m_mu2);
    				
    				double angle = pmu1.Angle(pmu2.Vect()); //Angle between mu+ and mu- 
					double mass = (pmu1 + pmu2).M(); // Invarient mass calculation 
					
					
					hist1d_pT_mu_pos -> Fill(pt_mu2);
					hist1d_pT_mu_neg -> Fill(pt_mu1);
					hist1d_eta_mu_pos -> Fill(eta_mu2);
					hist1d_eta_mu_neg -> Fill(eta_mu1);
					hist1d_y_mu_pos -> Fill(y_mu2);
					hist1d_y_mu_neg -> Fill(y_mu1);
					hist1d_theta_mu_pos -> Fill(theta_mu2);
					hist1d_theta_mu_neg -> Fill(theta_mu1);
					hist1d_phi_mu_pos -> Fill(phi_mu2);
					hist1d_phi_mu_neg -> Fill(phi_mu1);
								
					hist1d_angle_mumu -> Fill(angle); 
					hist1d_mass_mumu -> Fill(mass);
    				break;
    			}
    			else
    			{	 				
 					j = mot1;	
    			}
    		}
    	
    	}
    	
    	else if(pythia.event[dau1].id() == 11 and pythia.event[dau2].id() == -11)
    	// Sort out the process : Z0 -> e-(11) + e+(-11)
    	{
    		
    		int j = i;
    		
    		while(j>0)
    		{	
    			int mot1 = pythia.event[j].mother1();
    			
    			if (pythia.event[mot1].id() != 23)
    			{	
    				//cout << pythia.event[j].name()<<"->" <<pythia.event[dau1].name() << " + " << pythia.event[dau2].name() << endl;
					
					TLorentzVector pe1; // for electron
    				double pt_e1 = pythia.event[dau1].pT();
    				double eta_e1 = pythia.event[dau1].eta();
    				double phi_e1 = pythia.event[dau1].phi();
    				double m_e1 = pythia.event[dau1].mCalc(); 
    				double y_e1 = pythia.event[dau1].y();
    				double theta_e1 = pythia.event[dau1].theta(); 
    				pe1.SetPtEtaPhiM(pt_e1, eta_e1, phi_e1, m_e1);
    				
    				TLorentzVector pe2; // for positron
    				double pt_e2 = pythia.event[dau2].pT();
    				double eta_e2 = pythia.event[dau2].eta();
    				double phi_e2 = pythia.event[dau2].phi();
    				double m_e2 = pythia.event[dau2].mCalc();
    				double y_e2 = pythia.event[dau2].y();
    				double theta_e2 = pythia.event[dau2].theta();  
    				pe2.SetPtEtaPhiM(pt_e2, eta_e2, phi_e2, m_e2);    			
    				
    				double angle = pe1.Angle(pe2.Vect()); //Angle between e+ and e- 
					double mass = (pe1 + pe2).M(); // invarient mass calculation
					
					
					hist1d_pT_e_pos -> Fill(pt_e2);
					hist1d_pT_e_neg -> Fill(pt_e1);
					hist1d_eta_e_pos -> Fill(eta_e2);
					hist1d_eta_e_neg -> Fill(eta_e1);
					hist1d_y_e_pos -> Fill(y_e2);
					hist1d_y_e_neg -> Fill(y_e1);
					hist1d_theta_e_pos -> Fill(theta_e2);
					hist1d_theta_e_neg -> Fill(theta_e1);
					hist1d_phi_e_pos -> Fill(phi_e2);
					hist1d_phi_e_neg -> Fill(phi_e1);
					
    				hist1d_angle_ee -> Fill(angle);
    				hist1d_mass_ee -> Fill(mass);
    				break;
    			}
    			
    			else
    			{	 				
 					j = mot1;	
    			}
    			
    		}
    		
    	}
    }
    
    
    else if(pythia.event[i].id() == 25) 

    {    
    	int b = pythia.event[i].daughter1();
    	int bbar = pythia.event[i].daughter2();
    	
    	if (pythia.event[b].id() == 5 and pythia.event[bbar].id() == -5) 
    	// Sort out the process  H -> b(5) + bbar(-5)
    	
    	{			
    		int j = i;
    	
    		while(j>0)
    		{	
    			int mot1 = pythia.event[j].mother1();
    				
    			if (pythia.event[mot1].id() != 25)
    			//actual parent 
    			
    			{	
    				
    				TLorentzVector pb1; // for b 
    				double pt_b = pythia.event[b].pT();
    				double eta_b = pythia.event[b].eta();
    				double phi_b = pythia.event[b].phi();
    				double m_b = pythia.event[b].mCalc();
    				double y_b = pythia.event[b].y();
    				double theta_b = pythia.event[b].theta();
    				pb1.SetPtEtaPhiM(pt_b, eta_b, phi_b, m_b);
    				
    				TLorentzVector pb2; // for bbar
    				double pt_bbar = pythia.event[bbar].pT();
    				double eta_bbar = pythia.event[bbar].eta();
    				double phi_bbar = pythia.event[bbar].phi();
    				double m_bbar = pythia.event[bbar].mCalc(); 
    				double y_bbar = pythia.event[bbar].y();
    				double theta_bbar = pythia.event[bbar].theta();
    				pb2.SetPtEtaPhiM(pt_bbar, eta_bbar, phi_bbar, m_bbar);    
    				

    				double angle = pb1.Angle(pb2.Vect()); //Angle between e+ and e- 
					double mass = (pb1 + pb2).M(); // invarient mass calculation    				
 			
					hist1d_pT_b -> Fill(pt_b);    				    			
					hist1d_pT_bbar -> Fill(pt_bbar);
					hist1d_eta_b -> Fill(eta_b);
					hist1d_eta_bbar -> Fill(eta_bbar);
					hist1d_y_b -> Fill(y_b);
					hist1d_y_bbar -> Fill(y_bbar);
					hist1d_theta_b -> Fill(theta_b);
					hist1d_theta_bbar -> Fill(theta_bbar);
					hist1d_phi_b -> Fill(phi_b);
					hist1d_phi_bbar -> Fill(phi_bbar);
					
					hist1d_angle_bbbar -> Fill(angle);
					hist1d_mass_bbbar -> Fill(mass);
			
					break;
    			}
    			
    			else
    			{	 				
 						j = mot1;	
    			}
    		}
    	
    	}
    	
    	else continue;
    }
    
     
    else if(pythia.event[i].id() == 11) 
    {
    	
    	int higgs = pythia.event[i].daughter1();
    	int zbose = pythia.event[i].daughter2();
    	
    	if (pythia.event[higgs].id() == 25 and pythia.event[zbose].id() == 23)
    	//Sort out the process:  e- + e+ -> Z + H 
		
		
		{
					//cout <<"->" <<pythia.event[dau1].name() << " + " << pythia.event[dau2].name() << endl;
    				TLorentzVector pH;
    				double pt_pH = pythia.event[higgs].pT();
    				double eta_pH = pythia.event[higgs].eta();
    				double phi_pH = pythia.event[higgs].phi();
    				double m_pH = pythia.event[higgs].mCalc();
					double y_pH = pythia.event[higgs].y(); 
					double theta_pH = pythia.event[higgs].theta();
    				pH.SetPtEtaPhiM(pt_pH, eta_pH, phi_pH, m_pH);
    				
    				TLorentzVector pZ;
    				double pt_pZ = pythia.event[zbose].pT();
    				double eta_pZ = pythia.event[zbose].eta();
    				double phi_pZ = pythia.event[zbose].phi();
    				double m_pZ = pythia.event[zbose].mCalc();
    				double y_pZ = pythia.event[zbose].y();
    				double theta_pZ = pythia.event[zbose].theta();
    				pZ.SetPtEtaPhiM(pt_pZ, eta_pZ, phi_pZ, m_pZ);    			
    				
    	
    				double angle = pH.Angle(pZ.Vect()); // Angle between H and Zboson
    				
    				hist1d_pT_Z -> Fill(pt_pZ);
					hist1d_pT_H -> Fill(pt_pH);
					hist1d_eta_Z -> Fill(eta_pZ);
					hist1d_eta_H -> Fill(eta_pH);
    				hist1d_y_Z -> Fill(y_pZ);
					hist1d_y_H -> Fill(y_pH);
					hist1d_theta_Z -> Fill(theta_pZ);
					hist1d_theta_H -> Fill(theta_pH);
					hist1d_phi_Z -> Fill(phi_pZ);
					hist1d_phi_H -> Fill(phi_pH);
    				
    				hist1d_angle_ZH -> Fill(angle); // Angle between Z and Higgs
    				hist1d_pT_ZH -> Fill(pt_pH - pt_pZ); // change in tranverse momenta
		}
	
		else continue;	
    }
    
	else continue;
	
	}

}
  
//pythia.stat();
out_file->Write();
out_file->Close();
return 0;
}


