#include <iostream>
using namespace::std;
#include "TCanvas.h"
#include "TF1.h"
#include "TH1.h"
#include "TH1F.h"
#include "TMath.h"
#include "math.h"
#include "TRandom.h"

// #include "/home/bottle/Documents/PhD/energyCalibration.cpp"
#include "/home/l_mayilyan/Documents/PhD/workspace/sls_detectors_package/slsDetectorCalibration/energyCalibration.cpp"

#define SCALING_SIZE		2		// Bin size of the signal local histograms

double Eff = 0;
double en = 0;
int pulse_count = 0;

// TH1F *Q_hist;
// TH1F *pos_hist;
// TH1F *Q_trap;
// TH1F *Q_trap1;
// TH2F *Q_trap2;
TH1F *Q_wnoise;
TH1F *en_wnoise;


// TH1F *en_hist = new TH1F("en_hist", "Energy data", 2000, 0, 2000);

///////////////////////////////////////////////////
// Version with three stripes

/*
Double_t get_Q(Double_t pos)
{
	// double Q_0 = 3000;	// Number of electrons generated


	if (((pos >=d) && (pos <= p)) || ((pos >=d+p) && (pos <= 2*p)) || ((pos >=d+2*p) && (pos <= 3*p)))
	{
		// Flat part
		Eff = Q_0;
		// cout << pos << " Middle region; Q = " << Eff << endl;
	}
	else if (pos < d)
	{
		// Left slope /
		Eff = Q_0 * pos/d;
		// cout << pos << " Left slope region; Q = " << Eff << endl;
	}
	else if ((pos > 3*p) && (pos < (3*p + d)))
	{
		// Right slope \
		Eff = Q_0 * (3*p+d-pos)/d;
		// cout << pos << " Right slope region; Q = " << Eff << endl;
		// cout << (3*p+d-pos)/d << " !!!!!!!!\n" << endl;
	}
	else if ((pos > p) && (pos < p+d))
	{
		Eff = Q_0 * (p+d-p101010os)/d + Q_0 * (pos-p)/d;
		// cout << pos << " First sharing region; Q = " << Eff << endl;
	}
	else if ((pos > 2*p) && (pos < 2*p + d))
	{
		Eff =  Q_0 * (2*p+d-pos)/d + Q_0 * (pos-2*p)/d;
		// cout << pos << " Second sharing region; Q = " << Eff << endl;
	}
	else
		Eff = 0;


	return Eff;
}
*/
///////////////////////////////////////////////////


double get_Q(Double_t pos, Double_t Q_0, Double_t p, Double_t d)
{

	if ((pos > -(p-d)/2) && (pos < (p-d)/2))
	{
		// flat part
		Eff = Q_0;
	}
	else if ((pos > -(p+d)/2) && (pos < -(p-d)/2))
	{
		// left slope
		Eff = Q_0 * (2*pos+p+d)/(2*d);
	}
	else if ((pos > (p-d)/2) && (pos < (p+d)/2))
	{
		// right slope
		Eff = Q_0 * (p+d-2*pos)/(2*d);
	}
	else
		Eff = 0;


	// Q_hist->Fill(Eff);
	// pos_hist->Fill(pos);

	return Eff;
	// return Q_0;
}
	


Double_t exp_func(Double_t *x, Double_t *par)
{
	Float_t xx =x[0];

	Double_t f = TMath::E()/par[0] * xx * exp(-xx/par[0]);

	return f;
}

TH1F *get_exp(int sig_dec_time)
{
	TH1F* loc_hist = new TH1F("loc_hist", "Single pulse", 1500/SCALING_SIZE, 0, 1500);
	double loc_hist_val;

	TF1 *expo = new TF1("Exp curve", exp_func, 0, 1500, 1);
	expo->SetParameter(0,sig_dec_time);
	// Creating local histograms from the function
	for (int ibin = 1; ibin < loc_hist->GetNbinsX()+1; ibin++)
	{
		loc_hist_val = expo->Eval(loc_hist->GetBinCenter(ibin));
		loc_hist->SetBinContent(ibin,loc_hist_val);
	}

	expo->SetParName(0,"sig_dec_time");
	// en_hist->Fill(loc_hist->GetMaximum());

	delete expo;
	return loc_hist;
}

TH1F *get_s_curve(TH1F *l_hist)
{
	// TCanvas *c2 = new TCanvas("S Curve");

	// TH1F* s_curve = new TH1F();

	double hist_val_old = 0;
	double hist_val_cur = 0;
	double hist_val_new = 0;

	int hist_max_val = l_hist->GetBinContent(l_hist->GetMaximumBin());

	TH1F* s_curve = new TH1F("s_curve", "S Curve", hist_max_val, 0, hist_max_val);

	int X_COUNT = l_hist->GetNbinsX();

	double *min = (double*)malloc(X_COUNT*sizeof(int));
	if(!min)
	{
		cout << "Error! memory not allocated.";
		exit(0);
	}

	double *max = (double*)malloc(X_COUNT*sizeof(int));
	if(!max)
	{
		cout << "Error! memory not allocated.";
		exit(0);
	}

	int count = 0;
	int j = 0;
	for (int i = 1; i <= X_COUNT; i++)
	{
		hist_val_old = l_hist->GetBinContent(i-1);

		j = i;
		while((l_hist->GetBinContent(j) == hist_val_old) && (j < X_COUNT))
			j++;
		if (X_COUNT == j)
			break;
		hist_val_cur = l_hist->GetBinContent(j);

		i = j;

		// find next level
		j++;
		while((l_hist->GetBinContent(j) == hist_val_cur) && (j < X_COUNT))
			j++;
		if (X_COUNT == j)
			break;
		hist_val_new = l_hist->GetBinContent(j);

		// suppose first comes a peak
		if ((hist_val_old < hist_val_cur) && (hist_val_new < hist_val_cur))		// peak
			max[count] = hist_val_cur;
		else
		if ((hist_val_cur < hist_val_old) && (hist_val_cur < hist_val_new))		// min
		{
			min[count] = hist_val_cur;
			count++;
		}
	}


	cout << "Now will count photons at given threshold\n";

	// Number of counts for certain threshold
	int step = 1;
	// for (int i = 0; i <= hist_max_val; i+=step)
	// {
		int min_count = 0;
		int max_count = 0;
		for (int j = 0; j < count; j++)
		{
			if (min[j] >= 7000)
				min_count++;
			if (max[j] >= 7000)
				max_count++;
		}

		cout << "Number of counts at given threshold " << (max_count - min_count)/step << endl;

		// s_curve->SetBinContent(s_curve->FindBin(i/step), (max_count - min_count)/step);
	// }


	int step = 1;
	for (int i = 0; i <= hist_max_val; i+=step)
	{
		int min_count = 0;
		int max_count = 0;
		for (int j = 0; j < count; j++)
		{
			if (min[j] >= i)
				min_count++;
			if (max[j] >= i)
				max_count++;
		}
		s_curve->SetBinContent(s_curve->FindBin(i/step), (max_count - min_count)/step);
	}

	if (min)
		free(min);
	if (max)
		free(max);

	return s_curve;
}

double get_half_en_count(TH1F *l_hist, int E_0)
{
	// TCanvas *c2 = new TCanvas("S Curve");

	// TH1F* s_curve = new TH1F();

	double hist_val_old = 0;
	double hist_val_cur = 0;
	double hist_val_new = 0;

	int hist_max_val = l_hist->GetBinContent(l_hist->GetMaximumBin());

	int X_COUNT = l_hist->GetNbinsX();

	double *min = (double*)malloc(X_COUNT*sizeof(int));
	if(!min)
	{
		cout << "Error! memory not allocated.";
		exit(0);
	}

	double *max = (double*)malloc(X_COUNT*sizeof(int));
	if(!max)
	{
		cout << "Error! memory not allocated.";
		exit(0);
	}

	int count = 0;
	int j = 0;
	for (int i = 1; i <= X_COUNT; i++)
	{
		hist_val_old = l_hist->GetBinContent(i-1);

		j = i;
		while((l_hist->GetBinContent(j) == hist_val_old) && (j < X_COUNT))
			j++;
		if (X_COUNT == j)
			break;
		hist_val_cur = l_hist->GetBinContent(j);

		i = j;

		// find next level
		j++;
		while((l_hist->GetBinContent(j) == hist_val_cur) && (j < X_COUNT))
			j++;
		if (X_COUNT == j)
			break;
		hist_val_new = l_hist->GetBinContent(j);

		// suppose first comes a peak
		if ((hist_val_old < hist_val_cur) && (hist_val_new < hist_val_cur))		// peak
			max[count] = hist_val_cur;
		else
		if ((hist_val_cur < hist_val_old) && (hist_val_cur < hist_val_new))		// min
		{
			min[count] = hist_val_cur;
			count++;
		}
	}


	cout << "Now will count photons at given threshold\n";


	cout << E_0*1000/2 << " E_0*1000/2\n";
	double th_lev = E_0*1000/2;
	// Number of counts for certain threshold
	int step = 1;
	int min_count = 0;
	int max_count = 0;
	for (int j = 0; j < count; j++)
	{
		if (min[j] >= th_lev)
			min_count++;
		if (max[j] >= th_lev)
			max_count++;
	}

	cout << "Number of counts at given threshold " << (max_count - min_count)/step << endl;

	if (min)
		free(min);
	if (max)
		free(max);

	return (max_count - min_count)/step;
}



TH1F *waveform(int E_0, double gauss_rms, int sig_dec_time, int count_rate, int num_pulses, double d)
{
	// double gen_electrons = E_0/3.6;	// Number of electrons generated

	// Average time between photons in ns (1E3 for having ns)
	Double_t exp_time_tau = 1E6/count_rate;
	Double_t gate=num_pulses*exp_time_tau;
	TH1F *hist = new TH1F("hist", "Signal data", gate/SCALING_SIZE, 0, gate);
	// TH1F *h0 = new TH1F("h0", "Photons", gate/SCALING_SIZE, 0, gate);
	TH1F *loc_hist;

	// new version
	// for (int i = en_wnoise->GetXaxis()->GetXmin(); i < en_wnoise->GetXaxis()->GetXmax(); ++i)
	// {
	// 	hist->SetBinContent(i,0);
	// }

	loc_hist = get_exp(sig_dec_time);
	Double_t gaus_noise;
	
	double p = 50;
	int SP=0;
	int T;
	double charge_pos;
	double Q;
	int start_pos;
	int i = 0;
	for (i = 0; i < 1.5*num_pulses; i++)
	{
		TRandom2 r2(0);

		T = r2.Exp(exp_time_tau);
		SP += T;

		if (SP >= gate)
			break;

		// Version for 3 strips
		// charge_pos = gRandom->Rndm() * (3*p+d);

		// Version with one strip
		charge_pos = r2.Uniform(-(p+d)/2, (p+d)/2);

		// Gauss rms is converted from eletrons to enenrgy units 
		gaus_noise = r2.Gaus(0,gauss_rms);
		// cout << "gauss " << gaus_noise << endl;
		// Getting the charge colection efficiency
		Q = get_Q(charge_pos, E_0, p, d) + gaus_noise;
		// cout << "Q = " << Q << endl;
		// Q_trap1->SetBinContent(charge_pos+(p+d)/2, Q/140);
		// Q_trap2->Fill(charge_pos,Q,1);

		// new version
		// Q_wnoise->Fill(Q);

		// Getting rid of negative amplitudes coming from gaus_noise negative values
		if (Q < 0)
			Q = 0;

		// h0->Fill(SP,Q);
		// Filling the main histogram from small ones
		for (int h_start = 0; h_start < loc_hist->GetNbinsX(); h_start++)
		{
			start_pos = SP+h_start*SCALING_SIZE;
			if (start_pos >= gate)
				break;
			else
			{
				// Filling the bins	
				hist->Fill(start_pos, Q*loc_hist->GetBinContent(h_start+1));
			}
		}
	}
	delete loc_hist;

	pulse_count = i;
	cout << "Generated number of photons i: " << pulse_count << endl;

	// Drawing the energy distribution
	// en_hist->Draw();

	// cout << gate << endl;
	// hist->Draw();

	// TFile hfile("rate_scan.root","RECREATE","Rate scan for 14 KeV");
	// int rn = 5;
	// TTree *Tr = new TTree("Tr","Tree of waveforms for 14KeV");

	// // TH1F *hp;
	// // Tr->Branch("hist", "TH1F", &hp, 32000, 0);
	// // hp=hist;

	// Tr->Branch("rn", &rn, 16000);
	// Tr->Branch("hist", "TH1F", &hist, 32000, 0);

	// // rn =


	// Tr->Fill();

	// Tr->Print();
	// hfile.Write();
	// // create_wf(hist, count_rate);

	// return h0;
	return hist;
}


void read_tree()
{
	TFile *f = new TFile("rate_scan.root");
	TTree *Tr = (TTree*)f->Get("Tr");
	TH1F *hist = 0;
	int rn;

	Tr->SetBranchAddress("rn",&rn);
	Tr->SetBranchAddress("hist",&hist);
	Tr->GetEntry(0);
	cout << rn << endl;
	hist->Draw();
}

// create_wf(TH1F *hist, int count_rate)
// {

// }


TH1F *main_func(int E_0 = 14, double gauss_rms = 1244,
				int sig_dec_time = 85, int count_rate = 10,
				int num_pulses = 10000, double d = 23.6151)
{
	cout << "\n////////////////////////////////////////////////////////////////////////////////\n\n";
	cout << "Energy of the photons in keV:\t" << E_0 << endl;
	cout << "RMS of gaussian noise:\t\t" << gauss_rms << endl;
	cout << "Analog signal decay time in ns:\t" << sig_dec_time << endl;
	cout << "Counting rate in kHz:\t\t" << count_rate << endl;
	cout << "Number of pulses:\t\t" << num_pulses << endl;
	cout << "Charge sharing distance:\t" << d << endl;
	cout << "\n////////////////////////////////////////////////////////////////////////////////\n\n";

	// Q_hist = new TH1F("Q_hist", "Q_hist", 25000, 0, 25000);
	// pos_hist = new TH1F("pos_hist", "pos_hist", 900, -45, 45);
	// Q_trap = new TH1F("Q_trap", "Q_trap", 674, -33.7, 33.7);
	// Q_trap1 = new TH1F("Q_trap1", "Q_trap1", 70, -35, 35);
	// Q_trap2 = new TH2F("Q_trap2", "Q_trap2", 70, -35, 35, 200, 0, 20000);

	// new version
	// Q_wnoise = new TH1F("Q_wnoise", "Q_wnoise", 30000, 0, 30000);

	// new version
	// en_wnoise = new TH1F("en_wnoise", "en_wnoise", 30000, 0, 30000);

	// comment Q_trap filling when getting pos_hist

	// for (int i = -350; i <= 350; i++)
	// {
	// 	// cout << get_Q(i,14,32,17) << endl;
	// 	Q_trap->SetBinContent(i+337,get_Q(i/10,14,32,17));
	// }

	// comment Q_trap filling when getting pos_hist

	// saving waveforms to a root file
	// create_wf();



	// E_0 is being converted to eV
	TH1F *main_histogram = waveform(E_0*1000, gauss_rms, sig_dec_time, count_rate, num_pulses, d);
	// Q_trap2->SetMarkerStyle(20);
	// Q_trap2->Draw("");
	// return main_histogram;


	// return 0;
	// old version (obtaining s_curve)
	// TH1F *s_hist = get_s_curve(main_histogram);

	// Number of counts request at half energy
	double noc_half_en = get_half_en_count(main_histogram, E_0);
	cout << noc_half_en << " noc_half_en\n";
	// delete hist;

	// return 0;

	//////////////////////////////////////////////////////////////////////////////////////////
	// new version

	// for (int i = Q_wnoise->GetNbinsX(); i > 0; i--)
	// {
	// 	en += Q_wnoise->GetBinContent(i);
	// 	// cout << en << endl;
	// 	en_wnoise->SetBinContent(i,en);
	// }

	// // en = 0;

	// double noph = en_wnoise->GetMaximum()*0.6;
	// double infp_val = en_wnoise->GetMaximum()*0.22;
	// double infp_bin;


	// // finding the inflection point bin
	// for (int j = en_wnoise->GetXaxis()->GetXmax(); j > en_wnoise->GetXaxis()->GetXmin(); j--)
	// {
	// 	double bin_content = en_wnoise->GetBinContent(j);

	// 	if (bin_content > infp_val)
	// 	{
	// 		infp_bin = j;
	// 		break;
	// 	}
	// }

	// // !
	// // cout << "\n////////////////////////////////////////\n";
	// // cout << "Dealing with simulation data\n";
	// // cout << "Parameters set as an initial guess\n\n";
	// // cout << "Number of photons:\t" << noph << endl;
	// // cout << "Inflection point bin:\t" << infp_bin << endl;
	// // cout << "\n////////////////////////////////////////\n\n";
	// // !

	// // Fitting simulation curve
	// Double_t mypar[6]={0,0,infp_bin,500,noph,0.01};
	// Double_t emypar[6]={0.1,0.1,0.1,0.1,0.1};




	// energyCalibrationFunctions* cal1 = new  energyCalibrationFunctions();
	// cal1->setScanSign(1);
	// TF1* s_sim_func  = new TF1("s_sim_func", cal1, 
	// 		&energyCalibrationFunctions::erfFunctionChargeSharing,
	// 		en_wnoise->GetXaxis()->GetXmin(),
	// 		en_wnoise->GetXaxis()->GetXmax(), 6,
	// 		"energyCalibrationFunctions",
	// 		"erfFunctionChargeSharing");

	// s_sim_func->SetParNames("Background Offset","Background Slope",
	//							"Inflection Point","Noise RMS",
	//							"Number of Photons","Charge Sharing Slope");
	// s_sim_func->SetParameters(mypar);

	// en_wnoise->Fit("s_sim_func", "R");

	// // double par_inf_point = s_sim_func->GetParameter(2);
	// // double par_num_ph = s_sim_func->GetParameter(4);

	// // double noc_sim = s_sim_func->Eval(par_inf_point/2);

	// // wrting rate data to a file
	// // FILE *f = fopen("rate.dat", "a");

	// // // fputs ("fopen example",f);
	// // fprintf(f, "%d\t%f\t%f\t%f\n", count_rate, par_num_ph, noc_sim, par_num_ph/noc_sim);

	// // fclose(f);

	//////////////////////////////////////////////////////////////////////////////////////////

	///////////////////////////////////////////////////

	// energyCalibration *cal = new energyCalibration();
	// cal->setScanSign(1);
	// TF1 *s_sim_func = cal->fitSCurve(en_wnoise, mypar, emypar);
	// // cal->fixParameter(0,0);
	// cal->fixParameter(1,0);
	// cal->setPlotFlag(0);

	// en_wnoise->Draw();
	// en_wnoise->Fit("fscurve");

	///////////////////////////////////////////////////




	////////////////////////////////////
	// double inf_point = s_sim_func->GetParameter(2);
	// double num_ph = s_sim_func->GetParameter(4);

	// double noc_sim = s_sim_func->Eval(inf_point/2);

	// cout << num_ph << " num_ph\n";

	// // cout << inf_point << "\t" << s_sim_func->Eval(inf_point/2) << " Eval val" << endl;

	// // int noc_bin_exp = hist_ch_new_reqen->FindBin((int)(par_inf_point))/2;
	// // cout << "noc_bin_exp " << noc_bin_exp << endl;
	// // double noc_exp = hist_ch_new_reqen->GetBinContent(noc_bin_exp);
	// // cout << "@@@ NOC_EXP " << noc_exp << " @@@\n";

	// d = (noc_sim-num_ph)/noc_sim*50;


	// cout << "\nInput charge sharing distance\n" << d << endl;

	// // return 0;

	// ////////////////////////////////////


	// for (int parnum = 0; parnum <= 5; parnum++)
	// 	cout << s_sim_func->GetParameter(parnum) << endl;


	// root file writing
	// TFile f("22_25.root", "recreate");
	// en_wnoise->Write();
	// s_sim_func->Write();
	// f.Close();

	// //////////////////////////////////////////////////////////////////////////////////////////


	// old version
	// double noph = s_hist->GetMaximum()*0.6;
	// double infp_val = s_hist->GetMaximum()*0.22;

	// double infp_bin;
	// // finding the inflection point bin
	// for (int j = s_hist->GetXaxis()->GetXmax(); j > s_hist->GetXaxis()->GetXmin(); j--)
	// {
	// 	double bin_content = s_hist->GetBinContent(j);

	// 	if (bin_content > infp_val)
	// 	{
	// 		infp_bin = j;
	// 		break;
	// 	}
	// }

	// return 0;


	// cout << "\n////////////////////////////////////////\n";
	// cout << "Dealing with simulation data\n";
	// cout << "Parameters set as an initial guess\n\n";
	// cout << "Number of photons:\t" << noph << endl;
	// cout << "Inflection point bin:\t" << infp_bin << endl;
	// cout << "\n////////////////////////////////////////\n\n";

	// old version
	// Fitting simulation curve
	// Double_t mypar[6]={0.1,0.001,infp_bin,1000,noph,0.01};
	// Double_t emypar[6]={0.1,0.1,0.1,0.1,0.1};

/*
	energyCalibration *cal = new energyCalibration();
	cal->setScanSign(1);
	TF1 *s_sim_func = cal->fitSCurve(s_hist, mypar, emypar);
	// cal->fixParameter(0,0);
	// cal->fixParameter(1,0);
	cal->setPlotFlag(0);
*/

	// for (int parnum = 0; parnum <= 5; parnum++)
	// 	cout << mypar[parnum] << endl;


	// old version
	// energyCalibrationFunctions* cal1 = new  energyCalibrationFunctions();
	// cal1->setScanSign(1);
	// TF1* s_sim_func  = new TF1("s_sim_func", cal1, 
	// 		&energyCalibrationFunctions::erfFunctionChargeSharing,
	// 		s_hist->GetXaxis()->GetXmin(),
	// 		s_hist->GetXaxis()->GetXmax(), 6,
	// 		"energyCalibrationFunctions",
	// 		"erfFunctionChargeSharing");

	// s_sim_func->SetParNames("Background Offset","Background Slope",
	// 						"Inflection Point","Noise RMS",
	// 						"Number of Photons","Charge Sharing Slope");
	// s_sim_func->SetParameters(mypar);

	// s_hist->Fit("s_sim_func", "R");



	// double par_inf_point = s_sim_func->GetParameter(2);
	// double par_num_ph = s_sim_func->GetParameter(4);

	// double noc_sim = s_sim_func->Eval(par_inf_point/2);


	// file writing
	
	// writing rate data to a file
	// FILE *f = fopen("rate_68.dat", "a");
	// fprintf(f, "%d\t%d\t%f\n", count_rate, pulse_count, noc_half_en);
	// fclose(f);


	//////////////////////////////////////////////////////////////////////////////////////////

	// old version
	// s_hist->SetLineColor(kGreen);
	// s_hist->Draw("SAME");

	// new version
	// en_wnoise->SetLineColor(kGreen);
	// en_wnoise->Draw("SAME");

	// main_histogram->Draw();

	// return en_wnoise;
	// return s_hist;
	// return main_histogram;
	return 0;
}