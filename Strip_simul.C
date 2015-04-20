// #define EXP_MAX		1.0 	//  Landau curve surrounding box
#define E_0			10000		// Energy of the photon in eV
#define Q_0			E_0/3.6		// Number of electrons generated
#define tau			30/*1000*/	// Given in nm
#define GATE		10000	// Gate size
#define gaus_rms	240			// RMS of gaussian noise spread

#define COUNT_RATE	1		// Counting rate in kHz

#include <iostream>
#include "TF1"

double p = 50;
double d = 17;
double Eff = 0;

double Q = 0;

TH1F *hist = new TH1F("hist", "Signal data", GATE/2, 0, GATE);

Double_t get_Q(Double_t pos)
{
	// double Q_0 = 3000;	// Number of electrons generated


	if (((pos >=d) && (pos <= p)) || ((pos >=d+p) && (pos <= 2*p)) || ((pos >=d+2*p) && (pos <= 3*p)))
	{
		Q = Q_0;
		Eff = Q_0;
		// cout << pos << " Middle region; Q = " << Eff << endl;
	}
	else if (pos < d)
	{
		// /
		Eff = Q_0 * pos/d;
		// cout << pos << " Left slope region; Q = " << Eff << endl;
	}
	else if ((pos > 3*p) && (pos < (3*p + d)))
	{
		Q = Q_0 * (3*p+d-pos)/d;
		Eff = Q_0 * (3*p+d-pos)/d;
		// cout << pos << " Right slope region; Q = " << Eff << endl;
		// cout << (3*p+d-pos)/d << " !!!!!!!!\n" << endl;
	}
	else if ((pos > p) && (pos < p+d))
	{
		Q = Q_0 * (p+d-pos)/d + Q_0 * (pos-p)/d;
		Eff = Q_0 * (p+d-pos)/d + Q_0 * (pos-p)/d;
		// cout << pos << " First sharing region; Q = " << Eff << endl;
	}
	else if ((pos > 2*p) && (pos < 2*p + d))
	{
		Q =  Q_0 * (2*p+d-pos)/d + Q_0 * (pos-2*p)/d;
		Eff =  Q_0 * (2*p+d-pos)/d + Q_0 * (pos-2*p)/d;
		// cout << pos << " Second sharing region; Q = " << Eff << endl;
	}
	else
		Eff = 0;

	return Eff;
}

Double_t exp_func(Double_t *x, Double_t *par)
{
	Float_t xx =x[0];

	Double_t f = (par[1] + par[2]) * xx * exp(-xx/par[0]);

	return f;
}

TH1F * get_exp(Double_t tau, Double_t Eff)
{
	TH1F* loc_hist = new TH1F("loc_hist", "Single pulse", 1000/2, 0, 1000);

	// cout << start_point << "\t" << tau << "\t" << Q << endl;

	// Double_t gaus_noise = gRandom->Gaus(0,gaus_rms);
	// Double_t f = (par[1] + gaus_noise) * xx * exp(-xx/par[0]);
	Double_t gaus_noise = gRandom->Gaus(0,gaus_rms);

	TF1 *expo = new TF1("Exp curve", exp_func, 0, 1000, 3);
	expo -> SetParameters(tau, Eff, gaus_noise);
	for (int ibin = 1; ibin < loc_hist->GetNbinsX()+1; ibin++)
	{
		loc_hist->SetBinContent(ibin,expo->Eval(loc_hist->GetBinCenter(ibin)));
		// cout << expo->Eval(loc_hist->GetBinCenter(ibin)) << endl;
	}

	expo -> SetParNames("tau", "Q", "Gauss noise");
	// loc_hist -> Draw();

	return loc_hist;
}

void get_s_curve()
{
	double max = hist->GetMaximum();
	double max_next = hist->GetMaximum(max);
	double maxbin = hist->GetMaximumBin();
	// cout << "MAX " << max << " max_next " << max_next << " bin " << maxbin << endl;
	int layer_old[5000] = {0,};
	int layer[5000] = {0,};
	int column_old[100000] = {0,};
	int column[100000] = {0,};

	double peakcount[100000] = {0,};


	TH1F* s_curve = new TH1F("s_curve", "S Curve", 100000/2, 0, 100000);

	double hist_val_old = 0.;
	double hist_val = 0.;

	hist_val_old = hist->GetBinContent(0);
	for (int i = 0; i < hist->GetNbinsX(); i++)
	{
		hist_val = hist->GetBinContent(i);
		// cout << i << " hist_val " << hist_val << endl;
		// cout << i << " hist_val_old " << hist_val_old << endl;
		// hist_val_old = hist_val;

		if (hist_val >= hist_val_old)
		{
			hist_val_old = hist_val;
			// cout << "! " << i << endl;
			continue;
		}
		else
		{
			hist_val_old = hist_val;
			for (int j = hist_val-1; j >= 0; j--)
			{
				peakcount[j] += 1;
				double s_content = s_curve->GetBinContent(j)+1;
				s_curve->SetBinContent(j,s_content);
				// cout << "s_content " << s_content << endl;
			}
		}

	}

	s_curve->Draw();

/*
	for (int i = 0; i < hist->GetNbinsX(); i++)
	{
		for (int j = 0; j < 10000; j++)
		{
			if (hist->GetBinContent(i) > j)
			{
				column[j] += 1;	
				if (j!=0)
					if (column[j-1] == 1)
						column[j-1] = 0;
			}

			double s_content = s_curve->GetBinContent(j)+column[j];
			// cout << "REOIREOIJOU" << s_content << endl;
			s_curve->SetBinContent(j,s_content);
		}
	}

	s_curve -> Draw();
*/

	// // bool under_curve = true;
	// for (double i = max; i > max-400; i--)
	// {
	// 	for (int j = 0; j < hist->GetNbinsX(); j++)
	// 	{
	// 		// Not sensitive to the case when first bin is non-zero
	// 		double hist_val[j] = hist->GetBinContent(j);
	// 		double hist_val_prev[j] = {0,};

	// 		if (hist_val[j] > i-1)
	// 		{
	// 			if (hist_val_prev[j] >= hist)
	// 			cout << "Քաքը մեծ ա" << endl;
	// 		}
	// 		// Filling layer with peak positions
	// 		// cout << hist->GetNbinsX() << " &&&&&&&&&&&" << endl;
	// 		// if (hist_val > i-1)
	// 		// {
	// 		// 	cout << "պուպուլ " << j*2 << " val " << hist_val << endl;
	// 		// 	// under_curve = true;
	// 		// 	layer[j] = 1;
	// 		// 	if (i == max)
	// 		// 		layer_old[j] = layer[j];

	// 		// 	if (j != 0)
	// 		// 		if ((layer[j] == 1) && (layer[j-1] == 0))
	// 		// 		{
	// 		// 			peakcount++;
	// 		// 			cout << "j " << 2*j << endl;
	// 		// 		}
				
	// 		// }
	// 		// else
	// 		// {
	// 		// 	// under_curve = false;
	// 		// 	layer[j] = 0;
	// 		// 	if (i == max)
	// 		// 		layer_old[j] = layer[j];
	// 		// 	// cout << "chka " << j << endl;
	// 		// }

	// 		// cout << layer[j] << "\t" << layer_old[j] << endl;

	// 		// if ((layer_old[j] <= layer[j]) && (layer_old[j] >= layer[j]))
	// 		// 	if ((layer[j-1] == 0) && (layer[j] == 1))

	// 	}

	// 	// for (int j = 0; j < hist->GetNbinsX(); j++)
	// 	// {
	// 	// 	layer_old[j] = layer[j];

	// 	// 	if (layer_old[j] == layer[j])
	// 	// 		continue;
	// 	// 	else
	// 	// 		peakcount++;
	// 	// }
	// 	cout << peakcount << " #####################" << endl;
	// }
	
}

void main()
{
	double x_graph[GATE] = {0,};
	double y[GATE] = {0,};

	for (int i = 0; i < GATE; ++i)
		x_graph[i] = i;

	for (int i = 0; i < GATE; ++i)
		y[i] = 0;

	// cout << "-(p+d)/2 " << -(p+d)/2 << "\t" << "-(p-d)/2 " << -(p-d)/2 << endl;
	// cout << "-(p-d)/2 " << -(p-d)/2 << "\t" << "(p-d)/2 " << (p-d)/2 << endl;
	// cout << "(p-d)/2 " << (p-d)/2 << "\t" << "(p+d)/2 " << (p+d)/2 << endl << endl;

	double time_spread_mean = GATE / COUNT_RATE / 1000; // 1000 is for conversion to Hz
	// cout << "Time spread mean " << time_spread_mean << endl << endl;

	// TCanvas *c1 = new TCanvas();

	for (int i = 0; i < 50; i++)
	{
		// double T = gRandom->Poisson(10) /** GATE*/;
		double Charge_pos = 0;
		double T = time_spread_mean - gRandom->Poisson(time_spread_mean);
		int RP = gRandom->Rndm() * GATE/2;
		int SP = RP + T;
		// cout << RP << "\t" << T << endl;

		// double T = gRandom->Rndm() * GATE;

		Charge_pos = gRandom->Rndm() * (3*p+d);
		// Charge_pos = gRandom->Rndm() * 3*p+4*d - (3*p+4*d)/2;

		// cout << "\nSP " << SP*2 << "\t";
		double Q = get_Q(Charge_pos);

		// get_exp(SP, 1000, 1000/*tau*/, Q, m);
		// Start time, number of points, tau, Number of phot_el, histogram

		get_exp(tau, Q);

		// cout << "loc_hist " << loc_hist->GetNbinsX() << endl;
		// cout << "hist " << hist->GetNbinsX() << endl;
		// cout << "SP " << SP << endl;
		for (int h_start = 0; h_start < loc_hist->GetNbinsX(); h_start++)
		{
			start_pos = SP+h_start;
			if (start_pos >= GATE)
				continue;
			else
			{
				// Filling the bins with integers
				hist->SetBinContent(start_pos, (int)hist->GetBinContent(start_pos) + (int)loc_hist->GetBinContent(h_start));
				// cout << h_start << "\t" << loc_hist->GetBinContent(h_ssart) << endl;
			}
		}
		
	}

	hist->Draw();
	c1->SaveAs("graph_shot.root");

	// Analysis of the data
	// Making the S-curve

	get_s_curve();


}