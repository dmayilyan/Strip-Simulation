#include <iostream>
using namespace::std;
#include "TCanvas.h"
#include "TF1.h"
#include "TH1.h"
#include "TH1F.h"
#include "TMath.h"
#include "math.h"
#include "TRandom.h"


#define E_0			10000			// Energy of the photon in eV
#define Q_0			E_0/3.6			// Number of electrons generated
#define TAU			80				// Given in ns

#define GAUS_RMS	240				// RMS of gaussian noise spread

// either n-photons or 
#define GATE		1000000000		// Gate size in ns
#define NUM_PULSES	1000			// Number of pulses

#define EXP_TAU		100				// Decay time in ns
#define COUNT_RATE	1				// Counting rate in kHz


double p = 50;
double d = 17;
double Eff = 0;

double Q = 0;

TH1I *hist = new TH1I("hist", "Signal data", GATE/2, 0, GATE);
TH1F *en_hist = new TH1F("en_hist", "Energy data", 2000, 0, 2000);

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


Double_t get_Q(Double_t pos)
{

	if ((pos >=d) && (pos <= p))
	{
		// Flat part
		Eff = Q_0;
	}
	else if (pos < d)
	{
		// Left slope
		Eff = Q_0 * pos/d;
	}
	else if ((pos > p) && (pos < p+d))
	{
		// Right slope
		Eff = Q_0 * (p+d-pos)/d;
	}
	else
		Eff = 0;

	return Eff;
}



Double_t exp_func(Double_t *x, Double_t *par)
{
	Float_t xx =x[0];

	Double_t f = (par[1] + par[2])/par[0] * xx * exp(-xx/par[0]);

	return f;
}

TH1I* get_exp(int tau_l, Double_t Eff)
{
	TH1I* loc_hist = new TH1I("loc_hist", "Single pulse", 1000/2, 0, 1000);

	// Double_t f = (par[1] + gaus_noise) * xx * exp(-xx/par[0]);
	Double_t gaus_noise = gRandom->Gaus(0,GAUS_RMS);

	TF1 *expo = new TF1("Exp curve", exp_func, 0, 1000, 3);
	expo->SetParameters(tau_l, Eff, gaus_noise);
	for (int ibin = 1; ibin < loc_hist->GetNbinsX()+1; ibin++)
		loc_hist->SetBinContent(ibin,expo->Eval(loc_hist->GetBinCenter(ibin)));

	expo->SetParNames("tau_l", "Q", "Gauss noise");
	en_hist->Fill(loc_hist->GetMaximum());

	return loc_hist;
}

void get_s_curve()
{
	TCanvas *c2 = new TCanvas("S Curve");


	TH1F* s_curve = new TH1F("s_curve", "S Curve", 10000/2, 0, 10000);

	int hist_val_old = 0;
	int hist_val_cur = 0;
	int hist_val_new = 0;


	int hist_max_val = hist->GetBinContent(hist->GetMaximumBin());

	cout << hist->GetMaximum() << " That's the max val" << endl;

	int *g_gauge_count = (int*)malloc(hist_max_val*sizeof(int));
	if(g_gauge_count==NULL)
	{
		printf("Error! memory not allocated.");
		exit(0);
	}

	cout << "hist_max_val " << hist_max_val << endl;

	int X_COUNT = hist->GetNbinsX();

	int *g_min = (int*)malloc(X_COUNT*sizeof(int));
	if(!g_gauge_count)
	{
		cout << "Error! memory not allocated.";
		exit(0);
	}

	int *g_max = (int*)malloc(X_COUNT*sizeof(int));
	if(!g_gauge_count)
	{
		cout << "Error! memory not allocated.";
		exit(0);
	}

	int min_count_new = 0;
	int max_count_new = 0;
	int g_count = 0;

	int j = 0;
	for (int i = 1; i < X_COUNT; i++)
	{

		hist_val_old = hist->GetBinContent(i-1);

		j = i;
		while((hist->GetBinContent(j) == hist_val_old) && (j < X_COUNT))
			j++;
		if (X_COUNT == j)
			break;
		hist_val_cur = hist->GetBinContent(j);

		i = j;

		// find next level
		j++;
		while((hist->GetBinContent(j) == hist_val_cur) && (j < X_COUNT))
			j++;
		if (X_COUNT == j)
			break;
		hist_val_new = hist->GetBinContent(j);


		// suppose first comes a peak
		if ((hist_val_old < hist_val_cur) && (hist_val_new < hist_val_cur) /*&& ((hist_val_old > 0) && (hist_val_new > 0))*/ )		// peak
		{
			g_max[g_count] = hist_val_cur;
			max_count_new++;
			// cout << "max found at bin " << hist->FindBin(i) << endl;
		}
		else
		if ((hist_val_cur < hist_val_old) && (hist_val_cur < hist_val_new) /*&& ((hist_val_old > 0) && (hist_val_new > 0))*/ )		// min
		{
			g_min[g_count] = hist_val_cur;
			g_count++;
			min_count_new++;
			// cout << i << " HIST VAL " << hist_val_new << endl;
			// cout << "min found at bin " << hist->FindBin(i) << endl;
		}

	}

	// cout << g_count << " @*&($#@!&*^$$)(*#@$(*&" << endl;
	// cout << "min_count_new " << min_count_new << " max_count_new " << max_count_new << endl;

	for (int i = 1; i < hist_max_val; i++)
	{
		int min_count = 0;
		int max_count = 0;
		for (int j = 0; j < g_count; j++)
		{
			if (g_min[j] >= i)
				min_count++;
			if (g_max[j] >= i)
				max_count++;
		}

		// g_gauge_count[i] = max_count - min_count;
		// cout << min_count << "\t" << max_count << endl;
		s_curve->SetBinContent(s_curve->FindBin(i), max_count - min_count);

	}


	s_curve->Draw();

	if (g_gauge_count)
		free(g_gauge_count);
	if (g_min)
		free(g_min);
	if (g_max)
		free(g_max);

}

int main()
{

	TH1I *loc_hist;

	// double time_spread_mean = GATE / COUNT_RATE / 1000; // 1000 is for conversion to Hz

	for (int i = 0; i < NUM_PULSES; i++)
	{
		double Charge_pos = 0;
		// double T = time_spread_mean - gRandom->Poisson(time_spread_mean);
		double T = /*gRandom->*/ TMath::Exp(EXP_TAU);
		int RP = gRandom->Rndm() * GATE/2;
		int SP = RP + T;
		// cout << RP << "\t" << T << endl;

		// double T = gRandom->Rndm() * GATE;

		// Version for 3 strips
		// Charge_pos = gRandom->Rndm() * (3*p+d);

		// Version with one strip
		Charge_pos = gRandom->Rndm() * (p+d);

		// cout << "\nSP " << SP*2 << "\t";
		double Q = get_Q(Charge_pos);

		loc_hist = get_exp(TAU, Q);

		for (int h_start = 0; h_start < loc_hist->GetNbinsX(); h_start++)
		{
			double start_pos = SP+h_start;
			if (start_pos >= GATE)
				continue;
			else
			{
				// Filling the bins with integers
				hist->SetBinContent(start_pos, (int)hist->GetBinContent(start_pos) + (int)loc_hist->GetBinContent(h_start));
			}
		}
		delete loc_hist;
		
	}

	// Uncumment to see the signal
	// hist->Draw();
	// Expanding signal out
	// c1->SaveAs("graph_shot.root");

	// Drawing the energy
	// en_hist->Draw();
	// Analysis of the data
	// Making the S-curve

	get_s_curve();

	return 0;

}
