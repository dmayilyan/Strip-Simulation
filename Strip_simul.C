// #define EXP_MAX		1.0 	//  Landau curve surrounding box
#define Q_0			3000	// Number of electrons generated
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
		cout << pos << " Middle region; Q = " << Eff << endl;
	}
	else if (pos < d)
	{
		// /
		Eff = Q_0 * pos/d;
		cout << pos << " Left slope region; Q = " << Eff << endl;
	}
	else if ((pos > 3*p) && (pos < (3*p + d)))
	{
		Q = Q_0 * (3*p+d-pos)/d;
		Eff = Q_0 * (3*p+d-pos)/d;
		// cout << "!!!!!!!!!!! " << Eff << " !!!!!!!!!!!!!!" << endl;
		cout << pos << " Right slope region; Q = " << Eff << endl;
		cout << (3*p+d-pos)/d << " !!!!!!!!\n" << endl;
	}
	else if ((pos > p) && (pos < p+d))
	{
		Q = Q_0 * (p+d-pos)/d + Q_0 * (pos-p)/d;
		Eff = Q_0 * (p+d-pos)/d + Q_0 * (pos-p)/d;
		// cout << "!!!!!!!!!!! " << Eff << " !!!!!!!!!!!!!!" << endl;
		cout << pos << " First sharing region; Q = " << Eff << endl;
	}
	else if ((pos > 2*p) && (pos < 2*p + d))
	{
		Q =  Q_0 * (2*p+d-pos)/d + Q_0 * (pos-2*p)/d;
		Eff =  Q_0 * (2*p+d-pos)/d + Q_0 * (pos-2*p)/d;
		cout << pos << " Second sharing region; Q = " << Eff << endl;
	}


/*	else if (((pos > d/2) && (pos < d)) || ((pos > d/2 + p) && (pos < d + p)) || ((pos > d/2 + 2*p) && (pos < d + 2*p)))
	{
		Q = Q_0 * (2*pos+p+d)/(2*d);
	}
	else if ((pos > (p-d)/2+p) && (pos < (p+d)/2+p))
	{
		Q = Q_0 * (p+d-2*pos)/(2*d);
		cout << "pos1 " << (p-d)/2 << endl;
		cout << "pos2 " << -(p-d)/2+p << endl;
	}
	else if ((pos > (p-d)/2-p) && (pos < -(p-d)/2))
	{
		Q = Q_0 * (2*pos+p+d)/(2*d) + Q_0 * (p+d-2*pos)/(2*d);

		cout << "Contribution from first stripe " << Q_0 * (2*pos+p+d)/(2*d) << endl;
		cout << "Contribution from second stripe " << Q_0 * (p+d-2*pos)/(2*d) << endl;
		cout << "Overall contribution " << Q << endl;

		cout << "First stripe val " << (2*pos+p+d)/(2*d) << endl;
		cout << "Second stripe val " << (p+d-2*pos)/(2*d) << endl;
	}
	else if ((pos > (p-d)/2) && (pos < -(p-d)/2+p))
	{
		// Q = Q_0 * (p+d-2*pos)/(2*d);
		Q = Q_0 * (2*pos+p+d)/(2*d) + Q_0 * (p+d-2*pos)/(2*d);

		cout << "Contribution from second stripe " << Q_0 * (p+d-2*pos)/(2*d) << endl;
		cout << "Contribution from third stripe " << Q_0 * (2*pos+p+d)/(2*d) << endl;
		cout << "Overall contribution " << Q << endl;
	}*/
	else
		Eff = 0;

	return Eff;
}

Double_t qwe(Double_t *x, Double_t *par)
{
	Float_t xx =x[0];

	// cout << "gn = " << par[2] <<  endl;
	Double_t f = (par[1] + par[2]) * xx * exp(-xx/par[0]);

	// Double_t val = ()
	return f;
}

TH1F * get_exp(Double_t tau, Double_t Eff)
{
	TH1F* loc_hist = new TH1F("loc_hist", "Single pulse", 1000/2, 0, 1000);

	// cout << start_point << "\t" << tau << "\t" << Q << endl;

	// Double_t gaus_noise = gRandom->Gaus(0,gaus_rms);
	// Double_t f = (par[1] + gaus_noise) * xx * exp(-xx/par[0]);
	Double_t gaus_noise = gRandom->Gaus(0,gaus_rms);

	TF1 *expo = new TF1("Exp curve", qwe, 0, 1000, 3);
	expo -> SetParameters(tau, Eff, gaus_noise);
	for (int ibin = 1; ibin < loc_hist->GetNbinsX()+1; ibin++)
	{
		loc_hist->SetBinContent(ibin,expo->Eval(loc_hist->GetBinCenter(ibin)));
		// cout << expo->Eval(loc_hist->GetBinCenter(ibin)) << endl;
	}

	expo -> SetParNames("tau", "Q");
	// loc_hist -> Draw();

	// (expo->GetHistogram())->Draw();


	// double gaus_noise = gRandom->Gaus(0,gaus_rms);
	// cout << "Start point " << start_point << " gn " << gaus_noise << " Eff. + Noise " << Eff + gaus_noise << endl;
	// for (int i = 0; i < 1000; ++i)
	// {
	// 	if (start_point + i > GATE-1)
	// 		continue;
	// 	else
	// 	{
	// 		y[start_point+i] += (Eff + gaus_noise) * i * exp(-i/tau);
	// 	// cout << y[i] << endl;
	// 	}
	// }

	return loc_hist;
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

	// TH1F* m=new TH1F("Hist data","Exp signal",GATE/2,0,GATE); // 2 ns sampling
	// TH1F* l=new TH1F("qwe", "qwe",100,-50,50);

	// for (int i = 0; i < 100; ++i)
	// {
	// 	// TRandom *R = new TRandom(time(0));  // create a pointer to a new instance of TRandom
	// 	l->Fill(10 - gRandom->Poisson(10));
	// 	cout << 10 - gRandom->Poisson(10) << endl;
	// }

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

		cout << "loc_hist " << loc_hist->GetNbinsX() << endl;
		cout << "hist " << hist->GetNbinsX() << endl;
		cout << "SP " << SP << endl;
		for (int h_start = 0; h_start < loc_hist->GetNbinsX(); h_start++)
		{
			start_pos = SP+h_start;
			if (start_pos >= GATE)
				continue;
			else
			{
				hist->SetBinContent(start_pos, hist->GetBinContent(start_pos) + loc_hist->GetBinContent(h_start));
				// cout << h_start << "\t" << loc_hist->GetBinContent(h_ssart) << endl;
			}
		}
		// cout << "START POIINT " << SP << endl;

		// delete loc_hist;

		// Y values, Start time, Decay time, Charge collection efficiency
		
	}

	hist->Draw();
	c1->SaveAs("graph_shot.root");

	// Analysis of the data
	


}