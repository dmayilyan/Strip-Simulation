// #define EXP_MAX		1.0 	//  Landau curve surrounding box
#define Q_0			3000	// Number of electrons generated
#define tau		30/*1000*/	// Given in nm
#define GATE		10000	// Gate size
#define gaus_rms	240			// RMS of gaussian noise spread

#define COUNT_RATE	1		// Counting rate in kHz

// #define p/*STRIP_PITCH*/	50
// #define d/*DISTANCE*/	17		// Distance between the stripes

#include <iostream>

double p = 50;
double d = 17;

Double_t get_Q(Double_t pos)
{
	double Q;

	if ((pos >= -(p-d)/2) && (pos <= (p-d)/2)) || ((pos >= -(p-d)/2-p) && (pos <= (p-d)/2-p)) || ((pos >= -(p-d)/2+p) && (pos <= (p-d)/2+p))
	{
		Q = Q_0;
	}
	else if ((pos < -(p-d)/2-p) && (pos > -(p+d)/2-p))
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
	}
	else
		Q = 0;

	return Q;
}

// Double_t get_exp(Double_t start_point, Int_t event_count, Double_t tau, Double_t Q, TH1F* m)
// {
// 	Double_t t, y, x;
// 	Int_t count = 0;
// 	Int_t event_no = 0;
// 	// Big number to avoid infinite loop

// 	cout << "Start point " << start_point << endl;

// 	TF1 *exp_func = new TF1("exp_func","[0]*x*exp(-x/[1])",start_point,start_point+tau);
// 	TF1 *exp_1 = new TF1("exp_func","x*exp(-x/30)",0,100);
// 	TF1 *exp_2 = new TF1("exp_func","x*exp(-x/30)",50,150);

// 	exp_func->SetParNames("Q", "tau");
// 	exp_func->SetParameters(Q, tau);

// 	TH1F* histo=new TH1F("Hist data","Exp signal",GATE/2,0,GATE);
// 	// histo->Draw();
// 	// exp_func->Draw("same");

// 	// exp_1->Draw();
// 	// exp_2->Draw("+");


// 	// while ((event_no < event_count) && (count < 1000000))
// 	// {
// 	// 	t = gRandom->Rndm()* 1000/*(GATE-1)*/;	// Random point in the gate
// 	// 	// y = gRandom->Rndm()*EXP_MAX; 	// Random point in the Box
// 	// 	x = start_point + t;
// 	// 	// tau *= 2;
// 	// 	// Checking if the point is under the curve
// 	// 	// cout << "!!!! tau " << tau << endl;

// 	// 	if (y <= Q*t*exp(-t/tau))
// 	// 	{
// 	// 		event_no++;
// 	// 		if ((x < GATE) && (x >= 0))
// 	// 			m->Fill(x);				// FillingE histogram
// 	// 	}
// 	// 	count++;
// 	// }

// }

Double_t get_exp(Double_t y[], Double_t start_point, Double_t tau, Double_t Q)
{
	// cout << start_point << "\t" << tau << "\t" << Q << endl;

	for (int i = 0; i < 1000; ++i)
	{
		if (start_point + i > GATE-1)
			continue;
		else
		{
			double gaus_noise = gRandom->Gaus(0,gaus_rms);
			// cout << "gn " << gaus_noise << endl;
			y[start_point+i] += Q * i * exp(-i/tau) + gaus_noise;
		// cout << y[i] << endl;
		}
	}

	return y;
}

void main()
{
	double x[GATE] = {0,};
	double y[GATE] = {0,};

	for (int i = 0; i < GATE; ++i)
		x[i] = i;

	for (int i = 0; i < GATE; ++i)
		y[i] = 0;

	cout << "-(p+d)/2 " << -(p+d)/2 << "\t" << "-(p-d)/2 " << -(p-d)/2 << endl;
	cout << "-(p-d)/2 " << -(p-d)/2 << "\t" << "(p-d)/2 " << (p-d)/2 << endl;
	cout << "(p-d)/2 " << (p-d)/2 << "\t" << "(p+d)/2 " << (p+d)/2 << endl << endl;

	double time_spread_mean = GATE / COUNT_RATE / 1000; // 1000 is for conversion to Hz
	cout << "Time spread mean " << time_spread_mean << endl << endl;

	// TCanvas *c1 = new TCanvas();

	TH1F* m=new TH1F("Hist data","Exp signal",GATE/2,0,GATE); // 2 ns sampling
	// TH1F* l=new TH1F("qwe", "qwe",100,-50,50);

	// for (int i = 0; i < 100; ++i)
	// {
	// 	// TRandom *R = new TRandom(time(0));  // create a pointer to a new instance of TRandom
	// 	l->Fill(10 - gRandom->Poisson(10));
	// 	cout << 10 - gRandom->Poisson(10) << endl;
	// }

	for (int i = 0; i < 10; i++)
	{
		// double T = gRandom->Poisson(10) /** GATE*/;
		double Charge_pos = 0;
		double T = time_spread_mean - gRandom->Poisson(time_spread_mean) /** GATE*/;
		double RP = gRandom->Rndm() * GATE;
		double SP = RP + T;
		// cout << RP << "\t" << T << endl;

		// double T = gRandom->Rndm() * GATE;


		Charge_pos = gRandom->Rndm() * 3*(p+d)-(p+d)*3/2;
		// Charge_pos = gRandom->Rndm() * 3*p+4*d - (3*p+4*d)/2;


		double Q = get_Q(Charge_pos);

		// cout << "Charge_pos " << Charge_pos;
		// cout << " Q " << Q << endl;


		// get_exp(SP, 1000, 1000/*tau*/, Q, m);
		// Start time, number of points, tau, Number of phot_el, histogram

		get_exp(y, SP, tau, Q);
		// X coordinate, Y values, Start time, Decay time, Charge collection efficiency
		
	}


	gr = new TGraph(GATE,x,y);

	gr->Draw();

	gr->GetXaxis()->SetRangeUser(0,GATE);
	// gr->GetYaxis()->SetRangeUser(0,35000);

	// m->Draw();
	// l->Draw();


}