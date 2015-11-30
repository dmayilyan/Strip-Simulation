#define AV_POINTS 3		//  Should be odd


#include "filter.C"

void ef_scan3_scan(double thr0 = 0.5, double thr1 = 1.5, double thr_low = 2.5, double thr_high = 3.5, double step = 0.1)
{
	int E_0 = 14;
	// double x[500] = {0,};
	// double y[500] = {0,};
	double x,y;

	double thr[500] = {0,};
	// thr[0] = thr0;
	// thr[1] = thr1;

	double x_sum[500] = {0,};
	double y_sum[500] = {0,};

	// making the root filename fot given energy
	char filename[32];
	snprintf(filename, sizeof(filename), "rate_scan_%d_cad.root", E_0);


	// opening the rate_scan tree
	TFile *f = new TFile(filename);
	TTree *Tr = (TTree*)f->Get("Tr");
	int rate, n_phot;
	TH1F *s_hist;
	Tr->SetBranchAddress("rate_c",&rate);
	Tr->SetBranchAddress("n_phot_c",&n_phot);
	Tr->SetBranchAddress("s_hist_c",&s_hist);

	int ch_sh;
	if (E_0 == 14)
		ch_sh = 0;
	else if (E_0 == 18)
		ch_sh = 1;
	else if (E_0 == 22)
		ch_sh = 2;


	double charge_sh[3] = {24.08,23.55,23.71};
	// Collected charge correction
	double coll_corr = 50./(50.+charge_sh[ch_sh]);

	TMultiGraph *mg = new TMultiGraph();


	int num_thrs = 0;
	//filling the threshold array
	for (double th_i = thr_low; th_i < thr_high+0.01; th_i += step)
	{
		thr[num_thrs] = th_i;
		cout << num_thrs << "\t!!! " << th_i << " !!!\t" << thr[num_thrs] << "\tqweqwe\t" << thr_high << endl;
		// cout << thr[0] << "\t" << thr[1] << "\t" << thr[2] << "\t" << thr[3] << "\t" << thr[4] << "\t" << thr[5] << "\t" << thr[6] << "\n"; 
		++num_thrs;
	}

	// cout << num_thrs << "\tthr_count !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";

	// Making graph array
	char graphname1[16];
	char graphname2[16];
	TGraphErrors *ef_graph[50];
	TGraphErrors *ef_sum[50];

	for (int i = 0; i <= num_thrs; i++)
	{
		sprintf(graphname1,"ef_graph_%d",i);
		ef_graph[i] = new TGraphErrors();

		sprintf(graphname2,"ef_sum_%d",i);
		ef_sum[i] = new TGraphErrors();
	}

	char all_name[16];
	sprintf(all_name,"th =");
	char sum_name[128];
	snprintf(sum_name, sizeof(sum_name), "th = %1.2f, %1.2f, %1.2f-%1.2f, in %1.2f step", thr0, thr1, thr_low, thr_high, step);
	Int_t nentries = (Int_t)Tr->GetEntries();
	for (int entry_num = 0; entry_num < nentries; entry_num++)
	{
		// cout << thr[entry_num] << " 213123123123" << endl;
		///////////////////////////////////////
		// First threshold

		// finding efficiency for thr0(0.5) to not repeat it on every threshold pair call
		Tr->GetEntry(entry_num);

		// getting number of counts at given threshold
		// double thr0 = 0.5;
		int noc_bin = s_hist->FindBin((int)(E_0*1000))*thr0;
		double noc = s_hist->GetBinContent(noc_bin);

		n_phot = n_phot*coll_corr;

		double y0 = noc/n_phot;

		x = rate*coll_corr;

		///////////////////////////////////////
		// Second threshold

		// finding efficiency for thr0(0.5) to not repeat it on every threshold pair call
		Tr->GetEntry(entry_num);

		// getting number of counts at given threshold
		// double thr1 = 1.5;
		int noc_bin = s_hist->FindBin((int)(E_0*1000))*thr1;
		double noc = s_hist->GetBinContent(noc_bin);

		n_phot = n_phot*coll_corr;

		double y1 = noc/n_phot;


		Tr->GetEntry(entry_num);
		// cout << "rate " << rate << endl;

		n_phot = n_phot*coll_corr;

		for (int t = 0; t < num_thrs; t++)
		{
			// cout << "Scanning threshold " << thr0 << ", " << thr1 << ", and " << thr[t] << " from range [" << thr_low << "," << thr_high << "] in step of " << step << ".\n";
			// cout << "I AM HERE\n\n\n";

			
			// getting number of counts at given threshold
			int noc_bin = s_hist->FindBin((int)(E_0*1000))*thr[t];
			double noc = s_hist->GetBinContent(noc_bin);

			// cout << "Number of count for " << rate << " is " << noc << " for threshold of " << thr[t] << endl;

			x = rate*coll_corr*1000;
			y = noc/n_phot;

			// if (entry_num == nentries-1)
			cout << t+1 << "\t" << thr[t] << "\t" << rate << "\ty0:\t" << y0 << " y1:\t" << y1 << " y+y0+y1:\t" << y+y0+y1 << endl;
			ef_graph[t]->SetPoint(entry_num,x,y+y0+y1);

			mg->Add(ef_graph[t]);
		}

	}

	for (int i = 0; i <= num_thrs; i++)
	{
		if (i == 0)
		{
			ef_graph[i]->Draw("A*");
			ef_graph[i]->GetXaxis()->SetTitle("Rate [Hz]");
			ef_graph[i]->GetYaxis()->SetTitle("Efficiency");
			ef_graph[i]->GetXaxis()->SetTitleOffset(1.2);
			ef_graph[i]->GetYaxis()->SetTitleOffset(1.1);
			ef_graph[i]->SetTitle(sum_name);


		}

		ef_graph[i]->Draw("SAME *");
		if (i < 3)
			ef_graph[i]->SetMarkerColor(kBlue-i);
		else if (i < 6)
			ef_graph[i]->SetMarkerColor(kGreen-i%3);
		else if (i < 9)
			ef_graph[i]->SetMarkerColor(kRed-i%3);
		else if (i < 12)
			ef_graph[i]->SetMarkerColor(kCyan-i%3);

	}

	c1->SetLogx();
	c1->SetGridx();
	c1->SetGridy();

	// ef_graph[0]->GetYaxis()->SetRangeUser(0,1.075);

/*	TGraphErrors *exp_eff = run();

	// mg->Draw("a* fb l3d");

	TCanvas *c2 = new TCanvas("Sum of Effs");
	ef_sum[0]->Draw("A*");
	ef_graph[0]->Draw("SAME *");
	ef_sum[1]->Draw("SAME *");
	ef_sum[2]->Draw("SAME *");
	exp_eff->Draw("SAME *");

	ef_sum[0]->GetXaxis()->SetTitle("Rate");


	ef_sum[0]->SetMarkerColor(kBlue+1);
	ef_graph[0]->SetMarkerColor(kRed+1);
	ef_sum[1]->SetMarkerColor(kGreen+1);
	ef_sum[2]->SetMarkerColor(kBlack);
	exp_eff->SetMarkerColor(kOrange+7);

 	c2->SetLogx();
	c2->SetGridx();
	c2->SetGridy();
	ef_sum[0]->SetTitle(sum_name);
	ef_sum[0]->GetYaxis()->SetRangeUser(0,1.1);
	// ef_sum[0]->GetYaxis()->SetRangeUser(0.8,1.1);

	mg->Draw();*/
	// mg->Draw();
}


void ef_scan3(double thr0 = 0.5, double thr1 = 1.5, double thr2 = 2.5)
{
	int E_0 = 14;
	double thr[3] = {thr0,thr1,thr2};
	double y0[500] = {0,};

	double x_sum[500] = {0,};
	double y_sum[500] = {0,};

	// making the root filename fot given energy
	char filename[32];
	snprintf(filename, sizeof(filename), "rate_scan_%d_cad.root", E_0);


	// opening the rate_scan tree
	TFile *f = new TFile(filename);
	TTree *Tr = (TTree*)f->Get("Tr");
	int rate, n_phot;
	TH1F *s_hist;
	Tr->SetBranchAddress("rate_c",&rate);
	Tr->SetBranchAddress("n_phot_c",&n_phot);
	Tr->SetBranchAddress("s_hist_c",&s_hist);

	int ch_sh;
	if (E_0 == 14)
		ch_sh = 0;
	else if (E_0 == 18)
		ch_sh = 1;
	else if (E_0 == 22)
		ch_sh = 2;


	double charge_sh[3] = {24.08,23.55,23.71};
	// Collected charge correction
	double coll_corr = 50./(50.+charge_sh[ch_sh]);

	TMultiGraph *mg = new TMultiGraph();

	// Making graph array
	char graphname1[16];
	char graphname2[16];
	TGraphErrors *ef_graph[10];
	TGraphErrors *ef_sum[10];
	for (int i = 0; i < 3; i++)
	{
		sprintf(graphname1,"ef_graph_%d",i);
		ef_graph[i] = new TGraphErrors();

		sprintf(graphname2,"ef_sum_%d",i);
		ef_sum[i] = new TGraphErrors();
	}


	char all_name[16];
	sprintf(all_name,"th =");
	char sum_name[32];
	sprintf(sum_name,"th = 0.5, 0.5 +");
	Int_t nentries = (Int_t)Tr->GetEntries();
	for (int t = 0; t < 3; t++)
	{
		for (int i = 0; i < nentries; i++)
		{
			double x[500] = {0,};
			double y[500] = {0,};
			double x_err[500] = {0,};
			double y_err[500] = {0,};

			Tr->GetEntry(i);
			
			// getting number of counts at given threshold
			int noc_bin = s_hist->FindBin((int)(E_0*1000))*thr[t];
			double noc = s_hist->GetBinContent(noc_bin);

			// cout << "Number of count for " << rate << " is " << noc << " for threshold of " << thr[t] << endl;


			n_phot = n_phot*coll_corr;

			x[i] = rate*coll_corr*1000;
			x_sum[i] = x[i];
			y[i] = noc/n_phot;
			// x_err[i] = sqrt(rate)*coll_corr;
			y_err[i] = sqrt(1/noc+1/n_phot);

			// if (y[i] == 0) y_err[i] = 0.2;
			// if (t == 3) cout << y[i] << endl;

			ef_graph[t]->SetPoint(i,x[i],y[i]);
			ef_graph[t]->SetPointError(i,0,y_err[i]);

			if (t != 0)
				ef_sum[t-1]->SetPoint(i,x[i],y0[i] + y[i]);
			else
			{
				y0[i] = y[i];
				y_sum[i] = y[i];
			}

			if (t == 1)
				y_sum[i] += y[i];
			else if (t == 2)
			{
				y_sum[i] += y[i];
				ef_sum[t]->SetPoint(i,x_sum[i], y_sum[i]);
			}
		}
		// Generating graph name
		char loc_name[16];
		sprintf(loc_name,"th = %1.2f",thr[t]);

		if (t != 2)
			sprintf(all_name + strlen(all_name)," %1.2f,",thr[t]);
		else
			sprintf(all_name + strlen(all_name)," %1.2f",thr[t]);

		if (t == 1)
			sprintf(sum_name + strlen(sum_name)," %1.2f,",thr[t]);
		if (t == 2)
			sprintf(sum_name + strlen(sum_name)," %1.2f",thr[t]);

		mg->Add(ef_graph[t]);
		ef_graph[t]->SetTitle(loc_name);
	}

	TGraphErrors *exp_eff = run();

/*
	ef_graph[0]->Draw("A*");
	ef_graph[1]->Draw("SAME *");
	ef_graph[2]->Draw("SAME *");

	ef_graph[0]->GetXaxis()->SetTitle("Rate");
	ef_graph[0]->GetYaxis()->SetTitle("Efficiency");

	ef_graph[0]->SetMarkerColor(kRed+1);
	ef_graph[1]->SetMarkerColor(kBlue+1);
	ef_graph[2]->SetMarkerColor(kGreen+1);

	c1->SetLogx();
	c1->SetGridx();
	c1->SetGridy();
	ef_graph[0]->SetTitle(all_name);
	ef_graph[0]->GetYaxis()->SetRangeUser(0,1.075);
*/
	// mg->Draw("a* fb l3d");

	TCanvas *c2 = new TCanvas("Sum of Effs");
	ef_sum[0]->Draw("A*");
	ef_graph[0]->Draw("SAME *");
	ef_sum[1]->Draw("SAME *");
	ef_sum[2]->Draw("SAME *");
	exp_eff->Draw("SAME *");

	ef_sum[0]->GetXaxis()->SetTitle("Rate");


	ef_sum[0]->SetMarkerColor(kBlue+1);
	ef_graph[0]->SetMarkerColor(kRed+1);
	ef_sum[1]->SetMarkerColor(kGreen+1);
	ef_sum[2]->SetMarkerColor(kBlack);
	exp_eff->SetMarkerColor(kOrange+7);

	c2->SetLogx();
	c2->SetGridx();
	c2->SetGridy();
	ef_sum[0]->SetTitle(sum_name);
	ef_sum[0]->GetYaxis()->SetRangeUser(0,1.1);
	// ef_sum[0]->GetYaxis()->SetRangeUser(0.8,1.1);
}

void ef_scan4(double thr0 = 0.5, double thr1 = 1.5, double thr2 = 2.5, double thr3 = 3.5)
{
	int E_0 = 14;
	double thr[4] = {thr0,thr1,thr2,thr3};
	double y0[500] = {0,};

	double x_sum[500] = {0,};
	double y_sum[500] = {0,};

	// making the root filename fot given energy
	char filename[32];
	snprintf(filename, sizeof(filename), "rate_scan_%d.root", E_0);


	// opening the rate_scan tree
	TFile *f = new TFile(filename);
	TTree *Tr = (TTree*)f->Get("Tr");
	int rate, n_phot;
	TH1F *s_hist;
	Tr->SetBranchAddress("rate",&rate);
	Tr->SetBranchAddress("n_phot",&n_phot);
	Tr->SetBranchAddress("s_hist",&s_hist);

	int ch_sh;
	if (E_0 == 14)
		ch_sh = 0;
	else if (E_0 == 18)
		ch_sh = 1;
	else if (E_0 == 22)
		ch_sh = 2;


	double charge_sh[3] = {24.08,23.55,23.71};
	// Collected charge correction
	double coll_corr = 50./(50.+charge_sh[ch_sh]);

	TMultiGraph *mg = new TMultiGraph();

	// Making graph array
	char graphname1[16];
	char graphname2[16];
	TGraphErrors *ef_graph[10];
	TGraphErrors *ef_sum[10];
	for (int i = 0; i < 4; i++)
	{
		sprintf(graphname1,"ef_graph_%d",i);
		ef_graph[i] = new TGraphErrors();

		sprintf(graphname2,"ef_sum_%d",i);
		ef_sum[i] = new TGraphErrors();
	}


	char all_name[16];
	sprintf(all_name,"th =");
	char sum_name[32];
	sprintf(sum_name,"th = 0.5, 0.5 +");
	Int_t nentries = (Int_t)Tr->GetEntries();
	for (int t = 0; t < 4; t++)
	{
		for (int i = 0; i < nentries; i++)
		{
			double x[500] = {0,};
			double y[500] = {0,};
			double x_err[500] = {0,};
			double y_err[500] = {0,};

			Tr->GetEntry(i);
			
			// getting number of counts at given threshold
			int noc_bin = s_hist->FindBin((int)(E_0*1000))*thr[t];
			double noc = s_hist->GetBinContent(noc_bin);

			// cout << "Number of count for " << rate << " is " << noc << " for threshold of " << thr[t] << endl;


			n_phot = n_phot*coll_corr;

			x[i] = rate*coll_corr*1000;
			x_sum[i] = x[i];
			y[i] = noc/n_phot;
			// x_err[i] = sqrt(rate)*coll_corr;
			y_err[i] = sqrt(1/noc+1/n_phot);

			// if (y[i] == 0) y_err[i] = 0.2;
			// if (t == 3) cout << y[i] << endl;

			ef_graph[t]->SetPoint(i,x[i],y[i]);
			ef_graph[t]->SetPointError(i,0,y_err[i]);

			if (t != 0)
				ef_sum[t-1]->SetPoint(i,x[i],y0[i] + y[i]);
			else
			{
				y0[i] = y[i];
				y_sum[i] = y[i];
			}

			if ((t == 1) || (t == 2))
				y_sum[i] += y[i];
			else if (t == 3)
			{
				y_sum[i] += y[i];
				ef_sum[t]->SetPoint(i,x_sum[i], y_sum[i]);
			}
		}
		// Generating graph name
		char loc_name[16];
		sprintf(loc_name,"th = %1.1f",thr[t]);

		if (t != 3)
			sprintf(all_name + strlen(all_name)," %1.1f,",thr[t]);
		else
			sprintf(all_name + strlen(all_name)," %1.1f",thr[t]);

		if ((t == 1) || (t == 2))
			sprintf(sum_name + strlen(sum_name)," %1.1f,",thr[t]);
		if (t == 3)
			sprintf(sum_name + strlen(sum_name)," %1.1f",thr[t]);

		mg->Add(ef_graph[t]);
		ef_graph[t]->SetTitle(loc_name);
	}

	// can be implemented with a loop for 4 thresholds
	// ef_sum1 = new TGraphErrors();
	// ef_sum2 = new TGraphErrors();
	// for (int i = 0; i < nentries; i++)
	// {
	// 	cout << ef_graph[0]->GetPoint(i,i) << endl;
	// 	ef_sim1 = SetPoint(i,ef_graph[0]->GetPoint(i,i)+ef_graph[1]->GetPoint(i,i));
	// 	// ef_sim2 = SetPoint(i,ef_graph[0]->GetPoint(i));
	// }


	ef_graph[0]->Draw("A*");
	ef_graph[1]->Draw("SAME *");
	ef_graph[2]->Draw("SAME *");
	ef_graph[3]->Draw("SAME *");

	ef_graph[0]->GetXaxis()->SetTitle("Rate");
	ef_graph[0]->GetYaxis()->SetTitle("Efficiency");

	ef_graph[0]->SetMarkerColor(kRed+1);
	ef_graph[1]->SetMarkerColor(kBlue+1);
	ef_graph[2]->SetMarkerColor(kGreen+1);
	ef_graph[3]->SetMarkerColor(kMagenta+1);

	c1->SetLogx();
	c1->SetGridx();
	c1->SetGridy();
	ef_graph[0]->SetTitle(all_name);
	ef_graph[0]->GetYaxis()->SetRangeUser(0,1.075);

	// c1->Modified();
	// c1->Clear();
	// mg->Draw("a* fb l3d");

	TCanvas *c2 = new TCanvas("Sum of Effs");
	ef_sum[0]->Draw("A*");
	ef_graph[0]->Draw("SAME *");
	ef_sum[1]->Draw("SAME *");
	ef_sum[2]->Draw("SAME *");
	ef_sum[3]->Draw("SAME *");

	ef_sum[0]->GetXaxis()->SetTitle("Rate");

	ef_sum[0]->SetMarkerColor(kBlue+1);
	ef_sum[1]->SetMarkerColor(kGreen+1);
	ef_sum[2]->SetMarkerColor(kMagenta-3);
	// ef_sum[3]->SetMarkerColor(kCyan+1);

 	c2->SetLogx();
	c2->SetGridx();
	c2->SetGridy();
	ef_sum[0]->SetTitle(sum_name);
	ef_sum[0]->GetYaxis()->SetRangeUser(0,1.1);
	// ef_sum[0]->GetYaxis()->SetRangeUser(0,3.1);

	// c2->Modified();
}

///

void get_dev3()
{
	int E_0 = 14;
	gr = new TGraphErrors();
	double x[500] = {0.,};
	double y[500] = {0.,};
	double xy[50][50] ={0.,}{0.,};

	h2 = new TH2F("h2","thresholds", 32, 0.5, 3.6, 32, 0.5, 3.6);

	// making the root filename fot given energy
	char filename[32];
	snprintf(filename, sizeof(filename), "rate_scan_%d_cad.root", E_0);


	// opening the rate_scan tree
	TFile *f = new TFile(filename);
	TTree *Tr = (TTree*)f->Get("Tr");
	int rate, n_phot;
	TH1F *s_hist;
	Tr->SetBranchAddress("rate_c",&rate);
	Tr->SetBranchAddress("n_phot_c",&n_phot);
	Tr->SetBranchAddress("s_hist_c",&s_hist);

	int ch_sh;
	if (E_0 == 14)
		ch_sh = 0;
	else if (E_0 == 18)
		ch_sh = 1;
	else if (E_0 == 22)
		ch_sh = 2;


	double charge_sh[3] = {24.08,23.55,23.71};
	// Collected charge correction
	double coll_corr = 50./(50.+charge_sh[ch_sh]);

	// TMultiGraph *mg = new TMultiGraph();

	double thr1_low_edge = 0.5;
	double thr1_high_edge = 3.6;
	double thr2_low_edge = 0.5;
	double thr2_high_edge = 3.6;


	double thr0 = 0.5;
	Int_t nentries = (Int_t)Tr->GetEntries();
	for (int entry_num = 0; entry_num <= 107/*nentries*/; entry_num++)
	{
		// finding efficiency for thr0(0.5) to not repeat it on every threshold pair call
		Tr->GetEntry(entry_num);

		// getting number of counts at given threshold
		double thr0 = 0.5;
		int noc_bin = s_hist->FindBin((int)(E_0*1000))*thr0;
		double noc = s_hist->GetBinContent(noc_bin);

		n_phot = n_phot*coll_corr;

		double y0 = noc/n_phot;

		x[entry_num] = rate*coll_corr;

		for (double thr1 = thr1_low_edge; thr1 < thr1_high_edge; thr1+=0.1)
			for (double thr2 = thr2_low_edge; thr2 < thr2_high_edge; thr2+=0.1)
			{
				// buffer for deviation calculation
				double thr_buffer = 0;
				// generating array with given thresholds
				double thr[2] = {thr1,thr2};

				int i1 = thr1*10-4;
				int i2 = thr2*10-4;

				for (int t = 0; t < 2; t++)
				{
					// getting number of counts at given threshold
					int noc_bin = s_hist->FindBin((int)(E_0*1000))*thr[t];
					double noc = s_hist->GetBinContent(noc_bin);

					thr_buffer += noc/n_phot;
				}
				// adding threshold counts at 0.5 energy
				thr_buffer += y0;

				// cout << i1 << "\t" << i2 << "\t" << xy[i1][i2] << endl;
				// cout << i1 << "\t" << i2 << "\t" << thr_buffer << "\t" << (1 - thr_buffer)*(1 - thr_buffer) << "\n\n";
				// // Calculating and summing deviation
				xy[i1][i2] += (1-thr_buffer)*(1-thr_buffer);
				// if ((thr_buffer > 1)&&(thr_buffer < 0.9))
				// 	xy[i1][i2] = 101;
			}

	}

	// filling the 2D histogram
	for (double thr1 = thr1_low_edge; thr1 < thr1_high_edge; thr1+=0.1)
		for (double thr2 = thr2_low_edge; thr2 < thr2_high_edge; thr2+=0.1)
		{
			int i1 = thr1*10-4;
			int i2 = thr2*10-4;
			h2->SetBinContent(h2->GetXaxis()->FindBin(thr1),h2->GetYaxis()->FindBin(thr2), xy[i1][i2]);
		}

	// ef_sum[0]->Draw("A*");
	// ef_graph[0]->Draw("SAME *");
	// ef_sum[1]->Draw("SAME *");
	// ef_sum[2]->Draw("SAME *");
	// ef_sum[3]->Draw("SAME *");

	// ef_sum[0]->GetXaxis()->SetTitle("Rate");

	// ef_sum[0]->SetMarkerColor(kBlue+1);
	// ef_sum[1]->SetMarkerColor(kGreen+1);
	// ef_sum[2]->SetMarkerColor(kMagenta-3);
	// // ef_sum[3]->SetMarkerColor(kCyan+1);

 // 	c1->SetLogx();
	// c1->SetGridx();
	// c1->SetGridy();
	// ef_sum[0]->SetTitle(sum_name);
	// ef_sum[0]->GetYaxis()->SetRangeUser(0,1.1);

	// TCanvas *c2 = new TCanvas("Scan");
	h2->Draw("CONT1");
	h2->GetXaxis()->SetRangeUser(5,3.5);
	h2->GetYaxis()->SetRangeUser(5,3.5);

	h2->GetXaxis()->SetTitle("thr1");
	h2->GetYaxis()->SetTitle("thr2");
	h2->GetZaxis()->SetTitle("Efficiency");
	h2->GetXaxis()->SetTitleOffset(1.5);
	h2->GetYaxis()->SetTitleOffset(1.5);
	h2->GetZaxis()->SetTitleOffset(1.5);
}

///

void get_dev4(int r = 2700)
{
	int E_0 = 14;
	h2 = new TH2F("h2","thresholds", 31, 0.5, 3.5, 31, 0.5, 3.5);
	// h2 = new TH2F("h2","thresholds");


	// making the root filename fot given energy
	char filename[32];
	snprintf(filename, sizeof(filename), "rate_scan_%d.root", E_0);


	// opening the rate_scan tree
	TFile *f = new TFile(filename);
	TTree *Tr = (TTree*)f->Get("Tr");
	int rate, n_phot;
	TH1F *s_hist;
	Tr->SetBranchAddress("rate",&rate);
	Tr->SetBranchAddress("n_phot",&n_phot);
	Tr->SetBranchAddress("s_hist",&s_hist);

	int ch_sh;
	if (E_0 == 14)
		ch_sh = 0;
	else if (E_0 == 18)
		ch_sh = 1;
	else if (E_0 == 22)
		ch_sh = 2;


	double charge_sh[3] = {24.08,23.55,23.71};
	// Collected charge correction
	double coll_corr = 50./(50.+charge_sh[ch_sh]);

	// TMultiGraph *mg = new TMultiGraph();

	// Making graph array
	char graphname1[16];
	char graphname2[16];
	TGraphErrors *ef_graph[10];
	TGraphErrors *ef_sum[10];
	for (int i = 0; i < 4; i++)
	{
		sprintf(graphname1,"ef_graph_%d",i);
		ef_graph[i] = new TGraphErrors();

		sprintf(graphname2,"ef_sum_%d",i);
		ef_sum[i] = new TGraphErrors();
	}


	char all_name[16];
	sprintf(all_name,"th =");
	char sum_name[32];
	sprintf(sum_name,"th = 0.5, 0.5 +");
	Int_t nentries = (Int_t)Tr->GetEntries();

	double thr0 = 0.5;
	double thr1 = 1.5;
	int scan_thr = 0;
	// for (double thr1 = 0.5; thr1 <= 3.5; thr1+=0.5)
		for (double thr2 = 0.5; thr2 <= 3.5; thr2+=0.1)
			for (double thr3 = 0.5; thr3 <= 3.5; thr3+=0.1)
			{
				// cout << thr2*10-4 << endl;
		
		// getting number of counts at given threshold
		int noc_bin = s_hist->FindBin((int)(E_0*1000))*thr[t];
		double noc = s_hist->GetBinContent(noc_bin);
				// generating array with given thresholds
				double thr[4] = {0.,};
				double thr[4] = {thr0,thr1,thr2,thr3};

				double x_sum[500] = {0,};
				double y_sum[500] = {0,};

				// double y_dev = 0;
				for (int t = 0; t < 4; t++)
				{
					double x[500] = {0,};
					double y[500] = {0,};
					double x_err[500] = {0,};
					double y_err[500] = {0,};


					// if (t > 1)
					// cout << thr[t]/20. << endl;
					for (int i = 0; i < nentries; i++)
					{

						Tr->GetEntry(i);
						
						// getting number of counts at given threshold
						int noc_bin = s_hist->FindBin((int)(E_0*1000))*thr[t];
						double noc = s_hist->GetBinContent(noc_bin);

						// cout << "Number of count for " << rate << " is " << noc << " for threshold of " << thr[t]/20. << endl;


						n_phot = n_phot*coll_corr;

						x[i] = rate*coll_corr*1000;
						// if (i != 0) cout << x[i] - x[i-1] << endl;
						x_sum[i] = x[i];
						y[i] = noc/n_phot;
						x_err[i] = sqrt(rate)*coll_corr;
						y_err[i] = sqrt(1/noc+1/n_phot);

						if ((t == 0) && (x[i] < r*1000))
							scan_thr = i;


						// ef_graph[t]->SetPoint(i,x[i],y[i]);
						// ef_graph[t]->SetPointError(i,x_err[i],y_err[i]);

						y_sum[i] += y[i];

					}


					// // Generating graph name
					// char loc_name[16];
					// sprintf(loc_name,"th = %1.1f",thr[t]);

					// if (t != 3)
					// 	sprintf(all_name + strlen(all_name)," %1.1f,",thr[t]);
					// else
					// 	sprintf(all_name + strlen(all_name)," %1.1f",thr[t]);

					// if ((t == 1) || (t == 2))
					// 	sprintf(sum_name + strlen(sum_name)," %1.1f,",thr[t]);
					// if (t == 3)
					// 	sprintf(sum_name + strlen(sum_name)," %1.1f",thr[t]);

					// // mg->Add(ef_graph[t]);
					// ef_graph[t]->SetTitle(loc_name);
				}
				// cout << y_sum[87] << endl;

				// int asd = h2->FindBin(qwe1,qwe2);
				h2->SetBinContent(h2->GetXaxis()->FindBin(thr2),h2->GetYaxis()->FindBin(thr3), y_sum[scan_thr]);
				cout << thr2 << " " << h2->GetXaxis()->FindBin(thr2) << "\t" << thr3 << " " << h2->GetYaxis()->FindBin(thr3) << "\t" << y_sum[scan_thr] << endl;
				// cout << thr2/2-4 << " " << thr2/20. << "\t" << thr3/2-4 << " " << thr3/20. << "\t" << y_sum[scan_thr] << endl;

				// cout << thr2/20. << "\t" << thr3/20. << "\t" << y_sum[scan_thr] << endl;
				// cout << "!!!\n";

				// h2->SetBinContent((int)(thr2*10-5), (int)(thr3*10-5),y_sum[scan_thr]);
				// h2->SetBinContent(thr2*10-5, thr3*10-5,y_sum[scan_thr]);
				// cout << thr2 << "\t" << (thr2*10-4) << "\t" << thr3 << "\t" << (thr3*10-4) << "\t" << y_sum[scan_thr] << "\t" << h2->GetBinContent(thr2*10, thr3*10) << endl;

				// for (int i = 0; i < nentries; i++)
					// y_dev += (y_sum[i] - 1) * (y_sum[i] - 1);

				// y_dev = sqrt(y_dev/nentries);

				// h2->SetBinContent(h2->GetXaxis()->FindBin(thr2)-1, h2->GetYaxis()->FindBin(thr3)-1,y_dev);
				
				// cout << h2->GetXaxis()->FindBin(thr2)-1 << "\t" << h2->GetYaxis()->FindBin(thr3)-1 << "\t" << y_dev << endl;
				// h2->SetBinContent(h2->GetXaxis()->FindBin(thr2)-1, h2->GetYaxis()->FindBin(thr3)-1,y_sum[scan_thr]);

				// cout << "x " << h2->GetXaxis()->FindBin(thr2)-1 << "\ny " << h2->GetYaxis()->FindBin(thr3)-1 << "\nval: " << y_sum[scan_thr] << endl;
			}
	

	// ef_sum[0]->Draw("A*");
	// ef_graph[0]->Draw("SAME *");
	// ef_sum[1]->Draw("SAME *");
	// ef_sum[2]->Draw("SAME *");
	// ef_sum[3]->Draw("SAME *");

	// ef_sum[0]->GetXaxis()->SetTitle("Rate");

	// ef_sum[0]->SetMarkerColor(kBlue+1);
	// ef_sum[1]->SetMarkerColor(kGreen+1);
	// ef_sum[2]->SetMarkerColor(kMagenta-3);
	// // ef_sum[3]->SetMarkerColor(kCyan+1);

 // 	c1->SetLogx();
	// c1->SetGridx();
	// c1->SetGridy();
	// ef_sum[0]->SetTitle(sum_name);
	// ef_sum[0]->GetYaxis()->SetRangeUser(0,1.1);

	TCanvas *c2 = new TCanvas("Scan");
	h2->Draw("CONT1");
	h2->GetXaxis()->SetRangeUser(5,35);
	h2->GetYaxis()->SetRangeUser(5,35);

	h2->GetXaxis()->SetTitle("thr2");
	h2->GetYaxis()->SetTitle("thr3");
	h2->GetZaxis()->SetTitle("Efficiency");
	h2->GetXaxis()->SetTitleOffset(1.5);
	h2->GetYaxis()->SetTitleOffset(1.5);
	h2->GetZaxis()->SetTitleOffset(1.5);
}

void s_scan(int r1 = 0, int r2 = 20, int r3 = 50)
{
	int E_0 = 14;
	double r[3] = {r1,r2,r3};
	double slope[200000] = {0,};
	slope[0] = 1;					//to avoid having a peak on 0 point

	// making the root filename fot given energy
	char filename[32];
	snprintf(filename, sizeof(filename), "rate_scan_%d.root", E_0);


	// opening the rate_scan tree
	TFile *f = new TFile(filename);
	TTree *Tr = (TTree*)f->Get("Tr");
	int rate, n_phot;
	TH1F *s_hist;
	Tr->SetBranchAddress("rate",&rate);
	Tr->SetBranchAddress("n_phot",&n_phot);
	Tr->SetBranchAddress("s_hist",&s_hist);

	TMultiGraph *mg = new TMultiGraph();

	// THStack *hs = new THStack("hs","Stacked histograms");

	// Tr->GetEntry(0);
	// cout << s_hist->GetNbinsX() << "!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
	// // int s_size = s_hist->GetNbinsX();
	// int s_size = s_hist->GetXaxis()->GetXmax();
	// // Making graph array
	// char histname[16];
	// TH1F *s_hists[10];
	// for (int i = 0; i < 3; ++i)
	// {
	// 	sprintf(histname,"s_hist_%d",i);
	// 	s_hists[i] = new TH1F(histname, "qwe", s_size, 0, s_size);
	// }

	// for (int h = 0; h < 3; h++)
	// {
	// 	Tr->GetEntry(h);

	// 	s_hists[h] = s_hist;
	// 	hs->Add(s_hist);

	// }

	Tr->GetEntry(r1);
	cout << "Plotted graph for the rate of " << rate << endl;
	s_hist->DrawCopy();
	Tr->GetEntry(r2);
	cout << "Plotted graph for the rate of " << rate << endl;
	s_hist->DrawCopy("SAME");
	Tr->GetEntry(r3);
	cout << "Plotted graph for the rate of " << rate << endl;
	s_hist->DrawCopy("SAME");


	// Making graph array
	char graphname[16];
	TGraphErrors *slope_graph[10];
	for (int i = 0; i < 3; i++)
	{
		sprintf(graphname,"slope_graph_%d",i);
		slope_graph[i] = new TGraphErrors();
	}



	double values[4] = {0,};
	for (int t = 0; t < 3; t++)
	{
		Tr->GetEntry(r[t]);

		int s_size = s_hist->GetNbinsX();

		// Slope calculation
		for (int i = AV_POINTS/2; i <= s_size-1 - AV_POINTS/2; i++)
		{
			// Emptying the buffers
			for (int iii = 0; iii <= 3; iii++)
				values[iii] = 0;

			for (int ii = i - AV_POINTS/2; ii <= i + AV_POINTS/2; ii++)
			{
				values[0] += ii;								//Calculating sum(x)
				values[1] += s_hist->GetBinContent(ii);			//Calculating sum(y)
				values[2] += ii * s_hist->GetBinContent(ii);	//Calculating sum(x*y)
				values[3] += ii * ii;							//Calculating sum(x^2)

			}
			slope[i] = -(values[2] - values[0] * values[1] / AV_POINTS)/(values[3] - values[0] * values[0]/AV_POINTS);
			slope_graph[t]->SetPoint(i,i,slope[i]);
			mg->Add(slope_graph[t]);
		}

	}

	slope_graph[0]->Draw("A*");
	// mg->Draw("a* fb 3d");
	slope_graph[0]->GetYaxis()->SetRangeUser(-4,16);


	// for (int h = 0; h < 3; h++)
	// {
	// 	Tr->GetEntry(r[h]);

	// 	if (h == 0)
	// 		s_hist->Draw("A*");
	// 	else
	// 		s_hist->Draw("SAME *");

	// }

	// s_hists[0]->DrawCopy();
	// s_hists[1]->DrawCopy("SAME");
	// s_hists[2]->DrawCopy("SAME");

	// ef_graph[0]->GetXaxis()->SetTitle("Rate");
	// ef_graph[0]->GetYaxis()->SetTitle("Efficiency");

	// s_hists[0]->SetMarkerColor(kRed+1);
	// s_hists[1]->SetMarkerColor(kBlue+1);
	// s_hists[2]->SetMarkerColor(kGreen+1);

	// ef_graph[0]->GetYaxis()->SetRangeUser(0,1.075);

	// mg->Draw("a* fb l3d");
}

void comp(int r = 10)
{
	gr1 = new TGraph();
	gr2 = new TGraph();
	gr3 = new TGraph();
	double x[500] = {0.,};
	double y1[500] = {0.,};
	double y2[500] = {0.,};
	double y3[500] = {0.,};

	// reading exponent data
	TFile *f = new TFile("rate_scan_14.root");
	TTree *Tr = (TTree*)f->Get("Tr");
	int rate, n_phot;
	TH1F *s_hist;
	Tr->SetBranchAddress("rate",&rate);
	Tr->SetBranchAddress("n_phot",&n_phot);
	Tr->SetBranchAddress("s_hist",&s_hist);

	// reading cadence data
	TFile *f_c = new TFile("rate_scan_14_cad.root");
	TTree *Tr_c = (TTree*)f_c->Get("Tr");
	int rate_c, n_phot_c;
	TH1F *s_hist_c;
	Tr_c->SetBranchAddress("rate_c",&rate_c);
	Tr_c->SetBranchAddress("n_phot_c",&n_phot_c);
	Tr_c->SetBranchAddress("s_hist_c",&s_hist_c);

	// reading experimental data
	TFile *g = new TFile("14_exp.root");
	TTree *Tr_exp = (TTree*)g->Get("Tr_exp");
	// Tr_exp->Scan();
	TH1F *s_hist_exp;
	Tr_exp->SetBranchAddress("s_hist_exp",&s_hist_exp);


	Int_t nentries = (Int_t)Tr->GetEntries();
	for (int i = 0; i <= nentries; i++)
	{
		Tr_c->GetEntry(i);
		Tr->GetEntry(i);

		x[i] = rate;
		if (s_hist->GetBinContent(s_hist->FindBin(7000)) != 0)
			y1[i] = s_hist_c->GetBinContent(s_hist_c->FindBin(7000))/s_hist->GetBinContent(s_hist->FindBin(7000));
		if (s_hist->GetBinContent(s_hist->FindBin(10500)) != 0)
			y2[i] = s_hist_c->GetBinContent(s_hist_c->FindBin(10500))/s_hist->GetBinContent(s_hist->FindBin(10500));
		if (s_hist->GetBinContent(s_hist->FindBin(17500)) != 0)
			y3[i] = s_hist_c->GetBinContent(s_hist_c->FindBin(17500))/s_hist->GetBinContent(s_hist->FindBin(17500));

		gr1->SetPoint(i,x[i],y1[i]);
		gr2->SetPoint(i,x[i],y2[i]);
		gr3->SetPoint(i,x[i],y3[i]);
		
	}

	// TCanvas *c2 = new TCanvas("Ratio of cadence and exponent");

	gr1->Draw("A*");
	gr2->Draw("SAME *");
	gr3->Draw("SAME *");
	

	gr1->SetTitle("Ratio of cadence and exponent");
	gr1->GetXaxis()->SetTitle("Rate [kHz]");
	gr1->GetYaxis()->SetTitle("cadence/exponent");
	gr1->GetXaxis()->SetTitleOffset(1.2);
	gr1->SetMarkerColor(kBlue+1);
	gr2->SetMarkerColor(kGreen+1);
	gr3->SetMarkerColor(kRed+1);
	c1->SetGridx();
	c1->SetGridy();


	// plotting cadence and exponent pulses' s-curves
	TCanvas *c2 = new TCanvas("Ratio of cadence and exponent at given thresholds");

	//Cadence and exponent pulses' s-curves
	int entry_num;
	Int_t nentries = (Int_t)Tr->GetEntries();
	// Finding entry number
	for (int i = 0; i < nentries; i++)
	{
		Tr->GetEntry(i);
		if (rate == r)
			entry_num = i;

	}

	Tr_c->GetEntry(entry_num);
	Tr->GetEntry(entry_num);
	cout << "Data at rate " << r << "\nEntry number is " << entry_num << endl;
	s_hist_c->Draw();
	s_hist->Draw("SAME");

	s_hist_c->SetTitle("Cadence and exponent pulses' s-curves");
	s_hist_c->GetXaxis()->SetTitle("threshold");
	s_hist_c->GetYaxis()->SetTitle("counts");
	s_hist_c->GetYaxis()->SetTitleOffset(1.4);
	s_hist_c->SetLineColor(kRed+1);
	s_hist->SetLineColor(kBlue+1);
	c2->SetGridx();
	c2->SetGridy();

	// cout << s_hist->GetBinContent(s_hist->FindBin(14000)) << "\t" << s_hist_c->GetBinContent(s_hist_c->FindBin(14000)) << endl;

	// cout << rate << "\t" << rate_c << endl;
	// cout << n_phot << "\t" << n_phot_c << endl;

}

read_tree()
{
	TFile *f = new TFile("rate_scan_14_4000.root");
	TTree *Tr = (TTree*)f->Get("Tr");
	int rate, n_phot;
	TH1F *s_hist;
	Tr->SetBranchAddress("rate",&rate);
	Tr->SetBranchAddress("n_phot",&n_phot);
	Tr->SetBranchAddress("s_hist",&s_hist);
	return 0;
}
