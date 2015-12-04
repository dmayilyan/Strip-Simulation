#include "NewMythenMacros.C"

void get_dec_time()
{
	gr_dec = new TGraphErrors();
	gr_dec_exp = new TGraphErrors();
	char line[1024] = {0};
	ifstream f;

	f.open("./filter/dec_time_beh.dat", ios::in);
	if (!f.is_open())
		cout << "\nCouldn't open the file!\n\n";

	f.getline(line,2048);

	int i = 0;
	while(!f.eof() && strlen(line) > 0)
	{
		double r, exp_0, exp_0_err,
				sim_0, sim_0_err;

		sscanf(line, "%le\t%le\t%le\t%le\t%le", &r, &exp_0, &exp_0_err,
													&sim_0, &sim_0_err);

		gr_dec->SetPoint(i,r,sim_0);
		gr_dec->SetPointError(i,0,sim_0_err);

		// experimental with points
		gr_dec_exp->SetPoint(i,r,exp_0);
		gr_dec_exp->SetPointError(i,0,exp_0_err);

		cout << "Decay time for " << r << " is " << (sim_0 - exp_0)/exp_0_err << endl;

		f.getline(line,2048);
		i++;
	}

	TF1 *exp_lin = new TF1("exp_lin","pol0",2,10);
	exp_lin->FixParameter(0,exp_0);

	gr_dec->Draw("A*");
	// exp_lin->DrawCopy("SAME");
	// exp_lin->Draw();
	gr_dec_exp->Draw("SAME L3");



	gr_dec->GetXaxis()->SetTitle("Decay time");
	gr_dec->GetYaxis()->SetTitle("Exp fit argument");
	gr_dec->GetXaxis()->SetTitleOffset(1.3);
	gr_dec->GetYaxis()->SetTitleOffset(1.2);
	gr_dec->SetMarkerColor(kBlue+1);

	gr_dec_exp->SetMarkerColor(kGreen+1);
	gr_dec_exp->SetFillColor(kRed);

}

TGraphErrors *run(int dec_time_ns = 67, bool outside_call = false)
{
	char filename[16];
	snprintf(filename, sizeof(filename), "rate_%d.dat", dec_time_ns);

	cout << "\n////////////////////////////////////////////////////////////////////////////////\n\n";
	cout << "Opening filtering files!\n";
	cout << "\n////////////////////////////////////////////////////////////////////////////////\n\n\n";
	TH2F *hist_2d = createScan("/home/l_mayilyan/Documents/PhD/workspace/Mythen/filter/filter_scan_12400eV_%d.raw",0,40,1,1280*4);

	int last_item = 41;
	double av_x[41] = {0.,};
	double av_x_err[41] = {0.,};
	for (int channel = 1; channel <= last_item; channel++)
	{
		for (int Xrow = 20; Xrow < 120; Xrow++)
		{
			double cell_content = hist_2d->GetBinContent(Xrow,channel);
			av_x[channel-1] += cell_content;
			av_x_err[channel-1] += cell_content*cell_content;
		}

		av_x[channel-1] /= 100; // 120-20
		av_x_err[channel-1] = sqrt(av_x[channel-1]) / 100;

		// cout << hist_2d->GetBinContent(20,channel) << "\t" << hist_2d->GetBinContent(599,channel) << endl;
		// cout << hist_2d->GetBinContent(20,channel) << "\t" << hist_2d->GetBinContent(598,channel) << "\t" << av_x[channel-1] << endl;

	}

	cout << "\n\n";

	gr1 = new TGraph();
	gr2 = new TGraphErrors();
	gr3 = new TGraphErrors();

	double x[500] = {0,};
	double y[500] = {0,};
	double x_err[500] = {0,};
	double y_err[500] = {0,};

	double x_eval[500] = {0,};
	double y_eval[500] = {0,};
	// double x_eval_err[500] = {0,};
	double y_eval_err[500] = {0,};

	double x_ratio[500] = {0,};
	double y_ratio[500] = {0,};
	// double x_ratio_err[500] = {0,};
	double y_ratio_err[500] = {0,};

	int shift = 0;
	int fit_lev = 0;
	bool thresh_found = false;
	double max_ratio = 0.;
	for (int chan = 498-1; chan <= 598-1; chan++)
	{
		gr = new TGraphErrors();
		TH1F *hist_ch = getCh(hist_2d, chan);
		for (int i = 1+shift; i <= last_item+shift; i++)
		{
			// cout << (i-1)%42 << endl;
			x[(i-1)%41] = av_x[(i-1)%41];
			y[(i-1)%41] = hist_ch->GetBinContent((i-1)%41 + 1);
			x_err[(i-1)%41] = av_x_err[(i-1)%41];
			y_err[(i-1)%41] = sqrt(y[(i-1)%41]);

			// cout << (i-1)%41 << "\t" << x[(i-1)%41] << "\t" << y[(i-1)%41] << endl;

			gr->SetPoint(i-1,x[(i-1)%41],y[(i-1)%41]);
			gr->SetPointError(i-1,x_err[(i-1)%41],y_err[(i-1)%41]);

			// finding level for linear fit
			if ((y[(i-1)%41]<=50000) && (thresh_found == false))
			{
				fit_lev = (i-1)%41;
				// cout << "fit_lev " << fit_lev << endl;
				thresh_found = true;
			}

			if ((i-1)%41==40)
			{
				// checking if the call is from outside
				if (outside_call == false)
				{
					gr->Draw("A*");
					c1->SetGridx();
					c1->SetGridy();					
				}

				// cout << x[last_item] << " x[last_item]" << endl;

				TF1 *f1 = new TF1("f1","pol1",x[fit_lev],x[last_item-1]);
				f1->FixParameter(0,0);
				// f1->SetParLimits(0,-0.1,0.1);
				gr->Fit("f1", "RCQN");
	
				f1->SetRange(x[40],x[0]);
				// checking if the call is from outside
				if (outside_call == false)
				{
					f1->SetLineColor(kRed);
					f1->Draw("SAME R");
				}
					
				for (int j = 0; j < last_item; j++)
				{
					x_eval[j] = f1->Eval(x[j]);
					y_eval[j] = y[j];
					y_eval_err[j] = y_err[j];


					// cout << y[(i-1)%41] << " y[(i-1)%41]\n";

					gr2->SetPoint(j+shift,x_eval[j],y_eval[j]);
					gr2->SetPointError(j+shift,0,y_eval_err[j]);
					// cout << x_eval[j] << "\t" << y_eval[j] << endl;

					x_ratio[j] = x_eval[j];
					y_ratio[j] = y_eval[j]/x_eval[j];
					y_ratio_err[j] = y_eval_err[j]/x_eval[j];
					gr3->SetPoint(j+shift, x_ratio[j], y_ratio[j]);
					gr3->SetPointError(j+shift, 0, y_ratio_err[j]);

					if (x_ratio[j]>=max_ratio)
						max_ratio = x_ratio[j];
				}
	
			// c1->SaveAs("qwe.png");
			// int q;
			// cin >> q;
				
			}

		}
		thresh_found = false;
		shift += last_item;
	}

	cout << max_ratio << " max_ratio\n";

	// plotting rate ratio
	// checking if the call is from outside
	if (outside_call == false)
	{
		TCanvas *c3 = new TCanvas("Ratio");
		gr3->Draw("A*");
	}

	TF1 *exp_fit = new TF1("exp_fit","exp(-[0]*x)",x_ratio[0],x_ratio[40]);
	if (outside_call == false)
	{
		exp_fit->SetLineColor(kRed);
		gr3->Fit("exp_fit", "R");
	}
	else
		gr3->Fit("exp_fit", "RN");

	double exp_par0 = exp_fit->GetParameter(0);
	double exp_par0_err = exp_fit->GetParError(0);


	// TCanvas *c2 = new TCanvas("decay");
	// gr1->Draw("A*");
	// TF1 *exp_expo = new TF1("f2","pol1",x_ratio[0],x_ratio[40]);

	// return 0;


	//////////////////////////////////////////////////////////////////////////////////////////
	// analysing theoretical data

	// reading theoretical data


	double x2[500] = {0,};
	double y2[500] = {0,};
	double x2_err[500] = {0,};
	double y2_err[500] = {0,};

	double x_sim_ratio[500] = {0,};
	double y_sim_ratio[500] = {0,};
	double x_sim_ratio_err[500] = {0,};
	double y_sim_ratio_err[500] = {0,};
	
	// checking if the call is from outside
	if (outside_call == false)
		TCanvas *c3 = new TCanvas("Simulation data");

	gr4 = new TGraphErrors();
	gr5 = new TGraphErrors();

	char line[1024] = {0};
	ifstream f;

	// ifstream f("rate_" + std::to_string(dec_time_ns) + ".dat");

	f.open(filename, ios::in);
	if (!f.is_open())
		cout << "\nCouldn't open the file!\n\n";

	int dc = 0;
	int max_r = 0;

	// finding maximum rate
	// Processing the file
	f.getline(line,1024);
	while(!f.eof() && strlen(line) > 0)
	{
		double r, tc, c;
		sscanf(line, "%lf\t%lf\t%lf", &r, &tc, &c);
		// printf("a = %lf, \tb = %lf, \tc = %lf\n", r, tc, c);

		x2[dc] = r;
		if (dc != 0)
			if (x2[dc] > x2[dc-1])
				max_r = x2[dc];

		dc++;

		f.getline(line, 1024);
	}

	f.clear();
	f.seekg(0, ios::beg);

	double charge_sh[3] = {24.08,23.55,23.71};
	double coll_corr = 50./(50.+charge_sh[0]);

	dc = 0;
	f.getline(line,1024);
	while(!f.eof() && strlen(line) > 0)
	{
		double r, tc, c;
		sscanf(line, "%lf\t%lf\t%lf", &r, &tc, &c);


		// cout << r << "\t" << tc << "\t" << c << " QWE\n";
		tc = tc*coll_corr;

		x2[dc] = r*coll_corr*1000;
		y2[dc] = c/tc;
		x2_err[dc] = sqrt(r)*coll_corr;
		y2_err[dc] = sqrt(1/c+1/tc);

		gr4->SetPoint(dc,x2[dc],y2[dc]);
		gr4->SetPointError(dc,x2_err[dc],y2_err[dc]);

		dc++;

		f.getline(line, 1024);
	}

	f.close();

	int fit_lev = 0;
	for (int i = 0; i <= dc; i++)
	{
		if (x2[i]>=max_ratio)
		{
			fit_lev = i;
			cout << "fit_lev " << fit_lev << endl;
			break;
		}

	}

	// TCanvas *c3 = new TCanvas("Ratio sim");
	// checking if the call is from outside
	if (outside_call == false)
		gr4->Draw("A*");

	TF1 *sim_fit = new TF1("sim_fit","exp(-[0]*x)",x2[0],x2[fit_lev]);
	sim_fit->SetParameter(0,1.55e-7);
	if (outside_call == false)
	{
		sim_fit->SetLineColor(kRed);
		gr4->Fit("sim_fit", "R");
	}
	else
		gr4->Fit("sim_fit", "RN");

	double sim_par0 = sim_fit->GetParameter(0);
	double sim_par0_err = sim_fit->GetParError(0);
	gr4->GetXaxis()->SetRangeUser(0,x2[fit_lev]);


	printf("\n////////////////////////////////////////////////////////////////////////////////\n\n");
	printf("exp_0\texp_0_err\tsim_0\tsim_0_err\n");
	printf("%e\t%e\t%e\t%e\n", exp_par0,exp_par0_err,sim_par0,sim_par0_err);
	printf("\n////////////////////////////////////////////////////////////////////////////////\n\n");


	return gr3;

}