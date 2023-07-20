////////////////////////////////////////////////////////////////////////////
//
//	Program do analizy plikow zawierajacych cale sekwencje pomiarow przeprowadzonych przy uzyciu ukladu
//  Caen DT5730 oraz programu Compass. Skrypt ten w glownej mierze sluzy tylko do pomiarow z 03.2023
//  Jednakze rdzen sortowania eventow i dopasowywania funckji jest ten dla wszystkich pomiarow
//  w ktorych uzywany byl uklad firmy Caen
//
//	Autor: Przemyslaw Sekowski
//
//	email: przemyslaw.sekowski@fuw.edu.pl
//
///////////////////////////////////////////////////////////////////////////
#include "tspectrum_coincidence2.hpp"
using namespace std;

//Glowna czesc skryptu
void tspectrum_coincidence_const_count()
{

///////////////////////////////////////////         Do zmiany            //////////////////////////////////////////////////////////////////////
	Double_t ile_sigma = 5;
	Int_t liczba_pomiarow = 50;
	Int_t min_count_per_spectrum = 1000;
	Char_t nazwa_wejsciowego[200] = "SDataR_sio2_1.root";
	Char_t nazwa_pliku_wyjsciowego[200] = "Analysis_coincidence_const_count_SDataR_sio2_1.root";
	Int_t para_detektorow = -1; // -1 - Wszystkie detektory; 0 - detektory 0 i 1; 1 - detektory 2 i 3; 2 - detektory 4 i 5;
////////////////////////////////////////////////////////////////////////////////////////////////////////////////




	czas_koincydencji = { 25000,
		20000,
		30000
	};


	Int_t liczba_par_det, start_det, stop_det, liczba_det;

	if (para_detektorow == -1)
	{
		liczba_par_det = 3;
		liczba_det = 6;
		start_det = 0;
		stop_det = 6;
	}
	else
	{
		liczba_par_det = 1;
		liczba_det = 2;
		start_det = 2 * para_detektorow;
		stop_det = start_det + 2;
	}

	Int_t numberofbins = 950;	//liczba binow w widme energergetycznym; zmniejszajac jej wartosc otrzymujemy lepsza statystykew binach ale gorsza rozdzielczosc (mniej pewna wartosc dopasowanej centroidy)
	Float_t minimum = 0, maksimum = 17000;	//zakres widma energetycznego wyrazone w kanalach
	Double_t zliczenia_amplituda = maksimum / numberofbins;
	TH1F *h[liczba_det][liczba_pomiarow];
	TH1F *total_h[liczba_det];
	TH2F *h_2d[liczba_par_det][liczba_pomiarow];
	TH2F *total_h_2d[liczba_par_det];
	TH1F *h_delta_time[liczba_par_det];
	TF1 *dopasowanie[liczba_det][liczba_pomiarow];	//tworzona jest macierz funkcji dla kazdego widma z kazdej pozycji i detektora
	//plik zrodlowy z danymi
	auto *f_1 = new TFile(nazwa_wejsciowego);	//tworzone sa tu plik wyjsciowy i wszystkie histogramy, a takze do zmiennych energia, czas i channel przypisywane sa wartosci eventow
	auto *t_1 = (TTree*) f_1->Get("Data_R");
	t_1->SetBranchAddress("Energy", &energia);
	t_1->SetBranchAddress("Timestamp", &czas);
	t_1->SetBranchAddress("Channel", &channel);
	auto *f_output = new TFile(nazwa_pliku_wyjsciowego, "RECREATE");
	vector<vector < Double_t>> zliczenia(liczba_det, vector<Double_t> (liczba_pomiarow));
	vector<vector < Double_t>> blad_zliczenia(liczba_det, vector<Double_t> (liczba_pomiarow));
	vector<vector < Double_t>> blad_zliczenia_fixed(liczba_det, vector<Double_t> (liczba_pomiarow));
	vector<vector< ULong64_t>> wektor_timestamp(liczba_par_det, vector <ULong64_t> (czasy));
	vector<vector< Int_t>> wektor_entry(liczba_par_det, vector<Int_t> (liczba_pomiarow));
	vector<Int_t> counts (liczba_det);
	vector<UShort_t> centroidy (liczba_det);
	vector<UShort_t> sigmy (liczba_det);
	std::vector<Double_t> wektor_czasu(300);


	zakres_energii[0][0] = 2600, zakres_energii[0][1] = 4500, zakres_energii[0][2] = 2400, zakres_energii[0][3] = 1700, zakres_energii[0][4] = 1400, zakres_energii[0][5] = 3700;
	zakres_energii[1][0] = 3400, zakres_energii[1][1] = 5600, zakres_energii[1][2] = 3050, zakres_energii[1][3] = 2300, zakres_energii[1][4] = 2000, zakres_energii[1][5] = 4400;
	cout << "Tworzenie histogramow" << endl;
	char name[20];
	char title[100];

	for (Int_t i = 0; i < liczba_det; i++)
	{
		for (Int_t m = 0; m < liczba_pomiarow; m++)
		{
			sprintf(name, "spek_%d_det_%d", m + 1, i);
			sprintf(title, "Spektrum pozycji h%d det %d", m + 1, i);
			h[i][m] = new TH1F(name, title, numberofbins, minimum, maksimum);
			sprintf(name, "spek_2D_%d_det_%d", m + 1, i);
			sprintf(title, "Spektrum 2D pozycji h%d det %d", m + 1, i);
			if (i % 2 == 0) h_2d[i / 2][m] = new TH2F(name, title, 100, 1000, 6000, 150, 1000, 6000);
		}

		cout << "total det: " << i << endl;
		sprintf(name, "spek_total_det_%i", i);
		sprintf(title, "Spektrum calkowite det %i", i);
		total_h[i] = new TH1F(name, title, 4000, minimum, maksimum);
		cout << "total det: " << i << endl;
		sprintf(name, "spek_delta_det_%i", i);
		sprintf(title, "Spektrum delta time det %i", i);
		h_delta_time[i] = new TH1F(name, title, 5e3, 0, 1e5);
		cout << "Utworzono histogram: name: " << name << " title: " << title << " obiekt " << &total_h[i] << endl;
	}

	cout << "Przygotowywanie widm czastkowych 2d" << endl;
	for (Int_t i = 0; i < liczba_par_det; i++)
	{
		sprintf(name, "spek_total_2d_det_%d", i);
		sprintf(title, "Spektrum calkowite 2D det %d", i);
		total_h_2d[i] = new TH2F(name, title, 400, minimum, maksimum, 400, minimum, maksimum);
	}

	cout << "Zakonczono tworzenie histogramow" << endl;

	nentries = (Int_t) t_1->GetEntries();

	cout << "Wypelnianie histogramow" << endl;
	for (Int_t i = 0; i < nentries; i++)
	{
		t_1->GetEntry(i);
		if (w_zakresie_511kev(energia, channel)) continue;
		if (channel == 7) h_step_time->Fill(czas / 1e12);	// widmo krokow silnika
		if (channel > 5) continue;
		if (channel < start_det || channel >= stop_det) continue;
		channel %= liczba_det;
		total_h[channel]->Fill(energia);
	}

	for (Int_t det = 0; det < liczba_det; det++)
	{
		pre_dopasowanie[det] = new TF1("dopasowanie", gausswithlinearbkg, zakres_energii[0][start_det + det], zakres_energii[1][start_det + det], 5);
		pre_dopasowanie[det]->SetParameters(20000, (zakres_energii[0][start_det + det] + zakres_energii[1][start_det + det]) / 2, 25, -0.00001, 10);
		pre_dopasowanie[det]->SetParNames("Amplituda", "Srednia", "Sigma", "A", "B");
		pre_dopasowanie[det]->SetParLimits(0, 20 *zliczenia_amplituda, 200000 *zliczenia_amplituda);
		pre_dopasowanie[det]->SetParLimits(1, zakres_energii[0][start_det + det], zakres_energii[1][start_det + det]);
		pre_dopasowanie[det]->SetParLimits(2, 25, 100);
		pre_dopasowanie[det]->SetParLimits(3, -0.1, 0.1);
		pre_dopasowanie[det]->SetParLimits(4, -10000, 10000);
		total_h[det]->Fit("dopasowanie", "LMQR0", "", zakres_energii[0][start_det + det], zakres_energii[1][start_det + det]);
		// total_h[det] -> Draw();
	}	
	
	for (Int_t det = 0; det < liczba_det; det++)
	{
		centroidy[det] = pre_dopasowanie[det]->GetParameter(1);
		sigmy[det] = pre_dopasowanie[det]->GetParameter(2);
		counts[det] = 0;
	}

	cout << "Tworzenie wektorow koincydencji" << endl;
	for (Int_t i = 0; i < nentries; i++)
	{
		t_1->GetEntry(i);
		if (w_zakresie_sigma(energia, channel, centroidy, sigmy, ile_sigma)) continue;fa
		if (channel > 5) continue;
		if (channel % 2 == 0) continue;
		if (channel < start_det || channel >= stop_det) continue;
		channel %= liczba_det;
		wektor_timestamp[(channel - 1) / 2].push_back(czas);	// wypelniany jest wektor z czasem do koincydencji
		wektor_entry[(channel - 1) / 2].push_back(i);	// wypelniane jest widmo energetyczne w zaleznosci od kanalu i pozycji
		// if (i%10000==0) cout<<"wrzucono do "<<channel/2<<" oraz pozycji "<<pozycja<<" wartosci "<<czas<<" oraz "<<i<<endl;
	}

	cout << "stworzyla sie macierz koincydencji" << endl;

	cout << "Wypelnianie widm kazdego pomiaru" << endl;
	Int_t pozycja = 0;
	for (Int_t i = 0; i < nentries; i++)
	{
		t_1->GetEntry(i);
		if (w_zakresie_sigma(energia, channel, centroidy, sigmy, ile_sigma)) continue;
		if (channel > 5) continue;
		if (channel % 2 == 1) continue;	// jesli event jest zebrany na kanale innym niz do ktorych byly podlaczone detektory, kod przechodzi do kolejnego eventu
		if (channel < start_det || channel >= stop_det) continue;
		channel %= liczba_det;
		auto closest_higher_time_index = closest(wektor_timestamp[(channel) / 2], czas);
		auto closest_lower_time_index = (closest_higher_time_index - 1) *((closest_higher_time_index - 1) > 0);
		if (closest_higher_time_index > wektor_timestamp[(channel) / 2].size() - 1 ||
			closest_lower_time_index > wektor_timestamp[(channel) / 2].size() - 1 ||
			closest_higher_time_index > wektor_entry[(channel) / 2].size() - 1 ||
			closest_lower_time_index > wektor_entry[(channel) / 2].size() - 1) continue;
		closest_time_high = wektor_timestamp[(channel) / 2][closest_higher_time_index];
		closest_time_low = wektor_timestamp[(channel) / 2][closest_lower_time_index];
		delta_time_high = max(closest_time_high, czas) - min(closest_time_high, czas);
		delta_time_low = max(closest_time_low, czas) - min(closest_time_low, czas);
		auto closest_time_index = closest_higher_time_index;
		delta_time = delta_time_high;
		if (delta_time_high > delta_time_low)
		{
			closest_time_index = closest_lower_time_index;
			delta_time = delta_time_low;
		}

		h_delta_time[(channel) / 2]->Fill(delta_time);
		// cout<<wektor_timestamp[(channel-1)/2][pozycja][closest_time_index]<<" najblizszy czas do "<<czas<<endl;
		// cout << "kanal "<<channel<<endl;
		// cout << "delta time "<<delta_time<<endl;
		if (i % 10000 == 0) cout << i << endl;
		if (delta_time < czas_koincydencji[(channel) / 2])
		{
			energia_other = energia;
			t_1->GetEntry(wektor_entry[(channel) / 2][closest_time_index]);
			if (channel < start_det || channel >= stop_det) continue;
			channel %= liczba_det;
			// if (w_zakresie_elipsy(energia_other, energia, channel, pre_dopasowanie, ile_sigma)) continue;
			h[channel - 1][pozycja]->Fill(energia_other);	// wypelniane jest widmo energetyczne w zaleznosci od kanalu i pozycji
			h[channel][pozycja]->Fill(energia);	// wypelniane jest widmo energetyczne w zaleznosci od kanalu i pozycji
			h_2d[(channel) / 2][pozycja]->Fill(energia_other, energia);
			total_h_2d[(channel) / 2]->Fill(energia_other, energia);
			n_entried_entries++;	//zwiekszana jest liczba oznaczajaca eventy ktore przeszly analize
			counts[(channel) / 2]++;
			if (counts[0] > min_count_per_spectrum)
			{
				wektor_czasu[pozycja] = czas;
				pozycja = pozycja + 1;
				for (Int_t det = 0; det < liczba_det; det++){
					counts[det] = 0;
				}
			}
		}
		else continue;
	}

	cout << "Eksport danych" << endl;
	wektor_czasu.insert(wektor_czasu.end(), czas);	// na koniec wektora czasu dorzucany jest koniec pomiaru.

	for (Int_t det = 0; det < liczba_det; det++)
	{
		for (Int_t i = 0; i < liczba_pomiarow; i++)
		{
			// h[channel-1][pozycja] -> Draw
			zliczenia[det][i] = (h[det][i]->IntegralAndError(h[det][i]->FindBin(centroidy[det]-ile_sigma*sigmy[det]), h[det][i]->FindBin(centroidy[det]+ile_sigma*sigmy[det]), blad_zliczenia[det][i]));
		}
	}

	char plik[200];
	wektor_czasu.insert(wektor_czasu.begin(), 0);
	for (Int_t det = 0; det < liczba_det; det++)
	{
		cout << "Liczba zliczen w pomiarze - detektor " << det << endl;
		sprintf(plik, "zliczenia_coinc_%d", det);
		ofstream aktywnosci(plik);
		cout << "Liczba zliczen w pomiarze - detektor " << start_det + det << endl;
		for (Int_t i = 0; i < liczba_pomiarow; i++)
		{
			aktywnosci <<wektor_czasu[i]/1e12<< "\t\t" << wektor_czasu[i+1]/1e12<< "\t\t" << zliczenia[det][i];
			aktywnosci << "\t\t" << blad_zliczenia[det][i] << endl;
		}

		aktywnosci.close();
	}

	f_output->Write();
}