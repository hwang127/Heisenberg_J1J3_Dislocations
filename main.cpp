//#include"stdafx.h"


#include"search_for.h"
#include <vector>
#include<stdio.h>
#include<iostream>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <string>
#include <fstream>

using namespace std;
int str_to_int(string str)
{
	return stoi(str);
}
long long int str_to_long_long(string str)
{
	return stoll(str);
}
double str_to_d(string str)
{
	return stod(str);
}
double ourrandom()
{
	return rand() / (RAND_MAX + 1.);
}
void read_file(string infilename, int &L, double &J,double &alpha, int &num,double &Jthird, double &T, long long int &nsamples,
	int &seed, string &outfilename)
{
	string str_ret;
	bool found;

	search_for(string("L"), infilename, str_ret, found);
	if (found) { L = str_to_int(str_ret); }
	else { L = 4; }

	search_for(string("num"), infilename, str_ret, found);
	if (found) { num = str_to_int(str_ret); }
	else { num = 100; }


	//cout<<str_ret<<L<<endl;
	search_for(string("J"), infilename, str_ret, found);
	if (found) { J = str_to_d(str_ret); }
	else { J = 0.0; }

	search_for(string("alpha"), infilename, str_ret, found);
	if (found) { alpha = str_to_d(str_ret); }
	else { alpha = 0.0; }

	//cout<<"J"<<str_ret<<J<<endl;
	search_for(string("Jthird"), infilename, str_ret, found);
	if (found) { Jthird = str_to_d(str_ret); }
	else { Jthird = 0.0; }
	//cout<<"Jthird"<<str_ret<<Jthird<<endl;
	search_for(string("T"), infilename, str_ret, found);
	if (found) { T = str_to_d(str_ret); }
	else { cout << "T not found" << endl; T = 0.0; }

	search_for(string("nsamples"), infilename, str_ret, found);
	if (found) { nsamples = str_to_long_long(str_ret); }
	else { nsamples = 1000; }


	search_for(string("seed"), infilename, str_ret, found);
	if (found) { seed = str_to_int(str_ret); }
	else { seed = 0; }


	search_for(string("outfilename"), infilename, outfilename, found);
	if (!found) { outfilename = "C.txt"; }
}

int randint(int nsites)
{
	return rand() % nsites;
}

vector<vector<double>> create_config(int L)
{
	vector<vector<double>> config(L*L);

	for (int i = 0; i < L*L; i++) {
		double x1, x2;
		double new_x, new_y, new_z;
		while (true) {
			x1 = (ourrandom() - 0.5) *2.0;
			x2 = (ourrandom() - 0.5) *2.0;
			if (x1*x1 + x2*x2 < 1.0) break;
		}
		new_x = 2.0 * x1 * sqrt(1.0 - x1*x1 - x2*x2);
		new_y = 2.0 * x2 * sqrt(1.0 - x1*x1 - x2*x2);
		new_z = 1.0 - 2.0 * (1.0 - x1*x1 - x2*x2);
		config[i] = { new_x,new_y,new_z };
	} 
	return config;
}
int Right(int i,int L){return (i + 1) % L + i / L*L;}
int Left(int i, int L){return (i - 1 + L) % L + i / L*L;}
int TopLeft(int i, int L) { return (i + L) % (L*L); }
int TopRight(int i,int L){return (TopLeft(i,L) + 1) % L + TopLeft(i,L) / L*L;}
int TopRightRight(int i, int L){	return (TopLeft(i,L) + 2) % L + TopLeft(i,L) / L*L;}
int BotRight(int i, int L) { return (i - L + L*L) % (L*L); }
int BotRightRight(int i, int L){ return (BotRight(i,L) + 1) % L + BotRight(i,L) / L*L; }
int TopTopLeft(int i, int L){ return(i + 2 * L) % (L*L); }
int TopTop(int i, int L){ return(TopTopLeft(i,L) + 1) % L + TopTopLeft(i,L) / L*L; }

void make_nn(int &L, vector< vector<int> > &neighbors)
{
	vector<int>    nnentry;

	//////////////////////////////////////////////
	//             Triangular lattice
	//////////////////////////////////////////////
	// Make sites on a quadrilateral with length L on both sides.

	for (int i = 0; i < L*L; i++) {

		///////////////////////////////////////////////////////////////////////////////////////////////
		int right = (i + 1) % L+i/L*L;
		int left = (i - 1 + L) % L+i/L*L;
		int topLeft = (i + L) % (L*L);
		int topRight = (topLeft + 1) % L+topLeft/L*L;
		int botRight = (i - L + L*L) % (L*L);
		int botLeft = (botRight - 1 + L) % L+botRight/L*L;
		nnentry.push_back(right);
		nnentry.push_back(left);
		nnentry.push_back(topLeft);
		nnentry.push_back(topRight);
		nnentry.push_back(botRight);
		nnentry.push_back(botLeft);
		neighbors.push_back(nnentry);
		nnentry.clear();
	}


}

void make_3rdn(int &L, vector< vector<int> > &third_neighbors)
{
	vector<int>    nnentry;

	//////////////////////////////////////////////
	//             Pyrochlore lattice
	//////////////////////////////////////////////
	// Make sites on cube of dimensions L,L,L
	// Associate 16 sites with every n1,n2,n3
	//
	// Tetrahedron coords (used in Lucile/Ross paper)
	// r0 = 1/8 (+1,+1,+1)
	// r1 = 1/8 (+1,-1,-1)
	// r2 = 1/8 (-1,+1,-1)
	// r3 = 1/8 (-1,-1,+1)
	for (int i = 0; i < L*L; i++) {

		///////////////////////////////////////////////////////////////////////////////////////////////
		int right = (i + 2) % L + i / L*L;
		int left = (i - 2 + L) % L + i / L*L;
		int topLeft = (i + 2*L) % (L*L);
		int topRight = (topLeft + 2) % L + topLeft / L*L;
		int botRight = (i - 2*L + L*L) % (L*L);
		int botLeft = (botRight - 2 + L) % L + botRight / L*L;
		nnentry.push_back(right);
		nnentry.push_back(left);
		nnentry.push_back(topLeft);
		nnentry.push_back(topRight);
		nnentry.push_back(botLeft);
		nnentry.push_back(botRight);
		third_neighbors.push_back(nnentry);
		nnentry.clear();
	}


}


vector<vector<double>> disorder(int&L, int&num, int &nsites, double &alpha,vector<vector<int>> &neighbors)
{
	int dislo =nsites / num;
	vector<vector<double>> nn = { { 1,0 },{ -1,0 },{ -0.5,sqrt(3)*0.5 },{ 0.5,sqrt(3)*0.5 },{ 0.5,-sqrt(3)*0.5 },{ -0.5,-sqrt(3)*0.5 } };
	vector < vector<double >> bergus_lst (dislo,{ 0,0 } );
	vector<vector<double>> core_lst (dislo,{ 0,0} );
	for (int i = 0; i < dislo; i++) {
		bergus_lst[i] = nn[randint(6)];
		int site = randint(nsites);
		int s_i = site / L;
		int s_j = site%L;
		core_lst[i] = { -0.5*s_i + s_j,sqrt(3)*0.5*s_i };
	}
//	bergus_lst[0]=nn[0];
//	bergus_lst[1]=nn[1];
//	core_lst[0]={-0.5*L/2.0+L/4.0,sqrt(3)*0.5*L/2.0};
//	core_lst[1]={-0.5*L/2.0+3*L/4.0,sqrt(3)*0.5*L/2.0};
	
	
	vector<vector<double>> result(nsites, { 0,0,0,0,0,0 });
	for (int i = 0; i < nsites; i++) {
		for (int j = 0; j < 6; j++) {
			double dj = 0;
			for (int k = 0; k < dislo; k++) {
				vector<double> bergus = bergus_lst[k];
				vector<double> core = core_lst[k];
				int n = neighbors[i][j];
				int i_i = i / L;  //i,j of the first site i
				int i_j = i%L;
				int n_i = n / L;  //i,j of the second site n
				int n_j = n%L;
				vector<double> r_i = { -0.5*i_i + i_j,sqrt(3) / 2.0* i_i }; //position vector of site i
				vector<double> r_n = { -0.5*n_i + n_j,sqrt(3) / 2.0* n_i };  //position vector of site n
				vector<double> r_avg = { 0.5*(r_i[0] + r_n[0]),0.5*(r_i[1] + r_n[1]) }; //position vector of the center of the bond between site i and site n
				vector<double> dr = { r_avg[0] - core[0],r_avg[1] - core[1] };  //displacement vector between the core and the bond center
				double d = sqrt(dr[0] * dr[0] + dr[1] * dr[1]);
				double cos = (dr[0] * bergus[0] + dr[1] * bergus[1]) / d;
				double sin = (bergus[0] * dr[1] - bergus[1] * dr[0]) / d;
				double e = cos / (2 * 3.14159265*d);//we make dj not decaying.
				double exx = e*sin*cos;
				double eyy = -exx;
				double exy = e*(2 * cos*cos - 1) / 2.0;
				double nx = nn[j][0];
				double ny = nn[j][1];
				dj  =dj+ exx*nx*nx + eyy*ny*ny + 2 * exy*nx*ny;//this is dl;
			}
			result[i][j] = alpha*dj;  //choose J'=2;
		}
	}
	return result;
}
double nnenergy(int&nsites, int&L, double &J, vector<vector<double>> &config, vector< vector<int> > &neighbors, vector< vector<double> > &djm)
{


	double energycalc = 0.0;

	for (int i = 0; i<nsites; i++)
	{

		for (int k = 0; k < neighbors[i].size(); k++)
		{
			//cout <<i<<"	"<< neighbors[i][k] << endl;
			energycalc = energycalc + (djm[i][k] + J)*(double(config[i][0]) * double(config[neighbors[i][k]][0]) +
				double(config[i][1]) * double(config[neighbors[i][k]][1]) +
				double(config[i][2]) * double(config[neighbors[i][k]][2]));
			//cout<< double(config[neighbors[i][k]])<<endl;		
		}
	}
	return energycalc / (2.0);
}
double third_energy(int&nsites, int&L, double &Jthird, vector<vector<double>> &config, vector< vector<int> > &third_neighbors)
{


	double energycalc = 0.0;

	for (int i = 0; i<nsites; i++)
	{

		for (int k = 0; k < third_neighbors[i].size(); k++)
		{
			energycalc = energycalc + (Jthird)*(double(config[i][0]) * double(config[third_neighbors[i][k]][0]) +
				double(config[i][1]) * double(config[third_neighbors[i][k]][1]) +
				double(config[i][2]) * double(config[third_neighbors[i][k]][2]));

			//cout<< double(config[neighbors[i][k]])<<endl;		
		}
	}
	return energycalc / (2.0);
}

double diffenergy_nn(double &J, vector<vector<double>> &config,
	vector< vector<int> > &neighbors, int &chosensite, vector<double> new_old, vector< vector<double> > &djm)
{


	double energycalc = 0.0;
	for (int k = 0; k < neighbors[chosensite].size(); k++)
	{
		energycalc = energycalc + (djm[chosensite][k] + J)*(double(new_old[0]) * double(config[neighbors[chosensite][k]][0]) +
			double(new_old[1]) * double(config[neighbors[chosensite][k]][1]) +
			double(new_old[2]) * double(config[neighbors[chosensite][k]][2]));
	}
	return energycalc;
}
double diffenergy_third(double &Jthird, vector<vector<double>> &config,
	vector< vector<int> > &third_neighbors, int &chosensite, vector<double> new_old)
{


	double energycalc = 0.0;
	for (int k = 0; k < third_neighbors[chosensite].size(); k++)
	{
		energycalc = energycalc + Jthird*(double(new_old[0]) * double(config[third_neighbors[chosensite][k]][0]) +
			double(new_old[1]) * double(config[third_neighbors[chosensite][k]][1]) +
			double(new_old[2]) * double(config[third_neighbors[chosensite][k]][2]));
	}
	return energycalc;
}



int main(int argc, char *argv[])
{

	double J = 0.0;
	double Jthird = 0.0;
	double alpha = 0.0;
	long long int nsamples = 10;

	double T;
	int seed;
	int L;
	int num;


	string outfilename;
	string infilename;

	infilename = argv[1];

	read_file(infilename, L, J,alpha ,num, Jthird, T, nsamples, seed, outfilename);
	std::srand(seed);
//	cout<<L<<J<<alpha<<Jthird<<T<<endl;
	nsamples = (long long int)400000 * (long long int)(L*L);
	long long int nwates =0; //(long long int)200000 * (long long int)(L*L);
	int nsites = L*L;
	int pointNumber = 0;
	
	T = T*Jthird;
	double t=0.5*Jthird;	

	std::vector< std::vector<int> > neighbors;

	std::vector< std::vector<int> > third_neighbors;

	make_nn(L, neighbors);

	make_3rdn(L, third_neighbors);


	vector<vector<double>> djm = disorder(L,num, nsites,alpha, neighbors);
	ofstream out;
	string name = outfilename+to_string(seed);
	out.open(name);
	for (int i = 0; i < nsites; i++) { out << i << "=" << djm[i][0] << "=" << djm[i][1] << "=" << djm[i][2] << "=" << djm[i][3] << "=" << djm[i][4] << "=" << djm[i][5] << endl; }
	out.close();
	std::srand(seed);
	int count = 20;
	ofstream outfile[count];


	double energySum = 0.0;
	double energySquareSum = 0.0;


	//// Create configuration and measure energy

	vector<vector<double>> config = create_config(L);
	//for (int i=0;i<nsites;i++) cout<<config[i];

	//cout<<config.size()<<endl;
	//	cout << "start Metropolis" << endl;
	//	cout<<nsamples<<endl;
	// Run Metropolis
	double energycurrent = nnenergy(nsites, L, J, config, neighbors,djm) + third_energy(nsites, L, Jthird, config, third_neighbors);
	//cout << energycurrent/nsites << endl;
	
	for (long long int n = 0; n < nsamples; n++)

	{	if (n == 0) { T = t; }
		if (n == nsamples / 10) { T = t*0.9; }
		if (n == nsamples / 10*2) { T = t*0.8; }
		if (n == nsamples / 10*3) { T = t*0.7; }
		if (n == nsamples / 10*4) { T = t*0.6; }
		if (n == nsamples / 10*5) { T = t*0.5; }
		if (n == nsamples / 10*6) { T = t*0.4; }
		if (n == nsamples / 10*7) { T = t*0.3; }
		if (n == nsamples / 10*8) { T = t*0.2; }
		if (n == nsamples / 10*9) { T = t*0.1; }
	//	if (n == nsamples ) { T = t*0.02; }

		int chosensite = randint(nsites);
		double x1, x2;
		double new_x, new_y, new_z;
		while (true) {
			x1 = (ourrandom() - 0.5) *2.0;
			x2 = (ourrandom() - 0.5) *2.0;
			if (x1*x1 + x2*x2 < 1.0) break;
		}
		new_x = 2.0 * x1 * sqrt(1.0 - x1*x1 - x2*x2);
		new_y = 2.0 * x2 * sqrt(1.0 - x1*x1 - x2*x2);
		new_z = 1.0 - 2.0 * (1.0 - x1*x1 - x2*x2);

		vector<double >new_old =
		{ new_x - config[chosensite][0], new_y - config[chosensite][1], new_z - config[chosensite][2] };

		double energydiff = diffenergy_nn(J, config, neighbors, chosensite, new_old,djm) +
			diffenergy_third(Jthird, config, third_neighbors, chosensite, new_old);
	//	cout << energydiff << endl;
		double boltzmann = exp(-energydiff / T);
		double r = ourrandom();

		if (r < boltzmann)
		{
		//	cout << "flip" << endl;
			energycurrent = energycurrent + energydiff;
			config[chosensite][0] = new_x;
			config[chosensite][1] = new_y;
			config[chosensite][2] = new_z;
			//accept += 1;
		}
		else
		{
		//	reject += 1;
		}
		//}

		//}
		if (n>=nwates && n % ((nsamples-nwates) / count) == 0)
		{
			int num = (n-nwates) / ((nsamples-nwates) / count);
			string name = outfilename + to_string(num) + ".txt";
			outfile[num].open(name);
			for (int i = 0; i<nsites; i++) { outfile[num] << i << "=" << config[i][0] << "=" << config[i][1] << "=" << config[i][2] << endl; }
			outfile[num].close();
		}

		

	}
	
	

	//cout << "program ended." << endl;
	//cin.get();
	return 0;
	//      cout <<"program finished"<<endl;
}
