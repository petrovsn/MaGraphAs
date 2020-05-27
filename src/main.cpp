#include <iostream>
#include <stdio.h>
#include <vector>
#include <string>
#include <fstream>
#include <cstdio>
#include <time.h>
#include <map>
#include <set>
#include <math.h>
#include <unordered_map>
#include <tuple>
#include <math.h>
#include <sstream>
#include <thread>
#include <future>


#include <ctime>
#include <cstdio>

#include<Graph2.h>
#include <Interface.h>

using namespace std;


int main(int argc, const char* argv[])
{
	string REF_IN = ""; 
	string READS_IN = "";  
	string OUTPUT = "out";
	int THREADS = 0; 

	bool verbose = false;
	float treshhold = 0.8;
	int HASHLEN = 13;
	for (int i = 1; i < argc; i++)
	{
		cout << i << '\t' << argv[i] << endl;
		string str_argi = argv[i];
		if (str_argi == "-g")
		{
			REF_IN = argv[i + 1];
		}

		if (str_argi == "-f")
		{
			READS_IN = argv[i + 1];
		}
		
		if (str_argi == "-t")
		{
			THREADS = stoi(argv[i+1]);
		}

		if (str_argi == "-v")
		{
			verbose = true;
		}

		if (str_argi == "-c")
		{
			treshhold = stof(argv[i + 1]);
		}
		if (str_argi == "-o")
		{
			OUTPUT = argv[i + 1];
		}

		if (str_argi == "-l")
		{
			HASHLEN = stoi(argv[i+1]);
		}

	}	
	cout << "Ref:   " << REF_IN << endl;
	cout << "Reads: " << READS_IN << endl;
	cout << "Theads:" << THREADS << endl;
	cout << "Tresh: " << treshhold << endl;

	Graph GMAIN;

	GMAIN.LoadFromGFA(REF_IN);

	/*int size = 0;
	int edges = 0;
	for (auto k : GMAIN.Body)
	{
		size = size + GMAIN.Body[k.first].str.size();
		edges = edges + GMAIN.Body[k.first].Next.size();
	}
	cout<<"edge:"<<edges<<endl;
	cout<<GMAIN.Body.size()<<endl;
	cout<<"size:"<<size<<endl;
	cout<<"edge:"<<edges / GMAIN.Body.size()<<endl;*/


	GMAIN.BuildIndex(5, HASHLEN);
	
	clock_t time_a = clock();
	Results res = Align2GMAIN(GMAIN, READS_IN, treshhold, THREADS, verbose);
	clock_t time_b = clock();
	
	//cout << "Aligned: " << res.sumAlignment.size() << endl;
	//cout << "Avgtime: " << double(time_b-time_a)/CLOCKS_PER_SEC << endl;

	ofstream fout = ofstream(OUTPUT+".vcf");
	for (auto var: res.VarMap)
	{
		fout << var.first.ID1 << '.' << var.first.pos1 << '.' << var.first.ID2 << '.' << var.first.pos2 << '\t' << var.first.alt<<'\t'<<var.second << '\n';
	}
	fout.close();


	fout = ofstream(OUTPUT+".log");
	fout << "Ref:"<< '\t' << REF_IN << endl;
	fout << "Reads:"<< '\t' << READS_IN << endl;
	fout << "Theads:"<< '\t' << THREADS << endl;
	fout << "Tresh:"<< '\t' << treshhold << endl;
	fout << "Aligned:"<< '\t' << res.sumAlignment.size() << endl;
	//fout << "Avgtime:"<< '\t' << res.time_avg << endl;
	fout.close();

}
