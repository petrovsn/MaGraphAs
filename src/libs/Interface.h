#pragma once
#include <iostream>
#include <stdio.h>
#include <vector>
#include <string>
#include <fstream>
#pragma once
#include <cstdio>
#include <time.h>
#include <map>
#include <set>
#include <math.h>
#include <unordered_map>
#include <tuple>
#include <math.h>
#include <sstream>
#include <Graph2.h>
#include <Fasta.h>
#include <thread>
#include <future>

using namespace std;

class Results
{
public:
	unordered_map<GVariation, vector<string>, MyHashFunction> VarMapStringCount;

	unordered_map<GVariation, int, MyHashFunction> VarMap;
	map<string, FAlignment> sumAlignment;

	clock_t time_avg = 0;


	Results() {};
	int nonaligned = 0;
	void Add(Results res2)
	{
		for (auto n : res2.VarMap)
		{
			VarMap[n.first] = VarMap[n.first] + n.second;
		}

		for (auto n : res2.sumAlignment)
		{
			sumAlignment[n.first] = n.second;
		}

		for (auto n : res2.VarMapStringCount)
		{
			for (int i = 0; i < n.second.size(); i++)
			{
				VarMapStringCount[n.first].push_back(n.second[i]);
			}
		}

		time_avg = (time_avg + res2.time_avg) / 2;
	}
};

class Worker
{
public:
	float treshold;
	int ID;
	future<Results> result;
	string curr_read;
	bool verbose_out;

	Worker(int nID, float cov, bool verbose)
	{
		treshold = cov;
		ID = nID;
		verbose_out = verbose;
	}

	void Run(Graph *GMAIN, ifstream *fin, mutex *mut)
	{
		result = async(launch::async, &Worker::ThreadTask, this, GMAIN, fin, mut);
	}

	Results ThreadTask(Graph *GMAIN, ifstream *fin, mutex *mut)
	{
		Results res;
		string forread = "";
		string header = "";


		mut->lock();
		
		getline(*fin, header);
		getline(*fin, forread);

		if (verbose_out)
		{
			if (forread == "")
			{
				cout << "Nullread: " << header << endl;
			}
			else
			{
				cout << "Thread:"<<ID<< header << endl;
				cout << forread << endl;
			}
		}


		if (fin->eof())
		{
			mut->unlock();
			return res;
		}
		mut->unlock();

		int read_counter = 1;
		clock_t avg_worktime = 0;

		do
		{

			string revread = reverse(forread);

			curr_read = forread;

			//clock_t time_a = clock();
			FAlignment aln = GMAIN->AlignHashSmWtmn(forread); 
			FAlignment aln_rev = GMAIN->AlignHashSmWtmn(revread);
			//clock_t time_b = clock();

			//avg_worktime = avg_worktime + (time_b - time_a);
			read_counter = read_counter + 1;

			if (aln_rev.score > aln.score)
			{
				aln = aln_rev; 
			}

			if (aln.tresh > treshold)
			{
				for (int j = 0; j < aln.vars.size(); j++)
				{
					GVariation key = aln.vars[j];
					res.VarMap[key] = res.VarMap[key] + 1;
					res.VarMapStringCount[key].push_back(header);
				}

				res.sumAlignment[header] = aln;
			}
			else
			{
				res.nonaligned = res.nonaligned + 1;
			}


			mut->lock();
			getline(*fin, header);
			getline(*fin, forread);

			if (verbose_out)
			{
				if (forread == "")
				{
					cout << "Thread:" << ID << "Nullread: " << header << endl;
					mut->unlock();
					continue;
				}
				else
				{
					cout << "Thread:" << ID << header << endl;
					cout << forread << endl;
				}
			}

			mut->unlock();

		} while (!fin->eof());

		//avg_worktime = avg_worktime / read_counter;
		res.time_avg = 0; //avg_worktime;

		return res;
	}

	Results Finalize()
	{
		future_status res_stat = result.wait_for(chrono::microseconds(5000));
		if (res_stat == future_status::ready)
		{
			cout<<"ID = "<<ID<<endl;
			cout<<"status: ready"<<endl;
			return result.get();
		}
		if (res_stat == future_status::timeout)
		{
			cout<<"ID = "<<ID<<endl;
			cout<<"status: timeout"<<endl;
			cout<<"Read: "<<curr_read<<endl;
			result.wait();
			return result.get();
		}
		if (res_stat == future_status::deferred)
		{
			cout<<"ID = "<<ID<<endl;
			cout<<"status: deferred"<<endl;
			Results res;
			return res;
		}
	}

};


class ThreadPool
{
public:
	ifstream fin;
	vector<Worker> threads;
	mutex mut;
	mutex outmut;

	ThreadPool(string input, float cov, int numthreads, bool verbose)
	{
		fin = ifstream(input);

		for (int i = 0; i < numthreads; i++)
		{
			threads.push_back(Worker(i,cov, verbose));
		}
		cout<<"TPool construct: numthreads = "<<threads.size()<<endl;
	}

	void Run(Graph &VG)
	{
		for (int i = 0; i < threads.size(); i++)
		{
			threads[i].Run(&VG, &fin, &mut);
		}
	}

	Results GetResults()
	{
		Results res;

		for (int i = 0; i < threads.size(); i++)
		{
			res.Add(threads[i].Finalize());
		}
		return res;
	}

	
};


Results Align2GMAIN(Graph &GRAPH, string input_reads, float treshhold, int numthreads, bool verbose)
{
	cout << "SRAAlign: numthreads " << numthreads << endl;
	ThreadPool TPool(input_reads, treshhold, numthreads, verbose);

	cout<<"SRAAlign: TPool created"<<endl;
	TPool.Run(GRAPH);
	cout<<"SRAAlign: TPool runned"<<endl;
	Results res = TPool.GetResults();

	return res;
}
