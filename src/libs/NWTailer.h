#pragma once


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
#include <algorithm>

#include <GenerHash.h>
#include <Node.h>
#include <BubbleIndex.h>

#include <NWAligner.h>


class NWTtracke_back
{
public:
	bool ready = false;
	bool loaded = false;
	bool terminal = false;

	NWTtracke_back(NWNode prot, int pos, string read)
	{

	}

	void NWTInit_back()
	{

	}
};


class NWTtrack
{
public:
	int pos1 = -1;
	int pos2 = -1;
	
	vector<pair<int, int>> RealCords;
	vector<vector<int>>NWTmatr;

	int NextID = -1;
	int i_start = -1;

	int currmax;
	int pos = -1;
	bool finished = false;
	bool error = false;

	string storedref = "";

	NWTtrack GetCopy(string read, int N_Id, int nextpos = 0)
	{
		NWTtrack tmp = NWTtrack(read, N_Id, nextpos);
		tmp.pos1 = pos1;
		tmp.pos2 = pos2;
		tmp.NWTmatr = NWTmatr;
		tmp.NextID = N_Id;
		tmp.storedref = storedref;
		tmp.pos = nextpos;
		tmp.currmax = currmax;
		tmp.RealCords = RealCords;
		return tmp;
	}


	NWTtrack(string &read, int Next, int nextpos = 0)
	{
		NWTmatr.push_back(vector<int>(read.size() + 1, 0));
		for (int j = 1; j < read.size() + 1; j++)
		{
			NWTmatr[0][j] = NWTmatr[0][j - 1] + WGlob.gap;
		}
		RealCords.push_back(pair<int, int>(-1, -1));
		currmax = NWTmatr[0].back();
		NextID = Next;
		this->pos = nextpos;
	}

	void Next_forward(vector<NWTtrack> &alltracks, Node &n, string &read)
	{
		int rlen = read.size();
		int reflen = n.str.size();

		int tmp_max = currmax;
		bool exit = false;

		int i = NWTmatr.size() - 1;
		int refpos = pos - 1;

		if (n.ID == -1)
		{
			finished = true;
			return;
		}
		do
		{
			if (NWTmatr.size() > rlen * 3)
			{
				error = true;
				return;

			}

			i = i + 1;
			refpos = refpos + 1;
			NWTmatr.push_back(vector<int>(rlen + 1, 0));
			NWTmatr[i][0] = NWTmatr[i - 1][0];

			RealCords.push_back(pair<int, int>(n.ID, refpos));

			for (int j = 1; j < rlen + 1; j++)
			{
				vector<int> vals = {
						NWTmatr[i - 1][j - 1] + WGlob.CharMatrix[(n.str[refpos])][read[j - 1]],
						NWTmatr[i - 1][j] + WGlob.gap ,
						NWTmatr[i][j - 1] + WGlob.gap };

				vector<int>::iterator res = max_element(vals.begin(), vals.end());
				int tmpidx = distance(vals.begin(), res);
				NWTmatr[i][j] = vals[tmpidx];
			}

			tmp_max = NWTmatr[i].back();

			if (tmp_max >= currmax)
			{
				currmax = tmp_max;
			}
			else
			{
				i_start = i - 1;
				refpos = refpos - 1;
				exit = true;
				finished = true;
				RealCords.pop_back();
				NWTmatr.pop_back();
			}

			if (refpos == reflen - 1)
			{
				exit = true;
			}
		} while (not exit);

		storedref = storedref + n.str.substr(pos, refpos - pos + 1);


		vector<NWTtrack> tmp_arr;

		if (not finished)
		{
			if (n.Next.size() == 1)
			{
				NextID = n.Next.begin()->second->ID;
				pos = 0;
			}


			else if (n.Next.size() > 1)
			{
				NextID = n.Next.begin()->second->ID;
				pos = 0;
				for (auto n_tmp : n.Next)
				{
					NWTtrack nw = GetCopy(read, n_tmp.second->ID, 0);
					if (nw.NextID != this->NextID)
					{
						tmp_arr.push_back(nw);
					}
				}
			}

			else if (n.Next.size() == 0)
			{
				finished = true;
			}
		}

		alltracks.insert(alltracks.end(), tmp_arr.begin(), tmp_arr.end());
	}

	void Next_backward(vector<NWTtrack>& alltracks, Node& n, string& read)
	{
		int rlen = read.size();
		int reflen = n.str.size();

		int tmp_max = currmax;
		bool exit = false;

		int i = NWTmatr.size() - 1;
		int refpos = pos+1;

		do
		{
			if (NWTmatr.size() > rlen * 3)
			{
				error = true;
				return;
			}
			i = i + 1;
			refpos = refpos - 1;
			NWTmatr.push_back(vector<int>(rlen + 1, 0));
			NWTmatr[i][0] = NWTmatr[i - 1][0] + WGlob.gap;

			RealCords.push_back(pair<int, int>(n.ID, refpos));

			for (int j = 1; j < rlen + 1; j++)
			{
				vector<int> vals = {
						NWTmatr[i - 1][j - 1] + WGlob.CharMatrix[(n.str[refpos])][read[rlen-j]],
						NWTmatr[i - 1][j] + WGlob.gap ,
						NWTmatr[i][j - 1] + WGlob.gap };

				vector<int>::iterator res = max_element(vals.begin(), vals.end());
				int tmpidx = distance(vals.begin(), res);
				NWTmatr[i][j] = vals[tmpidx];
			}

			tmp_max = NWTmatr[i].back();

			if (tmp_max >= currmax)
			{
				currmax = tmp_max;
			}
			else
			{
				i_start = i-1;
				refpos = refpos + 1;
				exit = true;
				finished = true;
				RealCords.pop_back();
				NWTmatr.pop_back();
			}

			if (refpos == 0)
			{
				exit = true;
			}
		} while (not exit);

		storedref = n.str.substr(refpos, pos - refpos +1) + storedref;
		vector<NWTtrack> tmp_arr;

		if (not finished)
		{
			if (n.Prev.size() == 1)
			{
				NextID = n.Prev.begin()->second->ID;
				pos = n.Prev.begin()->second->str.size()-1;
			}


			else if (n.Prev.size() > 1)
			{
				NextID = n.Prev.begin()->second->ID;
				pos = n.Prev.begin()->second->str.size()-1;
				for (auto n_tmp : n.Prev)
				{
					NWTtrack nw = GetCopy(read, n_tmp.second->ID, n_tmp.second->str.size()-1);
					if (nw.NextID != this->NextID)
					{
						tmp_arr.push_back(nw);
					}
				}
			}

			else if (n.Prev.size() == 0)
			{
				finished = true;
			}
		}
		alltracks.insert(alltracks.end(), tmp_arr.begin(), tmp_arr.end());
	}

	
	
	NWpart NWTtrace(string &read)
	{
		vector<GVariation> vars;
		int i_pos = i_start;
		int j_pos = read.size();
		if (i_start == -1)
		{
			return  NWpart();
		}
		//string alignedref = storedref.substr(i_pos - 1, 1);
		//string alignedread = read.substr(j_pos-1, 1);

		string alignedref = "";
		string alignedread = "";

	//	cout << "toalign:" << endl;
	//	cout << " ref: " << storedref << endl;
	//	cout << "read: " << read << endl;
		
		GVariation lastvar;

		NWpart result;

		result.ID1 = RealCords[1].first;
		result.pos1 = RealCords[1].second;

		result.ID2 = RealCords.back().first;
		result.pos2 = RealCords.back().second;

		int score = 0;
		int penalted_score = 0;

		int lasalignedpos = 0;
		do
		{
			vector<int> vals = { NWTmatr[i_pos - 1][j_pos - 1], //match
								 NWTmatr[i_pos - 1][j_pos], //left
								 NWTmatr[i_pos][j_pos - 1] }; //up
			int tmpidx = 0;


			if ((i_pos == 1) and (j_pos > 1))
			{
				tmpidx = 2;
			}
			else if ((i_pos > 1) and (j_pos == 1))
			{
				tmpidx = 1;
			}
			else
			{
				vector<int>::iterator res = max_element(vals.begin(), vals.end());
				tmpidx = distance(vals.begin(), res);
			}


			if (tmpidx == 0)
			{
				alignedref = storedref[i_pos - 1] + alignedref;
				alignedread = read[j_pos - 1] + alignedread;

				if (tolower(storedref[i_pos - 1]) == tolower(read[j_pos - 1]))
				{
					lasalignedpos = i_pos;
					score = score + 1;
					penalted_score = penalted_score + 1;
					if (lastvar.active) //if variation is over
					{
						lastvar.Finish(RealCords[i_pos].first, RealCords[i_pos].second);
						if (lastvar.ID2 != -1)
						{
							vars.push_back(lastvar);
						}
					}
				}
				else
				{
					penalted_score = penalted_score - 1;
					if (lastvar.active) //if variation active extende, else - start new and extend. 
					{
						lastvar.Extend(read[j_pos - 1]);
					}
					else
					{
						lastvar = GVariation(RealCords[lasalignedpos].first, RealCords[lasalignedpos].second);
						lastvar.Extend(read[j_pos - 1]);
					}
				}
				i_pos = i_pos - 1;
				j_pos = j_pos - 1;



			}
			else if (tmpidx == 1)
			{
				penalted_score = penalted_score - 1;
				if (lastvar.active) //if variation active extende, else - start new and extend. 
				{
					lastvar.Extend();
				}
				else
				{
					lastvar = GVariation(RealCords[lasalignedpos].first, RealCords[lasalignedpos].second);
					lastvar.Extend();
				}

				alignedref = storedref[i_pos - 1] + alignedref;
				alignedread = '-' + alignedread;

				i_pos = i_pos - 1;
			}
			else if (tmpidx == 2)
			{
				penalted_score = penalted_score - 1;
				if (lastvar.active) //if variation active extende, else - start new and extend. 
				{
					lastvar.Extend(read[j_pos - 1]);
				}
				else
				{
					lastvar = GVariation(RealCords[lasalignedpos].first, RealCords[lasalignedpos].second);
					lastvar.Extend(read[j_pos - 1]);
				}



				alignedread = read[j_pos - 1] + alignedread;
				alignedref = '-' + alignedref;
				j_pos = j_pos - 1;
			}


		} while ((i_pos != 0) and (j_pos != 0));
		
		result.vars = vars;
		result.score = score;
		result.penalted_score = penalted_score;
		result.alnread = alignedread;
		result.alnref = alignedref;
	//	cout << "score: " << score << endl;
	//	cout << "p.scr: " << penalted_score << endl;
	//	cout << "  ref: " << alignedref << endl;
	//	cout << " read: "<< alignedread << endl;


		return result;
	}

};



class NWTailer
{
public:
	//first letter in read is fixed
	NWTailer()
	{
		
	}


	NWpart NWTforw(int ID, int pos, map<int, Node>& Body, string read)
	{
		NWTtrack startnwt = NWTtrack(read, ID, pos);
		vector< NWTtrack> total_nwts = { startnwt };

		int Tots = total_nwts.size();
		bool exit = false;

		do
		{
			for (int i = 0; i < Tots; i++)
			{
				if (not total_nwts[i].finished)
				{
					total_nwts[i].Next_forward(total_nwts,
						Body[total_nwts[i].NextID],
						read);
				}
			}

			exit = true;
			for (int i = 0; i < total_nwts.size(); i++)
			{
				if ((not total_nwts[i].finished)and(not total_nwts[i].error))
				{
					exit = false;
				}
			}

			Tots = total_nwts.size();
		}
		while (not exit);
		

		int score = total_nwts[0].currmax;
		int bestid = 0;
		for (int i = 0; i < total_nwts.size(); i++)
		{
			if (total_nwts[i].currmax > score)
			{
				score = total_nwts[i].currmax;
				bestid = i;
			}
		}
			
		NWpart res = total_nwts[bestid].NWTtrace(read);

		return res;
	}

	NWpart NWTback(int ID, int pos, map<int, Node>& Body, string read)
	{
		NWTtrack startnwt = NWTtrack(read, ID, pos);
		vector< NWTtrack> total_nwts = { startnwt };

		int Tots = total_nwts.size();
		bool exit = false;

		do
		{
			for (int i = 0; i < Tots; i++)
			{
				if (not total_nwts[i].finished)
				{
					total_nwts[i].Next_backward(total_nwts,
						Body[total_nwts[i].NextID],
						read);
				}
			}

			exit = true;
			for (int i = 0; i < total_nwts.size(); i++)
			{
				if ((not total_nwts[i].finished) and (not total_nwts[i].error))
				{
					exit = false;
				}
			}

			Tots = total_nwts.size();
		} while (not exit);


		int score = total_nwts[0].currmax;
		int bestid = 0;
		for (int i = 0; i < total_nwts.size(); i++)
		{
			if (total_nwts[i].currmax > score)
			{
				score = total_nwts[i].currmax;
				bestid = i;
			}
		}



		string revread = read;
		reverse(revread.begin(), revread.end());

		reverse(total_nwts[bestid].storedref.begin(), total_nwts[bestid].storedref.end());
		NWpart res = total_nwts[bestid].NWTtrace(revread);

		int tmp;
		swap(res.ID1, res.ID2);
		swap(res.pos1, res.pos2);
		for (int i = 0; i < res.vars.size(); i++)
		{
			swap(res.vars[i].ID1, res.vars[i].ID2);
			swap(res.vars[i].pos1, res.vars[i].pos2);
			reverse(res.vars[i].alt.begin(), res.vars[i].alt.end());
		}

		return res;
	}
};