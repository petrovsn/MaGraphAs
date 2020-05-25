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




void printNWmart(string ref, string read, vector<vector<int>> NWMatrix1)
{
	cout << ref.size() << '/' << NWMatrix1.size() << endl;
	cout << read.size() << '/' << NWMatrix1[0].size() << endl;
	cout << '\t' << '\t';
	for (int i = 0; i < ref.size(); i++)
	{
		cout << ref[i] << '\t';
	}
	cout << endl;
	cout << '\t';
	for (int i = 0; i < NWMatrix1.size(); i++)
	{
		cout << NWMatrix1[i][0] << '\t';
	}
	cout << endl;
	for (int j = 1; j <= read.size(); j++)
	{
		cout << read[j - 1] << '\t';
		for (int i = 0; i < NWMatrix1.size(); i++)
		{
			cout << NWMatrix1[i][j] << '\t';
		}
		cout << endl;
	}
}

class Weights
{
public:
	map<char, map<char, int>> CharMatrix;
	int gap;

	Weights()
	{
		CharMatrix = InitCharMatrix();
		gap = -6;
	}

	map<char, map<char, int>> InitCharMatrix()
	{
		map<char, map<char, int>> CharMatrix;
		for (char c1 : { 'a', 'c', 'g', 't', 'A', 'C', 'G', 'T' })
		{
			for (char c2 : { 'a', 'c', 'g', 't', 'A', 'C', 'G', 'T' })
			{
				if (tolower(c1) == tolower(c2))
				{
					CharMatrix[c1][c2] = 0;
				}
				else
				{
					CharMatrix[c1][c2] = -2;
				}
			}
		}
		return CharMatrix;
	}
};

Weights WGlob = Weights();

class GVariation
{
public:

	int ID1;
	int pos1;
	int ID2;
	int pos2;

	string alt;
	bool active = false;

	GVariation()
	{

	}
	
	GVariation(int eID, int epos)
	{
		ID2 = eID;
		pos2 = epos;
		active = true;
	}

	void Extend(char c)
	{
		alt = c + alt;
	}

	void Extend()
	{
	}

	void Finish(int sID, int spos)
	{
		ID1 = sID;
		pos1 = spos;
		active = false;
	}

	bool operator==(const GVariation& other) const
	{
		return ((ID1 == other.ID1) and
			(pos1 == other.pos1) and
			(ID2 == other.ID2) and
			(pos2 == other.pos2) and
			(alt == other.alt));
	}
};


template <class T>
inline void hash_combine(std::size_t& s, const T& v)
{
	std::hash<T> h;
	s ^= h(v) + 0x9e3779b9 + (s << 6) + (s >> 2);
}


class MyHashFunction {
public:

	// Use sum of lengths of first and last names 
	// as hash function. 
	std::size_t operator()(const GVariation& k) const
	{
	/*	size_t res = 17;
		res = res * 31 + hash<int>()(k.ID1);
		res = res * 31 + hash<int>()(k.pos1);
		res = res * 31 + hash<int>()(k.ID2);
		res = res * 31 + hash<int>()(k.pos2);
		res = res * 31 + hash<string>()(k.alt);
		return res;*/


		std::size_t res = 0;
		hash_combine(res, k.ID1);
		hash_combine(res, k.pos1);
		hash_combine(res, k.ID2);
		hash_combine(res, k.pos2);
		hash_combine(res, k.alt);
		return res;
	}
};


struct NWpart
{
	vector<GVariation> vars;
	GVariation lastvar;

	string alnread = "";
	string alnref = "";

	int ID1 = -1;
	int pos1 = -1;
	int ID2 = -1;
	int pos2 = -1;

	int lastloc = -1;
	int score = 0;
	int penalted_score = 0;
	int ID = -1;

	int lastalignedID = -1;
	int lastalignedpos = -1;

	bool finished = false;
};




class NWNode
{
public:
	
	string ref;
	int ID;
	int startpos = 0;

	int max_path = 0;
	int min_path = 0;

	vector<int> NWNext;
	
	bool status = false;

	NWNode()
	{
		ID = -1;
		ref = "empty";
	}
	 
	NWNode(string &reference, int newID)
	{
		ID = newID;
		ref = reference;
	}


	bool CheckSize(int readlen, int Minlen, int Maxlen, int MaxInDels = 25)
	{
		bool check1 = false;

		if (readlen * 3 + MaxInDels < Maxlen)
		{
			check1 = true;
		}

		bool check2 = false;

		if (readlen * 2 + MaxInDels < Minlen)
		{
			check2 = true;
		}

		if (check1 or check2)
		{
			return false;
		}
		return true;
	}



	//            data vectors          |       ------positions------     |
	NWNode(Node &N, vector<int>& links, map<int, Node> &Body, int N1, int pos1, int N2, int pos2, int lenread, map<int, pair<int, int>> &MinMaxborders, map<int, vector<int>>& NextNodes, bool first = false)
	{
		bool cont = true;
		if (first)
		{
			startpos = pos1;
		}

		ID = N.ID;

		if ((ID == N1) and (ID == N2))//in one node
		{
			int size = pos2 - pos1 + 1;
			

			MinMaxborders[ID] = pair<int, int>(MinMaxborders[ID].first + size, MinMaxborders[ID].second + size);
			bool Check = CheckSize(lenread, MinMaxborders[ID].first, MinMaxborders[ID].second);
			if (!Check)
			{
				return;
			}

			ref = N.str.substr(pos1, size);
			status = true;
			cont = false;

		}
		else if ((ID == N1) and (ID != N2))//first node of subgraph
		{
			int size = N.str.size() - pos1;

			MinMaxborders[ID] = pair<int, int>(MinMaxborders[ID].first + size, MinMaxborders[ID].second + size);
			bool Check = CheckSize(lenread, MinMaxborders[ID].first, MinMaxborders[ID].second);
			if (!Check)
			{
				return;
			}


			ref = N.str.substr(pos1);
			status = true;
		}
		else if ((ID != N1) and (ID == N2))//last node of subraph
		{

			int size = pos2 + 1;
			MinMaxborders[ID] = pair<int, int>(MinMaxborders[ID].first + size, MinMaxborders[ID].second + size);
			bool Check = CheckSize(lenread, MinMaxborders[ID].first, MinMaxborders[ID].second);
			if (!Check)
			{
				return;
			}

			ref = N.str.substr(0,pos2+1);
			status = true;
			cont = false;
		}
		else
		{
			int size = N.str.size();
			MinMaxborders[ID] = pair<int, int>(MinMaxborders[ID].first + size, MinMaxborders[ID].second + size);
			bool Check = CheckSize(lenread, MinMaxborders[ID].first, MinMaxborders[ID].second);
			if (!Check)
			{
				return;
			}

			ref = N.str;
			status = true;
		}


		if (cont)
		{

			for (auto p1 : N.Next)
			{
				bool internode = (*(p1.second)).IsInside(N1, N2, NextNodes);
				bool finalnode = ((*(p1.second)).ID == N2);

				if (internode or finalnode)
				{
					links.push_back((*(p1.second)).ID);
					NWNext.push_back((*(p1.second)).ID);
					MinMaxborders[(*(p1.second)).ID] = MinMaxborders[ID];
				}
			}
		}
	}

	void Link()
	{

	}

	void MinMaxMark(map<int, NWNode> &NWBody, vector<int>& links)
	{
		max_path = max_path + ref.size();
		min_path = min_path + ref.size();

		for (int i = 0; i < NWNext.size(); i++)
		{
			if (NWBody[NWNext[i]].max_path  <  (max_path + NWBody[NWNext[i]].ref.size()))
			{
				NWBody[NWNext[i]].max_path = max_path;
			}


			if ((NWBody[NWNext[i]].min_path > (min_path + NWBody[NWNext[i]].ref.size())) and
				(NWBody[NWNext[i]].min_path != -1))
			{
				NWBody[NWNext[i]].min_path = min_path;
			}
			else if (NWBody[NWNext[i]].min_path == -1)
			{
				NWBody[NWNext[i]].min_path = min_path;
			}

			links.push_back(NWNext[i]);
		}
	}

	vector<int> FirstIds;
	vector<vector<int>> NWMatrix;
	bool NwMatrLoaded = false;

	void InitMatrix(string &read, map<int, vector<int>> prevcolmns = map<int, vector<int>>())
	{
		FirstIds = vector<int>();
		
		int Ilen = ref.size() + 1;
		int Jlen = read.size() + 1;

		NWMatrix = vector<vector<int>>(Ilen);

		for (int i = 0; i < NWMatrix.size(); i++)
		{
			NWMatrix[i] = vector<int>(Jlen);
			for (int j = 0; j < NWMatrix[i].size(); j++)
			{
				NWMatrix[i][j] = 0;
			}
		}

		//fill first column
		if (prevcolmns.size() > 0) 
		{
			vector<int> previds;

			for (map<int, vector<int>>::iterator it = prevcolmns.begin(); it != prevcolmns.end(); ++it)
			{
				previds.push_back(it->first);
			}

			for (int j = 0; j < read.size() + 1; j++)
			{
				int curr_id = previds[0];
				int curr_max = prevcolmns[curr_id][j];

				for (int k = 0; k < previds.size(); k++)
				{
					if (prevcolmns[previds[k]][j] > curr_max)
					{
						curr_id = previds[k];
						curr_max = prevcolmns[previds[k]][j];
					}
				}

				FirstIds.push_back(curr_id);
				NWMatrix[0][j] = curr_max;
			}
		}  //fill first
		else
		{
			NWMatrix[0][0] = 0;
			for (int j = 1; j < NWMatrix[0].size(); j++)
			{
				NWMatrix[0][j] = NWMatrix[0][j - 1] + WGlob.gap;
			}
		}

		//fill first raw
		for (int i = 1; i < NWMatrix.size(); i++)
		{
			NWMatrix[i][0]=NWMatrix[i-1][0]+ WGlob.gap;
		}

		//fill other matrix
		for (int i = 1; i < NWMatrix.size(); i++)
		{
			for (int j = 1; j < NWMatrix[i].size(); j++)
			{
				vector<int> vals = { NWMatrix[i - 1][j - 1] + WGlob.CharMatrix[ref[i]][read[j]],
					NWMatrix[i - 1][j] + WGlob.gap ,
					NWMatrix[i][j - 1] + WGlob.gap };

				vector<int>::iterator res = max_element(vals.begin(), vals.end());
				int tmpidx = distance(vals.begin(), res);
				NWMatrix[i][j] = vals[tmpidx];
			}
		}
		NwMatrLoaded = true;
	}

	void AddVar(char readsimb, NWpart &head)
	{

	}

	NWpart NWTrace(string &read, NWpart head)
	{
		int i_pos = NWMatrix.size() - 1;
		int j_pos = NWMatrix.back().size() - 1; 
		int score = head.score;
		int penalted_score = head.penalted_score;

		if (head.lastloc!=-1)
		{
			j_pos = head.lastloc;
		}
		
		bool first = false;

		if (head.lastalignedID == -1)
		{
			first = true;
			head.lastalignedID = ID;
			head.lastalignedpos = ref.size() - 1;
		}

		string alignedref = head.alnref;
		string alignedread = head.alnread;

		//printNWmart(ref, read, NWMatrix);
		if (NWMatrix.size() == 1)
		{
			
		}
		else
		{
			do
			{
				vector<int> vals = { NWMatrix[i_pos - 1][j_pos - 1], //match
									 NWMatrix[i_pos - 1][j_pos], //left
									 NWMatrix[i_pos][j_pos - 1] }; //up
				int tmpidx = 0;


				if ((first) and (i_pos == NWMatrix.size() - 1))
				{
					tmpidx = 0;
				}
				else if ((i_pos == 1) and (j_pos > 1))
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
					alignedref = ref[i_pos - 1] + alignedref;
					alignedread = read[j_pos - 1] + alignedread;

					if (tolower(ref[i_pos - 1]) == tolower(read[j_pos - 1]))
					{
						score = score + 1;
						penalted_score = penalted_score + 1;
						head.lastalignedID = ID;
						head.lastalignedpos = i_pos - 1 + startpos;

						if (head.lastvar.active) //if variation is over
						{
							head.lastvar.Finish(ID, i_pos - 1 + startpos);
							head.vars.push_back(head.lastvar);
						}
					}
					else
					{
						penalted_score = penalted_score - 1;
						if (head.lastvar.active) //if variation active extende, else - start new and extend. 
						{
							head.lastvar.Extend(read[j_pos - 1]);
						}
						else
						{
							head.lastvar = GVariation(head.lastalignedID, head.lastalignedpos);
							head.lastvar.Extend(read[j_pos - 1]);
						}
					}
					i_pos = i_pos - 1;
					j_pos = j_pos - 1;



				}
				else if (tmpidx == 1)
				{
					penalted_score = penalted_score - 1;
					if (head.lastvar.active) //if variation active extende, else - start new and extend. 
					{
						head.lastvar.Extend();
					}
					else
					{
						head.lastvar = GVariation(head.lastalignedID, head.lastalignedpos);
						head.lastvar.Extend();
					}

					alignedref = ref[i_pos - 1] + alignedref;
					alignedread = '-' + alignedread;

					i_pos = i_pos - 1;
				}
				else if (tmpidx == 2)
				{
					penalted_score = penalted_score - 1;
					if (head.lastvar.active) //if variation active extende, else - start new and extend. 
					{
						head.lastvar.Extend(read[j_pos - 1]);
					}
					else
					{
						head.lastvar = GVariation(head.lastalignedID, head.lastalignedpos);
						head.lastvar.Extend(read[j_pos - 1]);
					}



					alignedread = read[j_pos - 1] + alignedread;
					alignedref = '-' + alignedref;
					j_pos = j_pos - 1;
				}


			} while ((i_pos != 0) and (j_pos != 0));
		}

		if (j_pos == 0)
		{
			head.finished = true;
		}

		head.score = score;
		head.penalted_score = penalted_score;

		head.alnread = alignedread;
		head.alnref = alignedref;
		head.lastloc = j_pos;
		

		if (FirstIds.size() != 0)
		{
			head.ID = FirstIds[j_pos];
		}
		return head;
	}
	 

	
};

class NWAligner
{
public:
	map<int, NWNode> NWBody;
	int lastID = -1;
	bool loaded = false;

	int nwID1 = -1;
	int nwpos1 = -1;
	int nwID2 = -1;
	int nwpos2 = -1;

	NWAligner()
	{

	}
	void Init(int ID1, int pos1, int ID2, int pos2, map<int, Node>& Body, string& read, map<int, vector<int>>& NextNodes, int MaxInDels = 25)
	{
		loaded = false;
		NWBody.clear();
		lastID = ID2;

		//Construct subgraph;

		map<int, pair<int, int>> MinMaxborders;
		MinMaxborders[ID1] = pair<int, int>(0, 0);

		vector<int> links;
		NWNode startNode = NWNode(Body[ID1], links, Body, ID1, pos1, ID2, pos2, read.size(), MinMaxborders, NextNodes, true);
		map<int, vector<int>> prevIds;
		if (!startNode.status) return;
		

		for (int j = 0; j < startNode.NWNext.size(); j++)
		{
			int t = startNode.NWNext[j];
			prevIds[t].push_back(startNode.ID);
		}

		NWBody[ID1] = startNode;

		while (links.size() != 0)
		{
			vector<int> new_links;
			for (int i = 0; i < links.size(); i++)
			{
				NWNode tmpNode = NWNode(Body[links[i]], new_links, Body, ID1, pos1, ID2, pos2, read.size(), MinMaxborders, NextNodes);
				if (!tmpNode.status) return;
				for (int j = 0; j < tmpNode.NWNext.size(); j++)
				{
					int t = tmpNode.NWNext[j];
					prevIds[t].push_back(tmpNode.ID);
				}

				NWBody[tmpNode.ID] = tmpNode;
			}
			sort(new_links.begin(), new_links.end());
			new_links.erase(unique(new_links.begin(), new_links.end()), new_links.end());
			links = new_links;
		}




		//Evaluate length
		links.clear();
		NWBody[ID1].MinMaxMark(NWBody, links);
		while (links.size() != 0)
		{
			vector<int> new_links;
			for (int i = 0; i < links.size(); i++)
			{
				NWBody[links[i]].MinMaxMark(NWBody, new_links);
			}
			sort(new_links.begin(), new_links.end());
			new_links.erase(unique(new_links.begin(), new_links.end()), new_links.end());
			links = new_links;
		}

		int Mult = 3;

		int Maxlen = NWBody[ID2].max_path;
		int Minlen = NWBody[ID2].min_path;

		int readlen = read.size();

		bool check1 = false;

		if (readlen*3+ MaxInDels <Maxlen)
		{
			check1 = true;
		}

		bool check2 = false;

		if (readlen * 2 + MaxInDels < Minlen)
		{
			check2 = true;
		}

		if (check1 or check2)
		{
			return;
		}

		//If graph is appropriate length, run NWInit
		//ATTENTION! Recursive B-graphs works incorectly, will be crashed at row 568 at .back() call
		links = { ID1 };
		vector<int> new_links = {};
		do
		{
			for (int i = 0; i < links.size(); i++)
			{
				int nId = links[i];
				bool ready2proc = true;

				for (int j = 0; j < prevIds[nId].size(); j++)
				{
					if (not NWBody[prevIds[nId][j]].NwMatrLoaded)
					{
						ready2proc = false;
						new_links.push_back(nId);
					}
				}

				if (ready2proc)
				{
					map<int, vector<int>> prevcolmns;
					for (int j = 0; j < prevIds[nId].size(); j++)
					{
						prevcolmns[prevIds[nId][j]] = NWBody[prevIds[nId][j]].NWMatrix.back();
					}
					if (prevcolmns.size() == 0)
					{
						NWBody[nId].InitMatrix(read);
					}
					else
					{
						NWBody[nId].InitMatrix(read, prevcolmns);
					}
					new_links.insert(new_links.end(), NWBody[nId].NWNext.begin(), NWBody[nId].NWNext.end());
				}
				
			}

			sort(new_links.begin(), new_links.end());
			new_links.erase(unique(new_links.begin(), new_links.end()), new_links.end());
			links = new_links;
			new_links.clear();

		} 
		while (links.size() != 0);
		

		nwID1 = ID1;
		nwID2 = ID2;
		nwpos1 = pos1;	
		nwpos2 = pos2;
		loaded = true;
	}

	NWpart NWTrace(string &read)
	{
		NWpart res = NWBody[lastID].NWTrace(read, NWpart());
		while (res.finished != true)
		{
			res = NWBody[res.ID].NWTrace(read, res);
		}
		reverse(res.vars.begin(), res.vars.end());
		res.ID1 = nwID1;
		res.ID2 = nwID2;
		res.pos1 = nwpos1;
		res.pos2 = nwpos2;
		return res;
	}

};