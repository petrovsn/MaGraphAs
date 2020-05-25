#pragma once
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

#include <numeric>


#include <ctime>
#include <cstdio>

#include <GenerHash.h>
#include <Fasta.h>
#include <Node.h>
#include <WArray.h>
#include<NWAligner.h>
#include<NWTailer.h>

using namespace std;

class FPoint
{
public:
	int score;
	int coverage;

	vector<int> seedIDs;
	vector<NWpart> NWlinks;

	FPoint()
	{
		score = 0;
		coverage = 0;
	}
};

class  FAlignment
{
public:
	int pos1 = -1;
	int pos2 = -1;

	int read_pos1 = -1;
	int read_pos2 = -1;

	vector<int> IDs;
	vector<GVariation> vars;

	int score = 0;
	int coverage = 0;
	float tresh = 0;

	//NWAligner nwaln;
	FAlignment() {};
	FAlignment(FPoint& fmap, vector<TArray>& StartTracks, NWpart head, NWpart tail, string read, map<int, Node>& Body)
	{
		pos1 = StartTracks[fmap.seedIDs[0]].startPos;
		pos2 = StartTracks[fmap.seedIDs.back()].endPos;

		read_pos1 = StartTracks[fmap.seedIDs[0]].read_startPos;
		read_pos2 = StartTracks[fmap.seedIDs.back()].read_endPos;

		for (int i = 0; i < fmap.seedIDs.size(); i++)
		{
			int idx = fmap.seedIDs[i];
			for (int j = 0; j < StartTracks[idx].NodeIDs.size(); j++)
			{
				IDs.push_back(StartTracks[idx].NodeIDs[j]);
			}

			coverage = coverage + StartTracks[idx].read_endPos - StartTracks[idx].read_startPos + 1;
		}

		sort(IDs.begin(), IDs.end());
		IDs.erase(unique(IDs.begin(), IDs.end()), IDs.end());


		vars.insert(vars.end(), head.vars.begin(), head.vars.end());
		for (int i = 0; i < fmap.NWlinks.size(); i++)
		{
			for (int j = 0; j < fmap.NWlinks[i].vars.size(); j++)
			{
				vars.push_back(fmap.NWlinks[i].vars[j]);
			}
		}
		vars.insert(vars.end(), tail.vars.begin(), tail.vars.end());
		

		score = fmap.score +tail.penalted_score + head.penalted_score;
		coverage = fmap.coverage +tail.score + head.score;
		tresh = coverage * 1.0 / read.size();
	}


};

class Graph
{
public:
	map<int, Node> Body;
	map<int, vector<int>> NextNodes;
	
	map<unsigned long long, vector<WArray>> Hashtable;
	
	int hashBase = -1;
	int hashLen = -1;


	Graph() {}

	//node operations

	void Link(Node* N1, Node* N2)
	{
		(N1->Next).insert(pair<int, Node*>(N2->ID, N2));
		(N2->Prev).insert(pair<int, Node*>(N1->ID, N1));
	}

	void LoadReference(string ref)
	{
		Node startNode(ref, 0);
		startNode.B.insert(pair<int, int>(0, -1));
		Body.insert(pair<int, Node>(0, startNode));

	//	Node endNode("x", -1);
	//	endNode.End = true;
	//	Body.insert(pair<int, Node>(-1, endNode));

	//	Link(&Body[0], &Body[-1]);

	}

	void LoadFromGFA(string file)
	{
		ifstream fin = ifstream(file);
		string header;
		do
		{
			getline(fin, header, '\t');
			if (header == "S")
			{
				string id;
				getline(fin, id, '\t');
				int ID = stoi(id);
				string seq;
				getline(fin, seq,'\t');
				Node node = Node(seq, ID);
				Body[ID] = node;
				string revseq = reverse(seq);
				Node revnode = Node(revseq, -ID);
				Body[-1 * ID] = revnode;
			}
			else if (header == "L")
			{
				string ID1buf;
				getline(fin, ID1buf, '\t');
				int ID1 = stoi(ID1buf);
				string sign1; 
				getline(fin, sign1, '\t');
				if (sign1 == "-")
				{
					ID1 = -1 * ID1;
				}

				string ID2buf;
				getline(fin, ID2buf, '\t');
				int ID2 = stoi(ID2buf);
				string sign2;
				getline(fin, sign2, '\t');
				if (sign2 == "-")
				{
					ID2 = -1 * ID2;
				}

				Link(&Body[ID1], &Body[ID2]);
			}
			else if (header == "P")
			{
				
			}
			else if (header == "i")
			{
			
			}
			else
			{
				//cout << "Smth wrong with gfa header\n";
			}
			getline(fin,header);
		} while (!fin.eof());
		fin.close();
	}

	//all for hash generator

	void BuildIndex(int base, int len)
	{
		Hashtable.clear();
		hashBase = base;
		hashLen = len;
		for (auto t1 : Body)
		{
			if ((!t1.second.NaNbp) && (t1.second.str.size() != 0))
			{
				GenHashIndexForNode(t1.second);
			}
			else
			{
				cout << "try to hash -1 or NaNbp Node\n";
			}
		}
		//BubbleIndexBuild();
		BfgIdxBuild(200);
	}

	void GenHashIndexForNode(Node n1)
	{

		GenerHash GH = GenerHash(hashBase, hashLen);
		GH.ReInit(n1.ID, n1.str);

		int res = 0;

		WArray tmpWarray = WArray(n1.ID, 0);

		for (int i = 0; i < n1.str.length(); i++)
		{
			if (GH.Next(n1.str[i]) == 0)
			{
				tmpWarray.finish(i);
				AddHash(tmpWarray, GH.currhash);


				tmpWarray.move_right();
			}
			else
			{
				//	cout<<"oooops"<<endl;
			}
		}

		for (auto node : n1.Next)
		{
			CallNode(*(node.second), GH, tmpWarray);
		}
	}

	void CallNode(Node n1, GenerHash GH2, WArray wray)
	{
		if (n1.NaNbp) return;
		unsigned long long res = 0;

		wray.extend(n1.ID);
		for (int i = 0; i < n1.str.length(); i++)
		{
			res = GH2.Next(n1.str[i]);
			if (res == -3)//pointer to position at the end of initial node. return;
			{
				return;
			}
			else
			{
				wray.finish(i);
				AddHash(wray, GH2.currhash);

				wray.move_right();

			}
		}

		for (auto p1 : n1.Next)
		{
			CallNode(*(p1.second), GH2, wray);
		}
	}

	void AddHash(WArray warray, unsigned long long hash)
	{
		Hashtable[hash].push_back(warray);
	}

	//B-index build

	void BubbleIndexBuild()
	{
		vector<int> markedNode;

		Body[0].bidx = BubbleIndex(true);
		markedNode.push_back(0);

		for (int i = 0; i < markedNode.size(); i++)
		{

			bool bubble = (Body[markedNode[i]].Next.size() > 1);
			int c = 0;
			for (auto n : Body[markedNode[i]].Next)
			{
				Body[n.first].GenBIdx(Body[markedNode[i]].bidx, c, bubble);
				markedNode.push_back(n.first);
				c++;
			}
		}
	}

	void BubbleIndexBuildTest()
	{
		int Bs = Body.size();
		Body[0].bidx = BubbleIndex(true);
		int c = 0;
		for (auto t : Body[0].Next)
		{
			t.second->GenBIdx(Body[0].bidx, c, true);
			c++;
		}

		for (int i = 2; i < Bs; i++)
		{
			bool bubble = (Body[i].Next.size() > 1);
			c = 0;
			for (auto n : Body[i].Next)
			{
				n.second->GenBIdx(Body[i].bidx, c, bubble);
				c++;
			}
		}
	}

	//Ancectors Idx

	void BfgIdxBuild(int len)
	{
		for (auto t1 : Body)
		{
			CallNode4Relate(t1.second, len);
		}
	}

	void CallNode4Relate(Node n1, int maxlen)
	{
		NextNodes[n1.ID] = vector<int>();
		for (auto p1 : n1.Next)
		{
			CallNode4Relate2(*p1.second, 0, maxlen, n1.ID);
		}
	}

	void CallNode4Relate2(Node &n1, int currlen, int maxlen, int RefID)
	{
		int left = maxlen - currlen;
		NextNodes[RefID].push_back(n1.ID);

		if (n1.str.size() < left)
		{
			for (auto p1 : n1.Next)
			{
				CallNode4Relate2(*p1.second, currlen+ n1.str.size(), maxlen, RefID);
			}
		}
	}


	//alignment&SNPcalling

	vector<TArray> GenerateSeeds(string& read)
	{
		GenerHash GH2 = GenerHash(hashBase, hashLen);
		GH2.ReInit(-2, read);

		vector<unsigned long long> num_read;
		for (int i = 0; i < read.size(); i++)
		{
			int res_hash = GH2.Next(read[i]);
			if (res_hash == 0)
			{
				num_read.push_back(GH2.currhash);
			}
		}

		vector<vector<WArray>> marks_read;
		for (int i = 0; i < num_read.size(); i++)
		{
			auto it = Hashtable.find(num_read[i]);
			if (it != Hashtable.end())
			{
				marks_read.push_back(Hashtable[num_read[i]]);
			}
			else
			{
				marks_read.push_back(vector<WArray>());
			}
			
		}



		vector<TArray> StartTracks;
		
		int cnt = 0;
		for (int i = 0; i < marks_read.size(); i++)
		{
			if (marks_read[i].size() == 0)
			{
				cnt = cnt + 1;
			}
		}

		if (cnt + hashLen > read.size() * 0.1)
		{
			//return StartTracks;
		}


		for (int i = 0; i < marks_read[0].size(); i++)
		{
			StartTracks.push_back(TArray(marks_read[0][i], 0, hashLen - 1));
		}

		vector<TArray> FinalStartTracks;

		for (int i = 1; i < marks_read.size(); i += 1)
		{
			for (int j = 0; j < StartTracks.size(); j++)
			{
				for (int k = 0; k < marks_read[i].size(); k++)
				{
					int res = StartTracks[j].TryExtRight(marks_read[i][k], Body, i, i + hashLen - 1);
					if (res == 0)
					{
						break;
					}
				}

				if (!StartTracks[j].extended)
				{
					StartTracks[j].futextended = false;
				}
			}

			for (int j = 0; j < StartTracks.size(); j++)
			{
				if (StartTracks[j].futextended == false)
				{
					FinalStartTracks.push_back(StartTracks[j]);
					StartTracks.erase(StartTracks.begin() + j);
					j = j - 1;
				}
			}

			for (int k = 0; k < marks_read[i].size(); k++)
			{
				if (marks_read[i][k].assembled == false)
				{
					StartTracks.push_back(TArray(marks_read[i][k], i, i + hashLen - 1));
				}
			}
		}

		StartTracks.insert(StartTracks.end(), FinalStartTracks.begin(), FinalStartTracks.end());



	/*	for (int i = 0; i < StartTracks.size(); i++)
		{
			if ((StartTracks[i].endPos >= Body[StartTracks[i].NodeIDs.back()].str.size()) or
				(StartTracks[i].startPos >= Body[StartTracks[i].NodeIDs[0]].str.size()))
			{
				cout << "seedwrong";
			}
		}*/




		return StartTracks;
	}

	void TrimSeeds(vector<TArray>& StartTracks)
	{
		for (int i = 0; i < StartTracks.size(); i++)
		{
			for (int j = 0; j < StartTracks.size(); j++)
			{
				if (i != j)
				{
					bool comp = StartTracks[i].IntersecSeed(StartTracks[j], NextNodes);
					if (comp)
					{
						StartTracks[j].CutFromLeft(StartTracks[i], Body);
					}
				}
			}
		}
	}

	FPoint FillLine(vector<vector<int>>& WeightMatrix, vector<vector<NWpart>>& NWWeightMatrix, vector<TArray>& StartTracks, int startline, string& read)
	{
		vector<int> starts = { startline };
		vector<FPoint> node_weights = vector<FPoint>(WeightMatrix.size());


		node_weights[startline].seedIDs.push_back(startline);
		node_weights[startline].score = StartTracks[startline].read_endPos - StartTracks[startline].read_startPos;
		node_weights[startline].coverage = node_weights[startline].score;

		NWAligner nwa;
		while (starts.size() != 0)
		{
			vector<int> new_starts;
			for (int idx = 0; idx < starts.size(); idx++)//for each of the edges
			{
				int i = starts[idx];
				for (int j = 0; j < WeightMatrix[i].size(); j++)
				{

					if (WeightMatrix[i][j] > 0)
					{
						
						NWpart nw_res = NWWeightMatrix[i][j];

						int RelScore = WeightMatrix[i][j] + nw_res.penalted_score + node_weights[i].score;

						if (RelScore >= node_weights[j].score)
						{
							node_weights[j].seedIDs = node_weights[i].seedIDs;
							node_weights[j].seedIDs.push_back(j);
							node_weights[j].score = RelScore;
							node_weights[j].coverage = node_weights[i].coverage + WeightMatrix[i][j] + nw_res.score;
							node_weights[j].NWlinks = node_weights[i].NWlinks;
							node_weights[j].NWlinks.push_back(nw_res);
							new_starts.push_back(j);
						}
					}
				}
			}
			starts = new_starts;
		}

		int fcov = 0;
		int ilast = -1;
		FPoint res;

		
		for (int i = 0; i < node_weights.size(); i++)
		{
			if (node_weights[i].coverage > fcov)
			{
				fcov = node_weights[i].coverage;
				res = node_weights[i];
				ilast = i;
			}
		}

		return res;
	}

	FAlignment AlignHashSmWtmn(string read)
	{

		int tmprs = read.size();
		if (tmprs < hashLen)
		{
			return FAlignment();
		}
		vector<TArray> StartTracks = GenerateSeeds(read); //call uninterrupted seeds

		//BAD CODE!!!==========================================


		if (StartTracks.size() > 2500)
		{
			return FAlignment();
		}
			


			//BAD CODE!!!==========================================

 		TrimSeeds(StartTracks); //trim seeds for indels with repeat;

		vector<vector<int>> WeightMatrix = vector<vector<int>>(StartTracks.size());

		for (int i = 0; i < WeightMatrix.size(); i++)
		{
			WeightMatrix[i] = vector<int>(StartTracks.size(), 0);
		}


		vector<vector<NWpart>> NWWeightMatrix = vector<vector<NWpart>>(StartTracks.size());
		for (int i = 0; i < NWWeightMatrix.size(); i++)
		{
			NWWeightMatrix[i] = vector<NWpart>(NWWeightMatrix.size(), NWpart());
		}


		NWAligner nwa;

		int StSize = StartTracks.size();

		//plot adjance matrix
		for (int i = 0; i < StSize; i++)
		{
			for (int j = 0; j < StSize; j++)
			{
				if (i != j)
				{
					if (StartTracks[i].read_endPos < StartTracks[j].read_startPos)
					{
						int ID1 = StartTracks[i].NodeIDs.back();
						int pos1 = StartTracks[i].endPos;
						int ID2 = StartTracks[j].NodeIDs[0];
						int pos2 = StartTracks[j].startPos;

						bool comp = Node1LessNode2(ID1, pos1, ID2, pos2, NextNodes);

						if (comp)
						{
							int w_plus = StartTracks[j].read_endPos - StartTracks[j].read_startPos + 1;
							WeightMatrix[i][j] = w_plus;
						}
					}
				}
			}
		}

		//filter excessive edges

		for (int i = 0; i < StSize; i++)
		{
			for (int j = 0; j < StSize; j++)
			{
				if (WeightMatrix[i][j] > 0)
				{
					for (int k = 0; k < StSize; k++)
					{
						if ((WeightMatrix[k][j] > 0) and (WeightMatrix[i][k] > 0))
						{
							WeightMatrix[i][j] = 0;
						}
					}
				}
			}
		}


		//Try build all NW matrix
		for (int i = 0; i < StSize; i++)
		{
			for (int j = 0; j < StSize; j++)
			{
				if ((i != j) and (WeightMatrix[i][j] > 0))
				{
					int ID1 = StartTracks[i].NodeIDs.back();
					int pos1 = StartTracks[i].endPos;
					int ID2 = StartTracks[j].NodeIDs[0];
					int pos2 = StartTracks[j].startPos;
					int rcoord = StartTracks[i].read_endPos;
					int rlen = StartTracks[j].read_startPos - StartTracks[i].read_endPos + 1;
					string nwread = read.substr(rcoord, rlen);


					nwa.Init(ID1, pos1, ID2, pos2, Body, nwread, NextNodes);

					NWpart nw_res;

					if (nwa.loaded)
					{
						nw_res = nwa.NWTrace(nwread);
						NWWeightMatrix[i][j] = nw_res;
					}
					else
					{
						WeightMatrix[i][j] = 0;
					}
				}
			}
		}
	




		vector<int> idxs_start; //all possible start points of D-Run;

		for (int j = 0; j < WeightMatrix.size(); j++)
		{
			int r = 0;
			for (int i = 0; i < WeightMatrix.size(); i++)
			{
				r = r + WeightMatrix[i][j];
			}
			if (r == 0)
			{
				idxs_start.push_back(j);
			}
		}

		

		

		vector<FPoint> maps;
		for (int i = 0; i < idxs_start.size(); i++)
		{
			FPoint res = FillLine(WeightMatrix, NWWeightMatrix, StartTracks, idxs_start[i], read);
			maps.push_back(res);
		}



		vector<FAlignment> total_aln;
		vector<FAlignment> filtered_aln;

		int maxscore = -1;
		int maxmapid = -1;
		for (int i = 0; i < maps.size(); i++)
		{
			if (maps[i].score > maxscore)
			{
				maxscore = maps[i].score;
				maxmapid = i;
			}
		}

		FAlignment tmp;
		NWTailer nwt;

		if (maxscore != -1)
		{
			int rlastpos = read.length() - 1;

			int seedhead = maps[maxmapid].seedIDs[0];

			//checkhead
			int backID = StartTracks[seedhead].NodeIDs[0];
			int backpos = StartTracks[seedhead].startPos;
			int rpos = StartTracks[seedhead].read_startPos;
			NWpart res_back;

			if ((rpos + 1) * 1.0 < tmprs * 0.25)
			{

				if ((rpos != 0) and (not((backpos < rpos) and (backID == 0))))
				{
					string tmpread = read.substr(0, rpos + 1);
					res_back = nwt.NWTback(backID, backpos, Body, tmpread);
				}
			}


			//checktail

			seedhead = maps[maxmapid].seedIDs.back();
			int forID = StartTracks[seedhead].NodeIDs.back();
			int forpos = StartTracks[seedhead].endPos;
			rpos = StartTracks[seedhead].read_endPos;
			NWpart res_forw;

			if ((read.size() - rpos) * 1.0 < tmprs * 0.25)
			{
				if (rpos != rlastpos)
				{
					string tmpread = read.substr(rpos);
					res_forw = nwt.NWTforw(forID, forpos, Body, tmpread);
				}
			}




			tmp = FAlignment(maps[maxmapid], StartTracks, res_back, res_forw, read, Body);
		}
		return tmp;
	}

	void GetCoverage(FAlignment aln)
	{
		if (aln.IDs.size() == 1)
		{
			int idx = aln.IDs[0];
			for (int i = aln.pos1; i <= aln.pos2; i++)
			{
				Body[idx].coverage[i] = Body[idx].coverage[i] + 1;
			}
		}
		else if (aln.IDs.size() == 2)
		{
			int IDfirst = aln.IDs[0];
			for (int i = aln.pos1; i < Body[aln.IDs[0]].coverage.size(); i++)
			{
				Body[IDfirst].coverage[i] = Body[IDfirst].coverage[i] + 1;
			}

			int IDlast = aln.IDs.back();
			for (int i = 0; i < aln.pos2; i++)
			{
				Body[IDlast].coverage[i] = Body[IDlast].coverage[i] + 1;
			}
		}
		else if (aln.IDs.size() > 2)
		{
			int IDfirst = aln.IDs[0];
			for (int i = aln.pos1; i < Body[aln.IDs[0]].coverage.size(); i++)
			{
				Body[IDfirst].coverage[i] = Body[IDfirst].coverage[i] + 1;
			}

			int IDlast = aln.IDs.back();
			for (int i = 0; i < aln.pos2; i++)
			{
				Body[IDlast].coverage[i] = Body[IDlast].coverage[i] + 1;
			}

			int sz = aln.IDs.size() - 1;
			for (int i = 1; i < sz; i++)
			{
				for (int j = 0; j < Body[aln.IDs[i]].coverage.size(); j++)
				{
					Body[aln.IDs[i]].coverage[j] = Body[aln.IDs[i]].coverage[j] + 1;
				}
			}
		}
	}

	//GraphTransformation

	map<int, vector<int>> Old2NewIds;

	void InitOld2NewIds()
	{
		for (auto n : Body)
		{
			Old2NewIds[n.first] = { n.first };
		}
	}

	pair<int, int> GetNewIDandPosition(int ID, int pos)//no optimal, too many accesses to Body and summarization;
	{
		int IDnew;
		int posnew;

		int curpos = 0;
		for (int i = 0; i < Old2NewIds[ID].size(); i++)
		{
			if ((curpos < pos) and (pos < curpos + Body[Old2NewIds[ID][i]].str.size()))
			{
				IDnew = Old2NewIds[ID][i];
				posnew = pos - curpos + 1;
				return pair<int, int>(IDnew, posnew);
			}
			curpos = curpos + Body[Old2NewIds[ID][i]].str.size();
		}

	}

	void InsertVariation(GVariation var)
	{
		int MaxCurrId = Body.size();
		if (var.ID1 != var.ID2) return;//fill after Big Tuesday

		pair<int, int> cord1 = GetNewIDandPosition(var.ID1, var.pos1);
		pair<int, int> cord2 = GetNewIDandPosition(var.ID2, var.pos2);

		if (cord1.first != cord2.first) return; //force check

		int ID = cord1.first;

		Node tail2 = Body[ID].SplitNoLink(cord2.second - 1, MaxCurrId);
		MaxCurrId++;
		Node tail1 = Body[ID].SplitNoLink(cord1.second, MaxCurrId);
		MaxCurrId++;
		Node alt = Node(var.alt, MaxCurrId);

		Body[tail2.ID] = tail2;
		Body[tail1.ID] = tail1;
		Body[alt.ID] = alt;

		//links to bubble
		Link(&Body[ID], &Body[tail1.ID]);
		Link(&Body[ID], &Body[alt.ID]);

		//links from bubble
		Link(&Body[tail1.ID], &Body[tail2.ID]);
		Link(&Body[alt.ID], &Body[tail2.ID]);

		//refill track from old nodeý
		vector<int>::iterator it = find(Old2NewIds[var.ID1].begin(), Old2NewIds[var.ID1].end(), ID) + 1;;
		vector<int> toinsert = { tail1.ID , tail2.ID };
		Old2NewIds[var.ID1].insert(it, toinsert.begin(), toinsert.end());
	}

	void UpgradeGraph(vector<GVariation> vars)
	{
		InitOld2NewIds();
		for (int i = 0; i < vars.size(); i++)
		{
			InsertVariation(vars[i]);
		}
	}

	bool CheckGVar(GVariation var, int counts, float tresh)
	{
		int cov1 = Body[var.ID1].coverage[var.pos1];
		int cov2 = Body[var.ID2].coverage[var.pos2];
		float medcov = (cov1 + cov2) * 1.0 / 2;

		if ((counts / medcov) > tresh) return true;
		return false;
	}

	vector<GVariation> FilterVarMap(unordered_map<GVariation, int, MyHashFunction>& VarMap)
	{
		vector<GVariation> res;
		vector<GVariation> todelete;
		for (auto v : VarMap)
		{
			bool check = CheckGVar(v.first, v.second, 0.8);
			if (!check)
			{
				todelete.push_back(v.first);
			}
			else
			{
				res.push_back(v.first);
			}
		}

		for (int i = 0; i < todelete.size(); i++)
		{
			VarMap.erase(todelete[i]);
		}
		return res;
	}

	void ReInitGraph()
	{
		for (auto n : Body)
		{
			Body[n.first].coverage = vector<int>(n.second.str.size(), 0);
			Body[n.first].bidx = BubbleIndex(1); //not correct, temporal constructor to distinguish in the code
		}
		Hashtable.clear();
	}

	//save and load
	void Save(string file)
	{
		ofstream fout = ofstream(file);
		for (auto n : Body)
		{
			fout << n.first << '\t';
			fout << n.second.str << '\t';
			for (auto n2 : n.second.Next)
			{
				fout << n2.second->ID << ',';
			}
			fout << endl;
		}
		fout.close();
	}
}; 
