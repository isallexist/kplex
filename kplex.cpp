#include<iostream>
#include<fstream>
#include <sstream>
#include<string>

#include<vector>
#include<set>
#include<stack>
#include<map>

#include<algorithm>
#include<math.h>

#include<omp.h>

#include "gurobi_c++.h"

using namespace std;

typedef unsigned int uint;			// uint is the same as size_t
typedef unsigned long ulong;

typedef vector<string> strv;
typedef set<string> strs;
typedef vector<int> intv;			
typedef set<int> ints;				// define the type of int vector for one adjacency list
typedef vector<ints> adj;			// adjacency lists

const int GRASPTL = 1800;			// Time limit for grasp
const int NIIL = 100;				// Limit of nonincreasing iteration
const double IPBCTL = 600;		// Time limit for IPBC, precisely it is the BC time limit
const double BCMINT = 30;			// Minimum run time for BC1
const double MIPTL = 3600;			// Time limit for simple mip

const int TABUSIZE = 3;
const double TABU_S_P = 0.5;

//#define KAVOIDS

#define OBJ GRB_DoubleAttr_ObjVal
#define TIM GRB_DoubleAttr_Runtime
#define DVV GRB_DoubleAttr_X 
#define THR GRB_IntParam_Threads
#define OUT GRB_IntParam_OutputFlag
#define STAT GRB_IntAttr_Status
#define MAX GRB_MAXIMIZE 

template <typename T>
string num2str ( T Number )
{
	stringstream ss;
	ss << Number;
	return ss.str();
}
void splitstr(strv & splits, const string& str, const string& deli = " ")
{	// Note: splits are never cleared inside!!!
	size_t start, end;
	
	start = str.find_first_not_of(deli);
	while(start!=string::npos)
	{
		end = str.find_first_of(deli,start+1);

		if(end!=string::npos)
		{
			splits.push_back(str.substr(start,end-start));
			start = str.find_first_not_of(deli,end+1);
		}
		else 
		{
			splits.push_back(str.substr(start));	
			break;
		}
	}
}
GRBVar* GRBVarArray(int n, GRBModel &model, double lb = 0, double ub = GRB_INFINITY, 
     char type = GRB_CONTINUOUS, double obj_coeff = 1) 
{
	 // delete after use
	 GRBVar* vars = new GRBVar[n];

	 for(int i=0;i<n;i++)
		vars[i] = model.addVar(lb, ub, obj_coeff, type);

	 return vars;
}
GRBLinExpr GRBSum(int n, GRBVar* vars, double* coeff = NULL)
{
	int i;
	GRBLinExpr expr = 0;

	if(coeff!=NULL)
	{
		for(i=0;i<n;i++)
		expr += (vars[i]*coeff[i]);  
	}
	else
	{
		for(i=0;i<n;i++)
		expr += vars[i];
	}
	return expr;
}
//////////////////////////////////////////////////////////////
// Classes
class graph
{	//	Caution: if the graph size is above 2,147,483,647;
	//	the reader will be wrong
private:
	// type: 
	// 0:	edge list, delimited using tab; 
	// 1:	edge list, delmited using space;
	// 2:	adj
	// 3:	dimacs 2
	// 4:	dimacs 10
	int thre, sense, type;					
	map<string,strs> adjlist;
	map<string,int> ms2i;

	//-------------------------------------------------------------//
	void load_edgelist_to_stradj(const string &file)
	{
		string str;
		ifstream fs(file.c_str());

		if(!fs)
		{
			cout<<"Cannot open!"<<endl;
			return;
		}

		while(!fs.eof())
		{
			getline(fs,str);

			if(str.size()<1 || str[0]=='#')
				continue;

			strv sv;
			if(type==0)
				splitstr(sv,str,"	");
			else
				splitstr(sv,str);

			if(sv.size()>2)
			{
				int wt = atoi(sv[2].c_str());
				if(wt*sense >= thre*sense)
				{
					string a = sv[0];
					string b = sv[1];
					if(a!=b)
					{
						adjlist[a].insert(b);
						adjlist[b].insert(a);
					}
				}
			}
			else
			{
				string a = sv[0];
				string b = sv[1];
				if(a!=b)
				{
					adjlist[a].insert(b);
					adjlist[b].insert(a);
				}
			}
		}

		int vi = 0;	ecnt = 0;
		map<string,strs>::iterator i;
		for(i=adjlist.begin();i!=adjlist.end();i++)
		{
			ecnt += (int)i->second.size();
			ms2i[i->first] = vi;
			mi2s.push_back(i->first);
			vi++;
		}

		vcnt = (int)adjlist.size();
		ecnt /= 2;
		fs.close();
	}
	void load_stradj(const string &file)
	{
		string str;
		strv::iterator b;
		ifstream fs(file.c_str());

		if(!fs)
		{
			cout<<"Cannot open!"<<endl;
			return;
		}

		while(!fs.eof())
		{
			getline(fs,str);

			if(str.size()<1)
				continue;

			strv sv;			
			splitstr(sv,str);
	
			string a = sv[0];			
			for(b=sv.begin()+1;b!=sv.end();b++)
			{
				if(a!=*b)
				{
					adjlist[a].insert(*b);
					adjlist[*b].insert(a);
				}
			}
		}

		int vi = 0;	ecnt = 0;
		map<string,strs>::iterator i;
		for(i=adjlist.begin();i!=adjlist.end();i++)
		{
			ecnt += (int)i->second.size();
			ms2i[i->first] = vi;
			mi2s.push_back(i->first);
			vi++;
		}

		vcnt = (int)adjlist.size();
		ecnt /= 2;
		fs.close();
	}
	void load_dimacs2(const string &file)
	{
		string str;
		ifstream fs(file.c_str());

		if(!fs)
		{
			cout<<"Cannot open!"<<endl;
			return;
		}

		while(!fs.eof())
		{
			getline(fs,str);

			if(str.size()<1 || str[0]!='e')
				continue;
			
			strv sv;
			splitstr(sv,str);
			
			string a = sv[1];
			string b = sv[2];
			if(a!=b)
			{
				adjlist[a].insert(b);
				adjlist[b].insert(a);
			}
		}

		int vi = 0;	ecnt = 0;
		map<string,strs>::iterator i;
		for(i=adjlist.begin();i!=adjlist.end();i++)
		{
			ecnt += (int)i->second.size();
			ms2i[i->first] = vi;
			mi2s.push_back(i->first);
			vi++;
		}

		vcnt = (int)adjlist.size();
		ecnt /= 2;
		fs.close();
	}
	void load_dimacs10(const string &file)
	{
		string str;
		strv::iterator b;
		ifstream fs(file.c_str());

		if(!fs)
		{
			cout<<"Cannot open!"<<endl;
			return;
		}
		
		strv parv;
		getline(fs,str);	//skip the first line		
		splitstr(parv,str);
		int vtx = 0;		//vid will be 1
		bool hasweight = false;
		if(parv.size()>2 && parv[2]=="1"){
			hasweight = true;
		}

		while(!fs.eof())
		{
			getline(fs,str);

			if(str.size()<1 || str[0]=='%')
				continue;

			vtx++;

			strv sv;			
			splitstr(sv,str);
	
			string a = num2str(vtx);			
			for(b=sv.begin();b!=sv.end();)
			{
				if(a!=*b)
				{
					adjlist[a].insert(*b);
					adjlist[*b].insert(a);
				}
				
				if(hasweight){
					b += 2;
				}
				else{
					b++;
				}
			}
		}

		int vi = 0;	ecnt = 0;
		map<string,strs>::iterator i;
		for(i=adjlist.begin();i!=adjlist.end();i++)
		{
			ecnt += (int)i->second.size();
			ms2i[i->first] = vi;
			mi2s.push_back(i->first);
			vi++;
		}

		vcnt = (int)adjlist.size();
		ecnt /= 2;
		fs.close();
	}
	/////////////////////////////////////////////////////////////////
	void convert_to_simpleadj()
	{
		//indexed from 0 to n-1
		AdjList.resize(vcnt);
		map<string,strs>::iterator i;
		strs::iterator j;
		int a, b;
		for(i=adjlist.begin();i!=adjlist.end();i++)
		{
			for(j=i->second.begin();j!=i->second.end();j++)
			{
				a = ms2i[i->first];
				b = ms2i[*j];
				AdjList[a].insert(b);
			}
		}
	}
public:
	int vcnt, ecnt;
	adj AdjList;
	strv mi2s;

	graph(const string &file, int _type, int _sense=1, int _thre=-5) 
	{	
		//seems -5 can include every edge
		thre = _thre;
		sense = _sense;
		type = _type;
		if(type<2)
		{
			load_edgelist_to_stradj(file);
		}
		else if(type==2)
		{
			load_stradj(file);
		}
		else if(type==3)
		{
			load_dimacs2(file);
		}
		else
		{
			load_dimacs10(file);
		}

		convert_to_simpleadj();
	}
	~graph(){}

	bool IsAdjTo(int a, int b)
	{
		return (AdjList[a].find(b)!=AdjList[a].end());
	}	
};

class grasp
{
private:

	graph& g;
	adj &AdjList;

	int n, k, unex, tl;		// tl: the time limit; unex: unexplored

	intv DegCand;			// degrees of nodes in the subgraph induced by candidates
	intv nncnt;				// keeps count of NonNeighboutrs 'nn'
	intv satu;				// saturate list
	intv stat;				// 1 added; 0 unexplored; 2 can not be added; 3 is completely deleted	
	vector<ints> tabu;		// old solutions (not the best solution) are tabued

	void DegUpdate(int del_v)
	{
		ints::iterator i;
		for(i=AdjList[del_v].begin();i!=AdjList[del_v].end();i++)
			DegCand[*i]--;
	}
	void ResetTrack(intv &sol)
	{
		sol.clear();	
		satu.clear();
		nncnt.assign(n,0);

		unex = glo_unex; 
		stat = glo_stat;		
		DegCand = glo_deg;
	}
	void Kick()
	{
		int i;
		bool vio;
		intv::iterator j;	
		
		for(i=0;i<n;i++)
		{
			vio = false;
			if(stat[i]==0)
			{
				if(nncnt[i] > (k-1))  
					vio = true;
				else
				{
					for (j=satu.begin(); j!=satu.end(); j++)
					{
						if (!g.IsAdjTo(i,*j))
						{
							vio = true;
							break;							
						}
					}
				}
			}
			if(vio)
			{
				stat[i] = 2;
				unex--;
				DegUpdate(i);
			}
		}
	}

	//------------------------------------------------------------------------//
	// Pre-process, update unex, stat and DegCand !!!!
	void peel(int v, stack<int>& peel_stack)
	{		
		int top;
		ints::iterator i;

		glo_stat[v] = 3;
		glo_unex--;
		peel_stack.push(v);

		while(!peel_stack.empty())
		{
			top = peel_stack.top();
			peel_stack.pop();		
			
			for(i=AdjList[top].begin();i!=AdjList[top].end();i++)
			{
				glo_deg[*i]--;
				if( (glo_stat[*i]!=3) && (glo_deg[*i] <= (int)bestsol.size() - k) )
				{
					glo_stat[*i] = 3;
					glo_unex--;
					peel_stack.push(*i);
				}
			}
		}	
	}
	void PeelProc()
	{
		stack<int> peel_stack;

		for(int i=0;i<n;i++)
		{
			if( (glo_stat[i]!=3) && (glo_deg[i] <= (int)bestsol.size() - k) )
				peel(i,peel_stack);
		}
	}

	// used when k+1 solution is smaller than k solution
	void peel(int v, stack<int>& peel_stack, int spec_k)
	{		
		int top;
		ints::iterator i;

		glo_stat[v] = 3;
		glo_unex--;
		peel_stack.push(v);

		while(!peel_stack.empty())
		{
			top = peel_stack.top();
			peel_stack.pop();		
			
			for(i=AdjList[top].begin();i!=AdjList[top].end();i++)
			{
				glo_deg[*i]--;
				if( (glo_stat[*i]!=3) && (glo_deg[*i] <= (int)bestsol.size() - spec_k) )
				{
					glo_stat[*i] = 3;
					glo_unex--;
					peel_stack.push(*i);
				}
			}
		}	
	}
	void PeelProc(int spec_k)
	{
		stack<int> peel_stack;

		for(int i=0;i<n;i++)
		{
			if( (glo_stat[i]!=3) && (glo_deg[i] <= (int)bestsol.size() - spec_k) )
				peel(i,peel_stack,spec_k);
		}
	}
	////////////////////////////////////////////////////////////////////////////


	//-----------------------------------------------//
	//	Construct 
	bool in_tabu(int j)
	{
		bool rtn = false;
		vector<ints>::iterator i;
		for(i=tabu.begin();i!=tabu.end();i++)
		{
			if( i->find(j) != i->end() )
			{
				rtn = true;
				break;
			}
		}
		return rtn;
	}
	int HybSelect(double alpha, bool avoid_tabu=false)
	{	// avoid tabu for the first few vertices
		intv RCL;
		double deg_thre;
		int i, max = -1, min = n, nn_thre = n;; 

		if(avoid_tabu)
		{
			for (i=0;i<n;i++)
			{
				if( (stat[i] == 0) && (nncnt[i] < nn_thre) && (!in_tabu(i)) )
					nn_thre = nncnt[i];
			}

			for(i=0;i<n;i++)
			{
				if( (stat[i] == 0) && (nncnt[i] <= nn_thre) && (!in_tabu(i)) )
				{
					if(DegCand[i] > max)	
						max = DegCand[i];
						
 					if(DegCand[i] < min)
						min = DegCand[i];
				}
			}
		
			deg_thre = min + alpha*(max-min);

			for(i=0;i<n;i++) 
			{
				if ( (stat[i] == 0) && (nncnt[i] <= nn_thre) && (DegCand[i]>=deg_thre) && (!in_tabu(i)) ) 
					RCL.push_back(i);
			}
		}
		else
		{
			for (i=0;i<n;i++)
			{
				if( (stat[i] == 0) && (nncnt[i] < nn_thre) )
					nn_thre = nncnt[i];
			}

			for(i=0;i<n;i++)
			{
				if( (stat[i] == 0) && (nncnt[i] <= nn_thre) )
				{
					if(DegCand[i] > max)	
						max = DegCand[i];
						
 					if(DegCand[i] < min)
						min = DegCand[i];
				}
			}
		
			deg_thre = min + alpha*(max-min);

			for(i=0;i<n;i++) 
			{
				if ( (stat[i] == 0) && (nncnt[i] <= nn_thre) && (DegCand[i]>=deg_thre) ) 
					RCL.push_back(i);
			}
		}
		
		if(RCL.empty()) 
			return -1;
		else
			return RCL[(int)floor( rand() * RCL.size() / (RAND_MAX + 1.0) )]; //simplified
	}
	void NoTabuHybCons(double alpha, intv &sol)
	{
		int i, chosen; 

		while(unex!=0)
		{
			chosen = HybSelect(alpha);

			if(chosen==-1) break;

			stat[chosen] = 1;
			unex--;
			sol.push_back(chosen);
			DegUpdate(chosen);
			
			////////////////////////////////////////////////
			// Update nncnt and satu
			for(i=0;i<n;i++)
			{	
				if(stat[i]==3) continue;	// skip preprocessed nodes

				if( (i!=chosen) && (!g.IsAdjTo(chosen,i)) )
				{
					nncnt[i] ++;
					if( (stat[i]==1) && nncnt[i] == k - 1 )
						satu.push_back(i); // May be good to put outside
				}
			}
			if( nncnt[chosen] == k-1)
				satu.push_back(chosen);		
			////////////////////////////////////////////////		
			
			Kick();
		}
	}
	void HybCons(double alpha, intv &sol)
	{
		int i, chosen, select_cnt = 0; 

		while(unex!=0)
		{
			select_cnt++;

			// set when to avoid
		#ifdef KAVOIDS
			if( select_cnt<=k )
		#else
			if( select_cnt<=1 )
		#endif			
				chosen = HybSelect(alpha,true);
			else
				chosen = HybSelect(alpha);

			if(chosen==-1) break;

			stat[chosen] = 1;
			unex--;
			sol.push_back(chosen);
			DegUpdate(chosen);
			
			////////////////////////////////////////////////
			// Update nncnt and satu
			for(i=0;i<n;i++)
			{	
				if(stat[i]==3) continue;	// skip preprocessed nodes

				if( (i!=chosen) && (!g.IsAdjTo(chosen,i)) )
				{
					nncnt[i] ++;
					if( (stat[i]==1) && nncnt[i] == k - 1 )
						satu.push_back(i); // May be good to put outside
				}
			}
			if( nncnt[chosen] == k-1)
				satu.push_back(chosen);		
			////////////////////////////////////////////////		
			
			Kick();
		}
	}
	void HybCons(double alpha, bool in_rp_niicnt, intv &sol)
	{
		int i, chosen, select_cnt = 0; 

		while(unex!=0)
		{
			select_cnt++;

			// set when to avoid
		#ifdef KAVOIDS
			if( in_rp_niicnt && select_cnt<=k )
		#else
			if( in_rp_niicnt && select_cnt<=1 )
		#endif			
				chosen = HybSelect(alpha,true);
			else
				chosen = HybSelect(alpha);

			if(chosen==-1) break;

			stat[chosen] = 1;
			unex--;
			sol.push_back(chosen);
			DegUpdate(chosen);
			
			////////////////////////////////////////////////
			// Update nncnt and satu
			for(i=0;i<n;i++)
			{	
				if(stat[i]==3) continue;	// skip preprocessed nodes

				if( (i!=chosen) && (!g.IsAdjTo(chosen,i)) )
				{
					nncnt[i] ++;
					if( (stat[i]==1) && nncnt[i] == k - 1 )
						satu.push_back(i); // May be good to put outside
				}
			}
			if( nncnt[chosen] == k-1)
				satu.push_back(chosen);		
			////////////////////////////////////////////////		
			
			Kick();
		}
	}
	///////////////////////////////////////////////////


	//-----------------------------------------------//
	//	Local Search
	void LocS(intv &sol)
	{
		intv Newnncnt;
		intv addable;

		int i = 0, j, u, w; 
		intv::iterator l;

		bool vio;

		while(i<(int)sol.size())
		{
			w = sol[i];
			i++;

			Newnncnt = nncnt;

			for(j=0;j<n;j++)
			{
				if(stat[i]==3) continue;	//skip the preprocessed nodes

				if ( !g.IsAdjTo(w,j) && (j != w) )
					Newnncnt[j] --;	
			}

			//////////////////////////////////////////////////////////////////
			// Find addable
			for (u=0; u<n; u++)
			{
				if(stat[u]==1 || stat[u]==3) continue;

				vio = false;

				if(Newnncnt[u] > k-1 )	vio=true;
				else
				{
					for(l=sol.begin();l!=sol.end();l++)
					{
						j=*l;
						if((Newnncnt[j]>=k-1)&&(j!=w)&&(!g.IsAdjTo(u,j)))
						{
							vio=true;
							break;
						}
					}
				}
				if(!vio)
				{
					for(l=addable.begin();l!=addable.end();l++)
					{
						if( (Newnncnt[*l] >= k-1) && (!g.IsAdjTo(u,*l)) )
						{
							vio = true;
							break;
						}
					}
				}
				if(!vio)
				{
					for(j=0;j<n;j++)
					{
						if( (j != u) && (!g.IsAdjTo(u,j)) )
							Newnncnt[j]++;
					}					
					addable.push_back(u);				
				}
			}//end for 'u'
			////////////////////////////////////////////////////////////////////
			
			if(addable.size()>1)
			{
				sol.erase(remove(sol.begin(),sol.end(),w),sol.end());
				stat[w]=2;
				
				while(!addable.empty())	
				{
					sol.push_back(addable.back());
					stat[addable.back()]=1;
					addable.pop_back();
				}

				nncnt=Newnncnt;			
				i=0;
			}
			else addable.clear();
		}
	}
	///////////////////////////////////////////////////

public:
	grasp(graph& g_obj,int k_value,int time):g(g_obj),AdjList(g.AdjList),glo_deg(g.vcnt)
	{
		k = k_value;
		n = g.vcnt;
		tl = time;

		//----------------------------------//
		// Initialize global
		glo_unex = n;
		glo_stat.assign(n,0);

		for (int i=0;i<n;i++){ 
			glo_deg[i] = (int)AdjList[i].size();
		}
		/////////////////////////////////////

		//-----------------------------------------//
		// Initialize part of statistics other part
		// is let go, since must be assigned later
		ct = 0;
		lt = 0;

		lirate = 0;
		licnt = 0;

		iter_cnt = 0;

		pret = 0;
		pre_cnt = 0;

		opt_proof = false;
		//////////////////////////////////////////////	

		//srand(19840729);
	}
	~grasp(){}

	//---------------------------------------------------------------------//
	//	Statistics
	intv bestsol;		// storing the best sol vertices
	size_t fs;			// 1st solution size

	double ft, bt;		// 1st and best time
	double ct, lt;		// construction & local search time
	double tt;			// total time

	double lirate;		// LS increase rate for solutions
	ulong licnt;		// LS increase rate for # of times	
	
	ulong iter_cnt;		// Total times of iterations

	double pret;		// Preprocess time
	int pre_cnt;		// # of times that preprocess is performed
						// the same as larger solutions are found

	intv update_iter;	// updated iteration

	bool opt_proof;
	/////////////////////////////////////////////////////////////////////////

	//------------------------------------------//
	// For peeling, also the stat passed to ipbc
	int glo_unex;
	intv glo_stat, glo_deg;
	//////////////////////////////////////////////
	
	//----------------------------------------------//
	// use tabu
	void Run(int most_recent)
	{
		double st = omp_get_wtime();

		intv sol;
		bool IsFirst = true;

		size_t cs;

		double t; // temporary timer to record time for construct and local search
		
		int uninc_iter = 0;

		while(true)
		{
			iter_cnt++;
			uninc_iter++;

			///////////////////////////////////////////////
			// Generate random factors 
			double alpha = rand() / (RAND_MAX + 1.0);
			///////////////////////////////////////////////

			ResetTrack(sol);

			//---------------------------//
			// Construction phase
			t = omp_get_wtime();
			HybCons(alpha,sol);
			ct += ( omp_get_wtime() - t );
			///////////////////////////////
			
			cs = sol.size();

			//--------------------//
			// Local search phase
			t = omp_get_wtime();
			LocS(sol);
			lt += ( omp_get_wtime() - t );			
			////////////////////////

			//--------------------------------------------------------//
			// tabu most recent most_recent solutions
			if( tabu.size() < most_recent )
			{
				ints temp;
				for(intv::iterator i=sol.begin();i!=sol.end();i++){
					temp.insert(*i);
				}	
				tabu.push_back(temp);
			}
			else
			{
				tabu[0].clear();
				for(intv::iterator i=sol.begin();i!=sol.end();i++){
					tabu[0].insert(*i);
				}	
			}
			////////////////////////////////////////////////////////////

			if(sol.size()>cs)
			{
				licnt++;
				lirate += (sol.size() - (double)cs)/cs;
			}
			
			if ( sol.size() > bestsol.size() )
			{
				bt = omp_get_wtime() - st;
				bestsol = sol;

				update_iter.push_back(uninc_iter);
				uninc_iter = 0;

				//-------------------------------------//
				// Preprocess to peel
				t = omp_get_wtime();
				PeelProc();	
				pret += ( omp_get_wtime() - t );

				pre_cnt++;
				/////////////////////////////////////////
			}

			if( IsFirst )
			{
				ft = bt;
				fs = sol.size();
				IsFirst = false;
			}	

			if (glo_unex<=(int)bestsol.size())
			{
				opt_proof = true;
				break;
			}
			
			if( (uninc_iter >= NIIL) || (omp_get_wtime() >= st + tl) ) break;
		}
		
		tt = omp_get_wtime() - st;

		if (licnt > 0) 
			lirate = lirate / licnt;		
	}
	void Run(double niil_portion)
	{		
		double st = omp_get_wtime();

		intv sol;
		bool IsFirst = true;

		size_t cs;

		double t; // temporary timer to record time for construct and local search
		
		int uninc_iter = 0;

		while(true)
		{
			iter_cnt++;
			uninc_iter++;

			int niil_thre = (int)floor(NIIL*niil_portion);

			///////////////////////////////////////////////
			// Generate random factors 
			double alpha = rand() / (RAND_MAX + 1.0);
			///////////////////////////////////////////////

			ResetTrack(sol);

			//---------------------------//
			// Construction phase
			t = omp_get_wtime();
			HybCons(alpha, uninc_iter>niil_thre, sol);
			ct += ( omp_get_wtime() - t );
			///////////////////////////////
			
			cs = sol.size();

			//--------------------//
			// Local search phase
			t = omp_get_wtime();
			LocS(sol);
			lt += ( omp_get_wtime() - t );			
			////////////////////////

			//--------------------------------------------------------//
			// start tabu after ...
			if( uninc_iter >= niil_thre )
			{
				if(tabu.size()==0){
					tabu.resize(1);
				}
				for(intv::iterator i=sol.begin();i!=sol.end();i++){
					tabu[0].insert(*i);
				}	
			}

			////////////////////////////////////////////////////////////

			if(sol.size()>cs)
			{
				licnt++;
				lirate += (sol.size() - (double)cs)/cs;
			}
			
			if ( sol.size() > bestsol.size() )
			{
				bt = omp_get_wtime() - st;
				bestsol = sol;

				update_iter.push_back(uninc_iter);
				uninc_iter = 0;
				tabu.clear();

				//-------------------------------------//
				// Preprocess to peel
				t = omp_get_wtime();
				PeelProc();	
				pret += ( omp_get_wtime() - t );

				pre_cnt++;
				/////////////////////////////////////////
			}

			if( IsFirst )
			{
				ft = bt;
				fs = sol.size();
				IsFirst = false;
			}	

			if (glo_unex<=(int)bestsol.size())
			{
				opt_proof = true;
				break;
			}
			
			if( (uninc_iter >= NIIL) || (omp_get_wtime() >= st + tl) ) break;
		}
		
		tt = omp_get_wtime() - st;

		if (licnt > 0) 
			lirate = lirate / licnt;
	}
	//////////////////////////////////////////////////

	//-----------------------------------------//
	// no tabu
	void Run()
	{
		double st = omp_get_wtime();

		intv sol;
		bool IsFirst = true;

		size_t cs;

		double t; // temporary timer to record time for construct and local search
		
		int uninc_iter = 0;

		while(true)
		{
			iter_cnt++;
			uninc_iter++;

			///////////////////////////////////////////////
			// Generate random factors 
			double alpha = rand() / (RAND_MAX + 1.0);
			///////////////////////////////////////////////

			ResetTrack(sol);

			//---------------------------//
			// Construction phase
			t = omp_get_wtime();
			NoTabuHybCons(alpha,sol);
			ct += ( omp_get_wtime() - t );
			///////////////////////////////
			
			cs = sol.size();

			//--------------------//
			// Local search phase
			t = omp_get_wtime();
			LocS(sol);
			lt += ( omp_get_wtime() - t );			
			////////////////////////

			if(sol.size()>cs)
			{
				licnt++;
				lirate += (sol.size() - (double)cs)/cs;
			}
			
			if ( sol.size() > bestsol.size() )
			{
				bt = omp_get_wtime() - st;
				bestsol = sol;

				update_iter.push_back(uninc_iter);
				uninc_iter = 0;

				////-------------------------------------//
				//// Preprocess to peel
				//t = omp_get_wtime();
				//PeelProc();	
				//pret += ( omp_get_wtime() - t );

				//pre_cnt++;
				///////////////////////////////////////////
			}

			if( IsFirst )
			{
				ft = bt;
				fs = sol.size();
				IsFirst = false;
			}	

			if (glo_unex<=(int)bestsol.size())
			{
				opt_proof = true;
				break;
			}
			
			if( (uninc_iter >= NIIL) || (omp_get_wtime() >= st + tl) ) break;
		}
		
		tt = omp_get_wtime() - st;
		
		if (licnt > 0) 
			lirate = lirate / licnt;
	}
	/////////////////////////////////////////////

	void ksol_for_kplus1sol()
	{
		//----------------------------------//
		// reset global
		glo_unex = n;
		glo_stat.assign(n,0);

		for (int i=0;i<n;i++){ 
			glo_deg[i] = (int)AdjList[i].size();
		}

		opt_proof = false;
		/////////////////////////////////////

		PeelProc(k+1);	

		if (glo_unex<=(int)bestsol.size()){
			opt_proof = true;
		}
	}

	void load_peel(int kplx_size)
	{	// load from text and do preprocessing
		// call without run()

		bestsol.assign(kplx_size,0); // Caution: all 0 for convenience !!!!!!!!!

		//----------------------------------//
		// reset global
		glo_unex = n;
		glo_stat.assign(n,0);

		for (int i=0;i<n;i++){ 
			glo_deg[i] = (int)AdjList[i].size();
		}

		opt_proof = false;
		/////////////////////////////////////

		PeelProc();	

		if (glo_unex<=(int)bestsol.size()){
			opt_proof = true;
		}
	}
};
struct deg_greater
{
	const intv & degs;
	deg_greater(const intv &degrees):degs(degrees){}

	bool operator()(int a, int b)
	{
		return degs[a] > degs[b];
	}
};

class ipbc
{
private:
	graph& g;
	adj &AdjList;

	int n, k, BC_size, remain_cnt;				// BC_size: the size of the graph for BC

	intv ddv, degs;								// ddv: degree decreasing vertecies

	vector<bool> deleted, BC_Cand;				// track
	intv N2;

	//------------------------------------------------------------------------------------//
	// N2 related functions
	void FindN2(int v)
	{
		ints::iterator i, j;

		BC_Cand.assign(n,false);
		N2.clear();

		BC_Cand[v] = true;
		N2.push_back(v);

		for(i=AdjList[v].begin();i!=AdjList[v].end();i++)
		{
			if(!deleted[*i])
			{
				BC_Cand[*i] = true;
				N2.push_back(*i);
			}
		}

		for(i=AdjList[v].begin();i!=AdjList[v].end();i++)
		{
			if(deleted[*i]) continue;
			
			for(j=AdjList[*i].begin();j!=AdjList[*i].end();j++)
			{
				if( (!deleted[*j]) && (!BC_Cand[*j]) )
				{
					BC_Cand[*j] = true;
					N2.push_back(*j);			
				}
			}
		}
	}

	void N2DegUpdate(int v, intv& BC_Cand_degs)		
	{	// Update degrees BC_Cand_degs in accordence with N2
		// it will be later updated again in peeling
		ints::iterator i;
		for(i=AdjList[v].begin();i!=AdjList[v].end();i++)
			BC_Cand_degs[*i]--;
	}

	void peelN2(int v, intv& BC_Cand_degs, stack<int>& peel_stack)
	{
		int top;
		ints::iterator i;

		BC_Cand[v] = false;
		BC_size --;
		peel_stack.push(v);

		while(!peel_stack.empty())
		{
			top = peel_stack.top();
			peel_stack.pop();		
			
			for(i=AdjList[top].begin();i!=AdjList[top].end();i++)
			{
				BC_Cand_degs[*i]--;
				if( BC_Cand[*i] && BC_Cand_degs[*i]<= (int)sol.size() - k )
				{
					BC_Cand[*i] = false;
					BC_size --;
					peel_stack.push(*i);
				}
			}
		}	
	}

	void PeelN2Proc()
	{
		int i;
		intv::iterator j;	
		stack<int> peel_stack;
	
		intv BC_Cand_degs = degs;
		BC_size = (int)N2.size();

		for(i=0;i<n;i++)
		{
			if(!BC_Cand[i]&&!deleted[i])
				N2DegUpdate(i, BC_Cand_degs);
		}	
	
		for(j=N2.begin();j!=N2.end();j++)
		{
			if( BC_Cand[*j] && BC_Cand_degs[*j]<= (int)sol.size() - k )
				peelN2(*j, BC_Cand_degs, peel_stack);
		}
	}
	////////////////////////////////////////////////////////////////////////////////////////

	void BC1(int v_equals_1,intv& BC1_sol, GRBEnv &env, double time_limit)		// BC for peeled graph
	{
		intv vertices;
		intv::iterator l;
		int i, d, j, cnt, x_equals_1;

		BC1_sol.clear();	

		for(l=N2.begin();l!=N2.end();l++)
		{
			if(BC_Cand[*l])
			{
				if( (*l) == v_equals_1 )
					x_equals_1 = (int)vertices.size();

				vertices.push_back(*l);
			}
		}

		cnt = (int)vertices.size();

		try
		{
			cout<<"building the model..."<<endl;
			env.set(GRB_DoubleParam_Cutoff, (double)sol.size());
			env.set(GRB_DoubleParam_TimeLimit, time_limit);
			GRBModel model(env);
		
			GRBVar* x = GRBVarArray(cnt, model, 0.0, 1.0, GRB_BINARY);
			model.update();

			model.setObjective(GRBSum(cnt, x), MAX);

			model.addConstr(x[x_equals_1] == 1);
			model.addConstr(GRBSum(cnt, x) >= (int)sol.size()+1);

			for(i=0; i<cnt; i++) 
			{	
				d = 0;
				GRBLinExpr v = 0;
				for(j=0; j<cnt; j++) 
				{
					if (( !g.IsAdjTo(vertices[i],vertices[j]) )&&(i!=j))
					{
						d += 1;
						v += x[j];
					}					
				}

				v -= ((k-1)*x[i]) + (d*(1-x[i]));			
				model.addConstr(v <= 0);
			}

			cout<<"optimizing..."<<endl;
			model.optimize();

			int opt_status = model.get(STAT);
			if( opt_status==GRB_OPTIMAL || opt_status==GRB_TIME_LIMIT )
			{	
				double temp_ub = model.get(GRB_DoubleAttr_ObjBound);
				temp_ub = min(temp_ub,(double)BC_size);
				if(temp_ub > ub){
					ub = temp_ub;
				}
				//cout<<ub<<endl<<model.get(GRB_DoubleAttr_ObjVal)<<endl;

				if( model.get(GRB_DoubleAttr_ObjVal)>0 )
				{
					for(i=0;i<cnt;i++)
					{
						if( x[i].get(DVV) >= 0.6) 
							BC1_sol.push_back(vertices[i]);				
					}		
				}
			}

			delete[] x;
		}
		catch(GRBException e) 
		{
			cout << "Error code = " << e.getErrorCode() << endl;
			cout << e.getMessage() << endl;
			if(BC_size>ub){
				ub = BC_size;		// use BC size as the upper bound
			} 
		} 
		catch(...) 
		{
			cout << "Exception during optimization" << endl;
			if(BC_size>ub){
				ub = BC_size;		// use BC size as the upper bound
			} 
		}
	}
	void BC2()									
	{	// BC for <= 2k-2	
		//This is a function that may never be used
		cout<<"Less than or equal to 2k-2, not large enough !!!";
	}

	void peel(int v)
	{
		int top;
		ints::iterator i;
		stack<int> peel_stack;

		deleted[v] = true;
		remain_cnt--;
		peel_stack.push(v);

		while(!peel_stack.empty())
		{
			top = peel_stack.top();
			peel_stack.pop();		
			
			for(i=AdjList[top].begin();i!=AdjList[top].end();i++)
			{
				degs[*i]--;
				if( (!deleted[*i]) && (degs[*i] <= (int)sol.size() - k) )
				{
					deleted[*i] = true;
					remain_cnt--;
					peel_stack.push(*i);
				}
			}
		}	
	}

public:
	
	ipbc(graph& g_obj, int k_value, intv& ini_sol, intv& stat_from_grasp,intv& degs_from_grasp,int left_from_grasp)
		: n(g.vcnt),g(g_obj),AdjList(g.AdjList),ddv(n),deleted(n,false)
	{
		ini_time = omp_get_wtime();

		k = k_value;		

		//---------------------------------//
		// Mainly transfer grasp results
		sol = ini_sol;
		remain_cnt = left_from_grasp;
		degs = degs_from_grasp;	

		for(int i=0;i<n;i++)
		{
			ddv[i] = i;

			if(stat_from_grasp[i] == 3)
				deleted[i] = true;
		}
		//////////////////////////////////////

		sort(ddv.begin(),ddv.end(),deg_greater(degs));

		//-----------------//
		// Initialize stat
		BCtimes = 0;

		BC_time = 0;
		incum_time = 0;
		incum_bc = 0;

		bc_max = 0;
		bc_min = 10000000;
		bc_total = 0;

		ub = (double)ini_sol.size();

		IfBC2 = false;
		/////////////////////

		ini_time = omp_get_wtime() - ini_time;
	}
	~ipbc(void){};

	//------------------------------------------------------------//
	//Statistics
	int opt;				// The optimal size 
	int BCtimes;			// How many times BC is executed
	int StopBC;				// Which BC reach the time limit

	int bc_max;				// The max graph size for BC 
	int bc_min;				// The min graph size for BC
	int bc_total;			// sum of the bc sizes, for time limit

	double ini_time;		// Iinitialization time
	double BC_time;			// Total BC time
	double run_time;		// Runing time

	double ub;				// upper bound;

	double incum_time;		// The time that incumbent is found
	int incum_bc;			// The BC that incumbent is found

	intv sol;				// Incumbent solution
	bool IfBC2;				// If BC2 runs
	/////////////////////////////////////////////////////////////////

	void Run()
	{
		double st = omp_get_wtime(); // starting time of ipbc
		double bc_timer, bc_span, tl_left=IPBCTL;

		intv BC1_sol;
		intv::iterator i;

		GRBEnv env;		
		env.set(OUT, 0);		
		env.set(GRB_IntParam_Cuts, 0);
		env.set(GRB_DoubleParam_Heuristics, 0);
		// CutOff and time limit are set in BC1

		for(i=ddv.begin();i!=ddv.end();i++)
		{
			if( remain_cnt<=(int)sol.size() )
				break;
			
			if(deleted[*i]) continue;

			FindN2(*i);

			if(N2.size() > sol.size())
			{
				PeelN2Proc();

				if(BC_Cand[*i] && BC_size > (int)sol.size() )//BC_list.size() > sol.size())
				{
					if(BC_size > bc_max)
						bc_max = BC_size;
					if(BC_size < bc_min)	
						bc_min = BC_size;

					bc_total += BC_size;
					BCtimes++;

					bc_timer = omp_get_wtime();					
					if(tl_left<=0)
					{	// <=0, no time left
						BC1_sol.clear();	// not going to do BC anymore so empty
						if(BC_size>ub){
							ub = BC_size;		// use BC size as the upper bound
						}
					}
					else
					{
						StopBC = BCtimes;
						BC1(*i, BC1_sol, env, max(tl_left,BCMINT));
					}
					bc_span = omp_get_wtime() - bc_timer;

					BC_time += bc_span;
					tl_left -= bc_span;

					if(BC1_sol.size() > sol.size())
					{
						sol = BC1_sol;
						incum_time = omp_get_wtime()-st;
						incum_bc = BCtimes;	// If 0 means grasp found the best solution !!!
					}
				}
			}

			peel(*i);
		}

		if((int)sol.size() <= 2*k-2) 
		{
			IfBC2 = true;
			BC2();
		}

		run_time = omp_get_wtime()-st;
		opt = (int)sol.size();		
	}
};
//////////////////////////////////////////////////////////////

void GetInput(const string &input, strv &gnames)
{
	string name;
	string path = "input/" + input + ".txt";

	ifstream fin(path.c_str());

	if(!fin)
	{
		cout<<"Cannot open!"<<endl;
		return;
	}

	while (!fin.eof()) 
	{
		getline(fin,name);

		if(name.size()>0)
			gnames.push_back(name);
	}

	fin.close();
}

//---------------------------------------------------------------------------------------------------------------------------//
// find the best grasp result
int best_grasp(int k, int core_cnt, grasp *gs, ofstream &gout, size_t &best_size, bool &opt_proof)
{
	int rtn;
	best_size = 0;
	opt_proof = false;

	for(int i=0;i<core_cnt;i++)
	{
		gout<<k<<",";

		if(gs[i].bestsol.size()>best_size)
		{
			rtn = i;
			best_size = gs[i].bestsol.size();
			opt_proof = gs[i].opt_proof;
		}

		gout<<gs[i].fs<<","<<gs[i].ft<<","<<gs[i].bestsol.size()<<","<<gs[i].bt<<","
			<<gs[i].tt<<","<<gs[i].ct<<","<<gs[i].lt<<","<<gs[i].licnt<<","<<gs[i].iter_cnt<<","
			<<gs[i].lirate<<","<<gs[i].pret<<","<<gs[i].pre_cnt<<","<<gs[i].glo_unex<<","
			<<gs[i].opt_proof<<",{";		
		
		for(intv::iterator j=gs[i].update_iter.begin();j!=gs[i].update_iter.end();j++){
			gout<<*j<<" ";
		}
		
		gout<<"}"<<endl;
	}
	return rtn;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void ParGRASPNoTabu(graph &g, intv &ks, int grasp_timelimit, ofstream &gout)
{	
	intv::iterator k;
	for(k=ks.begin();k!=ks.end();k++)
	{
		cout<<"GRASP Runing for k = "<<*k<<endl;
		
		double par_grasp_time = omp_get_wtime();

		// for 8 cores
		grasp gs0(g,*k,grasp_timelimit);
		grasp gs1(g,*k,grasp_timelimit);
		grasp gs2(g,*k,grasp_timelimit);
		grasp gs3(g,*k,grasp_timelimit);
		grasp gs4(g,*k,grasp_timelimit);
		grasp gs5(g,*k,grasp_timelimit);
		grasp gs6(g,*k,grasp_timelimit);
		grasp gs7(g,*k,grasp_timelimit);
		
		#pragma omp parallel
		{
			#pragma omp sections
			{
				#pragma omp section					
					gs0.Run();

				#pragma omp section					
					gs1.Run();

				#pragma omp section					
					gs2.Run();

				#pragma omp section					
					gs3.Run();

				#pragma omp section					
					gs4.Run();

				#pragma omp section					
					gs5.Run();

				#pragma omp section					
					gs6.Run();

				#pragma omp section					
					gs7.Run();
			}
			#pragma omp barrier
		}

		par_grasp_time = omp_get_wtime() - par_grasp_time;

		cout<<"GRASP Running Done."<<endl;

		grasp gs_array[] = {gs0, gs1, gs2, gs3, gs4, gs5, gs6, gs7};
		size_t best_size; bool opt_proof;
		int best_gs = best_grasp(*k,8,gs_array,gout,best_size,opt_proof);

		gout<<best_size<<","<<par_grasp_time<<","<<opt_proof<<","<<best_gs+1<<endl<<endl;	
	}
}
void ParGRASPSizeTabu(graph &g, intv &ks, int grasp_timelimit, ofstream &gout)
{	
	intv::iterator k;
	for(k=ks.begin();k!=ks.end();k++)
	{
		cout<<"GRASP Runing for k = "<<*k<<endl;
		
		double par_grasp_time = omp_get_wtime();

		// for 8 cores
		grasp gs0(g,*k,grasp_timelimit);
		grasp gs1(g,*k,grasp_timelimit);
		grasp gs2(g,*k,grasp_timelimit);
		grasp gs3(g,*k,grasp_timelimit);
		grasp gs4(g,*k,grasp_timelimit);
		grasp gs5(g,*k,grasp_timelimit);
		grasp gs6(g,*k,grasp_timelimit);
		grasp gs7(g,*k,grasp_timelimit);
		
		#pragma omp parallel
		{
			#pragma omp sections
			{
				#pragma omp section					
					gs0.Run(TABUSIZE);

				#pragma omp section					
					gs1.Run(TABUSIZE);

				#pragma omp section					
					gs2.Run(TABUSIZE);

				#pragma omp section					
					gs3.Run(TABUSIZE);

				#pragma omp section					
					gs4.Run(TABUSIZE);

				#pragma omp section					
					gs5.Run(TABUSIZE);

				#pragma omp section					
					gs6.Run(TABUSIZE);

				#pragma omp section					
					gs7.Run(TABUSIZE);
			}
			#pragma omp barrier
		}

		par_grasp_time = omp_get_wtime() - par_grasp_time;

		cout<<"GRASP Running Done."<<endl;

		grasp gs_array[] = {gs0, gs1, gs2, gs3, gs4, gs5, gs6, gs7};
		size_t best_size; bool opt_proof;
		int best_gs = best_grasp(*k,8,gs_array,gout,best_size,opt_proof);

		gout<<best_size<<","<<par_grasp_time<<","<<opt_proof<<","<<best_gs+1<<endl<<endl;	
	}
}
void ParGRASPPorTabu(graph &g, intv &ks, int grasp_timelimit, ofstream &gout)
{	
	intv::iterator k;
	for(k=ks.begin();k!=ks.end();k++)
	{
		cout<<"GRASP Runing for k = "<<*k<<endl;
		
		double par_grasp_time = omp_get_wtime();

		// for 8 cores
		grasp gs0(g,*k,grasp_timelimit);
		grasp gs1(g,*k,grasp_timelimit);
		grasp gs2(g,*k,grasp_timelimit);
		grasp gs3(g,*k,grasp_timelimit);
		grasp gs4(g,*k,grasp_timelimit);
		grasp gs5(g,*k,grasp_timelimit);
		grasp gs6(g,*k,grasp_timelimit);
		grasp gs7(g,*k,grasp_timelimit);
		
		#pragma omp parallel
		{
			#pragma omp sections
			{
				#pragma omp section					
					gs0.Run(TABU_S_P);

				#pragma omp section					
					gs1.Run(TABU_S_P);

				#pragma omp section					
					gs2.Run(TABU_S_P);

				#pragma omp section					
					gs3.Run(TABU_S_P);

				#pragma omp section					
					gs4.Run(TABU_S_P);

				#pragma omp section					
					gs5.Run(TABU_S_P);

				#pragma omp section					
					gs6.Run(TABU_S_P);

				#pragma omp section					
					gs7.Run(TABU_S_P);
			}
			#pragma omp barrier
		}

		par_grasp_time = omp_get_wtime() - par_grasp_time;

		cout<<"GRASP Running Done."<<endl;

		grasp gs_array[] = {gs0, gs1, gs2, gs3, gs4, gs5, gs6, gs7};
		size_t best_size; bool opt_proof;
		int best_gs = best_grasp(*k,8,gs_array,gout,best_size,opt_proof);

		gout<<best_size<<","<<par_grasp_time<<","<<opt_proof<<","<<best_gs+1<<endl<<endl;	
	}
}
void SingleGRASP(graph &g, intv &ks, int grasp_timelimit, ofstream &gout, int tabu_type) // tabu_type: 0 no tabu; 1 size tabu; 2 portion tabu
{	
	intv::iterator k;
	for(k=ks.begin();k!=ks.end();k++)
	{
		cout<<"GRASP Runing for k = "<<*k<<endl;
		
		double par_grasp_time = omp_get_wtime();

		grasp gs0(g,*k,grasp_timelimit);
		if(tabu_type==0)
			gs0.Run();
		else if(tabu_type==1)
			gs0.Run(TABUSIZE);
		else
			gs0.Run(TABU_S_P);

		par_grasp_time = omp_get_wtime() - par_grasp_time;

		cout<<"GRASP Running Done."<<endl;

		grasp gs_array[] = {gs0};
		size_t best_size; bool opt_proof;
		int best_gs = best_grasp(*k,1,gs_array,gout,best_size,opt_proof);

		gout<<best_size<<","<<par_grasp_time<<","<<opt_proof<<","<<best_gs+1<<endl<<endl;	
	}
}
void GRASPSolveBatch(intv &ks, const string &in_fn, int graph_type=0, const string &ext=".txt", int wsign=1, int wt=-5, bool use8p=true,int tabu_type=0) //go: grasp out, io: ipbc out. wsing: 1 >=w, -1 <=w. tabu_type: 0 no tabu; 1 size tabu; 2 portion tabu
{
	strv gnames;
	strv::iterator i;

	GetInput(in_fn, gnames);

	string path = "output/grasp_"+in_fn+".csv";
	ofstream gout(path.c_str());

	gout<<"Graph/k,1st_sol,1st_time,best_sol,best_time,total_time,Cons_time,LS_time,licnt,total_iter,lirate,pre_time,pre_count,final_rest_nodes,optimal"<<endl;
	
	for(i=gnames.begin();i!=gnames.end();i++)
	{
		cout<<*i<<" reading ..."<<endl;
		gout<<*i<<",";
		
		graph g("input/"+in_fn+"/"+*i+ext,graph_type,wsign,wt);

		gout<<g.vcnt<<","<<g.ecnt<<endl;

		if(use8p)
		{
			if(tabu_type==0)
				ParGRASPNoTabu(g,ks,GRASPTL,gout);
			else if(tabu_type==1)
				ParGRASPSizeTabu(g,ks,GRASPTL,gout);
			else
				ParGRASPPorTabu(g,ks,GRASPTL,gout);
		}
		else
		{
			SingleGRASP(g,ks,GRASPTL,gout,tabu_type);
		}
		
		cout<<"All Done!"<<endl<<endl;

		gout<<endl;
	}

	gout.close();
}
void LoadGRASP(const string &file, map<string,intv> &gr)
{
	string path = "input/" + file;
	ifstream fin(path.c_str());

	if(!fin)
	{
		cout<<"Cannot open!"<<endl;
		return;
	}

	string str, id;	
	
	while (!fin.eof()) 
	{
		getline(fin,str);
		if(str.size()<=1 || str.substr(0,7)=="Graph/k")
			continue; //skip head
		
		strv splits;
		splitstr(splits,str,",");
		id = splits[0];

		int s = atoi(splits[18].c_str());
		gr[id].push_back(s);
	}

	fin.close();
}

//--------------------------------------------------------------------------------------------------------------------------//
void IPBCSolve(map<string,intv> &gr, graph &g, string &gname, intv &ks, ofstream &iout)
{
	intv::iterator k;
	for(k=ks.begin();k!=ks.end();k++)
	{
		grasp gs(g,*k,0);
		gs.load_peel(gr[gname][*k-1]);
		
		cout<<"IPBC Running for k = "<<*k<<endl;

		if(gs.opt_proof)
			iout<<*k<<","<<gs.glo_unex<<","<<gs.bestsol.size()<<",Optimized by GRASP"<<endl;
		else
		{
			ipbc is(g,*k,gs.bestsol,gs.glo_stat,gs.glo_deg,gs.glo_unex);
			is.Run();

			iout<<*k<<","<<gs.glo_unex<<","<<gs.bestsol.size()<<","<<is.opt<<","<<is.ub<<","
				<<is.StopBC<<","<<is.BCtimes<<","<<is.ini_time<<","<<is.BC_time<<","
				<<is.run_time<<","<<is.incum_time<<","<<is.incum_bc<<","<<is.bc_max<<","
				<<is.bc_min<<","<<(double)is.bc_total/is.BCtimes<<","<<is.IfBC2<<endl;
		}

		cout<<"IPBC Running Done."<<endl;
	}
}
void IPBCSolveBatch(map<string,intv> &gr, intv &ks, const string &in_fn, int graph_type=0, const string &ext=".txt", int wsign=1, int wt=-5) 
{
	strv gnames;
	strv::iterator i;

	GetInput(in_fn, gnames);

	string path = "output/ipbc_"+in_fn+".csv";
	ofstream iout(path.c_str());

	iout<<"Graph/k,Nodes# after pre,Ini LB,Best k-plex Size,UB,StopBC,BC times,Initial Time,BC Time,IPBC Time,incumbent time,incumbent bc,BC MAX,BC MIN,BC AVG,IsBC2"<<endl;


	for(i=gnames.begin();i!=gnames.end();i++)
	{
		cout<<*i<<" reading ..."<<endl;
		iout<<*i<<",";

		graph g("input/"+in_fn+"/"+*i+ext,graph_type,wsign,wt);

		iout<<g.vcnt<<","<<g.ecnt<<endl;

		IPBCSolve(gr,g,*i,ks,iout);
		cout<<"All Done!"<<endl<<endl;

		iout<<endl;
	}
	
	iout.close();
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class BestTimeCB: public GRBCallback
{
private: 
	double st, gsol;
	bool is_1st;
public:
	double first_time, first_sol, better_time, first_better_sol, best_time;
	 ;
	BestTimeCB(double _st, double _gsol) 
	{
		st = _st;
		gsol = _gsol;
		is_1st = true;
		first_better_sol = -1;
		better_time = -1;
	}
	~BestTimeCB(){}
 protected:
	void callback () 
	{
		try
		{
			if (where == GRB_CB_MIPSOL)
			{
				best_time = omp_get_wtime() - st;
				double cur_sol = getDoubleInfo(GRB_CB_MIPSOL_OBJ);
				if(is_1st)
				{
					is_1st = false;
					first_time = best_time;
					first_sol = cur_sol;
				}
				if(cur_sol>=gsol-1e-8)
				{
					 first_better_sol = cur_sol;
					 better_time = best_time;
				}
			} 
		}
		catch (GRBException e) 
		{
			cout << "Error number: " << e.getErrorCode() << endl;
			cout << e.getMessage() << endl;
		} 
		catch (...) 
		{
			cout << "Error during callback" << endl;
		}
	}
};
void simple_mip(map<string,intv> &gr, graph &g, string &gname, intv &ks, double time_limit, ofstream &fout)
{
	int i, j, d, vcnt = g.vcnt;
	intv::iterator k;
	double obj, ub;

	GRBEnv env;
	env.set(OUT, 0);
	//env.set(GRB_DoubleParam_NodefileStart,80);	// not helping
	env.set(GRB_DoubleParam_TimeLimit, time_limit);
	//env.set(GRB_IntParam_MIPFocus,1);

	for(k=ks.begin();k!=ks.end();k++)
	{
		cout<<"MIP Runing for k = "<<*k<<endl;

		double rt = omp_get_wtime();

		try
		{
			GRBModel model(env);
			
			GRBVar* x = GRBVarArray(vcnt, model, 0.0, 1.0, GRB_BINARY);
			model.update();

			model.setObjective(GRBSum(vcnt, x), MAX);

			for(i=0; i<vcnt; i++) 
			{	
				d = 0;
				GRBLinExpr v = 0;
				for(j=0; j<vcnt; j++) 
				{
					if (( !g.IsAdjTo(i,j) )&&(i!=j))
					{
						d += 1;
						v += x[j];
					}					
				}

				v -= ((*k-1)*x[i]) + (d*(1-x[i]));			
				model.addConstr(v <= 0);
			}

			BestTimeCB cb = BestTimeCB(rt,(double)gr[gname][*k-1]);
			model.setCallback(&cb);
			model.optimize();

			rt = omp_get_wtime() - rt;

			cout<<"MIP Running Done."<<endl;

			int opt_status = model.get(STAT);
			if( opt_status==GRB_OPTIMAL || opt_status==GRB_TIME_LIMIT )
			{	
				ub = model.get(GRB_DoubleAttr_ObjBound);
				obj = model.get(GRB_DoubleAttr_ObjVal);		
			}

			delete[] x;

			fout<<*k<<","<<gr[gname][*k-1]<<","<<cb.first_sol<<","<<cb.first_better_sol<<","
				<<obj<<","<<ub<<","<<cb.first_time<<","<<cb.better_time<<","<<cb.best_time<<","
				<<rt<<endl<<endl;
		}
		catch(GRBException e) 
		{
			cout << "Error code = " << e.getErrorCode() << endl;
			fout << e.getMessage() << endl;
		} 
		catch(...) 
		{
			fout << "Exception during optimization" << endl;
		}
	}
}
void mip_batch(map<string,intv> &gr, intv &ks, const string &in_fn, int type=0, const string &ext=".txt")
{
	strv gnames;
	strv::iterator i;

	GetInput(in_fn, gnames);

	string path = "output/mip_"+in_fn+".csv";
	ofstream fout(path.c_str());

	fout<<"Graph/k,grasp_sol,1st_sol,1st_better,obj,ub,1st_time,1st_better_time,best_time,totall_time"<<endl;


	for(i=gnames.begin();i!=gnames.end();i++)
	{
		cout<<*i<<" reading ..."<<endl;
		fout<<*i<<",";

		graph g("input/"+in_fn+"/"+*i+ext,type);

		fout<<g.vcnt<<","<<g.ecnt<<endl;

		simple_mip(gr,g,*i,ks,MIPTL,fout);
		cout<<"All Done!"<<endl<<endl;

		fout<<endl;
	}
	
	fout.close();
}

int main()
{	
	map<string,intv> gr;
	LoadGRASP("grasp_snap.csv_sum.csv", gr);
	
	int kvs[] = {1,2,3,4,5,6};
	intv ks( kvs, kvs + sizeof(kvs)/sizeof(int) );	
	
	IPBCSolveBatch(gr, ks,"Autonomous systems");
	IPBCSolveBatch(gr, ks,"Citation network");
	IPBCSolveBatch(gr, ks,"Collaboration networks");
	IPBCSolveBatch(gr, ks,"Communication networks");
	IPBCSolveBatch(gr, ks,"Higgs Twitter Dataset",1,".edgelist");	
	IPBCSolveBatch(gr, ks,"Internet peer-to-peer networks");
	IPBCSolveBatch(gr, ks,"Location-based online social networks");
	IPBCSolveBatch(gr, ks,"Networks with ground-truth communities");
	IPBCSolveBatch(gr, ks,"Product co-purchasing networks");
	IPBCSolveBatch(gr, ks,"Road networks");
	IPBCSolveBatch(gr, ks,"Signed networks",0,".txt",1,0);
	IPBCSolveBatch(gr, ks,"Social Network");		
	IPBCSolveBatch(gr, ks,"Social Network2",1); //for combined social networks	
	IPBCSolveBatch(gr, ks,"Web graphs");	

	//system("pause");
	return 0;
}

