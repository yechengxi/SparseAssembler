#include "iostream"
#include "stdio.h"
#include "string"
#include "vector"
#include "cstdlib"
#include "bitset"
#include <map>
#include <math.h>
#include "memory"
#include <algorithm>
#include "fstream"
#include "sstream"
#include "list"
#include "stdlib.h"
#include "time.h"
#include "BasicDataStructure.h"
#include "BuildContigs.h"
#include "GraphConstruction.h"
#include "GraphSimplification.h"
#include "ScaffoldingDataStructure.h"
#include "ReadsCorrection.h"
#include "ReadsOperation.h"

using namespace std;

int main(int argc, char* argv[])
{
	//cout<<sizeof(kmer_info)<<endl;
	// how to use:		
	cout<<"Command: "<<endl;
	cout<<"Programfile g GAP_VALUE k KMER_SIZE LD LOAD_SKG GS GENOME_SIZE TrimN TRIM_READS_WITH_N f INPUT_FILE1 f INPUT_FILE2 i1 INWARD_PAIR_END1 i2 INWARD_PAIR_END2 o1 OUTWARD_PAIR_END1 o2 OUTWARD_PAIR_END2"<<endl;
	cout<<endl;
	cout<<"Parameters:"<<endl;
	cout<<"k: kmer size, support 15~127."<<endl;
	cout<<"g: number of skipped intermediate k-mers, support 1-25."<<endl;
	cout<<"f: single end reads. Multiple inputs shall be independently imported with this parameter."<<endl;
	cout<<"i1 & i2 (or p1 & p2): inward paired-end reads."<<endl;
	cout<<"o1 & o2 (or l1 & l2): outward paired-end reads."<<endl;
	cout<<"GS: genome size estimation in bp (used for memory pre-allocation), suggest a large value if possible.(e.g. ~ 2x genome size)"<<endl;
	cout<<"LD: load a saved k-mer graph. "<<endl;
	cout<<"BC: 1: build contigs.0: don't build."<<endl;
	cout<<"KmerTable: 1 if you want to output the kmer table."<<endl;
	cout<<"NodeCovTh: coverage threshold for spurious k-mers, support 0-16. (default 1)"<<endl;
	cout<<"EdgeCovTh: coverage threshold for spurious links, support 0-16. (default 0)"<<endl;
	cout<<"PathCovTh: coverage threshold for spurious paths in the breadth-first search, support 0-100."<<endl;
	cout<<"TrimLen: trim long sequences to this length."<<endl;
	cout<<"TrimN: throw away reads with more than this number of Ns."<<endl;
	cout<<"TrimQual: trim off tails with quality scores lower than this."<<endl;
	cout<<"QualBase: lowest quality score value (in ASCII value) in the current fastq scoring system, default: '!'."<<endl;
	cout<<endl;
	cout<<"For error correction:"<<endl;
	cout<<"Denoise: use 1 to call the error correction module. (default 0)"<<endl;
	cout<<"H: hybrid mode. 0 (Default): reads will be trimmed at the ends to ensure denoising accuracy (*MUST* set 0 for the last round). 1: reads will not be trimmed at the ends; "<<endl;
	cout<<"CovTh: coverage threshold for an error. A k-mer with coverage < this value will be checked. Setting 0 will allow the program to choose a value based on the coverage histogram."<<endl;	
	cout<<"CorrTh: coverage threshold for a correct k-mer candidate. A k-mer with coverage >= this value will be considered a candidate for correction. Setting 0 will allow the program to choose a value based on the coverage histogram."<<endl;
	cout<<endl;
	cout<<"For scaffolding:"<<endl;
	cout<<"ExpCov: expected average k-mer coverage in a unique contig. Used for scaffolding."<<endl;
	cout<<"Scaffold: 1: scaffolding with paired reads. 0: single end assembly."<<endl;
	cout<<"LinkCovTh: coverage threshold for spurious paired-end links, support 0-100. (default 5)"<<endl;
	cout<<"Iter_Scaffold: 1: iterative scaffolding using the already built scaffolds (/super contigs). 0: one round scaffolding."<<endl;
	cout<<"For mate pair scaffolding:"<<endl;
	cout<<"InsertSize: estimated insert size of the current pair."<<endl;
	cout<<"i1_mp & i2_mp: inward mate paired reads (large insert sizes >10k, for shorter libraries omit \"_mp\")."<<endl;
	cout<<"o1_mp & o2_mp : outward paired-end reads (large insert sizes >10k, for shorter libraries omit \"_mp\")."<<endl;
	
	///cout<<endl<<"Miscellaneous: "<<endl;

	cout<<endl;

	//debug
	/*
	string reads_info_name="NonContainedReads_info.txt";
	reads_table reads_table0;
	ConstructReadsOverlaps( reads_info_name,&reads_table0);
	*/
	//debug


	// parameters to be updated
	size_t hashTableSZ=1000000;
	int gap=0;
	uint64_t GenomeSize=0,totReads=0;
	uint64_t tot_bases=0;
	int MaxReadLen=0;
	int K_size=31,InsertSize=0;
	char QS_base='!';
	bool BUILD_GRAPH=1;

	bool BFS=1,BUILD_SCAFFOLDS=0,BUILD_CONTIGS=1,ResumePE=0,LOAD_GRAPH=0,LOAD_GRAPH1=0,GET_BRANCHES_INFO=0,BUILD_CONTIG_GRAPH=0,RESUME=0,BUILD_LONG_CONTIGS=0,Resolving_Branches_PE=0,SC_Rev_Comp=0,BUILD_SUPER_CONTIGS=0,LOAD_PE_DIST=0,SINGLE=0,Minimizer=0,Compress=0,KmerTable=0;
	uint8_t SpuriousTh=0;
	int NodeCovTh=1,EdgeCovTh=0,LinkCovTh=5,PathCovTh=-1,PathSim=5,MaxDepth=15,UniqueLenTh=150,TrimN=10,TrimLen=1000000000,TrimQual=-1;
	bool OLC=0;
	bool DBG2OLC=0;
	int ExpCov=0;
	//parameters for the bloom filter
	bool Bloom=1;
	int Coverage=0,AvgLen=100,Bloom_d=-1;
	double ErrorRate=0.01;
	double BF_FalsePositive=0.001;
	double two=2;
	bool Iter_Scaffold=0,MP_Scaffold=0;
	bool CollectNonContainedReads=0;
	//parameters for error correction
	size_t CovTh=5,CorrTh=5;
	bool Hybrid=0;
	map<int,int> cov_hist;
	bool Denoise=0;

	bool MemoryEfficient=1;


	struct read_index read_index;
	read_index.repeat_cnt=0;
	read_index.repeat_maps.clear();

	vector<int> insert_sz_vt;
	vector<bool> OutwardLib,SingleOutwardLib,OutwardLibMP;
	string filename,ContigsForScaffolds="Contigs.txt",bfilename;
	ofstream o_LogBases("BasesCount.txt");
	vector<string > sp_filenames_vt,p_filenames_vt,p1_filenames_vt,p2_filenames_vt,in_filenames_vt,single_filenames_vt,mp1_filenames_vt,mp2_filenames_vt,mp_filenames_vt;
	

	struct contigs_info  contigs_info,scaffolds_info;
	string ContigPrefix="";

	// update the parameters using input arguments
	for(int i=1;i<argc;++i)
	{
		
		if(strcmp(argv[i],"g")==0)
		{
			i++;
			gap=atoi(argv[i]);
			continue;
		}
		if(strcmp(argv[i],"Denoise")==0)
		{
			i++;
			Denoise=atoi(argv[i]);
			continue;
		}
		if(strcmp(argv[i],"ME")==0)
		{
			i++;
			MemoryEfficient=atoi(argv[i]);
			continue;
		}
		if(strcmp(argv[i],"MP_Scaffold")==0)
		{
			i++;
			MP_Scaffold=atoi(argv[i]);
			continue;
		}
		if(strcmp(argv[i],"KmerTable")==0)
		{
			i++;
			KmerTable=atoi(argv[i]);
			continue;
		}
		if(strcmp(argv[i],"Iter_Scaffold")==0)
		{
			i++;
			Iter_Scaffold=atoi(argv[i]);
			continue;
		}
		if(strcmp(argv[i],"Compress")==0)
		{
			i++;
			Compress=atoi(argv[i]);
			continue;
		}
		if(strcmp(argv[i],"g")==0)
		{
			i++;
			gap=atoi(argv[i]);
			continue;
		}
		if(strcmp(argv[i],"H")==0)
		{
			i++;
			Hybrid=atoi(argv[i]);
			continue;
		}
		
		if(strcmp(argv[i],"InsertSize")==0)
		{
			i++;
			InsertSize=atoi(argv[i]);
			insert_sz_vt.push_back(InsertSize);
			continue;
		}

		if(strcmp(argv[i],"RS")==0)
		{
			i++;
			RESUME=atoi(argv[i]);
			
			continue;
		}
		if(strcmp(argv[i],"ExpCov")==0)
		{
			i++;
			ExpCov=atoi(argv[i]);
			continue;
		}
		if(strcmp(argv[i],"CovTh")==0)
		{
			i++;
			CovTh=atoi(argv[i]);
			continue;
		}
		if(strcmp(argv[i],"CorrTh")==0)
		{
			i++;
			CorrTh=atoi(argv[i]);
			continue;
		}
		if(strcmp(argv[i],"AvgLen")==0)
		{
			i++;
			AvgLen=atoi(argv[i]);
			continue;
		}
		if(strcmp(argv[i],"Coverage")==0)
		{
			i++;
			Coverage=atoi(argv[i]);
			continue;
		}
		if(strcmp(argv[i],"ErrorRate")==0)
		{
			i++;
			ErrorRate=atof(argv[i]);
			continue;
		}
		if(strcmp(argv[i],"BF_FalsePositive")==0)
		{
			i++;
			BF_FalsePositive=atof(argv[i]);
			continue;
		}


		if(strcmp(argv[i],"LD")==0)
		{
			i++;
			LOAD_GRAPH=atoi(argv[i]);
			continue;
		}
		if(strcmp(argv[i],"LD1")==0)
		{
			i++;
			LOAD_GRAPH1=atoi(argv[i]);
			continue;
		}
		if(strcmp(argv[i],"Single")==0)
		{
			i++;
			SINGLE=atoi(argv[i]);
			continue;
		}
		if(strcmp(argv[i],"LDPE")==0)
		{
			i++;
			LOAD_PE_DIST=atoi(argv[i]);
			continue;
		}
		if(strcmp(argv[i],"BC")==0)
		{
			i++;
			BUILD_CONTIGS=atoi(argv[i]);
			continue;
		}
		if(strcmp(argv[i],"BFS")==0)
		{
			i++;
			BFS=atoi(argv[i]);
			continue;
		}
		if(strcmp(argv[i],"PathSim")==0)
		{
			i++;
			PathSim=atoi(argv[i]);
			continue;
		}
		if(strcmp(argv[i],"DBG2OLC")==0)
		{
			i++;
			DBG2OLC=atoi(argv[i]);
			continue;
		}
		
		if(strcmp(argv[i],"k")==0)
		{
			i++;
			K_size=atoi(argv[i]);
			continue;
		}
		if(strcmp(argv[i],"TrimN")==0)
		{
			i++;
			TrimN=atoi(argv[i]);
			continue;
		}
		if(strcmp(argv[i],"TrimQual")==0)
		{
			i++;
			TrimQual=atoi(argv[i]);
			continue;
		}
		if(strcmp(argv[i],"TrimLen")==0)
		{
			i++;
			TrimLen=atoi(argv[i]);
			continue;
		}
		if(strcmp(argv[i],"QualBase")==0)
		{
			i++;
			QS_base=*(argv[i]);
			continue;
		}
		if(strcmp(argv[i],"NodeCovTh")==0)
		{
			i++;
			NodeCovTh=atoi(argv[i]);
			continue;
		}
		if(strcmp(argv[i],"PathCovTh")==0)
		{
			i++;
			PathCovTh=atoi(argv[i]);
			continue;
		}
		if(strcmp(argv[i],"UniqueLenTh")==0)
		{
			i++;
			UniqueLenTh=atoi(argv[i]);
			continue;
		}
		if(strcmp(argv[i],"EdgeCovTh")==0)
		{
			i++;
			EdgeCovTh=atoi(argv[i]);
			continue;
		}
		if(strcmp(argv[i],"RSPE")==0)
		{
			i++;
			ResumePE=atoi(argv[i]);
			continue;
		}
		if(strcmp(argv[i],"LinkCovTh")==0)
		{
			i++;
			LinkCovTh=atoi(argv[i]);
			continue;
		}
		if(strcmp(argv[i],"MaxDepth")==0)
		{
			i++;
			MaxDepth=atoi(argv[i]);
			continue;
		}
		if(strcmp(argv[i],"Minimizer")==0)
		{
			i++;
			Minimizer=atoi(argv[i]);
			continue;
		}
		if(strcmp(argv[i],"Bloom")==0)
		{
			i++;
			Bloom=atoi(argv[i]);
			continue;
		}
		if(strcmp(argv[i],"Bloom_d")==0)
		{
			i++;
			Bloom_d=atoi(argv[i]);
			continue;
		}
		if(strcmp(argv[i],"f")==0)
		{
			i++;
			filename=(argv[i]);
			single_filenames_vt.push_back(filename);
			in_filenames_vt.push_back(filename);
			SingleOutwardLib.push_back(0);

			continue;
		}
		if(strcmp(argv[i],"p1")==0||strcmp(argv[i],"i1")==0)
		{
			i++;
			filename=(argv[i]);
			p1_filenames_vt.push_back(filename);
			in_filenames_vt.push_back(filename);
			p_filenames_vt.push_back(filename);
		
			
			SingleOutwardLib.push_back(0);
			OutwardLib.push_back(0);

			continue;
		}


		if(strcmp(argv[i],"p2")==0||strcmp(argv[i],"i2")==0)
		{
			i++;
			filename=(argv[i]);
			p2_filenames_vt.push_back(filename);
			in_filenames_vt.push_back(filename);
			p_filenames_vt.push_back(filename);
			
			continue;
		}

		
		if(strcmp(argv[i],"i1_mp")==0)
		{
			i++;
			filename=(argv[i]);
			mp1_filenames_vt.push_back(filename);
			mp_filenames_vt.push_back(filename);
			OutwardLibMP.push_back(0);

			continue;
		}
		if(strcmp(argv[i],"i2_mp")==0)
		{
			i++;
			filename=(argv[i]);
			mp2_filenames_vt.push_back(filename);
			mp_filenames_vt.push_back(filename);
			OutwardLibMP.push_back(0);

			continue;
		}

		if(strcmp(argv[i],"l1")==0||strcmp(argv[i],"o1")==0)
		{
			i++;
			filename=(argv[i]);
			p1_filenames_vt.push_back(filename);
			in_filenames_vt.push_back(filename);
			p_filenames_vt.push_back(filename);
			
			
			OutwardLib.push_back(1);
			continue;
		}
		if(strcmp(argv[i],"l2")==0||strcmp(argv[i],"o2")==0)
		{
			i++;
			filename=(argv[i]);
			p2_filenames_vt.push_back(filename);
			in_filenames_vt.push_back(filename);
			p_filenames_vt.push_back(filename);
			
			continue;
		}

		if(strcmp(argv[i],"o1_mp")==0)
		{
			i++;
			filename=(argv[i]);
			mp1_filenames_vt.push_back(filename);
			OutwardLibMP.push_back(1);

			continue;
		}
		if(strcmp(argv[i],"o2_mp")==0)
		{
			i++;
			filename=(argv[i]);
			mp2_filenames_vt.push_back(filename);
			OutwardLibMP.push_back(1);

			continue;
		}


		if(strcmp(argv[i],"GS")==0)
		{
			i++;
			//GenomeSize=atoi(argv[i]);
			stringstream strValue;

			strValue << argv[i];
			strValue >> GenomeSize;

			//GenomeSize=atoi(argv[i]);
			
			continue;
		}

		if(strcmp(argv[i],"R")==0) //only use this number of reads for the assembly
		{
			i++;
			//sscanf(argv[i],"%Lu",&totReads);
			stringstream strValue;

			strValue << argv[i];
			strValue >> totReads;
			
			continue;
		}


		if((strcmp(argv[i],"ResolveBranchesPE")==0)||strcmp(argv[i],"Scaffold")==0)
		{
			i++;
			Resolving_Branches_PE=atoi(argv[i]);
			
			continue;
		}




		if(strcmp(argv[i],"ContigPrefix")==0)
		{
			i++;
			ContigPrefix=(argv[i]);
			continue;
		}


		if(strcmp(argv[i],"bf")==0)
		{
			i++;
			bfilename=(argv[i]);
			
			continue;
		}


	}

	vector<int> N50_lst,TotLength_lst,MaxLen_lst;
	if (NodeCovTh > 0 & PathCovTh < 0)
	{
		PathCovTh = NodeCovTh * 10;
	}
	


	int num_inputs=0;
	char in_idx[100],prefix[100];

	while(1)
	{
		if(!bfilename.empty())
		{

			num_inputs++;
			prefix[0]='S';
			prefix[1]='\0';
			//itoa(num_inputs,in_idx,10);
			sprintf(in_idx,"%d",num_inputs);
			in_idx[strlen(in_idx)+1]='\0';
			in_idx[strlen(in_idx)]='_';
			strcat(prefix,in_idx);
			ContigPrefix=prefix;
			ifstream reads_in;
			strcat(prefix,bfilename.c_str());
			reads_in.open(prefix);
			string str;
			if(!getline(reads_in,str))
			{break;}
			else
			{
				reads_in.close();
				in_filenames_vt.clear();
				in_filenames_vt.push_back(prefix);
			}
		}
		if(K_size%2==0)
		{
			K_size--;
		}
		//if g not set
		if(gap==0)
		{
			gap=11;
			if(K_size<21)
			{
				gap=5;
			}

			if(K_size<17)
			{
				gap=max(K_size-2,1);
			}
			if(K_size>31)
			{
				gap=16;
			}
	
		}


		int Kmer_arr_sz=K_size/32+1;
		int rem1=K_size%32;
		if(rem1==0)
		{Kmer_arr_sz--;}
	
		if((gap>25&&K_size>31)||gap>64)
		{
			cout<<"Gap value too large, set to 25."<<endl;
			gap=25;
		}

		if(TrimQual>0)
		{
			QS_base+=TrimQual;
		}

		//initialize the hashtable structure

		string FilePrefix;
		if(Denoise)
		{
			FilePrefix="Dn_";
		}
		else
		{
			FilePrefix="Asm_";
		}
		struct hashtable ht;
		struct hashtable2 ht2;
		struct hashtable3 ht3;
		struct hashtable4 ht4;
		struct hashtable0 ht0;
		struct key_table key_table;
		key_table.KeysPerBlock=100000;
		key_table.current_block=0;
		key_table.current_index=key_table.KeysPerBlock;
		key_table.pblocks.clear();

		

	
		time_t beg_time,read_time;
		time(&beg_time);
		int64_t bucket_count=0,edge_cnt=0;
	
		uint64_t KmerEstNo;
	
		struct BF_info BF_info;
		BF_info.Bloom=0;
		uint64_t numReads=0;
		//cout<<numReads<<endl;
		//bloom filter
		if(NodeCovTh==0)
		{Bloom=0;}

		if((!LOAD_GRAPH)&&BUILD_GRAPH)
		{
			ofstream o_Log("Assembly_Log.txt");

			if(!LOAD_GRAPH1)
			{
				if(GenomeSize==0)
				{cout<<"Error! Genome size not given."<<endl;return -1;}

				tot_bases=0;
				if(Coverage<=0)
				{
					ScanDataset(in_filenames_vt,&tot_bases,&numReads,&totReads,&MaxReadLen);
			
					cout<<"Scan finished."<<endl;
					//cout<<numReads<<endl;
				}
			
				if(numReads>0)
				{
					AvgLen=(tot_bases/numReads);

					int mod=tot_bases%numReads;
					if(mod>(numReads/2))
					{AvgLen++;}
					if(tot_bases/GenomeSize>Coverage)
					{Coverage=tot_bases/GenomeSize;}
			
					o_Log<<"AvgLen: "<<AvgLen<<endl;
					o_Log<<"tot_bases: "<<tot_bases<<endl;
					o_Log<<"Coverage: "<<Coverage<<endl;
				}
				if(Coverage<=0)
				{Coverage=10;}

		
				if(AvgLen<=K_size+gap)
				{
					cout<<"Error! K-mer size too large!"<<endl;
					return -1;
				}



				KmerEstNo=(uint64_t) (((double)(GenomeSize*Coverage*(AvgLen-K_size+1)/AvgLen))*(1-pow((1-ErrorRate),K_size))+GenomeSize);
	
				o_Log<<"KmerEstNo: "<<KmerEstNo<<endl;
			
				if (Bloom)
				{

					uint64_t KmerEstNo2=KmerEstNo;
					if(!Minimizer)
					{
						KmerEstNo2=KmerEstNo2;
					}
					BF_info.m=(uint64_t) (((double)KmerEstNo2/gap)*(-log(BF_FalsePositive)/log(two)/log(two))+100);		
					BF_info.BF_HT=(uint8_t *)calloc(BF_info.m/8+10,sizeof(uint8_t));
					BF_info.d=min(max((int)(7*(BF_info.m)/KmerEstNo/10+0.5),2),2);
					if(Bloom_d>0)
					{	
						BF_info.d=Bloom_d;
					}

					//cout<<BF_info.d<<endl;
					BF_info.Bloom=1;
		

				}
			}
			numReads=0;
			tot_bases=0;
			//build the graph
			if(GenomeSize==0)
			{cout<<"Error! Genome size is not given."<<endl;return -1;}
	
			uint64_t ht_sz_tmp= KmerEstNo/max(gap-2,5);
			if(Bloom==1)
			{ht_sz_tmp= GenomeSize*2/max(gap-2,5);}
			if(ht_sz_tmp>((size_t) -1) )
			{
				cout<<"Problem scale too large. Reduce genome size."<<endl;
				return -1;

			}

			hashTableSZ=(size_t) ht_sz_tmp;
		
			uint64_t TotalSamplings=0;	
			//allocate memory for the hashtable
			if(!LOAD_GRAPH1)
			{
				if(MemoryEfficient&(K_size<=32))//K_size<=32)
				{
					Init_HT(&ht,hashTableSZ);
				}
				else
				{


					Init_HT0(&ht0,hashTableSZ);	
					/*

					if(K_size>32&&K_size<=64)
					{
						Init_HT2(&ht2,hashTableSZ);	
					}
					else
					{
						if(K_size<=96)
						{
							Init_HT3(&ht3,hashTableSZ);	
						}
						else
						{
							if(K_size<=128)
							{
								Init_HT4(&ht4,hashTableSZ);	
							}
						}
					}
					*/
				}
		
			}
			else
			{
				//load a saved graph
				if(LOAD_GRAPH1)
				{
					cout<<endl<<"Loading the graph..."<<endl;
					if(MemoryEfficient&(K_size<=32))
					{
						LoadingSparseKmerGraph(&ht,FilePrefix);
					}
					else
					{
						LoadingSparseKmerGraph0(&ht0,&key_table,FilePrefix,Kmer_arr_sz);
						//to be completed
						/*
						if(K_size>32&&K_size<=64)
						{
							LoadingSparseKmerGraph2(&ht2,FilePrefix);
						}
						else
						{
							if(K_size<=96)
							{
								LoadingSparseKmerGraph3(&ht3,FilePrefix);
							}
							else
							{
								if(K_size<=128)
								{
									LoadingSparseKmerGraph4(&ht4,FilePrefix);
								}
							}
						}

						*/
					}
	
				}


		
			}
	
			struct read_t read;
			read.read_bits =(uint64_t*) malloc(MaxReadLen/4+100);

			int readLen;
			string seq_s,str,tag1,tag2;
			//2 rounds to build the graph
			numReads=0;
			for (int round=1;round<=2;++round)
			{
				//LOAD_GRAPH1 for high end usage, omit
				if(LOAD_GRAPH1)
				{
					round=2;
				}
				cout<<"Building the sparse k-mer graph, round: "<<round<<endl;
				//scan files of fasta or fastq format. The below is long and trivial, skip it to the major function. 
				for(size_t jj=0;jj<in_filenames_vt.size();++jj)
				{
					uint64_t nLines=0;
					int seq_sz=0;

					ifstream infile(in_filenames_vt[jj].c_str());
		
					cout<<jj+1<<"/"<<in_filenames_vt.size()<<" files."<<endl;
					cout<<"Processing file: "<<in_filenames_vt[jj]<<endl;
		
					seq_s.clear();

					bool fq_flag=0;


					getline(infile,str);
					if(fq_flag==0&&str[0]=='@')
					{
						fq_flag=1;	
					}
					infile.close();

					infile.clear();
					infile.open(in_filenames_vt[jj].c_str());

					bool read_success=0;

					read_success=1;

					string tag,qs,n_tag;
					string QS_s;

					while(read_success)
					{
						if(fq_flag)
						{
							read_success=get_a_fastq_read(infile,tag,seq_s,QS_s);
					
						}
						else
						{
							read_success=get_a_fasta_read(infile,tag,seq_s,n_tag);
			
						}	
						if(read_success==0)
						{break;}
				
						seq_sz=seq_s.size();
						if(TrimLen<seq_sz)
						{
							seq_s.resize(TrimLen);
							seq_sz=TrimLen;
						}
					
						if (seq_s.size()==0)
						{
							cout<<"Empty sequence!"<<endl;
							continue;
						}
						bool bad_flag=0;
						int numN=0;
						for(int i=0;i<seq_sz;++i)
						{
							if(seq_s[i]!='A'&&seq_s[i]!='C'&&seq_s[i]!='G'&&seq_s[i]!='T'&&seq_s[i]!='N')
							{
								bad_flag=1;
								break;
							}
							if(seq_s[i]=='N')
							{
								numN++;
							}
						}
						if(numN>TrimN)
						{bad_flag=1;}
						if(bad_flag)
						{continue;}
							

						bad_flag=0;
					
						numReads++;
						//cout<<numReads<<endl;

						char QS_seq[1000];
					
						if(fq_flag)
						{
						
							if(QS_s[QS_s.size()-1]=='\n'||QS_s[QS_s.size()-1]=='\r')
							{
								QS_s.resize(QS_s.size()-1);
							}
							strcpy(QS_seq,QS_s.c_str());
							if(QS_s.size()>seq_sz)
							{
								QS_s.resize(seq_sz);
							}
						}


						if(totReads!=0&&numReads>totReads)
						{
						
							numReads=0;
							break;
						
						}

						if(round==1)
						{
							tot_bases+=seq_sz;
						}


						if(fq_flag&&TrimQual>0)
						{
							int QS_sz=strlen(QS_seq);
							for(int i=0;i<QS_sz;++i)
							{
								if(QS_seq[i]<QS_base)
								{
									seq_s[i]='\0';
									seq_s.resize(i);
									seq_sz=i;
									break;
								}
							}
						}
				

						int nN=seq_sz-1,isN=-1;
						for(int i=0;i<seq_sz;++i)
						{
						
							if(seq_s[i]=='-'||seq_s[i]=='N')
							{
								if(i<=seq_sz/2)
								{
									isN=i;
									continue;
								}
								else
								{
									nN=i-1;
									break;
								}
							}
						}
						int s=0;
						if((nN-isN)<=seq_sz/2)
						{
							bad_flag=1;
						}
					
						if(bad_flag==1)
						{
							seq_s.clear();
							continue;
						}

						if(isN>=0)
						{
							for(int i=isN+1;i<=nN;++i)
							{
								seq_s[s]=seq_s[i];
								s++;
							}
							seq_s[s]='\0';
							seq_s.resize(s);
						}
					
				
						if(Compress)
						{
							string_shrinkage(seq_s);
							seq_sz=seq_s.size();
						}

						if(bad_flag==1||seq_sz<(K_size+1))
						{seq_s.clear();continue;}
		
					
						//cout<<numReads<<endl;
						Init_Read(seq_s,read);
						//cout<<numReads<<endl;
						//strcpy(read.c_seq,seq_s.c_str());
						read.readLen=seq_s.size();

						seq_s.clear();
						size_t b2=bucket_count;
						if(read.readLen<K_size+gap)
						{seq_s.clear();continue;}



						///jump to here, the graph construction begins:
						//if(numReads==4902)
						//cout<<numReads<<endl;
						int ref_pos=0;
						//read.read_idx=numReads;
						



						if(MemoryEfficient&(K_size<=32))
						{
							uint64_t OverlapKmers=read.readLen-K_size+1;
							TotalSamplings+=OverlapKmers;
							Sparse_Kmer_Graph_Construction(&read,&ht,&bucket_count,&edge_cnt, K_size,gap,&BF_info,round);
							//Sparse_Kmer_Index_Construction(&read,&ht,&bucket_count, K_size,gap,&BF_info,round,ref_pos,&read_index);
						}
						else
						{

							uint64_t OverlapKmers=read.readLen-K_size+1;
							TotalSamplings+=OverlapKmers;
							Sparse_Kmer_Graph_Construction0(&read,&ht0,&key_table,&bucket_count,&edge_cnt,K_size, gap,&BF_info,round);

						}
				
				
						if (numReads%10000000==0)
						//if (numReads%100000==0)
						{
							cout<<"Reading: "<<numReads<<endl;
							if(round==1)
							{
								if(MemoryEfficient&(K_size<=32))
								{	cout<<"Memory used: "<<(sizeof(struct bucket_r1)*bucket_count)/1024/1024<<" MB."<<endl;}
								else
								{
									cout<<"Memory used: "<<((sizeof(struct bucket0_r1)+sizeof(uint64_t)*Kmer_arr_sz)*bucket_count)/1024/1024<<" MB."<<endl;
									
								}
							}
							else
							{
								if(MemoryEfficient&(K_size<=32))
								cout<<"Memory used: "<<(sizeof(struct bucket)*bucket_count+sizeof(struct edge_node)*edge_cnt)/1024/1024<<" MB."<<endl;
								else
								{

					
								}
							}

							//cout<<"bucket_count"<<bucket_count<<endl;
							time(&read_time);
							cout<<"Reading time: "<<difftime(read_time,beg_time)<<" secs."<<endl;

						}
				
			

					}
					//fclose(infile);
					infile.close();
					infile.clear();
				 
				 

				}


			
				if(round==1)
				{

					cout<<"Total bases: "<<tot_bases<<endl;
					o_LogBases<<"Total bases: "<<tot_bases<<endl;
				}
				if(numReads>0&&round==1)
				{
					int rlen=tot_bases/numReads;
					int mod=tot_bases%numReads;
					if(mod>(numReads))
					{rlen++;}

					o_LogBases<<"Avg read length: "<<rlen<<endl;
				}
				numReads=0;

				if(bucket_count>0)
				{
					cout<<"Total nodes in round "<<round<<": "<<bucket_count<<endl;
					o_Log<<"Total nodes in round "<<round<<": "<<bucket_count<<endl;

				}


				if(round==2&&bucket_count>0)//round==2&&
				{
					cout<<"Total nodes: "<<bucket_count<<endl;
					cout<<"Total edges in round "<<round<<": "<<edge_cnt<<endl;
					o_LogBases<<"Total nodes: "<<bucket_count<<endl;
					o_Log<<"Total edges in round "<<round<<": "<<edge_cnt<<endl;
				}


				//remove weak links and nodes.

				if(round==1)
				{
					//free BF
					if(Bloom)
					{
						free(BF_info.BF_HT);
					}

					if(MemoryEfficient&(K_size<=32))
					{
						if(Bloom&&(NodeCovTh>0))
						{
							RemovingWeakNodes_r1(&ht,&ht2, K_size,NodeCovTh-1, &bucket_count);
						}
						else
						{
							RemovingWeakNodes_r1(&ht,&ht2, K_size,NodeCovTh, &bucket_count);
						}
						SwitchBuckets(&ht,&ht2,K_size);
					}
					else
					{

						if(Bloom&&(NodeCovTh>0))
						{
							RemovingWeakNodes0_r1(&ht0, K_size,NodeCovTh-1, &bucket_count);
						}
						else
						{
							RemovingWeakNodes0_r1(&ht0, K_size,NodeCovTh, &bucket_count);
						}
						SwitchBuckets0(&ht0,K_size);

					
					}
				}
				// save to disk
			
				if(MemoryEfficient&(K_size<=32))
				{
					
					SavingSparseKmerGraph(&ht,FilePrefix);
					if(KmerTable)
					{
						OutputSparseKmers(&ht,K_size,Bloom);
					}

				}
				else
				{

					SavingSparseKmerGraph0(&ht0,FilePrefix,Kmer_arr_sz);

				
				}

			
				if(round==1)
				{
					if(MemoryEfficient&(K_size<=32))
					{
						for(size_t i=0;i<ht.ht_sz;++i)
						{				
							struct bucket* bktptr=ht.store_pos[i];
							while(bktptr!=NULL)
							{
								bktptr->kmer_info.cov1=0;	
								bktptr=bktptr->nxt_bucket;
							}
					
						}
				
					}
					else
					{


						for(size_t i=0;i<ht0.ht_sz;++i)
						{				
							struct bucket0* bktptr=ht0.store_pos[i];
							while(bktptr!=NULL)
							{
								bktptr->kmer_info.cov1=0;	
								bktptr=bktptr->nxt_bucket;
							}
					
						}
					
					}
				}


				time(&read_time);
				cout<<"Reading time: "<<difftime(read_time,beg_time)<<" secs."<<endl;


			}
	

			//free(read.read_bits);

			uint64_t MeanCov=TotalSamplings/GenomeSize;
		
		
		}

		hashtable merge_ht;
		hashtable2 merge_ht2;
		hashtable3 merge_ht3;
		hashtable4 merge_ht4;		
		hashtable0 merge_ht0;
		merge_ht.ht_sz=0;
		merge_ht2.ht_sz=0;
		merge_ht3.ht_sz=0;
		merge_ht4.ht_sz=0;

		merge_ht0.ht_sz=0;
		struct key_table merge_key_table;
		merge_key_table.KeysPerBlock=100000;
		merge_key_table.current_block=0;
		merge_key_table.current_index=key_table.KeysPerBlock;
		merge_key_table.pblocks.clear();

			

		if(LOAD_GRAPH)
		{
			cout<<endl<<"Loading the graph..."<<endl;
			if(MemoryEfficient&(K_size<=32))
			{
			
				LoadingSparseKmerGraph(&ht,FilePrefix);
				
				if(!BFS)
				{
					LoadingMergeHT(&merge_ht);
				}
			}
			else
			{

				LoadingSparseKmerGraph0(&ht0,&key_table,FilePrefix,Kmer_arr_sz);
				
				if(!BFS)
				{
					LoadingMergeHT0(&merge_ht0,&merge_key_table,Kmer_arr_sz);
				}

				
			}
	
		}


		if((BFS||Denoise)&&(!Iter_Scaffold))
		{
			cout<<"Screening off weak nodes and edges..."<<endl;

		
			if(MemoryEfficient&(K_size<=32))
			{
				RemovingWeakNodesAndEdges(&ht, K_size,NodeCovTh, EdgeCovTh,&bucket_count, &edge_cnt);
			}
			else
			{
				RemovingWeakNodesAndEdges0(&ht0, K_size,NodeCovTh, EdgeCovTh,&bucket_count, &edge_cnt);
				
			}
		}

		if(Denoise)
		{
			//denoise



			ifstream in_cov("CovHist.txt");
			map<int,int >::iterator mit;
			int cov,cnt;
			//int cov_vt[1000];
			int cnt_vt[200];
			//memset(cov_vt,0,sizeof(cov_vt));
			memset(cnt_vt,0,sizeof(cnt_vt));

			int max_cnt=0,cov_max_cnt=0;
			while(in_cov>>cov>>cnt)
			{
				cov_hist[cov]=cnt;
			
				if(cov<200)
				{
					cnt_vt[cov]=cnt;
				}
			}
			bool start=0;
			int Th=2;
			if(NodeCovTh>=Th-1)
			{
				Th=NodeCovTh+1;
			}
			for(int i=Th;i<100;++i)
			{
				if(start==0&&(cnt_vt[i]<cnt_vt[i-1]||cnt_vt[i+1]<cnt_vt[i]))
				{
					continue;
				}
				else
				{
					start=1;
				}
				if(start==1&&cnt_vt[i]>max_cnt)
				{max_cnt=cnt_vt[i];cov_max_cnt=i;}

			}
			cov_max_cnt=min(cov_max_cnt,100);
			ofstream o_DnLog("DenoisingLog.txt",ios::app);
			o_DnLog<<"max diversity: "<<cov_max_cnt<<endl;
		

		
			for(int i=Th;i<cov_max_cnt;++i)
			{
				if((cnt_vt[i]<cnt_vt[i-1]||cnt_vt[i+1]<cnt_vt[i])||i<(cov_max_cnt/15))
				{
					continue;
				}
				else
				{
					o_DnLog<<"Suggested CovTh: "<<i-1<<endl;
					if(CovTh==0)
					{
						CovTh=i+1;
						CorrTh=i+1;
					}
					cout<<"CovTh used: "<< CovTh<<endl;
					break;
				}
			
			}
		


			time(&read_time);

			cout<<"Read time: "<<difftime(read_time,beg_time)<<" secs."<<endl;

			cout<<"Start denoising."<<endl;

			if(K_size<=32)
			{
				MarkBranches(&ht);
			}
			else
			{
				if(K_size>32&&K_size<=64)
				{
					MarkBranches2(&ht2);
	
				}
		
			}
			numReads=0;
			uint64_t correction_cnt=0;
			//string temp,temp1;
			struct read_t read1,read2;
			read1.read_bits =(uint64_t*) malloc(MaxReadLen/4+100);
			read2.read_bits =(uint64_t*) malloc(MaxReadLen/4+100);
			bucket_count=0;


			string seq1,seq2,str,tag1,tag2,seq_s;
	
			uint64_t nLines=0;
			int seq_sz=0;
			
			for(int jj=0;jj<single_filenames_vt.size();++jj)
			{
			
				ifstream infile(single_filenames_vt[jj].c_str());
				bool fq_flag=0;
				nLines=0;


				string t_filename;
			//	t_filename=p1_filenames_vt[jj];
				int s,t=0;
				for(s=(int) single_filenames_vt[jj].size()-1;s>=0;--s)
				{
					if((single_filenames_vt[jj][s]== '\\' )||( single_filenames_vt[jj][s]=='/'))
					{
						t=s+1;
						break;
					}
				}
				t_filename=single_filenames_vt[jj].substr(t,single_filenames_vt[jj].size());
				string denoised_name="Denoised_"+t_filename;
		

				ofstream o_denoised(denoised_name.c_str());
				cout<<jj+1<<"/"<<single_filenames_vt.size()<<" files."<<endl;
				cout<<"Processing file: "<<single_filenames_vt[jj]<<endl;
		

				getline(infile,str);


				if(fq_flag==0&&str[0]=='@')
				{
					fq_flag=1;	
				}
				infile.close();

				infile.clear();
				infile.open(in_filenames_vt[jj].c_str());

				bool read_success=0;

				read_success=1;

				string tag,qs,n_tag;
				string QS_s;

				while(read_success)
				{
					if(fq_flag)
					{
						read_success=get_a_fastq_read(infile,tag,seq_s,QS_s);
					
					}
					else
					{
						read_success=get_a_fasta_read(infile,tag,seq_s,n_tag);
			
					}	


					
					if(read_success==0)
					{break;}
						
						
				
					seq_sz=seq_s.size();
					if(TrimLen<seq_sz)
					{
						seq_s.resize(TrimLen);
						seq_sz=TrimLen;
					}
					
					if (seq_s.size()==0)
					{
						cout<<"Empty sequence!"<<endl;
						continue;
					}
					bool bad_flag=0;
					int numN=0;
					for(int i=0;i<seq_sz;++i)
					{
						if(seq_s[i]!='A'&&seq_s[i]!='C'&&seq_s[i]!='G'&&seq_s[i]!='T'&&seq_s[i]!='N')
						{
							bad_flag=1;
							break;
						}
						if(seq_s[i]=='N')
						{
							numN++;
						}
					}
					if(numN>TrimN)
					{bad_flag=1;}
					if(bad_flag)
					{continue;}
							

					bad_flag=0;
					
					numReads++;

					char QS_seq[1000];
					string QS_s;
					if(fq_flag)
					{
					
						//cout<<QS_s<<endl;
						if(QS_s[QS_s.size()-1]=='\n'||QS_s[QS_s.size()-1]=='\r')
						{
							QS_s.resize(QS_s.size()-1);
						}
						strcpy(QS_seq,QS_s.c_str());
						if(QS_s.size()>seq_sz)
						{
							QS_s.resize(seq_sz);
						}
					}


					if(totReads!=0&&numReads>totReads)
					{
						
						numReads=0;
						break;
						
					}

			
					if(fq_flag&&TrimQual>0)
					{
						int QS_sz=strlen(QS_seq);
						for(int i=0;i<QS_sz;++i)
						{
							if(QS_seq[i]<QS_base)
							{
								seq_s[i]='\0';
								seq_s.resize(i);
								seq_sz=i;
								break;
							}
						}
					}
				

						
					int nN=seq_sz-1,isN=-1;
					for(int i=0;i<seq_sz;++i)
					{
						
						if(seq_s[i]=='-'||seq_s[i]=='N')
						{
							if(i<=seq_sz/2)
							{
								isN=i;
								continue;
							}
							else
							{
								nN=i-1;
								break;
							}
						}
					}
					int s=0;
					if((nN-isN)<=seq_sz/2)
					{
						bad_flag=1;
					}
					
					if(bad_flag==1)
					{
						seq_s.clear();
						continue;
					}

					if(isN>=0)
					{
						for(int i=isN+1;i<=nN;++i)
						{
							seq_s[s]=seq_s[i];
							s++;
						}
						seq_s[s]='\0';
						seq_s.resize(s);
					}
					
				
					if(Compress)
					{
						string_shrinkage(seq_s);
						seq_sz=seq_s.size();
					}

					if(bad_flag==1||seq_sz<(K_size+1))
					{seq_s.clear();continue;}
		
					seq1=seq_s;
		
		
					Init_Read(seq1,read1);
					strcpy(read1.tag,tag1.c_str());
					seq1.clear();
					size_t b2=bucket_count;

					if(read1.readLen<K_size)
					{continue;}


					int success_dn=0;
					if(K_size<=32)
					{
						uint64_t OverlapKmers=read1.readLen-K_size+1;
				
						success_dn=Sparse_Denoising(&read1,&ht,&ht2,&CovTh,&CorrTh, K_size, gap,&correction_cnt,Hybrid);

					}
					else
					{
						uint64_t OverlapKmers=read1.readLen-K_size+1;				
						success_dn=Sparse_Denoising(&read1,&ht,&ht2,&CovTh,&CorrTh, K_size, gap,&correction_cnt,Hybrid);
	
					}
					bool start_print=0;
					if(success_dn)
					{
						read1.tag[0]='>';
						o_denoised<<read1.tag<<endl;
						for(int i=0;i<read1.readLen;++i)
						{

							if((read1.error_nt[i]==0&&start_print==0)||Hybrid)
							{
								start_print=1;
							}
							if(start_print&&read1.error_nt[i]&&(!Hybrid))
							{break;}

							if(start_print)
							{
								o_denoised<<read1.c_seq[i];	
		
							}

		
						}
						o_denoised<<endl;

				
					}
	//				numReads++;
	

				
					if (numReads%10000000==0)
					{
						cout<<"Denoised: "<<numReads<<endl;
			
						time(&read_time);
						cout<<"Time: "<<difftime(read_time,beg_time)<<" secs."<<endl;
				
					}
				
			

				}

				infile.close();
				infile.clear();
				
			}

			//free(read1.read_bits);
			//free(read2.read_bits);

			for(int jj=0;jj<p1_filenames_vt.size();jj+=1)
			{
			
				ifstream in_pair1(p1_filenames_vt[jj].c_str());
				ifstream in_pair2(p2_filenames_vt[jj].c_str());


				string t_filename;

		
				int s,t=0;
				for(s=(int) p1_filenames_vt[jj].size()-1;s>=0;--s)
				{
					if((p1_filenames_vt[jj][s]== '\\' )||( p1_filenames_vt[jj][s]=='/'))
					{
						t=s+1;
						break;
					}
				}
				t_filename=p1_filenames_vt[jj].substr(t,p1_filenames_vt[jj].size());
				string denoised_name1="Denoised_"+t_filename;
				string denoised_name3="Denoised_Single_"+t_filename;

			//	t_filename=p2_filenames_vt[jj];
				s,t=0;
				for(s=(int) p2_filenames_vt[jj].size()-1;s>=0;--s)
				{
					if((p2_filenames_vt[jj][s]== '\\' )||( p2_filenames_vt[jj][s]=='/'))
					{
						t=s+1;
						break;
					}
				}
				t_filename=p2_filenames_vt[jj].substr(t,p2_filenames_vt[jj].size());


				string denoised_name2="Denoised_"+t_filename;


				ofstream o_denoised1(denoised_name1.c_str()),o_denoised2(denoised_name2.c_str()),o_denoised3(denoised_name3.c_str());
		
				cout<<jj+1<<"/"<<p1_filenames_vt.size()<<" pairs."<<endl;
				cout<<"Processing file: "<<p1_filenames_vt[jj]<<" & "<<p2_filenames_vt[jj]<<endl;
				//while(getline(infile,seq)&&seq.size()!=0)

				bool tag_mismatch=0;
				char lastLineFC1='0',lastLineFC2='0';
		
				bool fq_flag=0;
				uint64_t nLines1=0,nLines2=0;
		

				string seq_s1,seq_s2,tag_s1,tag_s1n,tag_s2,tag_s2n,str1,str2;
				int seq1_sz,seq2_sz;
				string fq_tmp;
	
				struct read_t Read1,Read2;

				Read1.read_bits =(uint64_t*) malloc(MaxReadLen/4+100);
				Read2.read_bits =(uint64_t*) malloc(MaxReadLen/4+100);
			

				int readLen1,readLen2;
				nLines1=0,nLines2=0;









			
				getline(in_pair1,str);
				if(fq_flag==0&&str[0]=='@')
				{
					fq_flag=1;	
				}
				in_pair1.close();

				in_pair1.clear();
				in_pair1.open(in_filenames_vt[jj].c_str());

				bool read_success1=1,read_success2=1;


				string tag1,n_tag1,tag2,n_tag2;
				string QS_s1,QS_s2;

				while(read_success1&read_success2)
				{
					if(fq_flag)
					{
						read_success1=get_a_fastq_read(in_pair1,tag1,seq_s1,QS_s1);
						read_success2=get_a_fastq_read(in_pair2,tag2,seq_s2,QS_s2);

					
					}
					else
					{
						read_success1=get_a_fasta_read(in_pair1,tag1,seq_s1,n_tag1);
						read_success1=get_a_fasta_read(in_pair2,tag2,seq_s2,n_tag2);
			
					}	

					bool bad_flag1=0,bad_flag2=0;
					bool success_dn1=0,success_dn2=0;
			
				

				
					seq1_sz=seq_s1.size();
					if(TrimLen<seq1_sz)
					{
						seq_s1.resize(TrimLen);
						seq1_sz=TrimLen;
					}
					
					if (seq_s1.size()==0)
					{
						cout<<"Empty sequence!"<<endl;
						continue;
					}
					int numN=0;
					for(int i=0;i<seq1_sz;++i)
					{
						if(seq_s1[i]!='A'&&seq_s1[i]!='C'&&seq_s1[i]!='G'&&seq_s1[i]!='T'&&seq_s1[i]!='N')
						{
							bad_flag1=1;
							break;
						}
						if(seq_s1[i]=='N')
						{
							numN++;
						}
					}
					if(numN>TrimN)
					{bad_flag1=1;}
			//		if(bad_flag1)
		//			{continue;}
							
		
					numReads++;
					if(1)
					{
						char QS_seq1[1000];
						if(fq_flag)
						{

							//cout<<QS_s<<endl;
							if(QS_s1[QS_s1.size()-1]=='\n'||QS_s1[QS_s1.size()-1]=='\r')
							{
								QS_s1.resize(QS_s1.size()-1);
							}
							strcpy(QS_seq1,QS_s1.c_str());
							if(QS_s1.size()>seq1_sz)
							{
								QS_s1.resize(seq1_sz);
							}
						}



						if(totReads!=0&&numReads>totReads)
						{
						
							numReads=0;
							break;
						
						}

						if(fq_flag&&TrimQual>0)
						{
							int QS_sz=strlen(QS_seq1);
							for(int i=0;i<QS_sz;++i)
							{
								if(QS_seq1[i]<QS_base)
								{
									seq_s1[i]='\0';
									seq_s1.resize(i);
									seq1_sz=i;
									break;
								}
							}
						}
				

						
						int nN=seq1_sz-1,isN=-1;
						for(int i=0;i<seq1_sz;++i)
						{
						
							if(seq_s1[i]=='-'||seq_s1[i]=='N')
							{
								if(i<=seq1_sz/2)
								{
									isN=i;
									continue;
								}
								else
								{
									nN=i-1;
									break;
								}
							}
						}
						int s=0;
						if((nN-isN)<=seq1_sz/2)
						{
							bad_flag1=1;
						}
					


						if(!bad_flag1&&isN>=0)
						{
							for(int i=isN+1;i<=nN;++i)
							{
								seq_s1[s]=seq_s1[i];
								s++;
							}
							seq_s1[s]='\0';
							seq_s1.resize(s);
						}
					
				

					}
				




					seq2_sz=seq_s2.size();
					if(TrimLen<seq2_sz)
					{
						seq_s2.resize(TrimLen);
						seq2_sz=TrimLen;
					}
					
					if (seq_s2.size()==0)
					{
						cout<<"Empty sequence!"<<endl;
						continue;
					}
					numN=0;
					for(int i=0;i<seq2_sz;++i)
					{
						if(seq_s2[i]!='A'&&seq_s2[i]!='C'&&seq_s2[i]!='G'&&seq_s2[i]!='T'&&seq_s2[i]!='N')
						{
							bad_flag2=1;
							break;
						}
						if(seq_s2[i]=='N')
						{
							numN++;
						}
					}
					if(numN>TrimN)
					{bad_flag2=1;}	
		
					if(1)
					{
						char QS_seq2[1000];
						string QS_s2;
						if(fq_flag)
						{

							if(QS_s2[QS_s2.size()-1]=='\n'||QS_s2[QS_s2.size()-1]=='\r')
							{
								QS_s2.resize(QS_s2.size()-1);
							}
							strcpy(QS_seq2,QS_s2.c_str());
							if(QS_s2.size()>seq2_sz)
							{
								QS_s2.resize(seq2_sz);
							}
						}


						if(totReads!=0&&numReads>totReads)
						{
						
							numReads=0;
							break;
						
						}

						if(fq_flag&&TrimQual>0)
						{
							int QS_sz=strlen(QS_seq2);
							for(int i=0;i<QS_sz;++i)
							{
								if(QS_seq2[i]<QS_base)
								{
									seq_s2[i]='\0';
									seq_s2.resize(i);
									seq2_sz=i;
									break;
								}
							}
						}
				

						
						int nN=seq2_sz-1,isN=-1;
						for(int i=0;i<seq2_sz;++i)
						{
						
							if(seq_s2[i]=='-'||seq_s2[i]=='N')
							{
								if(i<=seq2_sz/2)
								{
									isN=i;
									continue;
								}
								else
								{
									nN=i-1;
									break;
								}
							}
						}
						int s=0;
						if((nN-isN)<=seq2_sz/2)
						{
							bad_flag2=1;
						}
					


						if(!bad_flag2&&isN>=0)
						{
							for(int i=isN+1;i<=nN;++i)
							{
								seq_s2[s]=seq_s2[i];
								s++;
							}
							seq_s2[s]='\0';
							seq_s2.resize(s);
						}
					
				

					}
					if(Compress)
					{
						string_shrinkage(seq_s1);
					
						string_shrinkage(seq_s2);
					
					}
					readLen1=seq_s1.size();
					readLen2=seq_s2.size();

					if(readLen1<K_size)
					{bad_flag1=1;}
					if(readLen2<K_size)
					{bad_flag2=1;}

					if(bad_flag1==0)
					{
					
						Init_Read(seq_s1,read1);
						strcpy(read1.tag,tag_s1.c_str());
						seq_s1.clear();
						size_t b2=bucket_count;
						if(K_size<=32)
						{
							uint64_t OverlapKmers=read1.readLen-K_size+1;
				
							success_dn1=Sparse_Denoising(&read1,&ht,&ht2,&CovTh,&CorrTh, K_size, gap,&correction_cnt,Hybrid);


						}
						else
						{
							uint64_t OverlapKmers=read1.readLen-K_size+1;
				
							success_dn1=Sparse_Denoising(&read1,&ht,&ht2,&CovTh,&CorrTh, K_size, gap,&correction_cnt,Hybrid);
					
						}
					}
		
			
					if(bad_flag2==0)
					{
			
			
					
						Init_Read(seq_s2,read2);
						strcpy(read2.tag,tag_s2.c_str());
						seq_s2.clear();
						size_t b2=bucket_count;

				
						if(K_size<=32)
						{
							uint64_t OverlapKmers=read2.readLen-K_size+1;
				
							success_dn2=Sparse_Denoising(&read2,&ht,&ht2,&CovTh,&CorrTh, K_size, gap,&correction_cnt,Hybrid);

						}
						else
						{
							uint64_t OverlapKmers=read2.readLen-K_size+1;
				
							success_dn2=Sparse_Denoising(&read2,&ht,&ht2,&CovTh,&CorrTh, K_size, gap,&correction_cnt,Hybrid);

					
						}
					}
		



					bool start_print=0;
					if(success_dn1&&success_dn2)
					{
						read1.tag[0]='>';
						o_denoised1<<read1.tag<<endl;
				
						for(int i=0;i<read1.readLen;++i)
						{
							if((read1.error_nt[i]==0&&start_print==0)||Hybrid)
							{
								start_print=1;
							}
							if(start_print&&read1.error_nt[i]&&(!Hybrid))
							{break;}

							if(start_print)
							{
								o_denoised1<<read1.c_seq[i];	
		
							}

		
						}
						o_denoised1<<endl;
						start_print=0;
						read2.tag[0]='>';
						o_denoised2<<read2.tag<<endl;
						for(int i=0;i<read2.readLen;++i)
						{
							if((read2.error_nt[i]==0&&start_print==0)||Hybrid)
							{
								start_print=1;
							}
							if(start_print&&read2.error_nt[i]&&(!Hybrid))
							{break;}

							if(start_print)
							{
								o_denoised2<<read2.c_seq[i];	
		
							}

		
						}
						o_denoised2<<endl;




					}
					else
					{
						if(success_dn1==1)
						{
							read1.tag[0]='>';
							o_denoised3<<read1.tag<<endl;
							for(int i=0;i<read1.readLen;++i)
							{
								if((read1.error_nt[i]==0&&start_print==0)||Hybrid)
								{
									start_print=1;
								}
								if(start_print&&read1.error_nt[i]&&(!Hybrid))
								{break;}

								if(start_print)
								{
									o_denoised3<<read1.c_seq[i];	
		
								}

		
							}
							o_denoised3<<endl;
						}
						start_print=0;
						if(success_dn2==1)
						{
							read2.tag[0]='>';
							o_denoised3<<read2.tag<<endl;
							for(int i=0;i<read2.readLen;++i)
							{
								if((read2.error_nt[i]==0&&start_print==0)||Hybrid)
								{
									start_print=1;
								}
								if(start_print&&read2.error_nt[i]&&(!Hybrid))
								{break;}

								if(start_print)
								{
									o_denoised3<<read2.c_seq[i];	
		
								}

		
							}
							o_denoised3<<endl;
						}

				
					}
			
				
					//numReads++;
				
				//	if (numReads==180)
				//	{cout<<numReads<<endl;}
				
					if (numReads%10000000==0)
					{
						cout<<"Denoised: "<<numReads<< " pairs."<<endl;
			
						time(&read_time);
						cout<<"Time: "<<difftime(read_time,beg_time)<<" secs."<<endl;
				
					}
				
				
			

				}

				//free(Read1.read_bits);
				//free(Read2.read_bits);
	
				in_pair1.close();
				in_pair2.close();
				in_pair1.clear();
				in_pair2.clear();
			}


	
			time(&read_time);

			cout<<"Time to now: "<<difftime(read_time,beg_time)<<" secs."<<endl;
	
	
			bool OUT_LOG=1;
			if(OUT_LOG==1)
			{
				ofstream o_DnLog("DenoisingLog.txt",ios::app);
				o_DnLog<<"Total corrections:"<<endl<<correction_cnt<<endl;;
			}
	


		}



		if(!Denoise)
		{
			//assembly

			//call bubble removal
			if(BFS&&(!Iter_Scaffold))
			{
				cout<<"Graph simplification..."<<endl;
				size_t TipsRemoved=0;
				int TipLenTh=100;


			
				string Contig_Filename="Contigs.txt";

				if(MemoryEfficient&(K_size<=32))
				{
					Init_HT(&merge_ht,hashTableSZ/50);
				}
				else
				{
					Init_HT0(&merge_ht0,hashTableSZ/50);	
					/*
					if(K_size>32&&K_size<=64)
					{
						Init_HT2(&merge_ht2,hashTableSZ/50);	
					}
					else
					{
						if(K_size<=96)
						{
							Init_HT3(&merge_ht3,hashTableSZ/50);	
						}
						else
						{
							if(K_size<=128)
							{
								Init_HT4(&merge_ht4,hashTableSZ/50);	
							}
						}
					}
					*/
				}




				MaxDepth=max(MaxDepth,3*K_size/gap);
				int MinDepth=max(3,(1*K_size/gap)+1);//MaxDepth;// 
				int jmp=max((MaxDepth-MinDepth)/3,1);
				cout<<"Min searching depth: "<<MinDepth<<endl;
				cout<<"Max searching depth: "<<MaxDepth<<endl;

				//increase search depths in each iteration

				for(int depth=MinDepth;depth<=MaxDepth;depth+=jmp)
				{
					cout<<"Current depth: "<<depth<<endl;
			
			
					if(1)
					{
						//use the below:
						if(MemoryEfficient&(K_size<=32))
						{
							GraphSimplification(&ht,&merge_ht,&ht2, &merge_ht2,K_size, gap,PathCovTh,depth,PathSim);
							
						}
						else
						{
							GraphSimplification0(&ht0,&merge_ht0,K_size, gap,PathCovTh,depth,PathSim);
							
							/*
							if(K_size<=96)
							{
								GraphSimplification3(&ht3,&merge_ht3, K_size, gap,PathCovTh,depth,PathSim);
							}
							else
							{
								if(K_size<=128)
								{
									GraphSimplification4(&ht4,&merge_ht4, K_size, gap,PathCovTh,depth,PathSim);	
								}
							}
							*/
						}
				
						if(MemoryEfficient&(K_size<=32))
						{
							
							RemovingWeakNodesAndEdges(&ht, K_size,0, 0,&bucket_count, &edge_cnt);
							
						}
						else
						{
							RemovingWeakNodesAndEdges0(&ht0, K_size,0, 0,&bucket_count, &edge_cnt);
							/*
							if(K_size>32&&K_size<=64)
							{
								RemovingWeakNodesAndEdges2(&ht2, K_size,0, 0,&bucket_count, &edge_cnt);	
							}
							else
							{
								if(K_size<=96)
								{
									RemovingWeakNodesAndEdges3(&ht3, K_size,0, 0,&bucket_count, &edge_cnt);	
								}
								else
								{
									if(K_size<=128)
									{
										RemovingWeakNodesAndEdges4(&ht4, K_size,0, 0,&bucket_count, &edge_cnt);	
									}
								}
							}
							*/
						}

				
					}
		
				}
	
				if(MemoryEfficient&(K_size<=32))
				{
			
					SavingMergeHT(&merge_ht);
			
				}
				else
				{
					SavingMergeHT0(&merge_ht0,Kmer_arr_sz);
				}
				/*
				if(K_size>32&&K_size<=64)
				{
					SavingMergeHT2(&merge_ht2);
				}
		
				if(K_size>64&&K_size<=96)
				{
					SavingMergeHT3(&merge_ht3);
				}
		
		
				if(K_size>96&&K_size<=128)
				{
					SavingMergeHT4(&merge_ht4);
				}
				*/
			}

			//mark the branches in the simplified graph

	
			if(MemoryEfficient&(K_size<=32))
			{
				MarkBranches(&ht);
			
			}
			else
			{
				MarkBranches0(&ht0);
				/*
				if(K_size<=64)
				{
					MarkBranches2(&ht2);
				}
				else
				{
					if(K_size<=96)
					{
						MarkBranches3(&ht3);
					}
					else
					{
						MarkBranches4(&ht4);
					}
				}
				*/
			}

		//	struct contigs_info  contigs_info;
			string Contig_Filename="Contigs.txt";
			Contig_Filename=ContigPrefix+Contig_Filename;
			// build contigs, single end assembly complete.
			if(BUILD_CONTIGS&&!Iter_Scaffold)
			{
				cout<<"Build contigs..."<<endl;
				bool ScreenOffTips=0;
				if(MemoryEfficient&(K_size<=32))
				{
					
					build_contigs(&ht,K_size, gap,Contig_Filename,ScreenOffTips);
					
				}
				else
				{

					build_contigs0(&ht0,&key_table,K_size, gap,Contig_Filename,ScreenOffTips);
					/*
					if(K_size<=64)
					{
						build_contigs2(&ht2,K_size, gap,Contig_Filename,ScreenOffTips);
					}
					else
					{
						if(K_size<=96)
						{
							build_contigs3(&ht3,K_size, gap,Contig_Filename,ScreenOffTips);
						}
						else
						{
							if(K_size<=128)
							{
								build_contigs4(&ht4,K_size, gap,Contig_Filename,ScreenOffTips);
							}
						}
					}
					*/
				}
				/*
				ifstream in_long_contigs(Contig_Filename.c_str());

				string tag,seq_s,n_tag;
				vector<size_t> contig_len_vt;
				size_t cont_len;
				bool contig_success=1;
				
				while(contig_success)
				{
					contig_success=get_a_fasta_read(in_long_contigs,tag,seq_s,n_tag);
					
					cont_len=seq_s.size();

					contig_len_vt.push_back(cont_len);
				}



				cout<<"Sorting..."<<endl;
				sort(contig_len_vt.begin(),contig_len_vt.end());

				uint64_t sum=0,sum1=0;
				for (size_t i=0;i<contig_len_vt.size();++i)
				{
					if(contig_len_vt[i]>100)
					{
						sum+=contig_len_vt[i];
					}
				}

				string o_file=Contig_Filename+"_Stats.txt";
				ofstream o_stats(o_file.c_str());
				sum1=0;
				uint64_t sum2=0;
				double r=0.5;
				o_stats<<"Total len: "<<sum<<endl;
				TotLength_lst.push_back(sum);
				
				o_stats<<"Longest contig: "<<contig_len_vt[contig_len_vt.size()-1]<<endl;
				MaxLen_lst.push_back(contig_len_vt[contig_len_vt.size()-1]);
				if(GenomeSize>0)
				{
					sum=GenomeSize;
				}
				for (int i=contig_len_vt.size()-1;i>=1;--i)
				{
					sum1+=contig_len_vt[i];


					if(sum1>double(sum)*r)
					{
						o_stats<<"N50: "<<contig_len_vt[i]<<endl;
						N50_lst.push_back(contig_len_vt[i]);
						break;
					}		
		
				}
				*/
 				if(K_size<=32)
				{
					ConstructContigGraph(&ht,&merge_ht, K_size, &contigs_info,Contig_Filename);
				}
				else
				{
					ConstructContigGraph0(&ht0,&merge_ht0, K_size, &contigs_info,Contig_Filename);
				}

			}

			
			if(p_filenames_vt.size()==0&&DBG2OLC)
			{

				cout<<"Collecting non-contained reads"<<endl;
				
				if(MemoryEfficient&(K_size<=32))
				{
					
					CollectingNonContainedReadsSlow(&ht,&merge_ht, &ht2, &merge_ht2,K_size, in_filenames_vt,&contigs_info,Contig_Filename);
					
				}
				else
				{

					
					CollectingNonContainedReadsSlow0(&ht0,&merge_ht0,K_size, in_filenames_vt,&contigs_info,Contig_Filename);
					
				}

			}
			else
			{
				if (DBG2OLC)
				{
					cout << "Collecting non-contained pairs" << endl;
					if (MemoryEfficient&(K_size <= 32))
					{

						CollectingNonContainedPairsSlow(&ht, &merge_ht, &ht2, &merge_ht2, K_size, p_filenames_vt, &contigs_info, Contig_Filename);

					}
					else
					{


						CollectingNonContainedPairsSlow0(&ht0, &merge_ht0, K_size, in_filenames_vt, &contigs_info, Contig_Filename);

					}

				}
				
			}
			

			if(Resolving_Branches_PE&&(!Iter_Scaffold))
			{
				sp_filenames_vt=p_filenames_vt;
				SingleOutwardLib=OutwardLib;

				for(int f=0;f<in_filenames_vt.size();++f)
				{
					sp_filenames_vt.push_back(in_filenames_vt[f]);
					sp_filenames_vt.push_back(in_filenames_vt[f]);
					SingleOutwardLib.push_back(0);
				}

				p_filenames_vt=sp_filenames_vt;
				OutwardLib=SingleOutwardLib;


				string ContigFilename="Contigs.txt";
				if(LOAD_PE_DIST==0)
				{
					bool MatePair=0;
					if(K_size<=64)
					{
					
						ContigGapEst(&ht,&merge_ht, &ht2,&merge_ht2, K_size,insert_sz_vt, p_filenames_vt,OutwardLib,&contigs_info,ContigFilename,ResumePE,totReads,MatePair);
					}
					else
					{
						if(K_size<=96)
						{
							ContigGapEst3(&ht3,&merge_ht3, K_size,insert_sz_vt, p_filenames_vt,OutwardLib,&contigs_info,ContigFilename,ResumePE,totReads,MatePair);
						}
						else
						{
							if(K_size<=128)
							{
								ContigGapEst4(&ht4,&merge_ht4, K_size,insert_sz_vt, p_filenames_vt,OutwardLib,&contigs_info,ContigFilename,ResumePE,totReads,MatePair);
							}
						}
					}
				}
				else
				{
					if(K_size<=64)
					{
						ContigsRemapping(&ht,&ht2, K_size, &contigs_info,ContigFilename,0);
						RemoveUnmappedNodes(&ht,&ht2, K_size);
						BuildContigAdjacency(&ht, &ht2, &contigs_info, K_size,ContigFilename);

					}
					else
					{
						if(K_size<=96)
						{
							ContigsRemapping3(&ht3, K_size, &contigs_info,ContigFilename,0);
							RemoveUnmappedNodes3(&ht3, K_size);
							BuildContigAdjacency3(&ht3, &contigs_info, K_size,ContigFilename);	
				
						}
						else
						{
							if(K_size<=128)
							{
								ContigsRemapping4(&ht4, K_size, &contigs_info,ContigFilename,0);
								RemoveUnmappedNodes4(&ht4, K_size);
								BuildContigAdjacency4(&ht4, &contigs_info, K_size,ContigFilename);		
					
							}
						}
					}
				}
		
				//FindConnectedComponents(&contigs_info);

				if(mp1_filenames_vt.size()==0)
				{
					if(K_size<=32)
					{
						FreeSparseKmerGraph(&ht);
					}
					if(K_size>32&&K_size<=64)
					{
						FreeSparseKmerGraph2(&ht2);
					}
					if(K_size>64&&K_size<=96)
					{
						FreeSparseKmerGraph3(&ht3);
					}
					if(K_size>96&&K_size<=128)
					{
						FreeSparseKmerGraph4(&ht4);
					}
				}
			
				if(ExpCov==0)
				{cout<<"Error! Expected contig coverage (ExpCov) is not given."<<endl;return -1;}
				ResolvingRepeatsPE( insert_sz_vt, p_filenames_vt, &contigs_info,ContigFilename,LinkCovTh,UniqueLenTh,ExpCov);

			}

			//below is for mate pair scaffolding, omit first:

			if((MP_Scaffold&&(mp1_filenames_vt.size()>0))||Iter_Scaffold)
			{
				bool MatePair=1;
				string ContigFilename="SuperContigs.txt";			
				if(mp1_filenames_vt.size()==0)
				{
					MatePair=0;
					mp_filenames_vt=p_filenames_vt;
					OutwardLibMP=OutwardLib;
				}

				if(LOAD_PE_DIST==0)
				{
				
					if(K_size<=64)
					{					
						ContigGapEst(&ht,&merge_ht, &ht2,&merge_ht2, K_size,insert_sz_vt, mp_filenames_vt,OutwardLibMP,&scaffolds_info,ContigFilename,ResumePE,totReads,MatePair);
					}
					else
					{
						if(K_size<=96)
						{
							ContigGapEst3(&ht3,&merge_ht3, K_size,insert_sz_vt, mp_filenames_vt,OutwardLibMP,&scaffolds_info,ContigFilename,ResumePE,totReads,MatePair);
						}
						else
						{
							if(K_size<=128)
							{
								ContigGapEst4(&ht4,&merge_ht4, K_size,insert_sz_vt, mp_filenames_vt,OutwardLibMP,&scaffolds_info,ContigFilename,ResumePE,totReads,MatePair);
							}
						}
					}
				}
				else
				{
					if(K_size<=64)
					{
					


						if(K_size<=32)
						{
						//	MarkUniqueKmers(&ht,K_size);
						}
						else
						{
						//	MarkUniqueKmers2(&ht2,K_size);
						}
						SuperContigsRemapping(&ht,&ht2, K_size, &scaffolds_info,ContigFilename,0);
						RemoveUnmappedNodes(&ht,&ht2, K_size);
						BuildContigAdjacency(&ht, &ht2, &scaffolds_info, K_size,ContigFilename);
					}
					else
					{
						if(K_size<=96)
						{
						//	MarkUniqueKmers3(&ht3,K_size);
							SuperContigsRemapping3(&ht3, K_size, &scaffolds_info,ContigFilename,0);
							RemoveUnmappedNodes3(&ht3, K_size);
							BuildContigAdjacency3(&ht3, &scaffolds_info, K_size,ContigFilename);	
				
						}
						else
						{
							if(K_size<=128)
							{
							//	MarkUniqueKmers4(&ht4,K_size);
								SuperContigsRemapping4(&ht4, K_size, &scaffolds_info,ContigFilename,0);
								RemoveUnmappedNodes4(&ht4, K_size);
								BuildContigAdjacency4(&ht4, &scaffolds_info, K_size,ContigFilename);		
					
							}
						}
					}
				}
				if(scaffolds_info.cov_vt.size()<scaffolds_info.total_contigs)
				{
					scaffolds_info.cov_vt.resize(scaffolds_info.total_contigs);
				}
		
				//FindConnectedComponents(&scaffolds_info);
		
				if(K_size<=32)
				{
					FreeSparseKmerGraph(&ht);
				}
				if(K_size>32&&K_size<=64)
				{
					FreeSparseKmerGraph2(&ht2);
				}
				if(K_size>64&&K_size<=96)
				{
					FreeSparseKmerGraph3(&ht3);
				}
				if(K_size>96&&K_size<=128)
				{
					FreeSparseKmerGraph4(&ht4);
				}

				if(ExpCov==0)
				{cout<<"Error! Expected contig coverage (ExpCov) is not given."<<endl;return -1;}
				ResolvingRepeatsPE( insert_sz_vt, mp_filenames_vt, &scaffolds_info,ContigFilename,LinkCovTh,UniqueLenTh,ExpCov);

			}


		}

		if(bfilename.empty())
		{break;}

	}

	if(N50_lst.size()>1)
	{
		ofstream o_N50,o_max_ctg,o_tot_len;
		o_N50.open("N50_lst.txt");
		o_max_ctg.open("MaxLen_lst.txt");
		o_tot_len.open("TotLen_lst.txt");
		for(int i=0;i<N50_lst.size();++i)
		{
			o_N50<<N50_lst[i]<<endl;
			o_max_ctg<<MaxLen_lst[i]<<endl;
			o_tot_len<<TotLength_lst[i]<<endl;
		}

	}
	return 0;
}

