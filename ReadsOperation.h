#ifndef __READS_OPERATION_H
#define __READS_OPERATION_H

#include <iostream>
#include <string>
#include <string.h>
#include <stdint.h>
#include <vector>
#include <map>
#include <list>
#include <algorithm>
#include <fstream>
#include "time.h"
#include "BasicDataStructure.h"
#include "GraphConstruction.h"
#include "GraphSimplification.h"
#include "ScaffoldingDataStructure.h"

using namespace std;


struct SuperRead_t
{
	int PathsCnt;
	uint64_t idx;
	int PathsLength;
	int WellID;
	char seq[10000];
	bool PathFound;
	uint16_t depth;
	int ListSize;
	string extension;

};


struct SuperRead_ctg
{
	int PathsCnt;
	uint64_t idx;
	int PathsLength;

	bool PathFound;
	uint16_t depth;
	vector< vector<int> > Paths;

};


struct search_info
{
	int pos1,pos2;
	bool RightSearch;
	bool Flip_End;
	int offset;
};


struct BFS_path_info_v2
{
	int cov;
	int depth;
	int len;
	bool BothEndsUsed;
	uint64_t last_kmer;
	struct edge_node* last_bkt_edge;
};
struct BFS_path_info_v3
{
	int cov;
	int depth;
	int len;
	//bool BothEndsUsed;
	int RepeatVisits;
	int last_ctg;
	
};

struct reads_overlap_info
{
	vector< map<int32_t, vector<int32_t> > > left_overlaps,right_overlaps;
	vector<int> cov_vt;
	vector<bool> contained_vt,used_vt;

};

struct KmerInContig
{
	uint32_t contig_no : 31, flip : 1;
	int pos;
};
struct LongReadContigIndex
{
	map<int, KmerInContig> LR2CTG;
	map<int, vector<int> > CTG2LR, CTG2LR_2;
	map<int, vector<int> > layout;
	int KmerCovTh;
	int nMatches;
	//	vector<int> matching_positions_LR;
	//vector<KmerInContig> matching_positions_ctg;
};



void SavingReadsTable(reads_table *reads_table)
{
	ofstream  o_reads_table,o_reads_table_info;
	string reads_table_name="ReadsTable_content",reads_table_info_name="ReadsTable_info.txt";
	o_reads_table.open(reads_table_name.c_str(),ios_base::out|ios_base::binary);
	o_reads_table_info.open(reads_table_info_name.c_str(),ios_base::out|ios_base::binary);
	o_reads_table_info<<reads_table->BytesPerBlock;
	o_reads_table_info<<reads_table->pblocks.size();
	o_reads_table_info<<reads_table->read_len_vt.size();
	
	while(reads_table->pblocks.size()>0)
	{
		uint64_t *block_ptr=(reads_table->pblocks).front();
		(reads_table->pblocks).pop_front();
		o_reads_table.write((char*) block_ptr,reads_table->BytesPerBlock);

	}
	for(int i=0;i<reads_table->read_len_vt.size();++i)
	{
		o_reads_table_info<<reads_table->read_len_vt[i];
	}
	

}

void LoadingReadsTable(reads_table *reads_table)
{
	ifstream  in_reads_table,in_reads_table_info;
	string reads_table_name="ReadsTable_content",reads_table_info_name="ReadsTable_info.txt";
	in_reads_table.open(reads_table_name.c_str(),ios_base::in|ios_base::binary);
	in_reads_table_info.open(reads_table_info_name.c_str(),ios_base::in|ios_base::binary);
	in_reads_table_info>>reads_table->BytesPerBlock;
	int total_blocks;
	in_reads_table_info>>total_blocks;
	int n_reads;
	in_reads_table_info>>n_reads;
	reads_table->read_len_vt.resize(n_reads);
	
	while(reads_table->pblocks.size()>0)
	{
		uint64_t *block_ptr=(uint64_t*) malloc(reads_table->BytesPerBlock);
		in_reads_table.read((char*) (*block_ptr),sizeof(reads_table->BytesPerBlock));
		(reads_table->pblocks).push_back(block_ptr);
	}
	reads_table->current_byte=0;
	if(reads_table->pblocks.size()==0)
	{
		cout<<"Loading error!"<<endl;
		return;
	}
	list<uint64_t*>::iterator new_block_itr=reads_table->pblocks.begin();
	reads_table->current_read=0;
	reads_table->current_byte=0;
	for(int i=0;i<n_reads;++i)
	{
		in_reads_table_info>>reads_table->read_len_vt[i];

		int Read_arr_sz=reads_table->read_len_vt[i]/32+1;
		int rem=reads_table->read_len_vt[i]%32;
		if(rem==0)
		{Read_arr_sz--;}

		if((reads_table->current_byte+Read_arr_sz*8)>=reads_table->BytesPerBlock)
		{
			new_block_itr++;
			reads_table->current_byte=0;
		}
		uint64_t *block_ptr=*new_block_itr;
		reads_table->current_read++;
		
		reads_table->read_ptr[reads_table->current_read]=block_ptr+reads_table->current_byte/8;
		reads_table->current_byte+=Read_arr_sz*8;
		
	}
	
}


void CollectingNonContainedReads(struct hashtable *ht1,struct hashtable *merge_ht1, struct hashtable2 *ht2, struct hashtable2 *merge_ht2,int K_size,vector<string>& filenames_vt, contigs_info * contigs_info,string ContigFilename)
{

	time_t beg_time,read_time;
	string in_fname=ContigFilename;
	ifstream Contigs_in(in_fname.c_str());
	size_t num_Contigs=0;
	string tag,s,kmer,str,seq_s;
	uint64_t f_seq,hv;
	struct kmer_t2 f_seq_t2;
	size_t hash_idx;
	bool found;
	bool flip_1,flip_2,flip_0;
	size_t ht_sz;
	size_t numReads=0;
	ofstream out_reads("reads_temp.txt");
	
	int boundary=0,removed=0,bridge=0;
	if(K_size<=32)
	{
		ht_sz=ht1->ht_sz;
	}
	else
	{
		ht_sz=ht2->ht_sz;
	}
	cout<<"Contigs remapping."<<endl;
	if(ContigFilename=="Contigs.txt")
	{
		ContigsRemapping(ht1,ht2, K_size, contigs_info,ContigFilename,0);
		
	}
	if(K_size<=32)
	{
		 AppendMergeHT(ht1, merge_ht1);
	}
	else
	{
		 AppendMergeHT2(ht2,merge_ht2);
	
	}

	//cout<<"Collecting informative reads."<<endl;
	time(&beg_time);

	int64_t dist_sum=0,p_cnt=0;
	int mean_dist=0;
	int lib_no=0;
	
	for(size_t ii=0;ii<filenames_vt.size();ii++)
	{
		lib_no++;
	
		dist_sum=0,p_cnt=0;
		cout<<"Processing library: "<<ii<<endl;
		ifstream in_reads;
		in_reads.open(filenames_vt[ii].c_str());
		
		struct read_t Read1,Read2;

		Read1.read_bits =(uint64_t*) malloc(1000000/4+100);
		Read2.read_bits =(uint64_t*) malloc(1000000/4+100);

		int nLines1=0,nLines2=0;


		bool fq_flag=0;


		getline(in_reads,str);
		if(fq_flag==0&&str[0]=='@')
		{
			fq_flag=1;	
		}
		in_reads.close();

		in_reads.clear();
		in_reads.open(filenames_vt[ii].c_str());

		bool read_success=0;

		read_success=1;

		string tag,qs,n_tag;
		string QS_s;

		while(read_success)
		{
			if(fq_flag)
			{
				read_success=get_a_fastq_read(in_reads,tag,seq_s,QS_s);
					
			}
			else
			{
				read_success=get_a_fasta_read(in_reads,tag,seq_s,n_tag);
			
			}	
			if(read_success==0)
			{break;}



				
			int seq_sz=seq_s.size();
			
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
			
			if(bad_flag)
			{continue;}
							

			bad_flag=0;
					
			numReads++;


						
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
					


			Init_Read(seq_s,Read1);
			string read_str=seq_s;
			seq_s.clear();
			

			uint64_t bits1;
			int contig_no=-1;
			bool output_current_read=0;

			for(int i=0;i<Read1.readLen-K_size+1;++i )
			{
				if(K_size<=32)
				{
					get_sub_arr(Read1.read_bits,Read1.readLen,i,K_size,&bits1);
					f_seq=get_rev_comp_seq(bits1,K_size);
					flip_0=0;
					if(bits1>f_seq)
					{
						bits1=f_seq;
						flip_0=1;
					}

					hv=MurmurHash64A(&bits1,sizeof(bits1),0);
					hash_idx=(size_t) (hv%ht_sz);
					struct bucket **ptr1;
					ptr1= &(ht1->store_pos[hash_idx]);
					found=look_up_in_a_list(bits1,&ptr1);
					if(found)
					{
	
						if((*ptr1)->kmer_info.removed==0&&(*ptr1)->kmer_info.contig_no>0)
						{
							/*

							int cont1=(*ptr1)->kmer_info.contig_no;
							if(contig_no<0&&cont1>0)
							{
								contig_no=cont1;
							}
							else
							{
								
								if(cont1!=contig_no)
								{
									output_current_read=1;
									bridge++;
								}
								if((*ptr1)->kmer_info.masked==1)
								{
									output_current_read=1;
									boundary++;
								}


								////////////////// non-contained
							
								

								break;
								
							}
							*/
							// this small number of 50 gives better results.
								output_current_read=1;
								if(flip_0^(*ptr1)->kmer_info.flip)
								{
									if(i+50<=(*ptr1)->kmer_info.cod     &&              (Read1.readLen-i+50)<=(contigs_info->contig_sz_vt[(*ptr1)->kmer_info.contig_no]-(*ptr1)->kmer_info.cod))
									{
										output_current_read=0;
										//contained...
									}
								
								}
								else
								{
									if( (Read1.readLen-(i+K_size)+50)<=(*ptr1)->kmer_info.cod   &&   (i+K_size+50)  <=(contigs_info->contig_sz_vt[(*ptr1)->kmer_info.contig_no]-(*ptr1)->kmer_info.cod))
									{
										output_current_read=0;
										//contained...
									}
								}
								//if(output_current_read1==0&&output_current_read==1)
								//{cout<<"";}
								break;
								
						}
						else
						{

							/*
							output_current_read=1;
							removed++;
							break;

							*/


						}

						
					}
				}
				else
				{
					if(K_size>32&&K_size<=64)
					{
						kmer_t2 bits1_t2;
						get_sub_arr(Read1.read_bits,Read1.readLen,i,K_size,bits1_t2.kmer);

						f_seq_t2=bits1_t2;
						get_rev_comp_seq_arr(f_seq_t2.kmer,K_size,2);
						flip_0=0;
						if(uint64_t_cmp(bits1_t2.kmer,f_seq_t2.kmer,2)>0)
						{
							bits1_t2=f_seq_t2;
							flip_0=1;
						}

						hv=MurmurHash64A(bits1_t2.kmer,sizeof(bits1_t2),0);

						hash_idx=(size_t) (hv%ht_sz);
						struct bucket2 **ptr1;
						ptr1= &(ht2->store_pos[hash_idx]);
						found=look_up_in_a_list2(&bits1_t2,&ptr1);
						if(found)
						{
							
							if((*ptr1)->kmer_info.removed==0)
							{
							
								output_current_read=1;
								if(flip_0^(*ptr1)->kmer_info.flip)
								{
									if(i<=(*ptr1)->kmer_info.cod     &&         (Read1.readLen-i+(*ptr1)->kmer_info.cod)<=(contigs_info->contig_sz_vt[(*ptr1)->kmer_info.contig_no]))
									{
										output_current_read=0;
										//contained...
									}
								
								}
								else
								{
									if( (Read1.readLen-(i+K_size))<=(*ptr1)->kmer_info.cod   &&   (i+K_size+(*ptr1)->kmer_info.cod)  <=(contigs_info->contig_sz_vt[(*ptr1)->kmer_info.contig_no]))
									{
										output_current_read=0;
										//contained...
									}
								}
								//if(output_current_read1==0&&output_current_read==1)
								//{cout<<"";}
								break;
							}
							else
							{
								output_current_read=1;
								break;

							}

						}
					}
				}

			}

			if(output_current_read)
			{
				if(tag[0]=='@')
				{tag[0]='>';}
				out_reads<<tag<<endl;
				out_reads<<read_str<<endl;
			}


		}



	}

	//cout<<"boundary"<< boundary<<endl;
	//cout<<"removed"<<removed<<endl;
	//cout<<"bridge"<<bridge<<endl;

}


void CollectingNonContainedReads0(struct hashtable0 *ht0,struct hashtable0 *merge_ht0, int K_size,vector<string>& filenames_vt, contigs_info * contigs_info,string ContigFilename)
{
	
	time_t beg_time,read_time;
	string in_fname=ContigFilename;
	ifstream Contigs_in(in_fname.c_str());
	size_t num_Contigs=0;
	string tag,s,kmer,str,seq_s;
	uint64_t hv;
	uint64_t f_seq[100];
	size_t hash_idx;
	bool found;
	bool flip_1,flip_2,flip_0;
	size_t ht_sz;
	size_t numReads=0;
	ofstream out_reads("reads_temp.txt");

	int boundary=0,removed=0,bridge=0;

	int Kmer_arr_sz=K_size/32+1;
	int rem1=K_size%32;
	if(rem1==0)
	{Kmer_arr_sz--;}

	ht_sz=ht0->ht_sz;
	
	cout<<"Contigs remapping."<<endl;
	if(ContigFilename=="Contigs.txt")
	{
		ContigsRemapping0(ht0, K_size, contigs_info,ContigFilename,0);
		
	}

	AppendMergeHT0(ht0, merge_ht0,Kmer_arr_sz);

	//cout<<"Collecting informative reads."<<endl;
	time(&beg_time);

	int64_t dist_sum=0,p_cnt=0;
	int mean_dist=0;
	int lib_no=0;
	
	for(size_t ii=0;ii<filenames_vt.size();ii++)
	{
		lib_no++;
	
		dist_sum=0,p_cnt=0;
		cout<<"Processing library: "<<ii<<endl;
		ifstream in_reads;
		in_reads.open(filenames_vt[ii].c_str());
		
		struct read_t Read1,Read2;

		Read1.read_bits =(uint64_t*) malloc(1000000/4+100);
		Read2.read_bits =(uint64_t*) malloc(1000000/4+100);

		int nLines1=0,nLines2=0;


		bool fq_flag=0;


		getline(in_reads,str);
		if(fq_flag==0&&str[0]=='@')
		{
			fq_flag=1;	
		}
		in_reads.close();

		in_reads.clear();
		in_reads.open(filenames_vt[ii].c_str());

		bool read_success=0;

		read_success=1;

		string tag,qs,n_tag;
		string QS_s;

		while(read_success)
		{
			if(fq_flag)
			{
				read_success=get_a_fastq_read(in_reads,tag,seq_s,QS_s);
					
			}
			else
			{
				read_success=get_a_fasta_read(in_reads,tag,seq_s,n_tag);
			
			}	
			if(read_success==0)
			{break;}



				
			int seq_sz=seq_s.size();
			
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
			
			if(bad_flag)
			{continue;}
							

			bad_flag=0;
					
			numReads++;


						
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
					


			Init_Read(seq_s,Read1);
			string read_str=seq_s;
			seq_s.clear();
			

			//uint64_t bits1;
			int contig_no=-1;
			bool output_current_read=0;

			for(int i=0;i<Read1.readLen-K_size+1;++i )
			{
				uint64_t bits1[100];
				get_sub_arr(Read1.read_bits,Read1.readLen,i,K_size,bits1);

				memcpy(f_seq,bits1,sizeof(uint64_t)*Kmer_arr_sz);
				get_rev_comp_seq_arr(f_seq,K_size,Kmer_arr_sz);
				flip_0=0;
				if(uint64_t_cmp(bits1,f_seq,Kmer_arr_sz)>0)
				{
					memcpy(bits1,f_seq,sizeof(uint64_t)*Kmer_arr_sz);
					flip_0=1;
				}

				hv=MurmurHash64A(bits1,sizeof(uint64_t)*Kmer_arr_sz,0);

				hash_idx=(size_t) (hv%ht_sz);
				struct bucket0 **ptr1;
				ptr1= &(ht0->store_pos[hash_idx]);
				found=look_up_in_a_list0(bits1,&ptr1,Kmer_arr_sz);
				if(found)
				{
					
					if((*ptr1)->kmer_info.removed==0)
					{
					
						output_current_read=1;
						if(flip_0^(*ptr1)->kmer_info.flip)
						{
							if(i+50<=(*ptr1)->kmer_info.cod     &&              (Read1.readLen-i+50)<=(contigs_info->contig_sz_vt[(*ptr1)->kmer_info.contig_no]-(*ptr1)->kmer_info.cod))
							{
								output_current_read=0;
								//contained...
							}
								
						}
						else
						{
							if( (Read1.readLen-(i+K_size)+50)<=(*ptr1)->kmer_info.cod   &&   (i+K_size+50)  <=(contigs_info->contig_sz_vt[(*ptr1)->kmer_info.contig_no]-(*ptr1)->kmer_info.cod))
							{
								output_current_read=0;
								//contained...
							}
						}
	
						break;
								
					}
					else
					{
						output_current_read=1;
						removed++;
						break;


					}
	
				}
			}

			if(output_current_read)
			{
				if(tag[0]=='@')
				{tag[0]='>';}
				out_reads<<tag<<endl;
				out_reads<<read_str<<endl;
			}
		}
	}


	//cout<<"boundary"<< boundary<<endl;
	//cout<<"removed"<<removed<<endl;
	//cout<<"bridge"<<bridge<<endl;

}


int BFSearchGapCloser_ctg(struct contigs_info *contigs_info,int max_depth,int max_dist,vector<int> &ctgs_in_current_read,SuperRead_ctg *SuperRead)
{

	
	
	SuperRead->Paths.clear();
	SuperRead->PathFound=0;
	SuperRead->PathsCnt=0;
	bool correct=1;
	for(int i=0;i<ctgs_in_current_read.size()-1;++i)
	{
		int current_ctg=ctgs_in_current_read[i];
		int next_ctg=ctgs_in_current_read[i+1];
		if(current_ctg>0)
		{
			if(contigs_info->contig_adjacency_right[current_ctg].count(next_ctg)==0)
			{
				correct=0;
				break;
			}
		}
		else
		{
			if(contigs_info->contig_adjacency_left[-current_ctg].count(-next_ctg)==0)
			{
				correct=0;
				break;
			}
		}
	}

	if(correct==1)
	{
		SuperRead->PathFound=1;
		SuperRead->Paths.push_back(ctgs_in_current_read);
		return 0;
	}
	int beg_ctg=ctgs_in_current_read[0];
	int end_ctg=ctgs_in_current_read[ctgs_in_current_read.size()-1];


	map<uint64_t,struct BFS_path_info_v3 > Visited_Path;
	map< int,int > stacked_nodes;
	
	int max_stack=500;
	int DepthTh=max_depth;//min(300/gap,20);
	int LenTh=max_dist;
	bool RIGHT=0;

	map<int , list<int> > dist_ctgs;//neighborset
	dist_ctgs[0].push_back(beg_ctg);
	int NBs=1;

	int dist_searched=0;
	int new_node=beg_ctg;
	
	stacked_nodes[new_node]=1;
	//search direction

	Visited_Path[new_node].depth=1;
	Visited_Path[new_node].len=contigs_info->contig_sz_vt[abs(new_node)];
	Visited_Path[new_node].last_ctg=0;
	Visited_Path[new_node].RepeatVisits=0;
	
	map<int , list<int> >::iterator NB_it=dist_ctgs.begin();
	while(1)
	{
		NB_it=dist_ctgs.begin();
		
		if(NB_it==dist_ctgs.end())
		{break;}
		if(NB_it->first>max_dist)
		{
			if(SuperRead->PathsCnt!=1)
			{return -10000;}
			else
			{break;}
		}
		if(NB_it->second.size()==0)
		{dist_ctgs.erase(NB_it->first);continue;}
		if(NBs>max_stack)
		{
			SuperRead->PathsCnt=0;
			return -10000;
			break;
		}
		new_node=NB_it->second.front();

		NB_it->second.pop_front();
		NBs--;
		
		if(NB_it->second.size()==0)
		{
			dist_ctgs.erase(NB_it->first);
		}		

		if(Visited_Path[new_node].depth>DepthTh||Visited_Path[new_node].len>LenTh)
		{continue;}
		
		int current_ctg=new_node;
		
		if(current_ctg>0)
		{
			RIGHT=1;
		}
		else
		{
			RIGHT=0;
		}
		
		if(RIGHT)
		{
			map<int,struct adjacent_contig_info>::iterator tmp_it;
			for(tmp_it=contigs_info->contig_adjacency_right[abs(current_ctg)].begin();tmp_it!=contigs_info->contig_adjacency_right[abs(current_ctg)].end();++tmp_it)
			{
				
				int next_ctg=tmp_it->first;
								
				if(next_ctg==end_ctg)
				{
					SuperRead->PathsCnt++;
					SuperRead->PathFound=1;
					SuperRead->PathsLength=(int)(Visited_Path[new_node].len+contigs_info->contig_sz_vt[abs(next_ctg)]);

					int  backtrack_ctg=end_ctg;
					int depth=Visited_Path[backtrack_ctg].depth+5;
					vector<int> path;
					path.push_back(end_ctg);
					for (int i=depth;i>0;--i)
					{
						if(Visited_Path[backtrack_ctg].depth==1)
						{
							break;
						}
						if(Visited_Path[backtrack_ctg].RepeatVisits>0)
						{
							SuperRead->PathsCnt+=Visited_Path[backtrack_ctg].RepeatVisits;
						}
						int last_ctg=Visited_Path[backtrack_ctg].last_ctg;
						path.push_back(last_ctg);
						backtrack_ctg=last_ctg;
					}
					
				}

				// not in stack
				if(stacked_nodes[next_ctg]==0)
				{
					//Visited_Path[*ptr].cov=(int)(Visited_Path[new_node].cov+edge_ptr->edge_cov);
					Visited_Path[next_ctg].depth=(Visited_Path[new_node].depth+1);
						
					int cum_len=(int)(Visited_Path[new_node].len+contigs_info->contig_sz_vt[abs(next_ctg)]);
					Visited_Path[next_ctg].len=cum_len;//(Visited_Path[new_node].len+edge_ptr->len+1);
					Visited_Path[next_ctg].last_ctg=new_node;
					Visited_Path[next_ctg].RepeatVisits=0;
					stacked_nodes[next_ctg]=1;
					dist_ctgs[cum_len].push_back(next_ctg);
					NBs++;

				}
				else
				{
					//multiple paths, should do something
					Visited_Path[next_ctg].RepeatVisits++;
					
				}

		
			}

			


		}
		else
		{


			map<int,struct adjacent_contig_info>::iterator tmp_it;
			for(tmp_it=contigs_info->contig_adjacency_left[abs(current_ctg)].begin();tmp_it!=contigs_info->contig_adjacency_left[abs(current_ctg)].end();++tmp_it)
			{
				
				int next_ctg=-tmp_it->first;
								
				if(next_ctg==end_ctg)
				{
					SuperRead->PathsCnt++;
					SuperRead->PathFound=1;
					SuperRead->PathsLength=(int)(Visited_Path[new_node].len+contigs_info->contig_sz_vt[abs(next_ctg)]);


					int  backtrack_ctg=end_ctg;
					int depth=Visited_Path[backtrack_ctg].depth+5;
					vector<int> path;
					path.push_back(end_ctg);
					for (int i=depth;i>0;--i)
					{
						if(Visited_Path[backtrack_ctg].depth==1)
						{
							break;
						}
						if(Visited_Path[backtrack_ctg].RepeatVisits>0)
						{
							SuperRead->PathsCnt+=Visited_Path[backtrack_ctg].RepeatVisits;
						}
						int last_ctg=Visited_Path[backtrack_ctg].last_ctg;
						path.push_back(last_ctg);
						backtrack_ctg=last_ctg;
					
					}
					

						
				}

				// not in stack
				if(stacked_nodes[next_ctg]==0)
				{
					//Visited_Path[*ptr].cov=(int)(Visited_Path[new_node].cov+edge_ptr->edge_cov);
					Visited_Path[next_ctg].depth=(Visited_Path[new_node].depth+1);
						
					int cum_len=(int)(Visited_Path[new_node].len+contigs_info->contig_sz_vt[abs(next_ctg)]);
					Visited_Path[next_ctg].len=cum_len;//(Visited_Path[new_node].len+edge_ptr->len+1);
					Visited_Path[next_ctg].last_ctg=new_node;
					Visited_Path[next_ctg].RepeatVisits=0;
					stacked_nodes[next_ctg]=1;
					dist_ctgs[cum_len].push_back(next_ctg);
					NBs++;

				}
				else
				{
					//multiple paths, should do something
					Visited_Path[next_ctg].RepeatVisits++;
					
				}

			}

		}

	}


	if(SuperRead->PathsCnt!=1)
	{return -10000;}
	else
	{
		int  backtrack_ctg=end_ctg;
		int depth=Visited_Path[backtrack_ctg].depth+5;
		vector<int> path;
		path.push_back(end_ctg);
		for (int i=depth;i>0;--i)
		{

			if(Visited_Path[backtrack_ctg].depth==1)
			{
				reverse(path.begin(),path.end());
				SuperRead->Paths.push_back(path);
				return SuperRead->PathsLength;		
			}

			int last_ctg=Visited_Path[backtrack_ctg].last_ctg;
			
			path.push_back(last_ctg);
			backtrack_ctg=last_ctg;

		}
		return SuperRead->PathsLength;
	}

}

void ConstuctRefinedContigGraph(struct hashtable *ht1, struct hashtable *merge_ht1,int K_size, vector<string>& filenames_vt, contigs_info * contigs_info, string ContigFilename)
{
	time_t beg_time, read_time;
	string in_fname = ContigFilename;
	ifstream Contigs_in(in_fname.c_str());
	size_t num_Contigs = 0;
	string tag, s, kmer, str, seq_s;
	uint64_t f_seq, hv;
	struct kmer_t2 f_seq_t2;
	size_t hash_idx;
	bool found;
	bool flip_1, flip_2, flip_0;
	size_t ht_sz;
	size_t numReads = 0;
	ofstream out_refined_graph("RefinedContigGraph.txt");
	map<vector<int>, map<string, int> > temp_map;
	int boundary = 0, removed = 0, bridge = 0;
	if (K_size <= 32)
	{
		ht_sz = ht1->ht_sz;
	}
	hashtable2 ht2;
	cout << "Contigs remapping." << endl;
	if (ContigFilename == "Contigs.txt")
	{
		ContigsRemapping(ht1, &ht2, K_size, contigs_info, ContigFilename, 0);

	}
	if (K_size <= 32)
	{
		//AppendMergeHT(ht1, merge_ht1);
	}
	else
	{
		// AppendMergeHT2(ht2,merge_ht2);

	}
	//BuildContigAdjacency(ht1, ht2, contigs_info, K_size,ContigFilename);	

	//cout<<"Collecting informative reads."<<endl;
	time(&beg_time);

	int64_t dist_sum = 0, p_cnt = 0;
	int mean_dist = 0;
	int lib_no = 0;

	for (size_t ii = 0; ii<filenames_vt.size(); ii++)
	{
		lib_no++;

		dist_sum = 0, p_cnt = 0;
		cout << "Processing library: " << ii << endl;
		ifstream in_reads;
		in_reads.open(filenames_vt[ii].c_str());

		struct read_t Read1, Read2;

		Read1.read_bits = (uint64_t*)malloc(1000000 / 4 + 100);
		Read2.read_bits = (uint64_t*)malloc(1000000 / 4 + 100);

		int nLines1 = 0, nLines2 = 0;


		bool fq_flag = 0;


		getline(in_reads, str);
		if (fq_flag == 0 && str[0] == '@')
		{
			fq_flag = 1;
		}
		in_reads.close();

		in_reads.clear();
		in_reads.open(filenames_vt[ii].c_str());

		bool read_success = 0;

		read_success = 1;

		string tag, qs, n_tag;
		string QS_s;

		while (read_success)
		{
			if (fq_flag)
			{
				read_success = get_a_fastq_read(in_reads, tag, seq_s, QS_s);

			}
			else
			{
				read_success = get_a_fasta_read(in_reads, tag, seq_s, n_tag);

			}
			if (read_success == 0)
			{
				break;
			}




			int seq_sz = seq_s.size();
			int readLen = seq_sz;
			if (seq_s.size() == 0)
			{
				cout << "Empty sequence!" << endl;
				continue;
			}
			bool bad_flag = 0;
			int numN = 0;
			for (int i = 0; i<seq_sz; ++i)
			{
				if (seq_s[i] != 'A'&&seq_s[i] != 'C'&&seq_s[i] != 'G'&&seq_s[i] != 'T'&&seq_s[i] != 'N')
				{
					bad_flag = 1;
					break;
				}
				if (seq_s[i] == 'N')
				{
					numN++;
				}
			}

			if (bad_flag)
			{
				continue;
			}


			bad_flag = 0;

			numReads++;



			int nN = seq_sz - 1, isN = -1;
			for (int i = 0; i<seq_sz; ++i)
			{

				if (seq_s[i] == '-' || seq_s[i] == 'N')
				{
					if (i <= seq_sz / 2)
					{
						isN = i;
						continue;
					}
					else
					{
						nN = i - 1;
						break;
					}
				}
			}
			int s = 0;
			if ((nN - isN) <= seq_sz / 2)
			{
				bad_flag = 1;
			}

			if (bad_flag == 1)
			{
				seq_s.clear();
				continue;
			}

			if (isN >= 0)
			{
				for (int i = isN + 1; i <= nN; ++i)
				{
					seq_s[s] = seq_s[i];
					s++;
				}
				seq_s[s] = '\0';
				seq_s.resize(s);
			}



			Init_Read(seq_s, Read1);
			string read_str = seq_s;
			seq_s.clear();
			LongReadContigIndex LongReadContigIndex;

			uint64_t bits1;
			int contig_no = -1;
			bool output_current_read = 0;
			vector<int> ctgs_vt, cod_vt, pos_vt, offset_vt;
			bool bubble_flag = 0;
			for (int i = 0; i<Read1.readLen - K_size + 1; ++i)
			{
				if (K_size <= 32)
				{
					get_sub_arr(Read1.read_bits, Read1.readLen, i, K_size, &bits1);
					f_seq = get_rev_comp_seq(bits1, K_size);
					flip_0 = 0;
					if (bits1>f_seq)
					{
						bits1 = f_seq;
						flip_0 = 1;
					}

					hv = MurmurHash64A(&bits1, sizeof(bits1), 0);
					hash_idx = (size_t)(hv%ht_sz);
					struct bucket **ptr1;
					ptr1 = &(ht1->store_pos[hash_idx]);
					found = look_up_in_a_list(bits1, &ptr1);
					if (found)
					{
						int cod = (*ptr1)->kmer_info.cod;
						int contig_no = (*ptr1)->kmer_info.contig_no;

						if (cod > readLen + 50 && cod < (contigs_info->contig_sz_vt[abs(contig_no)] - readLen - 50))
						{

							continue;//contained hit
						}

						LongReadContigIndex.LR2CTG[i].contig_no = contig_no;
						LongReadContigIndex.LR2CTG[i].pos = (*(ptr1))->kmer_info.cod;
						LongReadContigIndex.LR2CTG[i].flip = flip_0 ^ (*(ptr1))->kmer_info.flip;

						int offset = 0;
						if ((*ptr1)->kmer_info.removed == 0 && (*ptr1)->kmer_info.contig_no>0 && (*ptr1)->kmer_info.repeat == 0)
						{
							if ((flip_0 ^ (*ptr1)->kmer_info.flip) == 0)
							{
								ctgs_vt.push_back((*ptr1)->kmer_info.contig_no);
								offset = i - (*ptr1)->kmer_info.cod + contigs_info->contig_sz_vt[(*(ptr1))->kmer_info.contig_no] / 2;
								//cout << offset << endl;
								LongReadContigIndex.CTG2LR[(*(ptr1))->kmer_info.contig_no].push_back(i);
								LongReadContigIndex.CTG2LR_2[(*(ptr1))->kmer_info.contig_no].push_back(offset);

							}
							else
							{
								ctgs_vt.push_back(-(*ptr1)->kmer_info.contig_no);;
								offset = i + (*ptr1)->kmer_info.cod - contigs_info->contig_sz_vt[(*(ptr1))->kmer_info.contig_no] / 2 + (K_size + 1) / 2;
								//cout << offset << endl;
								LongReadContigIndex.CTG2LR[-(*(ptr1))->kmer_info.contig_no].push_back(i);
								LongReadContigIndex.CTG2LR_2[-(*(ptr1))->kmer_info.contig_no].push_back(offset);

							}
							cod_vt.push_back((*ptr1)->kmer_info.cod);
							pos_vt.push_back(i);
							offset_vt.push_back(offset);
						}

					}
				}
				

			}


			map<int, vector<int> >::iterator tmp_it, tmp_it2;
			tmp_it2 = LongReadContigIndex.CTG2LR_2.begin();
			for (tmp_it = LongReadContigIndex.CTG2LR.begin(); tmp_it != LongReadContigIndex.CTG2LR.end(); ++tmp_it)
			{
				sort(tmp_it->second.begin(), tmp_it->second.end());
				sort(tmp_it2->second.begin(), tmp_it2->second.end());
				//LongReadContigIndex.layout[tmp_it->second[tmp_it->second.size() / 2]].push_back(tmp_it->first);//ctg_no
				//LongReadContigIndex.layout[tmp_it->second[tmp_it->second.size() / 2]].push_back(tmp_it->second.size());//cov
				LongReadContigIndex.layout[tmp_it->second[tmp_it->second.size() / 2]].push_back(tmp_it->first);//ctg_no
				LongReadContigIndex.layout[tmp_it->second[tmp_it->second.size() / 2]].push_back(tmp_it->second.size());//cov
				LongReadContigIndex.layout[tmp_it->second[tmp_it->second.size() / 2]].push_back(tmp_it2->second[tmp_it2->second.size() / 2]);//coord2
				tmp_it2++;
			}

			int contig_matches = 0;

			for (tmp_it = LongReadContigIndex.layout.begin(); tmp_it != LongReadContigIndex.layout.end(); ++tmp_it)
			{

				if (tmp_it->second[1] > 0)
				{
					contig_matches++;
					//LongReadContigIndex_info << tmp_it->first << ", " << tmp_it->second[0] << ", " << tmp_it->second[1] << endl;
				}
			}


			vector<int> final_ctgs_vt, final_cod_vt, overlap_vt, final_offset_vt;
			if (ctgs_vt.size()>1)
			{
				final_ctgs_vt.push_back(ctgs_vt[0]);
				final_cod_vt.push_back(cod_vt[0]);
				final_offset_vt.push_back(offset_vt[0]);
				for (int i = 1; i<ctgs_vt.size(); ++i)
				{

					//split the cases for easier debugging
					if (ctgs_vt[i] != ctgs_vt[i - 1])
					{
						final_ctgs_vt.push_back(ctgs_vt[i]);
						final_cod_vt.push_back(cod_vt[i]);
						final_offset_vt.push_back(offset_vt[i]);
						overlap_vt.push_back(pos_vt[i] - pos_vt[i - 1]);//
					}
					else
					{
						//ctgs_vt[i]==ctgs_vt[i-1]
						if ((ctgs_vt[i]>0 && cod_vt[i]<cod_vt[i - 1]) || (ctgs_vt[i]<0 && cod_vt[i]>cod_vt[i - 1]))
						{
							final_ctgs_vt.push_back(ctgs_vt[i]);
							final_cod_vt.push_back(cod_vt[i]);
							final_offset_vt.push_back(offset_vt[i]);
							overlap_vt.push_back(pos_vt[i] - pos_vt[i - 1]);//
						}
					}

				}
			}

			if (final_ctgs_vt.size()>1)
			{
				output_current_read = 1;
			}



			map<int, KmerInContig>::iterator it1, it2;

			for (it1 = LongReadContigIndex.LR2CTG.begin(); it1 != LongReadContigIndex.LR2CTG.end(); ++it1)
			{

				it2 = it1;
				it2++;
				if (it2 == LongReadContigIndex.LR2CTG.end())
				{
					break;
				}
				if (it2->second.contig_no != it1->second.contig_no)
				{
					bool flip1, flip2;
					flip1 = it1->second.flip;
					flip2 = it2->second.flip;
					int pos1 = it1->second.pos;
					int pos2 = it2->second.pos;
					int extra_bases1, extra_bases2;
					int contig1 = it1->second.contig_no;
					int contig2 = it2->second.contig_no;

					if (flip1)
					{
						contig1 = -contig1;
					}
					if (flip2)
					{
						contig2 = -contig2;
					}
					int dist = it2->first - it1->first;
					if (flip1 == 0)
					{
						extra_bases1 = contigs_info->contig_sz_vt[abs(contig1)] - pos1 - 1;
					}
					else
					{
						extra_bases1 = pos1;
					}

					if (flip2 == 0)
					{
						extra_bases2 = pos2;
					}
					else
					{
						extra_bases2 = contigs_info->contig_sz_vt[abs(contig2)] - pos2 - 1;
					}
					int extra_bases = extra_bases1 + extra_bases2;
					dist = dist - extra_bases - 1;

					string bridge;
					if (dist > 0)
					{
						bridge = read_str.substr(it1->first + extra_bases1, dist);
					}

					vector<int> temp_vec1, temp_vec2;
					temp_vec1.push_back(contig1);
					temp_vec1.push_back(contig2);
					temp_vec1.push_back(dist);

					temp_vec2.push_back(-contig2);
					temp_vec2.push_back(-contig1);
					temp_vec2.push_back(dist);

					bool take_rc = 0;
					for (int c = 0; c < 2; ++c)
					{
						if (temp_vec2[c] < temp_vec1[c])
						{
							take_rc = 1;
							break;
						}
						if (temp_vec2[c] > temp_vec1[c])
						{
							break;
						}

					}
					if (take_rc)
					{
						temp_vec1 = temp_vec2;
						if (bridge.size() > 0)
						{
							reverse(bridge.begin(), bridge.end());
							complement_str(bridge);
						}
					}

					temp_map[temp_vec1][bridge]++;


				}

			}

		}

	}


	map<vector<int>, map<string, int> >::iterator it;
	for (it = temp_map.begin(); it != temp_map.end(); ++it)
	{
		int contig1 = it->first[0];
		int contig2 = it->first[1];

		if (contig1 != 0 && contig2 != 0)
		{
			string best_str;
			int max_cov = 0;
			for (map<string, int>::iterator tmp_it = it->second.begin(); tmp_it != it->second.end(); ++tmp_it)
			{
				if (tmp_it->second > max_cov)
				{
					max_cov = tmp_it->second;
					best_str = tmp_it->first;
				}
			}

			if (contig1 > 0)
			{
				out_refined_graph << contig1 << " " << contig2;
			}
			else
			{
				out_refined_graph << contig1 << " " << -contig2;
			}

			out_refined_graph << " " << it->first[2] << " " << max_cov << " " << best_str << endl;
		}
	}

}

void ConstuctRefinedContigGraph0(struct hashtable0 *ht, struct hashtable0 *merge_ht, int K_size, vector<string>& filenames_vt, contigs_info * contigs_info, string ContigFilename)
{

	time_t beg_time, read_time;
	string in_fname = ContigFilename;
	ifstream Contigs_in(in_fname.c_str());
	size_t num_Contigs = 0;
	string tag, s, kmer, str, seq_s;
	uint64_t f_seq, hv;
	struct kmer_t2 f_seq_t2;
	size_t hash_idx;
	bool found;
	bool flip_1, flip_2, flip_0;
	size_t ht_sz;
	size_t numReads = 0;
	ofstream out_refined_graph("RefinedContigGraph.txt");
	map<vector<int>, map<string, int> > temp_map;

	int Kmer_arr_sz = K_size / 32 + 1;
	int rem1 = K_size % 32;
	if (rem1 == 0)
	{
		Kmer_arr_sz--;
	}

	int boundary = 0, removed = 0, bridge = 0;

	ht_sz = ht->ht_sz;


	cout << "Contigs remapping." << endl;
	if (ContigFilename == "Contigs.txt")
	{
		ContigsRemapping0(ht, K_size, contigs_info, ContigFilename, 0);

	}

	// AppendMergeHT0(ht, merge_ht,Kmer_arr_sz);
	//cout<<"Collecting informative reads."<<endl;
	//BuildContigAdjacency0(ht,contigs_info, K_size,ContigFilename);

	time(&beg_time);

	int64_t dist_sum = 0, p_cnt = 0;
	int mean_dist = 0;
	int lib_no = 0;

	for (size_t ii = 0; ii<filenames_vt.size(); ii++)
	{
		lib_no++;

		dist_sum = 0, p_cnt = 0;
		cout << "Processing library: " << ii << endl;
		ifstream in_reads;
		in_reads.open(filenames_vt[ii].c_str());

		struct read_t Read1, Read2;

		Read1.read_bits = (uint64_t*)malloc(1000000 / 4 + 100);
		Read2.read_bits = (uint64_t*)malloc(1000000 / 4 + 100);

		int nLines1 = 0, nLines2 = 0;


		bool fq_flag = 0;


		getline(in_reads, str);
		if (fq_flag == 0 && str[0] == '@')
		{
			fq_flag = 1;
		}
		in_reads.close();

		in_reads.clear();
		in_reads.open(filenames_vt[ii].c_str());

		bool read_success = 0;

		read_success = 1;

		string tag, qs, n_tag;
		string QS_s;

		while (read_success)
		{
			if (fq_flag)
			{
				read_success = get_a_fastq_read(in_reads, tag, seq_s, QS_s);

			}
			else
			{
				read_success = get_a_fasta_read(in_reads, tag, seq_s, n_tag);

			}
			if (read_success == 0)
			{
				break;
			}




			int seq_sz = seq_s.size();
			int readLen = seq_sz;
			if (seq_s.size() == 0)
			{
				cout << "Empty sequence!" << endl;
				continue;
			}
			bool bad_flag = 0;
			int numN = 0;
			for (int i = 0; i<seq_sz; ++i)
			{
				if (seq_s[i] != 'A'&&seq_s[i] != 'C'&&seq_s[i] != 'G'&&seq_s[i] != 'T'&&seq_s[i] != 'N')
				{
					bad_flag = 1;
					break;
				}
				if (seq_s[i] == 'N')
				{
					numN++;
				}
			}

			if (bad_flag)
			{
				continue;
			}


			bad_flag = 0;

			numReads++;



			int nN = seq_sz - 1, isN = -1;
			for (int i = 0; i<seq_sz; ++i)
			{

				if (seq_s[i] == '-' || seq_s[i] == 'N')
				{
					if (i <= seq_sz / 2)
					{
						isN = i;
						continue;
					}
					else
					{
						nN = i - 1;
						break;
					}
				}
			}
			int s = 0;
			if ((nN - isN) <= seq_sz / 2)
			{
				bad_flag = 1;
			}

			if (bad_flag == 1)
			{
				seq_s.clear();
				continue;
			}

			if (isN >= 0)
			{
				for (int i = isN + 1; i <= nN; ++i)
				{
					seq_s[s] = seq_s[i];
					s++;
				}
				seq_s[s] = '\0';
				seq_s.resize(s);
			}



			Init_Read(seq_s, Read1);
			string read_str = seq_s;
			seq_s.clear();
			LongReadContigIndex LongReadContigIndex;

			int contig_no = -1;
			bool output_current_read = 0;
			vector<int> ctgs_vt, cod_vt, pos_vt, offset_vt;
			bool bubble_flag = 0;
			for (int i = 0; i<Read1.readLen - K_size + 1; ++i)
			{


				uint64_t bits1[100], f_seq[100];
				get_sub_arr(Read1.read_bits, Read1.readLen, i, K_size, bits1);
				memcpy(f_seq, bits1, sizeof(uint64_t)*Kmer_arr_sz);
				get_rev_comp_seq_arr(f_seq, K_size, Kmer_arr_sz);
				flip_0 = 0;
				if (uint64_t_cmp(bits1, f_seq, Kmer_arr_sz)>0)
				{
					memcpy(bits1, f_seq, sizeof(uint64_t)*Kmer_arr_sz);

					flip_0 = 1;
				}

				hv = MurmurHash64A(bits1, sizeof(uint64_t)*Kmer_arr_sz, 0);

				hash_idx = (size_t)(hv%ht_sz);
				struct bucket0 **ptr1;
				ptr1 = &(ht->store_pos[hash_idx]);
				found = look_up_in_a_list0(bits1, &ptr1, Kmer_arr_sz);
				if (found)
				{
					int cod = (*ptr1)->kmer_info.cod;
					int contig_no = (*ptr1)->kmer_info.contig_no;

					if (cod > readLen + 50 && cod < (contigs_info->contig_sz_vt[abs(contig_no)] - readLen - 50))
					{

						continue;//contained hit
					}

					LongReadContigIndex.LR2CTG[i].contig_no = contig_no;
					LongReadContigIndex.LR2CTG[i].pos = (*(ptr1))->kmer_info.cod;
					LongReadContigIndex.LR2CTG[i].flip = flip_0 ^ (*(ptr1))->kmer_info.flip;


					int offset = 0;
					if ((*ptr1)->kmer_info.removed == 0 && (*ptr1)->kmer_info.contig_no > 0 && (*ptr1)->kmer_info.repeat == 0)
					{
						if ((flip_0 ^ (*ptr1)->kmer_info.flip) == 0)
						{
							ctgs_vt.push_back((*ptr1)->kmer_info.contig_no);
							offset = i - (*ptr1)->kmer_info.cod + contigs_info->contig_sz_vt[(*(ptr1))->kmer_info.contig_no] / 2;
							LongReadContigIndex.CTG2LR_2[(*(ptr1))->kmer_info.contig_no].push_back(offset);
							LongReadContigIndex.CTG2LR[(*(ptr1))->kmer_info.contig_no].push_back(i);

						}
						else
						{
							ctgs_vt.push_back(-(*ptr1)->kmer_info.contig_no);;
							offset = i + (*ptr1)->kmer_info.cod - contigs_info->contig_sz_vt[(*(ptr1))->kmer_info.contig_no] / 2 + (K_size + 1) / 2;
							LongReadContigIndex.CTG2LR_2[-(*(ptr1))->kmer_info.contig_no].push_back(offset);
							LongReadContigIndex.CTG2LR[-(*(ptr1))->kmer_info.contig_no].push_back(i);

						}
						cod_vt.push_back((*ptr1)->kmer_info.cod);
						pos_vt.push_back(i);
						offset_vt.push_back(offset);
					}
				}



			}


			map<int, vector<int> >::iterator tmp_it, tmp_it2;
			tmp_it2 = LongReadContigIndex.CTG2LR_2.begin();
			for (tmp_it = LongReadContigIndex.CTG2LR.begin(); tmp_it != LongReadContigIndex.CTG2LR.end(); ++tmp_it)
			{
				sort(tmp_it->second.begin(), tmp_it->second.end());
				sort(tmp_it2->second.begin(), tmp_it2->second.end());
				//LongReadContigIndex.layout[tmp_it->second[tmp_it->second.size() / 2]].push_back(tmp_it->first);//ctg_no
				//LongReadContigIndex.layout[tmp_it->second[tmp_it->second.size() / 2]].push_back(tmp_it->second.size());//cov
				LongReadContigIndex.layout[tmp_it->second[tmp_it->second.size() / 2]].push_back(tmp_it->first);//ctg_no
				LongReadContigIndex.layout[tmp_it->second[tmp_it->second.size() / 2]].push_back(tmp_it->second.size());//cov
				LongReadContigIndex.layout[tmp_it->second[tmp_it->second.size() / 2]].push_back(tmp_it2->second[tmp_it2->second.size() / 2]);//coord2
				tmp_it2++;
			}

			int contig_matches = 0;

			for (tmp_it = LongReadContigIndex.layout.begin(); tmp_it != LongReadContigIndex.layout.end(); ++tmp_it)
			{

				if (tmp_it->second[1] > 0)
				{
					contig_matches++;
					//LongReadContigIndex_info << tmp_it->first << ", " << tmp_it->second[0] << ", " << tmp_it->second[1] << endl;
				}
			}


		
			vector<int> final_ctgs_vt, final_cod_vt, overlap_vt, final_offset_vt;
			if (ctgs_vt.size()>1)
			{
				final_ctgs_vt.push_back(ctgs_vt[0]);
				final_cod_vt.push_back(cod_vt[0]);
				final_offset_vt.push_back(offset_vt[0]);
				for (int i = 1; i<ctgs_vt.size(); ++i)
				{

					//split the cases for easier debugging
					if (ctgs_vt[i] != ctgs_vt[i - 1])
					{
						final_ctgs_vt.push_back(ctgs_vt[i]);
						final_cod_vt.push_back(cod_vt[i]);
						final_offset_vt.push_back(offset_vt[i]);
						overlap_vt.push_back(pos_vt[i] - pos_vt[i - 1]);//
					}
					else
					{
						//ctgs_vt[i]==ctgs_vt[i-1]
						if ((ctgs_vt[i]>0 && cod_vt[i]<cod_vt[i - 1]) || (ctgs_vt[i]<0 && cod_vt[i]>cod_vt[i - 1]))
						{
							final_ctgs_vt.push_back(ctgs_vt[i]);
							final_cod_vt.push_back(cod_vt[i]);
							final_offset_vt.push_back(offset_vt[i]);
							overlap_vt.push_back(pos_vt[i] - pos_vt[i - 1]);//
						}
					}

				}
			}

			if (final_ctgs_vt.size()>1)
			{
				output_current_read = 1;
			}

		


			map<int, KmerInContig>::iterator it1, it2;

			for (it1 = LongReadContigIndex.LR2CTG.begin(); it1 != LongReadContigIndex.LR2CTG.end(); ++it1)
			{

				it2 = it1;
				it2++;
				if (it2 == LongReadContigIndex.LR2CTG.end())
				{
					break;
				}
				if (it2->second.contig_no != it1->second.contig_no)
				{
					bool flip1, flip2;
					flip1 = it1->second.flip;
					flip2 = it2->second.flip;
					int pos1 = it1->second.pos;
					int pos2 = it2->second.pos;
					int extra_bases1, extra_bases2;
					int contig1 = it1->second.contig_no;
					int contig2 = it2->second.contig_no;

					if (flip1)
					{
						contig1 = -contig1;
					}
					if (flip2)
					{
						contig2 = -contig2;
					}
					int dist = it2->first - it1->first;
					if (flip1 == 0)
					{
						extra_bases1 = contigs_info->contig_sz_vt[abs(contig1)] - pos1 - 1;
					}
					else
					{
						extra_bases1 = pos1;
					}

					if (flip2 == 0)
					{
						extra_bases2 = pos2;
					}
					else
					{
						extra_bases2 = contigs_info->contig_sz_vt[abs(contig2)] - pos2 - 1;
					}
					int extra_bases = extra_bases1 + extra_bases2;
					dist = dist - extra_bases - 1;

					string bridge;
					if (dist > 0)
					{
						bridge = read_str.substr(it1->first + extra_bases1, dist);
					}


					vector<int> temp_vec1, temp_vec2;
					temp_vec1.push_back(contig1);
					temp_vec1.push_back(contig2);
					temp_vec1.push_back(dist);

					temp_vec2.push_back(-contig2);
					temp_vec2.push_back(-contig1);
					temp_vec2.push_back(dist);

					bool take_rc = 0;
					for (int c = 0; c < 2; ++c)
					{
						if (temp_vec2[c] < temp_vec1[c])
						{
							take_rc = 1;
							break;
						}
						if (temp_vec2[c] > temp_vec1[c])
						{
							break;
						}

					}
					if (take_rc)
					{
						temp_vec1 = temp_vec2;
						if (bridge.size() > 0)
						{
							reverse(bridge.begin(), bridge.end());
							complement_str(bridge);
						}
					}

					temp_map[temp_vec1][bridge]++;

				}

			}
		}

	}


	map<vector<int>, map<string, int> >::iterator it;
	for (it = temp_map.begin(); it != temp_map.end(); ++it)
	{
		int contig1 = it->first[0];
		int contig2 = it->first[1];

		if (contig1 != 0 && contig2 != 0)
		{
			string best_str;
			int max_cov = 0;
			for (map<string, int>::iterator tmp_it = it->second.begin(); tmp_it != it->second.end(); ++tmp_it)
			{
				if (tmp_it->second > max_cov)
				{
					max_cov = tmp_it->second;
					best_str = tmp_it->first;
				}
			}

			if (contig1 > 0)
			{
				out_refined_graph << contig1 << " " << contig2;
			}
			else
			{
				out_refined_graph << contig1 << " " << -contig2;
			}

			out_refined_graph << " " << it->first[2] << " " << max_cov << " " << best_str << endl;
		}
	}
}


void CollectingNonContainedReadsSlow(struct hashtable *ht1,struct hashtable *merge_ht1, struct hashtable2 *ht2, struct hashtable2 *merge_ht2,int K_size,vector<string>& filenames_vt, contigs_info * contigs_info,string ContigFilename)
{

	time_t beg_time,read_time;
	string in_fname=ContigFilename;
	ifstream Contigs_in(in_fname.c_str());
	size_t num_Contigs=0;
	string tag,s,kmer,str,seq_s;
	uint64_t f_seq,hv;
	struct kmer_t2 f_seq_t2;
	size_t hash_idx;
	bool found;
	bool flip_1,flip_2,flip_0;
	size_t ht_sz;
	size_t numReads=0;
	
	ofstream out_refined_graph("RefinedContigGraph.txt");
	map<vector<int>, map<string,int> > temp_map;

	int boundary=0,removed=0,bridge=0;
	if(K_size<=32)
	{
		ht_sz=ht1->ht_sz;
	}
	else
	{
		ht_sz=ht2->ht_sz;
	}
	cout<<"Contigs remapping."<<endl;
	if(ContigFilename=="Contigs.txt")
	{
		ContigsRemapping(ht1,ht2, K_size, contigs_info,ContigFilename,0);
		
	}
	if(K_size<=32)
	{
		 //AppendMergeHT(ht1, merge_ht1);
	}
	else
	{
		// AppendMergeHT2(ht2,merge_ht2);
	
	}
	//BuildContigAdjacency(ht1, ht2, contigs_info, K_size,ContigFilename);	
	
	//cout<<"Collecting informative reads."<<endl;
	time(&beg_time);

	int64_t dist_sum=0,p_cnt=0;
	int mean_dist=0;
	int lib_no=0;
	
	for(size_t ii=0;ii<filenames_vt.size();ii++)
	{
		lib_no++;
	
		dist_sum=0,p_cnt=0;
		cout<<"Processing library: "<<ii<<endl;
		ifstream in_reads;
		in_reads.open(filenames_vt[ii].c_str());
		string file = filenames_vt[ii];
		int file_beg = 0;
		for (int c = 0; c < file.size(); ++c)
		{
			if (file[c] == '\\'||file[c]=='/')
			{
				file_beg = c;
			}
		}
		file = file.substr(file_beg+1,file.size());

		string Name1 = "NonContainedReadsFrom_" + file;
		string Name2 = "ReadsInfoFrom_" + file;
		ofstream out_reads(Name1.c_str());
		ofstream out_reads_info(Name2.c_str());;
		cout << "Creating: " << Name1 << endl;
		struct read_t Read1,Read2;

		Read1.read_bits =(uint64_t*) malloc(1000000/4+100);
		Read2.read_bits =(uint64_t*) malloc(1000000/4+100);

		int nLines1=0,nLines2=0;


		bool fq_flag=0;


		getline(in_reads,str);
		if(fq_flag==0&&str[0]=='@')
		{
			fq_flag=1;	
		}
		in_reads.close();

		in_reads.clear();
		in_reads.open(filenames_vt[ii].c_str());

		bool read_success=0;

		read_success=1;

		string tag,qs,n_tag;
		string QS_s;

		while(read_success)
		{
			if(fq_flag)
			{
				read_success=get_a_fastq_read(in_reads,tag,seq_s,QS_s);
					
			}
			else
			{
				read_success=get_a_fasta_read(in_reads,tag,seq_s,n_tag);
			
			}	
			if(read_success==0)
			{break;}
			
				
			int seq_sz=seq_s.size();
			int readLen = seq_sz;
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
			
			if(bad_flag)
			{continue;}
							

			bad_flag=0;
					
			numReads++;

			
						
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
					


			Init_Read(seq_s,Read1);
			string read_str=seq_s;
			seq_s.clear();
			LongReadContigIndex LongReadContigIndex;

			uint64_t bits1;
			int contig_no=-1;
			bool output_current_read=0;
			vector<int> ctgs_vt,cod_vt,pos_vt,offset_vt;
			bool bubble_flag=0;
			for(int i=0;i<Read1.readLen-K_size+1;++i )
			{
				if(K_size<=32)
				{
					get_sub_arr(Read1.read_bits,Read1.readLen,i,K_size,&bits1);
					f_seq=get_rev_comp_seq(bits1,K_size);
					flip_0=0;
					if(bits1>f_seq)
					{
						bits1=f_seq;
						flip_0=1;
					}

					hv=MurmurHash64A(&bits1,sizeof(bits1),0);
					hash_idx=(size_t) (hv%ht_sz);
					struct bucket **ptr1;
					ptr1= &(ht1->store_pos[hash_idx]);
					found=look_up_in_a_list(bits1,&ptr1);
					if(found)
					{
						int cod = (*ptr1)->kmer_info.cod;
						int contig_no = (*ptr1)->kmer_info.contig_no;
						
						if (contig_no==0||(cod > readLen + 50 && cod < (contigs_info->contig_sz_vt[abs(contig_no)] - readLen - 50)))
						{
							
							continue;//contained hit
						}
	
						LongReadContigIndex.LR2CTG[i].contig_no = contig_no;
						LongReadContigIndex.LR2CTG[i].pos = (*(ptr1))->kmer_info.cod;
						LongReadContigIndex.LR2CTG[i].flip = flip_0 ^ (*(ptr1))->kmer_info.flip;

						int offset = 0;
						if ((*ptr1)->kmer_info.removed == 0 && (*ptr1)->kmer_info.contig_no>0 && (*ptr1)->kmer_info.repeat==0)
						{
							if((flip_0^(*ptr1)->kmer_info.flip)==0)
							{
								ctgs_vt.push_back((*ptr1)->kmer_info.contig_no);
								offset = i - (*ptr1)->kmer_info.cod + contigs_info->contig_sz_vt[(*(ptr1))->kmer_info.contig_no] / 2;
								//cout << offset << endl;
								LongReadContigIndex.CTG2LR[(*(ptr1))->kmer_info.contig_no].push_back(i);
								LongReadContigIndex.CTG2LR_2[(*(ptr1))->kmer_info.contig_no].push_back(offset);

							}
							else
							{
								ctgs_vt.push_back(-(*ptr1)->kmer_info.contig_no);;	
								offset = i + (*ptr1)->kmer_info.cod - contigs_info->contig_sz_vt[(*(ptr1))->kmer_info.contig_no] / 2 + (K_size + 1) / 2;
								//cout << offset << endl;
								LongReadContigIndex.CTG2LR[-(*(ptr1))->kmer_info.contig_no].push_back(i);
								LongReadContigIndex.CTG2LR_2[-(*(ptr1))->kmer_info.contig_no].push_back(offset);

							}
							cod_vt.push_back((*ptr1)->kmer_info.cod);	
							pos_vt.push_back(i);
							offset_vt.push_back(offset);
						}
						
					}
				}
				else
				{
					if(K_size>32&&K_size<=64)
					{
						kmer_t2 bits1_t2;
						get_sub_arr(Read1.read_bits,Read1.readLen,i,K_size,bits1_t2.kmer);

						f_seq_t2=bits1_t2;
						get_rev_comp_seq_arr(f_seq_t2.kmer,K_size,2);
						flip_0=0;
						if(uint64_t_cmp(bits1_t2.kmer,f_seq_t2.kmer,2)>0)
						{
							bits1_t2=f_seq_t2;
							flip_0=1;
						}

						hv=MurmurHash64A(bits1_t2.kmer,sizeof(bits1_t2),0);

						hash_idx=(size_t) (hv%ht_sz);
						struct bucket2 **ptr1;
						ptr1= &(ht2->store_pos[hash_idx]);
						found=look_up_in_a_list2(&bits1_t2,&ptr1);
						if(found)
						{
							int cod = (*ptr1)->kmer_info.cod;
							int contig_no = (*ptr1)->kmer_info.contig_no;

							if (contig_no == 0 || (cod > readLen + 50 && cod < (contigs_info->contig_sz_vt[abs(contig_no)] - readLen - 50)))
							{
								continue;//contained hit
							}

							LongReadContigIndex.LR2CTG[i].contig_no = contig_no;
							LongReadContigIndex.LR2CTG[i].pos = (*(ptr1))->kmer_info.cod;
							LongReadContigIndex.LR2CTG[i].flip = flip_0 ^ (*(ptr1))->kmer_info.flip;

							int offset = 0;
							if ((*ptr1)->kmer_info.removed == 0 && (*ptr1)->kmer_info.contig_no>0 && (*ptr1)->kmer_info.repeat == 0)
							{
								if ((flip_0 ^ (*ptr1)->kmer_info.flip) == 0)
								{
									ctgs_vt.push_back((*ptr1)->kmer_info.contig_no);
									offset = i - (*ptr1)->kmer_info.cod + contigs_info->contig_sz_vt[(*(ptr1))->kmer_info.contig_no] / 2;
									
									LongReadContigIndex.CTG2LR_2[(*(ptr1))->kmer_info.contig_no].push_back(offset);
									LongReadContigIndex.CTG2LR[(*(ptr1))->kmer_info.contig_no].push_back(i);

								}
								else
								{
									ctgs_vt.push_back(-(*ptr1)->kmer_info.contig_no);;
									offset = i + (*ptr1)->kmer_info.cod - contigs_info->contig_sz_vt[(*(ptr1))->kmer_info.contig_no] / 2 + (K_size + 1) / 2;
									LongReadContigIndex.CTG2LR_2[-(*(ptr1))->kmer_info.contig_no].push_back(offset);
									LongReadContigIndex.CTG2LR[-(*(ptr1))->kmer_info.contig_no].push_back(i);

								}
								cod_vt.push_back((*ptr1)->kmer_info.cod);
								pos_vt.push_back(i);
								offset_vt.push_back(offset);
							}

						}
					}
				}


			}
			

			map<int, vector<int> >::iterator tmp_it, tmp_it2;
			tmp_it2 = LongReadContigIndex.CTG2LR_2.begin();
			for (tmp_it = LongReadContigIndex.CTG2LR.begin(); tmp_it != LongReadContigIndex.CTG2LR.end(); ++tmp_it)
			{
				sort(tmp_it->second.begin(), tmp_it->second.end());
				sort(tmp_it2->second.begin(), tmp_it2->second.end());
				//LongReadContigIndex.layout[tmp_it->second[tmp_it->second.size() / 2]].push_back(tmp_it->first);//ctg_no
				//LongReadContigIndex.layout[tmp_it->second[tmp_it->second.size() / 2]].push_back(tmp_it->second.size());//cov
				LongReadContigIndex.layout[tmp_it->second[tmp_it->second.size() / 2]].push_back(tmp_it->first);//ctg_no
				LongReadContigIndex.layout[tmp_it->second[tmp_it->second.size() / 2]].push_back(tmp_it->second.size());//cov
				LongReadContigIndex.layout[tmp_it->second[tmp_it->second.size() / 2]].push_back(tmp_it2->second[tmp_it2->second.size() / 2]);//coord2
				tmp_it2++;
			}

			int contig_matches = 0;

			for (tmp_it = LongReadContigIndex.layout.begin(); tmp_it != LongReadContigIndex.layout.end(); ++tmp_it)
			{

				if (tmp_it->second[1] > 0)
				{
					contig_matches++;
					//LongReadContigIndex_info << tmp_it->first << ", " << tmp_it->second[0] << ", " << tmp_it->second[1] << endl;
				}
			}
			

			if (contig_matches > 1)
			{
				tag[0] = '>';
				out_reads_info << tag << endl;

				for (tmp_it = LongReadContigIndex.layout.begin(); tmp_it != LongReadContigIndex.layout.end(); ++tmp_it)
				{

					if (tmp_it->second[1] > 0)
					{
						out_reads_info << tmp_it->first << ", " << tmp_it->second[0] << ", " << tmp_it->second[1] << ", " << tmp_it->second[2] <<endl;
					}
				}
				out_reads << tag << endl;
				out_reads << read_str << endl;
			}
			
			vector<int> final_ctgs_vt,final_cod_vt,overlap_vt,final_offset_vt;
			if(ctgs_vt.size()>1)
			{
				final_ctgs_vt.push_back(ctgs_vt[0]);
				final_cod_vt.push_back(cod_vt[0]);
				final_offset_vt.push_back(offset_vt[0]);
				for (int i=1;i<ctgs_vt.size();++i)
				{

					//split the cases for easier debugging
					if(ctgs_vt[i]!=ctgs_vt[i-1])
					{
						final_ctgs_vt.push_back(ctgs_vt[i]);
						final_cod_vt.push_back(cod_vt[i]);
						final_offset_vt.push_back(offset_vt[i]);
						overlap_vt.push_back(pos_vt[i]-pos_vt[i-1]);//
					}
					else
					{
						//ctgs_vt[i]==ctgs_vt[i-1]
						if((ctgs_vt[i]>0&&cod_vt[i]<cod_vt[i-1])||(ctgs_vt[i]<0&&cod_vt[i]>cod_vt[i-1]))
						{
							final_ctgs_vt.push_back(ctgs_vt[i]);
							final_cod_vt.push_back(cod_vt[i]);
							final_offset_vt.push_back(offset_vt[i]);
							overlap_vt.push_back(pos_vt[i]-pos_vt[i-1]);//
						}
					}
					
				}
			}

			if(final_ctgs_vt.size()>1)
			{
				output_current_read=1;
			}

			if(0)//output_current_read)
			{
				tag[0] = '>';
				//out_reads << tag << endl;
				//out_reads << read_str << endl;
				out_reads_info << tag << endl;
				
				for (int i = 0; i<final_ctgs_vt.size(); ++i)
				{
					out_reads_info << final_offset_vt[i]<<", " << final_ctgs_vt[i] << "," << contigs_info->contig_sz_vt[abs(final_ctgs_vt[i])] << ", " << endl;
				}

				
			}


			map<int, KmerInContig>::iterator it1,it2;

			for (it1 = LongReadContigIndex.LR2CTG.begin(); it1 != LongReadContigIndex.LR2CTG.end(); ++it1)
			{
				
				it2 = it1;
				it2++;
				if (it2 == LongReadContigIndex.LR2CTG.end())
				{
					break;
				}
				if (it2->second.contig_no != it1->second.contig_no)
				{
					bool flip1, flip2;
					flip1 = it1->second.flip;
					flip2 = it2->second.flip;
					int pos1 = it1->second.pos;
					int pos2 = it2->second.pos;
					int extra_bases1, extra_bases2;
					int contig1 = it1->second.contig_no;
					int contig2 = it2->second.contig_no;
					
					if (flip1)
					{
						contig1 = -contig1;
					}
					if (flip2)
					{
						contig2 = -contig2;
					}
					int dist = it2->first - it1->first;
					if (flip1 == 0)
					{
						extra_bases1 = contigs_info->contig_sz_vt[abs(contig1)]-pos1-1;
					}
					else
					{
						extra_bases1 = pos1;
					}

					if (flip2==0)
					{
						extra_bases2 = pos2;
					}
					else
					{
						extra_bases2 = contigs_info->contig_sz_vt[abs(contig2)] - pos2-1;
					}
					int extra_bases = extra_bases1 + extra_bases2;
					dist = dist - extra_bases - 1;
					
					string bridge;
					if (dist > 0)
					{
						bridge = read_str.substr(it1->first+extra_bases1,dist);
					}

					vector<int> temp_vec1, temp_vec2;
					temp_vec1.push_back(contig1);
					temp_vec1.push_back(contig2);
					temp_vec1.push_back(dist);

					temp_vec2.push_back(-contig2);
					temp_vec2.push_back(-contig1);
					temp_vec2.push_back(dist);

					bool take_rc = 0;
					for (int c = 0; c < 2; ++c)
					{
						if (temp_vec2[c] < temp_vec1[c])
						{
							take_rc = 1;
							break;
						}
						if (temp_vec2[c] > temp_vec1[c])
						{
							break;
						}

					}
					if (take_rc)
					{
						temp_vec1 = temp_vec2;
						if (bridge.size() > 0)
						{
							reverse(bridge.begin(),bridge.end());
							complement_str(bridge);
						}
					}

					temp_map[temp_vec1][bridge]++;

					
				}

			}
		
		}

	}


	map<vector<int>, map<string,int> >::iterator it;
	for (it = temp_map.begin(); it != temp_map.end(); ++it)
	{
		int contig1 = it->first[0];
		int contig2 = it->first[1];

		if (contig1 != 0 && contig2 != 0)
		{
			string best_str;
			int max_cov = 0;
			for (map<string, int>::iterator tmp_it = it->second.begin(); tmp_it != it->second.end(); ++tmp_it)
			{
				if (tmp_it->second > max_cov)
				{
					max_cov = tmp_it->second;
					best_str = tmp_it->first;
				}
			}

			if (contig1 > 0)
			{
				out_refined_graph << contig1 << " " << contig2;
			}
			else
			{
				out_refined_graph << contig1 << " " << -contig2;
			}
			
			out_refined_graph << " " << it->first[2] << " " << max_cov << " " << best_str<< endl;
		}
	}

}

void CollectingNonContainedReadsSlow0(struct hashtable0 *ht,struct hashtable0 *merge_ht, int K_size,vector<string>& filenames_vt, contigs_info * contigs_info,string ContigFilename)
{

	time_t beg_time,read_time;
	string in_fname=ContigFilename;
	ifstream Contigs_in(in_fname.c_str());
	size_t num_Contigs=0;
	string tag,s,kmer,str,seq_s;
	uint64_t f_seq,hv;
	struct kmer_t2 f_seq_t2;
	size_t hash_idx;
	bool found;
	bool flip_1,flip_2,flip_0;
	size_t ht_sz;
	size_t numReads=0;

	ofstream out_refined_graph("RefinedContigGraph.txt");
	map<vector<int>, map<string,int> > temp_map;

	int Kmer_arr_sz=K_size/32+1;
	int rem1=K_size%32;
	if(rem1==0)
	{Kmer_arr_sz--;}

	int boundary=0,removed=0,bridge=0;

	ht_sz=ht->ht_sz;
	
	
	cout<<"Contigs remapping."<<endl;
	if(ContigFilename=="Contigs.txt")
	{
		ContigsRemapping0(ht, K_size, contigs_info,ContigFilename,0);
		
	}

	// AppendMergeHT0(ht, merge_ht,Kmer_arr_sz);
	//cout<<"Collecting informative reads."<<endl;
	//BuildContigAdjacency0(ht,contigs_info, K_size,ContigFilename);
	
	time(&beg_time);

	int64_t dist_sum=0,p_cnt=0;
	int mean_dist=0;
	int lib_no=0;
	
	for(size_t ii=0;ii<filenames_vt.size();ii++)
	{
		lib_no++;
	
		dist_sum=0,p_cnt=0;
		cout<<"Processing library: "<<ii<<endl;
		ifstream in_reads;
		in_reads.open(filenames_vt[ii].c_str());
		string file = filenames_vt[ii];
		int file_beg = 0;
		for (int c = 0; c < file.size(); ++c)
		{
			if (file[c] == '\\' || file[c] == '/')
			{
				file_beg = c;
			}
		}
		file = file.substr(file_beg + 1, file.size());

		string Name1 = "NonContainedReadsFrom_" + file;
		string Name2 = "ReadsInfoFrom_" + file;
		
		ofstream out_reads(Name1.c_str());
		ofstream out_reads_info(Name2.c_str());;
		cout <<"Creating: "<< Name1 << endl;
		
		struct read_t Read1,Read2;

		Read1.read_bits =(uint64_t*) malloc(1000000/4+100);
		Read2.read_bits =(uint64_t*) malloc(1000000/4+100);

		int nLines1=0,nLines2=0;


		bool fq_flag=0;


		getline(in_reads,str);
		if(fq_flag==0&&str[0]=='@')
		{
			fq_flag=1;	
		}
		in_reads.close();

		in_reads.clear();
		in_reads.open(filenames_vt[ii].c_str());

		bool read_success=0;

		read_success=1;

		string tag,qs,n_tag;
		string QS_s;

		while(read_success)
		{
			if(fq_flag)
			{
				read_success=get_a_fastq_read(in_reads,tag,seq_s,QS_s);
					
			}
			else
			{
				read_success=get_a_fasta_read(in_reads,tag,seq_s,n_tag);
			
			}	
			if(read_success==0)
			{break;}

			/*
			cout << tag << endl;
			if (tag == "@SRR826442.57 57 length = 150")
			{
				cout << tag << endl;
			}
			*/
				
			int seq_sz=seq_s.size();
			int readLen = seq_sz;
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
			
			if(bad_flag)
			{continue;}
							

			bad_flag=0;
					
			numReads++;


						
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
					


			Init_Read(seq_s,Read1);
			string read_str=seq_s;
			seq_s.clear();
			LongReadContigIndex LongReadContigIndex;

			int contig_no=-1;
			bool output_current_read=0;
			vector<int> ctgs_vt, cod_vt, pos_vt, offset_vt;
			bool bubble_flag=0;
			for(int i=0;i<Read1.readLen-K_size+1;++i )
			{
				
					
				uint64_t bits1[100],f_seq[100];
				get_sub_arr(Read1.read_bits,Read1.readLen,i,K_size,bits1);
				memcpy(f_seq,bits1,sizeof(uint64_t)*Kmer_arr_sz);
				get_rev_comp_seq_arr(f_seq,K_size,Kmer_arr_sz);
				flip_0=0;
				if(uint64_t_cmp(bits1,f_seq,Kmer_arr_sz)>0)
				{
					memcpy(bits1,f_seq,sizeof(uint64_t)*Kmer_arr_sz);
					
					flip_0=1;
				}
				
				hv=MurmurHash64A(bits1,sizeof(uint64_t)*Kmer_arr_sz,0);

				hash_idx=(size_t) (hv%ht_sz);
				struct bucket0 **ptr1;
				ptr1= &(ht->store_pos[hash_idx]);
				found=look_up_in_a_list0(bits1,&ptr1,Kmer_arr_sz);
				if(found)
				{
					int cod = (*ptr1)->kmer_info.cod;
					int contig_no = (*ptr1)->kmer_info.contig_no;

					if (contig_no == 0 || (cod > readLen + 50 && cod < (contigs_info->contig_sz_vt[abs(contig_no)] - readLen - 50)))
					{

						continue;//contained hit
					}

					LongReadContigIndex.LR2CTG[i].contig_no = contig_no;
					LongReadContigIndex.LR2CTG[i].pos = (*(ptr1))->kmer_info.cod;
					LongReadContigIndex.LR2CTG[i].flip = flip_0 ^ (*(ptr1))->kmer_info.flip;
					

					int offset = 0;
					if ((*ptr1)->kmer_info.removed == 0 && (*ptr1)->kmer_info.contig_no > 0 && (*ptr1)->kmer_info.repeat == 0)
					{
						if ((flip_0 ^ (*ptr1)->kmer_info.flip) == 0)
						{
							ctgs_vt.push_back((*ptr1)->kmer_info.contig_no);
							offset = i - (*ptr1)->kmer_info.cod + contigs_info->contig_sz_vt[(*(ptr1))->kmer_info.contig_no] / 2;
							LongReadContigIndex.CTG2LR_2[(*(ptr1))->kmer_info.contig_no].push_back(offset);
							LongReadContigIndex.CTG2LR[(*(ptr1))->kmer_info.contig_no].push_back(i);

						}
						else
						{
							ctgs_vt.push_back(-(*ptr1)->kmer_info.contig_no);;
							offset = i + (*ptr1)->kmer_info.cod - contigs_info->contig_sz_vt[(*(ptr1))->kmer_info.contig_no] / 2 + (K_size + 1) / 2;
							LongReadContigIndex.CTG2LR_2[-(*(ptr1))->kmer_info.contig_no].push_back(offset);
							LongReadContigIndex.CTG2LR[-(*(ptr1))->kmer_info.contig_no].push_back(i);

						}
						cod_vt.push_back((*ptr1)->kmer_info.cod);
						pos_vt.push_back(i);
						offset_vt.push_back(offset);
					}
				}
				


			}


			map<int, vector<int> >::iterator tmp_it, tmp_it2;
			tmp_it2 = LongReadContigIndex.CTG2LR_2.begin();
			for (tmp_it = LongReadContigIndex.CTG2LR.begin(); tmp_it != LongReadContigIndex.CTG2LR.end(); ++tmp_it)
			{
				sort(tmp_it->second.begin(), tmp_it->second.end());
				sort(tmp_it2->second.begin(), tmp_it2->second.end());
				//LongReadContigIndex.layout[tmp_it->second[tmp_it->second.size() / 2]].push_back(tmp_it->first);//ctg_no
				//LongReadContigIndex.layout[tmp_it->second[tmp_it->second.size() / 2]].push_back(tmp_it->second.size());//cov
				LongReadContigIndex.layout[tmp_it->second[tmp_it->second.size() / 2]].push_back(tmp_it->first);//ctg_no
				LongReadContigIndex.layout[tmp_it->second[tmp_it->second.size() / 2]].push_back(tmp_it->second.size());//cov
				LongReadContigIndex.layout[tmp_it->second[tmp_it->second.size() / 2]].push_back(tmp_it2->second[tmp_it2->second.size() / 2]);//coord2
				tmp_it2++;
			}

			int contig_matches = 0;

			for (tmp_it = LongReadContigIndex.layout.begin(); tmp_it != LongReadContigIndex.layout.end(); ++tmp_it)
			{

				if (tmp_it->second[1] > 0)
				{
					contig_matches++;
					//LongReadContigIndex_info << tmp_it->first << ", " << tmp_it->second[0] << ", " << tmp_it->second[1] << endl;
				}
			}


			if (contig_matches > 1)
			{
				tag[0] = '>';
				out_reads_info << tag << endl;

				for (tmp_it = LongReadContigIndex.layout.begin(); tmp_it != LongReadContigIndex.layout.end(); ++tmp_it)
				{

					if (tmp_it->second[1] > 0)
					{

						out_reads_info << tmp_it->first << ", " << tmp_it->second[0] << ", " << tmp_it->second[1] << ", " << tmp_it->second[2] << endl;
					}
				}
				out_reads << tag << endl;
				out_reads << read_str << endl;
			}
			vector<int> final_ctgs_vt, final_cod_vt, overlap_vt, final_offset_vt;
			if (ctgs_vt.size()>1)
			{
				final_ctgs_vt.push_back(ctgs_vt[0]);
				final_cod_vt.push_back(cod_vt[0]);
				final_offset_vt.push_back(offset_vt[0]);
				for (int i = 1; i<ctgs_vt.size(); ++i)
				{

					//split the cases for easier debugging
					if (ctgs_vt[i] != ctgs_vt[i - 1])
					{
						final_ctgs_vt.push_back(ctgs_vt[i]);
						final_cod_vt.push_back(cod_vt[i]);
						final_offset_vt.push_back(offset_vt[i]);
						overlap_vt.push_back(pos_vt[i] - pos_vt[i - 1]);//
					}
					else
					{
						//ctgs_vt[i]==ctgs_vt[i-1]
						if ((ctgs_vt[i]>0 && cod_vt[i]<cod_vt[i - 1]) || (ctgs_vt[i]<0 && cod_vt[i]>cod_vt[i - 1]))
						{
							final_ctgs_vt.push_back(ctgs_vt[i]);
							final_cod_vt.push_back(cod_vt[i]);
							final_offset_vt.push_back(offset_vt[i]);
							overlap_vt.push_back(pos_vt[i] - pos_vt[i - 1]);//
						}
					}

				}
			}

			if(final_ctgs_vt.size()>1)
			{
				output_current_read=1;
			}

			if(0)//output_current_read)
			{

				tag[0] = '>';
				//out_reads << tag << endl;
				//out_reads << read_str << endl;
				out_reads_info << tag << endl;

				for (int i = 0; i<final_ctgs_vt.size(); ++i)
				{
					out_reads_info << final_offset_vt[i] << ", " << final_ctgs_vt[i] << "," << contigs_info->contig_sz_vt[abs(final_ctgs_vt[i])] << ", " << endl;
				}
			}


			map<int, KmerInContig>::iterator it1, it2;

			for (it1 = LongReadContigIndex.LR2CTG.begin(); it1 != LongReadContigIndex.LR2CTG.end(); ++it1)
			{

				it2 = it1;
				it2++;
				if (it2 == LongReadContigIndex.LR2CTG.end())
				{
					break;
				}
				if (it2->second.contig_no != it1->second.contig_no)
				{
					bool flip1, flip2;
					flip1 = it1->second.flip;
					flip2 = it2->second.flip;
					int pos1 = it1->second.pos;
					int pos2 = it2->second.pos;
					int extra_bases1, extra_bases2;
					int contig1 = it1->second.contig_no;
					int contig2 = it2->second.contig_no;

					if (flip1)
					{
						contig1 = -contig1;
					}
					if (flip2)
					{
						contig2 = -contig2;
					}
					int dist = it2->first - it1->first;
					if (flip1 == 0)
					{
						extra_bases1 = contigs_info->contig_sz_vt[abs(contig1)] - pos1-1;
					}
					else
					{
						extra_bases1 = pos1;
					}

					if (flip2 == 0)
					{
						extra_bases2 = pos2;
					}
					else
					{
						extra_bases2 = contigs_info->contig_sz_vt[abs(contig2)] - pos2-1;
					}
					int extra_bases = extra_bases1 + extra_bases2;
					dist = dist - extra_bases - 1;

					string bridge;
					if (dist > 0)
					{
						bridge = read_str.substr(it1->first + extra_bases1, dist);
					}


					vector<int> temp_vec1,temp_vec2;
					temp_vec1.push_back(contig1);
					temp_vec1.push_back(contig2);
					temp_vec1.push_back(dist);
					
					temp_vec2.push_back(-contig2);
					temp_vec2.push_back(-contig1);
					temp_vec2.push_back(dist);

					bool take_rc = 0;
					for (int c = 0; c < 2; ++c)
					{
						if (temp_vec2[c] < temp_vec1[c])
						{
							take_rc = 1;
							break;
						}
						if (temp_vec2[c] > temp_vec1[c])
						{
							break;
						}

					}
					if (take_rc)
					{
						temp_vec1 = temp_vec2;
						if (bridge.size() > 0)
						{
							reverse(bridge.begin(), bridge.end());
							complement_str(bridge);
						}
					}

					temp_map[temp_vec1][bridge]++;
					
				}

			}
		}

	}


	map<vector<int>, map<string, int> >::iterator it;
	for (it = temp_map.begin(); it != temp_map.end(); ++it)
	{
		int contig1 = it->first[0];
		int contig2 = it->first[1];

		if (contig1 != 0 && contig2 != 0)
		{
			string best_str;
			int max_cov = 0;
			for (map<string, int>::iterator tmp_it = it->second.begin(); tmp_it != it->second.end(); ++tmp_it)
			{
				if (tmp_it->second > max_cov)
				{
					max_cov = tmp_it->second;
					best_str = tmp_it->first;
				}
			}

			if (contig1 > 0)
			{
				out_refined_graph << contig1 << " " << contig2;
			}
			else
			{
				out_refined_graph << contig1 << " " << -contig2;
			}

			out_refined_graph << " " << it->first[2] << " " << max_cov << " " << best_str << endl;
		}
	}
}



void CollectingNonContainedPairsSlow(struct hashtable *ht1, struct hashtable *merge_ht1, struct hashtable2 *ht2, struct hashtable2 *merge_ht2, int K_size, vector<string>& filenames_vt, contigs_info * contigs_info, string ContigFilename)
{

	time_t beg_time, read_time;
	string in_fname = ContigFilename;
	ifstream Contigs_in(in_fname.c_str());
	size_t num_Contigs = 0;

	uint64_t f_seq, hv;	
	size_t hash_idx;
	bool found;
	bool flip_1, flip_2, flip_0;
	size_t ht_sz;
	size_t numReads = 0;
	ofstream out_pairs("NonContainedPairs.txt"), out_pairs_info("NonContainedPairs_info.txt");;
	ofstream out_contained_pairs_info("ContainedPairs_info.txt");
	uint64_t insert_size_sum = 0;
	uint64_t contained_pairs_cnt=0;
	ofstream out_refined_graph("RefinedContigGraph.txt");
	map<vector<int>, map<string,int > > temp_map;
	
	ht_sz = ht1->ht_sz;
		
	cout << "Contigs remapping." << endl;
	
	ContigsRemapping(ht1, ht2, K_size, contigs_info, ContigFilename, 0);	
	time(&beg_time);

	int64_t dist_sum = 0, p_cnt = 0;
	int mean_dist = 0;
	int lib_no = 0;

	for (size_t ii = 0; ii<filenames_vt.size(); ii+=2)
	{
		lib_no++;
		dist_sum = 0, p_cnt = 0;
		cout << "Processing libraries: " << ii + 1 << " and " << ii + 2<< endl;
		ifstream in_pair1,in_pair2;
		in_pair1.open(filenames_vt[ii].c_str());
		in_pair2.open(filenames_vt[ii + 1].c_str());

		struct read_t Read1, Read2;
		Read1.read_bits = (uint64_t*)malloc(1000000 / 4 + 100);
		Read2.read_bits = (uint64_t*)malloc(1000000 / 4 + 100);

		int nLines1 = 0, nLines2 = 0;
		string tag1, qs1, n_tag1, tag2, qs2, n_tag2, seq_s1, seq_s2;
		string QS_s1, QS_s2;

		bool fq_flag = 0;

		getline(in_pair1, tag1);
		if (fq_flag == 0 && tag1[0] == '@')
		{
			fq_flag = 1;
		}
		in_pair1.close();

		in_pair1.clear();
		in_pair1.open(filenames_vt[ii].c_str());

		bool read_success1 = 1, read_success2 = 1;


		
		while (read_success1&&read_success2)
		{
			if (fq_flag)
			{
				read_success1 = get_a_fastq_read(in_pair1, tag1, seq_s1, QS_s1);
				read_success2 = get_a_fastq_read(in_pair2, tag2, seq_s2, QS_s2);

			}
			else
			{
				read_success1 = get_a_fasta_read(in_pair1, tag1, seq_s1, n_tag1);
				read_success2 = get_a_fasta_read(in_pair2, tag2, seq_s2, n_tag2);

			}
			if (read_success1 == 0 || read_success2 == 0)
			{
				break;
			}

			int readLen1 = seq_s1.size();

			int readLen2 = seq_s2.size();


			bool good_read1 = basic_quality_check(seq_s1);
			bool good_read2 = basic_quality_check(seq_s2);

		

			if (readLen1<K_size)
			{
				good_read1 = 0;
			}
			if (readLen2<K_size)
			{
				good_read2 = 0;
			}


			if (good_read1 == 0 || good_read2 == 0)
			{
				continue;
			}

			Init_Read(seq_s1, Read1);
			Init_Read(seq_s2, Read2);

		
			
			LongReadContigIndex LongReadContigIndex1, LongReadContigIndex2;

			uint64_t bits1;
			int contig_no = -1;
			bool output_current_pair= 0;
			
			for (int i = 0; i<Read1.readLen - K_size + 1; ++i)
			{
				
				get_sub_arr(Read1.read_bits, Read1.readLen, i, K_size, &bits1);
				f_seq = get_rev_comp_seq(bits1, K_size);
				flip_0 = 0;
				if (bits1>f_seq)
				{
					bits1 = f_seq;
					flip_0 = 1;
				}

				hv = MurmurHash64A(&bits1, sizeof(bits1), 0);
				hash_idx = (size_t)(hv%ht_sz);
				struct bucket **ptr1;
				ptr1 = &(ht1->store_pos[hash_idx]);
				found = look_up_in_a_list(bits1, &ptr1);
				if (found)
				{
					int cod = (*ptr1)->kmer_info.cod;
					int contig_no = (*ptr1)->kmer_info.contig_no;
					if (contig_no == 0)
					{
						continue;
					}
					if (cod > readLen1 + 50 && cod < (contigs_info->contig_sz_vt[abs(contig_no)] - readLen1 - 50))
					{
						//continue;//contained hit
					}

					LongReadContigIndex1.LR2CTG[i].contig_no = contig_no;
					LongReadContigIndex1.LR2CTG[i].pos = (*(ptr1))->kmer_info.cod;
					LongReadContigIndex1.LR2CTG[i].flip = flip_0 ^ (*(ptr1))->kmer_info.flip;

					int offset = 0;
					if ((*ptr1)->kmer_info.removed == 0 && (*ptr1)->kmer_info.contig_no>0 && (*ptr1)->kmer_info.repeat == 0)
					{
						if ((flip_0 ^ (*ptr1)->kmer_info.flip) == 0)
						{
							offset = i - (*ptr1)->kmer_info.cod + contigs_info->contig_sz_vt[(*(ptr1))->kmer_info.contig_no] / 2;
							//cout << offset << endl;
							LongReadContigIndex1.CTG2LR[(*(ptr1))->kmer_info.contig_no].push_back(i);
							LongReadContigIndex1.CTG2LR_2[(*(ptr1))->kmer_info.contig_no].push_back(offset);

						}
						else
						{
							offset = i + (*ptr1)->kmer_info.cod - contigs_info->contig_sz_vt[(*(ptr1))->kmer_info.contig_no] / 2 + (K_size + 1) / 2;
							//cout << offset << endl;
							LongReadContigIndex1.CTG2LR[-(*(ptr1))->kmer_info.contig_no].push_back(i);
							LongReadContigIndex1.CTG2LR_2[-(*(ptr1))->kmer_info.contig_no].push_back(offset);

						}
						
					}

				}
				
			}



			for (int i = 0; i<Read2.readLen - K_size + 1; ++i)
			{

				get_sub_arr(Read2.read_bits, Read2.readLen, i, K_size, &bits1);
				f_seq = get_rev_comp_seq(bits1, K_size);
				flip_0 = 0;
				if (bits1>f_seq)
				{
					bits1 = f_seq;
					flip_0 = 1;
				}

				hv = MurmurHash64A(&bits1, sizeof(bits1), 0);
				hash_idx = (size_t)(hv%ht_sz);
				struct bucket **ptr1;
				ptr1 = &(ht1->store_pos[hash_idx]);
				found = look_up_in_a_list(bits1, &ptr1);
				if (found)
				{
					int cod = (*ptr1)->kmer_info.cod;
					int contig_no = (*ptr1)->kmer_info.contig_no;
					if (contig_no == 0)
					{
						continue;
					}
					if (cod > readLen2 + 50 && cod < (contigs_info->contig_sz_vt[abs(contig_no)] - readLen2 - 50))
					{
						//continue;//contained hit
					}

					LongReadContigIndex2.LR2CTG[i].contig_no = contig_no;
					LongReadContigIndex2.LR2CTG[i].pos = (*(ptr1))->kmer_info.cod;
					LongReadContigIndex2.LR2CTG[i].flip = flip_0 ^ (*(ptr1))->kmer_info.flip;

					int offset = 0;
					if ((*ptr1)->kmer_info.removed == 0 && (*ptr1)->kmer_info.contig_no>0 && (*ptr1)->kmer_info.repeat == 0)
					{
						if ((flip_0 ^ (*ptr1)->kmer_info.flip) == 0)
						{
							offset = i - (*ptr1)->kmer_info.cod + contigs_info->contig_sz_vt[(*(ptr1))->kmer_info.contig_no] / 2;
							//cout << offset << endl;
							LongReadContigIndex2.CTG2LR[(*(ptr1))->kmer_info.contig_no].push_back(i);
							LongReadContigIndex2.CTG2LR_2[(*(ptr1))->kmer_info.contig_no].push_back(offset);

						}
						else
						{
							offset = i + (*ptr1)->kmer_info.cod - contigs_info->contig_sz_vt[(*(ptr1))->kmer_info.contig_no] / 2 + (K_size + 1) / 2;
							//cout << offset << endl;
							LongReadContigIndex2.CTG2LR[-(*(ptr1))->kmer_info.contig_no].push_back(i);
							LongReadContigIndex2.CTG2LR_2[-(*(ptr1))->kmer_info.contig_no].push_back(offset);

						}

					}

				}

			}



			map<int, vector<int> >::iterator tmp_it, tmp_it2;
			tmp_it2 = LongReadContigIndex1.CTG2LR_2.begin();
			for (tmp_it = LongReadContigIndex1.CTG2LR.begin(); tmp_it != LongReadContigIndex1.CTG2LR.end(); ++tmp_it)
			{
				sort(tmp_it->second.begin(), tmp_it->second.end());
				sort(tmp_it2->second.begin(), tmp_it2->second.end());
				LongReadContigIndex1.layout[tmp_it->second[tmp_it->second.size() / 2]].push_back(tmp_it->first);//ctg_no
				LongReadContigIndex1.layout[tmp_it->second[tmp_it->second.size() / 2]].push_back(tmp_it->second.size());//cov
				LongReadContigIndex1.layout[tmp_it->second[tmp_it->second.size() / 2]].push_back(tmp_it2->second[tmp_it2->second.size() / 2]);//coord2
				tmp_it2++;
			}


			tmp_it2 = LongReadContigIndex2.CTG2LR_2.begin();
			for (tmp_it = LongReadContigIndex2.CTG2LR.begin(); tmp_it != LongReadContigIndex2.CTG2LR.end(); ++tmp_it)
			{
				sort(tmp_it->second.begin(), tmp_it->second.end());
				sort(tmp_it2->second.begin(), tmp_it2->second.end());
				LongReadContigIndex2.layout[tmp_it->second[tmp_it->second.size() / 2]].push_back(tmp_it->first);//ctg_no
				LongReadContigIndex2.layout[tmp_it->second[tmp_it->second.size() / 2]].push_back(tmp_it->second.size());//cov
				LongReadContigIndex2.layout[tmp_it->second[tmp_it->second.size() / 2]].push_back(tmp_it2->second[tmp_it2->second.size() / 2]);//coord2
				tmp_it2++;
			}



			int contig_matches1 = 0, contig_matches2 = 0;
			int contig1=0, contig2 = 0;
			for (tmp_it = LongReadContigIndex1.layout.begin(); tmp_it != LongReadContigIndex1.layout.end(); ++tmp_it)
			{

				if (tmp_it->second[1] > 0)
				{
					contig_matches1++;

					contig1 = tmp_it->second[0];
				}
			}

			for (tmp_it = LongReadContigIndex2.layout.begin(); tmp_it != LongReadContigIndex2.layout.end(); ++tmp_it)
			{

				if (tmp_it->second[1] > 0)
				{
					contig_matches2++;
					contig2 = tmp_it->second[0];
				}
			}

			if (contig1 == 0 || contig2 == 0)
			{
				continue;
			}
			if (contig_matches1 == 1 && contig_matches2 == 1 && (abs(contig1) == abs(contig2)))//a bit buggy
			{
				if (contigs_info->contig_sz_vt[abs(contig1)] > 500)
				{
					insert_size_sum += abs(abs(LongReadContigIndex1.layout.begin()->second[2]) - abs(LongReadContigIndex2.layout.begin()->second[2]));
					contained_pairs_cnt++;
					if (contained_pairs_cnt % 1000 == 0)
					{
						out_contained_pairs_info << insert_size_sum / contained_pairs_cnt << endl;
					}
				}
//				out_contained_pairs_info << contig1 << " " << LongReadContigIndex1.layout.begin()->second[2];
	//			out_contained_pairs_info <<" "<< contig2 << " " << LongReadContigIndex2.layout.begin()->second[2];
		//		out_contained_pairs_info << endl;
			}

			if ((contig_matches1>0 && contig_matches2>0) && (contig_matches1 > 1 || contig_matches2>1 || (abs(contig1) != abs(contig2))))
			{
				tag1[0] = '>';
				out_pairs_info << tag1 << endl;

				for (tmp_it = LongReadContigIndex1.layout.begin(); tmp_it != LongReadContigIndex1.layout.end(); ++tmp_it)
				{

					if (tmp_it->second[1] > 0)
					{
						out_pairs_info << tmp_it->first << ", " << tmp_it->second[0] << ", " << tmp_it->second[1] << ", " << tmp_it->second[2] << endl;
					}
				}
				tag2[0] = '>';
				out_pairs_info << tag2 << endl;

				for (tmp_it = LongReadContigIndex2.layout.begin(); tmp_it != LongReadContigIndex2.layout.end(); ++tmp_it)
				{

					if (tmp_it->second[1] > 0)
					{
						out_pairs_info << tmp_it->first << ", " << tmp_it->second[0] << ", " << tmp_it->second[1] << ", " << tmp_it->second[2] << endl;
					}
				}
				out_pairs << tag1 << endl;
				out_pairs << seq_s1 << endl;
				out_pairs << tag2 << endl;
				out_pairs << seq_s2 << endl;

			}
			

			
			string read_str = seq_s1;
			
			for (int pair = 1; pair <= 2; ++pair)
			{
				if (pair == 2)
				{
					LongReadContigIndex1 = LongReadContigIndex2;
					read_str = seq_s2;
				}
				map<int, KmerInContig>::iterator it1, it2;

				for (it1 = LongReadContigIndex1.LR2CTG.begin(); it1 != LongReadContigIndex1.LR2CTG.end(); ++it1)
				{

					it2 = it1;
					it2++;
					if (it2 == LongReadContigIndex1.LR2CTG.end())
					{
						break;
					}
					if (it2->second.contig_no != it1->second.contig_no)
					{
						bool flip1, flip2;
						flip1 = it1->second.flip;
						flip2 = it2->second.flip;
						int pos1 = it1->second.pos;
						int pos2 = it2->second.pos;
						int extra_bases1, extra_bases2;
						int contig1 = it1->second.contig_no;
						int contig2 = it2->second.contig_no;

						if (flip1)
						{
							contig1 = -contig1;
						}
						if (flip2)
						{
							contig2 = -contig2;
						}
						int dist = it2->first - it1->first;
						if (flip1 == 0)
						{
							extra_bases1 = contigs_info->contig_sz_vt[abs(contig1)] - pos1-1;
						}
						else
						{
							extra_bases1 = pos1;
						}

						if (flip2 == 0)
						{
							extra_bases2 = pos2;
						}
						else
						{
							extra_bases2 = contigs_info->contig_sz_vt[abs(contig2)] - pos2-1;
						}
						int extra_bases = extra_bases1 + extra_bases2;
						dist = dist - extra_bases - 1;
						string bridge;
						if (dist > 0)
						{
							bridge = read_str.substr(it1->first + extra_bases1, dist);
						}
						vector<int> temp_vec1, temp_vec2;
						temp_vec1.push_back(contig1);
						temp_vec1.push_back(contig2);
						temp_vec1.push_back(dist);

						temp_vec2.push_back(-contig2);
						temp_vec2.push_back(-contig1);
						temp_vec2.push_back(dist);

						bool take_rc = 0;
						for (int c = 0; c < 2; ++c)
						{
							if (temp_vec2[c] < temp_vec1[c])
							{
								take_rc = 1;
								break;
							}
							if (temp_vec2[c] > temp_vec1[c])
							{
								break;
							}

						}
						if (take_rc)
						{
							temp_vec1 = temp_vec2;
							if (bridge.size() > 0)
							{
								reverse(bridge.begin(), bridge.end());
								complement_str(bridge);
							}
						}

						temp_map[temp_vec1][bridge]++;



					}

				}
			}

		}

	}


	map<vector<int>, map<string, int> >::iterator it;
	for (it = temp_map.begin(); it != temp_map.end(); ++it)
	{

		int contig1 = it->first[0];
		int contig2 = it->first[1];

		if (contig1 != 0 && contig2 != 0)
		{

			string best_str;
			int max_cov = 0;
			for (map<string, int>::iterator tmp_it = it->second.begin(); tmp_it != it->second.end(); ++tmp_it)
			{
				if (tmp_it->second > max_cov)
				{
					max_cov = tmp_it->second;
					best_str = tmp_it->first;
				}
			}

			if (contig1 > 0)
			{
				out_refined_graph << contig1 << " " << contig2;
			}
			else
			{
				out_refined_graph << contig1 << " " << -contig2;
				if (best_str.size() > 0)
				{
					reverse(best_str.begin(), best_str.end());
					complement_str(best_str);
				}
			}
			                               //dist                 //cov           //bridge
			out_refined_graph << " " << it->first[2] << " " << max_cov << " " << best_str << endl;
		}
	}

}

void CollectingNonContainedPairsSlow0(struct hashtable0 *ht, struct hashtable0 *merge_ht, int K_size, vector<string>& filenames_vt, contigs_info * contigs_info, string ContigFilename)
{

	time_t beg_time, read_time;
	string in_fname = ContigFilename;
	ifstream Contigs_in(in_fname.c_str());
	size_t num_Contigs = 0;

	uint64_t f_seq, hv;

	size_t hash_idx;
	bool found;
	bool flip_1, flip_2, flip_0;
	size_t ht_sz;
	size_t numReads = 0;
	ofstream out_pairs("NonContainedPairs.txt"), out_pairs_info("NonContainedPairs_info.txt");;
	ofstream out_contained_pairs_info("ContainedPairs_info.txt");
	uint64_t insert_size_sum = 0;
	uint64_t contained_pairs_cnt = 0;
	ofstream out_refined_graph("RefinedContigGraph.txt");
	map<vector<int>, map<string,int> > temp_map;

	int Kmer_arr_sz = K_size / 32 + 1;
	int rem1 = K_size % 32;
	if (rem1 == 0)
	{
		Kmer_arr_sz--;
	}

	int boundary = 0, removed = 0, bridge = 0;

	ht_sz = ht->ht_sz;


	cout << "Contigs remapping." << endl;

	ContigsRemapping0(ht, K_size, contigs_info, ContigFilename, 0);


	// AppendMergeHT0(ht, merge_ht,Kmer_arr_sz);
	//cout<<"Collecting informative reads."<<endl;
	//BuildContigAdjacency0(ht,contigs_info, K_size,ContigFilename);

	time(&beg_time);

	int64_t dist_sum = 0, p_cnt = 0;
	int mean_dist = 0;
	int lib_no = 0;

	for (size_t ii = 0; ii<filenames_vt.size(); ii+=2)
	{
		lib_no++;

		dist_sum = 0, p_cnt = 0;

		cout << "Processing libraries: " << ii + 1 << " and " << ii + 2 << endl;
		ifstream in_pair1, in_pair2;
		in_pair1.open(filenames_vt[ii].c_str());
		in_pair2.open(filenames_vt[ii + 1].c_str());

		struct read_t Read1, Read2;

		Read1.read_bits = (uint64_t*)malloc(1000000 / 4 + 100);
		Read2.read_bits = (uint64_t*)malloc(1000000 / 4 + 100);

		int nLines1 = 0, nLines2 = 0;


		bool fq_flag = 0;


		string tag1, qs1, n_tag1, tag2, qs2, n_tag2, seq_s1, seq_s2;
		string QS_s1, QS_s2;

		getline(in_pair1, tag1);
		if (fq_flag == 0 && tag1[0] == '@')
		{
			fq_flag = 1;
		}
		in_pair1.close();

		in_pair1.clear();
		in_pair1.open(filenames_vt[ii].c_str());

		bool read_success1 = 1, read_success2 = 1;



		while (read_success1&&read_success2)
		{
			if (fq_flag)
			{
				read_success1 = get_a_fastq_read(in_pair1, tag1, seq_s1, QS_s1);
				read_success2 = get_a_fastq_read(in_pair2, tag2, seq_s2, QS_s2);

			}
			else
			{
				read_success1 = get_a_fasta_read(in_pair1, tag1, seq_s1, n_tag1);
				read_success2 = get_a_fasta_read(in_pair2, tag2, seq_s2, n_tag2);

			}
			if (read_success1 == 0 || read_success2 == 0)
			{
				break;
			}

			int readLen1 = seq_s1.size();

			int readLen2 = seq_s2.size();


			bool good_read1 = basic_quality_check(seq_s1);
			bool good_read2 = basic_quality_check(seq_s2);



			if (readLen1<K_size)
			{
				good_read1 = 0;
			}
			if (readLen2<K_size)
			{
				good_read2 = 0;
			}


			if (good_read1 == 0 || good_read2 == 0)
			{
				continue;
			}

			Init_Read(seq_s1, Read1);
			Init_Read(seq_s2, Read2);



			LongReadContigIndex LongReadContigIndex1, LongReadContigIndex2;



			for (int i = 0; i<Read1.readLen - K_size + 1; ++i)
			{


				uint64_t bits1[100], f_seq[100];
				get_sub_arr(Read1.read_bits, Read1.readLen, i, K_size, bits1);
				memcpy(f_seq, bits1, sizeof(uint64_t)*Kmer_arr_sz);
				get_rev_comp_seq_arr(f_seq, K_size, Kmer_arr_sz);
				flip_0 = 0;
				if (uint64_t_cmp(bits1, f_seq, Kmer_arr_sz)>0)
				{
					memcpy(bits1, f_seq, sizeof(uint64_t)*Kmer_arr_sz);

					flip_0 = 1;
				}

				hv = MurmurHash64A(bits1, sizeof(uint64_t)*Kmer_arr_sz, 0);

				hash_idx = (size_t)(hv%ht_sz);
				struct bucket0 **ptr1;
				ptr1 = &(ht->store_pos[hash_idx]);
				found = look_up_in_a_list0(bits1, &ptr1, Kmer_arr_sz);
				if (found)
				{
					int cod = (*ptr1)->kmer_info.cod;
					int contig_no = (*ptr1)->kmer_info.contig_no;
					if (contig_no == 0)
					{
						continue;
					}
					if (cod > readLen1 + 50 && cod < (contigs_info->contig_sz_vt[abs(contig_no)] - readLen1 - 50))
					{

						//continue;//contained hit
					}

					LongReadContigIndex1.LR2CTG[i].contig_no = contig_no;
					LongReadContigIndex1.LR2CTG[i].pos = (*(ptr1))->kmer_info.cod;
					LongReadContigIndex1.LR2CTG[i].flip = flip_0 ^ (*(ptr1))->kmer_info.flip;


					int offset = 0;
					if ((*ptr1)->kmer_info.removed == 0 && (*ptr1)->kmer_info.contig_no > 0 && (*ptr1)->kmer_info.repeat == 0)
					{
						if ((flip_0 ^ (*ptr1)->kmer_info.flip) == 0)
						{
							offset = i - (*ptr1)->kmer_info.cod + contigs_info->contig_sz_vt[(*(ptr1))->kmer_info.contig_no] / 2;
							LongReadContigIndex1.CTG2LR_2[(*(ptr1))->kmer_info.contig_no].push_back(offset);
							LongReadContigIndex1.CTG2LR[(*(ptr1))->kmer_info.contig_no].push_back(i);

						}
						else
						{
							offset = i + (*ptr1)->kmer_info.cod - contigs_info->contig_sz_vt[(*(ptr1))->kmer_info.contig_no] / 2 + (K_size + 1) / 2;
							LongReadContigIndex1.CTG2LR_2[-(*(ptr1))->kmer_info.contig_no].push_back(offset);
							LongReadContigIndex1.CTG2LR[-(*(ptr1))->kmer_info.contig_no].push_back(i);

						}
						
					}
				}



			}





			for (int i = 0; i<Read2.readLen - K_size + 1; ++i)
			{


				uint64_t bits1[100], f_seq[100];
				get_sub_arr(Read2.read_bits, Read2.readLen, i, K_size, bits1);
				memcpy(f_seq, bits1, sizeof(uint64_t)*Kmer_arr_sz);
				get_rev_comp_seq_arr(f_seq, K_size, Kmer_arr_sz);
				flip_0 = 0;
				if (uint64_t_cmp(bits1, f_seq, Kmer_arr_sz)>0)
				{
					memcpy(bits1, f_seq, sizeof(uint64_t)*Kmer_arr_sz);

					flip_0 = 1;
				}

				hv = MurmurHash64A(bits1, sizeof(uint64_t)*Kmer_arr_sz, 0);

				hash_idx = (size_t)(hv%ht_sz);
				struct bucket0 **ptr1;
				ptr1 = &(ht->store_pos[hash_idx]);
				found = look_up_in_a_list0(bits1, &ptr1, Kmer_arr_sz);
				if (found)
				{
					int cod = (*ptr1)->kmer_info.cod;
					int contig_no = (*ptr1)->kmer_info.contig_no;
					if (contig_no == 0)
					{
						continue;
					}
					if (cod > readLen2 + 50 && cod < (contigs_info->contig_sz_vt[abs(contig_no)] - readLen2 - 50))
					{

						//continue;//contained hit
					}

					LongReadContigIndex2.LR2CTG[i].contig_no = contig_no;
					LongReadContigIndex2.LR2CTG[i].pos = (*(ptr1))->kmer_info.cod;
					LongReadContigIndex2.LR2CTG[i].flip = flip_0 ^ (*(ptr1))->kmer_info.flip;


					int offset = 0;
					if ((*ptr1)->kmer_info.removed == 0 && (*ptr1)->kmer_info.contig_no > 0 && (*ptr1)->kmer_info.repeat == 0)
					{
						if ((flip_0 ^ (*ptr1)->kmer_info.flip) == 0)
						{
							offset = i - (*ptr1)->kmer_info.cod + contigs_info->contig_sz_vt[(*(ptr1))->kmer_info.contig_no] / 2;
							LongReadContigIndex2.CTG2LR_2[(*(ptr1))->kmer_info.contig_no].push_back(offset);
							LongReadContigIndex2.CTG2LR[(*(ptr1))->kmer_info.contig_no].push_back(i);

						}
						else
						{
							offset = i + (*ptr1)->kmer_info.cod - contigs_info->contig_sz_vt[(*(ptr1))->kmer_info.contig_no] / 2 + (K_size + 1) / 2;
							LongReadContigIndex2.CTG2LR_2[-(*(ptr1))->kmer_info.contig_no].push_back(offset);
							LongReadContigIndex2.CTG2LR[-(*(ptr1))->kmer_info.contig_no].push_back(i);

						}

					}
				}



			}


			
			map<int, vector<int> >::iterator tmp_it, tmp_it2;
			tmp_it2 = LongReadContigIndex1.CTG2LR_2.begin();
			for (tmp_it = LongReadContigIndex1.CTG2LR.begin(); tmp_it != LongReadContigIndex1.CTG2LR.end(); ++tmp_it)
			{
				sort(tmp_it->second.begin(), tmp_it->second.end());
				sort(tmp_it2->second.begin(), tmp_it2->second.end());
				LongReadContigIndex1.layout[tmp_it->second[tmp_it->second.size() / 2]].push_back(tmp_it->first);//ctg_no
				LongReadContigIndex1.layout[tmp_it->second[tmp_it->second.size() / 2]].push_back(tmp_it->second.size());//cov
				LongReadContigIndex1.layout[tmp_it->second[tmp_it->second.size() / 2]].push_back(tmp_it2->second[tmp_it2->second.size() / 2]);//coord2
				tmp_it2++;
			}


			tmp_it2 = LongReadContigIndex2.CTG2LR_2.begin();
			for (tmp_it = LongReadContigIndex2.CTG2LR.begin(); tmp_it != LongReadContigIndex2.CTG2LR.end(); ++tmp_it)
			{
				sort(tmp_it->second.begin(), tmp_it->second.end());
				sort(tmp_it2->second.begin(), tmp_it2->second.end());
				LongReadContigIndex2.layout[tmp_it->second[tmp_it->second.size() / 2]].push_back(tmp_it->first);//ctg_no
				LongReadContigIndex2.layout[tmp_it->second[tmp_it->second.size() / 2]].push_back(tmp_it->second.size());//cov
				LongReadContigIndex2.layout[tmp_it->second[tmp_it->second.size() / 2]].push_back(tmp_it2->second[tmp_it2->second.size() / 2]);//coord2
				tmp_it2++;
			}



			int contig_matches1 = 0, contig_matches2 = 0;
			int contig1 = 0, contig2 = 0;
			for (tmp_it = LongReadContigIndex1.layout.begin(); tmp_it != LongReadContigIndex1.layout.end(); ++tmp_it)
			{

				if (tmp_it->second[1] > 0)
				{
					contig_matches1++;

					contig1 = tmp_it->second[0];
				}
			}

			for (tmp_it = LongReadContigIndex2.layout.begin(); tmp_it != LongReadContigIndex2.layout.end(); ++tmp_it)
			{

				if (tmp_it->second[1] > 0)
				{
					contig_matches2++;
					contig2 = tmp_it->second[0];
				}
			}

			if (contig1 == 0 || contig2 == 0)
			{
				continue;
			}
			if (contig_matches1 == 1 && contig_matches2 == 1 && (abs(contig1) == abs(contig2)))//a bit buggy
			{
				if (contigs_info->contig_sz_vt[abs(contig1)] > 500)
				{
					insert_size_sum += abs(abs(LongReadContigIndex1.layout.begin()->second[2]) - abs(LongReadContigIndex2.layout.begin()->second[2]));
					contained_pairs_cnt++;
					if (contained_pairs_cnt % 1000 == 0)
					{
						out_contained_pairs_info << insert_size_sum / contained_pairs_cnt << endl;
					}
				}
				//				out_contained_pairs_info << contig1 << " " << LongReadContigIndex1.layout.begin()->second[2];
				//			out_contained_pairs_info <<" "<< contig2 << " " << LongReadContigIndex2.layout.begin()->second[2];
				//		out_contained_pairs_info << endl;
			}

			if ((contig_matches1>0&&contig_matches2>0)&&(contig_matches1 > 1 || contig_matches2>1 || (abs(contig1) != abs(contig2))))//a bit buggy
			{
				tag1[0] = '>';
				out_pairs_info << tag1 << endl;

				for (tmp_it = LongReadContigIndex1.layout.begin(); tmp_it != LongReadContigIndex1.layout.end(); ++tmp_it)
				{

					if (tmp_it->second[1] > 0)
					{
						out_pairs_info << tmp_it->first << ", " << tmp_it->second[0] << ", " << tmp_it->second[1] << ", " << tmp_it->second[2] << endl;
					}
				}
				tag2[0] = '>';
				out_pairs_info << tag2 << endl;

				for (tmp_it = LongReadContigIndex2.layout.begin(); tmp_it != LongReadContigIndex2.layout.end(); ++tmp_it)
				{

					if (tmp_it->second[1] > 0)
					{
						out_pairs_info << tmp_it->first << ", " << tmp_it->second[0] << ", " << tmp_it->second[1] << ", " << tmp_it->second[2] << endl;
					}
				}
				out_pairs << tag1 << endl;
				out_pairs << seq_s1 << endl;
				out_pairs << tag2 << endl;
				out_pairs << seq_s2 << endl;

			}


			string read_str=seq_s1;


			for (int pair = 1; pair <= 2; ++pair)
			{
				if (pair == 2)
				{
					read_str = seq_s2;
					LongReadContigIndex1 = LongReadContigIndex2;
				}
				map<int, KmerInContig>::iterator it1, it2;

				for (it1 = LongReadContigIndex1.LR2CTG.begin(); it1 != LongReadContigIndex1.LR2CTG.end(); ++it1)
				{

					it2 = it1;
					it2++;
					if (it2 == LongReadContigIndex1.LR2CTG.end())
					{
						break;
					}
					if (it2->second.contig_no != it1->second.contig_no)
					{
						bool flip1, flip2;
						flip1 = it1->second.flip;
						flip2 = it2->second.flip;
						int pos1 = it1->second.pos;
						int pos2 = it2->second.pos;
						int extra_bases1, extra_bases2;
						int contig1 = it1->second.contig_no;
						int contig2 = it2->second.contig_no;

						if (flip1)
						{
							contig1 = -contig1;
						}
						if (flip2)
						{
							contig2 = -contig2;
						}
						int dist = it2->first - it1->first;
						if (flip1 == 0)
						{
							extra_bases1 = contigs_info->contig_sz_vt[abs(contig1)] - pos1-1;
						}
						else
						{
							extra_bases1 = pos1;
						}

						if (flip2 == 0)
						{
							extra_bases2 = pos2;
						}
						else
						{
							extra_bases2 = contigs_info->contig_sz_vt[abs(contig2)] - pos2-1;
						}
						int extra_bases = extra_bases1 + extra_bases2;
						dist = dist - extra_bases - 1;
						string bridge;
						if (dist > 0)
						{
							bridge = read_str.substr(it1->first + extra_bases1, dist);
						}
						vector<int> temp_vec1, temp_vec2;
						temp_vec1.push_back(contig1);
						temp_vec1.push_back(contig2);
						temp_vec1.push_back(dist);

						temp_vec2.push_back(-contig2);
						temp_vec2.push_back(-contig1);
						temp_vec2.push_back(dist);

						bool take_rc = 0;
						for (int c = 0; c < 2; ++c)
						{
							if (temp_vec2[c] < temp_vec1[c])
							{
								take_rc = 1;
								break;
							}
							if (temp_vec2[c] > temp_vec1[c])
							{
								break;
							}

						}
						if (take_rc)
						{
							temp_vec1 = temp_vec2;
							if (bridge.size() > 0)
							{
								reverse(bridge.begin(), bridge.end());
								complement_str(bridge);
							}
						}

						temp_map[temp_vec1][bridge]++;


					}

				}
			}

		}

	}


	map<vector<int>, map<string, int> >::iterator it;
	for (it = temp_map.begin(); it != temp_map.end(); ++it)
	{
		int contig1 = it->first[0];
		int contig2 = it->first[1];

		if (contig1 != 0 && contig2 != 0)
		{
			string best_str;
			int max_cov = 0;
			for (map<string, int>::iterator tmp_it = it->second.begin(); tmp_it != it->second.end(); ++tmp_it)
			{
				if (tmp_it->second > max_cov)
				{
					max_cov = tmp_it->second;
					best_str = tmp_it->first;
				}
			}

			if (contig1 > 0)
			{
				out_refined_graph << contig1 << " " << contig2;
			}
			else
			{
				out_refined_graph << contig1 << " " << -contig2;
			}

			out_refined_graph << " " << it->first[2] << " " << max_cov << " " << best_str << endl;
		}
	}
}



struct BFS_reads_info
{
	int cov;
	int depth;
	int len;
	int last_read;
	vector<int> edge;
};



bool isSimplePath_read(reads_overlap_info *reads_overlap_info,int current_read,map<int,struct BFS_reads_info > &Visited_Path , map<int, int > &stacked_nodes)
{
	//return 1;
	int last_read=current_read;
	int dep=Visited_Path[current_read].depth;
	int node=last_read;
	for(int l=dep;l>=2;--l)//l>2
	{
	
		node=abs(Visited_Path[node].last_read);
		if(node==0)
		{return 1;}
		
		if(abs(stacked_nodes[node])>2)
		{return 1;}// only backtrack to here so return 1.
		
			
		if(stacked_nodes[node]>0&&reads_overlap_info->left_overlaps[abs(node)].size()>1)
		{return 0;}
		if(stacked_nodes[node]<0&&reads_overlap_info->right_overlaps[abs(node)].size()>1)
		{return 0;}
	}
	return 1;
}

void BreakLinks_read( reads_overlap_info *reads_overlap_info, map<int, int > &stacked_nodes,int node1, int node2)
{
	node1=abs(node1);
	node2=abs(node2);
	if(stacked_nodes[node1]>0)
	{
		int temp=abs(node2);
		if(stacked_nodes[node2]<0)
		{
			temp=-temp;
		}
		if(reads_overlap_info->right_overlaps[abs(node1)].count(temp))
		{
			reads_overlap_info->right_overlaps[abs(node1)].erase(temp);
			if(stacked_nodes[node1]>1)
			{
				stacked_nodes[abs(node1)]--;
			}
		}

	}


	if(stacked_nodes[node1]<0)
	{
		int temp=abs(node2);
		if(stacked_nodes[node2]>0)
		{
			temp=-temp;
		}

		if(reads_overlap_info->left_overlaps[abs(node1)].count(temp))
		{
			reads_overlap_info->left_overlaps[abs(node1)].erase(temp);
			if(stacked_nodes[node1]<-1)
			{
				stacked_nodes[abs(node1)]++;		
			}
		}


	}

	
	if(stacked_nodes[node2]<0)
	{
		int temp=abs(node1);
		if(stacked_nodes[node1]>0)
		{
			temp=-temp;
		}
		if(reads_overlap_info->right_overlaps[abs(node2)].count(temp))
		{
			reads_overlap_info->right_overlaps[abs(node2)].erase(temp);
				
		}
	}

	if(stacked_nodes[node2]>0)
	{
		int temp=abs(node1);
		if(stacked_nodes[node1]<0)
		{
			temp=-temp;
		}
		if(reads_overlap_info->left_overlaps[abs(node2)].count(temp))
		{
			reads_overlap_info->left_overlaps[abs(node2)].erase(temp);
				
		}
	}

}

void BacktrackBubbleRemoval_read(reads_overlap_info *reads_overlap_info,int last_read,int beg_read,map<int,struct BFS_reads_info > & Visited_Path , map<int ,int > &stacked_nodes)
{
	beg_read=abs(beg_read);
	last_read=abs(last_read);	
	int current_read=last_read;
	int dep=Visited_Path[last_read].depth;
	for(int l=dep;l>1;--l)
	{
		int previous_read=current_read;

		current_read=Visited_Path[current_read].last_read;
		if(current_read==0)
		{return;}
		
		if(stacked_nodes[current_read]>=1||stacked_nodes[current_read]<=-1)
		{
			

			if(abs(stacked_nodes[current_read])>2)
			{
				BreakLinks_read(reads_overlap_info,stacked_nodes,current_read,previous_read);
				break;
			}
			//else
			
			BreakLinks_read(reads_overlap_info,stacked_nodes,current_read,previous_read);
			
			if(Visited_Path[current_read].last_read==NULL)
			{
				break;
			}

			int free_read=current_read;

			if(beg_read==free_read)
			{
				break;
			}

			if(free_read!=last_read)
			{

				stacked_nodes[free_read]=stacked_nodes[free_read]/abs(stacked_nodes[free_read]);
				//freebkt->kmer_info.removed=1;			
			}
		}
	}
}


void BFSearchBubbleRemoval_read(reads_overlap_info *reads_overlap_info,reads_table *reads_table,int beg_read,int max_depth,int max_dist)
{
	map<int,struct BFS_reads_info > Visited_Path;
	map<int, int > stacked_nodes;
	int max_stack=300;
	int DepthTh=max_depth;
	int LenTh=30;
	bool RIGHT=0;
	map<int, list<int> > dist_reads;//neighborset
	dist_reads[0].push_back(beg_read);
	int NBs=1;
	int dist_searched=0;

	int new_node=abs(beg_read);
	if(beg_read>0)
	{
		stacked_nodes[abs(beg_read)]=1;
	}
	else
	{
		stacked_nodes[abs(beg_read)]=-1;
	}
	
	Visited_Path[new_node].cov=0;
	Visited_Path[new_node].depth=1;
	Visited_Path[new_node].last_read=0;
	
	map<int , list<int> >::iterator NB_it=dist_reads.begin();

	while(1)
	{
		NB_it=dist_reads.begin();
		if(NB_it==dist_reads.end())
		{break;}
		if(NB_it->second.size()==0)
		{dist_reads.erase(NB_it->first);continue;}

		if(NBs>max_stack)
		{
			break;
		}
		new_node=NB_it->second.front();

		NB_it->second.pop_front();
		NBs--;
		if(NB_it->second.size()==0)
		{
			dist_reads.erase(NB_it->first);
		}
		if(new_node>0)
		{
			RIGHT=1;
		}
		else
		{
			RIGHT=0;
		}

		new_node=abs(new_node);
		//if(new_node==173)
		//cout<<new_node<<endl;
		if(Visited_Path[new_node].depth>DepthTh||Visited_Path[new_node].len>LenTh)
		{continue;}
		//if(new_node==239)
		//{cout<<"";}

		if(RIGHT)
		{
			int rb=reads_overlap_info->right_overlaps[new_node].size();			
			if(stacked_nodes[new_node]==1&&rb>0)
			{
				stacked_nodes[new_node]=1+rb;
			}
			if(rb==0)
			{
				stacked_nodes[new_node]=2;
				
				if(abs(new_node)==abs(beg_read))
				{					
					continue;
				}
				//tip end reached so backtrack to the branching position.
				if(!isSimplePath_read(reads_overlap_info,new_node,Visited_Path, stacked_nodes))
				{
					continue;
				}

				BacktrackBubbleRemoval_read(reads_overlap_info,new_node,beg_read,Visited_Path,stacked_nodes);
				stacked_nodes[new_node]=1;
				continue;

			}

			map<int32_t, vector<int32_t> >::iterator tmp_it,tmp_it_n;
			for(tmp_it=reads_overlap_info->right_overlaps[new_node].begin();tmp_it!=reads_overlap_info->right_overlaps[new_node].end();)
			{
				tmp_it_n=tmp_it;
				tmp_it_n++;
				int next_read=tmp_it->first;
				// not in stack
				if(stacked_nodes[abs(next_read)]==0)
				{
//					Visited_Path[*ptr].cov=(int)(Visited_Path[new_node].cov+edge_ptr->edge_cov);
					Visited_Path[abs(next_read)].cov=(int)(Visited_Path[new_node].cov+reads_overlap_info->cov_vt[abs(next_read)]);
					Visited_Path[abs(next_read)].depth=(Visited_Path[abs(new_node)].depth+1);
					//int cum_len=(int)(Visited_Path[abs(next_read)].len+1);
					int cum_len=Visited_Path[abs(next_read)].depth;
					//Visited_Path[abs(next_read)].len=cum_len;//(Visited_Path[new_node].len+edge_ptr->len+1);
					Visited_Path[abs(next_read)].last_read=new_node;
					
					if(next_read<0)
					{
						stacked_nodes[abs(next_read)]=-1;	
					}
					else
					{
						stacked_nodes[abs(next_read)]=1;
					}
					dist_reads[cum_len].push_back(next_read);
					NBs++;

				}
				else
				{
					if((stacked_nodes[abs(next_read)]>0&&next_read>0)||(stacked_nodes[abs(next_read)]<0&&next_read<0))
					{
						
						if((Visited_Path[abs(new_node)].cov+reads_overlap_info->cov_vt[abs(next_read)]<=Visited_Path[abs(next_read)].cov))//||(BackCheckLoop(*ptr,new_node,Visited_Path)==1)//loop
						{
							//backtrack if the same direction is found
							
							if(!isSimplePath_read(reads_overlap_info,new_node,Visited_Path, stacked_nodes))
							{
								tmp_it=tmp_it_n;
								continue;
							}
							
							//backtrack the current path, common operation in this search
							if(stacked_nodes[new_node]>2)
							{
								BreakLinks_read(reads_overlap_info,stacked_nodes, new_node, next_read);
								tmp_it=tmp_it_n;
								continue;
							}
							else
							{
								if(stacked_nodes[new_node]<-2)
								{
									BreakLinks_read(reads_overlap_info,stacked_nodes, new_node, next_read);
									
									tmp_it=tmp_it_n;
									continue;
								}
								else
								{
									int free_node=(new_node);

									if(free_node==beg_read)
									{
										BreakLinks_read(reads_overlap_info,stacked_nodes, beg_read, new_node);
									}
									
									//free the node and edge.
									BreakLinks_read(reads_overlap_info,stacked_nodes,new_node,next_read);
									BacktrackBubbleRemoval_read(reads_overlap_info,new_node,beg_read,Visited_Path,stacked_nodes);
									tmp_it=tmp_it_n;
									continue;
			//
								}
							}

						}
						else
						{
							//backtrack the original path, rare operation in this search, can lead to errors
							
							//BacktrackBubbleRemoval(ht,merge_ht,*ptr,beg_bkt,beg_bkt,Visited_Path , stacked_nodes,K_size);
							BacktrackBubbleRemoval_read(reads_overlap_info,next_read,beg_read,Visited_Path,stacked_nodes);
							
							//marginal 
							Visited_Path[abs(next_read)].cov=(int)(Visited_Path[abs(new_node)].cov+reads_overlap_info->cov_vt[abs(next_read)]);
							Visited_Path[abs(next_read)].depth=Visited_Path[abs(new_node)].depth+1;
							Visited_Path[abs(next_read)].len=Visited_Path[abs(new_node)].depth;
							//marginal 

							Visited_Path[abs(next_read)].last_read=abs(new_node);
							tmp_it=tmp_it_n;
							continue;
						}

					}
					else
					{

						//don't do anything,since both strands are visited.

					}

				}
				tmp_it=tmp_it_n;
			}
			
		}
		else
		{

			int lb=reads_overlap_info->left_overlaps[new_node].size();			
			if(stacked_nodes[new_node]==-1&&lb>0)
			{
				stacked_nodes[new_node]=-1-lb;
			}
			if(lb==0)
			{
				stacked_nodes[new_node]=-2;

				
				if(abs(new_node)==abs(beg_read))
				{					
					continue;
				}
				//tip end reached so backtrack to the branching position.
				if(!isSimplePath_read(reads_overlap_info,new_node,Visited_Path, stacked_nodes))
				{
					continue;
				}

				BacktrackBubbleRemoval_read(reads_overlap_info,new_node,beg_read,Visited_Path,stacked_nodes);
				stacked_nodes[new_node]=-1;
				continue;

			}

			map<int32_t, vector<int32_t> >::iterator tmp_it,tmp_it_n;
			for(tmp_it=reads_overlap_info->left_overlaps[new_node].begin();tmp_it!=reads_overlap_info->left_overlaps[new_node].end();)
			{
				tmp_it_n=tmp_it;
				tmp_it_n++;
				int next_read=tmp_it->first;
				// not in stack
				if(stacked_nodes[abs(next_read)]==0)
				{
//					Visited_Path[*ptr].cov=(int)(Visited_Path[new_node].cov+edge_ptr->edge_cov);
					Visited_Path[abs(next_read)].cov=(int)(Visited_Path[new_node].cov+reads_overlap_info->cov_vt[abs(next_read)]);
					Visited_Path[abs(next_read)].depth=(Visited_Path[abs(new_node)].depth+1);
					//int cum_len=(int)(Visited_Path[abs(next_read)].len+1);
					int cum_len=Visited_Path[abs(next_read)].depth;
					//Visited_Path[abs(next_read)].len=cum_len;//(Visited_Path[new_node].len+edge_ptr->len+1);
					Visited_Path[abs(next_read)].last_read=new_node;
					
					if(next_read<0)
					{
						stacked_nodes[abs(next_read)]=1;	
					}
					else
					{
						stacked_nodes[abs(next_read)]=-1;
					}
					dist_reads[cum_len].push_back(-next_read);
					NBs++;

				}
				else
				{
					if((stacked_nodes[abs(next_read)]<0&&next_read>0)||(stacked_nodes[abs(next_read)]>0&&next_read<0))
					{
						
						if((Visited_Path[abs(new_node)].cov+reads_overlap_info->cov_vt[abs(next_read)]<=Visited_Path[abs(next_read)].cov))//||(BackCheckLoop(*ptr,new_node,Visited_Path)==1)//loop
						{

							//backtrack if the same direction is found

							if(!isSimplePath_read(reads_overlap_info,new_node,Visited_Path, stacked_nodes))
							{
								tmp_it=tmp_it_n;
								continue;
							}
							
							//backtrack the current path, common operation in this search
							if(stacked_nodes[new_node]>2)
							{
								BreakLinks_read(reads_overlap_info,stacked_nodes, new_node, next_read);
								tmp_it=tmp_it_n;
								continue;
							}
							else
							{
								if(stacked_nodes[new_node]<-2)
								{
									BreakLinks_read(reads_overlap_info,stacked_nodes, new_node, next_read);
									tmp_it=tmp_it_n;
									continue;
								}
								else
								{
									int free_node=(new_node);

									if(free_node==beg_read)
									{
										BreakLinks_read(reads_overlap_info,stacked_nodes, beg_read, new_node);
									}
									
									//free the node and edge.
									BreakLinks_read(reads_overlap_info,stacked_nodes,new_node,next_read);
									BacktrackBubbleRemoval_read(reads_overlap_info,new_node,beg_read,Visited_Path,stacked_nodes);
									tmp_it=tmp_it_n;
									continue;
								}
							}
						}
						else
						{
							//backtrack the original path, rare operation in this search, can lead to errors
							BacktrackBubbleRemoval_read(reads_overlap_info,next_read,beg_read,Visited_Path,stacked_nodes);
							
							//marginal 
							Visited_Path[abs(next_read)].cov=(int)(Visited_Path[abs(new_node)].cov+reads_overlap_info->cov_vt[abs(next_read)]);
							Visited_Path[abs(next_read)].depth=Visited_Path[abs(new_node)].depth+1;
							Visited_Path[abs(next_read)].len=Visited_Path[abs(new_node)].depth;
							//marginal 

							Visited_Path[abs(next_read)].last_read=abs(new_node);
							tmp_it=tmp_it_n;
							continue;
						}

					}
					else
					{

						//don't do anything,since both strands are visited.

					}


				}

				
				tmp_it=tmp_it_n;




			}

		}

	}

}



void ConstructReadsOverlaps(string reads_info_name,reads_table* reads_table)
{
	// 
	ofstream out_extended_reads("Extended_reads.txt");
	reads_overlap_info reads_overlap_info;
	bool EXACT=1;
	ifstream in_reads_info;
	in_reads_info.open(reads_info_name.c_str());

	string tag;
	int32_t num_reads=0,n_selected_reads=0;
	int n_ctgs;
	vector<vector<int32_t> > read_contig_index;
	map<int32_t,vector<int32_t> > contig_in_reads;
	map<vector<int>,int> selected_reads;

	reads_overlap_info.cov_vt.push_back(0);

	while(in_reads_info>>tag)
	{
		if(tag.size()<=1)
		{
			break;
		}
		num_reads++;
		in_reads_info>>n_ctgs;
		vector<int> contigs_vt,contigs_vt_rc;
		int ctg;
		// get the contigs into a vector
		for (int i=0;i<n_ctgs;++i)
		{
			in_reads_info>>ctg;
			contigs_vt.push_back(ctg);
		}
		int cod;
		for (int i=0;i<n_ctgs;++i)
		{
			in_reads_info>>cod;
			
		}
		int offset;
		for (int i=0;i<n_ctgs-1;++i)
		{
			in_reads_info>>offset;
			
		}
		//take rc
		contigs_vt_rc=contigs_vt;
		reverse(contigs_vt_rc.begin(),contigs_vt_rc.end());
		
		// compare and take the smaller vector
		for (int i=0;i<n_ctgs;++i)
		{
			contigs_vt_rc[i]=-contigs_vt_rc[i];
		}
		bool take_rc=0;	
		for (int i=0;i<n_ctgs;++i)
		{
			if(contigs_vt_rc[i]>contigs_vt[i])
			{
				take_rc=0;
				break;
			}
			if(contigs_vt_rc[i]<contigs_vt[i])
			{
				take_rc=1;
				break;
			}
		}
		if(take_rc)
		{
			contigs_vt=contigs_vt_rc;
		}
		if(selected_reads[contigs_vt]!=0)
		{
			
			reads_overlap_info.cov_vt[abs(selected_reads[contigs_vt])]++;
			continue;
		}
		else
		{
			n_selected_reads++;
			if(take_rc)
			{
				selected_reads[contigs_vt]=-n_selected_reads;
				reads_overlap_info.cov_vt.push_back(1);
			}
			else
			{
				selected_reads[contigs_vt]=n_selected_reads;
				reads_overlap_info.cov_vt.push_back(1);
			}
		}


		// now in selected reads
		for (int i=0;i<n_ctgs;++i)
		{
			ctg=contigs_vt[i];	
			if(contig_in_reads[abs(ctg)].size()==0||contig_in_reads[abs(ctg)].back()!=n_selected_reads)//check and don't save duplicate values
			{
				contig_in_reads[abs(ctg)].push_back(n_selected_reads);
			}
		}

		read_contig_index.push_back(contigs_vt);
	}
	in_reads_info.close();
	
	reads_overlap_info.left_overlaps.resize(n_selected_reads+1);
	reads_overlap_info.right_overlaps.resize(n_selected_reads+1);
	reads_overlap_info.contained_vt.resize(n_selected_reads+1);
	reads_overlap_info.used_vt.resize(n_selected_reads+1);

	for (int i=0;i<read_contig_index.size();++i)
	{

		if(EXACT)
		{

			
			for (int round1=1;round1<=2;++round1)
			{
				//front/rear overlap
				vector<int32_t> current_read=read_contig_index[i];
				
				if(round1==2)
				{
					reverse(current_read.begin(),current_read.end());
					for (int j=0;j<current_read.size();++j)
					{
						current_read[j]=-current_read[j];
					}

				}

			
				for (int r=0;r<contig_in_reads[abs(current_read[0])].size();++r)
				{
					int read_idx=contig_in_reads[abs(current_read[0])][r];
					if(read_idx-1==i)
					{
						continue;
					}

					vector<int32_t> ctgs_vt=read_contig_index[read_idx-1];
					for (int round2=1;round2<=2;++round2)
					{
						if(round2==2)
						{
							reverse(ctgs_vt.begin(),ctgs_vt.end());
							for (int j=0;j<ctgs_vt.size();++j)
							{
								ctgs_vt[j]=-ctgs_vt[j];
							}
						}
						for (int j=0;j<ctgs_vt.size();++j)
						{
							int pos=0;
							bool overlap=1;
							for(int k=j;k!=ctgs_vt.size();++k)
							{
								if (current_read[pos]!=ctgs_vt[k])
								{
									overlap=0;
									break;
								}
								pos++;
								if(pos==current_read.size()&&k!=ctgs_vt.size()-1)
								{
									overlap=0;
									break;
								}
							}
							//check if we have a proper overlap.
							if(overlap==1)
							{
								vector<int> edge;
								if(ctgs_vt.size()-j==current_read.size())
								{
									reads_overlap_info.contained_vt[i+1]=1;
									cout<<"";//contained overlap
								}


								for(int k=0;k<j;++k)
								{
									edge.push_back(ctgs_vt[k]);
								}

								// record the overlap and break;
								if(round1==1)
								{
									//front overlap
									if(round2==1)
									{
										//
										reads_overlap_info.left_overlaps[i+1][read_idx]=edge;
									}
									else
									{
										//flip
										reverse(edge.begin(),edge.end());
										for(int k=0;k<edge.size();++k)
										{
											edge[k]=-edge[k];
										}
										reads_overlap_info.left_overlaps[i+1][-read_idx]=edge;
									}

								}
								else
								{
									//rear overlap
									if(round2==1)
									{
										//flip
										reverse(edge.begin(),edge.end());
										for(int k=0;k<edge.size();++k)
										{
											edge[k]=-edge[k];
										}
										reads_overlap_info.right_overlaps[i+1][-read_idx]=edge;
									}
									else
									{
										
										reads_overlap_info.right_overlaps[i+1][read_idx]=edge;
									}
								}

								break;

							}
					
						}
					
					}

				}

			}

		}

	}
	
	for(int i=0;i<reads_overlap_info.right_overlaps.size();++i)
	{
		map<int32_t, vector<int32_t> >::iterator temp_it1,temp_it2;
		for(temp_it1=reads_overlap_info.right_overlaps[i].begin();temp_it1!=reads_overlap_info.right_overlaps[i].end(); )
		{
			temp_it2=temp_it1;
			temp_it2++;
			if(reads_overlap_info.contained_vt[abs(temp_it1->first)])
			{
				reads_overlap_info.right_overlaps[i].erase(temp_it1);
			}
			temp_it1=temp_it2;
		}

		for(temp_it1=reads_overlap_info.left_overlaps[i].begin();temp_it1!=reads_overlap_info.left_overlaps[i].end(); )
		{
			temp_it2=temp_it1;
			temp_it2++;
			if(reads_overlap_info.contained_vt[abs(temp_it1->first)])
			{
				reads_overlap_info.left_overlaps[i].erase(temp_it1);
			}
			temp_it1=temp_it2;
		}

	}

	int max_dist=20;
	int max_depth=20;
	
	for(int i=1;i<reads_overlap_info.contained_vt.size();++i)
	{
		cout<<i<<endl;
		int beg_read=i;
		if(reads_overlap_info.contained_vt[i]==0)
		{
			if(reads_overlap_info.right_overlaps[i].size()>1)
			{
				BFSearchBubbleRemoval_read(&reads_overlap_info,reads_table, beg_read, max_depth, max_dist);
			}
			cout<<"l"<<endl;
			if(reads_overlap_info.left_overlaps[i].size()>1)
			{
				BFSearchBubbleRemoval_read(&reads_overlap_info,reads_table, -beg_read, max_depth, max_dist);
			}
		}
	}

	cout<<"Overlap graph constructed"<<endl;

	for(int read_idx=1;read_idx<reads_overlap_info.contained_vt.size();++read_idx)
	{
		vector<int> left_ext,right_ext,extended_read;
		if(reads_overlap_info.contained_vt[read_idx]==0&&reads_overlap_info.used_vt[read_idx]==0)
		{
			
			for (int it=1;it<=2;++it)
			{
				int current_read=read_idx;
				reads_overlap_info.used_vt[current_read]=1;
	
				bool Right;
				if(it==1)
				{
					Right=1;
				}
				else
				{
					Right=0;
				}

				while(1)
				{
					int next_read;
					if(Right==1)
					{		
						if(reads_overlap_info.right_overlaps[abs(current_read)].size()==1)
						{
							next_read=reads_overlap_info.right_overlaps[abs(current_read)].begin()->first;
							if(reads_overlap_info.used_vt[abs(next_read)]==1)
							{
								break;
							}
							vector<int> edge=reads_overlap_info.right_overlaps[abs(current_read)].begin()->second;
							for(int e=0;e<edge.size();++e)
							{
								if(it==1)
								{
									right_ext.push_back(edge[e]);
								}
								else
								{
									left_ext.push_back(-edge[e]);
								}
							}
						}
						else
						{break;}
					}
					else
					{
						if(Right==0)
						{
							if(reads_overlap_info.left_overlaps[abs(current_read)].size()==1)
							{
								next_read=reads_overlap_info.left_overlaps[abs(current_read)].begin()->first;
								if(reads_overlap_info.used_vt[abs(next_read)]==1)
								{
									break;
								}
								vector<int> edge=reads_overlap_info.left_overlaps[abs(current_read)].begin()->second;
									
								for(int e=0;e<edge.size();++e)
								{
									if(it==1)
									{
										right_ext.push_back(-edge[e]);
									}
									else
									{
										left_ext.push_back(edge[e]);
									}
								}
							}
							else
							{
								break;
							}
						}
					}
					reads_overlap_info.used_vt[abs(next_read)]=1;
					if(next_read<0)
					{
						Right=!Right;
					}
					current_read=abs(next_read);
					if(reads_overlap_info.used_vt[abs(current_read)])
					{break;}
				}
			}

			for(int jj=0;jj<left_ext.size();++jj)
			{
				extended_read.push_back(left_ext[jj]);
				out_extended_reads<<left_ext[jj]<<" ";
			}
			for(int jj=0;jj<read_contig_index[read_idx-1].size();++jj)
			{
				extended_read.push_back(read_contig_index[read_idx-1][jj]);
				out_extended_reads<<read_contig_index[read_idx-1][jj]<<" ";
			}
			for(int jj=0;jj<right_ext.size();++jj)
			{
				extended_read.push_back(right_ext[jj]);
				out_extended_reads<<right_ext[jj]<<" ";
			}
			out_extended_reads<<endl;

		}

	}
	cout<<"";
	/////////////finished

}






/*
void CollectingNonContainedPairs(struct hashtable *ht1,struct hashtable *merge_ht1, struct hashtable2 *ht2, struct hashtable2 *merge_ht2,int K_size,vector<string>& filenames_vt, contigs_info * contigs_info,string ContigFilename)
{

	time_t beg_time,read_time;
	string in_fname=ContigFilename;
	ifstream Contigs_in(in_fname.c_str());
	size_t num_Contigs=0;
	string tag,s,kmer,str,seq_s;
	uint64_t f_seq,hv;
	struct kmer_t2 f_seq_t2;
	size_t hash_idx;
	bool found;
	bool flip_1,flip_2,flip_0;
	size_t ht_sz;
	size_t numReads=0;
	ofstream out_reads("reads_temp.txt");
	
	int boundary=0,removed=0,bridge=0;
	if(K_size<=32)
	{
		ht_sz=ht1->ht_sz;
	}
	
	cout<<"Contigs remapping."<<endl;
	if(ContigFilename=="Contigs.txt")
	{
		ContigsRemapping(ht1,ht2, K_size, contigs_info,ContigFilename,0);
		
	}

	AppendMergeHT(ht1, merge_ht1);

	//cout<<"Collecting informative reads."<<endl;
	time(&beg_time);

	int64_t dist_sum=0,p_cnt=0;
	int mean_dist=0;
	int lib_no=0;
	
	for(size_t ii=0;ii<filenames_vt.size();ii+=2)
	{
		lib_no++;
	
		dist_sum=0,p_cnt=0;
		cout<<"Processing library: "<<ii<<endl;
		ifstream in_reads;
		in_reads.open(filenames_vt[ii].c_str());
		
		struct read_t Read1,Read2;

		Read1.read_bits =(uint64_t*) malloc(1000000/4+100);
		Read2.read_bits =(uint64_t*) malloc(1000000/4+100);

		int nLines1=0,nLines2=0;


		bool fq_flag=0;


		getline(in_reads,str);
		if(fq_flag==0&&str[0]=='@')
		{
			fq_flag=1;	
		}
		in_reads.close();

		in_reads.clear();
		in_reads.open(filenames_vt[ii].c_str());

		bool read_success=0;

		read_success=1;

		string tag,qs,n_tag;
		string QS_s;

		while(read_success)
		{
			if(fq_flag)
			{
				read_success=get_a_fastq_read(in_reads,tag,seq_s,QS_s);
					
			}
			else
			{
				read_success=get_a_fasta_read(in_reads,tag,seq_s,n_tag);
			
			}	
			if(read_success==0)
			{break;}



				
			int seq_sz=seq_s.size();
			
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
			
			if(bad_flag)
			{continue;}
							

			bad_flag=0;
					
			numReads++;


						
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
					


			Init_Read(seq_s,Read1);
			string read_str=seq_s;
			seq_s.clear();
			

			uint64_t bits1;
			int contig_no=-1;
			bool output_current_read=0;

			for(int i=0;i<Read1.readLen-K_size+1;++i )
			{
				
				get_sub_arr(Read1.read_bits,Read1.readLen,i,K_size,&bits1);
				f_seq=get_rev_comp_seq(bits1,K_size);
				flip_0=0;
				if(bits1>f_seq)
				{
					bits1=f_seq;
					flip_0=1;
				}

				hv=MurmurHash64A(&bits1,sizeof(bits1),0);
				hash_idx=(size_t) (hv%ht_sz);
				struct bucket **ptr1;
				ptr1= &(ht1->store_pos[hash_idx]);
				found=look_up_in_a_list(bits1,&ptr1);
				if(found)
				{
	
					if((*ptr1)->kmer_info.removed==0)
					{
							

						int cont1=(*ptr1)->kmer_info.contig_no;
						if(contig_no<0)
						{
							contig_no=cont1;
						}
						else
						{
							if(cont1!=contig_no)
							{
								output_current_read=1;
								bridge++;
							}
							if((*ptr1)->kmer_info.masked==1)
							{
								output_current_read=1;
								boundary++;
							}
							break;	
						}
								
					}
					else
					{
						output_current_read=1;
						removed++;
						break;


					}

						
				}
			}
				


			



			if(output_current_read)
			{
				if(tag[0]=='@')
				{tag[0]='>';}
				out_reads<<tag<<endl;
				out_reads<<read_str<<endl;
			}


		}



	}

	//cout<<"boundary"<< boundary<<endl;
	//cout<<"removed"<<removed<<endl;
	//cout<<"bridge"<<bridge<<endl;

}
*/


//rewrite...	
/*
int BFSearchGapCloser(struct hashtable* ht,struct hashtable* merge_ht, struct bucket* bktptr,struct bucket* obj_bktptr,int K_size, stacked_bucket &kmer_stack_beg,int max_depth,int max_dist,SuperRead_t *SuperRead, search_info *search_info)
{
	map<bucket*,struct BFS_path_info > Visited_Path;
	map<struct bucket* ,int > stacked_nodes;
	struct bucket *beg_bkt= bktptr;
	int max_stack=500;
	int DepthTh=max_depth;//min(300/gap,20);
	int LenTh=max_dist;
	bool RIGHT=0;
	struct stacked_bucket stacked_bkt=kmer_stack_beg;

	map<int , list<stacked_bucket> > dist_ctgs;//neighborset
	dist_ctgs[0].push_back(kmer_stack_beg);
	int NBs=1;

	int dist_searched=0;
	bucket* new_node=stacked_bkt.bktptr;
	uint64_t kmer,f_kmer;
	if(stacked_bkt.RightSearch)
	{
		stacked_nodes[new_node]=1;
	}
	else
	{
		stacked_nodes[new_node]=-1;
	}
	//search direction

	//Visited_Path[new_node].cov=0;
	Visited_Path[new_node].depth=1;
	Visited_Path[new_node].len=K_size;
	Visited_Path[new_node].last_bkt=NULL;
	Visited_Path[new_node].BothEndsUsed=0;

	//Visited_Path[new_node].last_bkt_edge=NULL;

	map<int , list<stacked_bucket> >::iterator NB_it=dist_ctgs.begin();
	while(1)
	{
		NB_it=dist_ctgs.begin();
		
		if(NB_it==dist_ctgs.end())
		{break;}
		if(NB_it->first>max_dist)
		{
			if(SuperRead->PathsCnt!=1)
			{return -10000;}
			else
			{break;}
		}
		if(NB_it->second.size()==0)
		{dist_ctgs.erase(NB_it->first);continue;}
		if(NBs>max_stack)
		{
			
			SuperRead->PathsCnt=0;
			return -10000;
			break;
		}
		stacked_bkt=NB_it->second.front();

		NB_it->second.pop_front();
		NBs--;
		if(NB_it->second.size()==0)
		{
			dist_ctgs.erase(NB_it->first);
			
		}		
		new_node=stacked_bkt.bktptr;

		RIGHT=stacked_bkt.RightSearch;
		if(Visited_Path[new_node].depth>DepthTh||Visited_Path[new_node].len>LenTh)
		{continue;}

		if(RIGHT)
		{
			bool r_found=0;
			struct edge_node *edge_ptr;
			edge_ptr=(new_node)->kmer_info.right;
			if(edge_ptr==NULL)
			{
				continue;
			}

			while(edge_ptr!=NULL)
			{
				kmer=(new_node)->kmer_t.kmer;
				f_kmer=kmer;
				f_kmer=get_rev_comp_seq(kmer,K_size);			
				int edge_len=(int)(edge_ptr->len);
				for(int g=edge_len;g>=0;--g)
				{
					uint64_t b=(edge_ptr->edge>>2*g)&0x3;
					kmer&=(~((uint64_t)0x3<<(2*(K_size-1))));
					kmer<<=2;
					kmer|=b;
				}
				bool r_flip=0;
				uint64_t kmer2,f_kmer2;
				kmer2=kmer;
				f_kmer2=get_rev_comp_seq(kmer2,K_size);

				if(kmer2>f_kmer2)
				{
					kmer2=f_kmer2;
					r_flip=1;
				}
				else
				{
					r_flip=0;
				}

				uint64_t hv=MurmurHash64A(&kmer2,sizeof(kmer2),0);
				uint64_t hash_idx=(size_t) (hv%ht->ht_sz);

				struct bucket** ptr;
				ptr= &(ht->store_pos[hash_idx]);
				r_found=look_up_in_a_list(kmer2,&ptr);

				if(r_found&&(*ptr)->kmer_info.removed==1)
				{
					uint64_t bits1=(*ptr)->kmer_t.kmer;
					
					uint64_t hv=MurmurHash64A(&bits1,sizeof(bits1),0);
					uint64_t hash_idx=(size_t) (hv%(merge_ht->ht_sz));

					struct bucket_rm ** ptr_rm;
					ptr_rm=(bucket_rm **) &(merge_ht->store_pos[hash_idx]);
					r_found=look_up_in_a_list_rm(bits1,&ptr_rm);
					if(r_found==1)
					{
						r_flip^=(*ptr_rm)->flip;
						bits1=(*ptr_rm)->merged_kmer.kmer;
					}
					hv=MurmurHash64A(&bits1,sizeof(bits1),0);
					hash_idx=(size_t) (hv%(ht->ht_sz));
					ptr= &(ht->store_pos[hash_idx]);
					r_found=look_up_in_a_list(bits1,&ptr);

				}



				if(r_found)
				{
					if((*ptr)==(obj_bktptr)&&(search_info->Flip_End==r_flip))
					{
						SuperRead->PathsCnt++;
						SuperRead->PathFound=1;

						int edge_len=edge_ptr->len;
						SuperRead->PathsLength=(int)(Visited_Path[new_node].len+edge_len+1);
						//return (int)(Visited_Path[new_node].len+edge_len+1);//found distance
					}

					// not in stack
					if(stacked_nodes[*ptr]==0)
					{
						//Visited_Path[*ptr].cov=(int)(Visited_Path[new_node].cov+edge_ptr->edge_cov);
						Visited_Path[*ptr].depth=(Visited_Path[new_node].depth+1);
						
						int edge_len=edge_ptr->len;
						int cum_len=(int)(Visited_Path[new_node].len+edge_len+1);
						Visited_Path[*ptr].len=cum_len;//(Visited_Path[new_node].len+edge_ptr->len+1);

						Visited_Path[*ptr].last_bkt=new_node;
						Visited_Path[*ptr].last_bkt_edge=edge_ptr;
						Visited_Path[*ptr].BothEndsUsed=0;		

						if(r_flip)
						{
							stacked_nodes[*ptr]=-1;
							stacked_bkt.bktptr=*ptr;
							stacked_bkt.RightSearch=0;
							//stacked_bkt.BothEndsUsed=0;
						
						}
						else
						{
							stacked_nodes[*ptr]=1;
							stacked_bkt.bktptr=*ptr;
							stacked_bkt.RightSearch=1;
							//stacked_bkt.BothEndsUsed=0;							
						}
						dist_ctgs[cum_len].push_back(stacked_bkt);
						NBs++;

					}
					else
					{
						if((stacked_nodes[*ptr]==1&&r_flip==0)||(stacked_nodes[*ptr]==-1&&r_flip==1)||stacked_nodes[*ptr]==3)
						{
							edge_ptr=edge_ptr->nxt_edge;
							continue;
						
						}
						else
						{
							int edge_len=edge_ptr->len;
							int cum_len=(int)(Visited_Path[new_node].len+edge_len+1);
							if(r_flip)
							{
								stacked_nodes[*ptr]=3;
								stacked_bkt.bktptr=*ptr;
								stacked_bkt.RightSearch=0;
								//stacked_bkt.BothEndsUsed=1;
								//stacked_bkt.BothSideSearch=1;
								Visited_Path[*ptr].BothEndsUsed=1;
								
							}
							else
							{
								stacked_nodes[*ptr]=3;
								stacked_bkt.bktptr=*ptr;
								stacked_bkt.RightSearch=1;
								//stacked_bkt.BothEndsUsed=1;
								//stacked_bkt.BothSideSearch=1;
								Visited_Path[*ptr].BothEndsUsed=1;
							}
					
							dist_ctgs[cum_len].push_back(stacked_bkt);

						}


					}

				}

				edge_ptr=edge_ptr->nxt_edge;
			}


		}
		else
		{

			bool l_found=0;
			struct edge_node *edge_ptr;
			edge_ptr=(new_node)->kmer_info.left;


			int lb=0;
			while(edge_ptr!=NULL)
			{
				lb++;
				edge_ptr=edge_ptr->nxt_edge;
			}
			edge_ptr=(new_node)->kmer_info.left;

			
			if(lb==0)
			{
				continue;
			}

			while(edge_ptr!=NULL)
			{
				kmer=(new_node)->kmer_t.kmer;
				f_kmer=kmer;
				f_kmer=get_rev_comp_seq(kmer,K_size);
				int edge_len=(int)(edge_ptr->len);

				for(int g=0;g<=edge_len;++g)
				{
					uint64_t b=(edge_ptr->edge>>2*g)&0x3;

					kmer>>=2;
					//R_shift_NB(kmer.kmer,2,2);
					kmer|=(b<<(2*(K_size-1)));
				}

				bool l_flip=0;
				uint64_t kmer2,f_kmer2;;

				kmer2=kmer;
				//f_kmer2=kmer;

				f_kmer2=get_rev_comp_seq(kmer2,K_size);


				if(kmer2>f_kmer2)
				{
					kmer2=f_kmer2;
					l_flip=1;
				}
				else
				{l_flip=0;}

				uint64_t hv=MurmurHash64A(&kmer2,sizeof(kmer2),0);

				uint64_t hash_idx=(size_t) (hv%ht->ht_sz);

				struct bucket** ptr;

				ptr= &(ht->store_pos[hash_idx]);
				l_found=look_up_in_a_list(kmer2,&ptr);

				if(l_found&&(*ptr)->kmer_info.removed==1)
				{
					//cout<<"WarningL"<<endl;

					uint64_t bits1=(*ptr)->kmer_t.kmer;
					
					uint64_t hv=MurmurHash64A(&bits1,sizeof(bits1),0);
					uint64_t hash_idx=(size_t) (hv%(merge_ht->ht_sz));

					struct bucket_rm ** ptr_rm;
					ptr_rm=(bucket_rm **) &(merge_ht->store_pos[hash_idx]);
					l_found=look_up_in_a_list_rm(bits1,&ptr_rm);
					if(l_found==1)
					{
						l_flip^=(*ptr_rm)->flip;
						bits1=(*ptr_rm)->merged_kmer.kmer;
					}
					hv=MurmurHash64A(&bits1,sizeof(bits1),0);
					hash_idx=(size_t) (hv%(ht->ht_sz));
					ptr= &(ht->store_pos[hash_idx]);
					l_found=look_up_in_a_list(bits1,&ptr);

					
				}


				if(l_found)
				{
					if((*ptr)==obj_bktptr&&(search_info->Flip_End==!l_flip))
					{
						SuperRead->PathsCnt++;
						SuperRead->PathFound=1;
						int edge_len=edge_ptr->len;
						SuperRead->PathsLength=(int)(Visited_Path[new_node].len+edge_len+1);
						//return (int)(Visited_Path[new_node].len+edge_len+1);//found distance
					}

					if(stacked_nodes[*ptr]==0)
					{
						//Visited_Path[*ptr].cov=(int)(Visited_Path[new_node].cov+edge_ptr->edge_cov);
						Visited_Path[*ptr].depth=(Visited_Path[new_node].depth+1);
						int edge_len=edge_ptr->len;
						int cum_len=(int)(Visited_Path[new_node].len+edge_len+1);
						Visited_Path[*ptr].len=cum_len;

						Visited_Path[*ptr].last_bkt=new_node;
						Visited_Path[*ptr].last_bkt_edge=edge_ptr;
						Visited_Path[*ptr].BothEndsUsed=0;
						if(l_flip)
						{
							stacked_nodes[*ptr]=1;
							stacked_bkt.bktptr=*ptr;
							stacked_bkt.RightSearch=1;
							//stacked_bkt.BothEndsUsed=0;
							//kmer_stack.push_back(stacked_bkt);
							
						}
						else
						{
							stacked_nodes[*ptr]=-1;
							stacked_bkt.bktptr=*ptr;
							stacked_bkt.RightSearch=0;
							//stacked_bkt.BothEndsUsed=0;
							//kmer_stack.push_back(stacked_bkt);
							
						}
						dist_ctgs[cum_len].push_back(stacked_bkt);
						NBs++;

					}
					else
					{

						if((stacked_nodes[*ptr]==1&&l_flip==1)||(stacked_nodes[*ptr]==-1&&l_flip==0)||stacked_nodes[*ptr]==3)
						{
							edge_ptr=edge_ptr->nxt_edge;
							continue;
						
						}
						else
						{
							int edge_len=edge_ptr->len;
							int cum_len=(int)(Visited_Path[new_node].len+edge_len+1);
							
							if(l_flip)
							{
								stacked_nodes[*ptr]=3;
								stacked_bkt.bktptr=*ptr;
								stacked_bkt.RightSearch=1;
								stacked_bkt.BothEndsUsed=1;
								Visited_Path[*ptr].BothEndsUsed=1;
							}
							else
							{
								stacked_nodes[*ptr]=3;
								stacked_bkt.bktptr=*ptr;
								stacked_bkt.RightSearch=0;
								stacked_bkt.BothEndsUsed=1;
								Visited_Path[*ptr].BothEndsUsed=1;
							}
							dist_ctgs[cum_len].push_back(stacked_bkt);
							NBs++;
							

						}
					}

				}

			}

		}

	}


	if(SuperRead->PathsCnt!=1)
	{return -10000;}
	else
	{
		bucket * backtrack_bkt=obj_bktptr;
		int depth=Visited_Path[backtrack_bkt].depth+5;
		string extension;
		for (int i=depth;i>0;--i)
		{

			if(Visited_Path[backtrack_bkt].last_bkt==NULL)
			{
				SuperRead->extension=extension;
				return SuperRead->PathsLength;		
			}
			edge_node* edge=Visited_Path[backtrack_bkt].last_bkt_edge;
			uint64_t edge_bits=edge->edge;
			if(stacked_nodes[Visited_Path[backtrack_bkt].last_bkt]<0)
			{
				edge_bits=get_rev_comp_seq(edge_bits,(edge->len+1));
			}
			char edge_cstr[1000];
			bitsarr2str(&edge_bits,(edge->len+1),edge_cstr,1);
			extension=edge_cstr+extension;
			
			backtrack_bkt=Visited_Path[backtrack_bkt].last_bkt;

			if(Visited_Path[backtrack_bkt].BothEndsUsed|i==0)
			{
				extension.clear();
				SuperRead->extension.clear();
				return -10000;
			}
			


		}
		return SuperRead->PathsLength;
	}

}
*/

int BFSearchGapCloser_v2(struct hashtable* ht,struct hashtable* merge_ht, uint64_t beg_kmer,uint64_t end_kmer,int K_size,int max_depth,int max_dist,SuperRead_t *SuperRead, search_info *search_info)
{

	map<uint64_t,struct BFS_path_info_v2 > Visited_Path;
	map< uint64_t,int > stacked_nodes;
	
	int max_stack=500;
	int DepthTh=max_depth;//min(300/gap,20);
	int LenTh=max_dist;
	bool RIGHT=0;
	
	map<int , list<uint64_t> > dist_kmers;//neighborset
	dist_kmers[0].push_back(beg_kmer);
	int NBs=1;

	int dist_searched=0;
	uint64_t new_node=beg_kmer;
	uint64_t kmer,f_kmer;
	stacked_nodes[new_node]=1;
	//search direction

	Visited_Path[new_node].depth=1;
	Visited_Path[new_node].len=K_size;
	Visited_Path[new_node].last_kmer=NULL;
	
	map<int , list<uint64_t> >::iterator NB_it=dist_kmers.begin();
	while(1)
	{
		NB_it=dist_kmers.begin();
		
		if(NB_it==dist_kmers.end())
		{break;}
		if(NB_it->first>max_dist)
		{
			if(SuperRead->PathsCnt!=1)
			{return -10000;}
			else
			{break;}
		}
		if(NB_it->second.size()==0)
		{dist_kmers.erase(NB_it->first);continue;}
		if(NBs>max_stack)
		{
			
			SuperRead->PathsCnt=0;
			return -10000;
			break;
		}
		new_node=NB_it->second.front();

		NB_it->second.pop_front();
		NBs--;
		if(NB_it->second.size()==0)
		{
			dist_kmers.erase(NB_it->first);
			
		}		

		if(Visited_Path[new_node].depth>DepthTh||Visited_Path[new_node].len>LenTh)
		{continue;}
		
		uint64_t current_kmer=new_node;
		uint64_t current_kmer_rc=get_rev_comp_seq(current_kmer,K_size);

		if(current_kmer<current_kmer_rc)
		{
			RIGHT=1;
		}
		else
		{
			current_kmer=current_kmer_rc;
			RIGHT=0;
		}
		uint64_t hv=MurmurHash64A(&current_kmer,sizeof(current_kmer),0);
		uint64_t hash_idx=(size_t) (hv%ht->ht_sz);
		struct bucket** current_ptr;
		current_ptr= &(ht->store_pos[hash_idx]);
		bool found=look_up_in_a_list(current_kmer,&current_ptr);
		if(found==0||(*current_ptr)->kmer_info.removed==1)
		{
			continue;
		}
		if(RIGHT)
		{
			bool r_found=0;
			struct edge_node *edge_ptr;
			edge_ptr=(*current_ptr)->kmer_info.right;
			if(edge_ptr==NULL)
			{
				continue;
			}

			while(edge_ptr!=NULL)
			{
				kmer=current_kmer;					
				int edge_len=(int)(edge_ptr->len);
				for(int g=edge_len;g>=0;--g)
				{
					uint64_t b=(edge_ptr->edge>>2*g)&0x3;
					kmer&=(~((uint64_t)0x3<<(2*(K_size-1))));
					kmer<<=2;
					kmer|=b;
				}


				uint64_t kmer2,f_kmer2,next_kmer;
				next_kmer=kmer;
								
				if(next_kmer==end_kmer)
				{
					SuperRead->PathsCnt++;
					SuperRead->PathFound=1;
					int edge_len=edge_ptr->len;
					SuperRead->PathsLength=(int)(Visited_Path[new_node].len+edge_len+1);
						
				}

				// not in stack
				if(stacked_nodes[next_kmer]==0)
				{
					//Visited_Path[*ptr].cov=(int)(Visited_Path[new_node].cov+edge_ptr->edge_cov);
					Visited_Path[next_kmer].depth=(Visited_Path[new_node].depth+1);
						
					int edge_len=edge_ptr->len;
					int cum_len=(int)(Visited_Path[new_node].len+edge_len+1);
					Visited_Path[next_kmer].len=cum_len;//(Visited_Path[new_node].len+edge_ptr->len+1);

					Visited_Path[next_kmer].last_kmer=new_node;
					Visited_Path[next_kmer].last_bkt_edge=edge_ptr;
						
					stacked_nodes[next_kmer]=1;
					dist_kmers[cum_len].push_back(next_kmer);
					NBs++;

				}
				else
				{
					
					edge_ptr=edge_ptr->nxt_edge;
					continue;
						
			
				}

				edge_ptr=edge_ptr->nxt_edge;
		
			}

			


		}
		else
		{


			bool l_found=0;
			struct edge_node *edge_ptr;
			edge_ptr=(*current_ptr)->kmer_info.left;
			if(edge_ptr==NULL)
			{
				continue;
			}




			while(edge_ptr!=NULL)
			{
				kmer=current_kmer;					
				int edge_len=(int)(edge_ptr->len);
				for(int g=0;g<=edge_len;++g)
				{
					uint64_t b=(edge_ptr->edge>>2*g)&0x3;
					kmer>>=2;
					kmer|=(b<<(2*(K_size-1)));
				}

				kmer=get_rev_comp_seq(kmer,K_size);//take rc
				
				uint64_t next_kmer=kmer;


								
				if(next_kmer==end_kmer)
				{
					SuperRead->PathsCnt++;
					SuperRead->PathFound=1;
					int edge_len=edge_ptr->len;
					SuperRead->PathsLength=(int)(Visited_Path[new_node].len+edge_len+1);
						
				}

				// not in stack
				if(stacked_nodes[next_kmer]==0)
				{
					//Visited_Path[*ptr].cov=(int)(Visited_Path[new_node].cov+edge_ptr->edge_cov);
					Visited_Path[next_kmer].depth=(Visited_Path[new_node].depth+1);
						
					int edge_len=edge_ptr->len;
					int cum_len=(int)(Visited_Path[new_node].len+edge_len+1);
					Visited_Path[next_kmer].len=cum_len;//(Visited_Path[new_node].len+edge_ptr->len+1);

					Visited_Path[next_kmer].last_kmer=new_node;
					Visited_Path[next_kmer].last_bkt_edge=edge_ptr;
						
					
					stacked_nodes[next_kmer]=1;
					dist_kmers[cum_len].push_back(next_kmer);
					NBs++;

				}
				else
				{
					
					edge_ptr=edge_ptr->nxt_edge;
					continue;
						
			
				}

				edge_ptr=edge_ptr->nxt_edge;
		
			}


			

		}

	}


	if(SuperRead->PathsCnt!=1)
	{return -10000;}
	else
	{
		uint64_t  backtrack_kmer=end_kmer;
		int depth=Visited_Path[backtrack_kmer].depth+5;
		string extension;
		for (int i=depth;i>0;--i)
		{

			if(Visited_Path[backtrack_kmer].depth==1)
			{
				SuperRead->extension=extension;
				return SuperRead->PathsLength;		
			}
			edge_node* edge=Visited_Path[backtrack_kmer].last_bkt_edge;
			uint64_t edge_bits=edge->edge;
			uint64_t last_kmer=Visited_Path[backtrack_kmer].last_kmer;
			uint64_t last_kmer_rc=get_rev_comp_seq(last_kmer,K_size);
			if(last_kmer_rc<last_kmer)
			{
				edge_bits=get_rev_comp_seq(edge_bits,(edge->len+1));
			}

			char edge_cstr[1000];
			bitsarr2str(&edge_bits,(edge->len+1),edge_cstr,1);
			extension=edge_cstr+extension;
			backtrack_kmer=last_kmer;

		}
		return SuperRead->PathsLength;
	}

}


void CollectingNonContainedPairs(struct hashtable *ht1,struct hashtable *merge_ht1, struct hashtable2 *ht2, struct hashtable2 *merge_ht2,int K_size,vector<string>& filenames_vt,struct contigs_info * contigs_info,string ContigFilename)
{

	ofstream o_log;
	int lib_cnt=1;
	string pe_name;
	pe_name="InsertSizeEst.txt";
	o_log.open(pe_name.c_str());
	time_t beg_time,read_time;
	int64_t scaffold_len=0;
	string in_fname=ContigFilename;
	ifstream Contigs_in(in_fname.c_str());
	size_t num_Contigs=0;
	string tag,s,kmer;
	uint64_t f_seq,hv;
	struct kmer_t2 f_seq_t2;
	size_t hash_idx;
	//vector<int> contig_sz;
	bool found;
	bool flip_1,flip_2,flip_0;
	size_t ht_sz;
	
	ht_sz=ht1->ht_sz;
	size_t closed_gaps=0;
	size_t pair_gaps=0;
	size_t ambiguous_pairs=0;

	bool FAST=1;
	map<int,vector<int> > hit_position;
	contigs_info->contig_sz_vt.clear();
	contigs_info->contig_sz_vt.push_back(0);
	ofstream o_closed_gaps("Closed_gaps.txt");

	if(ContigFilename=="Contigs.txt")
	{
		ContigsRemapping(ht1,ht2, K_size, contigs_info,ContigFilename,0);	
	}	

	RemoveUnmappedNodes(ht1,ht2, K_size);
	//BuildContigAdjacency(ht1, ht2, contigs_info, K_size,  ContigFilename.c_str());
	///////////////marking finished
	if(K_size<=32)
	{
		AppendMergeHT(ht1, merge_ht1);
	}
	else
	{
		AppendMergeHT2(ht2,merge_ht2);
	}

	cout<<"Collecting paired ends information."<<endl;
	time(&beg_time);

	int64_t dist_sum=0,p_cnt=0;
	int mean_dist=0;
	int lib_no=0;
	for(size_t ii=0;ii<filenames_vt.size();ii+=2)
	{
		
		vector<int> insert_sz_est_vt;
		lib_no++;
		char o_sc_r_n[300],o_sc_l_n[300],o_sc_pd[300],o_sc_inward[300],o_sc_outward[300],o_sc_log[300],pair1_nc_name[300],pair2_nc_name[300];
		if(1)//Iter_Scaffold==0)
		{
			sprintf(o_sc_r_n,"Pdist_R_lib_%d.txt",lib_no);
			sprintf(o_sc_l_n,"Pdist_L_lib_%d.txt",lib_no);
			sprintf(o_sc_pd,"Pdist_lib_%d.txt",lib_no);
			sprintf(o_sc_inward,"Pdist_inward_lib_%d.txt",lib_no);
			sprintf(o_sc_outward,"Pdist_outward_lib_%d.txt",lib_no);
			sprintf(o_sc_log,"Pdist_lib_%d_log.txt",lib_no);
			sprintf(pair1_nc_name,"Lib_%d_NonContained_P1.txt",lib_no);
			sprintf(pair2_nc_name,"Lib_%d_NonContained_P2.txt",lib_no);
			
		}

		ofstream o_Pdist_R(o_sc_r_n),o_Pdist_L(o_sc_l_n),o_Pdist(o_sc_pd),o_Pdist_in(o_sc_inward),o_Pdist_out(o_sc_outward),o_Pdist_log(o_sc_log),o_pair1_nc(pair1_nc_name),o_pair2_nc(pair2_nc_name);
		uint64_t inward_found=0,inward_not_found=0,outward_found=0,outward_not_found=0;
		dist_sum=0,p_cnt=0;
		cout<<"Processing library: "<<ii<<" & "<<ii+1<<endl;
		//scaffold_len=insert_sz_vt[ii/2];
		int flip_cnt=0,non_flip_cnt=0;
		bool flip_library=0,flip_determined=0;
		int inward_cnt=0,outward_cnt=0;
		bool inward_library=0,orientation_determined=0,outward_library=0;
		bool SINGLE_READ=0;
		size_t n_pairs=0;
		int stable_insert_size=-10000;

		struct bucket **ptr1,**ptr2;
		bool RightSearch=0;
		for (int round=1;round<=3;++round)
		{
			if(round==2)
			{
				if(flip_determined==0)
				{
					cout<<"Bad library."<<endl;
					break;
				}
				if(flip_library==0)
				{
					continue;
				}
			}
		
			ifstream in_pair1,in_pair2;
	
			if(filenames_vt[ii]==filenames_vt[ii+1])
			{
				in_pair1.open(filenames_vt[ii].c_str());
				SINGLE_READ=1;
			}
			else
			{
				in_pair1.open(filenames_vt[ii].c_str());
				in_pair2.open(filenames_vt[ii+1].c_str());
				SINGLE_READ=0;
			}
		
			string t1,t2;
			int cont1,cont2;
			int64_t cod1,cod2;
			int readLen1,readLen2;
			size_t num_Reads=0;
			string kmer1,kmer2,seq_s1,seq_s2,tag_s1,tag_s1n,tag_s2,tag_s2n;
			bool fq_flag=0;
			string fq_tmp;
	
			struct read_t Read1,Read2;

			Read1.read_bits =(uint64_t*) malloc(100000+100);
			Read2.read_bits =(uint64_t*) malloc(100000+100);

			int nLines1=0,nLines2=0;
			bool read_success1=1,read_success2=1;
			read_success1=(bool) getline(in_pair1,t1);
		
			if(read_success1==0)
			{
				in_pair1.close();
				in_pair1.clear();
				continue;
			}


			if(fq_flag==0&&t1[0]=='@')
			{
				fq_flag=1;
			}

			in_pair1.close();
			in_pair1.clear();
			in_pair1.open(filenames_vt[ii].c_str());
		
			while(read_success1)
			{
				
			
				bool good_read1=1,good_read2=1;
				string QS_s1,QS_s2,tag1,tag2,n_tag1,n_tag2;
			

				if(fq_flag)
				{
					read_success1=get_a_fastq_read(in_pair1,tag1,seq_s1,QS_s1);
					read_success2=get_a_fastq_read(in_pair2,tag2,seq_s2,QS_s2);
				}
				else
				{
					read_success1=get_a_fasta_read(in_pair1,tag1,seq_s1,n_tag1);
					read_success2=get_a_fasta_read(in_pair2,tag2,seq_s2,n_tag2);
				}	

				good_read1=basic_quality_check(seq_s1);
				good_read2=basic_quality_check(seq_s2);

				readLen1=seq_s1.size();
				readLen2=seq_s2.size();

				if(readLen1<K_size)
				{good_read1=0;}
				if(readLen2<K_size)
				{good_read2=0;}

			

				int64_t pos1,pos2;
				uint64_t bits1,bits2;
			
				num_Reads++;

				if(num_Reads%10000000==0)
				{
					time(&read_time);
					cout<<num_Reads<<" Pairs Searched."<<endl;
					cout<<"Time: "<<difftime(read_time,beg_time)<<" secs."<<endl;
					continue;
				}

				if(good_read1==0||good_read2==0)
				{
					continue;
				}
		
				Init_Read(seq_s1,Read1);
				//seq_s1.clear();
			
		
				Init_Read(seq_s2,Read2);
				//seq_s2.clear();
			

				//base round

				bool output_current_pair=0;
				cont1=0,cont2=0;
				if(round==1)//first round determine if a reverse complement operation is necessary.
				{
					

					if(num_Reads>1000000)
					{
						cout<<"Bad library."<<endl;
					}

					bool found1=0,found2=0;
					struct bucket **ptr1_d,**ptr2_d;			
					bool flip_1d;
					for(int i=0;i<Read1.readLen-K_size+1;++i )
					{

						get_sub_arr(Read1.read_bits,Read1.readLen,i,K_size,&bits1);
						f_seq=get_rev_comp_seq(bits1,K_size);
						flip_0=0;
						if(bits1>f_seq)
						{
							bits1=f_seq;
							flip_0=1;
						}

						hv=MurmurHash64A(&bits1,sizeof(bits1),0);
						hash_idx=(size_t) (hv%ht_sz);
						struct bucket **ptr1;
						ptr1= &(ht1->store_pos[hash_idx]);
						found=look_up_in_a_list(bits1,&ptr1);
						if(found)
						{
					
							found1=1;
							if((*ptr1)->kmer_info.removed==0)
							{
								ptr1_d=ptr1;

								flip_1d=flip_0;
								cod1=(*ptr1)->kmer_info.cod;
								cont1=(*ptr1)->kmer_info.contig_no;
								flip_1=(*ptr1)->kmer_info.flip^flip_0;
								pos1=i;

								break;		
							}
						
					
						}
				


					}
					if(Read1.readLen<K_size||found1==0)
					{
						continue;
					}

					for(int i=0;i<=Read2.readLen-K_size;++i )
					{
				
						get_sub_arr(Read2.read_bits,Read2.readLen,i,K_size,&bits2);
						f_seq=get_rev_comp_seq(bits2,K_size);
						flip_0=0;
						if(bits2>f_seq)
						{
							bits2=f_seq;
							flip_0=1;
						}

						hv=MurmurHash64A(&bits2,sizeof(bits2),0);
						hash_idx=(size_t) (hv%ht_sz);
						struct bucket **ptr2;
						ptr2= &(ht1->store_pos[hash_idx]);
						found=look_up_in_a_list(bits2,&ptr2);
						if(found)
						{
						
							found2=1;
							if((*ptr2)->kmer_info.removed==0)
							{
								ptr2_d=ptr2;
								cod2=(*ptr2)->kmer_info.cod;
								cont2=(*ptr2)->kmer_info.contig_no;
								flip_2=(*ptr2)->kmer_info.flip^flip_0;
								pos2=i;

								break;		
						
							}
						}
					}
					if(Read2.readLen<K_size||found2==0)
					{continue;}

					if(cont1==0||cont2==0)
					{
						continue;
					}
					if(cont1==cont2)
					{
						if(flip_1==flip_2)
						{
							non_flip_cnt++;
						}
						else
						{
							flip_cnt++;
						}
					
					}

					
					if((flip_cnt>1000|non_flip_cnt>1000)&((flip_cnt>non_flip_cnt/2)|(flip_cnt/2<non_flip_cnt)))
					{
						flip_library=0;
						if(flip_cnt>non_flip_cnt/2)
						{
							flip_library=1;
							flip_determined=1;
							
						}
						cout<<"Initial Scan completed."<<endl;
						break;// break the reading process
					
					}
				}



				
				if(round==2)
				{
					

					if(num_Reads>1000000)
					{
						cout<<"bad library."<<endl;
					}

					bool found1=0,found2=0;
					struct bucket **ptr1_d,**ptr2_d;			
					bool flip_1d;
					for(int i=0;i<Read1.readLen-K_size+1;++i )
					{

						get_sub_arr(Read1.read_bits,Read1.readLen,i,K_size,&bits1);
						f_seq=get_rev_comp_seq(bits1,K_size);
						flip_0=0;
						if(bits1>f_seq)
						{
							bits1=f_seq;
							flip_0=1;
						}

						hv=MurmurHash64A(&bits1,sizeof(bits1),0);
						hash_idx=(size_t) (hv%ht_sz);
						struct bucket **ptr1;
						ptr1= &(ht1->store_pos[hash_idx]);
						found=look_up_in_a_list(bits1,&ptr1);
						if(found)
						{
					
							found1=1;
							if((*ptr1)->kmer_info.removed==0)
							{
								ptr1_d=ptr1;

								flip_1d=flip_0;
								cod1=(*ptr1)->kmer_info.cod;
								cont1=(*ptr1)->kmer_info.contig_no;
								flip_1=(*ptr1)->kmer_info.flip^flip_0;
								pos1=i;

								break;		
							}
						
					
						}
				


					}
					if(Read1.readLen<K_size||found1==0)
					{
						continue;
					}

					for(int i=0;i<=Read2.readLen-K_size;++i )
					{
				
						get_sub_arr(Read2.read_bits,Read2.readLen,i,K_size,&bits2);
						f_seq=get_rev_comp_seq(bits2,K_size);
						flip_0=0;
						if(bits2>f_seq)
						{
							bits2=f_seq;
							flip_0=1;
						}

						hv=MurmurHash64A(&bits2,sizeof(bits2),0);
						hash_idx=(size_t) (hv%ht_sz);
						struct bucket **ptr2;
						ptr2= &(ht1->store_pos[hash_idx]);
						found=look_up_in_a_list(bits2,&ptr2);
						if(found)
						{
						
							found2=1;
							if((*ptr2)->kmer_info.removed==0)
							{
								ptr2_d=ptr2;
								cod2=(*ptr2)->kmer_info.cod;
								cont2=(*ptr2)->kmer_info.contig_no;
								flip_2=(*ptr2)->kmer_info.flip^flip_0;
								pos2=i;

								break;		
						
							}
						}
					}
					if(Read2.readLen<K_size||found2==0)
					{continue;}

					if(cont1==0||cont2==0)
					{
						continue;
					}
					if(cont1==cont2)
					{
						if(flip_1==0&&flip_2==1)
						{
							if(cod2>cod1)
							{
								inward_cnt++;
							}
							else
							{
								outward_cnt++;
							}
						}
						else
						{
							if(cod2>cod1)
							{
								outward_cnt++;
							}
							else
							{
								inward_cnt++;
							}
						}
					
					}

					
					if((inward_cnt>1000|outward_cnt>1000)&((outward_cnt>inward_cnt/2)|(outward_cnt/2<inward_cnt)))
					{
						inward_library=0;
						outward_library=0;
						if(inward_cnt>outward_cnt/2)
						{
							inward_library=1;
							orientation_determined=1;
							
						}
						else
						{
							outward_library=1;
							orientation_determined=1;
						}
						cout<<"Orientation determined: ";
						if(inward_library==1)
						{ cout<< "inward."<<endl;}
						else
						{
							cout<< "outward."<<endl;
						}
						break;// break the reading process
					
					}
				}







				if(round==3)
				{
					search_info search_info;
					uint64_t beg_kmer,end_kmer;
					//cout<<"Collecting non-contained pairs..."<<endl;
					
					bool found1=0,found2=0;

					struct bucket **ptr1_d,**ptr2_d;
			
					if(inward_library)
					{
						int Read_arr_sz=Read2.readLen/32+1;
						int rem1=Read2.readLen%32;
						if(rem1==0)
						{Read_arr_sz--;}
						get_rev_comp_seq_arr(Read2.read_bits,Read2.readLen,Read_arr_sz);
					}
					
					if(outward_library)
					{
						int Read_arr_sz=Read1.readLen/32+1;
						int rem1=Read1.readLen%32;
						if(rem1==0)
						{Read_arr_sz--;}
						get_rev_comp_seq_arr(Read1.read_bits,Read1.readLen,Read_arr_sz);
					}

					bool flip_1d;
					for(int i=0;i<Read1.readLen-K_size+1;++i )
					{

						get_sub_arr(Read1.read_bits,Read1.readLen,i,K_size,&bits1);
						beg_kmer=bits1;
						f_seq=get_rev_comp_seq(bits1,K_size);
						flip_0=0;
						RightSearch=1;
						if(bits1>f_seq)
						{
							bits1=f_seq;
							flip_0=1;
							RightSearch=0;
						}

						hv=MurmurHash64A(&bits1,sizeof(bits1),0);
						hash_idx=(size_t) (hv%ht_sz);
						ptr1= &(ht1->store_pos[hash_idx]);
						found=look_up_in_a_list(bits1,&ptr1);
						if(found)
						{
					
							found1=1;
							if((*ptr1)->kmer_info.removed==0)
							{
								ptr1_d=ptr1;

								flip_1d=flip_0;
								cod1=(*ptr1)->kmer_info.cod;
								cont1=(*ptr1)->kmer_info.contig_no;
								flip_1=(*ptr1)->kmer_info.flip^flip_0;
								pos1=i;
								search_info.pos1=pos1;
								
								break;		
							}
							

						
					
						}
				


					}
					if(Read1.readLen<K_size||found1==0)
					{
						continue;
					}
					
					for(int i=Read2.readLen-K_size;i>=0;--i )
					{
				
						get_sub_arr(Read2.read_bits,Read2.readLen,i,K_size,&bits2);
						end_kmer=bits2;
						f_seq=get_rev_comp_seq(bits2,K_size);
						flip_0=0;
						if(bits2>f_seq)
						{
							bits2=f_seq;
							flip_0=1;
						}

						hv=MurmurHash64A(&bits2,sizeof(bits2),0);
						hash_idx=(size_t) (hv%ht_sz);
						
						ptr2= &(ht1->store_pos[hash_idx]);
						found=look_up_in_a_list(bits2,&ptr2);
						if(found)
						{
				
							found2=1;
							if((*ptr2)->kmer_info.removed==0)
							{
								ptr2_d=ptr2;
								cod2=(*ptr2)->kmer_info.cod;
								cont2=(*ptr2)->kmer_info.contig_no;
								flip_2=(*ptr2)->kmer_info.flip^flip_0;
								pos2=i;
								search_info.pos2=pos2;
								search_info.Flip_End=flip_0;
								break;		
							}
							
						}
					}
					if(Read2.readLen<K_size||found2==0)
					{continue;}

					if(cont1==0||cont2==0)
					{
						continue;
					}
					bool ignore=0;

					if(cont1==cont2&&flip_1==flip_2)
					{
						
				
						int pdist=(int)(abs(cod2-cod1)+pos1+readLen2-pos2);
						dist_sum+=pdist;
						o_Pdist<<pdist<<endl;
						p_cnt++;
					
						if(stable_insert_size>0&(abs(pdist-stable_insert_size)>stable_insert_size/3))
						{output_current_pair=1;}
						if(p_cnt>1000)
						{
							stable_insert_size=dist_sum/p_cnt+1;
							if(stable_insert_size<0)
							{
								cout<<"Bad library."<<endl;
								break;//quit
							}
						}

					}

					if(cont1!=cont2)
					{
						output_current_pair=1;
				
						int pdist=-10000;

						int64_t d1,d2,d3,d4,dist;

						// build adjacency info for cont1&2
						if(flip_1==0)
						{
							//struct adjacent_contig_info adj_contig;
							int Tcontig_no=cont2;
							int Tflip=flip_2;
							int Tadjacent=0;
							if(flip_2==0)
							{
								//int64_t Tdist=(scaffold_len-readLen2+pos2)-((int)contigs_info->contig_sz_vt[cont1]-cod1)-cod2;////
								dist=d1=	(scaffold_len-(int64_t)readLen2+pos2-pos1)-((int64_t)contigs_info->contig_sz_vt[cont1]-cod1)-cod2;
						
								o_Pdist_R<<cont1 <<" "<<cont2<<" "<<d1<<endl;//" Library: "<<lib_no<<endl;
						//		o_Pdist_L<<cont2 <<" "<<cont1<<" "<<d1<<" Library: "<<ii<<endl;

							}
							else
							{
								int64_t r_cod2=(int64_t)(contigs_info->contig_sz_vt[cont2])-cod2-(int64_t)K_size;
								dist=d2=(scaffold_len+pos2-(int64_t)readLen2-pos1)-((int64_t)contigs_info->contig_sz_vt[cont1]-cod1)-r_cod2;
								//int64_t Tdist=(scaffold_len+pos2-readLen2)-((int)contigs_info->contig_sz_vt[cont1]-cod1)-r_cod2;////
						
								o_Pdist_R<<cont1 <<" "<<-cont2<<" "<<d2<<endl;//" Library: "<<lib_no<<endl;
							//	o_Pdist_R<<cont2 <<" "<<-cont1<<" "<<d2<<" Library: "<<lib_no<<endl;


							}


						}
						else
						{

							//struct adjacent_contig_info adj_contig;
							int Tcontig_no=cont2;
							bool Tflip=flip_2;

							if(flip_2==0)
							{

								//int64_t Tdist=(scaffold_len-readLen2+pos2)-cod2-cod1;////
								dist=d3=(scaffold_len-(int64_t)readLen2+pos2-pos1)-cod2-cod1-(int64_t)K_size;
								o_Pdist_L<<cont1 <<" "<<-cont2<<" "<<d3<<endl;//" Library: "<<lib_no<<endl;
								//o_Pdist_L<<cont2 <<" "<<-cont1<<" "<<d3<<" Library: "<<lib_no<<endl;

							}
							else
							{
								int64_t r_cod2=(int64_t)contigs_info->contig_sz_vt[cont2]-cod2-(int64_t)K_size;
								//int64_t Tdist=(scaffold_len-readLen2+pos2)-cod1-r_cod2;////
								dist=d4=(scaffold_len-(int64_t)readLen2+pos2-pos1)-cod1-(int64_t)K_size-r_cod2;

								o_Pdist_L<<cont1 <<" "<<cont2<<" "<<d4<<endl;//" Library: "<<lib_no<<endl;
								//o_Pdist_R<<cont2 <<" "<<cont1<<" "<<d4<<" Library: "<<lib_no<<endl;

							}

						}




					}



					if(output_current_pair)
					{
						pair_gaps++;
						char pair1_cstr[1000],pair2_cstr[1000];
						int Read1_arr_sz=Read1.readLen/32+1;
						int rem=Read1.readLen%32;
						if(rem==0)
						{Read1_arr_sz--;}
						bitsarr2str(Read1.read_bits,Read1.readLen,pair1_cstr,Read1_arr_sz);
						
						int Read2_arr_sz=Read2.readLen/32+1;
						rem=Read2.readLen%32;
						if(rem==0)
						{Read2_arr_sz--;}
						bitsarr2str(Read2.read_bits,Read2.readLen,pair2_cstr,Read2_arr_sz);
						tag1[0]='>';
						
						tag2[0]='>';
						o_pair1_nc<<tag1<<endl<<pair1_cstr<<endl;
						o_pair2_nc<<tag2<<endl<<pair2_cstr<<endl;

						stacked_bucket kmer_stack_beg;
						int max_depth=300;
						int max_dist=250;

						kmer_stack_beg.bktptr=*ptr1;
						kmer_stack_beg.RightSearch=RightSearch;
						if(!RightSearch)
						{
							//search_info.Flip_End=!search_info.Flip_End;
						}
						SuperRead_t SuperRead;
						SuperRead.PathsCnt=0;
						SuperRead.depth=0;
						SuperRead.PathFound=0;
						SuperRead.PathsLength=0;

						int offset=(search_info.pos1+readLen2-search_info.pos2);
						search_info.offset=offset;
						SuperRead.extension.clear();
						char temp_cstr[1000];
						uint64_t temp_kmer=(*ptr1)->kmer_t.kmer;
						if(!RightSearch)
						{
							temp_kmer=get_rev_comp_seq(temp_kmer,K_size);
						}

						bitsarr2str(&temp_kmer,K_size,temp_cstr,1);
						
						string SuperReadStr;
						//if(num_Reads==84)
						//cout<<num_Reads<<endl;
//						int gap_dist = BFSearchGapCloser(ht1,merge_ht1, *ptr1, *ptr2,K_size, kmer_stack_beg, max_depth,max_dist,&SuperRead,&search_info);
						
						int gap_dist = BFSearchGapCloser_v2(ht1,merge_ht1, beg_kmer, end_kmer,K_size, max_depth,max_dist,&SuperRead,&search_info);
						if(SuperRead.extension.size()>0)
						{
							SuperReadStr=seq_s1.substr(0,pos1)+temp_cstr+SuperRead.extension+seq_s2.substr(pos2+K_size,seq_s2.size());
							o_closed_gaps<<tag1<<endl<<SuperReadStr<<endl;
						}
						
						if(gap_dist>0)
						{
							closed_gaps++;
						}
						if(SuperRead.PathsCnt>1)
						{
							ambiguous_pairs++;
						}
						


					}
				
				}





			}
			

		}
		cout<<endl<<pair_gaps<< " gaps in total."<<endl;
		cout<<closed_gaps<<" closed gaps."<<endl;
		cout<<ambiguous_pairs<<" ambiguous pairs."<<endl;
		

		o_Pdist_log<<"Inward found:"<<inward_found<<endl;
		o_Pdist_log<<"Inward not found:"<<inward_not_found<<endl;		
		o_Pdist_log<<"Outward found:"<<outward_found<<endl;
		o_Pdist_log<<"Outward not found:"<<outward_not_found<<endl;

		if(p_cnt>0)
		{o_log<<"InsertSizeEst_"<<lib_no<<": "<<int((double)dist_sum/(double)p_cnt+0.5)<<endl;}
		else
		{
		o_log<<"InsertSizeEst_"<<lib_no<<": "<<-10000<<endl;
		}



	}
	
	

}

/*
struct overlapper_info
{
	int search_depth;
	bool forward_end_search;

};

void FindOverlapsForEachRead(reads_table * reads_table, hashtable *ht, int K_size,read_index *read_index)
{
	//1. check the good region of the read --good kmers
	//	find the first and last sparse k-mer, and allow small offsets.
	//2. check if the reads overlap.
	//3. output the overlap
	
	//scan through the end to find the first  last overlap between reads and check if there are other matching kmers

	int search_depth=100;
	reads_table->read_ptr[1];
	uint64_t *read_arr1,*read_arr2;
	uint64_t kmer1,kmer2;
	read_arr1=reads_table->read_ptr[1];
	int readLen1,readLen2,read_idx1,read_idx2;

	read_idx1=1;
	readLen1=reads_table->read_len_vt[read_idx1-1];//double check;

	int Read_arr_sz1=reads_table->read_len_vt[read_idx1]/32+1;
	int rem=reads_table->read_len_vt[read_idx1]%32;
	if(rem==0)
	{Read_arr_sz1--;}

	for (int i=0;i<search_depth;++i)
	{
		//check for prefix overlap
		get_sub_arr(read_arr1,Read_arr_sz1,i,K_size,&kmer1);

		uint64_t f_kmer1=get_rev_comp_seq(kmer1,K_size);
		bool flip_read1=0;
		if(kmer1>f_kmer1)
		{
			uint64_t t=kmer1;
			kmer1=f_kmer1;
			f_kmer1=t;
			flip_read1=1;
			//flip_beg=1;
		}
	
		uint64_t hv=MurmurHash64A(&kmer1,sizeof(kmer1),0);

		size_t hash_idx=(size_t) (hv%ht->ht_sz);

		bucket ** bktp2p= &(ht->store_pos[hash_idx]);

		bool found=look_up_in_a_list(kmer1,&bktp2p);
		if(found)
		{
			int idx_pos=(*bktp2p)->kmer_info.cod;
			//map_it=read_index->repeat_maps[idx_pos].end();++map_it
			//read_index->repeat_maps[idx_pos].begin()
			for(map<uint64_t,bool>::iterator map_it=read_index->repeat_maps[idx_pos].begin();map_it!=read_index->repeat_maps[idx_pos].end();++map_it)
			{
				int read_idx2=(uint32_t) map_it->first;
				int read_cod2=(map_it->first>>32);
				if(read_idx2==read_idx1)
				{continue;}

				//check each of the reads, if overlap, do something, if the current read is contained, throw it away from future analysis
				if(flip_read1^map_it->second)
				{
					
				}
				else
				{
					//...end2 <end1 & head2<head1 then read2 is an overlap candidate. it is a solid overlap if it doesn't go into another contig.
					
				}

			}
		}
		//look up in a table for the key...
		//get_sub_arr to check each kmer
		//get reads that match and put them into a list.

		//for each read, check the overlap relation with the current one.
		//1. overlap
		//2. do not overlap
		//3. contain
	}
}
*/

#endif
