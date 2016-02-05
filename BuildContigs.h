#ifndef __BUILD_CONTIGS_H
#define __BUILD_CONTIGS_H


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
using namespace std;

//mark the branches in the sparse kmer graph
void MarkBranches(hashtable *ht)
{
	for(size_t i=0;i<ht->ht_sz;++i)
	{


		struct bucket *bkt_ptr=ht->store_pos[i];

		while(bkt_ptr!=NULL)
		{

			bkt_ptr->kmer_info.used=0;
			bkt_ptr->kmer_info.marked=0;
			



			if(bkt_ptr->kmer_info.left!=NULL&&bkt_ptr->kmer_info.left->nxt_edge!=NULL)
			{
				bkt_ptr->kmer_info.split_left=1;
			}
			else
			{
				bkt_ptr->kmer_info.split_left=0;
			}

			if(bkt_ptr->kmer_info.right!=NULL&&bkt_ptr->kmer_info.right->nxt_edge!=NULL)
			{
				bkt_ptr->kmer_info.split_right=1;
			}
			else
			{
				bkt_ptr->kmer_info.split_right=0;
			}
			bkt_ptr=bkt_ptr->nxt_bucket;
		}
	}
}


void MarkBranches2(hashtable2 *ht)
{
	for(size_t i=0;i<ht->ht_sz;++i)
	{


		struct bucket2 *bkt_ptr=ht->store_pos[i];

		while(bkt_ptr!=NULL)
		{
			bkt_ptr->kmer_info.used=0;
			bkt_ptr->kmer_info.marked=0;
			

			if(bkt_ptr->kmer_info.left!=NULL&&bkt_ptr->kmer_info.left->nxt_edge!=NULL)
			{
				bkt_ptr->kmer_info.split_left=1;
			}
			else
			{
				bkt_ptr->kmer_info.split_left=0;
			}

			if(bkt_ptr->kmer_info.right!=NULL&&bkt_ptr->kmer_info.right->nxt_edge!=NULL)
			{
				bkt_ptr->kmer_info.split_right=1;
			}
			else
			{
				bkt_ptr->kmer_info.split_right=0;
			}
			bkt_ptr=bkt_ptr->nxt_bucket;
		}
	}
}

void MarkBranches3(hashtable3 *ht)
{
	for(size_t i=0;i<ht->ht_sz;++i)
	{

		struct bucket3 *bkt_ptr=ht->store_pos[i];

		while(bkt_ptr!=NULL)
		{
			bkt_ptr->kmer_info.used=0;
			bkt_ptr->kmer_info.marked=0;
			

			if(bkt_ptr->kmer_info.left!=NULL&&bkt_ptr->kmer_info.left->nxt_edge!=NULL)
			{
				bkt_ptr->kmer_info.split_left=1;
			}
			else
			{
				bkt_ptr->kmer_info.split_left=0;
			}

			if(bkt_ptr->kmer_info.right!=NULL&&bkt_ptr->kmer_info.right->nxt_edge!=NULL)
			{
				bkt_ptr->kmer_info.split_right=1;
			}
			else
			{
				bkt_ptr->kmer_info.split_right=0;
			}
			bkt_ptr=bkt_ptr->nxt_bucket;
		}
	}
}

void MarkBranches4(hashtable4 *ht)
{
	for(size_t i=0;i<ht->ht_sz;++i)
	{

		struct bucket4 *bkt_ptr=ht->store_pos[i];

		while(bkt_ptr!=NULL)
		{
			bkt_ptr->kmer_info.used=0;
			bkt_ptr->kmer_info.marked=0;
			

			if(bkt_ptr->kmer_info.left!=NULL&&bkt_ptr->kmer_info.left->nxt_edge!=NULL)
			{
				bkt_ptr->kmer_info.split_left=1;
			}
			else
			{
				bkt_ptr->kmer_info.split_left=0;
			}

			if(bkt_ptr->kmer_info.right!=NULL&&bkt_ptr->kmer_info.right->nxt_edge!=NULL)
			{
				bkt_ptr->kmer_info.split_right=1;
			}
			else
			{
				bkt_ptr->kmer_info.split_right=0;
			}
			bkt_ptr=bkt_ptr->nxt_bucket;
		}
	}
}


void MarkBranches0(hashtable0 *ht)
{
	for(size_t i=0;i<ht->ht_sz;++i)
	{

		struct bucket0 *bkt_ptr=ht->store_pos[i];

		while(bkt_ptr!=NULL)
		{
			bkt_ptr->kmer_info.used=0;
			bkt_ptr->kmer_info.marked=0;
			

			if(bkt_ptr->kmer_info.left!=NULL&&bkt_ptr->kmer_info.left->nxt_edge!=NULL)
			{
				bkt_ptr->kmer_info.split_left=1;
			}
			else
			{
				bkt_ptr->kmer_info.split_left=0;
			}

			if(bkt_ptr->kmer_info.right!=NULL&&bkt_ptr->kmer_info.right->nxt_edge!=NULL)
			{
				bkt_ptr->kmer_info.split_right=1;
			}
			else
			{
				bkt_ptr->kmer_info.split_right=0;
			}
			bkt_ptr=bkt_ptr->nxt_bucket;
		}
	}
}



// produce the contigs, single end assembly complete
void build_contigs(struct hashtable *ht,int K_size, int gap,string Contig_Filename,bool ScreenOffTips)
{
	//bool ScreenOffTips=1;
	int TipLenTh=100;
	int TipCovTh=2;
	string o_name=Contig_Filename;
	ofstream o_contigs(o_name.c_str());
	ofstream o_contigsLen("Contigs_Len.txt");
	ofstream o_contigsCov("Contigs_Cov.txt");
	o_name="nLB_"+Contig_Filename;
	ofstream o_nLB(o_name.c_str());
	o_name="nRB_"+Contig_Filename;
	ofstream o_nRB(o_name.c_str());
	o_name="Cov_Hist"+Contig_Filename;
	ofstream o_cov_hist(o_name.c_str());

	map<int,int> cov_hist;
	bool FE_left=0,FE_right=0;
	ofstream o_Log("ContigsLog.txt");//,ios_base::app);
	bool COVERAGE_STS=1;
	uint64_t coverage=0;
	uint64_t kmer_cnt=0;
	size_t i=0,ht_sz=ht->ht_sz;
	size_t num_contigs=0;
	string s_left,s_right;
	s_left.reserve(100000);
	s_right.reserve(100000);
	struct bucket *bkt_ptr ,*beg_bkt_ptr,*last_bkt_ptr;
	uint64_t beg_kmer,t_kmer,t,f_kmer,hv,hash_idx;
	o_name="Contigs_info.txt";
	ofstream o_contig_info(o_name.c_str());//"Contigs_info.txt");
	//char c_str[1000];
	bool found;
//	vector<string> contigs_vt;
	for(i=0;i<ht_sz;++i)
	{
		bkt_ptr=ht->store_pos[i];

		while(bkt_ptr!=NULL)
		{
			beg_bkt_ptr=bkt_ptr;
			if(bkt_ptr->kmer_info.used==0&&bkt_ptr->kmer_info.removed==0)
			{
				num_contigs++;
				//if(num_contigs==9&&bkt_ptr->kmer_t.kmer==3550978978893533640)
				//{cout<<"";}
			//	if(bkt_ptr->kmer_t.kmer==14437493651275778)
			//	{cout<<"";}
				FE_left=0;FE_right=0;

				if(COVERAGE_STS)
				{
					coverage=0;
					kmer_cnt=0;
					coverage+=bkt_ptr->kmer_info.cov1;
					kmer_cnt+=1;
				}

				string contig;
				beg_kmer=bkt_ptr->kmer_t.kmer;
				bkt_ptr->kmer_info.used=1;

				//bkt_ptr->kmer_info.strand_visited|=1;
				bkt_ptr->kmer_info.flip=0;
				t_kmer=beg_kmer;
				char c_seq[1000];

				s_left.clear();
				s_right.clear();
				
				int nLBranches=0,nRBranches=0;
				uint64_t left_bits,right_bits;

				bool Right=0,flag=0;
				
				
				for(int it=1;it<=2;++it)
				{
					string sum_str,t_str;
					
					bool Free_End=0;
					int nBranches=0;
					sum_str.clear();
					t_str.clear();
					//right search
				
					if(it==1)
					{
						Right=1;
					}
					else
					{
						Right=0;
					}
					bkt_ptr=beg_bkt_ptr;
					t_kmer=bkt_ptr->kmer_t.kmer;

					while(1)
					{

						sum_str+=t_str;
						t_str.clear();
						bkt_ptr->kmer_info.used=1;
						if(Right==1)
						{
						
							if(bkt_ptr->kmer_info.right==NULL)
							{
								Free_End=1;
								nBranches=0;
								break;
							}
							if(bkt_ptr->kmer_info.right->nxt_edge!=NULL)
							{
								nBranches=0;
								edge_node *edge_ptr=bkt_ptr->kmer_info.right;
								while(edge_ptr!=NULL)
								{
									nBranches++;
									edge_ptr=edge_ptr->nxt_edge;
								}
							
								break;
							}
							int edge_len=(bkt_ptr->kmer_info.right)->len;
							for(int j=edge_len;j>=0;--j)
							{

								right_bits=(bkt_ptr->kmer_info.right->edge);
								right_bits>>=(2*j);
								switch(right_bits&0x3)
								{
									case 0:
										t=3;
										t<<=(K_size-1)*2;
										t_kmer&=(~t);
										t_kmer<<=2;

										t_str.push_back('A');

										break;
									case 1:
										t=3;
										t<<=(K_size-1)*2;
										t_kmer&=(~t);
										t_kmer<<=2;
										t=1;

										t_kmer|=t;
										t_str.push_back('C');
										break;
									case 2:
										t=3;
										t<<=(K_size-1)*2;
										t_kmer&=(~t);
										t_kmer<<=2;
										t=2;

										t_kmer|=t;
										t_str.push_back('G');
										break;
									case 3:
										t=3;
										t<<=(K_size-1)*2;
										t_kmer&=(~t);
										t_kmer<<=2;
										t=3;

										t_kmer|=t;
										t_str.push_back('T');
										break;
								}
							}

						}
						else
						{
		
							if(bkt_ptr->kmer_info.left==NULL)
							{
								Free_End=1;
								nBranches=0;
								break;
							}
							if(bkt_ptr->kmer_info.left->nxt_edge!=NULL)
							{
								nBranches=0;
								edge_node *edge_ptr=bkt_ptr->kmer_info.left;
								while(edge_ptr!=NULL)
								{
									nBranches++;
									edge_ptr=edge_ptr->nxt_edge;
								}
								break;
							}
							int edge_len=bkt_ptr->kmer_info.left->len;
							for(int j=0;j<=edge_len;++j)
							{

								left_bits=(bkt_ptr->kmer_info.left)->edge;
								left_bits>>=(2*j);
								switch(left_bits&0x3)
								{
									case 0:
										t_kmer>>=2;
										t=((uint64_t)0)<<((K_size-1)*2);
										t_kmer|=t;
										t_str.push_back('T');

										break;
									case 1:
										t_kmer>>=2;
										t=((uint64_t)1)<<((K_size-1)*2);
										t_kmer|=t;
										t_str.push_back('G');
										break;
									case 2:
										t_kmer>>=2;
										t=((uint64_t)2)<<(K_size-1)*2;
										t_kmer|=t;
										t_str.push_back('C');
										break;
									case 3:
										t_kmer>>=2;
										t=((uint64_t)3)<<((K_size-1)*2);
										t_kmer|=t;
										t_str.push_back('A');
										break;

								}
							}

						}


						f_kmer=get_rev_comp_seq(t_kmer,K_size);
						if(t_kmer>f_kmer)
						{
							t_kmer=f_kmer;
							Right=!Right;
						}

						hv=MurmurHash64A(&t_kmer,sizeof(t_kmer),0);

						hash_idx=(size_t) (hv%ht_sz);

						struct bucket ** ptr;

						ptr= &(ht->store_pos[hash_idx]);
						last_bkt_ptr=bkt_ptr;
						found=look_up_in_a_list(t_kmer,&ptr);

						bkt_ptr=*ptr;
					
						if(found==0)
						{
							Free_End=1;
							
							flag=1;

							break;
						}


	////////////////////////////////////////////////////////////////////////////////////

						if((Right==1&&bkt_ptr->kmer_info.split_left==1)||(Right==0&&bkt_ptr->kmer_info.split_right==1))
						{
							flag=1;
							last_bkt_ptr->kmer_info.marked=1;
							//t_str=t_str.substr(0,t_str.size()-1);
							t_str.clear();
							sum_str+=t_str;
							t_str.clear();
							nBranches=1;

							break;
						}
	////////////////////////////////////////////////////////////////////////////////////


						if(bkt_ptr->kmer_info.used==1)
						{
							flag=1;
							break;
						}

						if(COVERAGE_STS)
						{
							coverage+=bkt_ptr->kmer_info.cov1;
							kmer_cnt++;

						}
					
						if((Right==0&&it==1)||(Right==1&&it==2))
						{
							bkt_ptr->kmer_info.flip=1;
						}
						else
						{
							bkt_ptr->kmer_info.flip=0;
						}


					}
					if(it==1)
					{
						nRBranches=nBranches;
						s_right=sum_str;
						t_str.clear();
						sum_str.clear();
					}
					if(it==2)
					{
						nLBranches=nBranches;
						s_left=sum_str;
						reverse(s_left.begin(),s_left.end());
						complement_str(s_left);
						t_str.clear();
						sum_str.clear();
						
					}

					
				}
				//left,rght search ended

				bitsarr2str(&beg_kmer, K_size, c_seq,1);
			
				contig=s_left+c_seq+s_right;
				
				if(ScreenOffTips)
				{
					if((int)(contig.size())<TipLenTh&&((double(coverage)/double(kmer_cnt))<TipCovTh))
					{
						num_contigs--;
						bkt_ptr=beg_bkt_ptr->nxt_bucket;
						continue;
					}
				}


				//contigs_vt.push_back(contig);
				o_contigs<<">Contig_"<<num_contigs;//<< " ";
				if(0)//COVERAGE_STS)
				{
					o_contigs<<"Avg_Cov: "<< (double(coverage)/double(kmer_cnt));
//					o_contigs<<" Kmer_cnt: "<<kmer_cnt;
					o_contigs<<" Contig_size: "<<contig.size();
				}
				o_contigsLen<<contig.size()<<endl;
				o_contigsCov<<(int)(double(coverage)/double(kmer_cnt)+0.5)<<endl;
				cov_hist[int((double(coverage)/double(kmer_cnt)+0.5))]++;


				o_contigs<<endl;
				o_contigs<<contig<<endl;
			
				o_nLB<<nLBranches<<endl;
				o_nRB<<nRBranches<<endl;

				o_contig_info<<">Contig "<<num_contigs<< " ";
				if(COVERAGE_STS)
				{

					o_contig_info<<"Avg_Cov: "<< (double(coverage)/double(kmer_cnt));
					o_contig_info<<" Kmer_cnt: "<<kmer_cnt;
					o_contig_info<<" Contig_size: "<<contig.size()<<endl;

				}
			}
			bkt_ptr=beg_bkt_ptr->nxt_bucket;
		}

	}
	if(1)
	{
		o_Log<<"Contigs count:"<<endl;
		o_Log<<num_contigs<<endl;
	}

	map<int,int>::iterator mit;
	for(mit=cov_hist.begin();mit!=cov_hist.end();++mit)
	{
		o_cov_hist<<mit->first<<" "<<mit->second<<endl;
	}

}

void build_contigs2(struct hashtable2 *ht,int K_size, int gap,string Contig_Filename,bool ScreenOffTips)
{
	//bool ScreenOffTips=0;
	int TipLenTh=K_size+1;
	int TipCovTh=1;
	ofstream o_contigsLen("Contigs_Len.txt");
	ofstream o_contigsCov("Contigs_Cov.txt");
	string o_name=Contig_Filename;
	ofstream o_contigs(o_name.c_str());
	//o_name=Contig_Filename+"_FE_left.txt";
	//ofstream o_FE_left(o_name.c_str());
	o_name="nLB_"+Contig_Filename;
	ofstream o_nLB(o_name.c_str());
	//o_name=Contig_Filename+"_FE_right.txt";
	//ofstream o_FE_right(o_name.c_str());
	o_name="nRB_"+Contig_Filename;
	ofstream o_nRB(o_name.c_str());
	o_name="Cov_Hist"+Contig_Filename;
	ofstream o_cov_hist(o_name.c_str());

	map<int,int> cov_hist;
	bool FE_left=0,FE_right=0;
	ofstream o_Log("ContigsLog.txt");//,ios_base::app);
	bool COVERAGE_STS=1;
	uint64_t coverage=0;
	uint64_t kmer_cnt=0,t,hv;
	size_t i=0,ht_sz=ht->ht_sz,hash_idx;
	size_t num_contigs=0;
	string s_left,s_right;
	s_left.reserve(100000);
	s_right.reserve(100000);
	struct bucket2 *bkt_ptr ,*beg_bkt_ptr,*last_bkt_ptr;
	kmer_t2 beg_kmer,t_kmer,f_kmer;
	o_name="Contigs_info.txt";
	ofstream o_contig_info(o_name.c_str());//"Contigs_info.txt");
	//char c_str[1000];
	bool found;
//	vector<string> contigs_vt;
	for(i=0;i<ht_sz;++i)
	{
		bkt_ptr=ht->store_pos[i];

		while(bkt_ptr!=NULL)
		{
			beg_bkt_ptr=bkt_ptr;
			if(bkt_ptr->kmer_info.used==0&&bkt_ptr->kmer_info.removed==0)
			{
				num_contigs++;
				FE_left=0;FE_right=0;

				if(COVERAGE_STS)
				{
					coverage=0;
					kmer_cnt=0;
					coverage+=bkt_ptr->kmer_info.cov1;
					kmer_cnt+=1;
				}

				string contig;
				beg_kmer=bkt_ptr->kmer_t2;
				bkt_ptr->kmer_info.used=1;

				//bkt_ptr->kmer_info.strand_visited|=1;
				bkt_ptr->kmer_info.flip=0;
				t_kmer=beg_kmer;


				char c_seq[1000];

				s_left.clear();
				s_right.clear();
				
				uint64_t left_bits,right_bits;
				int nLBranches=0,nRBranches=0;
				bool Right=0;
				

				for (int it=1;it<=2;++it)
				{
					//1:right search,2:left;
					bool flag=0;
					bool Free_End=0;
					bkt_ptr=beg_bkt_ptr;
					t_kmer=bkt_ptr->kmer_t2;
					string sum_str, t_str;
					
					int nBranches=0;
	
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
						sum_str+=t_str;
						t_str.clear();
						bkt_ptr->kmer_info.used=1;
						if(Right==1)
						{
							if(bkt_ptr->kmer_info.right==NULL)
							{
								Free_End=1;
								nBranches=0;
								break;
							}
						
							if(bkt_ptr->kmer_info.right->nxt_edge!=NULL)
							{
								edge_node *edge_ptr=bkt_ptr->kmer_info.right;
								nBranches=0;
								while(edge_ptr!=NULL)
								{
									nBranches++;
									edge_ptr=edge_ptr->nxt_edge;
								}
								break;
							}
							//no branch
							for(int j=(bkt_ptr->kmer_info.right)->len;j>=0;--j)
							{

								right_bits=(bkt_ptr->kmer_info.right->edge)>>2*j;

								switch(right_bits&0x3)
								{
									case 0:
										t=3;
										t<<=(K_size-1-32)*2;
										t_kmer.kmer[0]&=(~t);
										L_shift_NB(t_kmer.kmer,2,2);

										t_str.push_back('A');

										break;
									case 1:
										t=3;
										t<<=(K_size-1-32)*2;
										t_kmer.kmer[0]&=(~t);
										L_shift_NB(t_kmer.kmer,2,2);
										t=1;

										t_kmer.kmer[1]|=t;
										t_str.push_back('C');
										break;
									case 2:
										t=3;
										t<<=(K_size-1-32)*2;
										t_kmer.kmer[0]&=(~t);
										L_shift_NB(t_kmer.kmer,2,2);
										t=2;

										t_kmer.kmer[1]|=t;
										t_str.push_back('G');
										break;
									case 3:
										t=3;
										t<<=(K_size-1-32)*2;
										t_kmer.kmer[0]&=(~t);
										L_shift_NB(t_kmer.kmer,2,2);
										t=3;

										t_kmer.kmer[1]|=t;
										t_str.push_back('T');
										break;



								}
							}


							
						}
						else
						{
					//left search
							if(bkt_ptr->kmer_info.left==NULL)
							{
								Free_End=1;
								nBranches=0;
								break;
							}
							if(bkt_ptr->kmer_info.left->nxt_edge!=NULL)
							{
								edge_node *edge_ptr=bkt_ptr->kmer_info.left;
								nBranches=0;
								while(edge_ptr!=NULL)
								{
									nBranches++;
									edge_ptr=edge_ptr->nxt_edge;
								}
								break;//or left=...;
							}
							int edge_len=bkt_ptr->kmer_info.left->len;

							for(int j=0;j<=edge_len;++j)
							{

								left_bits=(bkt_ptr->kmer_info.left)->edge>>2*j;

								switch(left_bits&0x3)
								{
									case 0:
										R_shift_NB(t_kmer.kmer,2,2);

										t_str.push_back('T');

										break;
									case 1:
										R_shift_NB(t_kmer.kmer,2,2);

										t=((uint64_t)1)<<((K_size-1-32)*2);
										t_kmer.kmer[0]|=t;
										t_str.push_back('G');
										break;
									case 2:
										R_shift_NB(t_kmer.kmer,2,2);

										t=((uint64_t)2)<<((K_size-1-32)*2);
										t_kmer.kmer[0]|=t;
										t_str.push_back('C');
										break;
									case 3:
										R_shift_NB(t_kmer.kmer,2,2);

										t=((uint64_t)3)<<((K_size-1-32)*2);
										t_kmer.kmer[0]|=t;
										t_str.push_back('A');
										break;

								}
							}

						}


						f_kmer=t_kmer;
						get_rev_comp_seq_arr(f_kmer.kmer,K_size,2);
						if(uint64_t_cmp(t_kmer.kmer,f_kmer.kmer,2)>0)
						{
							t_kmer=f_kmer;
							Right=!Right;
						}


						hv=MurmurHash64A(&t_kmer,sizeof(t_kmer),0);

						hash_idx=(size_t) (hv%ht_sz);


						struct bucket2 ** ptr;

						ptr= &(ht->store_pos[hash_idx]);
						last_bkt_ptr=bkt_ptr;
						found=look_up_in_a_list2(&t_kmer,&ptr);

						bkt_ptr=*ptr;
					

						if(found==0)
						{
							Free_End=1;
							flag=1;

							break;
						}


	////////////////////////////////////////////////////////////////////////////////////

						if((Right==1&&bkt_ptr->kmer_info.split_left==1)||(Right==0&&bkt_ptr->kmer_info.split_right==1))
						{
							flag=1;
							last_bkt_ptr->kmer_info.marked=1;
							//t_str=t_str.substr(0,t_str.size()-1);
							t_str.clear();
							sum_str+=t_str;
							t_str.clear();
							nBranches=1;
							break;
						}
	////////////////////////////////////////////////////////////////////////////////////


						if(bkt_ptr->kmer_info.used==1)
						{
							flag=1;
							break;
						}

						if(COVERAGE_STS)
						{
							coverage+=bkt_ptr->kmer_info.cov1;
							kmer_cnt++;
						}
						bkt_ptr->kmer_info.used=1;

						if((it==1&&Right==0)||(it==2&&Right==1))
						{
							bkt_ptr->kmer_info.flip=1;
						}
						else
						{
							bkt_ptr->kmer_info.flip=0;
						}


					}
					if(it==1)
					{
						s_right=sum_str;
						t_str.clear();
						sum_str.clear();
					}
					else
					{
						s_left=sum_str;
						reverse(s_left.begin(),s_left.end());
						complement_str(s_left);
						t_str.clear();
						sum_str.clear();
					}

					//left,rght search ended
					
					if(it==1)
					{
						nRBranches=nBranches;
					}
					if(it==2)
					{
						nLBranches=nBranches;
					}


				}
				bitsarr2str(beg_kmer.kmer, K_size, c_seq,2);
				contig=s_left+c_seq+s_right;

				if(ScreenOffTips)
				{
					if((int)(contig.size())<TipLenTh&&((double(coverage)/double(kmer_cnt))<TipCovTh))
					{
						num_contigs--;
						
						bkt_ptr=beg_bkt_ptr->nxt_bucket;
						continue;
					}
				}
				//contigs_vt.push_back(contig);
				o_contigs<<">Contig_"<<num_contigs;//<< " ";
				if(0)//COVERAGE_STS)
				{

					o_contigs<<"Avg_Cov: "<< (double(coverage)/double(kmer_cnt));
//					o_contigs<<" Kmer_cnt: "<<kmer_cnt;
					o_contigs<<" Contig_size: "<<contig.size();


				}
				o_contigsLen<<contig.size()<<endl;
				o_contigsCov<<(int)(double(coverage)/double(kmer_cnt)+0.5)<<endl;
				cov_hist[int((double(coverage)/double(kmer_cnt)+0.5))]++;


				o_contigs<<endl;
				o_contigs<<contig<<endl;
				o_nLB<<nLBranches<<endl;//" ";
				o_nRB<<nRBranches<<endl;//" ";
				
				o_contig_info<<">Contig "<<num_contigs<< " ";
				if(COVERAGE_STS)
				{

					o_contig_info<<"Avg_Cov: "<< (double(coverage)/double(kmer_cnt));
//					o_contig_info<<" Kmer_cnt: "<<kmer_cnt;
					o_contig_info<<" Contig_size: "<<contig.size()<<endl;

				}

			}



			bkt_ptr=beg_bkt_ptr->nxt_bucket;


		}

	}
	if(1)//Contig_Filename=="Contigs")
	{
	o_Log<<"Contigs count:"<<endl;
	o_Log<<num_contigs<<endl;
	}

	map<int,int>::iterator mit;
	for(mit=cov_hist.begin();mit!=cov_hist.end();++mit)
	{
		o_cov_hist<<mit->first<<" "<<mit->second<<endl;
	}

}


void build_contigs3(struct hashtable3 *ht,int K_size, int gap,string Contig_Filename,bool ScreenOffTips)
{
	//bool ScreenOffTips=0;
	int TipLenTh=K_size+1;
	int TipCovTh=1;
	ofstream o_contigsLen("Contigs_Len.txt");
	ofstream o_contigsCov("Contigs_Cov.txt");
	string o_name=Contig_Filename;
	ofstream o_contigs(o_name.c_str());
	//o_name=Contig_Filename+"_FE_left.txt";
	//ofstream o_FE_left(o_name.c_str());
	o_name="nLB_"+Contig_Filename;
	ofstream o_nLB(o_name.c_str());
	//o_name=Contig_Filename+"_FE_right.txt";
	//ofstream o_FE_right(o_name.c_str());
	o_name="nRB_"+Contig_Filename;
	ofstream o_nRB(o_name.c_str());
	o_name="Cov_Hist"+Contig_Filename;
	ofstream o_cov_hist(o_name.c_str());

	map<int,int> cov_hist;
	bool FE_left=0,FE_right=0;
	ofstream o_Log("ContigsLog.txt");//,ios_base::app);
	bool COVERAGE_STS=1;
	uint64_t coverage=0;
	uint64_t kmer_cnt=0,t,hv;
	size_t i=0,ht_sz=ht->ht_sz,hash_idx;
	size_t num_contigs=0;
	string s_left,s_right;
	s_left.reserve(100000);
	s_right.reserve(100000);
	struct bucket3 *bkt_ptr ,*beg_bkt_ptr,*last_bkt_ptr;
	kmer_t3 beg_kmer,t_kmer,f_kmer;
	o_name="Contigs_info.txt";
	ofstream o_contig_info(o_name.c_str());//"Contigs_info.txt");
	//char c_str[1000];
	bool found;
//	vector<string> contigs_vt;
	for(i=0;i<ht_sz;++i)
	{
		bkt_ptr=ht->store_pos[i];

		while(bkt_ptr!=NULL)
		{
			beg_bkt_ptr=bkt_ptr;
			if(bkt_ptr->kmer_info.used==0&&bkt_ptr->kmer_info.removed==0)
			{
				num_contigs++;
				FE_left=0;FE_right=0;

				if(COVERAGE_STS)
				{
					coverage=0;
					kmer_cnt=0;
					coverage+=bkt_ptr->kmer_info.cov1;
					kmer_cnt+=1;
				}

				string contig;
				beg_kmer=bkt_ptr->kmer_t3;
				bkt_ptr->kmer_info.used=1;

				//bkt_ptr->kmer_info.strand_visited|=1;
				bkt_ptr->kmer_info.flip=0;
				t_kmer=beg_kmer;


				char c_seq[1000];

				s_left.clear();
				s_right.clear();
				
				uint64_t left_bits,right_bits;
				int nLBranches=0,nRBranches=0;
				bool Right=0;
				

				for (int it=1;it<=2;++it)
				{
					//1:right search,2:left;
					bool flag=0;
					bool Free_End=0;
					bkt_ptr=beg_bkt_ptr;
					t_kmer=bkt_ptr->kmer_t3;
					string sum_str, t_str;
					
					int nBranches=0;
	
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
						sum_str+=t_str;
						t_str.clear();
						bkt_ptr->kmer_info.used=1;
						if(Right==1)
						{
							if(bkt_ptr->kmer_info.right==NULL)
							{
								Free_End=1;
								nBranches=0;
								break;
							}
						
							if(bkt_ptr->kmer_info.right->nxt_edge!=NULL)
							{
								edge_node *edge_ptr=bkt_ptr->kmer_info.right;
								nBranches=0;
								while(edge_ptr!=NULL)
								{
									nBranches++;
									edge_ptr=edge_ptr->nxt_edge;
								}
								break;
							}
							//no branch
							for(int j=(bkt_ptr->kmer_info.right)->len;j>=0;--j)
							{

								right_bits=(bkt_ptr->kmer_info.right->edge)>>2*j;

								switch(right_bits&0x3)
								{
									case 0:
										t=3;
										t<<=(K_size-1-64)*2;
										t_kmer.kmer[0]&=(~t);
										L_shift_NB(t_kmer.kmer,2,3);

										t_str.push_back('A');

										break;
									case 1:
										t=3;
										t<<=(K_size-1-64)*2;
										t_kmer.kmer[0]&=(~t);
										L_shift_NB(t_kmer.kmer,2,3);
										t=1;

										t_kmer.kmer[2]|=t;
										t_str.push_back('C');
										break;
									case 2:
										t=3;
										t<<=(K_size-1-64)*2;
										t_kmer.kmer[0]&=(~t);
										L_shift_NB(t_kmer.kmer,2,3);
										t=2;

										t_kmer.kmer[2]|=t;
										t_str.push_back('G');
										break;
									case 3:
										t=3;
										t<<=(K_size-1-64)*2;
										t_kmer.kmer[0]&=(~t);
										L_shift_NB(t_kmer.kmer,2,3);
										t=3;

										t_kmer.kmer[2]|=t;
										t_str.push_back('T');
										break;



								}
							}



							kmer_t3 tt_kmer=t_kmer;
							f_kmer=tt_kmer;
							get_rev_comp_seq_arr(f_kmer.kmer,K_size,3);
							if(uint64_t_cmp(t_kmer.kmer,f_kmer.kmer,3)>0)
							{
								tt_kmer=f_kmer;

							}

							hv=MurmurHash64A(&tt_kmer,sizeof(tt_kmer),0);

							hash_idx=(size_t) (hv%ht_sz);
							//2424710//1561973882672337490
							struct bucket3 ** ptr;

							ptr= &(ht->store_pos[hash_idx]);


							found=look_up_in_a_list3(&tt_kmer,&ptr);


						}
						else
						{
					//left search
							if(bkt_ptr->kmer_info.left==NULL)
							{
								Free_End=1;
								nBranches=0;
								break;
							}
							if(bkt_ptr->kmer_info.left->nxt_edge!=NULL)
							{
								edge_node *edge_ptr=bkt_ptr->kmer_info.left;
								nBranches=0;
								while(edge_ptr!=NULL)
								{
									nBranches++;
									edge_ptr=edge_ptr->nxt_edge;
								}
								break;//or left=...;
							}
							int edge_len=bkt_ptr->kmer_info.left->len;
							for(int j=0;j<=edge_len;++j)
							{

								left_bits=(bkt_ptr->kmer_info.left)->edge>>2*j;

								switch(left_bits&0x3)
								{
									case 0:
										R_shift_NB(t_kmer.kmer,2,3);

										t_str.push_back('T');

										break;
									case 1:
										R_shift_NB(t_kmer.kmer,2,3);

										t=((uint64_t)1)<<((K_size-1-64)*2);
										t_kmer.kmer[0]|=t;
										t_str.push_back('G');
										break;
									case 2:
										R_shift_NB(t_kmer.kmer,2,3);

										t=((uint64_t)2)<<((K_size-1-64)*2);
										t_kmer.kmer[0]|=t;
										t_str.push_back('C');
										break;
									case 3:
										R_shift_NB(t_kmer.kmer,2,3);

										t=((uint64_t)3)<<((K_size-1-64)*2);
										t_kmer.kmer[0]|=t;
										t_str.push_back('A');
										break;

								}
							}

						}


						f_kmer=t_kmer;
						get_rev_comp_seq_arr(f_kmer.kmer,K_size,3);
						if(uint64_t_cmp(t_kmer.kmer,f_kmer.kmer,3)>0)
						{
							t_kmer=f_kmer;
							Right=!Right;
						}


						hv=MurmurHash64A(&t_kmer,sizeof(t_kmer),0);

						hash_idx=(size_t) (hv%ht_sz);


						struct bucket3 ** ptr;

						ptr= &(ht->store_pos[hash_idx]);
						last_bkt_ptr=bkt_ptr;
						found=look_up_in_a_list3(&t_kmer,&ptr);

						bkt_ptr=*ptr;
					

						if(found==0)
						{
							Free_End=1;
							flag=1;

							break;
						}


	////////////////////////////////////////////////////////////////////////////////////

						if((Right==1&&bkt_ptr->kmer_info.split_left==1)||(Right==0&&bkt_ptr->kmer_info.split_right==1))
						{
							flag=1;
							last_bkt_ptr->kmer_info.marked=1;
							//t_str=t_str.substr(0,t_str.size()-1);
							t_str.clear();
							sum_str+=t_str;
							t_str.clear();
							nBranches=1;
							break;
						}
	////////////////////////////////////////////////////////////////////////////////////


						if(bkt_ptr->kmer_info.used==1)
						{
							flag=1;
							break;
						}

						if(COVERAGE_STS)
						{
							coverage+=bkt_ptr->kmer_info.cov1;
							kmer_cnt++;
						}
						bkt_ptr->kmer_info.used=1;

						if((it==1&&Right==0)||(it==2&&Right==1))
						{
							bkt_ptr->kmer_info.flip=1;
						}
						else
						{
							bkt_ptr->kmer_info.flip=0;
						}


					}
					if(it==1)
					{
						s_right=sum_str;
						t_str.clear();
						sum_str.clear();
					}
					else
					{
						s_left=sum_str;
						reverse(s_left.begin(),s_left.end());
						complement_str(s_left);
						t_str.clear();
						sum_str.clear();
					}

					//left,rght search ended
					
					if(it==1)
					{
						nRBranches=nBranches;
					}
					if(it==2)
					{
						nLBranches=nBranches;
					}


				}
				bitsarr2str(beg_kmer.kmer, K_size, c_seq,3);
				contig=s_left+c_seq+s_right;

				if(ScreenOffTips)
				{
					if((int)(contig.size())<TipLenTh&&((double(coverage)/double(kmer_cnt))<TipCovTh))
					{
						num_contigs--;
						
						bkt_ptr=beg_bkt_ptr->nxt_bucket;
						continue;
					}
				}
				//contigs_vt.push_back(contig);
				o_contigs<<">Contig_"<<num_contigs;//<< " ";
				if(0)//COVERAGE_STS)
				{

					o_contigs<<"Avg_Cov: "<< (double(coverage)/double(kmer_cnt));
//					o_contigs<<" Kmer_cnt: "<<kmer_cnt;
					o_contigs<<" Contig_size: "<<contig.size();


				}
				o_contigsLen<<contig.size()<<endl;
				o_contigsCov<<(int)(double(coverage)/double(kmer_cnt)+0.5)<<endl;
				cov_hist[int((double(coverage)/double(kmer_cnt)+0.5))]++;


				o_contigs<<endl;
				o_contigs<<contig<<endl;
				o_nLB<<nLBranches<<" ";
				o_nRB<<nRBranches<<" ";
				
				o_contig_info<<">Contig "<<num_contigs<< " ";
				if(COVERAGE_STS)
				{

					o_contig_info<<"Avg_Cov: "<< (double(coverage)/double(kmer_cnt));
//					o_contig_info<<" Kmer_cnt: "<<kmer_cnt;
					o_contig_info<<" Contig_size: "<<contig.size()<<endl;

				}

			}



			bkt_ptr=beg_bkt_ptr->nxt_bucket;


		}

	}
	if(1)//Contig_Filename=="Contigs")
	{
	o_Log<<"Contigs count:"<<endl;
	o_Log<<num_contigs<<endl;
	}

	map<int,int>::iterator mit;
	for(mit=cov_hist.begin();mit!=cov_hist.end();++mit)
	{
		o_cov_hist<<mit->first<<" "<<mit->second<<endl;
	}

}

void build_contigs4(struct hashtable4 *ht,int K_size, int gap,string Contig_Filename,bool ScreenOffTips)
{
	//bool ScreenOffTips=0;
	int TipLenTh=K_size+1;
	int TipCovTh=1;
	ofstream o_contigsLen("Contigs_Len.txt");
	ofstream o_contigsCov("Contigs_Cov.txt");
	string o_name=Contig_Filename;
	ofstream o_contigs(o_name.c_str());
	//o_name=Contig_Filename+"_FE_left.txt";
	//ofstream o_FE_left(o_name.c_str());
	o_name="nLB_"+Contig_Filename;
	ofstream o_nLB(o_name.c_str());
	//o_name=Contig_Filename+"_FE_right.txt";
	//ofstream o_FE_right(o_name.c_str());
	o_name="nRB_"+Contig_Filename;
	ofstream o_nRB(o_name.c_str());
	o_name="Cov_Hist"+Contig_Filename;
	ofstream o_cov_hist(o_name.c_str());

	map<int,int> cov_hist;
	bool FE_left=0,FE_right=0;
	ofstream o_Log("ContigsLog.txt");//,ios_base::app);
	bool COVERAGE_STS=1;
	uint64_t coverage=0;
	uint64_t kmer_cnt=0,t,hv;
	size_t i=0,ht_sz=ht->ht_sz,hash_idx;
	size_t num_contigs=0;
	string s_left,s_right;
	s_left.reserve(100000);
	s_right.reserve(100000);
	struct bucket4 *bkt_ptr ,*beg_bkt_ptr,*last_bkt_ptr;
	kmer_t4 beg_kmer,t_kmer,f_kmer;
	o_name="Contigs_info.txt";
	ofstream o_contig_info(o_name.c_str());//"Contigs_info.txt");
	//char c_str[1000];
	bool found;
//	vector<string> contigs_vt;
	for(i=0;i<ht_sz;++i)
	{
		bkt_ptr=ht->store_pos[i];

		while(bkt_ptr!=NULL)
		{
			beg_bkt_ptr=bkt_ptr;
			if(bkt_ptr->kmer_info.used==0&&bkt_ptr->kmer_info.removed==0)
			{
				num_contigs++;
				FE_left=0;FE_right=0;

				if(COVERAGE_STS)
				{
					coverage=0;
					kmer_cnt=0;
					coverage+=bkt_ptr->kmer_info.cov1;
					kmer_cnt+=1;
				}

				string contig;
				beg_kmer=bkt_ptr->kmer_t4;
				bkt_ptr->kmer_info.used=1;

				//bkt_ptr->kmer_info.strand_visited|=1;
				bkt_ptr->kmer_info.flip=0;
				t_kmer=beg_kmer;


				char c_seq[1000];

				s_left.clear();
				s_right.clear();
				
				uint64_t left_bits,right_bits;
				int nLBranches=0,nRBranches=0;
				bool Right=0;
				

				for (int it=1;it<=2;++it)
				{
					//1:right search,2:left;
					bool flag=0;
					bool Free_End=0;
					bkt_ptr=beg_bkt_ptr;
					t_kmer=bkt_ptr->kmer_t4;
					string sum_str, t_str;
					
					int nBranches=0;
	
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
						sum_str+=t_str;
						t_str.clear();
						bkt_ptr->kmer_info.used=1;
						if(Right==1)
						{
							if(bkt_ptr->kmer_info.right==NULL)
							{
								Free_End=1;
								nBranches=0;
								break;
							}
						
							if(bkt_ptr->kmer_info.right->nxt_edge!=NULL)
							{
								edge_node *edge_ptr=bkt_ptr->kmer_info.right;
								nBranches=0;
								while(edge_ptr!=NULL)
								{
									nBranches++;
									edge_ptr=edge_ptr->nxt_edge;
								}
								break;
							}
							//no branch
							for(int j=(bkt_ptr->kmer_info.right)->len;j>=0;--j)
							{

								right_bits=(bkt_ptr->kmer_info.right->edge)>>2*j;

								switch(right_bits&0x3)
								{
									case 0:
										t=3;
										t<<=(K_size-1-96)*2;
										t_kmer.kmer[0]&=(~t);
										L_shift_NB(t_kmer.kmer,2,4);

										t_str.push_back('A');

										break;
									case 1:
										t=3;
										t<<=(K_size-1-96)*2;
										t_kmer.kmer[0]&=(~t);
										L_shift_NB(t_kmer.kmer,2,4);
										t=1;

										t_kmer.kmer[3]|=t;
										t_str.push_back('C');
										break;
									case 2:
										t=3;
										t<<=(K_size-1-96)*2;
										t_kmer.kmer[0]&=(~t);
										L_shift_NB(t_kmer.kmer,2,4);
										t=2;

										t_kmer.kmer[3]|=t;
										t_str.push_back('G');
										break;
									case 3:
										t=3;
										t<<=(K_size-1-96)*2;
										t_kmer.kmer[0]&=(~t);
										L_shift_NB(t_kmer.kmer,2,4);
										t=3;

										t_kmer.kmer[3]|=t;
										t_str.push_back('T');
										break;



								}
							}



							kmer_t4 tt_kmer=t_kmer;
							f_kmer=tt_kmer;
							get_rev_comp_seq_arr(f_kmer.kmer,K_size,4);
							if(uint64_t_cmp(t_kmer.kmer,f_kmer.kmer,4)>0)
							{
								tt_kmer=f_kmer;

							}

							hv=MurmurHash64A(&tt_kmer,sizeof(tt_kmer),0);

							hash_idx=(size_t) (hv%ht_sz);
							//2424710//1561973882672337490
							struct bucket4 ** ptr;

							ptr= &(ht->store_pos[hash_idx]);


							found=look_up_in_a_list4(&tt_kmer,&ptr);


						}
						else
						{
					//left search
							if(bkt_ptr->kmer_info.left==NULL)
							{
								Free_End=1;
								nBranches=0;
								break;
							}
							if(bkt_ptr->kmer_info.left->nxt_edge!=NULL)
							{
								edge_node *edge_ptr=bkt_ptr->kmer_info.left;
								nBranches=0;
								while(edge_ptr!=NULL)
								{
									nBranches++;
									edge_ptr=edge_ptr->nxt_edge;
								}
								break;//or left=...;
							}
							int edge_len=bkt_ptr->kmer_info.left->len;
							for(int j=0;j<=edge_len;++j)
							{

								left_bits=(bkt_ptr->kmer_info.left)->edge>>2*j;

								switch(left_bits&0x3)
								{
									case 0:
										R_shift_NB(t_kmer.kmer,2,4);

										t_str.push_back('T');

										break;
									case 1:
										R_shift_NB(t_kmer.kmer,2,4);

										t=((uint64_t)1)<<((K_size-1-96)*2);
										t_kmer.kmer[0]|=t;
										t_str.push_back('G');
										break;
									case 2:
										R_shift_NB(t_kmer.kmer,2,4);

										t=((uint64_t)2)<<((K_size-1-96)*2);
										t_kmer.kmer[0]|=t;
										t_str.push_back('C');
										break;
									case 3:
										R_shift_NB(t_kmer.kmer,2,4);

										t=((uint64_t)3)<<((K_size-1-96)*2);
										t_kmer.kmer[0]|=t;
										t_str.push_back('A');
										break;

								}
							}

						}


						f_kmer=t_kmer;
						get_rev_comp_seq_arr(f_kmer.kmer,K_size,4);
						if(uint64_t_cmp(t_kmer.kmer,f_kmer.kmer,4)>0)
						{
							t_kmer=f_kmer;
							Right=!Right;
						}


						hv=MurmurHash64A(&t_kmer,sizeof(t_kmer),0);

						hash_idx=(size_t) (hv%ht_sz);


						struct bucket4 ** ptr;

						ptr= &(ht->store_pos[hash_idx]);
						last_bkt_ptr=bkt_ptr;
						found=look_up_in_a_list4(&t_kmer,&ptr);

						bkt_ptr=*ptr;
					

						if(found==0)
						{
							Free_End=1;
							flag=1;

							break;
						}


	////////////////////////////////////////////////////////////////////////////////////

						if((Right==1&&bkt_ptr->kmer_info.split_left==1)||(Right==0&&bkt_ptr->kmer_info.split_right==1))
						{
							flag=1;
							last_bkt_ptr->kmer_info.marked=1;
							//t_str=t_str.substr(0,t_str.size()-1);
							t_str.clear(); 
							sum_str += t_str;
							t_str.clear();
							nBranches=1;
							break;
						}
	////////////////////////////////////////////////////////////////////////////////////


						if(bkt_ptr->kmer_info.used==1)
						{
							flag=1;
							break;
						}

						if(COVERAGE_STS)
						{
							coverage+=bkt_ptr->kmer_info.cov1;
							kmer_cnt++;
						}
						bkt_ptr->kmer_info.used=1;

						if((it==1&&Right==0)||(it==2&&Right==1))
						{
							bkt_ptr->kmer_info.flip=1;
						}
						else
						{
							bkt_ptr->kmer_info.flip=0;
						}


					}
					if(it==1)
					{
						s_right=sum_str;
						t_str.clear();
						sum_str.clear();
					}
					else
					{
						s_left=sum_str;
						reverse(s_left.begin(),s_left.end());
						complement_str(s_left);
						t_str.clear();
						sum_str.clear();
					}

					//left,rght search ended
					
					if(it==1)
					{
						nRBranches=nBranches;
					}
					if(it==2)
					{
						nLBranches=nBranches;
					}


				}
				bitsarr2str(beg_kmer.kmer, K_size, c_seq,4);
				contig=s_left+c_seq+s_right;

				if(ScreenOffTips)
				{
					if((int)(contig.size())<TipLenTh&&((double(coverage)/double(kmer_cnt))<TipCovTh))
					{
						num_contigs--;
						
						bkt_ptr=beg_bkt_ptr->nxt_bucket;
						continue;
					}
				}
				//contigs_vt.push_back(contig);
				o_contigs<<">Contig_"<<num_contigs;//<< " ";
				if(0)//COVERAGE_STS)
				{

					o_contigs<<"Avg_Cov: "<< (double(coverage)/double(kmer_cnt));
//					o_contigs<<" Kmer_cnt: "<<kmer_cnt;
					o_contigs<<" Contig_size: "<<contig.size();


				}
				o_contigsLen<<contig.size()<<endl;
				o_contigsCov<<(int)(double(coverage)/double(kmer_cnt)+0.5)<<endl;
				cov_hist[int((double(coverage)/double(kmer_cnt)+0.5))]++;


				o_contigs<<endl;
				o_contigs<<contig<<endl;
				o_nLB<<nLBranches<<" ";
				o_nRB<<nRBranches<<" ";
				
				o_contig_info<<">Contig "<<num_contigs<< " ";
				if(COVERAGE_STS)
				{

					o_contig_info<<"Avg_Cov: "<< (double(coverage)/double(kmer_cnt));
//					o_contig_info<<" Kmer_cnt: "<<kmer_cnt;
					o_contig_info<<" Contig_size: "<<contig.size()<<endl;

				}

			}



			bkt_ptr=beg_bkt_ptr->nxt_bucket;


		}

	}
	if(1)//Contig_Filename=="Contigs")
	{
	o_Log<<"Contigs count:"<<endl;
	o_Log<<num_contigs<<endl;
	}

	map<int,int>::iterator mit;
	for(mit=cov_hist.begin();mit!=cov_hist.end();++mit)
	{
		o_cov_hist<<mit->first<<" "<<mit->second<<endl;
	}

}

void build_contigs0(struct hashtable0 *ht, key_table *key_table, int K_size, int gap,string Contig_Filename,bool ScreenOffTips)
{
	//bool ScreenOffTips=1;
	int TipLenTh=100;
	int TipCovTh=2;
	char ACGT[5]="ACGT";
	string o_name=Contig_Filename;
	ofstream o_contigs(o_name.c_str());
	ofstream o_contigsLen("Contigs_Len.txt");
	ofstream o_contigsCov("Contigs_Cov.txt");
	o_name="nLB_"+Contig_Filename;
	ofstream o_nLB(o_name.c_str());
	o_name="nRB_"+Contig_Filename;
	ofstream o_nRB(o_name.c_str());
	o_name="Cov_Hist"+Contig_Filename;
	ofstream o_cov_hist(o_name.c_str());

	map<int,int> cov_hist;
	bool FE_left=0,FE_right=0;
	ofstream o_Log("ContigsLog.txt");//,ios_base::app);
	bool COVERAGE_STS=1;
	uint64_t coverage=0;
	uint64_t kmer_cnt=0;
	size_t i=0,ht_sz=ht->ht_sz;
	size_t num_contigs=0;
	string s_left,s_right;
	s_left.reserve(100000);
	s_right.reserve(100000);


	int Kmer_arr_sz=K_size/32+1;
	int rem1=K_size%32;
	if(rem1==0)
	{Kmer_arr_sz--;}



	struct bucket0 *bkt_ptr ,*beg_bkt_ptr,*last_bkt_ptr;
	uint64_t beg_kmer[100],t_kmer[100],t,f_kmer[100],tt_kmer[100],hv,hash_idx;
	o_name="Contigs_info.txt";
	ofstream o_contig_info(o_name.c_str());//"Contigs_info.txt");
	//char c_str[1000];
	bool found;
//	vector<string> contigs_vt;
	for(i=0;i<ht_sz;++i)
	{
		bkt_ptr=ht->store_pos[i];

		while(bkt_ptr!=NULL)
		{
			beg_bkt_ptr=bkt_ptr;
			if(bkt_ptr->kmer_info.used==0&&bkt_ptr->kmer_info.removed==0)
			{
				num_contigs++;
				
				FE_left=0;FE_right=0;

				if(COVERAGE_STS)
				{
					coverage=0;
					kmer_cnt=0;
					coverage+=bkt_ptr->kmer_info.cov1;
					kmer_cnt+=1;
				}

				string contig;

				memcpy(beg_kmer,bkt_ptr->kmer_t,sizeof(uint64_t)*Kmer_arr_sz);
				bkt_ptr->kmer_info.used=1;

				//bkt_ptr->kmer_info.strand_visited|=1;
				bkt_ptr->kmer_info.flip=0;
				

				memcpy(t_kmer,beg_kmer,sizeof(uint64_t)*Kmer_arr_sz);

				char c_seq[1000];

				s_left.clear();
				s_right.clear();
				
				int nLBranches=0,nRBranches=0;
				uint64_t left_bits,right_bits;

				bool Right=0,flag=0;
				
				
				for(int it=1;it<=2;++it)
				{
					string sum_str,t_str;
					
					bool Free_End=0;
					int nBranches=0;
					sum_str.clear();
					t_str.clear();
					//right search
				
					if(it==1)
					{Right=1;}
					else
					{
						Right=0;
					}
					bkt_ptr=beg_bkt_ptr;
					

					memcpy(t_kmer,bkt_ptr->kmer_t,sizeof(uint64_t)*Kmer_arr_sz);


					while(1)
					{

						sum_str+=t_str;
						t_str.clear();
						bkt_ptr->kmer_info.used=1;
						if(Right==1)
						{
						
							if(bkt_ptr->kmer_info.right==NULL)
							{
								Free_End=1;
								nBranches=0;
								break;
							}
							if(bkt_ptr->kmer_info.right->nxt_edge!=NULL)
							{
								nBranches=0;
								edge_node *edge_ptr=bkt_ptr->kmer_info.right;
								while(edge_ptr!=NULL)
								{
									nBranches++;
									edge_ptr=edge_ptr->nxt_edge;
								}
							
								break;
							}
							int edge_len=(bkt_ptr->kmer_info.right)->len;
							for(int j=edge_len;j>=0;--j)
							{

								right_bits=(bkt_ptr->kmer_info.right->edge);
								right_bits>>=(2*j);
								uint64_t b=right_bits&0x3;

								//set first bit 0;
								t=3;
								t<<=(((K_size-1)%32)*2);
								t_kmer[0]&=(~t);
								L_shift_NB(t_kmer,2,Kmer_arr_sz);
								//do the rest work
								t_kmer[Kmer_arr_sz-1]|=b;
								t_str.push_back(ACGT[b]);
								
								
							}

						}
						else
						{
		
							if(bkt_ptr->kmer_info.left==NULL)
							{
								Free_End=1;
								nBranches=0;
								break;
							}
							if(bkt_ptr->kmer_info.left->nxt_edge!=NULL)
							{
								nBranches=0;
								edge_node *edge_ptr=bkt_ptr->kmer_info.left;
								while(edge_ptr!=NULL)
								{
									nBranches++;
									edge_ptr=edge_ptr->nxt_edge;
								}
								break;
							}
							int edge_len=bkt_ptr->kmer_info.left->len;
							for(int j=0;j<=edge_len;++j)
							{

								left_bits=(bkt_ptr->kmer_info.left)->edge;
								left_bits>>=(2*j);

								uint64_t b=(left_bits&0x3);
								R_shift_NB(t_kmer,2,Kmer_arr_sz);
								t=b<<((K_size-1)%32*2);
								t_kmer[0]|=t;
								t_str.push_back(ACGT[3-b]);


							}


						}

					

						memcpy(f_kmer,t_kmer,sizeof(uint64_t)*Kmer_arr_sz);
						get_rev_comp_seq_arr(f_kmer,K_size,Kmer_arr_sz);

						if(uint64_t_cmp(t_kmer,f_kmer,Kmer_arr_sz)>0)
						{
							
							memcpy(t_kmer,f_kmer,sizeof(uint64_t)*Kmer_arr_sz);
							Right=!Right;
						}

						

						hv=MurmurHash64A(t_kmer,sizeof(uint64_t)*Kmer_arr_sz,0);

						hash_idx=(size_t) (hv%ht_sz);


						struct bucket0 ** ptr;

						ptr= &(ht->store_pos[hash_idx]);
						last_bkt_ptr=bkt_ptr;
						found=look_up_in_a_list0(t_kmer,&ptr,Kmer_arr_sz);

						bkt_ptr=*ptr;
					
						if(found==0)
						{
							Free_End=1;
							
							flag=1;

							break;
						}


	////////////////////////////////////////////////////////////////////////////////////

						if((Right==1&&bkt_ptr->kmer_info.split_left==1)||(Right==0&&bkt_ptr->kmer_info.split_right==1))
						{
							flag=1;
							last_bkt_ptr->kmer_info.marked=1;
							//t_str=t_str.substr(0,t_str.size()-1);
							t_str.clear();
							sum_str+=t_str;
							t_str.clear();
							nBranches=1;

							break;
						}
	////////////////////////////////////////////////////////////////////////////////////


						if(bkt_ptr->kmer_info.used==1)
						{
							flag=1;
							break;
						}

						if(COVERAGE_STS)
						{
							coverage+=bkt_ptr->kmer_info.cov1;
							kmer_cnt++;

						}
					
						if((Right==0&&it==1)||(Right==1&&it==2))
						{
							bkt_ptr->kmer_info.flip=1;
						}
						else
						{
							bkt_ptr->kmer_info.flip=0;
						}


					}
					if(it==1)
					{
						nRBranches=nBranches;
						s_right=sum_str;
						t_str.clear();
						sum_str.clear();
					}
					if(it==2)
					{
						nLBranches=nBranches;
						s_left=sum_str;
						reverse(s_left.begin(),s_left.end());
						complement_str(s_left);
						t_str.clear();
						sum_str.clear();
						
					}

					
				}
				//left,rght search ended

				bitsarr2str(beg_kmer, K_size, c_seq,Kmer_arr_sz);
			
				contig=s_left+c_seq+s_right;
				
				if(ScreenOffTips)
				{
					if((int)(contig.size())<TipLenTh&&((double(coverage)/double(kmer_cnt))<TipCovTh))
					{
						num_contigs--;
						bkt_ptr=beg_bkt_ptr->nxt_bucket;
						continue;
					}
				}


				//contigs_vt.push_back(contig);
				o_contigs<<">Contig_"<<num_contigs;//<< " ";
				if(0)//COVERAGE_STS)
				{
					o_contigs<<"Avg_Cov: "<< (double(coverage)/double(kmer_cnt));
//					o_contigs<<" Kmer_cnt: "<<kmer_cnt;
					o_contigs<<" Contig_size: "<<contig.size();
				}
				o_contigsLen<<contig.size()<<endl;
				o_contigsCov<<(int)(double(coverage)/double(kmer_cnt)+0.5)<<endl;
				cov_hist[int((double(coverage)/double(kmer_cnt)+0.5))]++;


				o_contigs<<endl;
				o_contigs<<contig<<endl;
			
				o_nLB<<nLBranches<<endl;
				o_nRB<<nRBranches<<endl;

				o_contig_info<<">Contig "<<num_contigs<< " ";
				if(COVERAGE_STS)
				{

					o_contig_info<<"Avg_Cov: "<< (double(coverage)/double(kmer_cnt));
					o_contig_info<<" Kmer_cnt: "<<kmer_cnt;
					o_contig_info<<" Contig_size: "<<contig.size()<<endl;

				}
			}
			bkt_ptr=beg_bkt_ptr->nxt_bucket;
		}

	}
	if(1)
	{
		o_Log<<"Contigs count:"<<endl;
		o_Log<<num_contigs<<endl;
	}

	map<int,int>::iterator mit;
	for(mit=cov_hist.begin();mit!=cov_hist.end();++mit)
	{
		o_cov_hist<<mit->first<<" "<<mit->second<<endl;
	}

}
















#endif

