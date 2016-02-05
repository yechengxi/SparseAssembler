#ifndef __SCAFFOLDING_DATA_STRUCTURE_H
#define __SCAFFOLDING_DATA_STRUCTURE_H

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


using namespace std;


//The below is for scaffolding, can be omitted first %{

struct contig_t
{
	uint64_t contig_bits[100000];
	int contigLen;
};


struct adjacent_contig_info
{
	int32_t dist_sum;
	int8_t cov;
	string bridge;

};


struct scaffold_contig_info
{
	int32_t dist_sum;
	int8_t cov;

};

struct c_info
{
	 
	uint8_t used:1,removed:1,unique:1,flip:1,marked:1,loc_unique:1,positive:1,rcomp:1;
};


//%} The above is for scaffolding, can be omitted first 



bool it_cmp( const vector<int>::iterator &a,const vector<int>::iterator &b)
{
	return (*a) > (*b);
}

void Init_Contig(string &seq,struct contig_t & contig)
{
	contig.contigLen=(int)seq.size();
	int Contig_arr_sz=contig.contigLen/32+1;
	int rem=contig.contigLen%32;
	if(rem==0)
	{Contig_arr_sz--;}

	str2bitsarr(seq.c_str(),contig.contigLen,contig.contig_bits,Contig_arr_sz);

}



void ContigsRemapping(struct hashtable *ht,struct hashtable2 *ht2, int K_size,struct contigs_info * contigs_info,string Contig_Filename,bool RecordKmerCnt)
{
	time_t beg_time,read_time;
	//int64_t scaffold_len;
	//int Kmer_arr_sz;
////////////////////////////////
	contigs_info->K_size=K_size;
	string in_fname=Contig_Filename;
	ifstream Contigs_in(in_fname.c_str());
	string ctg_cov_fname="Contigs_Cov.txt",ctg_hp_fname="Contigs_HP.txt";

	if(Contig_Filename=="SuperContigs.txt")
	{
		ctg_cov_fname="SuperContigs_Cov.txt";
		ctg_hp_fname="SuperContigs_HP.txt";
	}
	ifstream in_ContigsCov(ctg_cov_fname.c_str());
	ofstream o_hit_pos(ctg_hp_fname.c_str());
	size_t num_Contigs=0;
	string tag,s,kmer;
	uint64_t seq,f_seq;
	kmer_t2 seq_t2,f_seq_t2;
	uint64_t hv;
	size_t hash_idx;
	bool flip_c;
	size_t nKmer=0;
	size_t tot_Kmers=0;
	size_t boundary_kmers=0;
	size_t non_repeat_kmers=0,repeat_kmers=0;
	//char c_str[1000];
	//vector<int> contig_sz;
	//bool flip_1,flip_2,flip_0;

	size_t ht_sz;
	if(K_size<=32)
	{
	ht_sz=ht->ht_sz;
	}
	else
	{
		ht_sz=ht2->ht_sz;
	}
	if(K_size<=32)
	{
		for(size_t i=0;i<ht->ht_sz;++i)
		{
			size_t list_sz=0;

			struct bucket *bktptr=ht->store_pos[i];
			while(bktptr!=NULL)
			{

				bktptr->kmer_info.contig_no=0;
				bktptr->kmer_info.cod=0;
				bktptr->kmer_info.repeat=0;
				bktptr=bktptr->nxt_bucket;
				list_sz++;
			}

		}
	}
	else
	{
		if(K_size>32&&K_size<=64)
		{
			for(size_t i=0;i<ht2->ht_sz;++i)
			{
				size_t list_sz=0;
				struct bucket2 *bktptr=ht2->store_pos[i];
				while(bktptr!=NULL)
				{

					bktptr->kmer_info.contig_no=0;
					bktptr->kmer_info.cod=0;
					bktptr->kmer_info.repeat=0;
					bktptr=bktptr->nxt_bucket;
					list_sz++;
				}
			}
		}

	}

	bool found;
	cout<<"Marking kmers..."<<endl;
	time(&beg_time);

	contigs_info->contig_sz_vt.clear();
	contigs_info->contig_sz_vt.push_back(0);
	contigs_info->cov_vt.clear();
	contigs_info->cov_vt.push_back(0);
	
	if(RecordKmerCnt)
	{
		contigs_info->kmer_cnt_vt.clear();
		contigs_info->kmer_cnt_vt.push_back(0);
	}
	string s1,s2,s3,cont_s;
	int ContLen;
	//double AvgCov ;

	while(getline(Contigs_in,s1))//(Contigs_in>>s1>>ContigNo>>s2>>AvgCov>>s3>>ContLen)		//getline(inContigs,tag))
	{
		getline(Contigs_in,cont_s);
		int ctg_cov=0;
		in_ContigsCov>>ctg_cov;
		contigs_info->cov_vt.push_back(ctg_cov);
		if(cont_s[cont_s.size()-1]=='\r'||cont_s[cont_s.size()-1]=='\n')
		{cont_s.resize(cont_s.size()-1);}
		if(cont_s.size()==0)
		{
			getline(Contigs_in,cont_s);
			if(cont_s[cont_s.size()-1]=='\r'||cont_s[cont_s.size()-1]=='\n')
			{cont_s.resize(cont_s.size()-1);}
		}
		nKmer=0;
		num_Contigs++;

		ContLen=cont_s.size();
		contigs_info->contig_sz_vt.push_back(ContLen);
		//contigs_info->cov_vt.push_back(int(AvgCov+0.5));
		struct contig_t Contig;
		Init_Contig(cont_s,Contig);
	
		bool Begin=0;
		int hp=0;
		for(int i=0;i<Contig.contigLen-K_size+1;++i)
		{


			if(K_size<=32)
			{
				get_sub_arr(Contig.contig_bits,Contig.contigLen,i,K_size,&seq);
				f_seq=get_rev_comp_seq(seq,K_size);

				
				flip_c=0;
				if(seq>f_seq)
				{
					flip_c=1;
					seq=f_seq;
				}

				hv=MurmurHash64A(&seq,sizeof(seq),0);
				hash_idx=(size_t) (hv%ht_sz);
				struct bucket **ptr;
				ptr= &(ht->store_pos[hash_idx]);
				found=look_up_in_a_list(seq,&ptr);

				if(found)
				{

					hp=i;
					if(Begin==0)
					{
						o_hit_pos<<hp<<" ";
						Begin=1;
					}
					
					nKmer++;
					tot_Kmers++;
					if(RecordKmerCnt)
					{
						(*ptr)->kmer_info.cod=nKmer;
					}
					else
					{
						(*ptr)->kmer_info.cod=(int )i;
					}
					if((*ptr)->kmer_info.repeat==0&&(*ptr)->kmer_info.contig_no==0)
					{
						(*ptr)->kmer_info.contig_no=(int)num_Contigs;
						non_repeat_kmers++;
					}
					else
					{
						(*ptr)->kmer_info.repeat=1;
						(*ptr)->kmer_info.contig_no=0;
						repeat_kmers++;
					}
					(*ptr)->kmer_info.flip=flip_c;
					if(i<100||i>(Contig.contigLen-K_size-100))
					{
						(*ptr)->kmer_info.masked=1;
					}
					else
					{
						(*ptr)->kmer_info.masked=0;
					}
				}

			}
			else
			{
				get_sub_arr(Contig.contig_bits,Contig.contigLen,i,K_size,seq_t2.kmer);
				f_seq_t2=seq_t2;
				get_rev_comp_seq_arr(f_seq_t2.kmer,K_size,2);
				flip_c=0;
				if(uint64_t_cmp(seq_t2.kmer,f_seq_t2.kmer,2) > 0)
				{
					flip_c=1;
					seq_t2=f_seq_t2;
				}

				hv=MurmurHash64A(seq_t2.kmer,sizeof(seq_t2),0);
				hash_idx=(size_t) (hv%ht_sz);
				struct bucket2 **ptr;
				ptr= &(ht2->store_pos[hash_idx]);
				found=look_up_in_a_list2(&seq_t2,&ptr);

				if(found)
				{
					hp=i;
					if(Begin==0)
					{
						o_hit_pos<<hp<<" ";
						Begin=1;
					}

					nKmer++;
					tot_Kmers++;
					if(RecordKmerCnt)
					{
						(*ptr)->kmer_info.cod=nKmer;
					}
					else
					{
						(*ptr)->kmer_info.cod=(int )i;
					}

					if(((*ptr)->kmer_info.repeat==0)&&(*ptr)->kmer_info.contig_no==0)
					{
						(*ptr)->kmer_info.contig_no=(int)num_Contigs;
					}
					else
					{
						(*ptr)->kmer_info.repeat=1;
						(*ptr)->kmer_info.contig_no=0;
					}

					(*ptr)->kmer_info.flip=flip_c;
					if(i<100||i>(Contig.contigLen-K_size-100))
					{
						(*ptr)->kmer_info.masked=1;
					}
					else
					{
						(*ptr)->kmer_info.masked=0;
					}

				}
			}

			if(i==Contig.contigLen-K_size)
			o_hit_pos<<hp<<endl;;

		}

		if(RecordKmerCnt)
		{
			contigs_info->kmer_cnt_vt.push_back(nKmer);
		}

	}

	//cout<<"tot_Kmers "<<tot_Kmers<<endl;
	//cout<<"repeat_kmers "<<repeat_kmers<<endl;
	//cout<<"non_repeat_kmers "<<non_repeat_kmers<<endl;
	vector<int>::iterator vit=contigs_info->contig_sz_vt.begin();
	contigs_info->LengthRank.clear();
	num_Contigs=contigs_info->contig_sz_vt.size()-1;
	vit++;
	for(size_t i=0;i<num_Contigs;++i)
	{
		contigs_info->LengthRank.push_back(vit);
		vit++;
	}

	contigs_info->total_contigs=num_Contigs;
	sort(contigs_info->LengthRank.begin(),contigs_info->LengthRank.end(),it_cmp);

	time(&read_time);
	cout<<"Kmer marking time: "<<difftime(read_time,beg_time)<<" secs."<<endl;

}


void ContigsRemapping3(struct hashtable3 *ht,int K_size,struct contigs_info * contigs_info,string Contig_Filename,bool RecordKmerCnt)
{
	time_t beg_time,read_time;
	contigs_info->K_size=K_size;
	string in_fname=Contig_Filename;
	ifstream Contigs_in(in_fname.c_str());	string ctg_cov_fname="Contigs_Cov.txt",ctg_hp_fname="Contigs_HP.txt";

	if(Contig_Filename=="SuperContigs.txt")
	{
		ctg_cov_fname="SuperContigs_Cov.txt";
		ctg_hp_fname="SuperContigs_HP.txt";
	}
	ifstream in_ContigsCov(ctg_cov_fname.c_str());
	ofstream o_hit_pos(ctg_hp_fname.c_str());
	size_t num_Contigs=0;
	string tag,s,kmer;
	kmer_t3 seq_t3,f_seq_t3;
	uint64_t hv;
	size_t hash_idx;
	bool flip_c;
	size_t nKmer=0;
	
	size_t ht_sz;
	
	ht_sz=ht->ht_sz;
	
	for(size_t i=0;i<ht->ht_sz;++i)
	{
		size_t list_sz=0;
		struct bucket3 *bktptr=ht->store_pos[i];
		while(bktptr!=NULL)
		{
			bktptr->kmer_info.contig_no=0;
			bktptr->kmer_info.cod=0;
			bktptr->kmer_info.repeat=0;
			bktptr=bktptr->nxt_bucket;
			list_sz++;
		}
	}
		


	bool found;
	cout<<"Marking kmers..."<<endl;
	time(&beg_time);

	contigs_info->contig_sz_vt.clear();
	contigs_info->contig_sz_vt.push_back(0);
	contigs_info->cov_vt.clear();
	contigs_info->cov_vt.push_back(0);
	
	if(RecordKmerCnt)
	{
		contigs_info->kmer_cnt_vt.clear();
		contigs_info->kmer_cnt_vt.push_back(0);
	}
	string s1,s2,s3,cont_s;
	int ContLen;
	//double AvgCov ;

	while(getline(Contigs_in,s1))//(Contigs_in>>s1>>ContigNo>>s2>>AvgCov>>s3>>ContLen)		//getline(inContigs,tag))
	{
		getline(Contigs_in,cont_s);
		int ctg_cov=0;
		in_ContigsCov>>ctg_cov;
		contigs_info->cov_vt.push_back(ctg_cov);
		if(cont_s[cont_s.size()-1]=='\r'||cont_s[cont_s.size()-1]=='\n')
		{cont_s.resize(cont_s.size()-1);}
		if(cont_s.size()==0)
		{
			getline(Contigs_in,cont_s);
			if(cont_s[cont_s.size()-1]=='\r'||cont_s[cont_s.size()-1]=='\n')
			{cont_s.resize(cont_s.size()-1);}
		}
		nKmer=0;
		num_Contigs++;

		ContLen=cont_s.size();
		contigs_info->contig_sz_vt.push_back(ContLen);
		struct contig_t Contig;
		Init_Contig(cont_s,Contig);
		
		bool Begin=0;
		int hp=0;
		for(int i=0;i<Contig.contigLen-K_size+1;++i)
		{
			get_sub_arr(Contig.contig_bits,Contig.contigLen,i,K_size,seq_t3.kmer);
			f_seq_t3=seq_t3;
			get_rev_comp_seq_arr(f_seq_t3.kmer,K_size,3);
			flip_c=0;
			if(uint64_t_cmp(seq_t3.kmer,f_seq_t3.kmer,3) > 0)
			{
				flip_c=1;
				seq_t3=f_seq_t3;
			}

			hv=MurmurHash64A(seq_t3.kmer,sizeof(seq_t3),0);
			hash_idx=(size_t) (hv%ht_sz);
			struct bucket3 **ptr;
			ptr= &(ht->store_pos[hash_idx]);
			found=look_up_in_a_list3(&seq_t3,&ptr);

			if(found)
			{
				hp=i;
				if(Begin==0)
				{
					o_hit_pos<<hp<<" ";
					Begin=1;
				}

				nKmer++;
				if(RecordKmerCnt)
				{
					(*ptr)->kmer_info.cod=nKmer;
				}
				else
				{
					(*ptr)->kmer_info.cod=(int )i;
				}

				if(((*ptr)->kmer_info.repeat==0)&&(*ptr)->kmer_info.contig_no==0)
				{
					(*ptr)->kmer_info.contig_no=(int)num_Contigs;
				}
				else
				{
					(*ptr)->kmer_info.repeat=1;
					(*ptr)->kmer_info.contig_no=0;
				}

				(*ptr)->kmer_info.flip=flip_c;
			}
		

			if(i==Contig.contigLen-K_size)
			o_hit_pos<<hp<<endl;

		}

		if(RecordKmerCnt)
		{
			contigs_info->kmer_cnt_vt.push_back(nKmer);
		}

	}

	vector<int>::iterator vit=contigs_info->contig_sz_vt.begin();
	num_Contigs=contigs_info->contig_sz_vt.size()-1;
	contigs_info->LengthRank.clear();
	vit++;
	for(size_t i=0;i<num_Contigs;++i)
	{
		contigs_info->LengthRank.push_back(vit);
		vit++;
	}

	contigs_info->total_contigs=num_Contigs;
	sort(contigs_info->LengthRank.begin(),contigs_info->LengthRank.end(),it_cmp);

	time(&read_time);
	cout<<"Kmer marking time: "<<difftime(read_time,beg_time)<<" secs."<<endl;

}

void ContigsRemapping4(struct hashtable4 *ht,int K_size,struct contigs_info * contigs_info,string Contig_Filename,bool RecordKmerCnt)
{
	time_t beg_time,read_time;
	contigs_info->K_size=K_size;
	string in_fname=Contig_Filename;
	ifstream Contigs_in(in_fname.c_str());

	string ctg_cov_fname="Contigs_Cov.txt",ctg_hp_fname="Contigs_HP.txt";

	if(Contig_Filename=="SuperContigs.txt")
	{
		ctg_cov_fname="SuperContigs_Cov.txt";
		ctg_hp_fname="SuperContigs_HP.txt";
	}
	ifstream in_ContigsCov(ctg_cov_fname.c_str());
	ofstream o_hit_pos(ctg_hp_fname.c_str());

	size_t num_Contigs=0;
	string tag,s,kmer;
	kmer_t4 seq_t4,f_seq_t4;
	uint64_t hv;
	size_t hash_idx;
	bool flip_c;
	size_t nKmer=0;
	
	size_t ht_sz;
	
	ht_sz=ht->ht_sz;
	
	for(size_t i=0;i<ht->ht_sz;++i)
	{
		size_t list_sz=0;
		struct bucket4 *bktptr=ht->store_pos[i];
		while(bktptr!=NULL)
		{
			bktptr->kmer_info.contig_no=0;
			bktptr->kmer_info.cod=0;
			bktptr->kmer_info.repeat=0;
			bktptr=bktptr->nxt_bucket;
			list_sz++;
		}
	}
		


	bool found;
	cout<<"Marking kmers..."<<endl;
	time(&beg_time);

	contigs_info->contig_sz_vt.clear();
	contigs_info->contig_sz_vt.push_back(0);
	contigs_info->cov_vt.clear();
	contigs_info->cov_vt.push_back(0);
	
	if(RecordKmerCnt)
	{
		contigs_info->kmer_cnt_vt.clear();
		contigs_info->kmer_cnt_vt.push_back(0);
	}
	string s1,s2,s3,cont_s;
	int ContLen;
	//double AvgCov ;

	while(getline(Contigs_in,s1))//(Contigs_in>>s1>>ContigNo>>s2>>AvgCov>>s3>>ContLen)		//getline(inContigs,tag))
	{
		getline(Contigs_in,cont_s);
		int ctg_cov=0;
		in_ContigsCov>>ctg_cov;
		contigs_info->cov_vt.push_back(ctg_cov);
		if(cont_s[cont_s.size()-1]=='\r'||cont_s[cont_s.size()-1]=='\n')
		{cont_s.resize(cont_s.size()-1);}
		if(cont_s.size()==0)
		{
			getline(Contigs_in,cont_s);
			if(cont_s[cont_s.size()-1]=='\r'||cont_s[cont_s.size()-1]=='\n')
			{cont_s.resize(cont_s.size()-1);}
		}
		nKmer=0;
		num_Contigs++;

		ContLen=cont_s.size();
		contigs_info->contig_sz_vt.push_back(ContLen);
		struct contig_t Contig;
		Init_Contig(cont_s,Contig);
		
		bool Begin=0;
		int hp=0;
		for(int i=0;i<Contig.contigLen-K_size+1;++i)
		{
			get_sub_arr(Contig.contig_bits,Contig.contigLen,i,K_size,seq_t4.kmer);
			f_seq_t4=seq_t4;
			get_rev_comp_seq_arr(f_seq_t4.kmer,K_size,4);
			flip_c=0;
			if(uint64_t_cmp(seq_t4.kmer,f_seq_t4.kmer,4) > 0)
			{
				flip_c=1;
				seq_t4=f_seq_t4;
			}

			hv=MurmurHash64A(seq_t4.kmer,sizeof(seq_t4),0);
			hash_idx=(size_t) (hv%ht_sz);
			struct bucket4 **ptr;
			ptr= &(ht->store_pos[hash_idx]);
			found=look_up_in_a_list4(&seq_t4,&ptr);

			if(found)
			{
				hp=i;
				if(Begin==0)
				{
					o_hit_pos<<hp<<" ";
					Begin=1;
				}

				nKmer++;
				if(RecordKmerCnt)
				{
					(*ptr)->kmer_info.cod=nKmer;
				}
				else
				{
					(*ptr)->kmer_info.cod=(int )i;
				}

				if(((*ptr)->kmer_info.repeat==0)&&(*ptr)->kmer_info.contig_no==0)
				{
					(*ptr)->kmer_info.contig_no=(int)num_Contigs;
				}
				else
				{
					(*ptr)->kmer_info.repeat=1;
					(*ptr)->kmer_info.contig_no=0;
				}

				(*ptr)->kmer_info.flip=flip_c;
			}
		

			if(i==Contig.contigLen-K_size)
			o_hit_pos<<hp<<endl;

		}

		if(RecordKmerCnt)
		{
			contigs_info->kmer_cnt_vt.push_back(nKmer);
		}

	}

	vector<int>::iterator vit=contigs_info->contig_sz_vt.begin();
	num_Contigs=contigs_info->contig_sz_vt.size()-1;
	contigs_info->LengthRank.clear();
	vit++;
	for(size_t i=0;i<num_Contigs;++i)
	{
		contigs_info->LengthRank.push_back(vit);
		vit++;
	}

	contigs_info->total_contigs=num_Contigs;
	sort(contigs_info->LengthRank.begin(),contigs_info->LengthRank.end(),it_cmp);

	time(&read_time);
	cout<<"Kmer marking time: "<<difftime(read_time,beg_time)<<" secs."<<endl;

}


void ContigsRemapping0(struct hashtable0 *ht,int K_size,struct contigs_info * contigs_info,string Contig_Filename,bool RecordKmerCnt)
{
	time_t beg_time,read_time;
	contigs_info->K_size=K_size;
	string in_fname=Contig_Filename;
	ifstream Contigs_in(in_fname.c_str());	
	string ctg_cov_fname="Contigs_Cov.txt",ctg_hp_fname="Contigs_HP.txt";

	if(Contig_Filename=="SuperContigs.txt")
	{
		ctg_cov_fname="SuperContigs_Cov.txt";
		ctg_hp_fname="SuperContigs_HP.txt";
	}
	ifstream in_ContigsCov(ctg_cov_fname.c_str());
	ofstream o_hit_pos(ctg_hp_fname.c_str());
	size_t num_Contigs=0;
	string tag,s,kmer;
	uint64_t seq[100],f_seq[100];
	uint64_t hv;
	size_t hash_idx;
	bool flip_c;
	size_t nKmer=0;
	size_t tot_Kmers=0;
	size_t boundary_kmers=0;
	size_t non_repeat_kmers=0,repeat_kmers=0;
	int Kmer_arr_sz=K_size/32+1;
	int rem1=K_size%32;
	if(rem1==0)
	{Kmer_arr_sz--;}

	size_t ht_sz;
	
	ht_sz=ht->ht_sz;
	
	for(size_t i=0;i<ht->ht_sz;++i)
	{
		size_t list_sz=0;
		struct bucket0 *bktptr=ht->store_pos[i];
		while(bktptr!=NULL)
		{
			bktptr->kmer_info.contig_no=0;
			bktptr->kmer_info.cod=0;
			bktptr->kmer_info.repeat=0;
			bktptr=bktptr->nxt_bucket;
			list_sz++;
		}
	}
		


	bool found;
	cout<<"Marking kmers..."<<endl;
	time(&beg_time);

	contigs_info->contig_sz_vt.clear();
	contigs_info->contig_sz_vt.push_back(0);
	contigs_info->cov_vt.clear();
	contigs_info->cov_vt.push_back(0);
	
	if(RecordKmerCnt)
	{
		contigs_info->kmer_cnt_vt.clear();
		contigs_info->kmer_cnt_vt.push_back(0);
	}
	string s1,s2,s3,cont_s;
	int ContLen;
	//double AvgCov ;

	while(getline(Contigs_in,s1))//(Contigs_in>>s1>>ContigNo>>s2>>AvgCov>>s3>>ContLen)		//getline(inContigs,tag))
	{
		getline(Contigs_in,cont_s);
		int ctg_cov=0;
		in_ContigsCov>>ctg_cov;
		contigs_info->cov_vt.push_back(ctg_cov);
		if(cont_s[cont_s.size()-1]=='\r'||cont_s[cont_s.size()-1]=='\n')
		{cont_s.resize(cont_s.size()-1);}
		if(cont_s.size()==0)
		{
			getline(Contigs_in,cont_s);
			if(cont_s[cont_s.size()-1]=='\r'||cont_s[cont_s.size()-1]=='\n')
			{cont_s.resize(cont_s.size()-1);}
		}
		nKmer=0;

		num_Contigs++;

		ContLen=cont_s.size();
		contigs_info->contig_sz_vt.push_back(ContLen);
		struct contig_t Contig;
		Init_Contig(cont_s,Contig);
		
		bool Begin=0;
		int hp=0;
		for(int i=0;i<Contig.contigLen-K_size+1;++i)
		{
			
			get_sub_arr(Contig.contig_bits,Contig.contigLen,i,K_size,seq);
			memcpy(f_seq,seq,sizeof(uint64_t)*Kmer_arr_sz);
			

			get_rev_comp_seq_arr(f_seq,K_size,Kmer_arr_sz);


			
			flip_c=0;
			if(uint64_t_cmp(seq,f_seq,Kmer_arr_sz) > 0)
			{
				flip_c=1;
			
				memcpy(seq,f_seq,sizeof(uint64_t)*Kmer_arr_sz);
			}


			hv=MurmurHash64A(seq,sizeof(uint64_t)*Kmer_arr_sz,0);
			hash_idx=(size_t) (hv%ht_sz);
			struct bucket0 **ptr;
			ptr= &(ht->store_pos[hash_idx]);
			found=look_up_in_a_list0(seq,&ptr,Kmer_arr_sz);

			if(found)
			{
				hp=i;
				if(Begin==0)
				{
					o_hit_pos<<hp<<" ";
					Begin=1;
				}

				nKmer++;
				tot_Kmers++;
				if(RecordKmerCnt)
				{
					(*ptr)->kmer_info.cod=nKmer;
				}
				else
				{
					(*ptr)->kmer_info.cod=(int )i;
				}

				if(((*ptr)->kmer_info.repeat==0)&&(*ptr)->kmer_info.contig_no==0)
				{
					(*ptr)->kmer_info.contig_no=(int)num_Contigs;
					non_repeat_kmers++;
				}
				else
				{
					(*ptr)->kmer_info.repeat=1;
					(*ptr)->kmer_info.contig_no=0;
					repeat_kmers++;
				}

				(*ptr)->kmer_info.flip=flip_c;
				if(i<100||i>(Contig.contigLen-K_size-100))
				{
					(*ptr)->kmer_info.masked=1;
				}
				else
				{
					(*ptr)->kmer_info.masked=0;
				}
			}
		

			if(i==Contig.contigLen-K_size)
			o_hit_pos<<hp<<endl;

		}

		if(RecordKmerCnt)
		{
			contigs_info->kmer_cnt_vt.push_back(nKmer);
		}

	}
	//cout<<"tot_Kmers "<<tot_Kmers<<endl;
	//cout<<"repeat_kmers "<<repeat_kmers<<endl;
	//cout<<"non_repeat_kmers "<<non_repeat_kmers<<endl;
	vector<int>::iterator vit=contigs_info->contig_sz_vt.begin();
	num_Contigs=contigs_info->contig_sz_vt.size()-1;
	contigs_info->LengthRank.clear();
	vit++;
	for(size_t i=0;i<num_Contigs;++i)
	{
		contigs_info->LengthRank.push_back(vit);
		vit++;
	}

	contigs_info->total_contigs=num_Contigs;
	sort(contigs_info->LengthRank.begin(),contigs_info->LengthRank.end(),it_cmp);

	time(&read_time);
	cout<<"Kmer marking time: "<<difftime(read_time,beg_time)<<" secs."<<endl;

}


void SuperContigsRemapping(struct hashtable *ht,struct hashtable2 *ht2, int K_size,struct contigs_info * contigs_info,string Contig_Filename,bool RecordKmerCnt)
{
	time_t beg_time,read_time;
	//int64_t scaffold_len;
	//int Kmer_arr_sz;
////////////////////////////////
	contigs_info->K_size=K_size;
	string in_fname=Contig_Filename;
	ifstream Contigs_in(in_fname.c_str());
	//ifstream in_ContigsCov("Contigs_Cov.txt");
	ofstream o_hit_pos("SuperContigs_HP.txt");
	size_t num_Contigs=0;
	string tag,s,kmer;
	uint64_t seq,f_seq;
	kmer_t2 seq_t2,f_seq_t2;
	uint64_t hv;
	size_t hash_idx;
	bool flip_c;
	size_t nKmer=0;
	
	size_t ht_sz;
	if(K_size<=32)
	{
		ht_sz=ht->ht_sz;
	}
	else
	{
		ht_sz=ht2->ht_sz;
	}
	if(K_size<=32)
	{
		for(size_t i=0;i<ht->ht_sz;++i)
		{
			size_t list_sz=0;

			struct bucket *bktptr=ht->store_pos[i];
			while(bktptr!=NULL)
			{

				bktptr->kmer_info.contig_no=0;
				bktptr->kmer_info.cod=0;
				bktptr->kmer_info.repeat=0;
				bktptr=bktptr->nxt_bucket;
				list_sz++;
			}

		}
	}
	else
	{
		if(K_size>32&&K_size<=64)
		{
			for(size_t i=0;i<ht2->ht_sz;++i)
			{
				size_t list_sz=0;
				struct bucket2 *bktptr=ht2->store_pos[i];
				while(bktptr!=NULL)
				{

					bktptr->kmer_info.contig_no=0;
					bktptr->kmer_info.cod=0;
					bktptr->kmer_info.repeat=0;
					bktptr=bktptr->nxt_bucket;
					list_sz++;
				}
			}
		}

	}

	bool found;
	cout<<"Marking kmers..."<<endl;
	time(&beg_time);

	contigs_info->contig_sz_vt.clear();
	contigs_info->contig_sz_vt.push_back(0);
	contigs_info->cov_vt.clear();
	contigs_info->cov_vt.push_back(0);
	
	if(RecordKmerCnt)
	{
		contigs_info->kmer_cnt_vt.clear();
		contigs_info->kmer_cnt_vt.push_back(0);
	}
	string s1,s2,s3,cont_s;
	int ContLen;
	//double AvgCov ;

	while(getline(Contigs_in,s1))//(Contigs_in>>s1>>ContigNo>>s2>>AvgCov>>s3>>ContLen)		//getline(inContigs,tag))
	{
		getline(Contigs_in,cont_s);
		int ctg_cov=0;
		//in_ContigsCov>>ctg_cov;
		//contigs_info->cov_vt.push_back(ctg_cov);
		if(cont_s[cont_s.size()-1]=='\r'||cont_s[cont_s.size()-1]=='\n')
		{cont_s.resize(cont_s.size()-1);}
		if(cont_s.size()==0)
		{
			getline(Contigs_in,cont_s);
			if(cont_s[cont_s.size()-1]=='\r'||cont_s[cont_s.size()-1]=='\n')
			{cont_s.resize(cont_s.size()-1);}
		}
		nKmer=0;
		num_Contigs++;

		ContLen=cont_s.size();
		contigs_info->contig_sz_vt.push_back(ContLen);
		//contigs_info->cov_vt.push_back(int(AvgCov+0.5));

	
		bool Begin=0;
		int hp=0;
		int k=0;
		for(int i=0;i<ContLen-K_size+1;++i)
		{
			i=k;

			while(cont_s[i]=='N'&&(i<ContLen-K_size+1))
			{++i;}
			if(i==(ContLen-K_size+1))
			{break;}

			
			for (k=i+1;k<ContLen-K_size+1;++k)
			{
				if(cont_s[k]=='N')
				{break;}
			}
			
			
			string contig_subseq=cont_s.substr(i,k-i);

			struct contig_t Contig;
			Init_Contig(contig_subseq,Contig);
			
			if(K_size<=32)
			{
				for (int j=0;j<((int)contig_subseq.size())-K_size+1;++j)
				{

					get_sub_arr(Contig.contig_bits,Contig.contigLen,j,K_size,&seq);
					f_seq=get_rev_comp_seq(seq,K_size);
					flip_c=0;
					if(seq>f_seq)
					{
						flip_c=1;
						seq=f_seq;
					}

					hv=MurmurHash64A(&seq,sizeof(seq),0);
					hash_idx=(size_t) (hv%ht_sz);
					struct bucket **ptr;
					ptr= &(ht->store_pos[hash_idx]);
					found=look_up_in_a_list(seq,&ptr);

					if(found)
					{

						hp=i+j;
						if(Begin==0)
						{
							o_hit_pos<<hp<<" ";
							Begin=1;
						}
					
						nKmer++;
						if(RecordKmerCnt)
						{
							(*ptr)->kmer_info.cod=nKmer;
						}
						else
						{
							(*ptr)->kmer_info.cod=(int )i+j;
						}
						if((*ptr)->kmer_info.repeat==0&&(*ptr)->kmer_info.contig_no==0)
						{
							(*ptr)->kmer_info.contig_no=(int)num_Contigs;
						}
						else
						{
							(*ptr)->kmer_info.repeat=1;
							(*ptr)->kmer_info.contig_no=0;
						}
						(*ptr)->kmer_info.flip=flip_c;
					
					}
				}

			}
			else
			{

				for (int j=0;j<((int)contig_subseq.size())-K_size+1;++j)
				{

					get_sub_arr(Contig.contig_bits,Contig.contigLen,j,K_size,seq_t2.kmer);
					f_seq_t2=seq_t2;
					get_rev_comp_seq_arr(f_seq_t2.kmer,K_size,2);
					flip_c=0;
					if(uint64_t_cmp(seq_t2.kmer,f_seq_t2.kmer,2) > 0)
					{
						flip_c=1;
						seq_t2=f_seq_t2;
					}

					hv=MurmurHash64A(seq_t2.kmer,sizeof(seq_t2),0);
					hash_idx=(size_t) (hv%ht_sz);
					struct bucket2 **ptr;
					ptr= &(ht2->store_pos[hash_idx]);
					found=look_up_in_a_list2(&seq_t2,&ptr);

					if(found)
					{
						hp=i+j;
						if(Begin==0)
						{
							o_hit_pos<<hp<<" ";
							Begin=1;
						}

						nKmer++;
						if(RecordKmerCnt)
						{
							(*ptr)->kmer_info.cod=nKmer;
						}
						else
						{
							(*ptr)->kmer_info.cod=(int )i+j;
						}

						if(((*ptr)->kmer_info.repeat==0)&&(*ptr)->kmer_info.contig_no==0)
						{
							(*ptr)->kmer_info.contig_no=(int)num_Contigs;
						}
						else
						{
							(*ptr)->kmer_info.repeat=1;
							(*ptr)->kmer_info.contig_no=0;
						}

						(*ptr)->kmer_info.flip=flip_c;
					}
				}

			

			}

		
		}
		if(Begin==0)
		{o_hit_pos<<hp<<" ";}
		o_hit_pos<<hp<<endl;


		if(RecordKmerCnt)
		{
			contigs_info->kmer_cnt_vt.push_back(nKmer);
		}

	}

	vector<int>::iterator vit=contigs_info->contig_sz_vt.begin();
	num_Contigs=contigs_info->contig_sz_vt.size()-1;
	contigs_info->LengthRank.clear();
	vit++;
	for(size_t i=0;i<num_Contigs;++i)
	{
		contigs_info->LengthRank.push_back(vit);
		vit++;
	}

	contigs_info->total_contigs=num_Contigs;
	sort(contigs_info->LengthRank.begin(),contigs_info->LengthRank.end(),it_cmp);

	time(&read_time);
	cout<<"Kmer marking time: "<<difftime(read_time,beg_time)<<" secs."<<endl;

}


void SuperContigsRemapping3(struct hashtable3 *ht,int K_size,struct contigs_info * contigs_info,string Contig_Filename,bool RecordKmerCnt)
{
	time_t beg_time,read_time;
	contigs_info->K_size=K_size;
	string in_fname=Contig_Filename;
	ifstream Contigs_in(in_fname.c_str());
	//ifstream in_ContigsCov("Contigs_Cov.txt");
	ofstream o_hit_pos("SuperContigs_HP.txt");
	size_t num_Contigs=0;
	string tag,s,kmer;
	kmer_t3 seq_t3,f_seq_t3;
	uint64_t hv;
	size_t hash_idx;
	bool flip_c;
	size_t nKmer=0;
	
	size_t ht_sz;
	
	ht_sz=ht->ht_sz;
	
	for(size_t i=0;i<ht->ht_sz;++i)
	{
		size_t list_sz=0;
		struct bucket3 *bktptr=ht->store_pos[i];
		while(bktptr!=NULL)
		{
			bktptr->kmer_info.contig_no=0;
			bktptr->kmer_info.cod=0;
			bktptr->kmer_info.repeat=0;
			bktptr=bktptr->nxt_bucket;
			list_sz++;
		}
	}
		


	bool found;
	cout<<"Marking kmers..."<<endl;
	time(&beg_time);

	contigs_info->contig_sz_vt.clear();
	contigs_info->contig_sz_vt.push_back(0);
	contigs_info->cov_vt.clear();
	contigs_info->cov_vt.push_back(0);
	
	if(RecordKmerCnt)
	{
		contigs_info->kmer_cnt_vt.clear();
		contigs_info->kmer_cnt_vt.push_back(0);
	}
	string s1,s2,s3,cont_s;
	int ContLen;
	//double AvgCov ;

	while(getline(Contigs_in,s1))//(Contigs_in>>s1>>ContigNo>>s2>>AvgCov>>s3>>ContLen)		//getline(inContigs,tag))
	{
		getline(Contigs_in,cont_s);
		int ctg_cov=0;
		//in_ContigsCov>>ctg_cov;
		//contigs_info->cov_vt.push_back(ctg_cov);
		if(cont_s[cont_s.size()-1]=='\r'||cont_s[cont_s.size()-1]=='\n')
		{cont_s.resize(cont_s.size()-1);}
		if(cont_s.size()==0)
		{
			getline(Contigs_in,cont_s);
			if(cont_s[cont_s.size()-1]=='\r'||cont_s[cont_s.size()-1]=='\n')
			{cont_s.resize(cont_s.size()-1);}
		}
		nKmer=0;
		num_Contigs++;

		ContLen=cont_s.size();
		contigs_info->contig_sz_vt.push_back(ContLen);
		
		bool Begin=0;
		int hp=0;


		int k=0;
		for(int i=0;i<ContLen-K_size+1;++i)
		{
			i=k;

			while(cont_s[i]=='N'&&(i<ContLen-K_size+1))
			{++i;}
			if(i==(ContLen-K_size+1))
			{break;}

			
			for (k=i+1;k<ContLen-K_size+1;++k)
			{
				if(cont_s[k]=='N')
				{break;}
			}
			
			
			string contig_subseq=cont_s.substr(i,k-i);

			struct contig_t Contig;
			Init_Contig(contig_subseq,Contig);
			for (int j=0;j<((int)contig_subseq.size())-K_size+1;++j)
			{
				get_sub_arr(Contig.contig_bits,Contig.contigLen,j,K_size,seq_t3.kmer);
				f_seq_t3=seq_t3;
				get_rev_comp_seq_arr(f_seq_t3.kmer,K_size,3);
				flip_c=0;
				if(uint64_t_cmp(seq_t3.kmer,f_seq_t3.kmer,3) > 0)
				{
					flip_c=1;
					seq_t3=f_seq_t3;
				}

				hv=MurmurHash64A(seq_t3.kmer,sizeof(seq_t3),0);
				hash_idx=(size_t) (hv%ht_sz);
				struct bucket3 **ptr;
				ptr= &(ht->store_pos[hash_idx]);
				found=look_up_in_a_list3(&seq_t3,&ptr);

				if(found)
				{
					hp=i+j;
					if(Begin==0)
					{
						o_hit_pos<<hp<<" ";
						Begin=1;
					}

					nKmer++;
					if(RecordKmerCnt)
					{
						(*ptr)->kmer_info.cod=nKmer;
					}
					else
					{
						(*ptr)->kmer_info.cod=(int )i+j;
					}

					if(((*ptr)->kmer_info.repeat==0)&&(*ptr)->kmer_info.contig_no==0)
					{
						(*ptr)->kmer_info.contig_no=(int)num_Contigs;
					}
					else
					{
						(*ptr)->kmer_info.repeat=1;
						(*ptr)->kmer_info.contig_no=0;
					}

					(*ptr)->kmer_info.flip=flip_c;
				}
			}

			//if(i==Contig.contigLen-K_size)
			//o_hit_pos<<hp<<endl;

		}
		if(Begin==0)
		{o_hit_pos<<hp<<" ";}
		o_hit_pos<<hp<<endl;

		if(RecordKmerCnt)
		{
			contigs_info->kmer_cnt_vt.push_back(nKmer);
		}

	}

	vector<int>::iterator vit=contigs_info->contig_sz_vt.begin();
	num_Contigs=contigs_info->contig_sz_vt.size()-1;
	contigs_info->LengthRank.clear();
	vit++;
	for(size_t i=0;i<num_Contigs;++i)
	{
		contigs_info->LengthRank.push_back(vit);
		vit++;
	}

	contigs_info->total_contigs=num_Contigs;
	sort(contigs_info->LengthRank.begin(),contigs_info->LengthRank.end(),it_cmp);

	time(&read_time);
	cout<<"Kmer marking time: "<<difftime(read_time,beg_time)<<" secs."<<endl;

}

void SuperContigsRemapping4(struct hashtable4 *ht,int K_size,struct contigs_info * contigs_info,string Contig_Filename,bool RecordKmerCnt)
{
	time_t beg_time,read_time;
	contigs_info->K_size=K_size;
	string in_fname=Contig_Filename;
	ifstream Contigs_in(in_fname.c_str());
	//ifstream in_ContigsCov("Contigs_Cov.txt");
	ofstream o_hit_pos("SuperContigs_HP.txt");
	size_t num_Contigs=0;
	string tag,s,kmer;
	kmer_t4 seq_t4,f_seq_t4;
	uint64_t hv;
	size_t hash_idx;
	bool flip_c;
	size_t nKmer=0;
	
	size_t ht_sz;
	
	ht_sz=ht->ht_sz;
	
	for(size_t i=0;i<ht->ht_sz;++i)
	{
		size_t list_sz=0;
		struct bucket4 *bktptr=ht->store_pos[i];
		while(bktptr!=NULL)
		{
			bktptr->kmer_info.contig_no=0;
			bktptr->kmer_info.cod=0;
			bktptr->kmer_info.repeat=0;
			bktptr=bktptr->nxt_bucket;
			list_sz++;
		}
	}
		


	bool found;
	cout<<"Marking kmers..."<<endl;
	time(&beg_time);

	contigs_info->contig_sz_vt.clear();
	contigs_info->contig_sz_vt.push_back(0);
	contigs_info->cov_vt.clear();
	contigs_info->cov_vt.push_back(0);
	
	if(RecordKmerCnt)
	{
		contigs_info->kmer_cnt_vt.clear();
		contigs_info->kmer_cnt_vt.push_back(0);
	}
	string s1,s2,s3,cont_s;
	int ContLen;
	//double AvgCov ;

	while(getline(Contigs_in,s1))//(Contigs_in>>s1>>ContigNo>>s2>>AvgCov>>s3>>ContLen)		//getline(inContigs,tag))
	{
		getline(Contigs_in,cont_s);
		int ctg_cov=0;
		//in_ContigsCov>>ctg_cov;
		//contigs_info->cov_vt.push_back(ctg_cov);
		if(cont_s[cont_s.size()-1]=='\r'||cont_s[cont_s.size()-1]=='\n')
		{cont_s.resize(cont_s.size()-1);}
		if(cont_s.size()==0)
		{
			getline(Contigs_in,cont_s);
			if(cont_s[cont_s.size()-1]=='\r'||cont_s[cont_s.size()-1]=='\n')
			{cont_s.resize(cont_s.size()-1);}
		}
		nKmer=0;
		num_Contigs++;

		ContLen=cont_s.size();
		contigs_info->contig_sz_vt.push_back(ContLen);
		//struct contig_t Contig;
		//Init_Contig(cont_s,Contig);
		
		bool Begin=0;
		int hp=0;


		int k=0;
		for(int i=0;i<ContLen-K_size+1;++i)
		{
			i=k;

			while(cont_s[i]=='N'&&(i<ContLen-K_size+1))
			{++i;}
			if(i==(ContLen-K_size+1))
			{break;}

			
			for (k=i+1;k<ContLen-K_size+1;++k)
			{
				if(cont_s[k]=='N')
				{break;}
			}
			
			
			string contig_subseq=cont_s.substr(i,k-i);

			struct contig_t Contig;
			Init_Contig(contig_subseq,Contig);
			for (int j=0;j<((int)contig_subseq.size())-K_size+1;++j)
			{

				get_sub_arr(Contig.contig_bits,Contig.contigLen,j,K_size,seq_t4.kmer);
				f_seq_t4=seq_t4;
				get_rev_comp_seq_arr(f_seq_t4.kmer,K_size,4);
				flip_c=0;
				if(uint64_t_cmp(seq_t4.kmer,f_seq_t4.kmer,4) > 0)
				{
					flip_c=1;
					seq_t4=f_seq_t4;
				}

				hv=MurmurHash64A(seq_t4.kmer,sizeof(seq_t4),0);
				hash_idx=(size_t) (hv%ht_sz);
				struct bucket4 **ptr;
				ptr= &(ht->store_pos[hash_idx]);
				found=look_up_in_a_list4(&seq_t4,&ptr);

				if(found)
				{
					hp=i+j;
					if(Begin==0)
					{
						o_hit_pos<<hp<<" ";
						Begin=1;
					}

					nKmer++;
					if(RecordKmerCnt)
					{
						(*ptr)->kmer_info.cod=nKmer;
					}
					else
					{
						(*ptr)->kmer_info.cod=(int )i+j;
					}

					if(((*ptr)->kmer_info.repeat==0)&&(*ptr)->kmer_info.contig_no==0)
					{
						(*ptr)->kmer_info.contig_no=(int)num_Contigs;
					}
					else
					{
						(*ptr)->kmer_info.repeat=1;
						(*ptr)->kmer_info.contig_no=0;
					}

					(*ptr)->kmer_info.flip=flip_c;
				}
			}

	//		if(i==Contig.contigLen-K_size)
		//	o_hit_pos<<hp<<endl;

		}
			if(Begin==0)
		{o_hit_pos<<hp<<" ";}
		o_hit_pos<<hp<<endl;

		if(RecordKmerCnt)
		{
			contigs_info->kmer_cnt_vt.push_back(nKmer);
		}

	}

	vector<int>::iterator vit=contigs_info->contig_sz_vt.begin();
	num_Contigs=contigs_info->contig_sz_vt.size()-1;
	contigs_info->LengthRank.clear();
	vit++;
	for(size_t i=0;i<num_Contigs;++i)
	{
		contigs_info->LengthRank.push_back(vit);
		vit++;
	}

	contigs_info->total_contigs=num_Contigs;
	sort(contigs_info->LengthRank.begin(),contigs_info->LengthRank.end(),it_cmp);

	time(&read_time);
	cout<<"Kmer marking time: "<<difftime(read_time,beg_time)<<" secs."<<endl;

}


void BuildContigAdjacency(hashtable *ht1, hashtable2 *ht2, struct contigs_info *contigs_info,int K_size, string ContigFilename)
{
	ifstream in_ctg(ContigFilename.c_str());
	ofstream out_ctg_graph("ContigGraph.txt");
	string tag,contig_s;
	int contig_no=0;
	contigs_info->contig_adjacency_left.resize(contigs_info->total_contigs+1);
	contigs_info->contig_adjacency_right.resize(contigs_info->total_contigs+1);

	while(getline(in_ctg,tag))
	{
		contig_no++;
	//	cout<<contig_no<<endl;
	
		getline(in_ctg,contig_s);

		if(K_size<=32)
		{
			uint64_t seq,f_seq,hv,hash_idx;
			bool flip_c=0,found=0;
			struct bucket **ptr;
			int gMax=100;
			if((contig_s.size()-K_size)<100)
			{
				gMax=contig_s.size()-K_size;
			}
			for(int j=0;j<=gMax;++j)
			{
				string K_mer=contig_s.substr(j,K_size);
				str2bitsarr(K_mer.c_str(),K_size,&seq,1);

				f_seq=get_rev_comp_seq(seq,K_size);
				flip_c=0;
				if(seq>f_seq)
				{
					flip_c=1;
					seq=f_seq;
				}

				hv=MurmurHash64A(&seq,sizeof(seq),0);
				hash_idx=(size_t) (hv%(ht1->ht_sz));

				ptr= &(ht1->store_pos[hash_idx]);
				found=look_up_in_a_list(seq,&ptr);

				if(found==0)
				{

					//cout<<"searching error."<<endl;
					continue;
				}
				else
				{
					break;
				}
			}

			if(found==0)
			{

			//	cout<<"searching error."<<endl;
				return;
			}

			if((*ptr)->kmer_info.flip==0)
			{
				//search left, append to left
				edge_node* edge_ptr=(*ptr)->kmer_info.left;
				while(edge_ptr!=NULL)
				{
					uint64_t t_kmer=seq,f_kmer;
					uint64_t t=0;
					bool flip_nc=0;

					for(int j=0;j<=edge_ptr->len;++j)
					{
						int left=(int)((edge_ptr->edge)>>2*j);


						switch(left&0x3)
						{
							case 0:
								t_kmer>>=2;
								break;
							case 1:
								t_kmer>>=2;
								t=((uint64_t)1)<<((K_size-1)*2);
								t_kmer|=t;

								break;
							case 2:
								t_kmer>>=2;
								t=((uint64_t)2)<<((K_size-1)*2);
								t_kmer|=t;

								break;
							case 3:
								t_kmer>>=2;
								t=((uint64_t)3)<<((K_size-1)*2);
								t_kmer|=t;

								break;

						}


					}

					flip_nc=0;

					f_kmer=get_rev_comp_seq(t_kmer,K_size);
					if(t_kmer>f_kmer)
					{
						t_kmer=f_kmer;
						flip_nc=!flip_nc;
					}

					hv=MurmurHash64A(&t_kmer,sizeof(t_kmer),0);

					hash_idx=(size_t) (hv%ht1->ht_sz);


					struct bucket ** nc_ptr;

					nc_ptr= &(ht1->store_pos[hash_idx]);

					found=look_up_in_a_list(t_kmer,&nc_ptr);


					if(found==0)
					{
						edge_ptr=edge_ptr->nxt_edge;
						//cout<<"WarningL."<<endl;
						continue;
					}
					int n_contig_no=(*nc_ptr)->kmer_info.contig_no;
					
					if(n_contig_no==0)
					{
						edge_ptr=edge_ptr->nxt_edge;
						continue;
					}

					//if((*nc_ptr)->kmer_info.flip^flip_nc)
					bool t_flip=(*nc_ptr)->kmer_info.flip;

					char edge_cstr[200];
					uint64_t edgebits=edge_ptr->edge;
					int edge_len= edge_ptr->len;
					edge_len++;
					string edge_s = bitsarr2str(&edgebits,edge_len,edge_cstr,1);


					if((t_flip^flip_nc)==1)
					{

						if(contigs_info->contig_adjacency_left[contig_no].count(-n_contig_no)==0)
						{
							contigs_info->contig_adjacency_left[contig_no][-n_contig_no].cov=edge_ptr->edge_cov;//1;
							contigs_info->contig_adjacency_left[contig_no][-n_contig_no].bridge=edge_s;
						}
						else
						{
						//	contigs_info->contig_adjacency_left[contig_no][-n_contig_no].cov++;
						}

					}
					else
					{
						if(contigs_info->contig_adjacency_left[contig_no].count(n_contig_no)==0)
						{
							contigs_info->contig_adjacency_left[contig_no][n_contig_no].cov=edge_ptr->edge_cov;
							contigs_info->contig_adjacency_left[contig_no][n_contig_no].bridge=edge_s;

						}
						else
						{
							//contigs_info->contig_adjacency_left[contig_no][n_contig_no].cov++;
						}

					}




					edge_ptr=edge_ptr->nxt_edge;

				}

			}
			else
			{
				//search right, app to adj left


				edge_node* edge_ptr=(*ptr)->kmer_info.right;
				while(edge_ptr!=NULL)
				{
					uint64_t t_kmer=seq,f_kmer;
					uint64_t t=0;
					bool flip_nc=0;

					for(int j=edge_ptr->len;j>=0;--j)
					{
						int left=(int)((edge_ptr->edge)>>2*j);


						//left=(edge_ptr->edge)>>2*j;
						switch(left&0x3)
						{
							case 0:
								t=3;
								t<<=(K_size-1)*2;
								t_kmer&=(~t);
								t_kmer<<=2;//



								break;
							case 1:
								t=3;
								t<<=(K_size-1)*2;
								t_kmer&=(~t);
								t_kmer<<=2;//
								t=1;

								t_kmer|=t;

								break;
							case 2:
								t=3;
								t<<=(K_size-1)*2;
								t_kmer&=(~t);
								t_kmer<<=2;
								t=2;

								t_kmer|=t;

								break;
							case 3:
								t=3;
								t<<=(K_size-1)*2;
								t_kmer&=(~t);
								t_kmer<<=2;
								t=3;

								t_kmer|=t;

								break;

						}

					}

					flip_nc=0;

					f_kmer=get_rev_comp_seq(t_kmer,K_size);
					if(t_kmer>f_kmer)
					{
						t_kmer=f_kmer;
						flip_nc=!flip_nc;
					}

					hv=MurmurHash64A(&t_kmer,sizeof(t_kmer),0);

					hash_idx=(size_t) (hv%ht1->ht_sz);


					struct bucket ** nc_ptr;

					nc_ptr= &(ht1->store_pos[hash_idx]);

					found=look_up_in_a_list(t_kmer,&nc_ptr);


					if(found==0)
					{
						edge_ptr=edge_ptr->nxt_edge;
						//cout<<"WarningL."<<endl;
						continue;
					}
					int n_contig_no=(*nc_ptr)->kmer_info.contig_no;
					if(n_contig_no==0)
					{
						edge_ptr=edge_ptr->nxt_edge;
						continue;
					}

					bool t_flip=(*nc_ptr)->kmer_info.flip;
					char edge_cstr[200];
					uint64_t edgebits=edge_ptr->edge;
					int edge_len= edge_ptr->len;
					edge_len++;
					string edge_s = bitsarr2str(&edgebits,edge_len,edge_cstr,1);
					
					if((t_flip^flip_nc)==0)
					{



						if(contigs_info->contig_adjacency_left[contig_no].count(-n_contig_no)==0)
						{
							contigs_info->contig_adjacency_left[contig_no][-n_contig_no].cov=edge_ptr->edge_cov;
							contigs_info->contig_adjacency_left[contig_no][-n_contig_no].bridge=edge_s;

						}
						else
						{
							//contigs_info->contig_adjacency_left[contig_no][-n_contig_no].cov++;
						}

					}
					else
					{
						if(contigs_info->contig_adjacency_left[contig_no].count(n_contig_no)==0)
						{
							contigs_info->contig_adjacency_left[contig_no][n_contig_no].cov=edge_ptr->edge_cov;;
							contigs_info->contig_adjacency_left[contig_no][n_contig_no].bridge=edge_s;
						}
						else
						{
							//contigs_info->contig_adjacency_left[contig_no][n_contig_no].cov++;
						}

					}




					edge_ptr=edge_ptr->nxt_edge;

				}

			}


			//contig endpoint

			int gMin=contig_s.size()-K_size-100;
			if(gMin<0)
			{gMin=0;}
			for(int j=contig_s.size()-K_size;j>=gMin;--j)
			{
				string K_mer=contig_s.substr(j,K_size);
				str2bitsarr(K_mer.c_str(),K_size,&seq,1);

				f_seq=get_rev_comp_seq(seq,K_size);
				flip_c=0;
				if(seq>f_seq)
				{
					flip_c=1;
					seq=f_seq;
				}

				hv=MurmurHash64A(&seq,sizeof(seq),0);
				hash_idx=(size_t) (hv%(ht1->ht_sz));

				ptr= &(ht1->store_pos[hash_idx]);
				found=look_up_in_a_list(seq,&ptr);
				if(found==0)
				{
					continue;
				}
				else
				{
					break;
				}
			}


			if(found==0)
			{

			//	cout<<"searching error."<<endl;
				return ;
			}


			if((*ptr)->kmer_info.flip==0)
			{
				//search right, app to adj right

				edge_node* edge_ptr=(*ptr)->kmer_info.right;
				while(edge_ptr!=NULL)
				{
					uint64_t t_kmer=seq,f_kmer;
					uint64_t t=0;
					bool flip_nc=0;

					for(int j=edge_ptr->len;j>=0;--j)
					{
						int left=(int)((edge_ptr->edge)>>2*j);


					//	left=(edge_ptr->edge)>>2*j;
						switch(left&0x3)
						{
							case 0:
								t=3;
								t<<=(K_size-1)*2;
								t_kmer&=(~t);
								t_kmer<<=2;//



								break;
							case 1:
								t=3;
								t<<=(K_size-1)*2;
								t_kmer&=(~t);
								t_kmer<<=2;//
								t=1;

								t_kmer|=t;

								break;
							case 2:
								t=3;
								t<<=(K_size-1)*2;
								t_kmer&=(~t);
								t_kmer<<=2;
								t=2;

								t_kmer|=t;

								break;
							case 3:
								t=3;
								t<<=(K_size-1)*2;
								t_kmer&=(~t);
								t_kmer<<=2;
								t=3;

								t_kmer|=t;

								break;

						}

					}

					flip_nc=0;

					f_kmer=get_rev_comp_seq(t_kmer,K_size);
					if(t_kmer>f_kmer)
					{
						t_kmer=f_kmer;
						flip_nc=!flip_nc;
					}

					hv=MurmurHash64A(&t_kmer,sizeof(t_kmer),0);

					hash_idx=(size_t) (hv%ht1->ht_sz);


					struct bucket ** nc_ptr;

					nc_ptr= &(ht1->store_pos[hash_idx]);

					found=look_up_in_a_list(t_kmer,&nc_ptr);


					if(found==0)
					{
						edge_ptr=edge_ptr->nxt_edge;
						//cout<<"WarningR."<<endl;
						continue;
					}
					int n_contig_no=(*nc_ptr)->kmer_info.contig_no;
					if(n_contig_no==0)
					{
						edge_ptr=edge_ptr->nxt_edge;
						continue;
					}
					bool t_flip=(*nc_ptr)->kmer_info.flip;

					char edge_cstr[200];
					uint64_t edgebits=edge_ptr->edge;
					//edgebits=get_rev_comp_seq(edgebits,edge_ptr->len);
					int edge_len= edge_ptr->len;
					edge_len++;
					string edge_s = bitsarr2str(&edgebits,edge_len,edge_cstr,1);

					if((t_flip^flip_nc)==1)
					{

						if(contigs_info->contig_adjacency_right[contig_no].count(-n_contig_no)==0)
						{
							contigs_info->contig_adjacency_right[contig_no][-n_contig_no].cov=edge_ptr->edge_cov;
							contigs_info->contig_adjacency_right[contig_no][-n_contig_no].bridge=edge_s;

						}
						else
						{
						//	contigs_info->contig_adjacency_right[contig_no][-n_contig_no].cov++;
						}

					}
					else
					{
						if(contigs_info->contig_adjacency_right[contig_no].count(n_contig_no)==0)
						{
							contigs_info->contig_adjacency_right[contig_no][n_contig_no].cov=edge_ptr->edge_cov;
							contigs_info->contig_adjacency_right[contig_no][n_contig_no].bridge=edge_s;

						}
						else
						{
							//contigs_info->contig_adjacency_right[contig_no][n_contig_no].cov++;
						}

					}




					edge_ptr=edge_ptr->nxt_edge;

				}

			}
			else
			{
				//search left, app to adj right

				edge_node* edge_ptr=(*ptr)->kmer_info.left;
				while(edge_ptr!=NULL)
				{
					uint64_t t_kmer=seq,f_kmer;
					uint64_t t=0;
					bool flip_nc=0;

					for(int j=0;j<=edge_ptr->len;++j)
					{
						int left=(int)((edge_ptr->edge)>>2*j);


						switch(left&0x3)
						{
							case 0:
								t_kmer>>=2;
								break;
							case 1:
								t_kmer>>=2;
								t=((uint64_t)1)<<((K_size-1)*2);
								t_kmer|=t;

								break;
							case 2:
								t_kmer>>=2;
								t=((uint64_t)2)<<((K_size-1)*2);
								t_kmer|=t;

								break;
							case 3:
								t_kmer>>=2;
								t=((uint64_t)3)<<((K_size-1)*2);
								t_kmer|=t;

								break;

						}


					}

					flip_nc=0;

					f_kmer=get_rev_comp_seq(t_kmer,K_size);
					if(t_kmer>f_kmer)
					{
						t_kmer=f_kmer;
						flip_nc=!flip_nc;
					}

					hv=MurmurHash64A(&t_kmer,sizeof(t_kmer),0);

					hash_idx=(size_t) (hv%ht1->ht_sz);


					struct bucket ** nc_ptr;

					nc_ptr= &(ht1->store_pos[hash_idx]);

					found=look_up_in_a_list(t_kmer,&nc_ptr);


					if(found==0)
					{
						edge_ptr=edge_ptr->nxt_edge;
						//cout<<"WarningR."<<endl;
						continue;
					}
					int n_contig_no=(*nc_ptr)->kmer_info.contig_no;
					if(n_contig_no==0)
					{
						edge_ptr=edge_ptr->nxt_edge;
						continue;
					}
					bool t_flip=(*nc_ptr)->kmer_info.flip;

					char edge_cstr[200];
					uint64_t edgebits=edge_ptr->edge;
					int edge_len= edge_ptr->len;
					edge_len++;
					edgebits=get_rev_comp_seq(edgebits,edge_len);
					string edge_s = bitsarr2str(&edgebits,edge_len,edge_cstr,1);

					if((t_flip^flip_nc)==0)
					{

						if(contigs_info->contig_adjacency_right[contig_no].count(-n_contig_no)==0)
						{
							contigs_info->contig_adjacency_right[contig_no][-n_contig_no].cov=edge_ptr->edge_cov;
							contigs_info->contig_adjacency_right[contig_no][-n_contig_no].bridge=edge_s;

						}
						else
						{
							//contigs_info->contig_adjacency_right[contig_no][-n_contig_no].cov++;
						}

					}
					else
					{
						if(contigs_info->contig_adjacency_right[contig_no].count(n_contig_no)==0)
						{
							contigs_info->contig_adjacency_right[contig_no][n_contig_no].cov=edge_ptr->edge_cov;
							contigs_info->contig_adjacency_right[contig_no][n_contig_no].bridge=edge_s;

						}
						else
						{
						//	contigs_info->contig_adjacency_right[contig_no][n_contig_no].cov++;
						}

					}




					edge_ptr=edge_ptr->nxt_edge;

				}

			}


		}




		//64

		if(K_size>32&&K_size<=64)
		{


			kmer_t2 seq,f_seq;
			uint64_t hv;
			size_t hash_idx;
			bool flip_c=0,found=0;
			struct bucket2 **ptr;
			int gMax=100;
			if((contig_s.size()-K_size)<100)
			{
				gMax=contig_s.size()-K_size;
			}
			for(int j=0;j<=gMax;++j)
			{

				string K_mer=contig_s.substr(j,K_size);
				str2bitsarr(K_mer.c_str(),K_size,seq.kmer,2);
				f_seq=seq;
				get_rev_comp_seq_arr(f_seq.kmer,K_size,2);
				flip_c=0;
				if(uint64_t_cmp(seq.kmer,f_seq.kmer,2)>0)
				{
					flip_c=1;
					seq=f_seq;
				}

				hv=MurmurHash64A(&seq,sizeof(seq),0);
				hash_idx=(size_t) (hv%(ht2->ht_sz));

				ptr= &(ht2->store_pos[hash_idx]);
				found=look_up_in_a_list2(&seq,&ptr);
				if(found==0)
				{continue;}
				else
				{break;}
			}

			if(found==0)
			{

				//cout<<"searching error."<<endl;
				//return ;
			}
			if(found)
			{
				if((*ptr)->kmer_info.flip==0)
				{
					//search left, append to left
					edge_node* edge_ptr=(*ptr)->kmer_info.left;
			
					while(edge_ptr!=NULL)
					{
						kmer_t2 t_kmer=seq,f_kmer;
						uint64_t t=0;
						bool flip_nc=0;

						for(int j=0;j<=edge_ptr->len;++j)
						{
							int left=(int)((edge_ptr->edge)>>2*j);



							switch(left&0x3)
							{
								case 0:
									R_shift_NB(t_kmer.kmer,2,2);


									break;
								case 1:
									R_shift_NB(t_kmer.kmer,2,2);
									t=((uint64_t)1)<<((K_size-1-32)*2);
									t_kmer.kmer[0]|=t;

									break;
								case 2:
									R_shift_NB(t_kmer.kmer,2,2);
									t=((uint64_t)2)<<((K_size-1-32)*2);
									t_kmer.kmer[0]|=t;

									break;
								case 3:
									R_shift_NB(t_kmer.kmer,2,2);
									t=((uint64_t)3)<<((K_size-1-32)*2);
									t_kmer.kmer[0]|=t;

									break;


							}


						}

						flip_nc=0;
						f_kmer=t_kmer;
						get_rev_comp_seq_arr(f_kmer.kmer,K_size,2);
						if(uint64_t_cmp(t_kmer.kmer,f_kmer.kmer,2)>0)
						{
							t_kmer=f_kmer;
							flip_nc=!flip_nc;
						}

						hv=MurmurHash64A(&t_kmer,sizeof(t_kmer),0);

						hash_idx=(size_t) (hv%ht2->ht_sz);


						struct bucket2 ** nc_ptr;

						nc_ptr= &(ht2->store_pos[hash_idx]);

						found=look_up_in_a_list2(&t_kmer,&nc_ptr);


						if(found==0)
						{
							edge_ptr=edge_ptr->nxt_edge;
							//cout<<"WarningL."<<endl;
							continue;
						}
						int n_contig_no=(*nc_ptr)->kmer_info.contig_no;
						if(n_contig_no==0)
						{
							edge_ptr=edge_ptr->nxt_edge;
							continue;
						}
						//if((*nc_ptr)->kmer_info.flip^flip_nc)
						bool t_flip=(*nc_ptr)->kmer_info.flip;

						char edge_cstr[200];
						uint64_t edgebits=edge_ptr->edge;

							int edge_len= edge_ptr->len;
							edge_len++;
							string edge_s = bitsarr2str(&edgebits,edge_len,edge_cstr,1);
							

						if((t_flip^flip_nc)==1)
						{

							if(contigs_info->contig_adjacency_left[contig_no].count(-n_contig_no)==0)
							{
								contigs_info->contig_adjacency_left[contig_no][-n_contig_no].cov=edge_ptr->edge_cov;
								contigs_info->contig_adjacency_left[contig_no][-n_contig_no].bridge=edge_s;
					
							}
							else
							{
							//	contigs_info->contig_adjacency_left[contig_no][-n_contig_no].cov++;
							}

						}
						else
						{
							if(contigs_info->contig_adjacency_left[contig_no].count(n_contig_no)==0)
							{
								contigs_info->contig_adjacency_left[contig_no][n_contig_no].cov=edge_ptr->edge_cov;
								contigs_info->contig_adjacency_left[contig_no][n_contig_no].bridge=edge_s;

							}
							else
							{
								//contigs_info->contig_adjacency_left[contig_no][n_contig_no].cov++;
							}

						}




						edge_ptr=edge_ptr->nxt_edge;

					}

				}
				else
				{
					//search right, app to adj left


					edge_node* edge_ptr=(*ptr)->kmer_info.right;
					while(edge_ptr!=NULL)
					{
						kmer_t2 t_kmer=seq,f_kmer;
						uint64_t t=0;
						bool flip_nc=0;

						for(int j=edge_ptr->len;j>=0;--j)
						{
							int right=(int)((edge_ptr->edge)>>2*j);


						//	right=(edge_ptr->edge)>>2*j;



							switch(right&0x3)
							{
								case 0:
									t=3;
									t<<=(K_size-1-32)*2;
									t_kmer.kmer[0]&=(~t);
									L_shift_NB(t_kmer.kmer,2,2);//



									break;
								case 1:
									t=3;
									t<<=(K_size-1-32)*2;
									t_kmer.kmer[0]&=(~t);
									L_shift_NB(t_kmer.kmer,2,2);//
									t=1;

									t_kmer.kmer[1]|=t;

									break;
								case 2:
									t=3;
									t<<=(K_size-1-32)*2;
									t_kmer.kmer[0]&=(~t);
									L_shift_NB(t_kmer.kmer,2,2);//
									t=2;

									t_kmer.kmer[1]|=t;

									break;
								case 3:
									t=3;
									t<<=(K_size-1-32)*2;
									t_kmer.kmer[0]&=(~t);
									L_shift_NB(t_kmer.kmer,2,2);//
									t=3;

									t_kmer.kmer[1]|=t;

									break;

							}

						}

						flip_nc=0;

						f_kmer=t_kmer;
						get_rev_comp_seq_arr(f_kmer.kmer,K_size,2);
						if(uint64_t_cmp(t_kmer.kmer,f_kmer.kmer,2)>0)
						{
							t_kmer=f_kmer;
							flip_nc=!flip_nc;
						}

						hv=MurmurHash64A(&t_kmer,sizeof(t_kmer),0);

						hash_idx=(size_t) (hv%ht2->ht_sz);


						struct bucket2 ** nc_ptr;

						nc_ptr= &(ht2->store_pos[hash_idx]);

						found=look_up_in_a_list2(&t_kmer,&nc_ptr);


						if(found==0)
						{
							edge_ptr=edge_ptr->nxt_edge;
							//cout<<"WarningL."<<endl;
							continue;
						}
						int n_contig_no=(*nc_ptr)->kmer_info.contig_no;
						if(n_contig_no==0)
						{
							edge_ptr=edge_ptr->nxt_edge;
							continue;
						}
						//if((*nc_ptr)->kmer_info.flip^flip_nc==0)
						bool t_flip=(*nc_ptr)->kmer_info.flip;
						char edge_cstr[200];
						uint64_t edgebits=edge_ptr->edge;
						int edge_len= edge_ptr->len;
						edge_len++;
						edgebits=get_rev_comp_seq(edgebits,edge_len);
						string edge_s = bitsarr2str(&edgebits,edge_len,edge_cstr,1);
						
						if((t_flip^flip_nc)==0)
						{

							if(contigs_info->contig_adjacency_left[contig_no].count(-n_contig_no)==0)
							{
								contigs_info->contig_adjacency_left[contig_no][-n_contig_no].cov=edge_ptr->edge_cov;;
								contigs_info->contig_adjacency_left[contig_no][-n_contig_no].bridge=edge_s;

							}
							else
							{
								//contigs_info->contig_adjacency_left[contig_no][-n_contig_no].cov++;
							}

						}
						else
						{
							if(contigs_info->contig_adjacency_left[contig_no].count(n_contig_no)==0)
							{
								contigs_info->contig_adjacency_left[contig_no][n_contig_no].cov=edge_ptr->edge_cov;;
								contigs_info->contig_adjacency_left[contig_no][n_contig_no].bridge=edge_s;
					
							}
							else
							{
							//	contigs_info->contig_adjacency_left[contig_no][n_contig_no].cov++;
							}

						}




						edge_ptr=edge_ptr->nxt_edge;

					}

				}

			}
			//contig endpoint

			int gMin=contig_s.size()-K_size-100;
			if(gMin<0)
			{gMin=0;}

			for(int j=contig_s.size()-K_size;j>=gMin;--j)
			{
				string K_mer=contig_s.substr(j,K_size);
				str2bitsarr(K_mer.c_str(),K_size,seq.kmer,2);

				f_seq=seq;
				get_rev_comp_seq_arr(f_seq.kmer,K_size,2);
				flip_c=0;
				if(uint64_t_cmp(seq.kmer,f_seq.kmer,2)>0)
				{
					flip_c=1;
					seq=f_seq;
				}

				hv=MurmurHash64A(&seq,sizeof(seq),0);
				hash_idx=(size_t) (hv%(ht2->ht_sz));

				ptr= &(ht2->store_pos[hash_idx]);
				found=look_up_in_a_list2(&seq,&ptr);
				if(found==0)
				{
					continue;
				}
				else
				{
					break;
				}
			}

			if(found==0)
			{

			//	cout<<"searching error."<<endl;
			//	return ;
			}

			if(found)
			{
				if((*ptr)->kmer_info.flip==0)
				{
					//search right, app to adj right

					edge_node* edge_ptr=(*ptr)->kmer_info.right;
					while(edge_ptr!=NULL)
					{
						kmer_t2 t_kmer=seq,f_kmer;
						uint64_t t=0;
						bool flip_nc=0;

						for(int j=edge_ptr->len;j>=0;--j)
						{
							int right=(int)((edge_ptr->edge)>>2*j);


						



							switch(right&0x3)
							{
								case 0:
									t=3;
									t<<=(K_size-1-32)*2;
									t_kmer.kmer[0]&=(~t);
									L_shift_NB(t_kmer.kmer,2,2);//



									break;
								case 1:
									t=3;
									t<<=(K_size-1-32)*2;
									t_kmer.kmer[0]&=(~t);
									L_shift_NB(t_kmer.kmer,2,2);//
									t=1;

									t_kmer.kmer[1]|=t;

									break;
								case 2:
									t=3;
									t<<=(K_size-1-32)*2;
									t_kmer.kmer[0]&=(~t);
									L_shift_NB(t_kmer.kmer,2,2);//
									t=2;

									t_kmer.kmer[1]|=t;

									break;
								case 3:
									t=3;
									t<<=(K_size-1-32)*2;
									t_kmer.kmer[0]&=(~t);
									L_shift_NB(t_kmer.kmer,2,2);//
									t=3;

									t_kmer.kmer[1]|=t;

									break;

							}

						}

						flip_nc=0;

						f_kmer=t_kmer;
						get_rev_comp_seq_arr(f_kmer.kmer,K_size,2);
						if(uint64_t_cmp(t_kmer.kmer,f_kmer.kmer,2)>0)
						{
							t_kmer=f_kmer;
							flip_nc=!flip_nc;
						}

						hv=MurmurHash64A(&t_kmer,sizeof(t_kmer),0);

						hash_idx=(size_t) (hv%ht2->ht_sz);


						struct bucket2 ** nc_ptr;

						nc_ptr= &(ht2->store_pos[hash_idx]);

						found=look_up_in_a_list2(&t_kmer,&nc_ptr);


						if(found==0)
						{
							edge_ptr=edge_ptr->nxt_edge;
							//cout<<"WarningL."<<endl;
							continue;
						}
						int n_contig_no=(*nc_ptr)->kmer_info.contig_no;
						if(n_contig_no==0)
						{
							edge_ptr=edge_ptr->nxt_edge;
							continue;
						}
						//if((*nc_ptr)->kmer_info.flip^flip_nc==0)
						bool t_flip=(*nc_ptr)->kmer_info.flip;


						char edge_cstr[200];
						uint64_t edgebits=edge_ptr->edge;
						//edgebits=get_rev_comp_seq(edgebits,edge_ptr->len);
						int edge_len= edge_ptr->len;
						edge_len++;
						string edge_s = bitsarr2str(&edgebits,edge_len,edge_cstr,1);
						

						if((t_flip^flip_nc)==1)
						{

							if(contigs_info->contig_adjacency_right[contig_no].count(-n_contig_no)==0)
							{
								contigs_info->contig_adjacency_right[contig_no][-n_contig_no].cov=edge_ptr->edge_cov;
								contigs_info->contig_adjacency_right[contig_no][-n_contig_no].bridge=edge_s;

							}
							else
							{
							//	contigs_info->contig_adjacency_right[contig_no][-n_contig_no].cov++;
							}

						}
						else
						{
							if(contigs_info->contig_adjacency_right[contig_no].count(n_contig_no)==0)
							{
								contigs_info->contig_adjacency_right[contig_no][n_contig_no].cov=edge_ptr->edge_cov;
								contigs_info->contig_adjacency_right[contig_no][n_contig_no].bridge=edge_s;

							}
							else
							{
							//	contigs_info->contig_adjacency_right[contig_no][n_contig_no].cov++;
							}

						}




						edge_ptr=edge_ptr->nxt_edge;

					}

				}
				else
				{
					//search left, app to adj right

					edge_node* edge_ptr=(*ptr)->kmer_info.left;
					while(edge_ptr!=NULL)
					{
						kmer_t2 t_kmer=seq,f_kmer;
						uint64_t t=0;
						bool flip_nc=0;

						for(int j=0;j<=edge_ptr->len;++j)
						{
							int left=(int)((edge_ptr->edge)>>2*j);



							switch(left&0x3)
							{
								case 0:
									R_shift_NB(t_kmer.kmer,2,2);


									break;
								case 1:
									R_shift_NB(t_kmer.kmer,2,2);
									t=((uint64_t)1)<<((K_size-1-32)*2);
									t_kmer.kmer[0]|=t;

									break;
								case 2:
									R_shift_NB(t_kmer.kmer,2,2);
									t=((uint64_t)2)<<((K_size-1-32)*2);
									t_kmer.kmer[0]|=t;

									break;
								case 3:
									R_shift_NB(t_kmer.kmer,2,2);
									t=((uint64_t)3)<<((K_size-1-32)*2);
									t_kmer.kmer[0]|=t;

									break;


							}


						}

						flip_nc=0;
						f_kmer=t_kmer;
						get_rev_comp_seq_arr(f_kmer.kmer,K_size,2);
						if(uint64_t_cmp(t_kmer.kmer,f_kmer.kmer,2)>0)
						{
							t_kmer=f_kmer;
							flip_nc=!flip_nc;
						}

						hv=MurmurHash64A(&t_kmer,sizeof(t_kmer),0);

						hash_idx=(size_t) (hv%ht2->ht_sz);


						struct bucket2 ** nc_ptr;

						nc_ptr= &(ht2->store_pos[hash_idx]);

						found=look_up_in_a_list2(&t_kmer,&nc_ptr);


						if(found==0)
						{
							edge_ptr=edge_ptr->nxt_edge;
							//cout<<"WarningR."<<endl;
							continue;
						}
						int n_contig_no=(*nc_ptr)->kmer_info.contig_no;
						if(n_contig_no==0)
						{
							edge_ptr=edge_ptr->nxt_edge;
							continue;
						}
					//	if((*nc_ptr)->kmer_info.flip^flip_nc)
						bool t_flip=(*nc_ptr)->kmer_info.flip;

					
						char edge_cstr[200];
						uint64_t edgebits=edge_ptr->edge;

						int edge_len= edge_ptr->len;
						edge_len++;
						edgebits=get_rev_comp_seq(edgebits,edge_len);
						string edge_s = bitsarr2str(&edgebits,edge_len,edge_cstr,1);

						

						if((t_flip^flip_nc)==0)
						{

							if(contigs_info->contig_adjacency_right[contig_no].count(-n_contig_no)==0)
							{
								contigs_info->contig_adjacency_right[contig_no][-n_contig_no].cov=edge_ptr->edge_cov;
								contigs_info->contig_adjacency_right[contig_no][-n_contig_no].bridge=edge_s;

							}
							else
							{
							//	contigs_info->contig_adjacency_right[contig_no][-n_contig_no].cov++;
							}

						}
						else
						{
							if(contigs_info->contig_adjacency_right[contig_no].count(n_contig_no)==0)
							{
								contigs_info->contig_adjacency_right[contig_no][n_contig_no].cov=edge_ptr->edge_cov;
								contigs_info->contig_adjacency_right[contig_no][n_contig_no].bridge=edge_s;

							}
							else
							{
								//contigs_info->contig_adjacency_right[contig_no][n_contig_no].cov++;
							}

						}




						edge_ptr=edge_ptr->nxt_edge;

					}
				}
			}






		}



	}

	out_ctg_graph<<contigs_info->total_contigs<<endl;
	for (int i=1;i<=contigs_info->total_contigs;++i)
	{
		out_ctg_graph<<i<<endl;
		out_ctg_graph<<contigs_info->contig_adjacency_left[i].size()<<endl;
		if(contigs_info->contig_adjacency_left[i].size()>0)
		{
		
			map<int,struct adjacent_contig_info>::iterator map_beg=contigs_info->contig_adjacency_left[i].begin();
			map<int,struct adjacent_contig_info>::iterator map_end=contigs_info->contig_adjacency_left[i].end();
			map<int,struct adjacent_contig_info>::iterator map_it;
			for(map_it=map_beg;map_it!=map_end;++map_it)
			{
				out_ctg_graph << map_it->first << " " << K_size - map_it->second.bridge.size() << ", ";
			}
			out_ctg_graph<<endl;
		}

		out_ctg_graph<<contigs_info->contig_adjacency_right[i].size()<<endl;
		if(contigs_info->contig_adjacency_right[i].size()>0)
		{
		
			map<int,struct adjacent_contig_info>::iterator map_beg=contigs_info->contig_adjacency_right[i].begin();
			map<int,struct adjacent_contig_info>::iterator map_end=contigs_info->contig_adjacency_right[i].end();
			map<int,struct adjacent_contig_info>::iterator map_it;
			for(map_it=map_beg;map_it!=map_end;++map_it)
			{
				out_ctg_graph << map_it->first << " " << K_size - map_it->second.bridge.size() << ", ";
			}
			out_ctg_graph<<endl;
		}

	}


}



void BuildContigAdjacency3(hashtable3 *ht, struct contigs_info *contigs_info,int K_size, string ContigFilename)
{
	ifstream in_ctg(ContigFilename.c_str());
	ofstream out_ctg_graph("ContigGraph.txt");
	string tag,contig_s;
	int contig_no=0;
	contigs_info->contig_adjacency_left.resize(contigs_info->total_contigs+1);
	contigs_info->contig_adjacency_right.resize(contigs_info->total_contigs+1);

	while(getline(in_ctg,tag))
	{
		contig_no++;	
		getline(in_ctg,contig_s);
		kmer_t3 seq,f_seq;
		uint64_t hv;
		size_t hash_idx;
		bool flip_c=0,found=0;
		struct bucket3 **ptr;
		int gMax=100;
		if((contig_s.size()-K_size)<100)
		{
			gMax=contig_s.size()-K_size;
		}
		for(int j=0;j<=gMax;++j)
		{

			string K_mer=contig_s.substr(j,K_size);
			str2bitsarr(K_mer.c_str(),K_size,seq.kmer,3);
			f_seq=seq;
			get_rev_comp_seq_arr(f_seq.kmer,K_size,3);
			flip_c=0;
			if(uint64_t_cmp(seq.kmer,f_seq.kmer,3)>0)
			{
				flip_c=1;
				seq=f_seq;
			}

			hv=MurmurHash64A(&seq,sizeof(seq),0);
			hash_idx=(size_t) (hv%(ht->ht_sz));

			ptr= &(ht->store_pos[hash_idx]);
			found=look_up_in_a_list3(&seq,&ptr);
			if(found==0)
			{continue;}
			else
			{break;}
		}

		if(found==0)
		{

		//	cout<<"searching error."<<endl;
			//return ;
		}

		if(found)
		{
			if((*ptr)->kmer_info.flip==0)
			{
				//search left, append to left
				edge_node* edge_ptr=(*ptr)->kmer_info.left;
			
				while(edge_ptr!=NULL)
				{
					kmer_t3 t_kmer=seq,f_kmer;
					uint64_t t=0;
					bool flip_nc=0;

					for(int j=0;j<=edge_ptr->len;++j)
					{
						int left=(int)((edge_ptr->edge)>>2*j);

						switch(left&0x3)
						{
							case 0:
								R_shift_NB(t_kmer.kmer,2,3);


								break;
							case 1:
								R_shift_NB(t_kmer.kmer,2,3);
								t=((uint64_t)1)<<((K_size-1-64)*2);
								t_kmer.kmer[0]|=t;

								break;
							case 2:
								R_shift_NB(t_kmer.kmer,2,3);
								t=((uint64_t)2)<<((K_size-1-64)*2);
								t_kmer.kmer[0]|=t;

								break;
							case 3:
								R_shift_NB(t_kmer.kmer,2,3);
								t=((uint64_t)3)<<((K_size-1-64)*2);
								t_kmer.kmer[0]|=t;

								break;

						}


					}

					flip_nc=0;
					f_kmer=t_kmer;
					get_rev_comp_seq_arr(f_kmer.kmer,K_size,3);
					if(uint64_t_cmp(t_kmer.kmer,f_kmer.kmer,3)>0)
					{
						t_kmer=f_kmer;
						flip_nc=!flip_nc;
					}

					hv=MurmurHash64A(&t_kmer,sizeof(t_kmer),0);

					hash_idx=(size_t) (hv%ht->ht_sz);


					struct bucket3 ** nc_ptr;

					nc_ptr= &(ht->store_pos[hash_idx]);

					found=look_up_in_a_list3(&t_kmer,&nc_ptr);


					if(found==0)
					{
						edge_ptr=edge_ptr->nxt_edge;
						//cout<<"WarningL."<<endl;
						continue;
					}
					int n_contig_no=(*nc_ptr)->kmer_info.contig_no;
					if(n_contig_no==0)
					{
						edge_ptr=edge_ptr->nxt_edge;
						continue;
					}
					//if((*nc_ptr)->kmer_info.flip^flip_nc)
					bool t_flip=(*nc_ptr)->kmer_info.flip;

					char edge_cstr[200];
					uint64_t edgebits=edge_ptr->edge;
					int edge_len= edge_ptr->len;
					edge_len++;
					//edgebits=get_rev_comp_seq(edgebits,edge_len);
					string edge_s = bitsarr2str(&edgebits,edge_len,edge_cstr,1);

					//string edge_s = bitsarr2str(&edgebits,(int)(edge_ptr->len+1),edge_cstr,1);


					if((t_flip^flip_nc)==1)
					{

						if(contigs_info->contig_adjacency_left[contig_no].count(-n_contig_no)==0)
						{
							contigs_info->contig_adjacency_left[contig_no][-n_contig_no].cov=edge_ptr->edge_cov;
							contigs_info->contig_adjacency_left[contig_no][-n_contig_no].bridge=edge_s;
					
						}
						else
						{
						//	contigs_info->contig_adjacency_left[contig_no][-n_contig_no].cov++;
						}

					}
					else
					{
						if(contigs_info->contig_adjacency_left[contig_no].count(n_contig_no)==0)
						{
							contigs_info->contig_adjacency_left[contig_no][n_contig_no].cov=edge_ptr->edge_cov;
							contigs_info->contig_adjacency_left[contig_no][n_contig_no].bridge=edge_s;

						}
						else
						{
							//contigs_info->contig_adjacency_left[contig_no][n_contig_no].cov++;
						}

					}


					edge_ptr=edge_ptr->nxt_edge;

				}

			}
			else
			{
				//search right, app to adj left


				edge_node* edge_ptr=(*ptr)->kmer_info.right;
				while(edge_ptr!=NULL)
				{
					kmer_t3 t_kmer=seq,f_kmer;
					uint64_t t=0;
					bool flip_nc=0;

					for(int j=edge_ptr->len;j>=0;--j)
					{
						int right=(int)((edge_ptr->edge)>>2*j);





						switch(right&0x3)
						{
							case 0:
								t=3;
								t<<=(K_size-1-64)*2;
								t_kmer.kmer[0]&=(~t);
								L_shift_NB(t_kmer.kmer,2,3);//



								break;
							case 1:
								t=3;
								t<<=(K_size-1-64)*2;
								t_kmer.kmer[0]&=(~t);
								L_shift_NB(t_kmer.kmer,2,3);//
								t=1;

								t_kmer.kmer[2]|=t;

								break;
							case 2:
								t=3;
								t<<=(K_size-1-64)*2;
								t_kmer.kmer[0]&=(~t);
								L_shift_NB(t_kmer.kmer,2,3);//
								t=2;

								t_kmer.kmer[2]|=t;

								break;
							case 3:
								t=3;
								t<<=(K_size-1-64)*2;
								t_kmer.kmer[0]&=(~t);
								L_shift_NB(t_kmer.kmer,2,3);//
								t=3;

								t_kmer.kmer[2]|=t;

								break;

						}

					}

					flip_nc=0;

					f_kmer=t_kmer;
					get_rev_comp_seq_arr(f_kmer.kmer,K_size,3);
					if(uint64_t_cmp(t_kmer.kmer,f_kmer.kmer,3)>0)
					{
						t_kmer=f_kmer;
						flip_nc=!flip_nc;
					}

					hv=MurmurHash64A(&t_kmer,sizeof(t_kmer),0);

					hash_idx=(size_t) (hv%ht->ht_sz);


					struct bucket3 ** nc_ptr;

					nc_ptr= &(ht->store_pos[hash_idx]);

					found=look_up_in_a_list3(&t_kmer,&nc_ptr);


					if(found==0)
					{
						edge_ptr=edge_ptr->nxt_edge;
						//cout<<"WarningL."<<endl;
						continue;
					}
					int n_contig_no=(*nc_ptr)->kmer_info.contig_no;
					if(n_contig_no==0)
					{
						edge_ptr=edge_ptr->nxt_edge;
						continue;
					}
					//if((*nc_ptr)->kmer_info.flip^flip_nc==0)
					bool t_flip=(*nc_ptr)->kmer_info.flip;
					char edge_cstr[200];
					uint64_t edgebits=edge_ptr->edge;
					int edge_len= edge_ptr->len;
					edge_len++;
					edgebits=get_rev_comp_seq(edgebits,edge_len);
					string edge_s = bitsarr2str(&edgebits,edge_len,edge_cstr,1);
					//edgebits=get_rev_comp_seq(edgebits,(int)(edge_ptr->len+1));
					//string edge_s = bitsarr2str(&edgebits,(int)(edge_ptr->len+1),edge_cstr,1);

					if((t_flip^flip_nc)==0)
					{

						if(contigs_info->contig_adjacency_left[contig_no].count(-n_contig_no)==0)
						{
							contigs_info->contig_adjacency_left[contig_no][-n_contig_no].cov=edge_ptr->edge_cov;;
							contigs_info->contig_adjacency_left[contig_no][-n_contig_no].bridge=edge_s;

						}
						else
						{
							//contigs_info->contig_adjacency_left[contig_no][-n_contig_no].cov++;
						}

					}
					else
					{
						if(contigs_info->contig_adjacency_left[contig_no].count(n_contig_no)==0)
						{
							contigs_info->contig_adjacency_left[contig_no][n_contig_no].cov=edge_ptr->edge_cov;;
							contigs_info->contig_adjacency_left[contig_no][n_contig_no].bridge=edge_s;
					
						}
						else
						{
						//	contigs_info->contig_adjacency_left[contig_no][n_contig_no].cov++;
						}

					}




					edge_ptr=edge_ptr->nxt_edge;

				}

			}
		}

		//contig endpoint

		int gMin=contig_s.size()-K_size-100;
		if(gMin<0)
		{gMin=0;}

		for(int j=contig_s.size()-K_size;j>=gMin;--j)
		{
			string K_mer=contig_s.substr(j,K_size);
			str2bitsarr(K_mer.c_str(),K_size,seq.kmer,3);

			f_seq=seq;
			get_rev_comp_seq_arr(f_seq.kmer,K_size,3);
			flip_c=0;
			if(uint64_t_cmp(seq.kmer,f_seq.kmer,3)>0)
			{
				flip_c=1;
				seq=f_seq;
			}

			hv=MurmurHash64A(&seq,sizeof(seq),0);
			hash_idx=(size_t) (hv%(ht->ht_sz));

			ptr= &(ht->store_pos[hash_idx]);
			found=look_up_in_a_list3(&seq,&ptr);
			if(found==0)
			{
				continue;
			}
			else
			{
				break;
			}
		}

		if(found==0)
		{

//			cout<<"searching error."<<endl;
	//		return ;
		}

		if(found)
		{
			if((*ptr)->kmer_info.flip==0)
			{
				//search right, app to adj right

				edge_node* edge_ptr=(*ptr)->kmer_info.right;
				while(edge_ptr!=NULL)
				{
					kmer_t3 t_kmer=seq,f_kmer;
					uint64_t t=0;
					bool flip_nc=0;

					for(int j=edge_ptr->len;j>=0;--j)
					{
						int right=(int)((edge_ptr->edge)>>2*j);





						switch(right&0x3)
						{
							case 0:
								t=3;
								t<<=(K_size-1-64)*2;
								t_kmer.kmer[0]&=(~t);
								L_shift_NB(t_kmer.kmer,2,3);//



								break;
							case 1:
								t=3;
								t<<=(K_size-1-64)*2;
								t_kmer.kmer[0]&=(~t);
								L_shift_NB(t_kmer.kmer,2,3);//
								t=1;

								t_kmer.kmer[2]|=t;

								break;
							case 2:
								t=3;
								t<<=(K_size-1-64)*2;
								t_kmer.kmer[0]&=(~t);
								L_shift_NB(t_kmer.kmer,2,3);//
								t=2;

								t_kmer.kmer[2]|=t;

								break;
							case 3:
								t=3;
								t<<=(K_size-1-64)*2;
								t_kmer.kmer[0]&=(~t);
								L_shift_NB(t_kmer.kmer,2,3);//
								t=3;

								t_kmer.kmer[2]|=t;

								break;

						}

					}

					flip_nc=0;

					f_kmer=t_kmer;
					get_rev_comp_seq_arr(f_kmer.kmer,K_size,3);
					if(uint64_t_cmp(t_kmer.kmer,f_kmer.kmer,3)>0)
					{
						t_kmer=f_kmer;
						flip_nc=!flip_nc;
					}

					hv=MurmurHash64A(&t_kmer,sizeof(t_kmer),0);

					hash_idx=(size_t) (hv%ht->ht_sz);


					struct bucket3 ** nc_ptr;

					nc_ptr= &(ht->store_pos[hash_idx]);

					found=look_up_in_a_list3(&t_kmer,&nc_ptr);


					if(found==0)
					{
						edge_ptr=edge_ptr->nxt_edge;
						//cout<<"WarningL."<<endl;
						continue;
					}
					int n_contig_no=(*nc_ptr)->kmer_info.contig_no;
					if(n_contig_no==0)
					{
						edge_ptr=edge_ptr->nxt_edge;
						continue;
					}
					//if((*nc_ptr)->kmer_info.flip^flip_nc==0)
					bool t_flip=(*nc_ptr)->kmer_info.flip;


					char edge_cstr[200];
					uint64_t edgebits=edge_ptr->edge;
					//edgebits=get_rev_comp_seq(edgebits,edge_ptr->len);
					int edge_len= edge_ptr->len;
					edge_len++;
					//edgebits=get_rev_comp_seq(edgebits,edge_len);
					string edge_s = bitsarr2str(&edgebits,edge_len,edge_cstr,1);

					//string edge_s = bitsarr2str(&edgebits,(int)(edge_ptr->len+1),edge_cstr,1);


					if((t_flip^flip_nc)==1)
					{

						if(contigs_info->contig_adjacency_right[contig_no].count(-n_contig_no)==0)
						{
							contigs_info->contig_adjacency_right[contig_no][-n_contig_no].cov=edge_ptr->edge_cov;
							contigs_info->contig_adjacency_right[contig_no][-n_contig_no].bridge=edge_s;

						}
						else
						{
						//	contigs_info->contig_adjacency_right[contig_no][-n_contig_no].cov++;
						}

					}
					else
					{
						if(contigs_info->contig_adjacency_right[contig_no].count(n_contig_no)==0)
						{
							contigs_info->contig_adjacency_right[contig_no][n_contig_no].cov=edge_ptr->edge_cov;
							contigs_info->contig_adjacency_right[contig_no][n_contig_no].bridge=edge_s;

						}
						else
						{
						//	contigs_info->contig_adjacency_right[contig_no][n_contig_no].cov++;
						}

					}




					edge_ptr=edge_ptr->nxt_edge;

				}

			}
			else
			{
				//search left, app to adj right

				edge_node* edge_ptr=(*ptr)->kmer_info.left;
				while(edge_ptr!=NULL)
				{
					kmer_t3 t_kmer=seq,f_kmer;
					uint64_t t=0;
					bool flip_nc=0;

					for(int j=0;j<=edge_ptr->len;++j)
					{
						int left=(int)((edge_ptr->edge)>>2*j);



						switch(left&0x3)
						{
							case 0:
								R_shift_NB(t_kmer.kmer,2,3);


								break;
							case 1:
								R_shift_NB(t_kmer.kmer,2,3);
								t=((uint64_t)1)<<((K_size-1-64)*2);
								t_kmer.kmer[0]|=t;

								break;
							case 2:
								R_shift_NB(t_kmer.kmer,2,3);
								t=((uint64_t)2)<<((K_size-1-64)*2);
								t_kmer.kmer[0]|=t;

								break;
							case 3:
								R_shift_NB(t_kmer.kmer,2,3);
								t=((uint64_t)3)<<((K_size-1-64)*2);
								t_kmer.kmer[0]|=t;

								break;


						}


					}

					flip_nc=0;
					f_kmer=t_kmer;
					get_rev_comp_seq_arr(f_kmer.kmer,K_size,3);
					if(uint64_t_cmp(t_kmer.kmer,f_kmer.kmer,3)>0)
					{
						t_kmer=f_kmer;
						flip_nc=!flip_nc;
					}

					hv=MurmurHash64A(&t_kmer,sizeof(t_kmer),0);

					hash_idx=(size_t) (hv%ht->ht_sz);


					struct bucket3 ** nc_ptr;

					nc_ptr= &(ht->store_pos[hash_idx]);

					found=look_up_in_a_list3(&t_kmer,&nc_ptr);


					if(found==0)
					{
						edge_ptr=edge_ptr->nxt_edge;
						//cout<<"WarningR."<<endl;
						continue;
					}
					int n_contig_no=(*nc_ptr)->kmer_info.contig_no;
					if(n_contig_no==0)
					{
						edge_ptr=edge_ptr->nxt_edge;
						continue;
					}
				//	if((*nc_ptr)->kmer_info.flip^flip_nc)
					bool t_flip=(*nc_ptr)->kmer_info.flip;

					
					char edge_cstr[200];
					uint64_t edgebits=edge_ptr->edge;
					int edge_len= edge_ptr->len;
					edge_len++;
					edgebits=get_rev_comp_seq(edgebits,edge_len);
					string edge_s = bitsarr2str(&edgebits,edge_len,edge_cstr,1);
					//edgebits=get_rev_comp_seq(edgebits,(int)(edge_ptr->len+1));
					//string edge_s = bitsarr2str(&edgebits,(int)(edge_ptr->len+1),edge_cstr,1);


					if((t_flip^flip_nc)==0)
					{

						if(contigs_info->contig_adjacency_right[contig_no].count(-n_contig_no)==0)
						{
							contigs_info->contig_adjacency_right[contig_no][-n_contig_no].cov=edge_ptr->edge_cov;
							contigs_info->contig_adjacency_right[contig_no][-n_contig_no].bridge=edge_s;

						}
						else
						{
						//	contigs_info->contig_adjacency_right[contig_no][-n_contig_no].cov++;
						}

					}
					else
					{
						if(contigs_info->contig_adjacency_right[contig_no].count(n_contig_no)==0)
						{
							contigs_info->contig_adjacency_right[contig_no][n_contig_no].cov=edge_ptr->edge_cov;
							contigs_info->contig_adjacency_right[contig_no][n_contig_no].bridge=edge_s;

						}
						else
						{
							//contigs_info->contig_adjacency_right[contig_no][n_contig_no].cov++;
						}

					}




					edge_ptr=edge_ptr->nxt_edge;

				}
			}


		}




	}

	out_ctg_graph<<contigs_info->total_contigs<<endl;
	for (int i=1;i<=contigs_info->total_contigs;++i)
	{
		out_ctg_graph<<i<<endl;
		out_ctg_graph<<contigs_info->contig_adjacency_left[i].size()<<endl;
		if(contigs_info->contig_adjacency_left[i].size()>0)
		{
		
			map<int,struct adjacent_contig_info>::iterator map_beg=contigs_info->contig_adjacency_left[i].begin();
			map<int,struct adjacent_contig_info>::iterator map_end=contigs_info->contig_adjacency_left[i].end();
			map<int,struct adjacent_contig_info>::iterator map_it;
			for(map_it=map_beg;map_it!=map_end;++map_it)
			{
				out_ctg_graph << map_it->first << " " << K_size - map_it->second.bridge.size() << ", ";
			}
			out_ctg_graph<<endl;
		}

		out_ctg_graph<<contigs_info->contig_adjacency_right[i].size()<<endl;
		if(contigs_info->contig_adjacency_right[i].size()>0)
		{
		
			map<int,struct adjacent_contig_info>::iterator map_beg=contigs_info->contig_adjacency_right[i].begin();
			map<int,struct adjacent_contig_info>::iterator map_end=contigs_info->contig_adjacency_right[i].end();
			map<int,struct adjacent_contig_info>::iterator map_it;
			for(map_it=map_beg;map_it!=map_end;++map_it)
			{
				out_ctg_graph << map_it->first << " " << K_size - map_it->second.bridge.size() << ", ";
			}
			out_ctg_graph<<endl;
		}

	}

}

void BuildContigAdjacency4(hashtable4 *ht, struct contigs_info *contigs_info,int K_size, string ContigFilename)
{
	ifstream in_ctg(ContigFilename.c_str());
	ofstream out_ctg_graph("ContigGraph.txt");
	string tag,contig_s;
	int contig_no=0;
	contigs_info->contig_adjacency_left.resize(contigs_info->total_contigs+1);
	contigs_info->contig_adjacency_right.resize(contigs_info->total_contigs+1);

	while(getline(in_ctg,tag))
	{
		contig_no++;	
		getline(in_ctg,contig_s);
		kmer_t4 seq,f_seq;
		uint64_t hv;
		size_t hash_idx;
		bool flip_c=0,found=0;
		struct bucket4 **ptr;
		int gMax=100;
		if((contig_s.size()-K_size)<100)
		{
			gMax=contig_s.size()-K_size;
		}
		for(int j=0;j<=gMax;++j)
		{

			string K_mer=contig_s.substr(j,K_size);
			str2bitsarr(K_mer.c_str(),K_size,seq.kmer,4);
			f_seq=seq;
			get_rev_comp_seq_arr(f_seq.kmer,K_size,4);
			flip_c=0;
			if(uint64_t_cmp(seq.kmer,f_seq.kmer,4)>0)
			{
				flip_c=1;
				seq=f_seq;
			}

			hv=MurmurHash64A(&seq,sizeof(seq),0);
			hash_idx=(size_t) (hv%(ht->ht_sz));

			ptr= &(ht->store_pos[hash_idx]);
			found=look_up_in_a_list4(&seq,&ptr);
			if(found==0)
			{continue;}
			else
			{break;}
		}

		if(found==0)
		{

		//	cout<<"searching error."<<endl;
		//	return ;
		}
		if(found)
		{
			if((*ptr)->kmer_info.flip==0)
			{
				//search left, append to left
				edge_node* edge_ptr=(*ptr)->kmer_info.left;
			
				while(edge_ptr!=NULL)
				{
					kmer_t4 t_kmer=seq,f_kmer;
					uint64_t t=0;
					bool flip_nc=0;

					for(int j=0;j<=edge_ptr->len;++j)
					{
						int left=(int)((edge_ptr->edge)>>2*j);



						switch(left&0x3)
						{
							case 0:
								R_shift_NB(t_kmer.kmer,2,4);


								break;
							case 1:
								R_shift_NB(t_kmer.kmer,2,4);
								t=((uint64_t)1)<<((K_size-1-96)*2);
								t_kmer.kmer[0]|=t;

								break;
							case 2:
								R_shift_NB(t_kmer.kmer,2,4);
								t=((uint64_t)2)<<((K_size-1-96)*2);
								t_kmer.kmer[0]|=t;

								break;
							case 3:
								R_shift_NB(t_kmer.kmer,2,4);
								t=((uint64_t)3)<<((K_size-1-96)*2);
								t_kmer.kmer[0]|=t;

								break;


						}


					}

					flip_nc=0;
					f_kmer=t_kmer;
					get_rev_comp_seq_arr(f_kmer.kmer,K_size,4);
					if(uint64_t_cmp(t_kmer.kmer,f_kmer.kmer,4)>0)
					{
						t_kmer=f_kmer;
						flip_nc=!flip_nc;
					}

					hv=MurmurHash64A(&t_kmer,sizeof(t_kmer),0);

					hash_idx=(size_t) (hv%ht->ht_sz);


					struct bucket4 ** nc_ptr;

					nc_ptr= &(ht->store_pos[hash_idx]);

					found=look_up_in_a_list4(&t_kmer,&nc_ptr);


					if(found==0)
					{
						edge_ptr=edge_ptr->nxt_edge;
						//cout<<"WarningL."<<endl;
						continue;
					}
					int n_contig_no=(*nc_ptr)->kmer_info.contig_no;
					if(n_contig_no==0)
					{
						edge_ptr=edge_ptr->nxt_edge;
						continue;
					}
					//if((*nc_ptr)->kmer_info.flip^flip_nc)
					bool t_flip=(*nc_ptr)->kmer_info.flip;

					char edge_cstr[200];
					uint64_t edgebits=edge_ptr->edge;
					int edge_len= edge_ptr->len;
					edge_len++;
					//edgebits=get_rev_comp_seq(edgebits,edge_len);
					string edge_s = bitsarr2str(&edgebits,edge_len,edge_cstr,1);


					//string edge_s = bitsarr2str(&edgebits,(int)(edge_ptr->len+1),edge_cstr,1);


					if((t_flip^flip_nc)==1)
					{

						if(contigs_info->contig_adjacency_left[contig_no].count(-n_contig_no)==0)
						{
							contigs_info->contig_adjacency_left[contig_no][-n_contig_no].cov=edge_ptr->edge_cov;
							contigs_info->contig_adjacency_left[contig_no][-n_contig_no].bridge=edge_s;
					
						}
						else
						{
						//	contigs_info->contig_adjacency_left[contig_no][-n_contig_no].cov++;
						}

					}
					else
					{
						if(contigs_info->contig_adjacency_left[contig_no].count(n_contig_no)==0)
						{
							contigs_info->contig_adjacency_left[contig_no][n_contig_no].cov=edge_ptr->edge_cov;
							contigs_info->contig_adjacency_left[contig_no][n_contig_no].bridge=edge_s;

						}
						else
						{
							//contigs_info->contig_adjacency_left[contig_no][n_contig_no].cov++;
						}

					}


					edge_ptr=edge_ptr->nxt_edge;

				}

			}
			else
			{
				//search right, app to adj left


				edge_node* edge_ptr=(*ptr)->kmer_info.right;
				while(edge_ptr!=NULL)
				{
					kmer_t4 t_kmer=seq,f_kmer;
					uint64_t t=0;
					bool flip_nc=0;

					for(int j=edge_ptr->len;j>=0;--j)
					{
						int right=(int)((edge_ptr->edge)>>2*j);




						switch(right&0x3)
						{
							case 0:
								t=3;
								t<<=(K_size-1-96)*2;
								t_kmer.kmer[0]&=(~t);
								L_shift_NB(t_kmer.kmer,2,4);//



								break;
							case 1:
								t=3;
								t<<=(K_size-1-96)*2;
								t_kmer.kmer[0]&=(~t);
								L_shift_NB(t_kmer.kmer,2,4);//
								t=1;

								t_kmer.kmer[3]|=t;

								break;
							case 2:
								t=3;
								t<<=(K_size-1-96)*2;
								t_kmer.kmer[0]&=(~t);
								L_shift_NB(t_kmer.kmer,2,4);//
								t=2;

								t_kmer.kmer[3]|=t;

								break;
							case 3:
								t=3;
								t<<=(K_size-1-96)*2;
								t_kmer.kmer[0]&=(~t);
								L_shift_NB(t_kmer.kmer,2,4);//
								t=3;

								t_kmer.kmer[3]|=t;

								break;

						}

					}

					flip_nc=0;

					f_kmer=t_kmer;
					get_rev_comp_seq_arr(f_kmer.kmer,K_size,4);
					if(uint64_t_cmp(t_kmer.kmer,f_kmer.kmer,4)>0)
					{
						t_kmer=f_kmer;
						flip_nc=!flip_nc;
					}

					hv=MurmurHash64A(&t_kmer,sizeof(t_kmer),0);

					hash_idx=(size_t) (hv%ht->ht_sz);


					struct bucket4 ** nc_ptr;

					nc_ptr= &(ht->store_pos[hash_idx]);

					found=look_up_in_a_list4(&t_kmer,&nc_ptr);


					if(found==0)
					{
						edge_ptr=edge_ptr->nxt_edge;
						//cout<<"WarningL."<<endl;
						continue;
					}
					int n_contig_no=(*nc_ptr)->kmer_info.contig_no;
					if(n_contig_no==0)
					{
						edge_ptr=edge_ptr->nxt_edge;
						continue;
					}
					//if((*nc_ptr)->kmer_info.flip^flip_nc==0)
					bool t_flip=(*nc_ptr)->kmer_info.flip;
					char edge_cstr[200];
					uint64_t edgebits=edge_ptr->edge;
					int edge_len= edge_ptr->len;
					edge_len++;
					edgebits=get_rev_comp_seq(edgebits,edge_len);
					string edge_s = bitsarr2str(&edgebits,edge_len,edge_cstr,1);

					//edgebits=get_rev_comp_seq(edgebits,(int)(edge_ptr->len+1));
					//string edge_s = bitsarr2str(&edgebits,(int)(edge_ptr->len+1),edge_cstr,1);

					if((t_flip^flip_nc)==0)
					{

						if(contigs_info->contig_adjacency_left[contig_no].count(-n_contig_no)==0)
						{
							contigs_info->contig_adjacency_left[contig_no][-n_contig_no].cov=edge_ptr->edge_cov;;
							contigs_info->contig_adjacency_left[contig_no][-n_contig_no].bridge=edge_s;

						}
						else
						{
							//contigs_info->contig_adjacency_left[contig_no][-n_contig_no].cov++;
						}

					}
					else
					{
						if(contigs_info->contig_adjacency_left[contig_no].count(n_contig_no)==0)
						{
							contigs_info->contig_adjacency_left[contig_no][n_contig_no].cov=edge_ptr->edge_cov;;
							contigs_info->contig_adjacency_left[contig_no][n_contig_no].bridge=edge_s;
					
						}
						else
						{
						//	contigs_info->contig_adjacency_left[contig_no][n_contig_no].cov++;
						}

					}




					edge_ptr=edge_ptr->nxt_edge;

				}

			}
		}

		//contig endpoint

		int gMin=contig_s.size()-K_size-100;
		if(gMin<0)
		{gMin=0;}

		for(int j=contig_s.size()-K_size;j>=gMin;--j)
		{
			string K_mer=contig_s.substr(j,K_size);
			str2bitsarr(K_mer.c_str(),K_size,seq.kmer,4);

			f_seq=seq;
			get_rev_comp_seq_arr(f_seq.kmer,K_size,4);
			flip_c=0;
			if(uint64_t_cmp(seq.kmer,f_seq.kmer,4)>0)
			{
				flip_c=1;
				seq=f_seq;
			}

			hv=MurmurHash64A(&seq,sizeof(seq),0);
			hash_idx=(size_t) (hv%(ht->ht_sz));

			ptr= &(ht->store_pos[hash_idx]);
			found=look_up_in_a_list4(&seq,&ptr);
			if(found==0)
			{
				continue;
			}
			else
			{
				break;
			}
		}

		if(found==0)
		{

			//cout<<"searching error."<<endl;
			//return ;
		}

		if(found)
		{
			if((*ptr)->kmer_info.flip==0)
			{
				//search right, app to adj right

				edge_node* edge_ptr=(*ptr)->kmer_info.right;
				while(edge_ptr!=NULL)
				{
					kmer_t4 t_kmer=seq,f_kmer;
					uint64_t t=0;
					bool flip_nc=0;

					for(int j=edge_ptr->len;j>=0;--j)
					{
						int right=(int)((edge_ptr->edge)>>2*j);


						switch(right&0x3)
						{
							case 0:
								t=3;
								t<<=(K_size-1-96)*2;
								t_kmer.kmer[0]&=(~t);
								L_shift_NB(t_kmer.kmer,2,4);//



								break;
							case 1:
								t=3;
								t<<=(K_size-1-96)*2;
								t_kmer.kmer[0]&=(~t);
								L_shift_NB(t_kmer.kmer,2,4);//
								t=1;

								t_kmer.kmer[3]|=t;

								break;
							case 2:
								t=3;
								t<<=(K_size-1-96)*2;
								t_kmer.kmer[0]&=(~t);
								L_shift_NB(t_kmer.kmer,2,4);//
								t=2;

								t_kmer.kmer[3]|=t;

								break;
							case 3:
								t=3;
								t<<=(K_size-1-96)*2;
								t_kmer.kmer[0]&=(~t);
								L_shift_NB(t_kmer.kmer,2,4);//
								t=3;

								t_kmer.kmer[3]|=t;

								break;

						}

					}

					flip_nc=0;

					f_kmer=t_kmer;
					get_rev_comp_seq_arr(f_kmer.kmer,K_size,4);
					if(uint64_t_cmp(t_kmer.kmer,f_kmer.kmer,4)>0)
					{
						t_kmer=f_kmer;
						flip_nc=!flip_nc;
					}

					hv=MurmurHash64A(&t_kmer,sizeof(t_kmer),0);

					hash_idx=(size_t) (hv%ht->ht_sz);


					struct bucket4 ** nc_ptr;

					nc_ptr= &(ht->store_pos[hash_idx]);

					found=look_up_in_a_list4(&t_kmer,&nc_ptr);


					if(found==0)
					{
						edge_ptr=edge_ptr->nxt_edge;
						//cout<<"WarningL."<<endl;
						continue;
					}
					int n_contig_no=(*nc_ptr)->kmer_info.contig_no;

					if(n_contig_no==0)
					{
						edge_ptr=edge_ptr->nxt_edge;
						continue;
					}

					//if((*nc_ptr)->kmer_info.flip^flip_nc==0)
					bool t_flip=(*nc_ptr)->kmer_info.flip;


					char edge_cstr[200];
					uint64_t edgebits=edge_ptr->edge;
					//edgebits=get_rev_comp_seq(edgebits,edge_ptr->len);
					int edge_len= edge_ptr->len;
					edge_len++;
					//edgebits=get_rev_comp_seq(edgebits,edge_len);
					string edge_s = bitsarr2str(&edgebits,edge_len,edge_cstr,1);
					//string edge_s = bitsarr2str(&edgebits,(int)(edge_ptr->len+1),edge_cstr,1);


					if((t_flip^flip_nc)==1)
					{

						if(contigs_info->contig_adjacency_right[contig_no].count(-n_contig_no)==0)
						{
							contigs_info->contig_adjacency_right[contig_no][-n_contig_no].cov=edge_ptr->edge_cov;
							contigs_info->contig_adjacency_right[contig_no][-n_contig_no].bridge=edge_s;

						}
						else
						{
						//	contigs_info->contig_adjacency_right[contig_no][-n_contig_no].cov++;
						}

					}
					else
					{
						if(contigs_info->contig_adjacency_right[contig_no].count(n_contig_no)==0)
						{
							contigs_info->contig_adjacency_right[contig_no][n_contig_no].cov=edge_ptr->edge_cov;
							contigs_info->contig_adjacency_right[contig_no][n_contig_no].bridge=edge_s;

						}
						else
						{
						//	contigs_info->contig_adjacency_right[contig_no][n_contig_no].cov++;
						}

					}




					edge_ptr=edge_ptr->nxt_edge;

				}

			}
			else
			{
				//search left, app to adj right

				edge_node* edge_ptr=(*ptr)->kmer_info.left;
				while(edge_ptr!=NULL)
				{
					kmer_t4 t_kmer=seq,f_kmer;
					uint64_t t=0;
					bool flip_nc=0;

					for(int j=0;j<=edge_ptr->len;++j)
					{
						int left=(int)((edge_ptr->edge)>>2*j);



						switch(left&0x3)
						{
							case 0:
								R_shift_NB(t_kmer.kmer,2,4);


								break;
							case 1:
								R_shift_NB(t_kmer.kmer,2,4);
								t=((uint64_t)1)<<((K_size-1-96)*2);
								t_kmer.kmer[0]|=t;

								break;
							case 2:
								R_shift_NB(t_kmer.kmer,2,4);
								t=((uint64_t)2)<<((K_size-1-96)*2);
								t_kmer.kmer[0]|=t;

								break;
							case 3:
								R_shift_NB(t_kmer.kmer,2,4);
								t=((uint64_t)3)<<((K_size-1-96)*2);
								t_kmer.kmer[0]|=t;

								break;


						}


					}

					flip_nc=0;
					f_kmer=t_kmer;
					get_rev_comp_seq_arr(f_kmer.kmer,K_size,4);
					if(uint64_t_cmp(t_kmer.kmer,f_kmer.kmer,4)>0)
					{
						t_kmer=f_kmer;
						flip_nc=!flip_nc;
					}

					hv=MurmurHash64A(&t_kmer,sizeof(t_kmer),0);

					hash_idx=(size_t) (hv%ht->ht_sz);


					struct bucket4 ** nc_ptr;

					nc_ptr= &(ht->store_pos[hash_idx]);

					found=look_up_in_a_list4(&t_kmer,&nc_ptr);


					if(found==0)
					{
						edge_ptr=edge_ptr->nxt_edge;
						//cout<<"WarningR."<<endl;
						continue;
					}
					int n_contig_no=(*nc_ptr)->kmer_info.contig_no;
					if(n_contig_no==0)
					{
						edge_ptr=edge_ptr->nxt_edge;
						continue;
					}
				//	if((*nc_ptr)->kmer_info.flip^flip_nc)
					bool t_flip=(*nc_ptr)->kmer_info.flip;

					
					char edge_cstr[200];
					uint64_t edgebits=edge_ptr->edge;
					int edge_len= edge_ptr->len;
					edge_len++;
					edgebits=get_rev_comp_seq(edgebits,edge_len);
					string edge_s = bitsarr2str(&edgebits,edge_len,edge_cstr,1);
//					edgebits=get_rev_comp_seq(edgebits,(int)(edge_ptr->len+1));
	//				string edge_s = bitsarr2str(&edgebits,(int)(edge_ptr->len+1),edge_cstr,1);


					if((t_flip^flip_nc)==0)
					{

						if(contigs_info->contig_adjacency_right[contig_no].count(-n_contig_no)==0)
						{
							contigs_info->contig_adjacency_right[contig_no][-n_contig_no].cov=edge_ptr->edge_cov;
							contigs_info->contig_adjacency_right[contig_no][-n_contig_no].bridge=edge_s;

						}
						else
						{
						//	contigs_info->contig_adjacency_right[contig_no][-n_contig_no].cov++;
						}

					}
					else
					{
						if(contigs_info->contig_adjacency_right[contig_no].count(n_contig_no)==0)
						{
							contigs_info->contig_adjacency_right[contig_no][n_contig_no].cov=edge_ptr->edge_cov;
							contigs_info->contig_adjacency_right[contig_no][n_contig_no].bridge=edge_s;

						}
						else
						{
							//contigs_info->contig_adjacency_right[contig_no][n_contig_no].cov++;
						}

					}




					edge_ptr=edge_ptr->nxt_edge;

				}
			}


		}	




	}

	out_ctg_graph<<contigs_info->total_contigs<<endl;
	for (int i=1;i<=contigs_info->total_contigs;++i)
	{
		out_ctg_graph<<i<<endl;
		out_ctg_graph<<contigs_info->contig_adjacency_left[i].size()<<endl;
		if(contigs_info->contig_adjacency_left[i].size()>0)
		{
		
			map<int,struct adjacent_contig_info>::iterator map_beg=contigs_info->contig_adjacency_left[i].begin();
			map<int,struct adjacent_contig_info>::iterator map_end=contigs_info->contig_adjacency_left[i].end();
			map<int,struct adjacent_contig_info>::iterator map_it;
			for(map_it=map_beg;map_it!=map_end;++map_it)
			{
				out_ctg_graph<<map_it->first<<" "<<K_size-map_it->second.bridge.size()<<", ";
			}
			out_ctg_graph<<endl;
		}

		out_ctg_graph<<contigs_info->contig_adjacency_right[i].size()<<endl;
		if(contigs_info->contig_adjacency_right[i].size()>0)
		{
		
			map<int,struct adjacent_contig_info>::iterator map_beg=contigs_info->contig_adjacency_right[i].begin();
			map<int,struct adjacent_contig_info>::iterator map_end=contigs_info->contig_adjacency_right[i].end();
			map<int,struct adjacent_contig_info>::iterator map_it;
			for(map_it=map_beg;map_it!=map_end;++map_it)
			{
				out_ctg_graph<<map_it->first<<" "<<K_size-map_it->second.bridge.size()<<", ";
			}
			out_ctg_graph<<endl;
		}

	}

}



void BuildContigAdjacency0(hashtable0 *ht, struct contigs_info *contigs_info,int K_size, string ContigFilename)
{
	ifstream in_ctg(ContigFilename.c_str());
	ofstream out_ctg_graph("ContigGraph.txt");
	string tag,contig_s;

	int Kmer_arr_sz=K_size/32+1;
	int rem1=K_size%32;
	if(rem1==0)
	{Kmer_arr_sz--;}

	int contig_no=0;
	contigs_info->contig_adjacency_left.resize(contigs_info->total_contigs+1);
	contigs_info->contig_adjacency_right.resize(contigs_info->total_contigs+1);

	while(getline(in_ctg,tag))
	{
		contig_no++;	
		getline(in_ctg,contig_s);
		uint64_t seq[100],f_seq[100];
		uint64_t hv;
		size_t hash_idx;
		bool flip_c=0,found=0;
		struct bucket0 **ptr;
		int gMax=100;
		if((contig_s.size()-K_size)<100)
		{
			gMax=contig_s.size()-K_size;
		}
		for(int j=0;j<=gMax;++j)
		{

			string K_mer=contig_s.substr(j,K_size);
			str2bitsarr(K_mer.c_str(),K_size,seq,Kmer_arr_sz);
			
			memcpy(f_seq,seq,Kmer_arr_sz*sizeof(uint64_t));
			get_rev_comp_seq_arr(f_seq,K_size,Kmer_arr_sz);
			flip_c=0;
			if(uint64_t_cmp(seq,f_seq,Kmer_arr_sz)>0)
			{
				flip_c=1;
				memcpy(seq,f_seq,Kmer_arr_sz*sizeof(uint64_t));			
			}

			hv=MurmurHash64A(&seq,sizeof(uint64_t)*Kmer_arr_sz,0);
			hash_idx=(size_t) (hv%(ht->ht_sz));

			ptr= &(ht->store_pos[hash_idx]);
			
			found=look_up_in_a_list0(seq, &ptr ,Kmer_arr_sz);
			if(found==0)
			{continue;}
			else
			{break;}
		}

		if(found==0)
		{

		//	cout<<"searching error."<<endl;
			//return ;
		}

		if(found)
		{
			if((*ptr)->kmer_info.flip==0)
			{
				//search left, append to left
				edge_node* edge_ptr=(*ptr)->kmer_info.left;
			
				while(edge_ptr!=NULL)
				{
					uint64_t t_kmer[100],f_kmer[100];
					memcpy(t_kmer,seq,Kmer_arr_sz*sizeof(uint64_t));
					
					uint64_t t=0;
					bool flip_nc=0;
					for(int j=0;j<=edge_ptr->len;++j)
					{
						uint64_t left=(uint64_t)((edge_ptr->edge)>>2*j);

						left&=0x3;
						R_shift_NB(t_kmer,2,Kmer_arr_sz);
						uint64_t t=((uint64_t)1)<<(((K_size-1)%32)*2);
						t_kmer[0]|=t;

					}

					flip_nc=0;
					memcpy(f_kmer,t_kmer,Kmer_arr_sz*sizeof(uint64_t));
					get_rev_comp_seq_arr(f_kmer,K_size,Kmer_arr_sz);
					if(uint64_t_cmp(t_kmer,f_kmer,Kmer_arr_sz)>0)
					{
						memcpy(t_kmer,f_kmer,Kmer_arr_sz*sizeof(uint64_t));
						flip_nc=!flip_nc;
					}

					hv=MurmurHash64A(&t_kmer,sizeof(uint64_t)*Kmer_arr_sz,0);

					hash_idx=(size_t) (hv%ht->ht_sz);


					struct bucket0 ** nc_ptr;

					nc_ptr= &(ht->store_pos[hash_idx]);

					found=look_up_in_a_list0(t_kmer,&nc_ptr,Kmer_arr_sz);


					if(found==0)
					{
						edge_ptr=edge_ptr->nxt_edge;
						//cout<<"WarningL."<<endl;
						continue;
					}
					int n_contig_no=(*nc_ptr)->kmer_info.contig_no;
					if(n_contig_no==0)
					{
						edge_ptr=edge_ptr->nxt_edge;
						continue;
					}
					//if((*nc_ptr)->kmer_info.flip^flip_nc)
					bool t_flip=(*nc_ptr)->kmer_info.flip;

					char edge_cstr[200];
					uint64_t edgebits=edge_ptr->edge;
					int edge_len= edge_ptr->len;
					edge_len++;
					//edgebits=get_rev_comp_seq(edgebits,edge_len);
					string edge_s = bitsarr2str(&edgebits,edge_len,edge_cstr,1);

					//string edge_s = bitsarr2str(&edgebits,(int)(edge_ptr->len+1),edge_cstr,1);


					if((t_flip^flip_nc)==1)
					{

						if(contigs_info->contig_adjacency_left[contig_no].count(-n_contig_no)==0)
						{
							contigs_info->contig_adjacency_left[contig_no][-n_contig_no].cov=edge_ptr->edge_cov;
							contigs_info->contig_adjacency_left[contig_no][-n_contig_no].bridge=edge_s;
					
						}
						else
						{
						//	contigs_info->contig_adjacency_left[contig_no][-n_contig_no].cov++;
						}

					}
					else
					{
						if(contigs_info->contig_adjacency_left[contig_no].count(n_contig_no)==0)
						{
							contigs_info->contig_adjacency_left[contig_no][n_contig_no].cov=edge_ptr->edge_cov;
							contigs_info->contig_adjacency_left[contig_no][n_contig_no].bridge=edge_s;

						}
						else
						{
							//contigs_info->contig_adjacency_left[contig_no][n_contig_no].cov++;
						}

					}


					edge_ptr=edge_ptr->nxt_edge;

				}

			}
			else
			{
				//search right, app to adj left


				edge_node* edge_ptr=(*ptr)->kmer_info.right;
				while(edge_ptr!=NULL)
				{
					uint64_t t_kmer[100],f_kmer[100];
					memcpy(t_kmer,seq,Kmer_arr_sz*sizeof(uint64_t));
					
					uint64_t t=0;
					bool flip_nc=0;

					for(int j=edge_ptr->len;j>=0;--j)
					{


						uint64_t right=(uint64_t)((edge_ptr->edge)>>2*j);
						right&=0x3;
						t=3;
						t<<=((K_size-1)%32)*2;
						t_kmer[0]&=(~t);
						L_shift_NB(t_kmer,2,Kmer_arr_sz);//
						
						t_kmer[Kmer_arr_sz-1]|=right;

					}

					flip_nc=0;

					memcpy(f_kmer,t_kmer,Kmer_arr_sz*sizeof(uint64_t));
					
					get_rev_comp_seq_arr(f_kmer,K_size,Kmer_arr_sz);

					if(uint64_t_cmp(t_kmer,f_kmer,Kmer_arr_sz)>0)
					{
						
						memcpy(t_kmer,f_kmer,Kmer_arr_sz*sizeof(uint64_t));
					
						flip_nc=!flip_nc;
					}

					hv=MurmurHash64A(&t_kmer,Kmer_arr_sz*sizeof(uint64_t),0);

					hash_idx=(size_t) (hv%ht->ht_sz);


					struct bucket0 ** nc_ptr;

					nc_ptr= &(ht->store_pos[hash_idx]);

					found=look_up_in_a_list0(t_kmer,&nc_ptr,Kmer_arr_sz);


					if(found==0)
					{
						edge_ptr=edge_ptr->nxt_edge;
						//cout<<"WarningL."<<endl;
						continue;
					}
					int n_contig_no=(*nc_ptr)->kmer_info.contig_no;
					if(n_contig_no==0)
					{
						edge_ptr=edge_ptr->nxt_edge;
						continue;
					}
					//if((*nc_ptr)->kmer_info.flip^flip_nc==0)
					bool t_flip=(*nc_ptr)->kmer_info.flip;
					char edge_cstr[200];
					uint64_t edgebits=edge_ptr->edge;
					int edge_len= edge_ptr->len;
					edge_len++;
					edgebits=get_rev_comp_seq(edgebits,edge_len);
					string edge_s = bitsarr2str(&edgebits,edge_len,edge_cstr,1);
					//edgebits=get_rev_comp_seq(edgebits,(int)(edge_ptr->len+1));
					//string edge_s = bitsarr2str(&edgebits,(int)(edge_ptr->len+1),edge_cstr,1);

					if((t_flip^flip_nc)==0)
					{

						if(contigs_info->contig_adjacency_left[contig_no].count(-n_contig_no)==0)
						{
							contigs_info->contig_adjacency_left[contig_no][-n_contig_no].cov=edge_ptr->edge_cov;;
							contigs_info->contig_adjacency_left[contig_no][-n_contig_no].bridge=edge_s;

						}
						else
						{
							//contigs_info->contig_adjacency_left[contig_no][-n_contig_no].cov++;
						}

					}
					else
					{
						if(contigs_info->contig_adjacency_left[contig_no].count(n_contig_no)==0)
						{
							contigs_info->contig_adjacency_left[contig_no][n_contig_no].cov=edge_ptr->edge_cov;;
							contigs_info->contig_adjacency_left[contig_no][n_contig_no].bridge=edge_s;
					
						}
						else
						{
						//	contigs_info->contig_adjacency_left[contig_no][n_contig_no].cov++;
						}

					}




					edge_ptr=edge_ptr->nxt_edge;

				}

			}
		}

		//contig endpoint


		int gMin=contig_s.size()-K_size-100;
		if(gMin<0)
		{gMin=0;}

		for(int j=contig_s.size()-K_size;j>=gMin;--j)
		{
			string K_mer=contig_s.substr(j,K_size);


			
			str2bitsarr(K_mer.c_str(),K_size,seq,Kmer_arr_sz);
			
			memcpy(f_seq,seq,Kmer_arr_sz*sizeof(uint64_t));
			get_rev_comp_seq_arr(f_seq,K_size,Kmer_arr_sz);
			flip_c=0;
			if(uint64_t_cmp(seq,f_seq,Kmer_arr_sz)>0)
			{
				flip_c=1;
				memcpy(seq,f_seq,Kmer_arr_sz*sizeof(uint64_t));			
			}

			

			hv=MurmurHash64A(&seq,sizeof(uint64_t)*Kmer_arr_sz,0);
			hash_idx=(size_t) (hv%(ht->ht_sz));

			ptr= &(ht->store_pos[hash_idx]);
			found=look_up_in_a_list0(seq,&ptr,Kmer_arr_sz);
			if(found==0)
			{
				continue;
			}
			else
			{
				break;
			}
		}

		if(found==0)
		{

//			cout<<"searching error."<<endl;
	//		return ;
		}

		if(found)
		{
			if((*ptr)->kmer_info.flip==0)
			{
				//search right, app to adj right

				edge_node* edge_ptr=(*ptr)->kmer_info.right;
				while(edge_ptr!=NULL)
				{

					uint64_t t_kmer[100],f_kmer[100];
					memcpy(t_kmer,seq,Kmer_arr_sz*sizeof(uint64_t));
					
					uint64_t t=0;
					bool flip_nc=0;


					for(int j=edge_ptr->len;j>=0;--j)
					{



						uint64_t right=(uint64_t)((edge_ptr->edge)>>2*j);
						right&=0x3;
						t=3;
						t<<=((K_size-1)%32)*2;
						t_kmer[0]&=(~t);
						L_shift_NB(t_kmer,2,Kmer_arr_sz);//
						
						t_kmer[Kmer_arr_sz-1]|=right;


					}





					flip_nc=0;

					memcpy(f_kmer,t_kmer,Kmer_arr_sz*sizeof(uint64_t));
					
					get_rev_comp_seq_arr(f_kmer,K_size,Kmer_arr_sz);

					if(uint64_t_cmp(t_kmer,f_kmer,Kmer_arr_sz)>0)
					{
						
						memcpy(t_kmer,f_kmer,Kmer_arr_sz*sizeof(uint64_t));
					
						flip_nc=!flip_nc;
					}

					hv=MurmurHash64A(&t_kmer,Kmer_arr_sz*sizeof(uint64_t),0);


					hash_idx=(size_t) (hv%ht->ht_sz);


					struct bucket0 ** nc_ptr;

					nc_ptr= &(ht->store_pos[hash_idx]);

					found=look_up_in_a_list0(t_kmer,&nc_ptr,Kmer_arr_sz);


					if(found==0)
					{
						edge_ptr=edge_ptr->nxt_edge;
						//cout<<"WarningL."<<endl;
						continue;
					}
					int n_contig_no=(*nc_ptr)->kmer_info.contig_no;
					if(n_contig_no==0)
					{
						edge_ptr=edge_ptr->nxt_edge;
						continue;
					}
					//if((*nc_ptr)->kmer_info.flip^flip_nc==0)
					bool t_flip=(*nc_ptr)->kmer_info.flip;


					char edge_cstr[200];
					uint64_t edgebits=edge_ptr->edge;
					//edgebits=get_rev_comp_seq(edgebits,edge_ptr->len);
					int edge_len= edge_ptr->len;
					edge_len++;
					//edgebits=get_rev_comp_seq(edgebits,edge_len);
					string edge_s = bitsarr2str(&edgebits,edge_len,edge_cstr,1);

					//string edge_s = bitsarr2str(&edgebits,(int)(edge_ptr->len+1),edge_cstr,1);


					if((t_flip^flip_nc)==1)
					{

						if(contigs_info->contig_adjacency_right[contig_no].count(-n_contig_no)==0)
						{
							contigs_info->contig_adjacency_right[contig_no][-n_contig_no].cov=edge_ptr->edge_cov;
							contigs_info->contig_adjacency_right[contig_no][-n_contig_no].bridge=edge_s;

						}
						else
						{
						//	contigs_info->contig_adjacency_right[contig_no][-n_contig_no].cov++;
						}

					}
					else
					{
						if(contigs_info->contig_adjacency_right[contig_no].count(n_contig_no)==0)
						{
							contigs_info->contig_adjacency_right[contig_no][n_contig_no].cov=edge_ptr->edge_cov;
							contigs_info->contig_adjacency_right[contig_no][n_contig_no].bridge=edge_s;

						}
						else
						{
						//	contigs_info->contig_adjacency_right[contig_no][n_contig_no].cov++;
						}

					}




					edge_ptr=edge_ptr->nxt_edge;

				}

			}
			else
			{
				//search left, app to adj right

				edge_node* edge_ptr=(*ptr)->kmer_info.left;
				while(edge_ptr!=NULL)
				{




					uint64_t t_kmer[100],f_kmer[100];
					memcpy(t_kmer,seq,Kmer_arr_sz*sizeof(uint64_t));
					
					uint64_t t=0;
					bool flip_nc=0;
					for(int j=0;j<=edge_ptr->len;++j)
					{
						uint64_t left=(uint64_t)((edge_ptr->edge)>>2*j);

						left&=0x3;
						R_shift_NB(t_kmer,2,Kmer_arr_sz);
						uint64_t t=((uint64_t)1)<<(((K_size-1)%32)*2);
						t_kmer[0]|=t;

					}



					flip_nc=0;
					memcpy(f_kmer,t_kmer,Kmer_arr_sz*sizeof(uint64_t));
					get_rev_comp_seq_arr(f_kmer,K_size,Kmer_arr_sz);
					if(uint64_t_cmp(t_kmer,f_kmer,Kmer_arr_sz)>0)
					{
						memcpy(t_kmer,f_kmer,Kmer_arr_sz*sizeof(uint64_t));
						flip_nc=!flip_nc;
					}

					hv=MurmurHash64A(&t_kmer,sizeof(uint64_t)*Kmer_arr_sz,0);

					hash_idx=(size_t) (hv%ht->ht_sz);


					struct bucket0 ** nc_ptr;

					nc_ptr= &(ht->store_pos[hash_idx]);

					found=look_up_in_a_list0(t_kmer,&nc_ptr,Kmer_arr_sz);



					if(found==0)
					{
						edge_ptr=edge_ptr->nxt_edge;
						//cout<<"WarningR."<<endl;
						continue;
					}
					int n_contig_no=(*nc_ptr)->kmer_info.contig_no;
					if(n_contig_no==0)
					{
						edge_ptr=edge_ptr->nxt_edge;
						continue;
					}
				//	if((*nc_ptr)->kmer_info.flip^flip_nc)
					bool t_flip=(*nc_ptr)->kmer_info.flip;

					
					char edge_cstr[200];
					uint64_t edgebits=edge_ptr->edge;
					int edge_len= edge_ptr->len;
					edge_len++;
					edgebits=get_rev_comp_seq(edgebits,edge_len);
					string edge_s = bitsarr2str(&edgebits,edge_len,edge_cstr,1);
					//edgebits=get_rev_comp_seq(edgebits,(int)(edge_ptr->len+1));
					//string edge_s = bitsarr2str(&edgebits,(int)(edge_ptr->len+1),edge_cstr,1);


					if((t_flip^flip_nc)==0)
					{

						if(contigs_info->contig_adjacency_right[contig_no].count(-n_contig_no)==0)
						{
							contigs_info->contig_adjacency_right[contig_no][-n_contig_no].cov=edge_ptr->edge_cov;
							contigs_info->contig_adjacency_right[contig_no][-n_contig_no].bridge=edge_s;

						}
						else
						{
							//	contigs_info->contig_adjacency_right[contig_no][-n_contig_no].cov++;
						}

					}
					else
					{
						if(contigs_info->contig_adjacency_right[contig_no].count(n_contig_no)==0)
						{
							contigs_info->contig_adjacency_right[contig_no][n_contig_no].cov=edge_ptr->edge_cov;
							contigs_info->contig_adjacency_right[contig_no][n_contig_no].bridge=edge_s;

						}
						else
						{
							//contigs_info->contig_adjacency_right[contig_no][n_contig_no].cov++;
						}

					}


					edge_ptr=edge_ptr->nxt_edge;

				}
			}

		}

	}


	out_ctg_graph<<contigs_info->total_contigs<<endl;
	for (int i=1;i<=contigs_info->total_contigs;++i)
	{
		out_ctg_graph<<i<<endl;
		out_ctg_graph<<contigs_info->contig_adjacency_left[i].size()<<endl;
		if(contigs_info->contig_adjacency_left[i].size()>0)
		{
		
			map<int,struct adjacent_contig_info>::iterator map_beg=contigs_info->contig_adjacency_left[i].begin();
			map<int,struct adjacent_contig_info>::iterator map_end=contigs_info->contig_adjacency_left[i].end();
			map<int,struct adjacent_contig_info>::iterator map_it;
			for(map_it=map_beg;map_it!=map_end;++map_it)
			{
				out_ctg_graph<<map_it->first<<" "<<K_size-map_it->second.bridge.size()<<", ";
			}
			out_ctg_graph<<endl;
		}

		out_ctg_graph<<contigs_info->contig_adjacency_right[i].size()<<endl;
		if(contigs_info->contig_adjacency_right[i].size()>0)
		{
		
			map<int,struct adjacent_contig_info>::iterator map_beg=contigs_info->contig_adjacency_right[i].begin();
			map<int,struct adjacent_contig_info>::iterator map_end=contigs_info->contig_adjacency_right[i].end();
			map<int,struct adjacent_contig_info>::iterator map_it;
			for(map_it=map_beg;map_it!=map_end;++map_it)
			{
				out_ctg_graph<<map_it->first<<" "<<K_size-map_it->second.bridge.size()<<", ";
			}
			out_ctg_graph<<endl;
		}

	}

}


int BFSearchDist(struct hashtable* ht,struct hashtable* merge_ht, struct bucket* bktptr,struct bucket* obj_bktptr,int K_size, stacked_bucket &kmer_stack_beg,int max_depth,int max_dist)
{
	map<bucket*,struct BFS_path_info > Visited_Path;
	map<struct bucket* ,int > stacked_nodes;
	struct bucket *beg_bkt= bktptr;
	int max_stack=300;
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
	//Visited_Path[new_node].last_bkt_edge=NULL;

	map<int , list<stacked_bucket> >::iterator NB_it=dist_ctgs.begin();
	while(1)
	{
		NB_it=dist_ctgs.begin();
		
		if(NB_it==dist_ctgs.end())
		{break;}
		if(NB_it->first>max_dist)
		{return -10000;}
		if(NB_it->second.size()==0)
		{dist_ctgs.erase(NB_it->first);continue;}
		if(NBs>max_stack)
		{
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
				//tip end reached so continue.
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
				{r_flip=0;}

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
					//r_found=look_up_in_a_list_rm((*ptr)->kmer_t.kmer,&ptr_rm);
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
					if((*ptr)==(obj_bktptr))
					{
						int edge_len=edge_ptr->len;
						return (int)(Visited_Path[new_node].len+edge_len+1);//found distance
					}

					// not in stack
					if(stacked_nodes[*ptr]==0)
					{

						//Visited_Path[*ptr].cov=(int)(Visited_Path[new_node].cov+edge_ptr->edge_cov);
						Visited_Path[*ptr].depth=(Visited_Path[new_node].depth+1);
						
						int edge_len=edge_ptr->len;
						int cum_len=(int)(Visited_Path[new_node].len+edge_len+1);
						Visited_Path[*ptr].len=cum_len;//(Visited_Path[new_node].len+edge_ptr->len+1);

						//Visited_Path[*ptr].last_bkt=new_node;
						//Visited_Path[*ptr].last_bkt_edge=edge_ptr;
						if(r_flip)
						{
							stacked_nodes[*ptr]=-1;
							stacked_bkt.bktptr=*ptr;
							stacked_bkt.RightSearch=0;
						
						}
						else
						{
							stacked_nodes[*ptr]=1;
							stacked_bkt.bktptr=*ptr;
							stacked_bkt.RightSearch=1;
							
							
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
								//stacked_bkt.BothSideSearch=1;
								
							}
							else
							{
								stacked_nodes[*ptr]=3;
								stacked_bkt.bktptr=*ptr;
								stacked_bkt.RightSearch=1;
								//stacked_bkt.BothSideSearch=1;
							
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
					if((*ptr)==obj_bktptr)
					{
						int edge_len=edge_ptr->len;
						return (int)(Visited_Path[new_node].len+edge_len+1);//found
						
					}

					if(stacked_nodes[*ptr]==0)
					{
						//Visited_Path[*ptr].cov=(int)(Visited_Path[new_node].cov+edge_ptr->edge_cov);
						Visited_Path[*ptr].depth=(Visited_Path[new_node].depth+1);
						int edge_len=edge_ptr->len;
						int cum_len=(int)(Visited_Path[new_node].len+edge_len+1);
						Visited_Path[*ptr].len=cum_len;

						//Visited_Path[*ptr].last_bkt=new_node;
						//Visited_Path[*ptr].last_bkt_edge=edge_ptr;
						if(l_flip)
						{
							stacked_nodes[*ptr]=1;
							stacked_bkt.bktptr=*ptr;
							stacked_bkt.RightSearch=1;
							//kmer_stack.push_back(stacked_bkt);
						}
						else
						{
							stacked_nodes[*ptr]=-1;
							stacked_bkt.bktptr=*ptr;
							stacked_bkt.RightSearch=0;
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
								//kmer_stack.push_back(stacked_bkt);
							}
							else
							{
								stacked_nodes[*ptr]=3;
								stacked_bkt.bktptr=*ptr;
								stacked_bkt.RightSearch=0;
								//kmer_stack.push_back(stacked_bkt);
							}
							dist_ctgs[cum_len].push_back(stacked_bkt);
							NBs++;
							//don't do anything,since both strands are visited.

						

						}
					}

				}





			}




		}

	}

	return -10000;//return a large negative number.
}



int BFSearchDist2(struct hashtable2* ht,struct hashtable2* merge_ht, struct bucket2* bktptr,struct bucket2* obj_bktptr,int K_size, stacked_bucket2 &kmer_stack_beg,int max_depth,int max_dist)
{
	map<bucket2*,struct BFS_path_info > Visited_Path;
	map<struct bucket2* ,int > stacked_nodes;
	struct bucket2 *beg_bkt= bktptr;
	int max_stack=300;
	int DepthTh=max_depth;//min(300/gap,20);
	int LenTh=max_dist;
	bool RIGHT=0;
	struct stacked_bucket2 stacked_bkt=kmer_stack_beg;

	map<int , list<stacked_bucket2> > dist_ctgs;//neighborset
	dist_ctgs[0].push_back(kmer_stack_beg);
	int NBs=1;

	int dist_searched=0;
	bucket2* new_node=stacked_bkt.bktptr;
	kmer_t2 kmer,f_kmer;
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
	//Visited_Path[new_node].last_bkt_edge=NULL;

	map<int , list<stacked_bucket2> >::iterator NB_it=dist_ctgs.begin();
	while(1)
	{
		NB_it=dist_ctgs.begin();
		
		if(NB_it==dist_ctgs.end())
		{break;}
		if(NB_it->first>max_dist)
		{return -10000;}
		if(NB_it->second.size()==0)
		{dist_ctgs.erase(NB_it->first);continue;}
		if(NBs>max_stack)
		{
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
				//tip end reached so continue.
				continue;
				
			}

			while(edge_ptr!=NULL)
			{
							kmer=(new_node)->kmer_t2;

				f_kmer=kmer;
				get_rev_comp_seq_arr(f_kmer.kmer,K_size,2);
				int edge_len=(int) (edge_ptr->len);


				for(int g=edge_len;g>=0;--g)
				{
					uint64_t b=(edge_ptr->edge>>2*g)&0x3;
					kmer.kmer[0]&=(~((uint64_t)0x3<<(2*(K_size-1-32))));
					L_shift_NB(kmer.kmer,2,2);
					kmer.kmer[1]|=b;



				}
				bool r_flip=0;
				kmer_t2 kmer2,f_kmer2;
				kmer2=kmer;
				f_kmer2=kmer2;
				get_rev_comp_seq_arr(f_kmer2.kmer,K_size,2);
				if(uint64_t_cmp(kmer2.kmer,f_kmer2.kmer,2)>0)
				{
					kmer2=f_kmer2;
					r_flip=1;
				}
				else
				{r_flip=0;}

				uint64_t hv=MurmurHash64A(&kmer2,sizeof(kmer2),0);
				uint64_t hash_idx=(size_t) (hv%ht->ht_sz);

				struct bucket2** ptr;
				ptr= &(ht->store_pos[hash_idx]);
				r_found=look_up_in_a_list2(&kmer2,&ptr);



				if(r_found&&(*ptr)->kmer_info.removed==1)
				{

					kmer_t2 bits1=(*ptr)->kmer_t2;
					
					uint64_t hv=MurmurHash64A(&bits1,sizeof(bits1),0);
					uint64_t hash_idx=(size_t) (hv%(merge_ht->ht_sz));

					struct bucket_rm2 ** ptr_rm;
					ptr_rm=(bucket_rm2 **) &(merge_ht->store_pos[hash_idx]);
					r_found=look_up_in_a_list_rm2(&bits1,&ptr_rm);
					if(r_found==1)
					{
						r_flip^=(*ptr_rm)->flip;
						bits1=(*ptr_rm)->merged_kmer;
					}
					hv=MurmurHash64A(&bits1,sizeof(bits1),0);
					hash_idx=(size_t) (hv%(ht->ht_sz));
					ptr= &(ht->store_pos[hash_idx]);
					r_found=look_up_in_a_list2(&bits1,&ptr);

			
					
				}



				if(r_found)
				{
					if((*ptr)==(obj_bktptr))
					{
						int edge_len=edge_ptr->len;
						
						return (int)(Visited_Path[new_node].len+edge_len+1);//found distance
					}

					// not in stack
					if(stacked_nodes[*ptr]==0)
					{

						//Visited_Path[*ptr].cov=(int)(Visited_Path[new_node].cov+edge_ptr->edge_cov);
						Visited_Path[*ptr].depth=(Visited_Path[new_node].depth+1);
						int edge_len=edge_ptr->len;
						int cum_len=(int)(Visited_Path[new_node].len+edge_len+1);
						Visited_Path[*ptr].len=cum_len;//(Visited_Path[new_node].len+edge_ptr->len+1);

						//Visited_Path[*ptr].last_bkt=new_node;
						//Visited_Path[*ptr].last_bkt_edge=edge_ptr;
						if(r_flip)
						{
							stacked_nodes[*ptr]=-1;
							stacked_bkt.bktptr=*ptr;
							stacked_bkt.RightSearch=0;
						
						}
						else
						{
							stacked_nodes[*ptr]=1;
							stacked_bkt.bktptr=*ptr;
							stacked_bkt.RightSearch=1;
							
							
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
								//stacked_bkt.BothSideSearch=1;
								
							}
							else
							{
								stacked_nodes[*ptr]=3;
								stacked_bkt.bktptr=*ptr;
								stacked_bkt.RightSearch=1;
								//stacked_bkt.BothSideSearch=1;
							
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
				kmer=(new_node)->kmer_t2;
				f_kmer=kmer;
				get_rev_comp_seq_arr(f_kmer.kmer,K_size,2);
				int edge_len=(int)(edge_ptr->len);

				for(int g=0;g<=edge_len;++g)
				{
					uint64_t b=(edge_ptr->edge>>2*g)&0x3;

					//kmer.kmer[0]&=(~((uint64_t)0x3<<(K_size-1-32)));

					R_shift_NB(kmer.kmer,2,2);
					//R_shift_NB(kmer.kmer,2,2);
					kmer.kmer[0]|=(b<<(2*(K_size-1-32)));
				}

				bool l_flip=0;
				kmer_t2 kmer2,f_kmer2;;

				kmer2=kmer;
				f_kmer2=kmer2;

				get_rev_comp_seq_arr(f_kmer2.kmer,K_size,2);

				if(uint64_t_cmp(kmer2.kmer,f_kmer2.kmer,2)>0)
				{
					kmer2=f_kmer2;
					l_flip=1;
				}
				else
				{l_flip=0;}

				uint64_t hv=MurmurHash64A(&kmer2,sizeof(kmer2),0);

				uint64_t hash_idx=(size_t) (hv%ht->ht_sz);

				struct bucket2** ptr;

				ptr= &(ht->store_pos[hash_idx]);

				l_found=look_up_in_a_list2(&kmer2,&ptr);
			

				if(l_found&&(*ptr)->kmer_info.removed==1)
				{
					//cout<<"WarningL"<<endl;

					kmer_t2 bits1=(*ptr)->kmer_t2;
					
					uint64_t hv=MurmurHash64A(&bits1,sizeof(bits1),0);
					uint64_t hash_idx=(size_t) (hv%(merge_ht->ht_sz));

					struct bucket_rm2 ** ptr_rm;
					ptr_rm=(bucket_rm2 **) &(merge_ht->store_pos[hash_idx]);
					l_found=look_up_in_a_list_rm2(&bits1,&ptr_rm);
					if(l_found==1)
					{
						l_flip^=(*ptr_rm)->flip;
						bits1=(*ptr_rm)->merged_kmer;
					}
					hv=MurmurHash64A(&bits1,sizeof(bits1),0);
					hash_idx=(size_t) (hv%(ht->ht_sz));
					ptr= &(ht->store_pos[hash_idx]);
					l_found=look_up_in_a_list2(&bits1,&ptr);

					
				}


				if(l_found)
				{
					if((*ptr)==obj_bktptr)
					{
						int edge_len=edge_ptr->len;
						return (int)(Visited_Path[new_node].len+edge_len+1);//found
						
					}

					if(stacked_nodes[*ptr]==0)
					{
						//Visited_Path[*ptr].cov=(int)(Visited_Path[new_node].cov+edge_ptr->edge_cov);
						Visited_Path[*ptr].depth=(Visited_Path[new_node].depth+1);
						int edge_len=edge_ptr->len;
						int cum_len=(int)(Visited_Path[new_node].len+edge_len+1);
						Visited_Path[*ptr].len=cum_len;

						//Visited_Path[*ptr].last_bkt=new_node;
						//Visited_Path[*ptr].last_bkt_edge=edge_ptr;
						if(l_flip)
						{
							stacked_nodes[*ptr]=1;
							stacked_bkt.bktptr=*ptr;
							stacked_bkt.RightSearch=1;
							//kmer_stack.push_back(stacked_bkt);
						}
						else
						{
							stacked_nodes[*ptr]=-1;
							stacked_bkt.bktptr=*ptr;
							stacked_bkt.RightSearch=0;
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
								//kmer_stack.push_back(stacked_bkt);
							}
							else
							{
								stacked_nodes[*ptr]=3;
								stacked_bkt.bktptr=*ptr;
								stacked_bkt.RightSearch=0;
								//kmer_stack.push_back(stacked_bkt);
							}
							dist_ctgs[cum_len].push_back(stacked_bkt);
							NBs++;
							//don't do anything,since both strands are visited.

						

						}
					}

				}





			}




		}

	}

	return -10000;//return a large negative number.
}


int BFSearchDist3(struct hashtable3* ht,struct hashtable3* merge_ht, struct bucket3* bktptr,struct bucket3* obj_bktptr,int K_size, stacked_bucket3 &kmer_stack_beg,int max_depth,int max_dist)
{
	map<bucket3*,struct BFS_path_info > Visited_Path;
	map<struct bucket3* ,int > stacked_nodes;
	struct bucket3 *beg_bkt= bktptr;
	int max_stack=300;
	int DepthTh=max_depth;//min(300/gap,20);
	int LenTh=max_dist;
	bool RIGHT=0;
	struct stacked_bucket3 stacked_bkt=kmer_stack_beg;

	map<int , list<stacked_bucket3> > dist_ctgs;//neighborset
	dist_ctgs[0].push_back(kmer_stack_beg);
	int NBs=1;

	int dist_searched=0;
	bucket3* new_node=stacked_bkt.bktptr;
	kmer_t3 kmer,f_kmer;
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
	//Visited_Path[new_node].last_bkt_edge=NULL;

	map<int , list<stacked_bucket3> >::iterator NB_it=dist_ctgs.begin();
	while(1)
	{
		NB_it=dist_ctgs.begin();
		
		if(NB_it==dist_ctgs.end())
		{break;}
		if(NB_it->first>max_dist)
		{return -10000;}
		if(NB_it->second.size()==0)
		{dist_ctgs.erase(NB_it->first);continue;}
		if(NBs>max_stack)
		{
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
				//tip end reached so continue.
				continue;
				
			}

					while(edge_ptr!=NULL)
			{


				kmer=(new_node)->kmer_t3;

				f_kmer=kmer;
				get_rev_comp_seq_arr(f_kmer.kmer,K_size,3);
				int edge_len=(int)(edge_ptr->len);


				for(int g=edge_len;g>=0;--g)
				{
					uint64_t b=(edge_ptr->edge>>2*g)&0x3;
					kmer.kmer[0]&=(~((uint64_t)0x3<<(2*(K_size-1-64))));
					L_shift_NB(kmer.kmer,2,3);
					kmer.kmer[2]|=b;



				}
				bool r_flip=0;
				kmer_t3 kmer2,f_kmer2;
				kmer2=kmer;
				f_kmer2=kmer2;
				get_rev_comp_seq_arr(f_kmer2.kmer,K_size,3);
				if(uint64_t_cmp(kmer2.kmer,f_kmer2.kmer,3)>0)
				{
					kmer2=f_kmer2;
					r_flip=1;
				}
				else
				{r_flip=0;}

				uint64_t hv=MurmurHash64A(&kmer2,sizeof(kmer2),0);
				uint64_t hash_idx=(size_t) (hv%ht->ht_sz);

				struct bucket3** ptr;
				ptr= &(ht->store_pos[hash_idx]);
				r_found=look_up_in_a_list3(&kmer2,&ptr);



				if(r_found&&(*ptr)->kmer_info.removed==1)
				{

					kmer_t3 bits1=(*ptr)->kmer_t3;
					
					uint64_t hv=MurmurHash64A(&bits1,sizeof(bits1),0);
					uint64_t hash_idx=(size_t) (hv%(merge_ht->ht_sz));

					struct bucket_rm3 ** ptr_rm;
					ptr_rm=(bucket_rm3 **) &(merge_ht->store_pos[hash_idx]);
					r_found=look_up_in_a_list_rm3(&bits1,&ptr_rm);
					if(r_found==1)
					{
						r_flip^=(*ptr_rm)->flip;
						bits1=(*ptr_rm)->merged_kmer;
					}
					hv=MurmurHash64A(&bits1,sizeof(bits1),0);
					hash_idx=(size_t) (hv%(ht->ht_sz));
					ptr= &(ht->store_pos[hash_idx]);
					r_found=look_up_in_a_list3(&bits1,&ptr);

			
					
				}



				if(r_found)
				{
					if((*ptr)==(obj_bktptr))
					{
						int edge_len=edge_ptr->len;
						return (int)(Visited_Path[new_node].len+edge_len+1);//found distance
					}

					// not in stack
					if(stacked_nodes[*ptr]==0)
					{

						//Visited_Path[*ptr].cov=(int)(Visited_Path[new_node].cov+edge_ptr->edge_cov);
						Visited_Path[*ptr].depth=(Visited_Path[new_node].depth+1);
						int edge_len=edge_ptr->len;
						int cum_len=(int)(Visited_Path[new_node].len+edge_len+1);
						Visited_Path[*ptr].len=cum_len;//(Visited_Path[new_node].len+edge_ptr->len+1);

						//Visited_Path[*ptr].last_bkt=new_node;
						//Visited_Path[*ptr].last_bkt_edge=edge_ptr;
						if(r_flip)
						{
							stacked_nodes[*ptr]=-1;
							stacked_bkt.bktptr=*ptr;
							stacked_bkt.RightSearch=0;
						
						}
						else
						{
							stacked_nodes[*ptr]=1;
							stacked_bkt.bktptr=*ptr;
							stacked_bkt.RightSearch=1;
							
							
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
								//stacked_bkt.BothSideSearch=1;
								
							}
							else
							{
								stacked_nodes[*ptr]=3;
								stacked_bkt.bktptr=*ptr;
								stacked_bkt.RightSearch=1;
								//stacked_bkt.BothSideSearch=1;
							
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
				kmer=(new_node)->kmer_t3;
				f_kmer=kmer;
				get_rev_comp_seq_arr(f_kmer.kmer,K_size,3);
				int edge_len=(int)(edge_ptr->len);

				for(int g=0;g<=edge_len;++g)
				{
					uint64_t b=(edge_ptr->edge>>2*g)&0x3;

					//kmer.kmer[0]&=(~((uint64_t)0x3<<(K_size-1-32)));

					R_shift_NB(kmer.kmer,2,3);
					//R_shift_NB(kmer.kmer,2,2);
					kmer.kmer[0]|=(b<<(2*(K_size-1-64)));
				}

				bool l_flip=0;
				kmer_t3 kmer2,f_kmer2;;

				kmer2=kmer;
				f_kmer2=kmer2;

				get_rev_comp_seq_arr(f_kmer2.kmer,K_size,3);

				if(uint64_t_cmp(kmer2.kmer,f_kmer2.kmer,3)>0)
				{
					kmer2=f_kmer2;
					l_flip=1;
				}
				else
				{l_flip=0;}

				uint64_t hv=MurmurHash64A(&kmer2,sizeof(kmer2),0);

				uint64_t hash_idx=(size_t) (hv%ht->ht_sz);

				struct bucket3** ptr;

				ptr= &(ht->store_pos[hash_idx]);

				l_found=look_up_in_a_list3(&kmer2,&ptr);
			
				if(l_found&&(*ptr)->kmer_info.removed==1)
				{
					//cout<<"WarningL"<<endl;

					kmer_t3 bits1=(*ptr)->kmer_t3;
					
					uint64_t hv=MurmurHash64A(&bits1,sizeof(bits1),0);
					uint64_t hash_idx=(size_t) (hv%(merge_ht->ht_sz));

					struct bucket_rm3 ** ptr_rm;
					ptr_rm=(bucket_rm3 **) &(merge_ht->store_pos[hash_idx]);
					l_found=look_up_in_a_list_rm3(&bits1,&ptr_rm);
					if(l_found==1)
					{
						l_flip^=(*ptr_rm)->flip;
						bits1=(*ptr_rm)->merged_kmer;
					}
					hv=MurmurHash64A(&bits1,sizeof(bits1),0);
					hash_idx=(size_t) (hv%(ht->ht_sz));
					ptr= &(ht->store_pos[hash_idx]);
					l_found=look_up_in_a_list3(&bits1,&ptr);

					
				}


				if(l_found)
				{
					if((*ptr)==obj_bktptr)
					{
						int edge_len=edge_ptr->len;
						return (int)(Visited_Path[new_node].len+edge_len+1);//found
						
					}

					if(stacked_nodes[*ptr]==0)
					{
						//Visited_Path[*ptr].cov=(int)(Visited_Path[new_node].cov+edge_ptr->edge_cov);
						Visited_Path[*ptr].depth=(Visited_Path[new_node].depth+1);
						
						int edge_len=edge_ptr->len;
						int cum_len=(int)(Visited_Path[new_node].len+edge_len+1);
						Visited_Path[*ptr].len=cum_len;

						//Visited_Path[*ptr].last_bkt=new_node;
						//Visited_Path[*ptr].last_bkt_edge=edge_ptr;
						if(l_flip)
						{
							stacked_nodes[*ptr]=1;
							stacked_bkt.bktptr=*ptr;
							stacked_bkt.RightSearch=1;
							//kmer_stack.push_back(stacked_bkt);
						}
						else
						{
							stacked_nodes[*ptr]=-1;
							stacked_bkt.bktptr=*ptr;
							stacked_bkt.RightSearch=0;
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
								//kmer_stack.push_back(stacked_bkt);
							}
							else
							{
								stacked_nodes[*ptr]=3;
								stacked_bkt.bktptr=*ptr;
								stacked_bkt.RightSearch=0;
								//kmer_stack.push_back(stacked_bkt);
							}
							dist_ctgs[cum_len].push_back(stacked_bkt);
							NBs++;
							//don't do anything,since both strands are visited.

						

						}
					}

				}





			}




		}

	}

	return -10000;//return a large negative number.
}


int BFSearchDist4(struct hashtable4* ht,struct hashtable4* merge_ht, struct bucket4* bktptr,struct bucket4* obj_bktptr,int K_size, stacked_bucket4 &kmer_stack_beg,int max_depth,int max_dist)
{
	map<bucket4*,struct BFS_path_info > Visited_Path;
	map<struct bucket4* ,int > stacked_nodes;
	struct bucket4 *beg_bkt= bktptr;
	int max_stack=300;
	int DepthTh=max_depth;//min(300/gap,20);
	int LenTh=max_dist;
	bool RIGHT=0;
	struct stacked_bucket4 stacked_bkt=kmer_stack_beg;

	map<int , list<stacked_bucket4> > dist_ctgs;//neighborset
	dist_ctgs[0].push_back(kmer_stack_beg);
	int NBs=1;

	int dist_searched=0;
	bucket4* new_node=stacked_bkt.bktptr;
	kmer_t4 kmer,f_kmer;
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
	//Visited_Path[new_node].last_bkt_edge=NULL;

	map<int , list<stacked_bucket4> >::iterator NB_it=dist_ctgs.begin();
	while(1)
	{
		NB_it=dist_ctgs.begin();
		
		if(NB_it==dist_ctgs.end())
		{break;}
		if(NB_it->first>max_dist)
		{return -10000;}
		if(NB_it->second.size()==0)
		{dist_ctgs.erase(NB_it->first);continue;}
		if(NBs>max_stack)
		{
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
				//tip end reached so continue.
				continue;
				
			}

			while(edge_ptr!=NULL)
			{


				kmer=(new_node)->kmer_t4;

				f_kmer=kmer;
				get_rev_comp_seq_arr(f_kmer.kmer,K_size,4);
				int edge_len=(int)(edge_ptr->len);


				for(int g=edge_len;g>=0;--g)
				{
					uint64_t b=(edge_ptr->edge>>2*g)&0x3;
					kmer.kmer[0]&=(~((uint64_t)0x3<<(2*(K_size-1-96))));
					L_shift_NB(kmer.kmer,2,4);
					kmer.kmer[3]|=b;



				}
				bool r_flip=0;
				kmer_t4 kmer2,f_kmer2;
				kmer2=kmer;
				f_kmer2=kmer2;
				get_rev_comp_seq_arr(f_kmer2.kmer,K_size,4);
				if(uint64_t_cmp(kmer2.kmer,f_kmer2.kmer,4)>0)
				{
					kmer2=f_kmer2;
					r_flip=1;
				}
				else
				{r_flip=0;}

				uint64_t hv=MurmurHash64A(&kmer2,sizeof(kmer2),0);
				uint64_t hash_idx=(size_t) (hv%ht->ht_sz);

				struct bucket4** ptr;
				ptr= &(ht->store_pos[hash_idx]);
				r_found=look_up_in_a_list4(&kmer2,&ptr);

				if(r_found&&(*ptr)->kmer_info.removed==1)
				{

					kmer_t4 bits1=(*ptr)->kmer_t4;
					
					uint64_t hv=MurmurHash64A(&bits1,sizeof(bits1),0);
					uint64_t hash_idx=(size_t) (hv%(merge_ht->ht_sz));

					struct bucket_rm4 ** ptr_rm;
					ptr_rm=(bucket_rm4 **) &(merge_ht->store_pos[hash_idx]);
					r_found=look_up_in_a_list_rm4(&bits1,&ptr_rm);
					if(r_found==1)
					{
						r_flip^=(*ptr_rm)->flip;
						bits1=(*ptr_rm)->merged_kmer;
					}
					hv=MurmurHash64A(&bits1,sizeof(bits1),0);
					hash_idx=(size_t) (hv%(ht->ht_sz));
					ptr= &(ht->store_pos[hash_idx]);
					r_found=look_up_in_a_list4(&bits1,&ptr);

			
					
				}



				if(r_found)
				{
					if((*ptr)==(obj_bktptr))
					{
						
						int edge_len=edge_ptr->len;
						return (int)(Visited_Path[new_node].len+edge_len+1);//found distance
					}

					// not in stack
					if(stacked_nodes[*ptr]==0)
					{

						//Visited_Path[*ptr].cov=(int)(Visited_Path[new_node].cov+edge_ptr->edge_cov);
						Visited_Path[*ptr].depth=(Visited_Path[new_node].depth+1);
						
						int edge_len=edge_ptr->len;
						int cum_len=(int)(Visited_Path[new_node].len+edge_len+1);
						Visited_Path[*ptr].len=cum_len;//(Visited_Path[new_node].len+edge_ptr->len+1);

						//Visited_Path[*ptr].last_bkt=new_node;
						//Visited_Path[*ptr].last_bkt_edge=edge_ptr;
						if(r_flip)
						{
							stacked_nodes[*ptr]=-1;
							stacked_bkt.bktptr=*ptr;
							stacked_bkt.RightSearch=0;
						
						}
						else
						{
							stacked_nodes[*ptr]=1;
							stacked_bkt.bktptr=*ptr;
							stacked_bkt.RightSearch=1;
							
							
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
								//stacked_bkt.BothSideSearch=1;
								
							}
							else
							{
								stacked_nodes[*ptr]=3;
								stacked_bkt.bktptr=*ptr;
								stacked_bkt.RightSearch=1;
								//stacked_bkt.BothSideSearch=1;
							
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
				kmer=(new_node)->kmer_t4;
				f_kmer=kmer;
				get_rev_comp_seq_arr(f_kmer.kmer,K_size,4);
				int edge_len=(int)(edge_ptr->len);

				for(int g=0;g<=edge_len;++g)
				{
					uint64_t b=(edge_ptr->edge>>2*g)&0x3;

					//kmer.kmer[0]&=(~((uint64_t)0x3<<(K_size-1-32)));

					R_shift_NB(kmer.kmer,2,4);
					//R_shift_NB(kmer.kmer,2,2);
					kmer.kmer[0]|=(b<<(2*(K_size-1-96)));
				}

				bool l_flip=0;
				kmer_t4 kmer2,f_kmer2;;

				kmer2=kmer;
				f_kmer2=kmer2;

				get_rev_comp_seq_arr(f_kmer2.kmer,K_size,4);

				if(uint64_t_cmp(kmer2.kmer,f_kmer2.kmer,4)>0)
				{
					kmer2=f_kmer2;
					l_flip=1;
				}
				else
				{l_flip=0;}

				uint64_t hv=MurmurHash64A(&kmer2,sizeof(kmer2),0);

				uint64_t hash_idx=(size_t) (hv%ht->ht_sz);

				struct bucket4** ptr;

				ptr= &(ht->store_pos[hash_idx]);

				l_found=look_up_in_a_list4(&kmer2,&ptr);
			

				if(l_found&&(*ptr)->kmer_info.removed==1)
				{
					//cout<<"WarningL"<<endl;

					kmer_t4 bits1=(*ptr)->kmer_t4;
					
					uint64_t hv=MurmurHash64A(&bits1,sizeof(bits1),0);
					uint64_t hash_idx=(size_t) (hv%(merge_ht->ht_sz));

					struct bucket_rm4 ** ptr_rm;
					ptr_rm=(bucket_rm4 **) &(merge_ht->store_pos[hash_idx]);
					l_found=look_up_in_a_list_rm4(&bits1,&ptr_rm);
					if(l_found==1)
					{
						l_flip^=(*ptr_rm)->flip;
						bits1=(*ptr_rm)->merged_kmer;
					}
					hv=MurmurHash64A(&bits1,sizeof(bits1),0);
					hash_idx=(size_t) (hv%(ht->ht_sz));
					ptr= &(ht->store_pos[hash_idx]);
					l_found=look_up_in_a_list4(&bits1,&ptr);

					
				}


				if(l_found)
				{
					if((*ptr)==obj_bktptr)
					{
						
							int edge_len=edge_ptr->len;
						return (int)(Visited_Path[new_node].len+edge_len+1);//found
						
					}

					if(stacked_nodes[*ptr]==0)
					{
						//Visited_Path[*ptr].cov=(int)(Visited_Path[new_node].cov+edge_ptr->edge_cov);
						Visited_Path[*ptr].depth=(Visited_Path[new_node].depth+1);
						
						int edge_len=edge_ptr->len;
						int cum_len=(int)(Visited_Path[new_node].len+edge_len+1);
						Visited_Path[*ptr].len=cum_len;

						//Visited_Path[*ptr].last_bkt=new_node;
						//Visited_Path[*ptr].last_bkt_edge=edge_ptr;
						if(l_flip)
						{
							stacked_nodes[*ptr]=1;
							stacked_bkt.bktptr=*ptr;
							stacked_bkt.RightSearch=1;
							//kmer_stack.push_back(stacked_bkt);
						}
						else
						{
							stacked_nodes[*ptr]=-1;
							stacked_bkt.bktptr=*ptr;
							stacked_bkt.RightSearch=0;
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
								//kmer_stack.push_back(stacked_bkt);
							}
							else
							{
								stacked_nodes[*ptr]=3;
								stacked_bkt.bktptr=*ptr;
								stacked_bkt.RightSearch=0;
								//kmer_stack.push_back(stacked_bkt);
							}
							dist_ctgs[cum_len].push_back(stacked_bkt);
							NBs++;
							//don't do anything,since both strands are visited.

						

						}
					}

				}





			}




		}

	}

	return -10000;//return a large negative number.
}



void AppendMergeHT(hashtable *ht,hashtable *merge_ht)
{

	kmer_t kmer;
	struct bucket_rm ** ptr_rm;
	struct bucket ** ptr;
	
	for(size_t i=0;i<(merge_ht->ht_sz);++i)
	{
		ptr_rm=(bucket_rm **) &(merge_ht->store_pos[i]);
		while((*ptr_rm)!=NULL)
		{
			kmer=(*ptr_rm)->kmer_t;
			uint64_t hv=MurmurHash64A(&kmer,sizeof(kmer),0);
			uint64_t hash_idx=(size_t) (hv%(ht->ht_sz));
			ptr= &(ht->store_pos[hash_idx]);
			bool r_found=look_up_in_a_list(kmer.kmer,&ptr);
			if(r_found==0)
			{
				(*ptr)=(struct bucket*)malloc(sizeof(struct bucket));
				memset(*ptr,0,sizeof(struct bucket));
				((struct bucket*) *ptr)->kmer_t=kmer;
				((struct bucket*) *ptr)->kmer_info.removed=1;
			}

			

			(ptr_rm)=&((*ptr_rm)->nxt_bucket);
		}
	}

}

void AppendMergeHT2(hashtable2 *ht,hashtable2 *merge_ht)
{

	kmer_t2 kmer;
	struct bucket_rm2 ** ptr_rm;
	struct bucket2 ** ptr;
	
	for(size_t i=0;i<(merge_ht->ht_sz);++i)
	{
		ptr_rm=(bucket_rm2 **) &(merge_ht->store_pos[i]);
		while((*ptr_rm)!=NULL)
		{
			kmer=(*ptr_rm)->kmer_t2;
			uint64_t hv=MurmurHash64A(&kmer,sizeof(kmer),0);
			uint64_t hash_idx=(size_t) (hv%(ht->ht_sz));
			ptr= &(ht->store_pos[hash_idx]);
			bool r_found=look_up_in_a_list2(&kmer,&ptr);
			if(r_found==0)
			{
				(*ptr)=(struct bucket2*)malloc(sizeof(struct bucket2));
				memset(*ptr,0,sizeof(struct bucket2));
				((struct bucket2*) *ptr)->kmer_t2=kmer;
				((struct bucket2*) *ptr)->kmer_info.removed=1;
			}

			

			(ptr_rm)=&((*ptr_rm)->nxt_bucket);
		}
	}

}

void AppendMergeHT3(hashtable3 *ht,hashtable3 *merge_ht)
{

	kmer_t3 kmer;
	struct bucket_rm3 ** ptr_rm;
	struct bucket3 ** ptr;
	
	for(size_t i=0;i<(merge_ht->ht_sz);++i)
	{
		ptr_rm=(bucket_rm3 **) &(merge_ht->store_pos[i]);
		while((*ptr_rm)!=NULL)
		{
			kmer=(*ptr_rm)->kmer_t3;
			uint64_t hv=MurmurHash64A(&kmer,sizeof(kmer),0);
			uint64_t hash_idx=(size_t) (hv%(ht->ht_sz));
			ptr= &(ht->store_pos[hash_idx]);
			bool r_found=look_up_in_a_list3(&kmer,&ptr);
			if(r_found==0)
			{
				(*ptr)=(struct bucket3*)malloc(sizeof(struct bucket3));
				memset(*ptr,0,sizeof(struct bucket3));
				((struct bucket3*) *ptr)->kmer_t3=kmer;
				((struct bucket3*) *ptr)->kmer_info.removed=1;
			}

			

			(ptr_rm)=&((*ptr_rm)->nxt_bucket);
		}
	}

}

void AppendMergeHT4(hashtable4 *ht,hashtable4 *merge_ht)
{

	kmer_t4 kmer;
	struct bucket_rm4 ** ptr_rm;
	struct bucket4 ** ptr;
	
	for(size_t i=0;i<(merge_ht->ht_sz);++i)
	{
		ptr_rm=(bucket_rm4 **) &(merge_ht->store_pos[i]);
		while((*ptr_rm)!=NULL)
		{
			kmer=(*ptr_rm)->kmer_t4;
			uint64_t hv=MurmurHash64A(&kmer,sizeof(kmer),0);
			uint64_t hash_idx=(size_t) (hv%(ht->ht_sz));
			ptr= &(ht->store_pos[hash_idx]);
			bool r_found=look_up_in_a_list4(&kmer,&ptr);
			if(r_found==0)
			{
				(*ptr)=(struct bucket4*)malloc(sizeof(struct bucket4));
				memset(*ptr,0,sizeof(struct bucket4));
				((struct bucket4*) *ptr)->kmer_t4=kmer;
				((struct bucket4*) *ptr)->kmer_info.removed=1;
			}

			

			(ptr_rm)=&((*ptr_rm)->nxt_bucket);
		}
	}

}

void AppendMergeHT0(hashtable0 *ht,hashtable0 *merge_ht,int Kmer_arr_sz)
{

	uint64_t kmer[100];
	struct bucket_rm0 ** ptr_rm;
	struct bucket0 ** ptr;

	for(size_t i=0;i<(merge_ht->ht_sz);++i)
	{
		ptr_rm=(bucket_rm0 **) &(merge_ht->store_pos[i]);
		while((*ptr_rm)!=NULL)
		{
			memcpy(kmer,(*ptr_rm)->kmer_t,sizeof(uint64_t)*Kmer_arr_sz);
			
			uint64_t hv=MurmurHash64A(kmer,sizeof(uint64_t)*Kmer_arr_sz,0);
			uint64_t hash_idx=(size_t) (hv%(ht->ht_sz));
			ptr= &(ht->store_pos[hash_idx]);
			bool r_found=look_up_in_a_list0(kmer,&ptr,Kmer_arr_sz);
			if(r_found==0)
			{
				(*ptr)=(struct bucket0*)malloc(sizeof(struct bucket0));

				memset(*ptr,0,sizeof(struct bucket0));

				(*ptr)->kmer_t=kmer;
				(*ptr)->kmer_info.removed=1;
			}
			(ptr_rm)=&((*ptr_rm)->nxt_bucket);
		}
	}

}



void RemoveUnmappedNodes(hashtable *ht,hashtable2 *ht2,int K_size)
{
	
	if(K_size<=32)
	{
		bucket *bktptr;
		for(size_t i=0;i<ht->ht_sz;++i)
		{
			
			bktptr=ht->store_pos[i];
			while(bktptr!=NULL)
			{
				if(bktptr->kmer_info.contig_no==0)
				{
					bktptr->kmer_info.removed=1;
				}
				bktptr=bktptr->nxt_bucket;
			}
				
		}

	}
	else
	{
		if(K_size>32&&K_size<=64)
		{
		

			bucket2 *bktptr;
			for(size_t i=0;i<ht2->ht_sz;++i)
			{
				bktptr=ht2->store_pos[i];
				while(bktptr!=NULL)
				{
					if(bktptr->kmer_info.contig_no==0)
					{
						bktptr->kmer_info.removed=1;
					}
					bktptr=bktptr->nxt_bucket;
				}
			}

		}
	}

	int64_t bucket_count=0,edge_cnt=0;
	
	if(K_size<=32)
	{
		RemovingWeakNodesAndEdges(ht, K_size,0, 0,&bucket_count, &edge_cnt);
	}
	else
	{
		if(K_size>32&&K_size<=64)
		{
			RemovingWeakNodesAndEdges2(ht2, K_size,0, 0,&bucket_count, &edge_cnt);	
		}
	}
}


void RemoveUnmappedNodes3(hashtable3 *ht,int K_size)
{
		

	bucket3 *bktptr;
	for(size_t i=0;i<ht->ht_sz;++i)
	{
		bktptr=ht->store_pos[i];
		while(bktptr!=NULL)
		{
			if(bktptr->kmer_info.contig_no==0)
			{
				bktptr->kmer_info.removed=1;
			}
			bktptr=bktptr->nxt_bucket;
		}
	}

	int64_t bucket_count=0,edge_cnt=0;
	
	RemovingWeakNodesAndEdges3(ht, K_size,0, 0,&bucket_count, &edge_cnt);	
	
}


void RemoveUnmappedNodes4(hashtable4 *ht,int K_size)
{
		

	bucket4 *bktptr;
	for(size_t i=0;i<ht->ht_sz;++i)
	{
		bktptr=ht->store_pos[i];
		while(bktptr!=NULL)
		{
			if(bktptr->kmer_info.contig_no==0)
			{
				bktptr->kmer_info.removed=1;
			}
			bktptr=bktptr->nxt_bucket;
		}
	}

	int64_t bucket_count=0,edge_cnt=0;
	
	RemovingWeakNodesAndEdges4(ht, K_size,0, 0,&bucket_count, &edge_cnt);	
	
}


void ContigGapEst(struct hashtable *ht1,struct hashtable *merge_ht1, struct hashtable2 *ht2, struct hashtable2 *merge_ht2,int K_size,vector<int> &insert_sz_vt,vector<string>& filenames_vt,vector<bool> &LongLib,struct contigs_info * contigs_info,string ContigFilename,bool ResumePE,int64_t totReads,bool MatePair)
{
	bool Iter_Scaffold=0;
	if(ContigFilename=="SuperContigs.txt")
	{
		Iter_Scaffold=1;
	}

	int MaxDepth=300,MaxSearchLen=2000;
	
	ofstream o_log;
	int lib_cnt=1;
	string pe_name;
	if(Iter_Scaffold==0)
	{
		pe_name="InsertSizeEst.txt";
	}
	else
	{
		pe_name="InsertSizeEst.txt";//pe_name="InsertSizeEstIterated.txt";
	}
	if(0)//ResumePE)
	{
		ifstream in_pe_info(pe_name.c_str(),ios::app);
		string str,s;
		int ins_est;
		vector<int> ins_est_vt;
		while(in_pe_info>>str>>ins_est)
		{
			ins_est_vt.push_back(ins_est);
		}
		lib_cnt=ins_est_vt.size();
	}
	else
	{
		o_log.open(pe_name.c_str());
	}
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
	if(K_size<=32)
	{
		ht_sz=ht1->ht_sz;
	}
	else
	{
		ht_sz=ht2->ht_sz;
	}

	bool FAST=1;
	map<int,vector<int> > hit_position;

	contigs_info->contig_sz_vt.clear();
	contigs_info->contig_sz_vt.push_back(0);

	if(ContigFilename=="Contigs.txt")
	{
		ContigsRemapping(ht1,ht2, K_size, contigs_info,ContigFilename,0);
		
	}
	else
	{
		SuperContigsRemapping(ht1,ht2, K_size, contigs_info,ContigFilename,0);
	}
	RemoveUnmappedNodes(ht1,ht2, K_size);
	BuildContigAdjacency(ht1, ht2, contigs_info, K_size,  ContigFilename.c_str());
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
	if(0)//ResumePE)
	{
		lib_no=lib_cnt;
	}
	
	for(size_t ii=0;ii<filenames_vt.size();ii+=2)
	{
		int BaseInsertSize=-100000;
		vector<int> insert_sz_est_vt;
		if(insert_sz_vt.size()==(filenames_vt.size()/2))
		{
			BaseInsertSize=insert_sz_vt[ii/2];
			if(BaseInsertSize>10000)
			{MatePair=1;}
		}

		bool isLongLib=LongLib[ii/2];
		lib_no++;
		if(0)//lib_no>lib_cnt)
		{
			if(p_cnt>0)
			{o_log<<"InsertSizeEst_"<<lib_no<<": "<<int((double)dist_sum/(double)p_cnt+0.5)<<endl;}
			else
			{
			o_log<<"InsertSizeEst_"<<lib_no<<": "<<-10000<<endl;
			}	
		
		}
		char o_sc_r_n[300],o_sc_l_n[300],o_sc_pd[300],o_sc_inward[300],o_sc_outward[300],o_sc_log[300];
		if(1)//Iter_Scaffold==0)
		{
		sprintf(o_sc_r_n,"Pdist_R_lib_%d.txt",lib_no);
		sprintf(o_sc_l_n,"Pdist_L_lib_%d.txt",lib_no);
		sprintf(o_sc_pd,"Pdist_lib_%d.txt",lib_no);
		sprintf(o_sc_inward,"Pdist_inward_lib_%d.txt",lib_no);
		sprintf(o_sc_outward,"Pdist_outward_lib_%d.txt",lib_no);
		sprintf(o_sc_log,"Pdist_lib_%d_log.txt",lib_no);
		}
		else
		{
		sprintf(o_sc_r_n,"PdistIter_R_lib_%d.txt",lib_no);
		sprintf(o_sc_l_n,"PdistIter_L_lib_%d.txt",lib_no);
		sprintf(o_sc_pd,"PdistIter_lib_%d.txt",lib_no);
		sprintf(o_sc_inward,"PdistIter_inward_lib_%d.txt",lib_no);
		sprintf(o_sc_outward,"PdistIter_outward_lib_%d.txt",lib_no);
		sprintf(o_sc_log,"PdistIter_lib_%d_log.txt",lib_no);
		}
		ofstream o_Pdist_R(o_sc_r_n),o_Pdist_L(o_sc_l_n),o_Pdist(o_sc_pd),o_Pdist_in(o_sc_inward),o_Pdist_out(o_sc_outward),o_Pdist_log(o_sc_log);
		uint64_t inward_found=0,inward_not_found=0,outward_found=0,outward_not_found=0;
		dist_sum=0,p_cnt=0;
		cout<<"Processing library: "<<ii<<" & "<<ii+1<<endl;
		//scaffold_len=insert_sz_vt[ii/2];
		ifstream in_pair1,in_pair2;
		bool SINGLE_READ;

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
		while(getline(in_pair1,t1))
		{
			nLines1++;
			nLines1%=4;

			if(t1.size()==0)
			{
				nLines1--;
				continue;
			}
			if(t1[t1.size()-1]=='\n'||t1[t1.size()-1]=='\r')
			{
				t1.resize(t1.size()-1);
			}
			if(t1.size()==0)
			{
				nLines1--;
				continue;
			}

			if(fq_flag==0&&t1[0]=='@')
			{
				fq_flag=1;
				//fq_tmp=t1.substr(0,5);
			}
			
			bool bad_flag1=0,bad_flag2=0;
			int tag1_sz,tag2_sz;
			if((!fq_flag&&t1[0]=='>')||(fq_flag&&nLines1%4==1))
			{
				nLines1%=4;
				tag1_sz=t1.size();
				tag_s1=t1;

				

			}
			else
			{
				seq_s1=t1;
				if(tag_s1n.size()>0)
				{
					tag_s1=tag_s1n;
				}
				if(t1[0]!='N'&&t1[0]!='A'&&t1[0]!='C'&&t1[0]!='G'&&t1[0]!='T')
				{
					seq_s1.clear();
					continue;
				}
			}
			
			
			while(getline(in_pair1,t1))
			{
				nLines1++;
				nLines1%=4;
				if(t1.size()==0)
				{
					nLines1--;
					continue;
				}
				if(t1[t1.size()-1]=='\n'||t1[t1.size()-1]=='\r')
				{
					t1.resize(t1.size()-1);
				}
				if(t1.size()==0)
				{
					nLines1--;
					continue;
				}
				
			

				if(t1[0]!='N'&&t1[0]!='A'&&t1[0]!='C'&&t1[0]!='G'&&t1[0]!='T')
				{
				
					break;
				}
				seq_s1+=t1;
						
			}
			

			//while(t1.size()>0&&((fq_flag&&(t1.substr(0,5)!=fq_tmp))||((!fq_flag)&&t1[0]!='>')))
			while(t1.size()>0&&((fq_flag&&nLines1!=1)||((!fq_flag)&&t1[0]!='>')))
			{
				getline(in_pair1,t1);
				nLines1++;
				nLines1%=4;
	
			}
			//cout<<t1<<endl;
			if(t1.size()==0)
			{
				break;
			}
			tag_s1n=t1;
										
						
			readLen1=seq_s1.size();
					
			if (readLen1==0)
			{
				cout<<"Empty sequence!"<<endl;
				bad_flag1=1;
			}
			
			for(int i=0;i<readLen1;++i)
			{
				if(seq_s1[i]!='A'&&seq_s1[i]!='C'&&seq_s1[i]!='G'&&seq_s1[i]!='T'&&seq_s1[i]!='N')
				{
					bad_flag1=1;
					break;
				}
			}
			
			

			int nN=readLen1-1,isN=-1;

			for(int i=0;i<readLen1;++i)
			{
						
				if(seq_s1[i]=='-'||seq_s1[i]=='N')
				{
					if(i<=readLen1/2)
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
			if((nN-isN)<=readLen1/2)
			{
				bad_flag1=1;
			}
		

			if(isN>=0)
			{
				for(int i=isN+1;i<=nN;++i)
				{
					seq_s1[s]=seq_s1[i];
					s++;
				}
				seq_s1[s]='\0';
				seq_s1.resize(s);
			}

		
			
	
			if(1)
			{
				
				if(!SINGLE_READ)
				{
					getline(in_pair2,t2);
					nLines2++;
					nLines2%=4;
	
					if(t2.size()==0)
					{
					
						getline(in_pair2,t2);
					
					
					}
			

					if(t2[t2.size()-1]=='\n'||t2[t2.size()-1]=='\r')
					{
						t2.resize(t2.size()-1);
					}
					if(t2.size()==0)
					{
						getline(in_pair2,t2);
					}
			
			
			
					//if(((!fq_flag)&&t2[0]=='>')||((fq_flag)&&(t2.substr(0,5)==fq_tmp)))
					if(((!fq_flag)&&t2[0]=='>')||((fq_flag)&&(nLines2==1)))
					{
						tag2_sz=t2.size();
						tag_s2=t2;
					}
					else
					{
						seq_s2=t2;
						if(tag_s2n.size()>0)
						{
							tag_s2=tag_s2n;
						}
					
						if(t2[0]!='N'&&t2[0]!='A'&&t2[0]!='C'&&t2[0]!='G'&&t2[0]!='T')
						{
							seq_s2.clear();
						
							continue;
						}
					}
			
			
					while(getline(in_pair2,t2))
					{
						nLines2++;
						nLines2%=4;
						if(t2.size()==0)
						{
							nLines2--;
							continue;
						}
						if(t2[t2.size()-1]=='\n'||t2[t2.size()-1]=='\r')
						{
							t2.resize(t2.size()-1);
						}
						if(t2.size()==0)
						{
							nLines2--;
							continue;
						}
				
					

						if(t2[0]!='N'&&t2[0]!='A'&&t2[0]!='C'&&t2[0]!='G'&&t2[0]!='T')
						{
					
							break;
						}
						seq_s2+=t2;
						
					}
			

					//while(t2.size()>0&&((fq_flag&&(t2.substr(0,5)!=fq_tmp))||((!fq_flag)&&t2[0]!='>')))
					while(t2.size()>0&&((fq_flag&&nLines2!=1)||((!fq_flag)&&t2[0]!='>')))
					{
						getline(in_pair2,t2);
						nLines2++;
						nLines2%=4;
				
					}
					if(t2.size()==0)
					{
						break;
					}

					//cout<<t2<<endl;
				
					tag_s2n=t2;
				
										
						
					readLen2=seq_s2.size();
					
					if (readLen2==0)
					{
						cout<<"Empty sequence!"<<endl;
						bad_flag2=1;
					}
			
					for(int i=0;i<readLen2;++i)
					{
						if(seq_s2[i]!='A'&&seq_s2[i]!='C'&&seq_s2[i]!='G'&&seq_s2[i]!='T'&&seq_s2[i]!='N')
						{
							bad_flag2=1;
							break;
						}
					}
			
						
			
					

					int nN=readLen2-1,isN=-1;
					for(int i=0;i<readLen2;++i)
					{
						
						if(seq_s2[i]=='-'||seq_s2[i]=='N')
						{
							if(i<=readLen2/2)
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
					if((nN-isN)<=readLen2/2)
					{
						bad_flag2=1;
					}
		

					if(isN>=0)
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
				else
				{
					readLen2=readLen1;
					seq_s2=seq_s1;
					tag_s2=tag_s1;
				}
			}
			if(readLen1<K_size)
			{bad_flag1=1;}
			if(readLen2<K_size)
			{bad_flag2=1;}


			int64_t pos1,pos2;
			uint64_t bits1,bits2;
			kmer_t2 bits1_t2,bits2_t2;

			num_Reads++;

			if(totReads!=0&&num_Reads>totReads)
			{
				break;
			}
			if(num_Reads%10000000==0)
			{
				
				time(&read_time);
				cout<<num_Reads<<" Pairs Searched."<<endl;
				cout<<"Time: "<<difftime(read_time,beg_time)<<" secs."<<endl;
				continue;

			}
			//cout<<num_Reads<<endl;
		//if(num_Reads<2721)
		//{continue;}
				//cout<<"";

			/*
			if(bad_flag1==1||bad_flag2==1)
			{
				cout<<tag_s1<<endl;cout<<seq_s1<<endl;
				cout<<tag_s2<<endl;cout<<seq_s2<<endl<<endl;
				continue;
			}
			*/
		
			if(!SINGLE_READ)
			{
				if(isLongLib)
				{
					reverse(seq_s1.begin(),seq_s1.end());
					complement_str(seq_s1);
				}
				else
				{
					reverse(seq_s2.begin(),seq_s2.end());
					complement_str(seq_s2);
						
				}
			}
			Init_Read(seq_s1,Read1);
			seq_s1.clear();
			
		
			Init_Read(seq_s2,Read2);
			seq_s2.clear();
			
			tag1_sz=tag_s1.size();
			tag2_sz=tag_s2.size();

	
			int mismatch_cnt=0;
			bool tag_mismatch=0;
			/*
			if(abs(tag1_sz-tag2_sz)>2)
			{
				tag_mismatch=1;
				break;
			}

			int tag_sz_t=tag1_sz;
			if(tag2_sz<tag1_sz)
			{
				tag_sz_t=tag2_sz;
			}
			for(int t=0;t<tag_sz_t;++t)
			{
				if(tag_s1[t]!=tag_s2[t])
				{
					mismatch_cnt++;

					
				}
			}
			if(mismatch_cnt>3)
			{
				tag_mismatch=1;
				cout<<num_Reads<<endl;
				break;
			}

			*/

			bool found1=0,found2=0,LocalSearch=1;

			//for non adjacent relation

		
			struct bucket **ptr1_d,**ptr2_d;
			struct bucket2 **ptr1_t2d,**ptr2_t2d;
			bool flip_1d;
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
						
					//	cout<<"f"<<endl;

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
						else
						{
							bool flip_1=flip_0;
							
							while(1)//up search in the tree flip as needed,and if finally found, update the values.
							{
								uint64_t hv=MurmurHash64A(&bits1,sizeof(bits1),0);
								uint64_t hash_idx=(size_t) (hv%(merge_ht1->ht_sz));

								struct bucket_rm ** ptr;
								ptr=(bucket_rm **) &(merge_ht1->store_pos[hash_idx]);
								//r_found=look_up_in_a_list_rm((*ptr1)->kmer_t.kmer,&ptr);
								bool r_found=look_up_in_a_list_rm(bits1,&ptr);
								
								if(r_found==1)
								{
									flip_1^=(*ptr)->flip;
									bits1=(*ptr)->merged_kmer.kmer;

								}
								else
								{break;}
							}


							hv=MurmurHash64A(&bits1,sizeof(bits1),0);
							hash_idx=(size_t) (hv%ht_sz);
							struct bucket **ptr1;
							ptr1= &(ht1->store_pos[hash_idx]);
							found1=look_up_in_a_list(bits1,&ptr1);
							

							if(found1&&(!(*ptr1)->kmer_info.removed))
							{
								ptr1_d=ptr1;
								flip_1d=flip_0;
								cod1=(*ptr1)->kmer_info.cod;
								cont1=(*ptr1)->kmer_info.contig_no;
								flip_1=(*ptr1)->kmer_info.flip^flip_1;
								pos1=i;

							
								
							}
							else
							{
								found1=0;

							}
							break;


						}

						
					}
				}
				else
				{
					if(K_size>32&&K_size<=64)
					{
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
							
							found1=1;






							if((*ptr1)->kmer_info.removed==0)
							{
								ptr1_t2d=ptr1;
								flip_1d=flip_0;
								cod1=(*ptr1)->kmer_info.cod;
								cont1=(*ptr1)->kmer_info.contig_no;
								flip_1=(*ptr1)->kmer_info.flip^flip_0;
								pos1=i;

								break;		
							}
							else
							{
								bool flip_1=flip_0;
							
								while(1)//up search in the tree flip as needed,and if finally found, update the values.
								{
									uint64_t hv=MurmurHash64A(&bits1_t2,sizeof(bits1_t2),0);
									uint64_t hash_idx=(size_t) (hv%(merge_ht2->ht_sz));

									struct bucket_rm2 ** ptr;
									ptr=(bucket_rm2 **) &(merge_ht2->store_pos[hash_idx]);
									
									//r_found=look_up_in_a_list_rm2(&((*ptr1)->kmer_t2),&ptr);
									bool r_found=look_up_in_a_list_rm2(&bits1_t2,&ptr);
									
									if(r_found==1)
									{
										flip_1^=(*ptr)->flip;
										bits1_t2=(*ptr)->merged_kmer;
									}
									else
									{break;}
								}


								hv=MurmurHash64A(&bits1_t2,sizeof(bits1_t2),0);
								hash_idx=(size_t) (hv%ht_sz);
								struct bucket2 **ptr1;
								ptr1= &(ht2->store_pos[hash_idx]);
								found1=look_up_in_a_list2(&bits1_t2,&ptr1);

								if(found1&&(!(*ptr1)->kmer_info.removed))
								{
									ptr1_t2d=ptr1;
									cod1=(*ptr1)->kmer_info.cod;
									cont1=(*ptr1)->kmer_info.contig_no;
									flip_1=(*ptr1)->kmer_info.flip^flip_1;
									pos1=i;

							
								
								}
								else
								{
									found1=0;
								}
								break;


							}








						
						}
					}
				}


			}
			if(Read1.readLen<K_size||found1==0)
			{
				//cout<<"not found"<<endl;
				continue;
			}

			for(int i=Read2.readLen-K_size;i>=0;--i )
			{
				if(K_size<=32)
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
						
					//	cout<<"f"<<endl;

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
						else
						{
							bool flip_2=flip_0;
							
							while(1)//up search in the tree flip as needed,and if finally found, update the values.
							{
								uint64_t hv=MurmurHash64A(&bits2,sizeof(bits2),0);
								uint64_t hash_idx=(size_t) (hv%(merge_ht1->ht_sz));

								struct bucket_rm ** ptr;
								ptr=(bucket_rm **) &(merge_ht1->store_pos[hash_idx]);
								//r_found=look_up_in_a_list_rm((*ptr2)->kmer_t.kmer,&ptr);
								bool r_found=look_up_in_a_list_rm(bits2,&ptr);
								if(r_found==1)
								{
									flip_2^=(*ptr)->flip;
									bits2=(*ptr)->merged_kmer.kmer;
								}
								else
								{break;}
							}


							hv=MurmurHash64A(&bits2,sizeof(bits2),0);
							hash_idx=(size_t) (hv%ht_sz);
							struct bucket **ptr2;
							ptr2= &(ht1->store_pos[hash_idx]);
							found2=look_up_in_a_list(bits2,&ptr2);
							if(found2&&(!(*ptr2)->kmer_info.removed))
							{
								ptr2_d=ptr2;

								cod2=(*ptr2)->kmer_info.cod;
								cont2=(*ptr2)->kmer_info.contig_no;
								flip_2=(*ptr2)->kmer_info.flip^flip_2;
								pos2=i;
			
							}
							else
							{
								found2=0;

							}
							break;


						}

						/*

						found2=1;
						cod2=(*ptr2)->kmer_info.cod;
						cont2=(*ptr2)->kmer_info.contig_no;
						flip_2=(*ptr2)->kmer_info.flip^flip_0;
						pos2=i;

						break;

						*/
					}
				}
				else
				{
					if(K_size<=64&&K_size>32)
					{
						get_sub_arr(Read2.read_bits,Read2.readLen,i,K_size,bits2_t2.kmer);
						f_seq_t2=bits2_t2;
						get_rev_comp_seq_arr(f_seq_t2.kmer,K_size,2);
						flip_0=0;
						if(uint64_t_cmp(bits2_t2.kmer,f_seq_t2.kmer,2)>0)
						{
							bits2_t2=f_seq_t2;
							flip_0=1;
						}

						hv=MurmurHash64A(bits2_t2.kmer,sizeof(bits2_t2),0);

						hash_idx=(size_t) (hv%ht_sz);
						struct bucket2 **ptr2;
						ptr2= &(ht2->store_pos[hash_idx]);
						found=look_up_in_a_list2(&bits2_t2,&ptr2);
						if(found)
						{



							found2=1;




							if((*ptr2)->kmer_info.removed==0)
							{
								ptr2_t2d=ptr2;
								cod2=(*ptr2)->kmer_info.cod;
								cont2=(*ptr2)->kmer_info.contig_no;
								flip_2=(*ptr2)->kmer_info.flip^flip_0;
								pos2=i;

								break;		
							}
							else
							{
								bool flip_2=flip_0;
							
								while(1)//up search in the tree flip as needed,and if finally found, update the values.
								{
									uint64_t hv=MurmurHash64A(&bits2_t2,sizeof(bits2_t2),0);
									uint64_t hash_idx=(size_t) (hv%(merge_ht2->ht_sz));

									struct bucket_rm2 ** ptr;
									ptr=(bucket_rm2 **) &(merge_ht2->store_pos[hash_idx]);
									//r_found=look_up_in_a_list_rm2(&((*ptr2)->kmer_t2),&ptr);
									bool r_found=look_up_in_a_list_rm2(&bits2_t2,&ptr);
									if(r_found==1)
									{
										flip_2^=(*ptr)->flip;
										bits2_t2=(*ptr)->merged_kmer;
									}
									else
									{break;}
								}


								hv=MurmurHash64A(&bits2_t2,sizeof(bits2_t2),0);
								hash_idx=(size_t) (hv%ht_sz);
								struct bucket2 **ptr2;
								ptr2= &(ht2->store_pos[hash_idx]);
								found2=look_up_in_a_list2(&bits2_t2,&ptr2);
								if(found2&&(!(*ptr2)->kmer_info.removed))
								{
									cod2=(*ptr2)->kmer_info.cod;
									cont2=(*ptr2)->kmer_info.contig_no;
									flip_2=(*ptr2)->kmer_info.flip^flip_2;
									pos2=i;

							
								
								}
								else
								{
									found2=0;
								}
								break;


							}



							/*


							found2=1;
							cod2=(*ptr2)->kmer_info.cod;
							cont2=(*ptr2)->kmer_info.contig_no;
							if(cont2==0)
							{break;cout<<"Error Pair2."<<endl;}
							flip_2=(*ptr2)->kmer_info.flip^flip_0;
							pos2=i;

							break;*/
						}

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
				int pdist=-10000;

				
				if(MatePair&&LocalSearch)
				{

					if(K_size<=32)
					{
						struct stacked_bucket stacked_bucket;

						stacked_bucket.bktptr=*ptr1_d;
						if(flip_1d==0)
						{
							stacked_bucket.RightSearch=1;
						}
						else
						{
							stacked_bucket.RightSearch=0;
						}
						pdist= BFSearchDist(ht1, merge_ht1, *ptr1_d,*ptr2_d,K_size, stacked_bucket,MaxDepth,MaxSearchLen);
						
						if(pdist>0)
						{
							ignore=1;
							o_Pdist_in<<(*ptr1_d)->kmer_info.contig_no<<" "<<(*ptr2_d)->kmer_info.contig_no<<" "<<pdist<<endl;
							inward_found++;
						}
						else
						{
							inward_not_found++;
						}

						stacked_bucket.RightSearch=!stacked_bucket.RightSearch;
						pdist= BFSearchDist(ht1, merge_ht1, *ptr1_d,*ptr2_d,K_size, stacked_bucket,MaxDepth,MaxSearchLen);
						
						if(pdist>0)
						{
							ignore=1;
							o_Pdist_out<<(*ptr1_d)->kmer_info.contig_no<<" "<<(*ptr2_d)->kmer_info.contig_no<<" "<<pdist<<endl;
							outward_found++;
						}
						else
						{
							outward_not_found++;
						}
						
						
					}
					else
					{
						struct stacked_bucket2 stacked_bucket;

						stacked_bucket.bktptr=*ptr1_t2d;
						if(flip_1d==0)
						{
							stacked_bucket.RightSearch=1;
						}
						else
						{
							stacked_bucket.RightSearch=0;
						}
						pdist= BFSearchDist2(ht2, merge_ht2, *ptr1_t2d,*ptr2_t2d,K_size, stacked_bucket,MaxDepth,MaxSearchLen);
						
						if(pdist>0)
						{
							ignore=1;
							o_Pdist_in<<(*ptr1_t2d)->kmer_info.contig_no<<" "<<(*ptr2_t2d)->kmer_info.contig_no<<" "<<pdist<<endl;
							inward_found++;
						}
						else
						{
							inward_not_found++;
						}
						stacked_bucket.RightSearch=!stacked_bucket.RightSearch;
						pdist= BFSearchDist2(ht2, merge_ht2, *ptr1_t2d,*ptr2_t2d,K_size, stacked_bucket,MaxDepth,MaxSearchLen);
						
						if(pdist>0)
						{
							ignore=1;
							o_Pdist_out<<(*ptr1_t2d)->kmer_info.contig_no<<" "<<(*ptr2_t2d)->kmer_info.contig_no<<" "<<pdist<<endl;
							outward_found++;
						}
						else
						{
							outward_not_found++;
						}
					
					}
				}
				
				if(ignore==1)
				{
							
					continue;
				}
			

				if(flip_1==0&&(cod2>cod1))
				{
					int pdist=(int)((cod2-cod1)+pos1+readLen2-pos2);
					if(BaseInsertSize<-1000||(BaseInsertSize<=500&&((pdist<=(2*BaseInsertSize))&&(pdist>=(-100))))||((pdist<=(3*BaseInsertSize/2))&&(pdist>=(BaseInsertSize/2))))
					{
						dist_sum+=pdist;
						o_Pdist<<pdist<<endl;
						p_cnt++;
						if(insert_sz_est_vt.size()<=102)
						{
							insert_sz_est_vt.push_back(pdist);
							if(insert_sz_est_vt.size()==101)
							{
								sort(insert_sz_est_vt.begin(),insert_sz_est_vt.end());
								BaseInsertSize=insert_sz_est_vt[50];
							}
						}
					}
				}
				else
				{
					if(flip_1==1&&cod1>cod2)
					{
						int pdist=(int)((cod1-cod2)+pos1+readLen2-pos2);
						if(BaseInsertSize<-1000||(BaseInsertSize<=500&&((pdist<=(2*BaseInsertSize))&&(pdist>=(-100))))||((pdist<=(3*BaseInsertSize/2))&&(pdist>=(BaseInsertSize/2))))
						{
							dist_sum+=pdist;
							p_cnt++;
							o_Pdist<<pdist<<endl;
							if(insert_sz_est_vt.size()<=102)
							{
								insert_sz_est_vt.push_back(pdist);
								if(insert_sz_est_vt.size()==101)
								{
									sort(insert_sz_est_vt.begin(),insert_sz_est_vt.end());
									BaseInsertSize=insert_sz_est_vt[50];
								}
							}
						}

					}

				}
			}

			if(cont1!=cont2)
			{
				
				int pdist=-10000;

				
				if(MatePair&&LocalSearch)
				{

					if(K_size<=32)
					{
						struct stacked_bucket stacked_bucket;

						stacked_bucket.bktptr=*ptr1_d;
						if(flip_1d==0)
						{
							stacked_bucket.RightSearch=1;
						}
						else
						{
							stacked_bucket.RightSearch=0;
						}
						pdist= BFSearchDist(ht1, merge_ht1, *ptr1_d,*ptr2_d,K_size, stacked_bucket,MaxDepth,MaxSearchLen);
						
						if(pdist>0)
						{
							ignore=1;
							o_Pdist_in<<(*ptr1_d)->kmer_info.contig_no<<" "<<(*ptr2_d)->kmer_info.contig_no<<" "<<pdist<<endl;
							inward_found++;
						}
						else
						{
							inward_not_found++;
						}

						stacked_bucket.RightSearch=!stacked_bucket.RightSearch;
						pdist= BFSearchDist(ht1, merge_ht1, *ptr1_d,*ptr2_d,K_size, stacked_bucket,MaxDepth,MaxSearchLen);
						
						if(pdist>0)
						{
							ignore=1;
							o_Pdist_out<<(*ptr1_d)->kmer_info.contig_no<<" "<<(*ptr2_d)->kmer_info.contig_no<<" "<<pdist<<endl;
							outward_found++;
						}
						else
						{
							outward_not_found++;
						}
						
						
					}
					else
					{
						struct stacked_bucket2 stacked_bucket;

						stacked_bucket.bktptr=*ptr1_t2d;
						if(flip_1d==0)
						{
							stacked_bucket.RightSearch=1;
						}
						else
						{
							stacked_bucket.RightSearch=0;
						}
						pdist= BFSearchDist2(ht2, merge_ht2, *ptr1_t2d,*ptr2_t2d,K_size, stacked_bucket,MaxDepth,MaxSearchLen);
						
						if(pdist>0)
						{
							ignore=1;
							o_Pdist_in<<(*ptr1_t2d)->kmer_info.contig_no<<" "<<(*ptr2_t2d)->kmer_info.contig_no<<" "<<pdist<<endl;
							inward_found++;
						}
						else
						{
							inward_not_found++;
						}
						stacked_bucket.RightSearch=!stacked_bucket.RightSearch;
						pdist= BFSearchDist2(ht2, merge_ht2, *ptr1_t2d,*ptr2_t2d,K_size, stacked_bucket,MaxDepth,MaxSearchLen);
						
						if(pdist>0)
						{
							ignore=1;
							o_Pdist_out<<(*ptr1_t2d)->kmer_info.contig_no<<" "<<(*ptr2_t2d)->kmer_info.contig_no<<" "<<pdist<<endl;
							outward_found++;
						}
						else
						{
							outward_not_found++;
						}
					
					}
				}
				
				if(ignore==1)
				{
							
					continue;
				}
			


				//cout<<pdist<<endl;
				

				///////
              

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
					//	if(0)//if(Tdist<0)
					//	{
					//		Tdist=0;///////////////
					//	}

						o_Pdist_L<<cont1 <<" "<<cont2<<" "<<d4<<endl;//" Library: "<<lib_no<<endl;
						//o_Pdist_R<<cont2 <<" "<<cont1<<" "<<d4<<" Library: "<<lib_no<<endl;

					}

//					contigs_info->scaffold_adjacency_left[cont1].push_back(adj_contig);

				}




			}



		}

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
	
	
///end for all insert size


}



void ContigGapEst3(struct hashtable3 *ht,struct hashtable3 *merge_ht,int K_size,vector<int> &insert_sz_vt,vector<string>& filenames_vt,vector<bool> &LongLib,struct contigs_info * contigs_info,string ContigFilename,bool ResumePE,int64_t totReads,bool MatePair)
{
	bool Iter_Scaffold=0;
	if(ContigFilename=="SuperContigs.txt")
	{
		Iter_Scaffold=1;
	}
	int MaxDepth=300,MaxSearchLen=2000;
	string pe_name;
	if(MatePair==0)
	{
		pe_name="InsertSizeEst.txt";
	}
	else
	{
		pe_name="InsertSizeEst.txt";//pe_name="InsertSizeEstIterated.txt";
	}
	ofstream o_log;
	int lib_cnt=1;
	if(0)//ResumePE)
	{
		
		ifstream in_pe_info(pe_name.c_str(),ios::app);

		string str,s;
		int ins_est;
		vector<int> ins_est_vt;
		while(in_pe_info>>str>>ins_est)
		{
			ins_est_vt.push_back(ins_est);
		}
		lib_cnt=ins_est_vt.size();
		
	}
	else
	{
		o_log.open(pe_name.c_str());
	}
	time_t beg_time,read_time;
	int64_t scaffold_len=0;
	string in_fname=ContigFilename;
	ifstream Contigs_in(in_fname.c_str());
	size_t num_Contigs=0;
	string tag,s,kmer;
	uint64_t hv;
	struct kmer_t3 f_seq_t3;
	size_t hash_idx;
	//vector<int> contig_sz;
	bool found;
	bool flip_1,flip_2,flip_0;
	size_t ht_sz;
	
	ht_sz=ht->ht_sz;
	

	bool FAST=1;
	map<int,vector<int> > hit_position;

	contigs_info->contig_sz_vt.clear();
	contigs_info->contig_sz_vt.push_back(0);
	if(ContigFilename=="Contigs.txt")
	{
		ContigsRemapping3(ht, K_size, contigs_info,ContigFilename,0);
	
	}
	else
	{
		SuperContigsRemapping3(ht, K_size, contigs_info,ContigFilename,0);
	}
	RemoveUnmappedNodes3(ht, K_size);
	BuildContigAdjacency3(ht, contigs_info, K_size,  ContigFilename.c_str());
	///////////////marking finished
	AppendMergeHT3(ht,merge_ht);

	cout<<"Collecting paired ends information."<<endl;
	time(&beg_time);

	int64_t dist_sum=0,p_cnt=0;
	int lib_no=0;
	if(0)//ResumePE)
	{
		lib_no=lib_cnt;
	}
	
	
	for(size_t ii=0;ii<filenames_vt.size();ii+=2)
	{
		int BaseInsertSize=-100000;
		vector<int> insert_sz_est_vt;
		if(insert_sz_vt.size()==(filenames_vt.size()/2))
		{
			BaseInsertSize=insert_sz_vt[ii/2];
			if(BaseInsertSize>10000)
			{MatePair=1;}
		}
		bool isLongLib=LongLib[ii/2];
		lib_no++;
		if(0)//lib_no>lib_cnt)
		{
			if(p_cnt>0)
			{o_log<<"InsertSizeEst_"<<lib_no<<": "<<int((double)dist_sum/(double)p_cnt+0.5)<<endl;}
			else
			{
			o_log<<"InsertSizeEst_"<<lib_no<<": "<<-10000<<endl;
			}
		}

		char o_sc_r_n[300],o_sc_l_n[300],o_sc_pd[300],o_sc_inward[300],o_sc_outward[300],o_sc_log[300];
		
		if(1)//Iter_Scaffold==0)
		{
		sprintf(o_sc_r_n,"Pdist_R_lib_%d.txt",lib_no);
		sprintf(o_sc_l_n,"Pdist_L_lib_%d.txt",lib_no);
		sprintf(o_sc_pd,"Pdist_lib_%d.txt",lib_no);
		sprintf(o_sc_inward,"Pdist_inward_lib_%d.txt",lib_no);
		sprintf(o_sc_outward,"Pdist_outward_lib_%d.txt",lib_no);
		sprintf(o_sc_log,"Pdist_lib_%d_log.txt",lib_no);
		}
		else
		{
		sprintf(o_sc_r_n,"PdistIter_R_lib_%d.txt",lib_no);
		sprintf(o_sc_l_n,"PdistIter_L_lib_%d.txt",lib_no);
		sprintf(o_sc_pd,"PdistIter_lib_%d.txt",lib_no);
		sprintf(o_sc_inward,"PdistIter_inward_lib_%d.txt",lib_no);
		sprintf(o_sc_outward,"PdistIter_outward_lib_%d.txt",lib_no);
		sprintf(o_sc_log,"PdistIter_lib_%d_log.txt",lib_no);
		}
		ofstream o_Pdist_R(o_sc_r_n),o_Pdist_L(o_sc_l_n),o_Pdist(o_sc_pd),o_Pdist_in(o_sc_inward),o_Pdist_out(o_sc_outward),o_Pdist_log(o_sc_log);
		uint64_t inward_found=0,inward_not_found=0,outward_found=0,outward_not_found=0;
		
		dist_sum=0,p_cnt=0;
		cout<<"Processing library: "<<ii<<" & "<<ii+1<<endl;
		//scaffold_len=insert_sz_vt[ii/2];
		ifstream in_pair1,in_pair2;
		bool SINGLE_READ;

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
		while(getline(in_pair1,t1))
		{
			nLines1++;
			nLines1%=4;

			if(t1.size()==0)
			{
				nLines1--;
				continue;
			}
			if(t1[t1.size()-1]=='\n'||t1[t1.size()-1]=='\r')
			{
				t1.resize(t1.size()-1);
			}
			if(t1.size()==0)
			{
				nLines1--;
				continue;
			}

			if(fq_flag==0&&t1[0]=='@')
			{
				fq_flag=1;
				//fq_tmp=t1.substr(0,5);
			}
			
			bool bad_flag1=0,bad_flag2=0;
			int tag1_sz,tag2_sz;
			if((!fq_flag&&t1[0]=='>')||(fq_flag&&nLines1%4==1))
			{
				nLines1%=4;
				tag1_sz=t1.size();
				tag_s1=t1;

				

			}
			else
			{
				seq_s1=t1;
				if(tag_s1n.size()>0)
				{
					tag_s1=tag_s1n;
				}
				if(t1[0]!='N'&&t1[0]!='A'&&t1[0]!='C'&&t1[0]!='G'&&t1[0]!='T')
				{
					seq_s1.clear();
					continue;
				}
			}
			
			
			while(getline(in_pair1,t1))
			{
				nLines1++;
				nLines1%=4;
				if(t1.size()==0)
				{
					nLines1--;
					continue;
				}
				if(t1[t1.size()-1]=='\n'||t1[t1.size()-1]=='\r')
				{
					t1.resize(t1.size()-1);
				}
				if(t1.size()==0)
				{
					nLines1--;
					continue;
				}
				
			

				if(t1[0]!='N'&&t1[0]!='A'&&t1[0]!='C'&&t1[0]!='G'&&t1[0]!='T')
				{
				
					break;
				}
				seq_s1+=t1;
						
			}
			

			//while(t1.size()>0&&((fq_flag&&(t1.substr(0,5)!=fq_tmp))||((!fq_flag)&&t1[0]!='>')))
			while(t1.size()>0&&((fq_flag&&nLines1!=1)||((!fq_flag)&&t1[0]!='>')))
			{
				getline(in_pair1,t1);
				nLines1++;
				nLines1%=4;
	
			}
			//cout<<t1<<endl;
			if(t1.size()==0)
			{
				break;
			}
			tag_s1n=t1;
										
						
			readLen1=seq_s1.size();
					
			if (readLen1==0)
			{
				cout<<"Empty sequence!"<<endl;
				bad_flag1=1;
			}
			
			for(int i=0;i<readLen1;++i)
			{
				if(seq_s1[i]!='A'&&seq_s1[i]!='C'&&seq_s1[i]!='G'&&seq_s1[i]!='T'&&seq_s1[i]!='N')
				{
					bad_flag1=1;
					break;
				}
			}
			
			

			int nN=readLen1-1,isN=-1;

			for(int i=0;i<readLen1;++i)
			{
						
				if(seq_s1[i]=='-'||seq_s1[i]=='N')
				{
					if(i<=readLen1/2)
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
			if((nN-isN)<=readLen1/2)
			{
				bad_flag1=1;
			}
		

			if(isN>=0)
			{
				for(int i=isN+1;i<=nN;++i)
				{
					seq_s1[s]=seq_s1[i];
					s++;
				}
				seq_s1[s]='\0';
				seq_s1.resize(s);
			}

		
			
	
			if(1)
			{
				
				if(!SINGLE_READ)
				{
					getline(in_pair2,t2);
					nLines2++;
					nLines2%=4;
	
					if(t2.size()==0)
					{
					
						getline(in_pair2,t2);
					
					
					}
			

					if(t2[t2.size()-1]=='\n'||t2[t2.size()-1]=='\r')
					{
						t2.resize(t2.size()-1);
					}
					if(t2.size()==0)
					{
						getline(in_pair2,t2);
					}
			
			
			
					//if(((!fq_flag)&&t2[0]=='>')||((fq_flag)&&(t2.substr(0,5)==fq_tmp)))
					if(((!fq_flag)&&t2[0]=='>')||((fq_flag)&&(nLines2==1)))
					{
						tag2_sz=t2.size();
						tag_s2=t2;
					}
					else
					{
						seq_s2=t2;
						if(tag_s2n.size()>0)
						{
							tag_s2=tag_s2n;
						}
					
						if(t2[0]!='N'&&t2[0]!='A'&&t2[0]!='C'&&t2[0]!='G'&&t2[0]!='T')
						{
							seq_s2.clear();
						
							continue;
						}
					}
			
			
					while(getline(in_pair2,t2))
					{
						nLines2++;
						nLines2%=4;
						if(t2.size()==0)
						{
							nLines2--;
							continue;
						}
						if(t2[t2.size()-1]=='\n'||t2[t2.size()-1]=='\r')
						{
							t2.resize(t2.size()-1);
						}
						if(t2.size()==0)
						{
							nLines2--;
							continue;
						}
				
					

						if(t2[0]!='N'&&t2[0]!='A'&&t2[0]!='C'&&t2[0]!='G'&&t2[0]!='T')
						{
					
							break;
						}
						seq_s2+=t2;
						
					}
			

					//while(t2.size()>0&&((fq_flag&&(t2.substr(0,5)!=fq_tmp))||((!fq_flag)&&t2[0]!='>')))
					while(t2.size()>0&&((fq_flag&&nLines2!=1)||((!fq_flag)&&t2[0]!='>')))
					{
						getline(in_pair2,t2);
						nLines2++;
						nLines2%=4;
				
					}
					if(t2.size()==0)
					{
						break;
					}

					//cout<<t2<<endl;
				
					tag_s2n=t2;
				
										
						
					readLen2=seq_s2.size();
					
					if (readLen2==0)
					{
						cout<<"Empty sequence!"<<endl;
						bad_flag2=1;
					}
			
					for(int i=0;i<readLen2;++i)
					{
						if(seq_s2[i]!='A'&&seq_s2[i]!='C'&&seq_s2[i]!='G'&&seq_s2[i]!='T'&&seq_s2[i]!='N')
						{

							bad_flag2=1;
							break;
						}
					}
			
						
			
					
				


					int nN=readLen2-1,isN=-1;
					for(int i=0;i<readLen2;++i)
					{
						
						if(seq_s2[i]=='-'||seq_s2[i]=='N')
						{
							if(i<=readLen2/2)
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
					if((nN-isN)<=readLen2/2)
					{
						bad_flag2=1;
					}
		

					if(isN>=0)
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
				else
				{
					readLen2=readLen1;
					seq_s2=seq_s1;
					tag_s2=tag_s1;
				}
			}
			if(readLen1<K_size)
			{bad_flag1=1;}
			if(readLen2<K_size)
			{bad_flag2=1;}


			int64_t pos1,pos2;
			kmer_t3 bits1_t3,bits2_t3;

			num_Reads++;
			if(totReads!=0&&num_Reads>totReads)
			{
				break;
			}
			if(num_Reads%10000000==0)
			{
				
				time(&read_time);
				cout<<num_Reads<<" Pairs Searched."<<endl;
				cout<<"Time: "<<difftime(read_time,beg_time)<<" secs."<<endl;
				continue;

			}

			/*
			if(bad_flag1==1||bad_flag2==1)
			{
				cout<<tag_s1<<endl;cout<<seq_s1<<endl;
				cout<<tag_s2<<endl;cout<<seq_s2<<endl<<endl;
				continue;
			}
			*/
		
			if(!SINGLE_READ)
			{
				if(isLongLib)
				{
					reverse(seq_s1.begin(),seq_s1.end());
					complement_str(seq_s1);
				}
				else
				{
					reverse(seq_s2.begin(),seq_s2.end());
					complement_str(seq_s2);
						
				}
			}
			Init_Read(seq_s1,Read1);
			seq_s1.clear();
			
		
			Init_Read(seq_s2,Read2);
			seq_s2.clear();
			
			tag1_sz=tag_s1.size();
			tag2_sz=tag_s2.size();

	
			int mismatch_cnt=0;
			bool tag_mismatch=0;
			/*
			if(abs(tag1_sz-tag2_sz)>2)
			{
				tag_mismatch=1;
				break;
			}

			int tag_sz_t=tag1_sz;
			if(tag2_sz<tag1_sz)
			{
				tag_sz_t=tag2_sz;
			}
			for(int t=0;t<tag_sz_t;++t)
			{
				if(tag_s1[t]!=tag_s2[t])
				{
					mismatch_cnt++;

					
				}
			}
			if(mismatch_cnt>3)
			{
				tag_mismatch=1;
				cout<<num_Reads<<endl;
				break;
			}

			*/

			bool found1=0,found2=0,LocalSearch=1;

			
		
			struct bucket3 **ptr1_d,**ptr2_d;
			//struct bucket2 **ptr1_t2d,**ptr2_t2d;
			bool flip_1d;

			//for non adjacent relation

			//	cout<<num_Reads<<endl;
		
			for(int i=0;i<Read1.readLen-K_size+1;++i )
			{
				
				
				get_sub_arr(Read1.read_bits,Read1.readLen,i,K_size,bits1_t3.kmer);

				f_seq_t3=bits1_t3;
				get_rev_comp_seq_arr(f_seq_t3.kmer,K_size,3);
				flip_0=0;
				if(uint64_t_cmp(bits1_t3.kmer,f_seq_t3.kmer,3)>0)
				{
					bits1_t3=f_seq_t3;
					flip_0=1;
				}

				hv=MurmurHash64A(bits1_t3.kmer,sizeof(bits1_t3),0);

				hash_idx=(size_t) (hv%ht_sz);
				struct bucket3 **ptr1;
				ptr1= &(ht->store_pos[hash_idx]);
				found=look_up_in_a_list3(&bits1_t3,&ptr1);
				if(found)
				{
					found1=1;






					if((*ptr1)->kmer_info.removed==0)
					{
						flip_1d=flip_0;
						ptr1_d=ptr1;
						cod1=(*ptr1)->kmer_info.cod;
						cont1=(*ptr1)->kmer_info.contig_no;
						flip_1=(*ptr1)->kmer_info.flip^flip_0;
						pos1=i;

						break;		
					}
					else
					{
						bool flip_1=flip_0;
							
						while(1)//up search in the tree flip as needed,and if finally found, update the values.
						{
							uint64_t hv=MurmurHash64A(&bits1_t3,sizeof(bits1_t3),0);
							uint64_t hash_idx=(size_t) (hv%(merge_ht->ht_sz));

							struct bucket_rm3 ** ptr;
							ptr=(bucket_rm3 **) &(merge_ht->store_pos[hash_idx]);
							//r_found=look_up_in_a_list_rm3(&((*ptr1)->kmer_t3),&ptr);
							bool r_found=look_up_in_a_list_rm3(&bits1_t3,&ptr);
							if(r_found==1)
							{
								flip_1^=(*ptr)->flip;
								bits1_t3=(*ptr)->merged_kmer;
							}
							else
							{break;}
						}


						hv=MurmurHash64A(&bits1_t3,sizeof(bits1_t3),0);
						hash_idx=(size_t) (hv%ht_sz);
						struct bucket3 **ptr1;
						ptr1= &(ht->store_pos[hash_idx]);
						found1=look_up_in_a_list3(&bits1_t3,&ptr1);
						if(found1&&(!(*ptr1)->kmer_info.removed))
						{
							ptr1_d=ptr1;

							cod1=(*ptr1)->kmer_info.cod;
							cont1=(*ptr1)->kmer_info.contig_no;
							flip_1=(*ptr1)->kmer_info.flip^flip_1;
							pos1=i;

							
								
						}
						else
						{
							found1=0;
						}
						break;


					}








						
				}
			
			


			}
			if(Read1.readLen<K_size||found1==0)
			{
				//cout<<"not found"<<endl;
				continue;
			}

			for(int i=Read2.readLen-K_size;i>=0;--i )
			{
				get_sub_arr(Read2.read_bits,Read2.readLen,i,K_size,bits2_t3.kmer);
				f_seq_t3=bits2_t3;
				get_rev_comp_seq_arr(f_seq_t3.kmer,K_size,3);
				flip_0=0;
				if(uint64_t_cmp(bits2_t3.kmer,f_seq_t3.kmer,3)>0)
				{
					bits2_t3=f_seq_t3;
					flip_0=1;
				}

				hv=MurmurHash64A(bits2_t3.kmer,sizeof(bits2_t3),0);

				hash_idx=(size_t) (hv%ht_sz);
				struct bucket3 **ptr2;
				ptr2= &(ht->store_pos[hash_idx]);
				found=look_up_in_a_list3(&bits2_t3,&ptr2);
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
					else
					{
						bool flip_2=flip_0;
							
						while(1)//up search in the tree flip as needed,and if finally found, update the values.
						{
							uint64_t hv=MurmurHash64A(&bits2_t3,sizeof(bits2_t3),0);
							uint64_t hash_idx=(size_t) (hv%(merge_ht->ht_sz));

							struct bucket_rm3 ** ptr;
							ptr=(bucket_rm3 **) &(merge_ht->store_pos[hash_idx]);
							//r_found=look_up_in_a_list_rm3(&((*ptr2)->kmer_t3),&ptr);
							bool r_found=look_up_in_a_list_rm3(&bits2_t3,&ptr);
							if(r_found==1)
							{
								flip_2^=(*ptr)->flip;
								bits2_t3=(*ptr)->merged_kmer;
							}
							else
							{break;}
						}


						hv=MurmurHash64A(&bits2_t3,sizeof(bits2_t3),0);
						hash_idx=(size_t) (hv%ht_sz);
						struct bucket3 **ptr2;
						ptr2= &(ht->store_pos[hash_idx]);
						found2=look_up_in_a_list3(&bits2_t3,&ptr2);
						if(found2&&(!(*ptr2)->kmer_info.removed))
						{
							ptr2_d=ptr2;
							cod2=(*ptr2)->kmer_info.cod;
							cont2=(*ptr2)->kmer_info.contig_no;
							flip_2=(*ptr2)->kmer_info.flip^flip_2;
							pos2=i;
	
						}
						else
						{
							found2=0;
						}
						break;


					}



					/*


					found2=1;
					cod2=(*ptr2)->kmer_info.cod;
					cont2=(*ptr2)->kmer_info.contig_no;
					if(cont2==0)
					{break;cout<<"Error Pair2."<<endl;}
					flip_2=(*ptr2)->kmer_info.flip^flip_0;
					pos2=i;

					break;*/
				}


			}
			if(Read2.readLen<K_size||found2==0)
			{continue;}

			if(cont1==0||cont2==0)
			{
				continue;
			}

			if(cont1==cont2&&flip_1==flip_2)
			{
				if(flip_1==0&&(cod2>cod1))
				{
					
					int pdist=(int)((cod2-cod1)+pos1+readLen2-pos2);
					if(BaseInsertSize<-1000||(BaseInsertSize<=500&&((pdist<=(2*BaseInsertSize))&&(pdist>=(-100))))||((pdist<=(3*BaseInsertSize/2))&&(pdist>=(BaseInsertSize/2))))
					{
						dist_sum+=pdist;
						o_Pdist<<pdist<<endl;
						p_cnt++;
						if(insert_sz_est_vt.size()<=102)
						{
							insert_sz_est_vt.push_back(pdist);
							if(insert_sz_est_vt.size()==101)
							{
								sort(insert_sz_est_vt.begin(),insert_sz_est_vt.end());
								BaseInsertSize=insert_sz_est_vt[50];
							}
						}
					}
				}
				else
				{
					if(flip_1==1&&cod1>cod2)
					{
						int pdist=(int)((cod1-cod2)+pos1+readLen2-pos2);
						if(BaseInsertSize<-1000||(BaseInsertSize<=500&&((pdist<=(2*BaseInsertSize))&&(pdist>=(-100))))||((pdist<=(3*BaseInsertSize/2))&&(pdist>=(BaseInsertSize/2))))
						{

							dist_sum+=pdist;
							p_cnt++;
							o_Pdist<<pdist<<endl;
							if(insert_sz_est_vt.size()<=102)
							{
								insert_sz_est_vt.push_back(pdist);
								if(insert_sz_est_vt.size()==101)
								{
									sort(insert_sz_est_vt.begin(),insert_sz_est_vt.end());
									BaseInsertSize=insert_sz_est_vt[50];
								}
							}
						}
					}

				}
			}
			bool ignore=0;
			if(cont1!=cont2)
			{
				int pdist=-10000;

				
				if(MatePair&&LocalSearch)
				{
					struct stacked_bucket3 stacked_bucket;

					stacked_bucket.bktptr=*ptr1_d;
					if(flip_1d==0)
					{
						stacked_bucket.RightSearch=1;
					}
					else
					{
						stacked_bucket.RightSearch=0;
					}
					pdist= BFSearchDist3(ht, merge_ht, *ptr1_d,*ptr2_d,K_size, stacked_bucket,MaxDepth,MaxSearchLen);
					
					if(pdist>0)
					{
						ignore=1;
						o_Pdist_in<<(*ptr1_d)->kmer_info.contig_no<<" "<<(*ptr2_d)->kmer_info.contig_no<<" "<<pdist<<endl;
						inward_found++;
					}
					else
					{
						inward_not_found++;
					}
					stacked_bucket.RightSearch=!stacked_bucket.RightSearch;
					pdist= BFSearchDist3(ht, merge_ht, *ptr1_d,*ptr2_d,K_size, stacked_bucket,MaxDepth,MaxSearchLen);
					if(pdist>0)
					{
						ignore=1;
						o_Pdist_out<<(*ptr1_d)->kmer_info.contig_no<<" "<<(*ptr2_d)->kmer_info.contig_no<<" "<<pdist<<endl;
						outward_found++;
					}
					else
					{
						outward_not_found++;
					}

				}
				if(ignore==1)
				{continue;}
			
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
					//	if(0)//if(Tdist<0)
					//	{
					//		Tdist=0;///////////////
					//	}

						o_Pdist_L<<cont1 <<" "<<cont2<<" "<<d4<<endl;//" Library: "<<lib_no<<endl;
						//o_Pdist_R<<cont2 <<" "<<cont1<<" "<<d4<<" Library: "<<lib_no<<endl;

					}

//					contigs_info->scaffold_adjacency_left[cont1].push_back(adj_contig);

				}




			}



		}
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

void ContigGapEst4(struct hashtable4 *ht,struct hashtable4 *merge_ht,int K_size,vector<int> &insert_sz_vt,vector<string>& filenames_vt,vector<bool> &LongLib,struct contigs_info * contigs_info,string ContigFilename,bool ResumePE,int64_t totReads,bool MatePair)
{
	bool Iter_Scaffold=0;
	if(ContigFilename=="SuperContigs.txt")
	{
		Iter_Scaffold=1;
	}
	int MaxDepth=300,MaxSearchLen=2000;
	string pe_name;
	if(MatePair==0)
	{
		pe_name="InsertSizeEst.txt";
	}
	else
	{
		pe_name="InsertSizeEst.txt";//pe_name="InsertSizeEstIterated.txt";
	}
	ofstream o_log;
	int lib_cnt=1;
	if(0)//ResumePE)
	{
		ifstream in_pe_info(pe_name.c_str(),ios::app);
		string str,s;
		int ins_est;
		vector<int> ins_est_vt;
		while(in_pe_info>>str>>ins_est)
		{
			ins_est_vt.push_back(ins_est);
		}
		lib_cnt=ins_est_vt.size();
		
	}
	else
	{
		o_log.open(pe_name.c_str());
	}
	time_t beg_time,read_time;
	int64_t scaffold_len=0;
	string in_fname=ContigFilename;
	ifstream Contigs_in(in_fname.c_str());
	size_t num_Contigs=0;
	string tag,s,kmer;
	uint64_t hv;
	struct kmer_t4 f_seq_t4;
	size_t hash_idx;
	//vector<int> contig_sz;
	bool found;
	bool flip_1,flip_2,flip_0;
	size_t ht_sz;
	
	ht_sz=ht->ht_sz;
	

	bool FAST=1;
	map<int,vector<int> > hit_position;

	contigs_info->contig_sz_vt.clear();
	contigs_info->contig_sz_vt.push_back(0);

	if(ContigFilename=="Contigs.txt")
	{
		ContigsRemapping4(ht, K_size, contigs_info,ContigFilename,0);
	
	}
	else
	{
		SuperContigsRemapping4(ht, K_size, contigs_info,ContigFilename,0);
	}
	RemoveUnmappedNodes4(ht, K_size);
	BuildContigAdjacency4(ht, contigs_info, K_size,  ContigFilename.c_str());
	///////////////marking finished
	AppendMergeHT4(ht,merge_ht);

	cout<<"Collecting paired ends information."<<endl;
	time(&beg_time);

	int64_t dist_sum=0,p_cnt=0;
	int lib_no=0;
	if(0)//ResumePE)
	{
		lib_no=lib_cnt;
	}
	

	
	for(size_t ii=0;ii<filenames_vt.size();ii+=2)
	{
		int BaseInsertSize=-100000;
		vector<int> insert_sz_est_vt;
		if(insert_sz_vt.size()==(filenames_vt.size()/2))
		{
			BaseInsertSize=insert_sz_vt[ii/2];	
			if(BaseInsertSize>10000)
			{MatePair=1;}
		}
		bool isLongLib=LongLib[ii/2];
		lib_no++;
		if(0)//lib_no>lib_cnt)
		{
			if(p_cnt>0)
			{o_log<<"InsertSizeEst_"<<lib_no<<": "<<int((double)dist_sum/(double)p_cnt+0.5)<<endl;}
			else
			{
			o_log<<"InsertSizeEst_"<<lib_no<<": "<<-10000<<endl;
			}	
		
		}
		
		
		char o_sc_r_n[300],o_sc_l_n[300],o_sc_pd[300],o_sc_inward[300],o_sc_outward[300],o_sc_log[300];
		if(1)//Iter_Scaffold==0)
		{
		sprintf(o_sc_r_n,"Pdist_R_lib_%d.txt",lib_no);
		sprintf(o_sc_l_n,"Pdist_L_lib_%d.txt",lib_no);
		sprintf(o_sc_pd,"Pdist_lib_%d.txt",lib_no);
		sprintf(o_sc_inward,"Pdist_inward_lib_%d.txt",lib_no);
		sprintf(o_sc_outward,"Pdist_outward_lib_%d.txt",lib_no);
		sprintf(o_sc_log,"Pdist_lib_%d_log.txt",lib_no);
		}
		else
		{
		sprintf(o_sc_r_n,"PdistIter_R_lib_%d.txt",lib_no);
		sprintf(o_sc_l_n,"PdistIter_L_lib_%d.txt",lib_no);
		sprintf(o_sc_pd,"PdistIter_lib_%d.txt",lib_no);
		sprintf(o_sc_inward,"PdistIter_inward_lib_%d.txt",lib_no);
		sprintf(o_sc_outward,"PdistIter_outward_lib_%d.txt",lib_no);
		sprintf(o_sc_log,"PdistIter_lib_%d_log.txt",lib_no);
		}
		
		ofstream o_Pdist_R(o_sc_r_n),o_Pdist_L(o_sc_l_n),o_Pdist(o_sc_pd),o_Pdist_in(o_sc_inward),o_Pdist_out(o_sc_outward),o_Pdist_log(o_sc_log);
		uint64_t inward_found=0,inward_not_found=0,outward_found=0,outward_not_found=0;
		
		dist_sum=0,p_cnt=0;
		cout<<"Processing library: "<<ii<<" & "<<ii+1<<endl;
		//scaffold_len=insert_sz_vt[ii/2];
		ifstream in_pair1,in_pair2;
		bool SINGLE_READ;

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
		while(getline(in_pair1,t1))
		{
			nLines1++;
			nLines1%=4;

			if(t1.size()==0)
			{
				nLines1--;
				continue;
			}
			if(t1[t1.size()-1]=='\n'||t1[t1.size()-1]=='\r')
			{
				t1.resize(t1.size()-1);
			}
			if(t1.size()==0)
			{
				nLines1--;
				continue;
			}

			if(fq_flag==0&&t1[0]=='@')
			{
				fq_flag=1;
				//fq_tmp=t1.substr(0,5);
			}
			
			bool bad_flag1=0,bad_flag2=0;
			int tag1_sz,tag2_sz;
			if((!fq_flag&&t1[0]=='>')||(fq_flag&&nLines1%4==1))
			{
				nLines1%=4;
				tag1_sz=t1.size();
				tag_s1=t1;

				

			}
			else
			{
				seq_s1=t1;
				if(tag_s1n.size()>0)
				{
					tag_s1=tag_s1n;
				}
				if(t1[0]!='N'&&t1[0]!='A'&&t1[0]!='C'&&t1[0]!='G'&&t1[0]!='T')
				{
					seq_s1.clear();
					continue;
				}
			}
			
			
			while(getline(in_pair1,t1))
			{
				nLines1++;
				nLines1%=4;
				if(t1.size()==0)
				{
					nLines1--;
					continue;
				}
				if(t1[t1.size()-1]=='\n'||t1[t1.size()-1]=='\r')
				{
					t1.resize(t1.size()-1);
				}
				if(t1.size()==0)
				{
					nLines1--;
					continue;
				}
				
			

				if(t1[0]!='N'&&t1[0]!='A'&&t1[0]!='C'&&t1[0]!='G'&&t1[0]!='T')
				{
				
					break;
				}
				seq_s1+=t1;
						
			}
			

			//while(t1.size()>0&&((fq_flag&&(t1.substr(0,5)!=fq_tmp))||((!fq_flag)&&t1[0]!='>')))
			while(t1.size()>0&&((fq_flag&&nLines1!=1)||((!fq_flag)&&t1[0]!='>')))
			{
				getline(in_pair1,t1);
				nLines1++;
				nLines1%=4;
	
			}
			//cout<<t1<<endl;
			if(t1.size()==0)
			{
				break;
			}
			tag_s1n=t1;
										
						
			readLen1=seq_s1.size();
					
			if (readLen1==0)
			{
				cout<<"Empty sequence!"<<endl;
				bad_flag1=1;
			}
			
			for(int i=0;i<readLen1;++i)
			{
				if(seq_s1[i]!='A'&&seq_s1[i]!='C'&&seq_s1[i]!='G'&&seq_s1[i]!='T'&&seq_s1[i]!='N')
				{
					bad_flag1=1;
					break;
				}
			}
			
			int nN=readLen1-1,isN=-1;

			for(int i=0;i<readLen1;++i)
			{
						
				if(seq_s1[i]=='-'||seq_s1[i]=='N')
				{
					if(i<=readLen1/2)
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
			if((nN-isN)<=readLen1/2)
			{
				bad_flag1=1;
			}
		

			if(isN>=0)
			{
				for(int i=isN+1;i<=nN;++i)
				{
					seq_s1[s]=seq_s1[i];
					s++;
				}
				seq_s1[s]='\0';
				seq_s1.resize(s);
			}

		
			
	
			if(1)
			{
				
				if(!SINGLE_READ)
				{
					getline(in_pair2,t2);
					nLines2++;
					nLines2%=4;
	
					if(t2.size()==0)
					{
					
						getline(in_pair2,t2);
					
					
					}
			

					if(t2[t2.size()-1]=='\n'||t2[t2.size()-1]=='\r')
					{
						t2.resize(t2.size()-1);
					}
					if(t2.size()==0)
					{
						getline(in_pair2,t2);
					}
			
			
			
					//if(((!fq_flag)&&t2[0]=='>')||((fq_flag)&&(t2.substr(0,5)==fq_tmp)))
					if(((!fq_flag)&&t2[0]=='>')||((fq_flag)&&(nLines2==1)))
					{
						tag2_sz=t2.size();
						tag_s2=t2;
					}
					else
					{
						seq_s2=t2;
						if(tag_s2n.size()>0)
						{
							tag_s2=tag_s2n;
						}
					
						if(t2[0]!='N'&&t2[0]!='A'&&t2[0]!='C'&&t2[0]!='G'&&t2[0]!='T')
						{
							seq_s2.clear();
						
							continue;
						}
					}
			
			
					while(getline(in_pair2,t2))
					{
						nLines2++;
						nLines2%=4;
						if(t2.size()==0)
						{
							nLines2--;
							continue;
						}
						if(t2[t2.size()-1]=='\n'||t2[t2.size()-1]=='\r')
						{
							t2.resize(t2.size()-1);
						}
						if(t2.size()==0)
						{
							nLines2--;
							continue;
						}
				
					

						if(t2[0]!='N'&&t2[0]!='A'&&t2[0]!='C'&&t2[0]!='G'&&t2[0]!='T')
						{
					
							break;
						}
						seq_s2+=t2;
						
					}
			

					//while(t2.size()>0&&((fq_flag&&(t2.substr(0,5)!=fq_tmp))||((!fq_flag)&&t2[0]!='>')))
					while(t2.size()>0&&((fq_flag&&nLines2!=1)||((!fq_flag)&&t2[0]!='>')))
					{
						getline(in_pair2,t2);
						nLines2++;
						nLines2%=4;
				
					}
					if(t2.size()==0)
					{
						break;
					}

					//cout<<t2<<endl;
				
					tag_s2n=t2;
				
										
						
					readLen2=seq_s2.size();
					
					if (readLen2==0)
					{
						cout<<"Empty sequence!"<<endl;
						bad_flag2=1;
					}
			
					for(int i=0;i<readLen2;++i)
					{
						if(seq_s2[i]!='A'&&seq_s2[i]!='C'&&seq_s2[i]!='G'&&seq_s2[i]!='T'&&seq_s2[i]!='N')
						{

							bad_flag2=1;
							break;
						}
					}
			
						
			
					
				


					int nN=readLen2-1,isN=-1;
					for(int i=0;i<readLen2;++i)
					{
						
						if(seq_s2[i]=='-'||seq_s2[i]=='N')
						{
							if(i<=readLen2/2)
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
					if((nN-isN)<=readLen2/2)
					{
						bad_flag2=1;
					}
		

					if(isN>=0)
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
				else
				{
					readLen2=readLen1;
					seq_s2=seq_s1;
					tag_s2=tag_s1;
				}
			}
			if(readLen1<K_size)
			{bad_flag1=1;}
			if(readLen2<K_size)
			{bad_flag2=1;}


			int64_t pos1,pos2;
			kmer_t4 bits1_t4,bits2_t4;

			num_Reads++;
			if(totReads!=0&&num_Reads>totReads)
			{
				break;
			}
			if(num_Reads%10000000==0)
			{
				
				time(&read_time);
				cout<<num_Reads<<" Pairs Searched."<<endl;
				cout<<"Time: "<<difftime(read_time,beg_time)<<" secs."<<endl;
				continue;

			}

			/*
			if(bad_flag1==1||bad_flag2==1)
			{
				cout<<tag_s1<<endl;cout<<seq_s1<<endl;
				cout<<tag_s2<<endl;cout<<seq_s2<<endl<<endl;
				continue;
			}
			*/
		
			if(!SINGLE_READ)
			{
				if(isLongLib)
				{
					reverse(seq_s1.begin(),seq_s1.end());
					complement_str(seq_s1);
				}
				else
				{
					reverse(seq_s2.begin(),seq_s2.end());
					complement_str(seq_s2);		
				}
			}
			Init_Read(seq_s1,Read1);
			seq_s1.clear();
			
		
			Init_Read(seq_s2,Read2);
			seq_s2.clear();
			
			tag1_sz=tag_s1.size();
			tag2_sz=tag_s2.size();

	
			int mismatch_cnt=0;
			bool tag_mismatch=0;
			/*
			if(abs(tag1_sz-tag2_sz)>2)
			{
				tag_mismatch=1;
				break;
			}

			int tag_sz_t=tag1_sz;
			if(tag2_sz<tag1_sz)
			{
				tag_sz_t=tag2_sz;
			}
			for(int t=0;t<tag_sz_t;++t)
			{
				if(tag_s1[t]!=tag_s2[t])
				{
					mismatch_cnt++;

					
				}
			}
			if(mismatch_cnt>3)
			{
				tag_mismatch=1;
				cout<<num_Reads<<endl;
				break;
			}

			*/

			bool found1=0,found2=0,LocalSearch=1;

			//for non adjacent relation

			//	cout<<num_Reads<<endl;
			
		
			struct bucket4 **ptr1_d,**ptr2_d;
			bool flip_1d;
		
			for(int i=0;i<Read1.readLen-K_size+1;++i )
			{
				
				
				get_sub_arr(Read1.read_bits,Read1.readLen,i,K_size,bits1_t4.kmer);

				f_seq_t4=bits1_t4;
				get_rev_comp_seq_arr(f_seq_t4.kmer,K_size,4);
				flip_0=0;
				if(uint64_t_cmp(bits1_t4.kmer,f_seq_t4.kmer,4)>0)
				{
					bits1_t4=f_seq_t4;
					flip_0=1;
				}

				hv=MurmurHash64A(bits1_t4.kmer,sizeof(bits1_t4),0);

				hash_idx=(size_t) (hv%ht_sz);
				struct bucket4 **ptr1;
				ptr1= &(ht->store_pos[hash_idx]);
				found=look_up_in_a_list4(&bits1_t4,&ptr1);
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
					else
					{
						bool flip_1=flip_0;
							
						while(1)//up search in the tree flip as needed,and if finally found, update the values.
						{
							uint64_t hv=MurmurHash64A(&bits1_t4,sizeof(bits1_t4),0);
							uint64_t hash_idx=(size_t) (hv%(merge_ht->ht_sz));

							struct bucket_rm4 ** ptr;
							ptr=(bucket_rm4 **) &(merge_ht->store_pos[hash_idx]);
							//r_found=look_up_in_a_list_rm4(&((*ptr1)->kmer_t4),&ptr);
							bool r_found=look_up_in_a_list_rm4(&bits1_t4,&ptr);
							if(r_found==1)
							{
								flip_1^=(*ptr)->flip;
								bits1_t4=(*ptr)->merged_kmer;
							}
							else
							{break;}
						}


						hv=MurmurHash64A(&bits1_t4,sizeof(bits1_t4),0);
						hash_idx=(size_t) (hv%ht_sz);
						struct bucket4 **ptr1;
						ptr1= &(ht->store_pos[hash_idx]);
						found1=look_up_in_a_list4(&bits1_t4,&ptr1);
						if(found1&&(!(*ptr1)->kmer_info.removed))
						{
							ptr1_d=ptr1;
							cod1=(*ptr1)->kmer_info.cod;
							cont1=(*ptr1)->kmer_info.contig_no;
							flip_1=(*ptr1)->kmer_info.flip^flip_1;
							pos1=i;

							
								
						}
						else
						{
							found1=0;
						}
						break;


					}
	
				}
			
			


			}
			if(Read1.readLen<K_size||found1==0)
			{
				//cout<<"not found"<<endl;
				continue;
			}

			for(int i=Read2.readLen-K_size;i>=0;--i )
			{
				get_sub_arr(Read2.read_bits,Read2.readLen,i,K_size,bits2_t4.kmer);
				f_seq_t4=bits2_t4;
				get_rev_comp_seq_arr(f_seq_t4.kmer,K_size,4);
				flip_0=0;
				if(uint64_t_cmp(bits2_t4.kmer,f_seq_t4.kmer,4)>0)
				{
					bits2_t4=f_seq_t4;
					flip_0=1;
				}

				hv=MurmurHash64A(bits2_t4.kmer,sizeof(bits2_t4),0);

				hash_idx=(size_t) (hv%ht_sz);
				struct bucket4 **ptr2;
				ptr2= &(ht->store_pos[hash_idx]);
				found=look_up_in_a_list4(&bits2_t4,&ptr2);
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
					else
					{
						bool flip_2=flip_0;
							
						while(1)//up search in the tree flip as needed,and if finally found, update the values.
						{
							uint64_t hv=MurmurHash64A(&bits2_t4,sizeof(bits2_t4),0);
							uint64_t hash_idx=(size_t) (hv%(merge_ht->ht_sz));

							struct bucket_rm4 ** ptr;
							ptr=(bucket_rm4 **) &(merge_ht->store_pos[hash_idx]);
							//r_found=look_up_in_a_list_rm4(&((*ptr2)->kmer_t4),&ptr);
							bool r_found=look_up_in_a_list_rm4(&bits2_t4,&ptr);
							if(r_found==1)
							{
								flip_2^=(*ptr)->flip;
								bits2_t4=(*ptr)->merged_kmer;
							}
							else
							{break;}
						}


						hv=MurmurHash64A(&bits2_t4,sizeof(bits2_t4),0);
						hash_idx=(size_t) (hv%ht_sz);
						struct bucket4 **ptr2;
						ptr2= &(ht->store_pos[hash_idx]);
						found2=look_up_in_a_list4(&bits2_t4,&ptr2);
						if(found2&&(!(*ptr2)->kmer_info.removed))
						{
							ptr2_d=ptr2;
							cod2=(*ptr2)->kmer_info.cod;
							cont2=(*ptr2)->kmer_info.contig_no;
							flip_2=(*ptr2)->kmer_info.flip^flip_2;
							pos2=i;
	
						}
						else
						{
							found2=0;
						}
						break;


					}



					/*


					found2=1;
					cod2=(*ptr2)->kmer_info.cod;
					cont2=(*ptr2)->kmer_info.contig_no;
					if(cont2==0)
					{break;cout<<"Error Pair2."<<endl;}
					flip_2=(*ptr2)->kmer_info.flip^flip_0;
					pos2=i;

					break;*/
				}


			}
			if(Read2.readLen<K_size||found2==0)
			{continue;}

			if(cont1==0||cont2==0)
			{
				continue;
			}

			if(cont1==cont2&&flip_1==flip_2)
			{
				if(flip_1==0&&(cod2>cod1))
				{
					int pdist=(int)((cod2-cod1)+pos1+readLen2-pos2);
					
					if(BaseInsertSize<-1000||(BaseInsertSize<=500&&((pdist<=(2*BaseInsertSize))&&(pdist>=(-100))))||((pdist<=(3*BaseInsertSize/2))&&(pdist>=(BaseInsertSize/2))))
					{
						
						dist_sum+=pdist;
						o_Pdist<<pdist<<endl;
						p_cnt++;
						if(insert_sz_est_vt.size()<=102)
						{
							insert_sz_est_vt.push_back(pdist);
							if(insert_sz_est_vt.size()==101)
							{
								sort(insert_sz_est_vt.begin(),insert_sz_est_vt.end());
								BaseInsertSize=insert_sz_est_vt[50];
							}
						}
					}
				}
				else
				{
					
					if(flip_1==1&&cod1>cod2)
					{
						int pdist=(int)((cod1-cod2)+pos1+readLen2-pos2);
						if(BaseInsertSize<-1000||(BaseInsertSize<=500&&((pdist<=(2*BaseInsertSize))&&(pdist>=(-100))))||((pdist<=(3*BaseInsertSize/2))&&(pdist>=(BaseInsertSize/2))))
						{
							dist_sum+=pdist;
							p_cnt++;
							o_Pdist<<pdist<<endl;
							if(insert_sz_est_vt.size()<=102)
							{
								insert_sz_est_vt.push_back(pdist);
								if(insert_sz_est_vt.size()==101)
								{
									sort(insert_sz_est_vt.begin(),insert_sz_est_vt.end());
									BaseInsertSize=insert_sz_est_vt[50];
								}
							}

						}
					}
					else
					{
						//orientation errors
						if(flip_1==1&&(cod2>cod1))
						{
							
						}
						if(flip_2==0&&(cod2<cod1))
						{
							
						}
					
					}

				}
			}
			bool ignore=0;
			if(cont1!=cont2)
			{

				int pdist=-10000;

				
				if(MatePair&&LocalSearch)
				{
					struct stacked_bucket4 stacked_bucket;

					stacked_bucket.bktptr=*ptr1_d;
					if(flip_1d==0)
					{
						stacked_bucket.RightSearch=1;
					}
					else
					{
						stacked_bucket.RightSearch=0;
					}
					pdist= BFSearchDist4(ht, merge_ht, *ptr1_d,*ptr2_d,K_size, stacked_bucket,MaxDepth,MaxSearchLen);
					if(pdist>0)
					{
						ignore=1;
						o_Pdist_in<<(*ptr1_d)->kmer_info.contig_no<<" "<<(*ptr2_d)->kmer_info.contig_no<<" "<<pdist<<endl;
						inward_found++;
					}
					else
					{
						inward_not_found++;
					}
					stacked_bucket.RightSearch=!stacked_bucket.RightSearch;
					pdist= BFSearchDist4(ht, merge_ht, *ptr1_d,*ptr2_d,K_size, stacked_bucket,MaxDepth,MaxSearchLen);
					if(pdist>0)
					{
						ignore=1;
						o_Pdist_out<<(*ptr1_d)->kmer_info.contig_no<<" "<<(*ptr2_d)->kmer_info.contig_no<<" "<<pdist<<endl;
						outward_found++;
					}
					else
					{
						outward_not_found++;
					}

				}
				if(ignore==1)
				{continue;}
			
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
					//	if(0)//if(Tdist<0)
					//	{
					//		Tdist=0;///////////////
					//	}

						o_Pdist_L<<cont1 <<" "<<cont2<<" "<<d4<<endl;//" Library: "<<lib_no<<endl;
						//o_Pdist_R<<cont2 <<" "<<cont1<<" "<<d4<<" Library: "<<lib_no<<endl;

					}

//					contigs_info->scaffold_adjacency_left[cont1].push_back(adj_contig);

				}




			}



		}
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






bool BackCheckLoop_ctg(int new_ctg,int end_node, map<int,struct BFS_path_info_ctg > & Visited_Path )
{
	new_ctg=abs(new_ctg);
	end_node=abs(end_node);
	int cur_node=abs(end_node);
	int dep=Visited_Path[cur_node].depth;
	int pre_node;
	for(int l=dep;l>=1;--l)
	{

		pre_node=abs(Visited_Path[cur_node].last_ctg);
		
		if(new_ctg==cur_node)
		{
			return 1;
		}

		cur_node=pre_node;
	}

	return 0;

}


void BFSearchPathFinder(struct contigs_info *contigs_info,list<int> ctg_stack,map<int,list<int> > &dist_ctg,map<int,int > & unitig_dist,int dist_searched,map<int,vector<int> > &node_cov,int max_depth,bool GapClosingMode,bool MatePair)//
{	
		
	map<int, struct BFS_path_info_ctg > Visited_Path;
	map<int, int > stacked_nodes;
//	map<int, vector<int> > node_cov;

	int beg_node=ctg_stack.front();
	beg_node=abs(beg_node);
	int max_stack=300;
	int DepthTh=max_depth;//min(300/gap,20);
	//int LenTh=1000;
	int TipLenTh=100;
	
	int dist_searched_t=0;
	int DistVar=1000;
	int DistVarFactor=2;
	bool RIGHT=0;
	int new_node=ctg_stack.front();
	
	if(new_node>0)
	{
		stacked_nodes[new_node]=1;
	}
	else
	{
		stacked_nodes[abs(new_node)]=-1;
	}
	Visited_Path[abs(new_node)].cov=0;
	Visited_Path[abs(new_node)].depth=1;
	Visited_Path[abs(new_node)].len=dist_searched;
	Visited_Path[abs(new_node)].last_ctg=0;

	while(ctg_stack.size()>0)
	{
		if(((int)ctg_stack.size())>max_stack)
		{
			break;
		}

		new_node=ctg_stack.front();
		ctg_stack.pop_front();
	
		
		if(new_node>0)
		{
			RIGHT=1;
		}
		else
		{
			RIGHT=0;
		}

		if(Visited_Path.count(abs(new_node)))
		{
			if(Visited_Path[abs(new_node)].depth>DepthTh)//||Visited_Path[abs(new_node)].len>LenTh)
			{					
				continue;
			}
		}

		if(RIGHT)
		{				
			map<int, struct adjacent_contig_info>::iterator ctg_it;
			int nRB=0;
			for(ctg_it=contigs_info->contig_adjacency_right[new_node].begin();ctg_it!=contigs_info->contig_adjacency_right[new_node].end();)
			{
				map<int, struct adjacent_contig_info>::iterator n_ctg_it=ctg_it;
				n_ctg_it++;
				int ctg_no=ctg_it->first;
				nRB++;				
				ctg_it=n_ctg_it;			
			}
			stacked_nodes[new_node]=1+nRB;
			
			
			if(nRB==0)
			{
				stacked_nodes[new_node]=2;
				
				//tip end reached
				continue;

			}
			
			//so at least one ctg on the right

			dist_searched_t=Visited_Path[abs(new_node)].len;
			map<int,list<int> >::iterator dist_ctg_it,dist_ctg_it_n;
			map<int,int> marked_ctgs;



			for(dist_ctg_it=dist_ctg.begin();dist_ctg_it!=dist_ctg.end();)
			{
				dist_ctg_it_n=dist_ctg_it;
				dist_ctg_it_n++;
						
						
				list<int>::iterator lit,lit2,n_lit;
				for(lit=dist_ctg_it->second.begin();lit!=dist_ctg_it->second.end();)
				{
					lit2=lit;
					lit2++;
					n_lit=lit2;
					n_lit++;
					int dist0=*lit2;
						
					int DistVar2=DistVar;
													
					if(MatePair)//dist0>3000)
					{
						DistVar2=20000;//dist0;//dist0/DistVarFactor;
					}


					if(dist_ctg_it->first<=(dist_searched-DistVar2))
					{
					//	dist_ctg.erase(dist_ctg_it);
					//	dist_ctg_it->second.erase(lit2);
						//dist_ctg_it->second.erase(lit);
						lit=n_lit;
						continue;
					}

							
					if((dist_ctg_it->first>dist_searched-DistVar2)&&(dist_ctg_it->first<dist_searched+DistVar2))
					{
						marked_ctgs[abs(*lit)]=dist_ctg_it->first-dist_searched;
						lit=n_lit;
						continue;
					}
					lit=n_lit;
						
				}
				if(dist_ctg_it->second.size()==0)
				{
					dist_ctg.erase(dist_ctg_it);
				}
		
				dist_ctg_it=dist_ctg_it_n;

					
			}

			
			for(ctg_it=contigs_info->contig_adjacency_right[new_node].begin();ctg_it!=contigs_info->contig_adjacency_right[new_node].end(); )
			{
				
				map<int, struct adjacent_contig_info>::iterator n_ctg_it=ctg_it;
				n_ctg_it++;
				
				int nxt_ctg=ctg_it->first;

				int is_marked=0,is_unitig;

				if(marked_ctgs.count(abs(nxt_ctg)))
				{
					is_marked=1;
				}

				if(1)//
				{
					// not in stack
					if(stacked_nodes.count(abs(nxt_ctg))==0)//&&contigs_info->c_info_vt[abs(nxt_ctg)].removed==0)
					{
						
						Visited_Path[abs(nxt_ctg)].cov=(Visited_Path[abs(new_node)].cov+is_marked);
						
						Visited_Path[abs(nxt_ctg)].depth=(Visited_Path[abs(new_node)].depth+1);
						Visited_Path[abs(nxt_ctg)].len=(Visited_Path[abs(new_node)].len+contigs_info->contig_sz_vt[abs(nxt_ctg)]);
						Visited_Path[abs(nxt_ctg)].last_ctg=abs(new_node);
			
						if(unitig_dist.count(abs(nxt_ctg)))//is_mark
						{
							if(Visited_Path[abs(new_node)].cov==0)//GapClosingMode&&
							{
								node_cov[(nxt_ctg)].clear();
								node_cov[(nxt_ctg)].push_back(1);//cov
								int dep=Visited_Path[abs(nxt_ctg)].depth;

								int p_ctg=nxt_ctg;
								while(p_ctg!=0)
								{
//									
									dep--;
									if(dep<0)
									{
										//cout<<"w"<<endl;
										//node_cov[(nxt_ctg)].clear();
										node_cov.erase(nxt_ctg);
										ctg_it=n_ctg_it;
										break;
									}
									
									node_cov[(nxt_ctg)].push_back(p_ctg);
									p_ctg=Visited_Path[abs(p_ctg)].last_ctg;
									if(stacked_nodes[abs(p_ctg)]<0)
									{
										p_ctg=-p_ctg;
									}
								}
								if(dep>=0)
								{
									reverse((node_cov[(nxt_ctg)].begin()+1),node_cov[(nxt_ctg)].end());
									ctg_it=n_ctg_it;
									continue;
								}
							}
							if((GapClosingMode==0)&&Visited_Path[abs(new_node)].cov>=1)
							{
								node_cov[(nxt_ctg)].clear();
								node_cov[(nxt_ctg)].push_back(Visited_Path[abs(new_node)].cov+1);
								node_cov.erase(new_node);
															
								int p_ctg=nxt_ctg;

								int dep=Visited_Path[abs(nxt_ctg)].depth;

								while(p_ctg!=0)
								{

									dep--;
									if(dep<0)
									{
										//cout<<"w"<<endl;
										
										node_cov.erase(nxt_ctg);
										ctg_it=n_ctg_it;
										break;
									}
									
									node_cov[(nxt_ctg)].push_back(p_ctg);
									p_ctg=Visited_Path[abs(p_ctg)].last_ctg;
									if(stacked_nodes[abs(p_ctg)]<0)
									{
										p_ctg=-p_ctg;
									}
								}
								if(dep>=0)
								{
									reverse(node_cov[(nxt_ctg)].begin()+1,node_cov[(nxt_ctg)].end());

									ctg_it=n_ctg_it;
									continue;
								}
							}
							
						}

						if(nxt_ctg<0)
						{
							stacked_nodes[abs(nxt_ctg)]=-1;
						}
						else
						{
							stacked_nodes[abs(nxt_ctg)]=1;
						}

						ctg_stack.push_back(nxt_ctg);
						if(((int)ctg_stack.size())>max_stack)
						{
							break;
						}

						ctg_it=n_ctg_it;
						continue;

					}
					else
					{
						if((stacked_nodes[abs(nxt_ctg)]>0&&nxt_ctg>0)||(stacked_nodes[abs(nxt_ctg)]<0&&nxt_ctg<0))
						{
							//backtrack if the same direction is found
							if(((Visited_Path[new_node].cov+is_marked)<=Visited_Path[abs(nxt_ctg)].cov)||(BackCheckLoop_ctg( abs(nxt_ctg),abs(new_node),Visited_Path)==1))
							{
								ctg_it=n_ctg_it;
								continue;
								
							}
							else
							{
								//backtrack the original path
							
								Visited_Path[abs(nxt_ctg)].cov=(Visited_Path[abs(new_node)].cov+is_marked);
								Visited_Path[abs(nxt_ctg)].depth=Visited_Path[abs(new_node)].depth+1;
								Visited_Path[abs(nxt_ctg)].len=(Visited_Path[abs(new_node)].len+contigs_info->contig_sz_vt[abs(nxt_ctg)]);
								Visited_Path[abs(nxt_ctg)].last_ctg=abs(new_node);


								if(unitig_dist.count(abs(nxt_ctg)))//is_mark
								{
									if(Visited_Path[abs(new_node)].cov==0)//GapClosingMode&&
									{
										node_cov[(nxt_ctg)].clear();
										node_cov[(nxt_ctg)].push_back(1);
										
										int p_ctg=nxt_ctg;
										int dep=Visited_Path[abs(nxt_ctg)].depth;

										while(p_ctg!=0)
										{
											
										//	if(sc_cnt==4866)
											//{cout<<"q";}
											dep--;
											if(dep<0)
											{
												//cout<<"w"<<endl;
												//node_cov[(nxt_ctg)].clear();
												node_cov.erase(nxt_ctg);
												ctg_it=n_ctg_it;
												break;
											}
									
											node_cov[(nxt_ctg)].push_back(p_ctg);
											p_ctg=Visited_Path[abs(p_ctg)].last_ctg;
											if(stacked_nodes[abs(p_ctg)]<0)
											{
												p_ctg=-p_ctg;
											}
										}
										if(dep>=0)
										{
											reverse(node_cov[(nxt_ctg)].begin()+1,node_cov[(nxt_ctg)].end());	
											ctg_it=n_ctg_it;
											continue;
										}
									}
									if((GapClosingMode==0)&&Visited_Path[abs(new_node)].cov>=1)
									{
										node_cov[(nxt_ctg)].clear();
										node_cov[(nxt_ctg)].push_back(Visited_Path[abs(new_node)].cov+1);
										node_cov.erase(new_node);

										int p_ctg=nxt_ctg;
										int dep=Visited_Path[abs(nxt_ctg)].depth;

										while(p_ctg!=0)
										{
										//	if(sc_cnt==4866)
											//{cout<<"w";}
											dep--;
											if(dep<0)
											{
												//cout<<"w"<<endl;
												//node_cov[(nxt_ctg)].clear();
												node_cov.erase(nxt_ctg);
												ctg_it=n_ctg_it;
												break;;
											}
									
											node_cov[(nxt_ctg)].push_back(p_ctg);
											p_ctg=Visited_Path[abs(p_ctg)].last_ctg;
											if(stacked_nodes[abs(p_ctg)]<0)
											{
												p_ctg=-p_ctg;
											}
										}
										if(dep>=0)
										{
											reverse(node_cov[(nxt_ctg)].begin()+1,node_cov[(nxt_ctg)].end());
											ctg_it=n_ctg_it;
											continue;
										}
									}
									//node_cov[abs(nxt_ctg)]=Visited_Path[abs(new_node)].cov+is_marked;
								}
								
							}

						}
						else
						{

							//don't do anything,since both strands are visited.

						}


					}

				}
				ctg_it=n_ctg_it;
				
			}
		}
		else
		{
			new_node=abs(new_node);

			map<int, struct adjacent_contig_info>::iterator ctg_it;
			int nLB=0;
			for(ctg_it=contigs_info->contig_adjacency_left[new_node].begin();ctg_it!=contigs_info->contig_adjacency_left[new_node].end();)
			{
				
				map<int, struct adjacent_contig_info>::iterator n_ctg_it=ctg_it;
				n_ctg_it++;
				int ctg_no=ctg_it->first;
				
				nLB++;
				
				ctg_it=n_ctg_it;
			}
			stacked_nodes[new_node]=-1-nLB;
			
			if(nLB==0)
			{

				stacked_nodes[abs(new_node)]=-2;
			
				continue;
			}



			
			dist_searched_t=Visited_Path[abs(new_node)].len;
			map<int,list<int> >::iterator dist_ctg_it,dist_ctg_it_n;
			map<int,int> marked_ctgs;

			for(dist_ctg_it=dist_ctg.begin();dist_ctg_it!=dist_ctg.end();)
			{
				dist_ctg_it_n=dist_ctg_it;
				dist_ctg_it_n++;
						
				
						
				list<int>::iterator lit,lit2,n_lit;
				for(lit=dist_ctg_it->second.begin();lit!=dist_ctg_it->second.end();)
				{
					lit2=lit;
					lit2++;
					n_lit=lit2;
					n_lit++;
					int dist0=*lit2;
						
					int DistVar2=DistVar;
													
					if(MatePair)//dist0>3000)
					{
						DistVar2=20000;//dist0;//dist0/DistVarFactor;
					}


					if(dist_ctg_it->first<=(dist_searched-DistVar2))
					{
					//	dist_ctg.erase(dist_ctg_it);
					//	dist_ctg_it->second.erase(lit2);
					//	dist_ctg_it->second.erase(lit);
						lit=n_lit;
						continue;
					}

							
					if((dist_ctg_it->first>dist_searched-DistVar2)&&(dist_ctg_it->first<dist_searched+DistVar2))
					{
						marked_ctgs[abs(*lit)]=dist_ctg_it->first-dist_searched;
						lit=n_lit;
						continue;
					}
					lit=n_lit;
						
				}
				if(dist_ctg_it->second.size()==0)
				{
					dist_ctg.erase(dist_ctg_it);
				}



						
				dist_ctg_it=dist_ctg_it_n;

					
			}




			for(ctg_it=contigs_info->contig_adjacency_left[abs(new_node)].begin();ctg_it!=contigs_info->contig_adjacency_left[abs(new_node)].end();)
			{

				map<int, struct adjacent_contig_info>::iterator n_ctg_it=ctg_it;
				n_ctg_it++;
				int nxt_ctg=ctg_it->first;
				int is_marked=0;

				if(marked_ctgs.count(abs(nxt_ctg)))
				{
					is_marked=1;
				}
				
				if(1)//contigs_info->c_info_vt[abs(nxt_ctg)].removed==0)
				{
					
					if(stacked_nodes.count(abs(nxt_ctg))==0)//&&contigs_info->c_info_vt[abs(nxt_ctg)].removed==0)
					{
						
						Visited_Path[abs(nxt_ctg)].cov=(Visited_Path[abs(new_node)].cov+is_marked);						
						Visited_Path[abs(nxt_ctg)].depth=(Visited_Path[abs(new_node)].depth+1);
						Visited_Path[abs(nxt_ctg)].len=(Visited_Path[abs(new_node)].len+contigs_info->contig_sz_vt[abs(nxt_ctg)]);
						Visited_Path[abs(nxt_ctg)].last_ctg=abs(new_node);//abs

						if(unitig_dist.count(abs(nxt_ctg)))//is_mark
						{
							//node_cov[abs(nxt_ctg)]=Visited_Path[abs(new_node)].cov+is_marked;
							//if(GapClosingMode&&Visited_Path[abs(new_node)].cov==0)
							if(GapClosingMode&&node_cov.count(-(nxt_ctg))==0)
							{
								node_cov[-(nxt_ctg)].clear();
								node_cov[-(nxt_ctg)].push_back(1);
																
								int p_ctg=-nxt_ctg;
								int dep=Visited_Path[abs(nxt_ctg)].depth;

								while(p_ctg!=0)
								{
											
//									if(sc_cnt==4866)
	//										{cout<<"e";}
									dep--;
									if(dep<0)
									{
										//cout<<"w"<<endl;
										//node_cov[(-nxt_ctg)].clear();
										node_cov.erase(-nxt_ctg);
										ctg_it=n_ctg_it;
										break;
									}
									
									node_cov[-(nxt_ctg)].push_back(p_ctg);
									p_ctg=Visited_Path[abs(p_ctg)].last_ctg;
									if(stacked_nodes[abs(p_ctg)]<0)
									{
										p_ctg=-p_ctg;
									}
								}
								if(dep>=0)
								{
									reverse(node_cov[-(nxt_ctg)].begin()+1,node_cov[-(nxt_ctg)].end());
									ctg_it=n_ctg_it;
									continue;
								}
								
							}
							if((GapClosingMode==0)&&Visited_Path[abs(new_node)].cov>=1)
							{
								node_cov.erase(-new_node);
								node_cov[-(nxt_ctg)].clear();
								node_cov[-(nxt_ctg)].push_back(Visited_Path[abs(new_node)].cov+1);

								
								int p_ctg=-nxt_ctg;
								
								int dep=Visited_Path[abs(nxt_ctg)].depth;

								while(p_ctg!=0)
								{
//									if(sc_cnt==4866)
	//										{cout<<"r";}
											
									dep--;
									if(dep<0)
									{
										//cout<<"w"<<endl;
										//node_cov[-(nxt_ctg)].clear();
										node_cov.erase(-nxt_ctg);
										ctg_it=n_ctg_it;
										break;
									}
									node_cov[-(nxt_ctg)].push_back(p_ctg);
									p_ctg=Visited_Path[abs(p_ctg)].last_ctg;
									if(stacked_nodes[abs(p_ctg)]<0)
									{
										p_ctg=-p_ctg;
									}
								}
								if(dep>=0)
								{
									reverse(node_cov[-(nxt_ctg)].begin()+1,node_cov[-(nxt_ctg)].end());
									ctg_it=n_ctg_it;
									continue;
								}
								
							}
							
						}
						
						if(nxt_ctg<0)
						{
							stacked_nodes[abs(nxt_ctg)]=1;
						}
						else
						{
							stacked_nodes[abs(nxt_ctg)]=-1;
						}
						ctg_stack.push_back(-(nxt_ctg));
						if(((int)ctg_stack.size())>max_stack)
						{
							break;
						}
						ctg_it=n_ctg_it;
						continue;

					}
					else
					{
						
						if((stacked_nodes[abs(nxt_ctg)]<0&&nxt_ctg>0)||(stacked_nodes[abs(nxt_ctg)]>0&&nxt_ctg<0))
						{
							//backtrack if the same direction is found
							if((Visited_Path[abs(new_node)].cov+is_marked<=Visited_Path[abs(nxt_ctg)].cov)||(BackCheckLoop_ctg( abs(nxt_ctg),abs(new_node),Visited_Path)==1))//
							{
								
								//backtrack the current path
								
								ctg_it=n_ctg_it;
								continue;
								
							
							}
							else
							{
								//backtrack the original path
						
								Visited_Path[abs(nxt_ctg)].cov=(Visited_Path[abs(new_node)].cov+is_marked);
	
								Visited_Path[abs(nxt_ctg)].depth=Visited_Path[abs(new_node)].depth+1;
								Visited_Path[abs(nxt_ctg)].len=(Visited_Path[abs(new_node)].len+contigs_info->contig_sz_vt[abs(nxt_ctg)]);

								Visited_Path[abs(nxt_ctg)].last_ctg=abs(new_node);


								if(unitig_dist.count(abs(nxt_ctg)))//is_mark
								{
							//		node_cov[abs(nxt_ctg)]=Visited_Path[abs(new_node)].cov+is_marked;
									//if(GapClosingMode&&Visited_Path[abs(new_node)].cov==0)
									if(GapClosingMode&&node_cov.count(-(nxt_ctg))==0)
									{
										node_cov[-(nxt_ctg)].clear();
										node_cov[-(nxt_ctg)].push_back(1);
								
										
										int p_ctg=-nxt_ctg;	
										int dep=Visited_Path[abs(nxt_ctg)].depth;

										while(p_ctg!=0)
										{
//											if(sc_cnt==4866)
	//										{cout<<"t";}
											
											dep--;
											if(dep<0)
											{
												//cout<<"w"<<endl;
												//node_cov[-(nxt_ctg)].clear();
												node_cov.erase(-nxt_ctg);
												ctg_it=n_ctg_it;
												break;
											}
											
											node_cov[-(nxt_ctg)].push_back(p_ctg);
											p_ctg=Visited_Path[abs(p_ctg)].last_ctg;
											if(stacked_nodes[abs(p_ctg)]<0)
											{
												p_ctg=-p_ctg;
											}
										}
										if(dep>=0)
										{
											reverse(node_cov[-(nxt_ctg)].begin()+1,node_cov[-(nxt_ctg)].end());
											ctg_it=n_ctg_it;
											continue;
										}
									}
									if((GapClosingMode==0)&&Visited_Path[abs(new_node)].cov>=1)
									{
									
										node_cov.erase(-new_node);
										node_cov[-(nxt_ctg)].clear();
										node_cov[-(nxt_ctg)].push_back(Visited_Path[abs(new_node)].cov+1);

										
										int p_ctg=-nxt_ctg;
										int dep=Visited_Path[abs(nxt_ctg)].depth;

										while(p_ctg!=0)
										{

											dep--;
											if(dep<0)
											{
												//cout<<"w"<<endl;
												//node_cov[-(nxt_ctg)].clear();
												node_cov.erase(-nxt_ctg);
												ctg_it=n_ctg_it;
												break;
											}
										
											
											node_cov[-(nxt_ctg)].push_back(p_ctg);
											p_ctg=Visited_Path[abs(p_ctg)].last_ctg;
											if(stacked_nodes[abs(p_ctg)]<0)
											{
												p_ctg=-p_ctg;
											}
										}
										if(dep>=0)
										{
											reverse(node_cov[-(nxt_ctg)].begin()+1,node_cov[-(nxt_ctg)].end());
											ctg_it=n_ctg_it;
											continue;
										}
								
									}
							
								}


								ctg_it=n_ctg_it;
								continue;

							}

						}
						else
						{

							//don't do anything,since both strands are visited.

						}




					}

				}
				ctg_it=n_ctg_it;

			}

		}

	}

}


void ResolvingRepeatsPE(vector<int> &insert_sz_vt,vector<string>& filenames_vt,struct contigs_info * contigs_info,string ContigFilename,int LinkCovTh,int UniqueLenTh,int ExpCov)
{
	bool MatePair=0;
	bool Iter_Scaffold=0;
	map<int,int> Unitigs;
	int LinkCovTh0=1;
	int LinkCovTh1=max(LinkCovTh/2,2);

	string pe_name="InsertSizeEst.txt";
	string uniqueness_name="UniqueContigs.txt";
	if(ContigFilename=="SuperContigs.txt")
	{
		pe_name="InsertSizeEst.txt";//Iterated
		uniqueness_name="UniqueSuperContigs.txt";
		Iter_Scaffold=1;
	}
	if(insert_sz_vt.size()>0)
	{
		MatePair=1;
		for(int i=0;i!=insert_sz_vt.size();++i)
		{
			if(insert_sz_vt[i]<=10000)
			{MatePair=0;}
		}
		
	}
	ifstream in_pe_info(pe_name.c_str());
	string str,s;
	int ins_est;
	int max_dep=9;
	int DistVar=500;

	int DistVar0=500;

	if(MatePair)
	{
		DistVar=10000;
		DistVar0=5000;
	}
	int DistVarFactor=2;
	
	vector<int> ins_est_vt;
	while(in_pe_info>>str>>ins_est)
	{
		ins_est_vt.push_back(ins_est);
	}
	int lib_cnt=ins_est_vt.size();

	for(int i=1;i<=lib_cnt;++i)
	{
		if(ins_est_vt[i-1]<0)
		{
			continue;
		}
		char sc_name[1000];
		if(1)//Iter_Scaffold==0)
		{
			sprintf(sc_name,"Pdist_L_lib_%d.txt",i);
		}
		else
		{
			sprintf(sc_name,"PdistIter_L_lib_%d.txt",i);
		}
		ifstream in_sc_l(sc_name);
		contigs_info->scaffold_adjacency_left.resize(contigs_info->total_contigs+1);
		contigs_info->scaffold_adjacency_right.resize(contigs_info->total_contigs+1);
		contigs_info->c_info_vt.resize(contigs_info->total_contigs+1);

		cout<<"Loading adjacent info. left. "<<i<<endl;
		int cont1,cont2,dist;
		int nr=0;
		
		while(in_sc_l>>cont1>>cont2>>dist)
		{
			nr++;
			if(nr%1000000==0)
			{
				cout<<nr<<endl;
			}
		//	cout<<nr<<endl;
			dist+=ins_est_vt[i-1];
			if(dist<-200||(dist>(ins_est_vt[i-1]*3/2)))//double check
			{continue;}

	
			if(1)//contigs_info->contig_sz_vt[cont1]>=UniqueLenTh)
			{

				if(contigs_info->scaffold_adjacency_left[cont1].count(cont2)==0)
				{
					contigs_info->scaffold_adjacency_left[cont1][cont2].cov=0;
					contigs_info->scaffold_adjacency_left[cont1][cont2].dist_sum=0;

				}

				if( contigs_info->scaffold_adjacency_left[cont1][cont2].cov<20)//contigs_info->contig_sz_vt[abs(cont1)]>LengthTh&&contigs_info->contig_sz_vt[abs(cont2)]>LengthTh&&
				{
					contigs_info->scaffold_adjacency_left[cont1][cont2].dist_sum+=dist;
					contigs_info->scaffold_adjacency_left[cont1][cont2].cov+=1;
				
				}
			}

			if(1)//contigs_info->contig_sz_vt[abs(cont2)]>=UniqueLenTh)
			{

				if(cont2>0)
				{
					if(contigs_info->scaffold_adjacency_right[cont2].count(cont1)==0)
					{
						contigs_info->scaffold_adjacency_right[cont2][cont1].cov=0;
						contigs_info->scaffold_adjacency_right[cont2][cont1].dist_sum=0;
					}

					if(contigs_info->scaffold_adjacency_right[cont2][cont1].cov<20)
					{

						contigs_info->scaffold_adjacency_right[cont2][cont1].dist_sum+=dist;
						contigs_info->scaffold_adjacency_right[cont2][cont1].cov+=1;
					}
				}
				else
				{
					if(contigs_info->scaffold_adjacency_left[-cont2].count(-cont1)==0)
					{
						contigs_info->scaffold_adjacency_left[-cont2][-cont1].cov=0;
						contigs_info->scaffold_adjacency_left[-cont2][-cont1].dist_sum=0;
					}

					if(contigs_info->scaffold_adjacency_left[-cont2][-cont1].cov<20)
					{
					contigs_info->scaffold_adjacency_left[-cont2][-cont1].dist_sum+=dist;
					contigs_info->scaffold_adjacency_left[-cont2][-cont1].cov+=1;
					}
				}
			}

		}
		in_sc_l.clear();
		in_sc_l.close();



		cout<<"Loading adjacent info. right."<<i<<endl;

		if(MatePair==0)
		{
		sprintf(sc_name,"Pdist_R_lib_%d.txt",i);
		}
		else
		{
			sprintf(sc_name,"PdistIter_R_lib_%d.txt",i);
		}
		ifstream in_sc_r(sc_name);
		 nr=0;
		while(in_sc_r>>cont1>>cont2>>dist)
		{


		//	if((abs(cont1)==16&&abs(cont2)==9)||(abs(cont1)==9&&abs(cont2)==16))
			//{cout<<"";}
			nr++;
			if(nr%1000000==0)
			{
				cout<<nr<<endl;
			}

			dist+=ins_est_vt[i-1];
			
			if(dist<-200||(dist>(ins_est_vt[i-1]*3/2)))//double check
			{continue;}

			if(1)//contigs_info->contig_sz_vt[cont1]>=UniqueLenTh)
			{
				if(contigs_info->scaffold_adjacency_right[cont1].count(cont2)==0)
				{
					contigs_info->scaffold_adjacency_right[cont1][cont2].cov=0;
					contigs_info->scaffold_adjacency_right[cont1][cont2].dist_sum=0;

				}

				if(contigs_info->scaffold_adjacency_right[cont1][cont2].cov<20)		//contigs_info->contig_sz_vt[abs(cont1)]>LengthTh&&contigs_info->contig_sz_vt[abs(cont2)]>LengthTh&&
				{
					contigs_info->scaffold_adjacency_right[cont1][cont2].dist_sum+=dist;
					contigs_info->scaffold_adjacency_right[cont1][cont2].cov+=1;
					
				}

			}


			if(1)//contigs_info->contig_sz_vt[abs(cont2)]>=UniqueLenTh)
			{
				if(cont2>0)
				{
					if(contigs_info->scaffold_adjacency_left[cont2].count(cont1)==0)
					{
						contigs_info->scaffold_adjacency_left[cont2][cont1].cov=0;
						contigs_info->scaffold_adjacency_left[cont2][cont1].dist_sum=0;
					}
					if(contigs_info->scaffold_adjacency_left[cont2][cont1].cov<20)		//contigs_info->contig_sz_vt[abs(cont1)]>LengthTh&&contigs_info->contig_sz_vt[abs(cont2)]>LengthTh&&
					{
						contigs_info->scaffold_adjacency_left[cont2][cont1].dist_sum+=dist;
						contigs_info->scaffold_adjacency_left[cont2][cont1].cov+=1;
					}
				}
				else
				{
					if(contigs_info->scaffold_adjacency_right[-cont2].count(-cont1)==0)
					{
						contigs_info->scaffold_adjacency_right[-cont2][-cont1].cov=0;
						contigs_info->scaffold_adjacency_right[-cont2][-cont1].dist_sum=0;
					}
					if(contigs_info->scaffold_adjacency_right[-cont2][-cont1].cov<20)
					{
						contigs_info->scaffold_adjacency_right[-cont2][-cont1].dist_sum+=dist;
						contigs_info->scaffold_adjacency_right[-cont2][-cont1].cov+=1;
					}

				}
			}
		}
		in_sc_r.clear();
		in_sc_r.close();



	}

	cout<<"Removing weak links..."<<endl;
	
	for(int i=1;i<=contigs_info->total_contigs;++i)
	{

		map<int,scaffold_contig_info>::iterator s_adj_it,s_adj_it_n;
		map<int,scaffold_contig_info> temp_map=contigs_info->scaffold_adjacency_left[i];
		for(s_adj_it=temp_map.begin();s_adj_it!=temp_map.end();)
		{
			s_adj_it_n=s_adj_it;
			s_adj_it_n++;
			if(s_adj_it->second.cov<=LinkCovTh0)
			{
				temp_map.erase(s_adj_it);
			}
			s_adj_it=s_adj_it_n;
		}

		contigs_info->scaffold_adjacency_left[i].clear();
		for(s_adj_it=temp_map.begin();s_adj_it!=temp_map.end();++s_adj_it)
		{
			contigs_info->scaffold_adjacency_left[i][s_adj_it->first]=s_adj_it->second;
		}

		temp_map=contigs_info->scaffold_adjacency_right[i];
		for(s_adj_it=temp_map.begin();s_adj_it!=temp_map.end();)
		{
			s_adj_it_n=s_adj_it;
			s_adj_it_n++;
			if(s_adj_it->second.cov<=LinkCovTh0)
			{
				temp_map.erase(s_adj_it);
			}
			s_adj_it=s_adj_it_n;
		}

		contigs_info->scaffold_adjacency_right[i].clear();
		for(s_adj_it=temp_map.begin();s_adj_it!=temp_map.end();++s_adj_it)
		{
			contigs_info->scaffold_adjacency_right[i][s_adj_it->first]=s_adj_it->second;
		}



	}
	
	cout<<"Done."<<endl;
	
	string ContigsHPName="Contigs_HP.txt";
	if(Iter_Scaffold==1)
	{
		ContigsHPName="SuperContigs_HP.txt";
	}

	ifstream in_contigs(ContigFilename.c_str()),in_contigs_HP(ContigsHPName.c_str());
	string s1,cont_s;
	int codB,codE;
	contigs_info->contigs_str.resize(contigs_info->total_contigs+1);
	contigs_info->contigs_hp_b.resize(contigs_info->total_contigs+1);
	contigs_info->contigs_hp_e.resize(contigs_info->total_contigs+1);

	int n_ctgs=0;
	while(getline(in_contigs,s1))
	{
		n_ctgs++;
		getline(in_contigs,cont_s);
		if(cont_s[cont_s.size()-1]=='\r'||cont_s[cont_s.size()-1]=='\n')
		{cont_s.resize(cont_s.size()-1);}
		in_contigs_HP>>codB>>codE;
		contigs_info->contigs_hp_b[n_ctgs]=codB;
		contigs_info->contigs_hp_e[n_ctgs]=codE;
		contigs_info->contigs_str[n_ctgs]=cont_s;
		
	}

	int K_size=contigs_info->K_size;
	
	string sc_info_name="SuperContigs_info.txt",sc_name="SuperContigs.txt";

	if(Iter_Scaffold)
	{
		sc_info_name="SuperContigs_info.txt",sc_name="SuperContigs.txt";//sc_info_name="IteratedSuperContigs_info.txt",sc_name="IteratedSuperContigs.txt";
	}
	/*
	if(MatePair==1)
	{
		sc_info_name="Scaffolds_info.txt";
		sc_name="Scaffolds.txt";
	}
	*/

	ofstream o_sc_info(sc_info_name.c_str());
	ofstream o_sc(sc_name.c_str());
	string sc_cov_name="SuperContigs_Cov.txt";
	ofstream osc_cov;
	if(Iter_Scaffold)
	{
		sc_cov_name="SuperContigs_Cov.txt";//sc_cov_name="IteratedSuperContigs_Cov.txt";
	}
	
	osc_cov.open(sc_cov_name.c_str());

	int sc_cnt=0;

	int tot_ctgs=contigs_info->total_contigs;
//	contigs_info->c_info_vt.resize(tot_ctgs+1);
	cout<<"Building super contigs..."<<endl;
	//ofstream o_ContigGap_info("ContigGap_info.txt");
	//contigs_info->scaffolds.resize(1);
	//contigs_info->gaps_in_scaffolds.resize(1);
		
	int break_points=0;
	int ambiguous_points=0;
	int gap_points=0;

	ofstream o_unique(uniqueness_name.c_str());

	for(int i=1;i<=tot_ctgs;++i)
	{

		/////////////////////////////////////////////////////////////////////////////////Record gaps
		//...
		/////////////////////////////////////////////////////////////////////////////////Record gaps


		//record branches
		int EdgeCovTh=max(1,LinkCovTh/10);
		int LeftBranches=0,RightBranches=0;
		 map<int,struct adjacent_contig_info>::iterator it;
		if(contigs_info->contig_adjacency_left[i].size()>0)
		for(it=contigs_info->contig_adjacency_left[i].begin();it!=contigs_info->contig_adjacency_left[i].end();++it)
		{
			if(it->second.cov>EdgeCovTh)
			{
				LeftBranches++;
			}
		}
		if(contigs_info->contig_adjacency_right[i].size()>0)
		for(it=contigs_info->contig_adjacency_right[i].begin();it!=contigs_info->contig_adjacency_right[i].end();++it)
		{
			if(it->second.cov>EdgeCovTh)
			{
				RightBranches++;
			}
		}
		////.record branches finished



		////critical importance, repeat mask


		if(LeftBranches<=1&&RightBranches<=1&&contigs_info->contig_sz_vt[i]>=UniqueLenTh&&contigs_info->cov_vt[i]<(ExpCov*3/2))
		{
			contigs_info->c_info_vt[i].unique=1;
			contigs_info->c_info_vt[i].loc_unique=1;
		}
		else
		{
			contigs_info->c_info_vt[i].unique=0;
			contigs_info->c_info_vt[i].loc_unique=0;
		}

	}

	////critical importance, repeat mask strategy2

	for(int i=1;i<=tot_ctgs;++i)
	{
	
		bool Unique=1;
		map<int,int> dist_ctg;
		map<int,struct scaffold_contig_info>::iterator ctg_info_it1;//ctg_info_it2;
		for(ctg_info_it1=contigs_info->scaffold_adjacency_left[i].begin();ctg_info_it1!=contigs_info->scaffold_adjacency_left[i].end();++ctg_info_it1)
		{
			if(contigs_info->c_info_vt[abs(ctg_info_it1->first)].unique==0||ctg_info_it1->second.cov<=LinkCovTh)
			{continue;}
			int dist1=ctg_info_it1->second.dist_sum/ctg_info_it1->second.cov;
			if(dist1>1000)
			{continue;}
			if(dist_ctg.count(dist1)==0)
			{
				dist_ctg[dist1]=abs(ctg_info_it1->first);
			}
			else
			{Unique=0;break;}
		}


		map<int,int>::iterator dist_it1,dist_it2;
		
		for(dist_it1=dist_ctg.begin();dist_it1!=dist_ctg.end();++dist_it1)
		{
			dist_it2=dist_it1;
			dist_it2++;
			if(dist_it2==dist_ctg.end())
			{
				break;
			}

			if((dist_it1->first+contigs_info->contig_sz_vt[dist_it1->second]-DistVar0)>dist_it2->first)
			{
				Unique=0;
			}

		}

		dist_ctg.clear();

		for(ctg_info_it1=contigs_info->scaffold_adjacency_right[i].begin();ctg_info_it1!=contigs_info->scaffold_adjacency_right[i].end();++ctg_info_it1)
		{
			if(contigs_info->c_info_vt[abs(ctg_info_it1->first)].unique==0||ctg_info_it1->second.cov<=LinkCovTh)
			{continue;}
			int dist1=ctg_info_it1->second.dist_sum/ctg_info_it1->second.cov;
			if(dist1>1000)
			{continue;}
			if(dist_ctg.count(dist1)==0)
			{
				dist_ctg[dist1]=abs(ctg_info_it1->first);
			}
			else
			{Unique=0;break;}
		}


		
		for(dist_it1=dist_ctg.begin();dist_it1!=dist_ctg.end();++dist_it1)
		{
			dist_it2=dist_it1;
			dist_it2++;
			if(dist_it2==dist_ctg.end())
			{
				break;
			}

			if((dist_it1->first+contigs_info->contig_sz_vt[dist_it1->second]-DistVar0)>dist_it2->first)
			{
				Unique=0;
			}



		}
		if(Unique==0)
		{
			contigs_info->c_info_vt[i].unique=0;
			contigs_info->c_info_vt[i].loc_unique=0;
		}


		if(contigs_info->c_info_vt[i].unique)
		{
			o_unique<<i<<endl;
			Unitigs[i]=1;
		}

//////////////////////////////////////////////////////

	}


	
	for(int i=0;i<tot_ctgs;++i )
	{
		int contig_no=(contigs_info->LengthRank[i])-contigs_info->contig_sz_vt.begin();

		if(contigs_info->c_info_vt[contig_no].used==1||contigs_info->c_info_vt[contig_no].unique!=1)
		{
			continue;
		}

		
		sc_cnt++;
		//cout<<sc_cnt<<endl;


		map<int,int> unitig_dist;

		unitig_dist[contig_no]=0;
		vector<int> left_ctgs,right_ctgs,t_left_ctgs,t_right_ctgs;	
		map<int,int> ctgs_for_second_round;		
		//make it 2 rounds to close the gaps 
		for(int round=1;round<=1;++round)
		{
			unitig_dist.clear();
			contigs_info->c_info_vt[contig_no].used=1;
			left_ctgs.clear();
			right_ctgs.clear();
			t_left_ctgs.clear();
			t_right_ctgs.clear();
			for(int it=1;it<=2;++it)
			{
				int beg_ctg=contig_no;
				int current_ctg=beg_ctg;
				int last_unique_ctg=beg_ctg;
				int nxt_ctg;
				map<int,int> searched_ctgs;
				nxt_ctg=beg_ctg;
				searched_ctgs[abs(nxt_ctg)]++;
				bool UniqueNxt=1,UniqueCurrent=1;
				bool Exist_Next_Unitig=0;
				int NextUnitig=0,Dist2NextUnitig=0;

				//right search
				bool Right,RightUnitig;
				if(it==1)
				{
					Right=1;
				}
				else
				{
					Right=0;
				}

				map<int, list<int> > dist_ctg,dist_ctg_rep,dist_unitig;
				
				map<int,int >::iterator ctgs_it;
				for(ctgs_it=ctgs_for_second_round.begin();ctgs_it!=ctgs_for_second_round.end();++ctgs_it)
				{
					if(it==1&&ctgs_it->second>0)
					{
						dist_ctg[ctgs_it->second].push_back(ctgs_it->first);
						dist_ctg[ctgs_it->second].push_back(ctgs_it->second);
					}
					if(it==2&&ctgs_it->second<0)
					{
						dist_ctg[-ctgs_it->second].push_back(ctgs_it->first);
						dist_ctg[-ctgs_it->second].push_back(ctgs_it->second);
					}
				}
				
				map<int,int> marked_ctgs;
				int dist_searched=0;
			
				while(1)
				{
				
					if(Right)
					{
						UniqueCurrent=UniqueNxt;
						current_ctg=abs(nxt_ctg);
					
						if(UniqueCurrent)
						{
							if(unitig_dist.count(current_ctg)==0)
							{
								if(it==1)
								{
									unitig_dist[abs(current_ctg)]=dist_searched;
								}
								else
								{
									unitig_dist[abs(current_ctg)]=-dist_searched;
								}
							}
							Exist_Next_Unitig=0;
							NextUnitig=0;
							RightUnitig=1;
						
							dist_searched+=contigs_info->contig_sz_vt[abs(current_ctg)];//
							marked_ctgs.clear();
							map<int,scaffold_contig_info>::iterator adj_it,adj_it_n;
							//check for unitig below
						
							for(adj_it=contigs_info->scaffold_adjacency_right[current_ctg].begin();adj_it!=contigs_info->scaffold_adjacency_right[current_ctg].end();++adj_it)
							{
							
								if(adj_it->second.cov>=LinkCovTh)
								{
									//record a next unitig
									if(contigs_info->c_info_vt[(abs(adj_it->first))].unique)
									{
										Exist_Next_Unitig=1;
										if(NextUnitig==0)
										{
										
											NextUnitig=adj_it->first;
											Dist2NextUnitig=adj_it->second.dist_sum/adj_it->second.cov;
										
										}
									}

									int dist0=adj_it->second.dist_sum/adj_it->second.cov;
									int dist=dist0+dist_searched;///
									if(it==1)
									{
										dist_ctg[dist].push_back(adj_it->first);
										dist_ctg[dist].push_back(dist0);
										dist_ctg_rep[dist].push_back(adj_it->first);
										dist_ctg_rep[dist].push_back(dist0);
									}
									else
									{
										dist_ctg[dist].push_back(-adj_it->first);
										dist_ctg[dist].push_back(dist0);
										dist_ctg_rep[dist].push_back(-adj_it->first);
										dist_ctg_rep[dist].push_back(dist0);

									}
								}
							}
							last_unique_ctg=current_ctg;

						}
						else
						{
							dist_searched+=contigs_info->contig_sz_vt[abs(current_ctg)];

							marked_ctgs.clear();
												
							//newly added repeat solver

							map<int,scaffold_contig_info>::iterator adj_it,adj_it_n;

							for(adj_it=contigs_info->scaffold_adjacency_right[current_ctg].begin();adj_it!=contigs_info->scaffold_adjacency_right[current_ctg].end();++adj_it)
							{
								if(adj_it->second.cov>=LinkCovTh)
								{
									int dist0=adj_it->second.dist_sum/adj_it->second.cov;
									int dist=dist0+dist_searched;///
									if(it==1)
									{
										dist_ctg_rep[dist].push_back(adj_it->first);
										dist_ctg_rep[dist].push_back(dist0);
									}
									else
									{
										dist_ctg_rep[dist].push_back(-adj_it->first);
										dist_ctg_rep[dist].push_back(dist0);
									}
								}
							}

							//newly added repeat solver
						}


						if(Exist_Next_Unitig==0)
						{
							//check if there is a unitig in the dist_ctg map
							bool UnitigFound=0;

							map<int,list<int> >::iterator dist_ctg_it,dist_ctg_it_n;
							for(dist_ctg_it=dist_ctg.begin();dist_ctg_it!=dist_ctg.end();)
							{
								dist_ctg_it_n=dist_ctg_it;
								dist_ctg_it_n++;
						
								list<int>::iterator lit,lit2,n_lit;
								for(lit=dist_ctg_it->second.begin();lit!=dist_ctg_it->second.end();)
								{
									lit2=lit;
									lit2++;
									n_lit=lit2;
									n_lit++;
									int dist0=*lit2;
						
									int DistVar2=DistVar;
													
									if(MatePair)//dist0>3000)
									{
										DistVar2=10000;//dist0;
									}
									if(contigs_info->c_info_vt[abs(*lit)].unique&&(!contigs_info->c_info_vt[abs(*lit)].used))
									if(dist_ctg_it->first>dist_searched-DistVar2)
									{
									
										if(unitig_dist.count(abs(*lit))==0)
										{
											if(it==1)
											{
												right_ctgs.push_back(*lit);
												
											}
											else
											{
												left_ctgs.push_back(*lit);
												
											}
											if(it==1)
											{
												unitig_dist[abs(*lit)]=dist_ctg_it->first;
											}
											else
											{
												unitig_dist[abs(*lit)]=-(dist_ctg_it->first);
											}
											//unitig_dist[abs(*lit)]=dist_searched;
											contigs_info->c_info_vt[abs(*lit)].used=1;

											UniqueNxt=1;
											nxt_ctg=*lit;
											t_right_ctgs.clear();
											t_left_ctgs.clear();
											contigs_info->c_info_vt[abs(*lit)].used=1;							
											NextUnitig=0;														
											dist_searched+=contigs_info->contig_sz_vt[abs(*lit)];
											UnitigFound=1;
											break;

										}
									
									}
									lit=n_lit;

						
								}
							
								if(UnitigFound)
								{break;}

								
								dist_ctg_it=dist_ctg_it_n;
								continue;



							}


							if(UnitigFound)
							{continue;}



							break;
						}
						map<int,list<int> >::iterator dist_ctg_it,dist_ctg_it_n;
						for(dist_ctg_it=dist_ctg.begin();dist_ctg_it!=dist_ctg.end();)
						{
							dist_ctg_it_n=dist_ctg_it;
							dist_ctg_it_n++;
						
							list<int>::iterator lit,lit2,n_lit;
							for(lit=dist_ctg_it->second.begin();lit!=dist_ctg_it->second.end();)
							{
								lit2=lit;
								lit2++;
								n_lit=lit2;
								n_lit++;
								int dist0=*lit2;
						
								int DistVar2=DistVar;
													
								if(MatePair)//dist0>3000)
								{
									DistVar2=20000;//dist0;
								}


								if(dist_ctg_it->first<=(dist_searched-DistVar2))
								{
								//	dist_ctg.erase(dist_ctg_it);
									dist_ctg_it->second.erase(lit2);
									dist_ctg_it->second.erase(lit);
									lit=n_lit;
									continue;
								}

							
								if((dist_ctg_it->first>dist_searched-DistVar2)&&(dist_ctg_it->first<dist_searched+DistVar2))
								{
									marked_ctgs[abs(*lit)]=dist_ctg_it->first-dist_searched;
									lit=n_lit;
									continue;
								}
								lit=n_lit;
						
							}
							if(dist_ctg_it->second.size()==0)
							{
								dist_ctg.erase(dist_ctg_it);
							}

							dist_ctg_it=dist_ctg_it_n;
							continue;



						}


						int marked_cnt=0;
						vector<int> nxt_ctg_list;
						vector<int> nxt_ctg_dist_list;
						map<int,adjacent_contig_info>::iterator adj_it;

						for(adj_it=contigs_info->contig_adjacency_right[current_ctg].begin();adj_it!=contigs_info->contig_adjacency_right[current_ctg].end();++adj_it)
						{
							if(marked_ctgs.count(abs(adj_it->first)))
							{
								marked_cnt++;
								nxt_ctg=adj_it->first;
								nxt_ctg_list.push_back(nxt_ctg);
								nxt_ctg_dist_list.push_back(marked_ctgs[abs(adj_it->first)]);
							}
						}

						bool rep_solver=0;

						if(marked_cnt==0)
						{
							map<int,int> key_unitigs;
							t_right_ctgs.clear();
							t_left_ctgs.clear();
							if(Exist_Next_Unitig)
							{
								//insert the unitigs before this 


								map<int,list<int> >::iterator dist_ctg_it,dist_ctg_it_n;
							
								for(dist_ctg_it=dist_ctg.begin();dist_ctg_it!=dist_ctg.end();)
								{
									dist_ctg_it_n=dist_ctg_it;
									dist_ctg_it_n++;
						
									list<int>::iterator lit,lit2,n_lit;
									for(lit=dist_ctg_it->second.begin();lit!=dist_ctg_it->second.end();)
									{
										lit2=lit;
										lit2++;
										n_lit=lit2;
										n_lit++;
										int dist0=*lit2;
						
										int DistVar2=DistVar;
													
										if(MatePair)//dist0>3000)
										{
											DistVar2=10000;//dist0;
										}

										if(contigs_info->c_info_vt[abs(*lit)].unique&&(!contigs_info->c_info_vt[abs(*lit)].used))
										if((dist_ctg_it->first>dist_searched-DistVar2)&&((dist_ctg_it->first+contigs_info->contig_sz_vt[abs(*lit)])<(dist_searched+Dist2NextUnitig+DistVar2)))
										{
											if(abs(*lit)==abs(NextUnitig))
											{
												break;
											}
											if(unitig_dist.count(abs(*lit))==0)
											{
												if(it==1)
												{
													t_right_ctgs.push_back(*lit);
												
												}
												else
												{
													t_left_ctgs.push_back(*lit);
												
												}
												if(it==1)
												{
													unitig_dist[abs(*lit)]=dist_ctg_it->first;
												}
												else
												{
													unitig_dist[abs(*lit)]=-(dist_ctg_it->first);
												}

												//unitig_dist[abs(*lit)]=dist_searched+Dist2NextUnitig/2;
												contigs_info->c_info_vt[abs(*lit)].used=1;
												key_unitigs[abs(*lit)]=dist_searched+Dist2NextUnitig/2;
											}
									
										}
										lit=n_lit;
						
									}
								
									dist_ctg_it=dist_ctg_it_n;
									continue;



								}












								Right=0;
								if((NextUnitig>0&&RightUnitig)||(NextUnitig<0&&RightUnitig==0))
								{
									Right=1;
								}
								UniqueNxt=1;
								nxt_ctg=NextUnitig;
								if(contigs_info->c_info_vt[abs(NextUnitig)].used==1)
								{break;}


								key_unitigs[abs(NextUnitig)]=dist_searched+Dist2NextUnitig;
							
								map<int,vector<int> > node_cov;
								if(it==1)
								{


									if(RightUnitig)
									{
										t_right_ctgs.push_back(NextUnitig);
									
										list<int> ctg_stack;
										ctg_stack.push_back(current_ctg);
									
										for(int d=3;d<max_dep;d+=2)
										{

											BFSearchPathFinder(contigs_info, ctg_stack,dist_ctg,key_unitigs,dist_searched,node_cov,d,1,MatePair);//,sc_cnt
											if(node_cov.size()==1)
											{
											
												break;
										
											}	
										
									

										}
									
							
									}
									else
									{
										t_right_ctgs.push_back(-NextUnitig);

										list<int> ctg_stack;
										ctg_stack.push_back(-current_ctg);
									
								
										for(int d=3;d<max_dep;d+=2)
										{

											BFSearchPathFinder(contigs_info, ctg_stack,dist_ctg,key_unitigs,dist_searched,node_cov,d,1,MatePair);//,sc_cnt
											if(node_cov.size()==1)
											{
											
												break;
										
											}
										}
								

									}
								


									if(node_cov.size()==1)
									{
										for(int ii=0;ii<t_right_ctgs.size();++ii)
										{
											if(t_right_ctgs[ii]==((*node_cov.begin()).first))
											{
												vector<int> tmp_vt=(*node_cov.begin()).second;
												for(int jj=1;jj<tmp_vt.size()-1;++jj)
												{
													if(UniqueCurrent&&jj==1)
													{continue;}
													right_ctgs.push_back(tmp_vt[jj]);
												}
												for(int jj=ii;jj<t_right_ctgs.size();++jj)
												{
													right_ctgs.push_back(t_right_ctgs[jj]);
												}
												break;
											}
											if(t_right_ctgs[ii]==-((*node_cov.begin()).first))
											{
												vector<int> tmp_vt=(*node_cov.begin()).second;
												for(int jj=1;jj<tmp_vt.size()-1;++jj)
												{
													if(UniqueCurrent&&jj==1)
													{continue;}
													right_ctgs.push_back(-tmp_vt[jj]);
												}
												for(int jj=ii;jj<t_right_ctgs.size();++jj)
												{
													right_ctgs.push_back(t_right_ctgs[jj]);
												}
												break;
											}
										}
										//break;
										
									}	
									else
									{
										for(int jj=0;jj<t_right_ctgs.size();++jj)
										{
											right_ctgs.push_back(t_right_ctgs[jj]);
										}
									}
										


									t_right_ctgs.clear();
								}
								else
								{
									if(RightUnitig)
									{
										t_left_ctgs.push_back(-NextUnitig);

											list<int> ctg_stack;
										ctg_stack.push_back(current_ctg);
									
										for(int d=3;d<max_dep;d+=2)
										{

											BFSearchPathFinder(contigs_info, ctg_stack,dist_ctg,key_unitigs,dist_searched,node_cov,d,1,MatePair);//,sc_cnt
											if(node_cov.size()==1)
											{
											
												break;
										
											}	
										
										}
								
									}
									else
									{
										t_left_ctgs.push_back(NextUnitig);


											list<int> ctg_stack;
										ctg_stack.push_back(-current_ctg);
									
										for(int d=3;d<max_dep;d+=2)
										{

											BFSearchPathFinder(contigs_info, ctg_stack,dist_ctg,key_unitigs,dist_searched,node_cov,d,1,MatePair);//,sc_cnt
											if(node_cov.size()==1)
											{
											
												break;
										
											}	
										
										}
								
									}





									if(node_cov.size()==1)
									{
										for(int ii=0;ii<t_left_ctgs.size();++ii)
										{
											if(t_left_ctgs[ii]==((*node_cov.begin()).first))
											{
												vector<int> tmp_vt=(*node_cov.begin()).second;
												for(int jj=1;jj<tmp_vt.size()-1;++jj)
												{
													if(UniqueCurrent&&jj==1)
													{continue;}
													left_ctgs.push_back(tmp_vt[jj]);

												}
												for(int jj=ii;jj<t_left_ctgs.size();++jj)
												{
													left_ctgs.push_back(t_left_ctgs[jj]);
												}
												break;
											}
											if(t_left_ctgs[ii]==-((*node_cov.begin()).first))
											{
												vector<int> tmp_vt=(*node_cov.begin()).second;
												for(int jj=1;jj<tmp_vt.size()-1;++jj)
												{
													if(UniqueCurrent&&jj==1)
													{continue;}
													left_ctgs.push_back(-tmp_vt[jj]);
												}
												for(int jj=ii;jj<t_left_ctgs.size();++jj)
												{
													left_ctgs.push_back(t_left_ctgs[jj]);
												}
												break;
											}
										}
									//	break;
										
									}
									else
									{
										for(int jj=0;jj<t_left_ctgs.size();++jj)
										{
											left_ctgs.push_back(t_left_ctgs[jj]);
										}
									}
								

								
									t_left_ctgs.clear();
								}

								contigs_info->c_info_vt[abs(NextUnitig)].used=1;
								NextUnitig=0;
								//key_unitigs[abs(NextUnitig)]=dist_searched+Dist2NextUnitig;
								dist_searched+=Dist2NextUnitig;

								continue;
							}

							break;
						
						
						}


						//newly added repeat solver


						if(marked_cnt>1)
						{
							map<int,int> key_unitigs;
							t_right_ctgs.clear();
							t_left_ctgs.clear();
							if(Exist_Next_Unitig)
							{



								map<int,list<int> >::iterator dist_ctg_it,dist_ctg_it_n;
							
								for(dist_ctg_it=dist_ctg.begin();dist_ctg_it!=dist_ctg.end();)
								{
									dist_ctg_it_n=dist_ctg_it;
									dist_ctg_it_n++;
						
									list<int>::iterator lit,lit2,n_lit;
									for(lit=dist_ctg_it->second.begin();lit!=dist_ctg_it->second.end();)
			 						{
										lit2=lit;
										lit2++;
										n_lit=lit2;
										n_lit++;
										int dist0=*lit2;
						
										int DistVar2=DistVar;
													
										if(MatePair)//dist0>3000)
										{
											DistVar2=10000;//dist0;
										}
										if(contigs_info->c_info_vt[abs(*lit)].unique&&(!contigs_info->c_info_vt[abs(*lit)].used))
										if((dist_ctg_it->first>dist_searched-DistVar2)&&((dist_ctg_it->first+contigs_info->contig_sz_vt[abs(*lit)])<(dist_searched+Dist2NextUnitig+DistVar2)))
										{
											if(abs(*lit)==abs(NextUnitig))
											{
												break;
											}
											if(unitig_dist.count(abs(*lit))==0)
											{
												if(it==1)
												{
													t_right_ctgs.push_back(*lit);
												
												}
												else
												{
													t_left_ctgs.push_back(*lit);
												
												}

												if(it==1)
												{
													unitig_dist[abs(*lit)]=dist_ctg_it->first;
												}
												else
												{
													unitig_dist[abs(*lit)]=-(dist_ctg_it->first);
												}

												//unitig_dist[abs(*lit)]=dist_searched+Dist2NextUnitig/2;
												contigs_info->c_info_vt[abs(*lit)].used=1;
												key_unitigs[abs(*lit)]=dist_searched+Dist2NextUnitig/2;
											}
										
										}
										lit=n_lit;
						
									}
								
									dist_ctg_it=dist_ctg_it_n;
									continue;



								}










								Right=0;
								if((NextUnitig>0&&RightUnitig)||(NextUnitig<0&&RightUnitig==0))
								{
									Right=1;
								}
								UniqueNxt=1;
								nxt_ctg=NextUnitig;
								if(contigs_info->c_info_vt[abs(NextUnitig)].used==1)
								{break;}

								key_unitigs[abs(NextUnitig)]=dist_searched+Dist2NextUnitig;
							
								map<int,vector<int> > node_cov;
								
								if(it==1)
								{


									if(RightUnitig)
									{
										t_right_ctgs.push_back(NextUnitig);
									
										list<int> ctg_stack;
										ctg_stack.push_back(current_ctg);
								
										for(int d=3;d<max_dep;d+=2)
										{

											BFSearchPathFinder(contigs_info, ctg_stack,dist_ctg,key_unitigs,dist_searched,node_cov,d,1,MatePair);//,sc_cnt
											if(node_cov.size()==1)
											{
											
												break;
										
											}	
										
									
										}
									
							
									}
									else
									{
										t_right_ctgs.push_back(-NextUnitig);

										list<int> ctg_stack;
										ctg_stack.push_back(-current_ctg);
								
										for(int d=3;d<max_dep;d+=2)
										{

											BFSearchPathFinder(contigs_info, ctg_stack,dist_ctg,key_unitigs,dist_searched,node_cov,d,1,MatePair);//,sc_cnt
											if(node_cov.size()==1)
											{
											
												break;
										
											}	
										
									
										}
								

									}
								









									if(node_cov.size()==1)
									{
										for(int ii=0;ii<t_right_ctgs.size();++ii)
										{
											if(t_right_ctgs[ii]==((*node_cov.begin()).first))
											{
												vector<int> tmp_vt=(*node_cov.begin()).second;
												for(int jj=1;jj<tmp_vt.size()-1;++jj)
												{
													if(UniqueCurrent&&jj==1)
													{continue;}
													right_ctgs.push_back(tmp_vt[jj]);
												}
												for(int jj=ii;jj<t_right_ctgs.size();++jj)
												{
													right_ctgs.push_back(t_right_ctgs[jj]);
												}
												break;
											}
											if(t_right_ctgs[ii]==-((*node_cov.begin()).first))
											{
												vector<int> tmp_vt=(*node_cov.begin()).second;
												for(int jj=1;jj<tmp_vt.size()-1;++jj)
												{
													if(UniqueCurrent&&jj==1)
													{continue;}
													right_ctgs.push_back(-tmp_vt[jj]);
												}
												for(int jj=ii;jj<t_right_ctgs.size();++jj)
												{
													right_ctgs.push_back(t_right_ctgs[jj]);
												}
												break;
											}
										}
										//break;
										
									}	
									else
									{
										for(int jj=0;jj<t_right_ctgs.size();++jj)
										{
											right_ctgs.push_back(t_right_ctgs[jj]);
										}
									}		

									t_right_ctgs.clear();
								}
								else
								{
									if(RightUnitig)
									{
										t_left_ctgs.push_back(-NextUnitig);

											list<int> ctg_stack;
										ctg_stack.push_back(current_ctg);
									
										for(int d=3;d<max_dep;d+=2)
										{

											BFSearchPathFinder(contigs_info, ctg_stack,dist_ctg,key_unitigs,dist_searched,node_cov,d,1,MatePair);//,sc_cnt
											if(node_cov.size()==1)
											{
											
												break;
										
											}	
										
										}
								
									}
									else
									{
										t_left_ctgs.push_back(NextUnitig);


											list<int> ctg_stack;
										ctg_stack.push_back(-current_ctg);
									
										for(int d=3;d<max_dep;d+=2)
										{

											BFSearchPathFinder(contigs_info, ctg_stack,dist_ctg,key_unitigs,dist_searched,node_cov,d,1,MatePair);//,sc_cnt
											if(node_cov.size()==1)
											{
											
												break;
										
											}	
										
										}
								
									}
								








									if(node_cov.size()==1)
									{
										for(int ii=0;ii<t_left_ctgs.size();++ii)
										{
											if(t_left_ctgs[ii]==((*node_cov.begin()).first))
											{
												vector<int> tmp_vt=(*node_cov.begin()).second;
												for(int jj=1;jj<tmp_vt.size()-1;++jj)
												{
													if(UniqueCurrent&&jj==1)
													{continue;}
													left_ctgs.push_back(tmp_vt[jj]);
												}
												for(int jj=ii;jj<t_left_ctgs.size();++jj)
												{
													left_ctgs.push_back(t_left_ctgs[jj]);
												}
												break;
											}
											if(t_left_ctgs[ii]==-((*node_cov.begin()).first))
											{
												vector<int> tmp_vt=(*node_cov.begin()).second;
												for(int jj=1;jj<tmp_vt.size()-1;++jj)
												{
													if(UniqueCurrent&&jj==1)
													{continue;}
													left_ctgs.push_back(-tmp_vt[jj]);
												}
												for(int jj=ii;jj<t_left_ctgs.size();++jj)
												{
													left_ctgs.push_back(t_left_ctgs[jj]);
												}
												break;
											}
										}
									//	break;
										
									}
									else
									{
										for(int jj=0;jj<t_left_ctgs.size();++jj)
										{
											left_ctgs.push_back(t_left_ctgs[jj]);
										}
									}
								

								


									t_left_ctgs.clear();
								}

							
								contigs_info->c_info_vt[abs(NextUnitig)].used=1;
								NextUnitig=0;
								//key_unitigs[abs(NextUnitig)]=dist_searched+Dist2NextUnitig;
								dist_searched+=Dist2NextUnitig;
								continue;
							}

					
							break;
						
						}
						if(marked_cnt==1)
						{
							//dist_searched-=marked_ctgs[abs(nxt_ctg)]/2;//smoothing
							if(searched_ctgs.count(abs(nxt_ctg))==0)
							{
								searched_ctgs[abs(nxt_ctg)]++;
							}
							else
							{
								searched_ctgs[abs(nxt_ctg)]++;
								if(searched_ctgs[abs(nxt_ctg)]>4)
								{
									map<int,int> key_unitigs;
									t_right_ctgs.clear();
									t_left_ctgs.clear();
									if(Exist_Next_Unitig)
									{

										map<int,list<int> >::iterator dist_ctg_it,dist_ctg_it_n;
									
										for(dist_ctg_it=dist_ctg.begin();dist_ctg_it!=dist_ctg.end();)
										{
											dist_ctg_it_n=dist_ctg_it;
											dist_ctg_it_n++;
						
											list<int>::iterator lit,lit2,n_lit;
											for(lit=dist_ctg_it->second.begin();lit!=dist_ctg_it->second.end();)
											{
												lit2=lit;
												lit2++;
												n_lit=lit2;
												n_lit++;
												int dist0=*lit2;
						
												int DistVar2=DistVar;
													
												if(MatePair)//dist0>3000)
												{
													DistVar2=10000;//dist0;
												}
												if(contigs_info->c_info_vt[abs(*lit)].unique&&(!contigs_info->c_info_vt[abs(*lit)].used))
												if((dist_ctg_it->first>dist_searched-DistVar2)&&((dist_ctg_it->first+contigs_info->contig_sz_vt[abs(*lit)])<(dist_searched+Dist2NextUnitig+DistVar2)))
												{
													if(abs(*lit)==abs(NextUnitig))
													{
														break;
													}
													if(unitig_dist.count(abs(*lit))==0)
													{
														if(it==1)
														{
															t_right_ctgs.push_back(*lit);
														
														}
														else
														{
															t_left_ctgs.push_back(*lit);
														
														}

														if(it==1)
														{
															unitig_dist[abs(*lit)]=dist_ctg_it->first;
														}
														else
														{
															unitig_dist[abs(*lit)]=-(dist_ctg_it->first);
														}
														contigs_info->c_info_vt[abs(*lit)].used=1;
														//unitig_dist[abs(*lit)]=dist_searched+Dist2NextUnitig/2;
														key_unitigs[abs(*lit)]=dist_searched+Dist2NextUnitig/2;
													}
												
												}
												lit=n_lit;
						
											}
								
											dist_ctg_it=dist_ctg_it_n;
											continue;



										}









										Right=0;
										if(NextUnitig>0)
										{
											Right=1;
										}
										UniqueNxt=1;
										nxt_ctg=NextUnitig;
										if(contigs_info->c_info_vt[abs(NextUnitig)].used==1)
										{break;}

										key_unitigs[abs(NextUnitig)]=dist_searched+Dist2NextUnitig;
										map<int,vector<int> > node_cov;
								
										if(it==1)
										{


											if(RightUnitig)
											{
												t_right_ctgs.push_back(NextUnitig);
									
												list<int> ctg_stack;
												ctg_stack.push_back(current_ctg);
											
												for(int d=3;d<max_dep;d+=2)
												{

													BFSearchPathFinder(contigs_info, ctg_stack,dist_ctg,key_unitigs,dist_searched,node_cov,d,1,MatePair);//,sc_cnt
													if(node_cov.size()==1)
													{
													
														break;
										
													}	
										
									
												}
									
							
											}
											else
											{
												t_right_ctgs.push_back(-NextUnitig);

												list<int> ctg_stack;
												ctg_stack.push_back(-current_ctg);
											
												for(int d=3;d<max_dep;d+=2)
												{

													BFSearchPathFinder(contigs_info, ctg_stack,dist_ctg,key_unitigs,dist_searched,node_cov,d,1,MatePair);//,sc_cnt
													if(node_cov.size()==1)
													{
														break;
										
													}	
										
									
												}
								

											}
								










											if(node_cov.size()==1)
											{
												for(int ii=0;ii<t_right_ctgs.size();++ii)
												{
													if(t_right_ctgs[ii]==((*node_cov.begin()).first))
													{
														vector<int> tmp_vt=(*node_cov.begin()).second;
														for(int jj=1;jj<tmp_vt.size()-1;++jj)
														{
															if(UniqueCurrent&&jj==1)
															{continue;}
															right_ctgs.push_back(tmp_vt[jj]);
														}
														for(int jj=ii;jj<t_right_ctgs.size();++jj)
														{
															right_ctgs.push_back(t_right_ctgs[jj]);
														}
														break;
													}
													if(t_right_ctgs[ii]==-((*node_cov.begin()).first))
													{
														vector<int> tmp_vt=(*node_cov.begin()).second;
														for(int jj=1;jj<tmp_vt.size()-1;++jj)
														{
															if(UniqueCurrent&&jj==1)
															{continue;}
															right_ctgs.push_back(-tmp_vt[jj]);
														}
														for(int jj=ii;jj<t_right_ctgs.size();++jj)
														{
															right_ctgs.push_back(t_right_ctgs[jj]);
														}
														break;
													}
												}
											//	break;
										
											}
											else
											{
												for(int jj=0;jj<t_right_ctgs.size();++jj)
												{
													right_ctgs.push_back(t_right_ctgs[jj]);
												}
											}
										

											t_right_ctgs.clear();
										}
										else
										{
											if(RightUnitig)
											{
												t_left_ctgs.push_back(-NextUnitig);

												list<int> ctg_stack;
												ctg_stack.push_back(current_ctg);
											
												for(int d=3;d<max_dep;d+=2)
												{

													BFSearchPathFinder(contigs_info, ctg_stack,dist_ctg,key_unitigs,dist_searched,node_cov,d,1,MatePair);//,sc_cnt
													if(node_cov.size()==1)
													{
													
														break;
										
													}	
										
												}
								
											}
											else
											{
												t_left_ctgs.push_back(NextUnitig);


													list<int> ctg_stack;
												ctg_stack.push_back(-current_ctg);
								
												for(int d=3;d<max_dep;d+=2)
												{

													BFSearchPathFinder(contigs_info, ctg_stack,dist_ctg,key_unitigs,dist_searched,node_cov,d,1,MatePair);//,sc_cnt
													if(node_cov.size()==1)
													{
													
														break;
										
													}	
										
												}
								
											}



											if(node_cov.size()==1)
											{
												for(int ii=0;ii<t_left_ctgs.size();++ii)
												{
													if(t_left_ctgs[ii]==((*node_cov.begin()).first))
													{
														vector<int> tmp_vt=(*node_cov.begin()).second;
														for(int jj=1;jj<tmp_vt.size()-1;++jj)
														{
															if(UniqueCurrent&&jj==1)
															{continue;}
															left_ctgs.push_back(tmp_vt[jj]);
														}
														for(int jj=ii;jj<t_left_ctgs.size();++jj)
														{
															left_ctgs.push_back(t_left_ctgs[jj]);
														}
														break;
													}
													if(t_left_ctgs[ii]==-((*node_cov.begin()).first))
													{
														vector<int> tmp_vt=(*node_cov.begin()).second;
														for(int jj=1;jj<tmp_vt.size()-1;++jj)
														{
															if(UniqueCurrent&&jj==1)
															{continue;}
															left_ctgs.push_back(-tmp_vt[jj]);
														}
														for(int jj=ii;jj<t_left_ctgs.size();++jj)
														{
															left_ctgs.push_back(t_left_ctgs[jj]);
														}
														break;
													}
												}
											//	break;
										
											}
											else
											{
												for(int jj=0;jj<t_left_ctgs.size();++jj)
												{
													left_ctgs.push_back(t_left_ctgs[jj]);
												}
											}
											t_left_ctgs.clear();
										}

								
										contigs_info->c_info_vt[abs(NextUnitig)].used=1;
										NextUnitig=0;
										//key_unitigs[abs(NextUnitig)]=dist_searched+Dist2NextUnitig;
										dist_searched+=Dist2NextUnitig;
										continue;
									}

									break;
								}
							}



							if(nxt_ctg>0)
							{
								UniqueNxt=contigs_info->c_info_vt[nxt_ctg].unique;

								if(UniqueNxt&&contigs_info->c_info_vt[nxt_ctg].used)
								{
									break;
									contigs_info->c_info_vt[nxt_ctg].unique=0;
									//break;
								}

								Right=1;
								if(it==1)
								{
									t_right_ctgs.push_back(nxt_ctg);
								}
								else
								{
									t_left_ctgs.push_back(-nxt_ctg);
								}
								//contigs_info->c_info_vt[abs(nxt_ctg)].used=1;
								//right search next
							}
							else
							{
								UniqueNxt=contigs_info->c_info_vt[-nxt_ctg].unique;

								if(UniqueNxt&&contigs_info->c_info_vt[-nxt_ctg].used)
								{
								
									UniqueNxt=contigs_info->c_info_vt[-nxt_ctg].unique=0;
								//	break;
								}

								Right=0;
								if(it==1)
								{
									t_right_ctgs.push_back(nxt_ctg);
								}
								else
								{
									t_left_ctgs.push_back(-nxt_ctg);
								}
								//contigs_info->c_info_vt[abs(nxt_ctg)].used=1;
								//left search next

							}


							if(UniqueNxt)
							{
								if(it==1)
								{
									for (int ii=0;ii<t_right_ctgs.size();++ii)
									{
										right_ctgs.push_back(t_right_ctgs[ii]);
										contigs_info->c_info_vt[abs(t_right_ctgs[ii])].used=1;
									}
									t_right_ctgs.clear();
								}
								else
								{
									for (int ii=0;ii<t_left_ctgs.size();++ii)
									{
										left_ctgs.push_back(t_left_ctgs[ii]);
										contigs_info->c_info_vt[abs(t_left_ctgs[ii])].used=1;
									}
									t_left_ctgs.clear();
							
								}
							}


						}
						else
						{
							map<int,int> key_unitigs;

							t_right_ctgs.clear();
							t_left_ctgs.clear();
							if(Exist_Next_Unitig)
							{

								map<int,list<int> >::iterator dist_ctg_it,dist_ctg_it_n;
							
								for(dist_ctg_it=dist_ctg.begin();dist_ctg_it!=dist_ctg.end();)
								{
									dist_ctg_it_n=dist_ctg_it;
									dist_ctg_it_n++;
						
									list<int>::iterator lit,lit2,n_lit;
									for(lit=dist_ctg_it->second.begin();lit!=dist_ctg_it->second.end();)
									{
										lit2=lit;
										lit2++;
										n_lit=lit2;
										n_lit++;
										int dist0=*lit2;
						
										int DistVar2=DistVar;
													
										if(MatePair)//dist0>3000)
										{
											DistVar2=10000;//dist0;
										}
										if(contigs_info->c_info_vt[abs(*lit)].unique&&(!contigs_info->c_info_vt[abs(*lit)].used))
										if((dist_ctg_it->first>dist_searched-DistVar2)&&((dist_ctg_it->first+contigs_info->contig_sz_vt[abs(*lit)])<(dist_searched+Dist2NextUnitig+DistVar2)))
										{
											if(abs(*lit)==abs(NextUnitig))
											{
												break;
											}
											if(unitig_dist.count(abs(*lit))==0)
											{
												if(it==1)
												{
													t_right_ctgs.push_back(*lit);
												}
												else
												{
													t_left_ctgs.push_back(*lit);
											
												}


												if(it==1)
												{
													unitig_dist[abs(*lit)]=dist_ctg_it->first;
												}
												else
												{
													unitig_dist[abs(*lit)]=-(dist_ctg_it->first);
												}

												contigs_info->c_info_vt[abs(*lit)].used=1;
												//unitig_dist[abs(*lit)]=dist_searched+Dist2NextUnitig/2;
												key_unitigs[abs(*lit)]=dist_searched+Dist2NextUnitig/2;
											}

										
										}
										lit=n_lit;
						
									}
								
									dist_ctg_it=dist_ctg_it_n;
									continue;



								}








								Right=0;
								if((NextUnitig>0&&RightUnitig)||(NextUnitig<0&&RightUnitig==0))
								{
									Right=1;
								}
								UniqueNxt=1;
								nxt_ctg=NextUnitig;
								if(contigs_info->c_info_vt[abs(NextUnitig)].used==1)
								{break;}

								key_unitigs[abs(NextUnitig)]=dist_searched+Dist2NextUnitig;
								map<int,vector<int> > node_cov;
								


								if(it==1)
								{


									if(RightUnitig)
									{
										t_right_ctgs.push_back(NextUnitig);
									
										list<int> ctg_stack;
										ctg_stack.push_back(current_ctg);
									
										for(int d=3;d<max_dep;d+=2)
										{

											BFSearchPathFinder(contigs_info, ctg_stack,dist_ctg,key_unitigs,dist_searched,node_cov,d,1,MatePair);//,sc_cnt
											if(node_cov.size()==1)
											{
											
												break;
										
											}	
										
										}
									
							
									}
									else
									{
										t_right_ctgs.push_back(-NextUnitig);

										list<int> ctg_stack;
										ctg_stack.push_back(-current_ctg);
									
										for(int d=3;d<max_dep;d+=2)
										{

											BFSearchPathFinder(contigs_info, ctg_stack,dist_ctg,key_unitigs,dist_searched,node_cov,d,1,MatePair);//,sc_cnt
											if(node_cov.size()==1)
											{
											
												break;
										
											}	
										
									

										}
								

									}



								



									if(node_cov.size()==1)
									{
										for(int ii=0;ii<t_right_ctgs.size();++ii)
										{
											if(t_right_ctgs[ii]==((*node_cov.begin()).first))
											{
												vector<int> tmp_vt=(*node_cov.begin()).second;
												for(int jj=1;jj<tmp_vt.size()-1;++jj)
												{
													if(UniqueCurrent&&jj==1)
													{continue;}
													right_ctgs.push_back(tmp_vt[jj]);
												}
												for(int jj=ii;jj<t_right_ctgs.size();++jj)
												{
													right_ctgs.push_back(t_right_ctgs[jj]);
												}
												break;
											}
											if(t_right_ctgs[ii]==-((*node_cov.begin()).first))
											{
												vector<int> tmp_vt=(*node_cov.begin()).second;
												for(int jj=1;jj<tmp_vt.size()-1;++jj)
												{
													if(UniqueCurrent&&jj==1)
													{continue;}
													
													right_ctgs.push_back(-tmp_vt[jj]);
												}
												for(int jj=ii;jj<t_right_ctgs.size();++jj)
												{
													right_ctgs.push_back(t_right_ctgs[jj]);
												}
												break;
											}
										}
									//	break;
										
									}	
									else
									{
										for(int jj=0;jj<t_right_ctgs.size();++jj)
										{
											right_ctgs.push_back(t_right_ctgs[jj]);
										}
									}		

									t_right_ctgs.clear();
								}
								else
								{
									if(RightUnitig)
									{
										t_left_ctgs.push_back(-NextUnitig);

											list<int> ctg_stack;
										ctg_stack.push_back(current_ctg);
								
										for(int d=3;d<max_dep;d+=2)
										{

											BFSearchPathFinder(contigs_info, ctg_stack,dist_ctg,key_unitigs,dist_searched,node_cov,d,1,MatePair);//,sc_cnt
											if(node_cov.size()==1)
											{
										
												break;


											}
										}
								
									}
									else
									{
										t_left_ctgs.push_back(NextUnitig);


											list<int> ctg_stack;
										ctg_stack.push_back(-current_ctg);
									
										for(int d=3;d<max_dep;d+=2)
										{

											BFSearchPathFinder(contigs_info, ctg_stack,dist_ctg,key_unitigs,dist_searched,node_cov,d,1,MatePair);//,sc_cnt
											if(node_cov.size()==1)
											{
											
												break;

											}
										}
								
									}



									if(node_cov.size()==1)
									{
										for(int ii=0;ii<t_left_ctgs.size();++ii)
										{
											if(t_left_ctgs[ii]==((*node_cov.begin()).first))
											{
												vector<int> tmp_vt=(*node_cov.begin()).second;
												for(int jj=1;jj<tmp_vt.size()-1;++jj)
												{
													if(UniqueCurrent&&jj==1)
													{continue;}
													
													left_ctgs.push_back(tmp_vt[jj]);
												}
												for(int jj=ii;jj<t_left_ctgs.size();++jj)
												{
													left_ctgs.push_back(t_left_ctgs[jj]);
												}
												break;
											}
											if(t_left_ctgs[ii]==-((*node_cov.begin()).first))
											{
												vector<int> tmp_vt=(*node_cov.begin()).second;
												for(int jj=1;jj<tmp_vt.size()-1;++jj)
												{
													if(UniqueCurrent&&jj==1)
													{continue;}
													
													left_ctgs.push_back(-tmp_vt[jj]);
												}
												for(int jj=ii;jj<t_left_ctgs.size();++jj)
												{
													left_ctgs.push_back(t_left_ctgs[jj]);
												}
												break;
											}
										}
									//	break;
										
									}
									else
									{
										for(int jj=0;jj<t_left_ctgs.size();++jj)
										{
											left_ctgs.push_back(t_left_ctgs[jj]);
										}
									}
								
								
									t_left_ctgs.clear();
								}

							
								contigs_info->c_info_vt[abs(NextUnitig)].used=1;
								NextUnitig=0;
								dist_searched+=Dist2NextUnitig;
								continue;
							}

							break;
						


						}
					







					}
					else
					{
						//left

						UniqueCurrent=UniqueNxt;
						current_ctg=abs(nxt_ctg);
						if(UniqueCurrent)
						{
							if(unitig_dist.count(current_ctg)==0)
							{
								//unitig_dist[current_ctg]=dist_searched;
								if(it==1)
								{
									unitig_dist[abs(current_ctg)]=dist_searched;
								}
								else
								{
									unitig_dist[abs(current_ctg)]=-dist_searched;
								}
							}
							Exist_Next_Unitig=0;
							NextUnitig=0;
							RightUnitig=0;
							//dist_searched=0;
						//	dist_ctg.clear();
							marked_ctgs.clear();
							dist_searched+=contigs_info->contig_sz_vt[abs(current_ctg)];///
						
							map<int,scaffold_contig_info>::iterator adj_it;

							for(adj_it=contigs_info->scaffold_adjacency_left[current_ctg].begin();adj_it!=contigs_info->scaffold_adjacency_left[current_ctg].end();++adj_it)
							{

								if(adj_it->second.cov>=LinkCovTh)
								{
									if(contigs_info->c_info_vt[(abs(adj_it->first))].unique)
									{
										Exist_Next_Unitig=1;
										if(NextUnitig==0)
										{
										
											NextUnitig=adj_it->first;
											Dist2NextUnitig=adj_it->second.dist_sum/adj_it->second.cov;
										
										
										}
									}
									int dist0=adj_it->second.dist_sum/adj_it->second.cov;
									int dist=dist0+dist_searched;///
									if(it==1)
									{
										dist_ctg[dist].push_back(-(adj_it->first));//////////////////////////////
										dist_ctg[dist].push_back(dist0);
										dist_ctg_rep[dist].push_back(-(adj_it->first));
										dist_ctg_rep[dist].push_back(dist0);
									}
									else
									{
										dist_ctg[dist].push_back((adj_it->first));	
										dist_ctg[dist].push_back(dist0);
										dist_ctg_rep[dist].push_back(adj_it->first);
										dist_ctg_rep[dist].push_back(dist0);
									}
								}
							}


						}
						else
						{
							dist_searched+=contigs_info->contig_sz_vt[abs(current_ctg)];

							marked_ctgs.clear();


							//newly added repeat solver
							map<int,scaffold_contig_info>::iterator adj_it;

							for(adj_it=contigs_info->scaffold_adjacency_left[current_ctg].begin();adj_it!=contigs_info->scaffold_adjacency_left[current_ctg].end();++adj_it)
							{
								if(adj_it->second.cov>=LinkCovTh)
								{
									int dist0=adj_it->second.dist_sum/adj_it->second.cov;
									int dist=dist0+dist_searched;///
									if(it==1)
									{
										dist_ctg_rep[dist].push_back(-(adj_it->first));//////////////////////////////
										dist_ctg_rep[dist].push_back(dist0);
									}
									else
									{
										dist_ctg_rep[dist].push_back((adj_it->first));								
										dist_ctg_rep[dist].push_back(dist0);
									}
								}
							}
						
							//newly added repeat solver

						}

						if(Exist_Next_Unitig==0)
						{
							//check ....
							bool UnitigFound=0;

							map<int,list<int> >::iterator dist_ctg_it,dist_ctg_it_n;
							for(dist_ctg_it=dist_ctg.begin();dist_ctg_it!=dist_ctg.end();)
							{
								dist_ctg_it_n=dist_ctg_it;
								dist_ctg_it_n++;
						
								list<int>::iterator lit,lit2,n_lit;
								for(lit=dist_ctg_it->second.begin();lit!=dist_ctg_it->second.end();)
								{
									lit2=lit;
									lit2++;
									n_lit=lit2;
									n_lit++;
									int dist0=*lit2;
						
									int DistVar2=DistVar;
													
									if(MatePair)//dist0>3000)
									{
										DistVar2=10000;//dist0;
									}
									if(contigs_info->c_info_vt[abs(*lit)].unique&&(!contigs_info->c_info_vt[abs(*lit)].used))										
									if(dist_ctg_it->first>dist_searched-DistVar2)
									{
									
										if(unitig_dist.count(abs(*lit))==0)
										{
											if(it==1)
											{
												right_ctgs.push_back(*lit);
												
											}
											else
											{
												left_ctgs.push_back(*lit);
												
											}
											if(it==1)
											{
												unitig_dist[abs(*lit)]=dist_ctg_it->first;
											}
											else
											{
												unitig_dist[abs(*lit)]=-(dist_ctg_it->first);
											}
											//unitig_dist[abs(*lit)]=dist_searched;
											contigs_info->c_info_vt[abs(*lit)].used=1;

											UniqueNxt=1;
											nxt_ctg=*lit;
											t_right_ctgs.clear();
											t_left_ctgs.clear();
											contigs_info->c_info_vt[abs(*lit)].used=1;							
											NextUnitig=0;														
											dist_searched+=contigs_info->contig_sz_vt[abs(*lit)];
											UnitigFound=1;
											break;

										}
									
									}
									lit=n_lit;
						
								}
								
								if(UnitigFound)
								{break;}

								dist_ctg_it=dist_ctg_it_n;
								continue;



							}


							if(UnitigFound)
							{continue;}

							break;
						}

						map<int,list<int> >::iterator dist_ctg_it,dist_ctg_it_n;
						for(dist_ctg_it=dist_ctg.begin();dist_ctg_it!=dist_ctg.end();)
						{
							dist_ctg_it_n=dist_ctg_it;
							dist_ctg_it_n++;
						


						
							list<int>::iterator lit,lit2,n_lit;
							for(lit=dist_ctg_it->second.begin();lit!=dist_ctg_it->second.end();)
							{
								lit2=lit;
								lit2++;
								n_lit=lit2;
								n_lit++;
								int dist0=*lit2;
						
								int DistVar2=DistVar;
													
								if(MatePair)//dist0>3000)
								{
									DistVar2=20000;//dist0;//dist0/DistVarFactor;;
								}


								if(dist_ctg_it->first<=(dist_searched-DistVar2))
								{
								//	dist_ctg.erase(dist_ctg_it);
									dist_ctg_it->second.erase(lit2);
									dist_ctg_it->second.erase(lit);
									lit=n_lit;
									continue;
								}

							
								if((dist_ctg_it->first>dist_searched-DistVar2)&&(dist_ctg_it->first<dist_searched+DistVar2))
								{
									marked_ctgs[abs(*lit)]=dist_ctg_it->first-dist_searched;
									lit=n_lit;
									continue;
								}
								lit=n_lit;
						
							}
							if(dist_ctg_it->second.size()==0)
							{
								dist_ctg.erase(dist_ctg_it);
							}


						
							dist_ctg_it=dist_ctg_it_n;

					
						}




						int marked_cnt=0;
						vector<int> nxt_ctg_list;
						vector<int> nxt_ctg_dist_list;
						map<int,adjacent_contig_info>::iterator adj_it;

						for(adj_it=contigs_info->contig_adjacency_left[current_ctg].begin();adj_it!=contigs_info->contig_adjacency_left[current_ctg].end();++adj_it)
						{
							if(marked_ctgs.count(abs(adj_it->first)))
							{
								marked_cnt++;
								nxt_ctg=adj_it->first;
								nxt_ctg_list.push_back(nxt_ctg);
								nxt_ctg_dist_list.push_back(marked_ctgs[abs(adj_it->first)]);
						
							}
						}

						bool rep_solver=0;
						if(marked_cnt==0)
						{
							map<int,int> key_unitigs;
							t_right_ctgs.clear();
							t_left_ctgs.clear();
							if(Exist_Next_Unitig)
							{


								map<int,list<int> >::iterator dist_ctg_it,dist_ctg_it_n;
							
								for(dist_ctg_it=dist_ctg.begin();dist_ctg_it!=dist_ctg.end();)
								{
									dist_ctg_it_n=dist_ctg_it;
									dist_ctg_it_n++;
						
									list<int>::iterator lit,lit2,n_lit;
									for(lit=dist_ctg_it->second.begin();lit!=dist_ctg_it->second.end();)
									{
										lit2=lit;
										lit2++;
										n_lit=lit2;
										n_lit++;
										int dist0=*lit2;
						
										int DistVar2=DistVar;
													
										if(MatePair)//dist0>3000)
										{
											DistVar2=10000;//dist0;
										}
										if(contigs_info->c_info_vt[abs(*lit)].unique&&(!contigs_info->c_info_vt[abs(*lit)].used))
										if((dist_ctg_it->first>dist_searched-DistVar2)&&((dist_ctg_it->first+contigs_info->contig_sz_vt[abs(*lit)])<(dist_searched+Dist2NextUnitig+DistVar2)))
										{
											if(abs(*lit)==abs(NextUnitig))
											{
												break;
											}
											if(unitig_dist.count(abs(*lit))==0)
											{
												if(it==1)
												{
													t_right_ctgs.push_back(*lit);
												
												}
												else
												{
													t_left_ctgs.push_back(*lit);
												
												}
												if(it==1)
												{
													unitig_dist[abs(*lit)]=dist_ctg_it->first;
												}
												else
												{
													unitig_dist[abs(*lit)]=-(dist_ctg_it->first);
												}
												contigs_info->c_info_vt[abs(*lit)].used=1;
												//unitig_dist[abs(*lit)]=dist_searched+Dist2NextUnitig/2;
												key_unitigs[abs(*lit)]=dist_searched+Dist2NextUnitig/2;
											}
										
										}
										lit=n_lit;
						
									}
								
									dist_ctg_it=dist_ctg_it_n;
									continue;



								}








								Right=0;
								if((NextUnitig>0&&RightUnitig)||(NextUnitig<0&&RightUnitig==0))
								{
									Right=1;
								}
								UniqueNxt=1;
								nxt_ctg=NextUnitig;
								if(contigs_info->c_info_vt[abs(NextUnitig)].used==1)
								{break;}
							
								key_unitigs[abs(NextUnitig)]=dist_searched+Dist2NextUnitig;
								map<int,vector<int> > node_cov;
								

								if(it==1)
								{


									if(RightUnitig)
									{
										t_right_ctgs.push_back(NextUnitig);
									
										list<int> ctg_stack;
										ctg_stack.push_back(current_ctg);
									
										for(int d=3;d<max_dep;d+=2)
										{

											BFSearchPathFinder(contigs_info, ctg_stack,dist_ctg,key_unitigs,dist_searched,node_cov,d,1,MatePair);//,sc_cnt
										
											if(node_cov.size()==1)
											{
											
												break;
										
											}	
										
									

										}
									
							
									}
									else
									{
										t_right_ctgs.push_back(-NextUnitig);

										list<int> ctg_stack;
										ctg_stack.push_back(-current_ctg);
									
										for(int d=3;d<max_dep;d+=2)
										{

											BFSearchPathFinder(contigs_info, ctg_stack,dist_ctg,key_unitigs,dist_searched,node_cov,d,1,MatePair);//,sc_cnt
											if(node_cov.size()==1)
											{
											
												break;
										
											}	
										
									


										}
								

									}

								



									if(node_cov.size()==1)
									{
										for(int ii=0;ii<t_right_ctgs.size();++ii)
										{
											if(t_right_ctgs[ii]==((*node_cov.begin()).first))
											{
												vector<int> tmp_vt=(*node_cov.begin()).second;
												for(int jj=1;jj<tmp_vt.size()-1;++jj)
												{
													if(UniqueCurrent&&jj==1)
													{continue;}
													
													right_ctgs.push_back(tmp_vt[jj]);
												}
												for(int jj=ii;jj<t_right_ctgs.size();++jj)
												{
													right_ctgs.push_back(t_right_ctgs[jj]);
												}
												break;
											}
											if(t_right_ctgs[ii]==-((*node_cov.begin()).first))
											{
												vector<int> tmp_vt=(*node_cov.begin()).second;
												for(int jj=1;jj<tmp_vt.size()-1;++jj)
												{
													if(UniqueCurrent&&jj==1)
													{continue;}
													
													right_ctgs.push_back(-tmp_vt[jj]);
												}
												for(int jj=ii;jj<t_right_ctgs.size();++jj)
												{
													right_ctgs.push_back(t_right_ctgs[jj]);
												}
												break;
											}
										}
									//	break;
										
									}
									else
									{
										for(int jj=0;jj<t_right_ctgs.size();++jj)
										{
											right_ctgs.push_back(t_right_ctgs[jj]);
										}
									}
										
									t_right_ctgs.clear();
								}
								else
								{
									if(RightUnitig)
									{
										t_left_ctgs.push_back(-NextUnitig);

											list<int> ctg_stack;
										ctg_stack.push_back(current_ctg);
								
										for(int d=3;d<max_dep;d+=2)
										{

											BFSearchPathFinder(contigs_info, ctg_stack,dist_ctg,key_unitigs,dist_searched,node_cov,d,1,MatePair);//,sc_cnt
											if(node_cov.size()==1)
											{
											
												break;
										
											}
										}
								
									}
									else
									{
										t_left_ctgs.push_back(NextUnitig);


											list<int> ctg_stack;
										ctg_stack.push_back(-current_ctg);
								
										for(int d=3;d<max_dep;d+=2)
										{

											BFSearchPathFinder(contigs_info, ctg_stack,dist_ctg,key_unitigs,dist_searched,node_cov,d,1,MatePair);//,sc_cnt
											if(node_cov.size()==1)
											{
											
												break;
										
											}	
										
									
										}
								
									}
								

									if(node_cov.size()==1)
									{
										for(int ii=0;ii<t_left_ctgs.size();++ii)
										{
											if(t_left_ctgs[ii]==((*node_cov.begin()).first))
											{
												vector<int> tmp_vt=(*node_cov.begin()).second;
												for(int jj=1;jj<tmp_vt.size()-1;++jj)
												{
													if(UniqueCurrent&&jj==1)
													{continue;}
													
													left_ctgs.push_back(tmp_vt[jj]);
												}
												for(int jj=ii;jj<t_left_ctgs.size();++jj)
												{
													left_ctgs.push_back(t_left_ctgs[jj]);
												}
												break;
											}
											if(t_left_ctgs[ii]==-((*node_cov.begin()).first))
											{
												vector<int> tmp_vt=(*node_cov.begin()).second;
												for(int jj=1;jj<tmp_vt.size()-1;++jj)
												{
													if(UniqueCurrent&&jj==1)
													{continue;}
													
													left_ctgs.push_back(-tmp_vt[jj]);
												}
												for(int jj=ii;jj<t_left_ctgs.size();++jj)
												{
													left_ctgs.push_back(t_left_ctgs[jj]);
												}
												break;
											}
										}
									//	break;
										
									}
									else
									{
										for(int jj=0;jj<t_left_ctgs.size();++jj)
										{
											left_ctgs.push_back(t_left_ctgs[jj]);
										}
									}
								
									t_left_ctgs.clear();
								}

							
								contigs_info->c_info_vt[abs(NextUnitig)].used=1;
								NextUnitig=0;
								dist_searched+=Dist2NextUnitig;
								continue;
							}
							break;
						
					
						}


						if(marked_cnt>1)
						{


						
							map<int,int> key_unitigs;
							t_right_ctgs.clear();
							t_left_ctgs.clear();
							if(Exist_Next_Unitig)
							{


								map<int,list<int> >::iterator dist_ctg_it,dist_ctg_it_n;
							
								for(dist_ctg_it=dist_ctg.begin();dist_ctg_it!=dist_ctg.end();)
								{
									dist_ctg_it_n=dist_ctg_it;
									dist_ctg_it_n++;
						
									list<int>::iterator lit,lit2,n_lit;
									for(lit=dist_ctg_it->second.begin();lit!=dist_ctg_it->second.end();)
									{
										lit2=lit;
										lit2++;
										n_lit=lit2;
										n_lit++;
										int dist0=*lit2;
						
										int DistVar2=DistVar;
													
										if(MatePair)//dist0>3000)
										{
											DistVar2=10000;//dist0;
										}
										if(contigs_info->c_info_vt[abs(*lit)].unique&&(!contigs_info->c_info_vt[abs(*lit)].used))
										if((dist_ctg_it->first>dist_searched-DistVar2)&&((dist_ctg_it->first+contigs_info->contig_sz_vt[abs(*lit)])<(dist_searched+Dist2NextUnitig+DistVar2)))
										{
											if(abs(*lit)==abs(NextUnitig))
											{
												break;
											}
											if(unitig_dist.count(abs(*lit))==0)
											{
												if(it==1)
												{
													t_right_ctgs.push_back(*lit);
														
												}
												else
												{
													t_left_ctgs.push_back(*lit);
														
												}
												if(it==1)
												{
													unitig_dist[abs(*lit)]=dist_ctg_it->first;
												}
												else
												{
													unitig_dist[abs(*lit)]=-(dist_ctg_it->first);
												}
												contigs_info->c_info_vt[abs(*lit)].used=1;
												//unitig_dist[abs(*lit)]=dist_searched+Dist2NextUnitig/2;
												key_unitigs[abs(*lit)]=dist_searched+Dist2NextUnitig/2;
											}
												
										}
										lit=n_lit;
						
									}
								
									dist_ctg_it=dist_ctg_it_n;
									continue;



								}







								Right=0;
								if((NextUnitig>0&&RightUnitig)||(NextUnitig<0&&RightUnitig==0))
								{
									Right=1;
								}
								UniqueNxt=1;
								nxt_ctg=NextUnitig;
								if(contigs_info->c_info_vt[abs(NextUnitig)].used==1)
								{break;}
								key_unitigs[abs(NextUnitig)]=dist_searched+Dist2NextUnitig;
								map<int,vector<int> > node_cov;
								

								if(it==1)
								{


									if(RightUnitig)
									{
										t_right_ctgs.push_back(NextUnitig);
									
										list<int> ctg_stack;
										ctg_stack.push_back(current_ctg);
									
										for(int d=3;d<max_dep;d+=2)
										{

											BFSearchPathFinder(contigs_info, ctg_stack,dist_ctg,key_unitigs,dist_searched,node_cov,d,1,MatePair);//,sc_cnt
											if(node_cov.size()==1)
											{
											
												break;
										
											}
										}
									
							
									}
									else
									{
										t_right_ctgs.push_back(-NextUnitig);

										list<int> ctg_stack;
										ctg_stack.push_back(-current_ctg);
									
										for(int d=3;d<max_dep;d+=2)
										{

											BFSearchPathFinder(contigs_info, ctg_stack,dist_ctg,key_unitigs,dist_searched,node_cov,d,1,MatePair);//,sc_cnt
											if(node_cov.size()==1)
											{
										
												break;
										
											}
										}
								

									}



								



									if(node_cov.size()==1)
									{
										for(int ii=0;ii<t_right_ctgs.size();++ii)
										{
											if(t_right_ctgs[ii]==((*node_cov.begin()).first))
											{
												vector<int> tmp_vt=(*node_cov.begin()).second;
												for(int jj=1;jj<tmp_vt.size()-1;++jj)
												{
													if(UniqueCurrent&&jj==1)
													{continue;}
													
													right_ctgs.push_back(tmp_vt[jj]);
												}
												for(int jj=ii;jj<t_right_ctgs.size();++jj)
												{
													right_ctgs.push_back(t_right_ctgs[jj]);
												}
												break;
											}
											if(t_right_ctgs[ii]==-((*node_cov.begin()).first))
											{
												vector<int> tmp_vt=(*node_cov.begin()).second;
												for(int jj=1;jj<tmp_vt.size()-1;++jj)
												{
													if(UniqueCurrent&&jj==1)
													{continue;}
													
													right_ctgs.push_back(-tmp_vt[jj]);
												}
												for(int jj=ii;jj<t_right_ctgs.size();++jj)
												{
													right_ctgs.push_back(t_right_ctgs[jj]);
												}
												break;
											}
										}
									//	break;
										
									}	
									else
									{
										for(int jj=0;jj<t_right_ctgs.size();++jj)
										{
											right_ctgs.push_back(t_right_ctgs[jj]);
										}
									}		
									t_right_ctgs.clear();
								}
								else
								{
									if(RightUnitig)
									{
										t_left_ctgs.push_back(-NextUnitig);

											list<int> ctg_stack;
										ctg_stack.push_back(current_ctg);
									
										for(int d=3;d<max_dep;d+=2)
										{

											BFSearchPathFinder(contigs_info, ctg_stack,dist_ctg,key_unitigs,dist_searched,node_cov,d,1,MatePair);//,sc_cnt
											if(node_cov.size()==1)
											{
											
												break;
										
											}	
										
										}
								
									}
									else
									{
										t_left_ctgs.push_back(NextUnitig);


											list<int> ctg_stack;
										ctg_stack.push_back(-current_ctg);
									
										for(int d=3;d<max_dep;d+=2)
										{

											BFSearchPathFinder(contigs_info, ctg_stack,dist_ctg,key_unitigs,dist_searched,node_cov,d,1,MatePair);//,sc_cnt
											if(node_cov.size()==1)
											{
											
												break;
										
											}	
										
										}
								
									}



									if(node_cov.size()==1)
									{
										for(int ii=0;ii<t_left_ctgs.size();++ii)
										{
											if(t_left_ctgs[ii]==((*node_cov.begin()).first))
											{
												vector<int> tmp_vt=(*node_cov.begin()).second;
												for(int jj=1;jj<tmp_vt.size()-1;++jj)
												{
													if(UniqueCurrent&&jj==1)
													{continue;}
													
													left_ctgs.push_back(tmp_vt[jj]);
												}
												for(int jj=ii;jj<t_left_ctgs.size();++jj)
												{
													left_ctgs.push_back(t_left_ctgs[jj]);
												}
												break;
											}
											if(t_left_ctgs[ii]==-((*node_cov.begin()).first))
											{
												vector<int> tmp_vt=(*node_cov.begin()).second;
												for(int jj=1;jj<tmp_vt.size()-1;++jj)
												{
													if(UniqueCurrent&&jj==1)
													{continue;}
													
													left_ctgs.push_back(-tmp_vt[jj]);
												}
												for(int jj=ii;jj<t_left_ctgs.size();++jj)
												{
													left_ctgs.push_back(t_left_ctgs[jj]);
												}
												break;
											}
										}
									//	break;
										
									}
									else
									{
										for(int jj=0;jj<t_left_ctgs.size();++jj)
										{
											left_ctgs.push_back(t_left_ctgs[jj]);
										}
								
									}
								
								
									t_left_ctgs.clear();
								}


								contigs_info->c_info_vt[abs(NextUnitig)].used=1;
								//key_unitigs[abs(NextUnitig)]=dist_searched+Dist2NextUnitig;
								NextUnitig=0;
								dist_searched+=Dist2NextUnitig;
								continue;
							}
							break;
						
						}

						if(marked_cnt==1)
						{
						
							//dist_searched-=marked_ctgs[abs(nxt_ctg)]/2;//smoothing
							if(searched_ctgs.count(abs(nxt_ctg))==0)
							{
								searched_ctgs[abs(nxt_ctg)]++;
							
							}
							else
							{
								searched_ctgs[abs(nxt_ctg)]++;
								if(searched_ctgs[abs(nxt_ctg)]>4)
								{		
									map<int,int> key_unitigs;
									t_right_ctgs.clear();
									t_left_ctgs.clear();
									if(Exist_Next_Unitig)
									{


										map<int,list<int> >::iterator dist_ctg_it,dist_ctg_it_n;
									
										for(dist_ctg_it=dist_ctg.begin();dist_ctg_it!=dist_ctg.end();)
										{
											dist_ctg_it_n=dist_ctg_it;
											dist_ctg_it_n++;
						
											list<int>::iterator lit,lit2,n_lit;
											for(lit=dist_ctg_it->second.begin();lit!=dist_ctg_it->second.end();)
											{
												lit2=lit;
												lit2++;
												n_lit=lit2;
												n_lit++;
												int dist0=*lit2;
						
												int DistVar2=DistVar;
													
												if(MatePair)//dist0>3000)
												{
													DistVar2=10000;//dist0;
												}
												if(contigs_info->c_info_vt[abs(*lit)].unique&&(!contigs_info->c_info_vt[abs(*lit)].used))
												if((dist_ctg_it->first>dist_searched-DistVar2)&&((dist_ctg_it->first+contigs_info->contig_sz_vt[abs(*lit)])<(dist_searched+Dist2NextUnitig+DistVar2)))
												{
													if(abs(*lit)==abs(NextUnitig))
													{
														break;
													}
													if(unitig_dist.count(abs(*lit))==0)
													{
														if(it==1)
														{
															right_ctgs.push_back(*lit);
														
														}
														else
														{
															left_ctgs.push_back(*lit);
														
														}
														if(it==1)
														{
															unitig_dist[abs(*lit)]=dist_ctg_it->first;
														}
														else
														{
															unitig_dist[abs(*lit)]=-(dist_ctg_it->first);
														}
														contigs_info->c_info_vt[abs(*lit)].used=1;
														//unitig_dist[abs(*lit)]=dist_searched+Dist2NextUnitig/2;
														key_unitigs[abs(*lit)]=dist_searched+Dist2NextUnitig/2;
													}
												
												}
												lit=n_lit;
						
											}
								
											dist_ctg_it=dist_ctg_it_n;
											continue;



										}







										Right=0;
										if((NextUnitig>0&&RightUnitig)||(NextUnitig<0&&RightUnitig==0))
										{
											Right=1;
										}
										UniqueNxt=1;
										nxt_ctg=NextUnitig;
										if(contigs_info->c_info_vt[abs(NextUnitig)].used==1)
										{break;}

										key_unitigs[abs(NextUnitig)]=dist_searched+Dist2NextUnitig;
										map<int,vector<int> > node_cov;
								

										if(it==1)
										{


											if(RightUnitig)
											{
												t_right_ctgs.push_back(NextUnitig);
									
												list<int> ctg_stack;
												ctg_stack.push_back(current_ctg);
											
												for(int d=3;d<max_dep;d+=2)
												{

													BFSearchPathFinder(contigs_info, ctg_stack,dist_ctg,key_unitigs,dist_searched,node_cov,d,1,MatePair);//,sc_cnt
													if(node_cov.size()==1)
													{
												
														break;


													}
												}
									
							
											}
											else
											{
												t_right_ctgs.push_back(-NextUnitig);

												list<int> ctg_stack;
												ctg_stack.push_back(-current_ctg);
											
												for(int d=3;d<max_dep;d+=2)
												{

													BFSearchPathFinder(contigs_info, ctg_stack,dist_ctg,key_unitigs,dist_searched,node_cov,d,1,MatePair);//,sc_cnt
													if(node_cov.size()==1)
													{
												
														break;
												
													}
												}
								

											}



										



											if(node_cov.size()==1)
											{
												for(int ii=0;ii<t_right_ctgs.size();++ii)
												{
													if(t_right_ctgs[ii]==((*node_cov.begin()).first))
													{
														vector<int> tmp_vt=(*node_cov.begin()).second;
														for(int jj=1;jj<tmp_vt.size()-1;++jj)
														{
															if(UniqueCurrent&&jj==1)
															{continue;}
													
															right_ctgs.push_back(tmp_vt[jj]);
														}
														for(int jj=ii;jj<t_right_ctgs.size();++jj)
														{
															right_ctgs.push_back(t_right_ctgs[jj]);
														}
														break;
													}
													if(t_right_ctgs[ii]==-((*node_cov.begin()).first))
													{
														vector<int> tmp_vt=(*node_cov.begin()).second;
														for(int jj=1;jj<tmp_vt.size()-1;++jj)
														{
															if(UniqueCurrent&&jj==1)
															{continue;}
													
															right_ctgs.push_back(-tmp_vt[jj]);
														}
														for(int jj=ii;jj<t_right_ctgs.size();++jj)
														{
															right_ctgs.push_back(t_right_ctgs[jj]);
														}
														break;
													}
												}
											//	break;
										
											}
											else
											{
												for(int jj=0;jj<t_right_ctgs.size();++jj)
												{
													right_ctgs.push_back(t_right_ctgs[jj]);
												}
											}
										
								
											t_right_ctgs.clear();
										}
										else
										{
											if(RightUnitig)
											{
												t_left_ctgs.push_back(-NextUnitig);

													list<int> ctg_stack;
												ctg_stack.push_back(current_ctg);
											
												for(int d=3;d<max_dep;d+=2)
												{

													BFSearchPathFinder(contigs_info, ctg_stack,dist_ctg,key_unitigs,dist_searched,node_cov,d,1,MatePair);//,sc_cnt
													if(node_cov.size()==1)
													{
													
														break;
										
													}	
										
												}
								
											}
											else
											{
												t_left_ctgs.push_back(NextUnitig);


													list<int> ctg_stack;
												ctg_stack.push_back(-current_ctg);
											
												for(int d=3;d<max_dep;d+=2)
												{

													BFSearchPathFinder(contigs_info, ctg_stack,dist_ctg,key_unitigs,dist_searched,node_cov,d,1,MatePair);//,sc_cnt
													if(node_cov.size()==1)
													{
													
														break;
										
													}	
										
												}
								
											}

											if(node_cov.size()==1)
											{
												for(int ii=0;ii<t_left_ctgs.size();++ii)
												{
													if(t_left_ctgs[ii]==((*node_cov.begin()).first))
													{
														vector<int> tmp_vt=(*node_cov.begin()).second;
														for(int jj=1;jj<tmp_vt.size()-1;++jj)
														{
															if(UniqueCurrent&&jj==1)
															{continue;}
													
															left_ctgs.push_back(tmp_vt[jj]);
														}
														for(int jj=ii;jj<t_left_ctgs.size();++jj)
														{
															left_ctgs.push_back(t_left_ctgs[jj]);
														}
														break;
													}
													if(t_left_ctgs[ii]==-((*node_cov.begin()).first))
													{
														vector<int> tmp_vt=(*node_cov.begin()).second;
														for(int jj=1;jj<tmp_vt.size()-1;++jj)
														{
															if(UniqueCurrent&&jj==1)
															{continue;}
													
															left_ctgs.push_back(-tmp_vt[jj]);
														}
														for(int jj=ii;jj<t_left_ctgs.size();++jj)
														{
															left_ctgs.push_back(t_left_ctgs[jj]);
														}
														break;
													}
												}
											//	break;
										
											}
											else
											{
												for(int jj=0;jj<t_left_ctgs.size();++jj)
												{
													left_ctgs.push_back(t_left_ctgs[jj]);
												}
											}
								
											t_left_ctgs.clear();
										}
									
									
										contigs_info->c_info_vt[abs(NextUnitig)].used=1;
										//key_unitigs[abs(NextUnitig)]=dist_searched+Dist2NextUnitig;
										NextUnitig=0;
										dist_searched+=Dist2NextUnitig;
										continue;
									}
									break;
								}
							}



							if(nxt_ctg>0)
							{
								UniqueNxt=contigs_info->c_info_vt[nxt_ctg].unique;

								if(UniqueNxt&&contigs_info->c_info_vt[nxt_ctg].used)
								{
									break;
									contigs_info->c_info_vt[nxt_ctg].unique=0;
								//	break;
								}

								Right=0;
								if(it==1)
								{
									t_right_ctgs.push_back(-nxt_ctg);
								}
								else
								{
									t_left_ctgs.push_back(nxt_ctg);
								}
								//contigs_info->c_info_vt[abs(nxt_ctg)].used=1;
								//left search next
							}
							else
							{
								UniqueNxt=contigs_info->c_info_vt[-nxt_ctg].unique;
								if(UniqueNxt&&contigs_info->c_info_vt[-nxt_ctg].used)
								{
									contigs_info->c_info_vt[-nxt_ctg].unique=0;
									break;
								}
								Right=1;
								if(it==1)
								{
									t_right_ctgs.push_back(-nxt_ctg);
								}
								else
								{
									t_left_ctgs.push_back(nxt_ctg);
								}
								//contigs_info->c_info_vt[abs(nxt_ctg)].used=1;
								//right search next

							}

							if(UniqueNxt)
							{
								if(it==1)
								{
									for (int ii=0;ii<t_right_ctgs.size();++ii)
									{
										right_ctgs.push_back(t_right_ctgs[ii]);
										contigs_info->c_info_vt[abs(t_right_ctgs[ii])].used=1;
									}
									t_right_ctgs.clear();
								}
								else
								{
									for (int ii=0;ii<t_left_ctgs.size();++ii)
									{
										left_ctgs.push_back(t_left_ctgs[ii]);
										contigs_info->c_info_vt[abs(t_left_ctgs[ii])].used=1;
									}
									t_left_ctgs.clear();
							
								}
							}


						}
						else
						{
							map<int,int> key_unitigs;
							t_right_ctgs.clear();
							t_left_ctgs.clear();


							if(Exist_Next_Unitig)
							{


								map<int,list<int> >::iterator dist_ctg_it,dist_ctg_it_n;
							
								for(dist_ctg_it=dist_ctg.begin();dist_ctg_it!=dist_ctg.end();)
								{
									dist_ctg_it_n=dist_ctg_it;
									dist_ctg_it_n++;
						
									list<int>::iterator lit,lit2,n_lit;
									for(lit=dist_ctg_it->second.begin();lit!=dist_ctg_it->second.end();)
									{
										lit2=lit;
										lit2++;
										n_lit=lit2;
										n_lit++;
										int dist0=*lit2;
						
										int DistVar2=DistVar;
													
										if(MatePair)//dist0>3000)
										{
											DistVar2=10000;//dist0;
										}
										if(contigs_info->c_info_vt[abs(*lit)].unique&&(!contigs_info->c_info_vt[abs(*lit)].used))
										if((dist_ctg_it->first>dist_searched-DistVar2)&&((dist_ctg_it->first+contigs_info->contig_sz_vt[abs(*lit)])<(dist_searched+Dist2NextUnitig+DistVar2)))
										{
											if(abs(*lit)==abs(NextUnitig))
											{
												break;
											}
											if(unitig_dist.count(abs(*lit))==0)
											{
												if(it==1)
												{
													t_right_ctgs.push_back(*lit);
												
												}
												else
												{
													t_left_ctgs.push_back(*lit);
												
												}
												if(it==1)
												{
													unitig_dist[abs(*lit)]=dist_ctg_it->first;
												}
												else
												{
													unitig_dist[abs(*lit)]=-(dist_ctg_it->first);
												}
												contigs_info->c_info_vt[abs(*lit)].used=1;
												//unitig_dist[abs(*lit)]=dist_searched+Dist2NextUnitig/2;
												key_unitigs[abs(*lit)]=dist_searched+Dist2NextUnitig/2;
											}
										
										}
										lit=n_lit;
						
									}
								
									dist_ctg_it=dist_ctg_it_n;
									continue;



								}







								Right=0;
								if((NextUnitig>0&&RightUnitig)||(NextUnitig<0&&RightUnitig==0))
								{
									Right=1;
								}
								UniqueNxt=1;
								nxt_ctg=NextUnitig;
								if(contigs_info->c_info_vt[abs(NextUnitig)].used==1)
								{break;}

								key_unitigs[abs(NextUnitig)]=dist_searched+Dist2NextUnitig;
								map<int,vector<int> > node_cov;
								

								if(it==1)
								{


									if(RightUnitig)
									{
										t_right_ctgs.push_back(NextUnitig);
									
										list<int> ctg_stack;
										ctg_stack.push_back(current_ctg);
									
										for(int d=3;d<max_dep;d+=2)
										{

											BFSearchPathFinder(contigs_info, ctg_stack,dist_ctg,key_unitigs,dist_searched,node_cov,d,1,MatePair);//,sc_cnt
											if(node_cov.size()==1)
											{
										
												break;

											}
										}
									
							
									}
									else
									{
										t_right_ctgs.push_back(-NextUnitig);

										list<int> ctg_stack;
										ctg_stack.push_back(-current_ctg);
									
										for(int d=3;d<max_dep;d+=2)
										{

											BFSearchPathFinder(contigs_info, ctg_stack,dist_ctg,key_unitigs,dist_searched,node_cov,d,1,MatePair);//,sc_cnt
											if(node_cov.size()==1)
											{
										
												break;

											}
										}
								

									}
								



								



									if(node_cov.size()==1)
									{
										for(int ii=0;ii<t_right_ctgs.size();++ii)
										{
											if(t_right_ctgs[ii]==((*node_cov.begin()).first))
											{
												vector<int> tmp_vt=(*node_cov.begin()).second;
												for(int jj=1;jj<tmp_vt.size()-1;++jj)
												{
													if(UniqueCurrent&&jj==1)
													{continue;}
													right_ctgs.push_back(tmp_vt[jj]);
												}
												for(int jj=ii;jj<t_right_ctgs.size();++jj)
												{
													right_ctgs.push_back(t_right_ctgs[jj]);
												}
												break;
											}
											if(t_right_ctgs[ii]==-((*node_cov.begin()).first))
											{
												vector<int> tmp_vt=(*node_cov.begin()).second;
												for(int jj=1;jj<tmp_vt.size()-1;++jj)
												{
													if(UniqueCurrent&&jj==1)
													{continue;}
													right_ctgs.push_back(-tmp_vt[jj]);
												}
												for(int jj=ii;jj<t_right_ctgs.size();++jj)
												{
													right_ctgs.push_back(t_right_ctgs[jj]);
												}
												break;
											}
										}
									//	break;
										
									}	
									else
									{
										for(int jj=0;jj<t_right_ctgs.size();++jj)
										{
											right_ctgs.push_back(t_right_ctgs[jj]);
										}
									}
										
									t_right_ctgs.clear();
								}
								else
								{
									if(RightUnitig)
									{
										t_left_ctgs.push_back(-NextUnitig);

											list<int> ctg_stack;
										ctg_stack.push_back(current_ctg);
									
										for(int d=3;d<max_dep;d+=2)
										{

											BFSearchPathFinder(contigs_info, ctg_stack,dist_ctg,key_unitigs,dist_searched,node_cov,d,1,MatePair);//,sc_cnt
											if(node_cov.size()==1)
											{
											
												break;
										
											}	
										
										}
								
									}
									else
									{
										t_left_ctgs.push_back(NextUnitig);


											list<int> ctg_stack;
										ctg_stack.push_back(-current_ctg);
								
										for(int d=3;d<max_dep;d+=2)
										{

											BFSearchPathFinder(contigs_info, ctg_stack,dist_ctg,key_unitigs,dist_searched,node_cov,d,1,MatePair);//,sc_cnt
											if(node_cov.size()==1)
											{
											
												break;
										
											}	
										
										}
								
									}


									if(node_cov.size()==1)
									{
										for(int ii=0;ii<t_left_ctgs.size();++ii)
										{
											if(t_left_ctgs[ii]==((*node_cov.begin()).first))
											{
												vector<int> tmp_vt=(*node_cov.begin()).second;
												for(int jj=1;jj<tmp_vt.size()-1;++jj)
												{
													if(UniqueCurrent&&jj==1)
													{continue;}
													left_ctgs.push_back(tmp_vt[jj]);
												}
												for(int jj=ii;jj<t_left_ctgs.size();++jj)
												{
													left_ctgs.push_back(t_left_ctgs[jj]);
												}
												break;
											}
											if(t_left_ctgs[ii]==-((*node_cov.begin()).first))
											{
												vector<int> tmp_vt=(*node_cov.begin()).second;
												for(int jj=1;jj<tmp_vt.size()-1;++jj)
												{
													if(UniqueCurrent&&jj==1)
													{continue;}
													left_ctgs.push_back(-tmp_vt[jj]);
												}
												for(int jj=ii;jj<t_left_ctgs.size();++jj)
												{
													left_ctgs.push_back(t_left_ctgs[jj]);
												}
												break;
											}
										}
									//	break;
										
									}
									else
									{
										for(int jj=0;jj<t_left_ctgs.size();++jj)
										{
											left_ctgs.push_back(t_left_ctgs[jj]);
										}
									}
								
								
									t_left_ctgs.clear();
								}
								//key_unitigs[abs(NextUnitig)]=dist_searched+Dist2NextUnitig;
								contigs_info->c_info_vt[abs(NextUnitig)].used=1;
								NextUnitig=0;
								dist_searched+=Dist2NextUnitig;
								continue;
							}

							break;
						
					

						}



					}


				}

			}


			if(0)//round==1)
			{
				//cleaning stuffs;
			//	vector<int> left_ctgs,right_ctgs;
				vector<int> unitigs_vt;
				map<int,vector<int> > ctgs_dist_cov;

				reverse(left_ctgs.begin(),left_ctgs.end());

				for(int j=0;j<(int)left_ctgs.size();++j)
				{
					int ctg_no=left_ctgs[j];
					if(contigs_info->c_info_vt[abs(ctg_no)].unique)
					{
						unitigs_vt.push_back(ctg_no);
					}
					contigs_info->c_info_vt[abs(ctg_no)].used=0;
					

				}
				if(contigs_info->c_info_vt[abs(contig_no)].unique)
				{
					unitigs_vt.push_back(contig_no);
				}
				contigs_info->c_info_vt[abs(contig_no)].used=0;


				for(int j=0;j<(int)right_ctgs.size();++j)
				{
					int ctg_no=right_ctgs[j];
					if(contigs_info->c_info_vt[abs(ctg_no)].unique)
					{
						unitigs_vt.push_back(ctg_no);
					}
					contigs_info->c_info_vt[abs(ctg_no)].used=0;

				}

				//add weak contigs.
				//cout<<"adding weak ctgs."<<endl;
				for(int j=0;j<unitigs_vt.size();++j)
				{
					int ctg_no=unitigs_vt[j];

					map<int,struct scaffold_contig_info>::iterator scf_it;
					for(scf_it=contigs_info->scaffold_adjacency_left[abs(ctg_no)].begin();scf_it!=contigs_info->scaffold_adjacency_left[abs(ctg_no)].end();++scf_it)
					{
						if(ctg_no>0)
						{
							if(ctgs_dist_cov.count(scf_it->first)==0)
							{
								ctgs_dist_cov[scf_it->first].resize(2);
							}
							
							if(unitig_dist[abs(ctg_no)]<=0)
							{
								ctgs_dist_cov[scf_it->first][0]-=(contigs_info->contig_sz_vt[abs(ctg_no)])*scf_it->second.cov;
								ctgs_dist_cov[scf_it->first][0]-=(((abs(unitig_dist[abs(ctg_no)]))*scf_it->second.cov)-scf_it->second.dist_sum);
						
							}
							else
							{
								ctgs_dist_cov[scf_it->first][0]-=(contigs_info->contig_sz_vt[abs(scf_it->first)])*scf_it->second.cov;
								ctgs_dist_cov[scf_it->first][0]+=(((abs(unitig_dist[abs(ctg_no)]))*scf_it->second.cov)+scf_it->second.dist_sum);
				
							}
				
								
							ctgs_dist_cov[scf_it->first][1]+=scf_it->second.cov;
						
							
						}
						else
						{
							if(ctgs_dist_cov.count(-scf_it->first)==0)
							{
								ctgs_dist_cov[-scf_it->first].resize(2);
							}
							if(unitig_dist[abs(ctg_no)]>=0)
							{
								ctgs_dist_cov[-scf_it->first][0]+=(contigs_info->contig_sz_vt[abs(ctg_no)])*scf_it->second.cov;
								ctgs_dist_cov[-scf_it->first][0]+=(((abs(unitig_dist[abs(ctg_no)]))*scf_it->second.cov)+scf_it->second.dist_sum);
						
							}
							else
							{
								ctgs_dist_cov[-scf_it->first][0]+=(contigs_info->contig_sz_vt[abs(scf_it->first)])*scf_it->second.cov;
								ctgs_dist_cov[-scf_it->first][0]-=(((abs(unitig_dist[abs(ctg_no)]))*scf_it->second.cov)+scf_it->second.dist_sum);
				
							}
				
								
							ctgs_dist_cov[-scf_it->first][1]+=scf_it->second.cov;
							
						}
					}


					for(scf_it=contigs_info->scaffold_adjacency_right[abs(ctg_no)].begin();scf_it!=contigs_info->scaffold_adjacency_right[abs(ctg_no)].end();++scf_it)
					{
						if(ctg_no<0)
						{
							if(ctgs_dist_cov.count(-scf_it->first)==0)
							{
								ctgs_dist_cov[-scf_it->first].resize(2);
							}
							
							if(unitig_dist[abs(ctg_no)]<=0)
							{
								ctgs_dist_cov[-scf_it->first][0]-=(contigs_info->contig_sz_vt[abs(ctg_no)])*scf_it->second.cov;
								ctgs_dist_cov[-scf_it->first][0]-=(((abs(unitig_dist[abs(ctg_no)]))*scf_it->second.cov)+scf_it->second.dist_sum);
				
										
							}
							else
							{
								ctgs_dist_cov[-scf_it->first][0]-=(contigs_info->contig_sz_vt[abs(scf_it->first)])*scf_it->second.cov;
								ctgs_dist_cov[-scf_it->first][0]+=(((abs(unitig_dist[abs(ctg_no)]))*scf_it->second.cov)+scf_it->second.dist_sum);
				
							}
				
								
							ctgs_dist_cov[-scf_it->first][1]+=scf_it->second.cov;
						
							
						}
						else
						{
							if(ctgs_dist_cov.count(scf_it->first)==0)
							{
								ctgs_dist_cov[scf_it->first].resize(2);
							}
							if(unitig_dist[abs(ctg_no)]>=0)
							{
								ctgs_dist_cov[scf_it->first][0]+=(contigs_info->contig_sz_vt[abs(ctg_no)])*scf_it->second.cov;
								ctgs_dist_cov[scf_it->first][0]+=(((abs(unitig_dist[abs(ctg_no)]))*scf_it->second.cov)+scf_it->second.dist_sum);
						
							}
							else
							{
								ctgs_dist_cov[scf_it->first][0]+=(contigs_info->contig_sz_vt[abs(scf_it->first)])*scf_it->second.cov;
								ctgs_dist_cov[scf_it->first][0]-=(((abs(unitig_dist[abs(ctg_no)]))*scf_it->second.cov)+scf_it->second.dist_sum);
				
							}
				
								
							ctgs_dist_cov[scf_it->first][1]+=scf_it->second.cov;
							
						}
					}


				}
				
				//cout<<"adding missing contigs ."<<endl;
				map<int,vector<int> >::iterator ctgs_it,ctgs_it_n;
				for(ctgs_it=ctgs_dist_cov.begin();ctgs_it!=ctgs_dist_cov.end();++ctgs_it)
				{

					if(ctgs_it->second[1]>=LinkCovTh)
					{
						ctgs_for_second_round[ctgs_it->first]=ctgs_it->second[0]/ctgs_it->second[1];
					}
				}


				
				//cout<<"adding missing unitigs ."<<endl;

				map<int,int>::iterator unitig_it;
				//for(unitig_it=unitig_dist.begin();unitig_it!=unitig_dist.end();++unitig_it)
				for(int j=0;j<unitigs_vt.size();++j)
				{
					ctgs_for_second_round[unitigs_vt[j]]=unitig_dist[abs(unitigs_vt[j])];

				}

				map<int, int > dist_ctg;
				
				
				//cout<<"dist_ctgs ."<<endl;
				map<int,int>::iterator ctgs_it2;
				for(ctgs_it2=ctgs_for_second_round.begin();ctgs_it2!=ctgs_for_second_round.end();++ctgs_it2)
				{
					int dist=ctgs_it2->second;
					dist_ctg[dist]=(ctgs_it2->first);
					
				}
				
				for(ctgs_it=ctgs_dist_cov.begin();ctgs_it!=ctgs_dist_cov.end();++ctgs_it)
				{
					
					
					if(ctgs_it->second[1]<LinkCovTh&&(ctgs_it->second[1]>LinkCovTh0&&contigs_info->c_info_vt[abs(ctgs_it->first)].unique))
					{
						int dist=ctgs_it->second[0]/ctgs_it->second[1];
						dist_ctg[dist]=ctgs_it->first;
					}
					
				}

				
				//cout<<"conflicts ."<<endl;

				map<int,int>::iterator dist_ctg_it_p,dist_ctg_it,dist_ctg_it_n;
				dist_ctg_it_p=dist_ctg.begin();
				dist_ctg_it=dist_ctg_it_p;
				cout<<dist_ctg.size()<<endl;
				
				if(dist_ctg.size()>0)
				{
				
					dist_ctg_it++;
				
					for (;dist_ctg_it!=dist_ctg.end();++dist_ctg_it)
					{
						
						dist_ctg_it_n=dist_ctg_it;
						dist_ctg_it_n++;
						if(dist_ctg_it_n==dist_ctg.end())
						{break;}

						if(ctgs_dist_cov.count(dist_ctg_it->second)==1)
						{
						
							if(ctgs_dist_cov[(dist_ctg_it->second)][1]<LinkCovTh&&ctgs_dist_cov[(dist_ctg_it->second)][1]>LinkCovTh0)
							{
								if(dist_ctg_it->first>0)
								{
									if(((dist_ctg_it_p->first)+contigs_info->contig_sz_vt[abs(dist_ctg_it_p->second)]-DistVar0<(dist_ctg_it->first))&&((dist_ctg_it->first)+contigs_info->contig_sz_vt[abs(dist_ctg_it->second)]-DistVar0<(dist_ctg_it_n->first)))
									{
										ctgs_for_second_round[dist_ctg_it->second]=dist_ctg_it->first;
									}
								}
								if(dist_ctg_it->first<0)
								{
									if((dist_ctg_it->first-contigs_info->contig_sz_vt[abs(dist_ctg_it->second)]+DistVar0>(dist_ctg_it_p->first))&&((dist_ctg_it_n->first)-contigs_info->contig_sz_vt[abs(dist_ctg_it_n->second)]+DistVar0>(dist_ctg_it->first)))
									{
										ctgs_for_second_round[dist_ctg_it->second]=dist_ctg_it->first;
									}
								}
							
							}
						}
						dist_ctg_it_p=dist_ctg_it;
					}
			

					//cout<<"conflicts end."<<endl;
		
				}
			}






		}
	//	ofstream osc_BP_info("Scaffold_Break_Points_info.txt");
	//	ofstream osc_cov("SuperContigs_Cov.txt");
	//	osc_BP_info<<"Break points: "<<break_points<<endl;
	//	osc_BP_info<<"Gap break points: "<<gap_points<<endl;
	//	osc_BP_info<<"Ambiguous break points: "<<ambiguous_points<<endl;




		o_sc_info<<">SuperContig_"<<sc_cnt<<endl;

		o_sc<<">SuperContig_"<<sc_cnt<<endl;
		
		
		vector<int> ctgs_vt;
		reverse(left_ctgs.begin(),left_ctgs.end());

		
		for(int j=0;j<(int)left_ctgs.size();++j)
		{
		//	o_sc_info<<left_ctgs[j]<<" ";
			ctgs_vt.push_back(left_ctgs[j]);
		}

		//o_sc_info<<contig_no<<" ";//endl;
		ctgs_vt.push_back(contig_no);

		for(int j=0;j<(int)right_ctgs.size();++j)
		{
			//o_sc_info<<right_ctgs[j]<<" ";
			ctgs_vt.push_back(right_ctgs[j]);

		}
		//o_sc_info<<endl;
		

		//output coverage		

		uint64_t tot_cov=0,uni_cov=0,tot_len=0,uni_len=0;
		for(int j=0;j<(int)ctgs_vt.size();++j)
		{
			if(contigs_info->c_info_vt[abs(ctgs_vt[j])].unique)
			{
				tot_cov+=contigs_info->contig_sz_vt[abs(ctgs_vt[j])]*contigs_info->cov_vt[abs(ctgs_vt[j])];
				tot_len+=contigs_info->contig_sz_vt[abs(ctgs_vt[j])];
				uni_cov+=contigs_info->contig_sz_vt[abs(ctgs_vt[j])]*contigs_info->cov_vt[abs(ctgs_vt[j])];
				uni_len+=contigs_info->contig_sz_vt[abs(ctgs_vt[j])];
			}
			tot_cov+=contigs_info->contig_sz_vt[abs(ctgs_vt[j])]*contigs_info->cov_vt[abs(ctgs_vt[j])];
			tot_len+=contigs_info->contig_sz_vt[abs(ctgs_vt[j])];

		}
		if(uni_len==0)
		{
			if(tot_len>0)
			{
				osc_cov<<(tot_cov/tot_len)<<endl;;
			}
			else
			{
				osc_cov<<"0"<<endl;
			}
		}
		else
		{
			osc_cov<<(uni_cov/uni_len)<<endl;;
		}


		
		int gap_closing_flag=0;
		int gap_offset=0;
		//o_sc<<sc_cnt<<endl;
		o_sc_info<<ctgs_vt.size()<<": ";

		//contigs_info->scaffolds.push_back(ctgs_vt);
		//vector<int> empty_vt;
		//empty_vt.resize(0);
		//contigs_info->gaps_in_scaffolds.push_back(empty_vt);
		int scf_len=0;
		for(int j=0;j<(int)ctgs_vt.size();++j)
		{
		//	if(1)//abs(ctgs_vt[j])==662)
		//	{cout<<ctgs_vt[j];}

			int ctg_no=ctgs_vt[j];
			o_sc_info<<ctg_no<<" ";//endl;
			/*
			if(contigs_info->c_info_vt[ctg_no].unique)
			{
				if(ctg_no>0)
				{
					contigs_info->ctg_in_scf[ctg_no].push_back(sc_cnt);
					contigs_info->ctg_in_scf[ctg_no].push_back(scf_len);
				}
				else
				{
					contigs_info->ctg_in_scf[-ctg_no].push_back(-sc_cnt);
					contigs_info->ctg_in_scf[-ctg_no].push_back(scf_len);
				}
			}
			*/
			if(ctgs_vt[j]>0)
			{
			//	cout<<"+"<<endl;
				int cur_ctg=ctgs_vt[j];
				int nxt_ctg=0;
				if(j!=ctgs_vt.size()-1)
				{
					nxt_ctg=ctgs_vt[j+1];
				}

				int codB;
				if(j==0)
				{
					codB=0;//contigs_info->contigs_hp_b[cur_ctg];
				}
				else
				{
					codB=contigs_info->contigs_hp_b[cur_ctg]+K_size;
				}
				int codE;
				if(nxt_ctg==0)
				{
					codE=contigs_info->contig_sz_vt[cur_ctg];
				}
				else
				{
					codE=contigs_info->contigs_hp_e[cur_ctg];
				}

				if(gap_closing_flag)
				{
					codB=gap_offset;
					gap_offset=0;
					gap_closing_flag=0;
				}

				o_sc<<contigs_info->contigs_str[cur_ctg].substr(codB,K_size+codE-codB);
				if(nxt_ctg!=0)
				{
					if(contigs_info->contig_adjacency_right[cur_ctg].count(nxt_ctg)>0)
					{
						//cout<<contigs_info->contig_adjacency_right[cur_ctg][nxt_ctg].bridge;
						string bridge=contigs_info->contig_adjacency_right[cur_ctg][nxt_ctg].bridge;
						//o_sc<<contigs_info->contig_adjacency_right[cur_ctg][nxt_ctg].bridge;
						
						o_sc<<bridge;

						int zero=0;
						o_sc_info<<"d: "<<zero<<" ";
						
						//contigs_info->gaps_in_scaffolds[sc_cnt].push_back(0);
					
					}
					else
					{
						if(contigs_info->scaffold_adjacency_right[cur_ctg].count(nxt_ctg)>0)
						{
							//o_sc<<contigs_info->contigs_str[cur_ctg].substr(..,..);
							
							int gap_sz=contigs_info->scaffold_adjacency_right[cur_ctg][nxt_ctg].dist_sum;
							int gap_cov=contigs_info->scaffold_adjacency_right[cur_ctg][nxt_ctg].cov;
							gap_sz=gap_sz/gap_cov;

							//cout<<gap_sz<<endl;
							if(gap_sz<150&&gap_sz>-150)
							{
							
								int start_pos=codB+K_size+codE-codB;
								if(start_pos<0)
								{
									start_pos=0;
								}
								if(start_pos>(contigs_info->contigs_str[cur_ctg].size()))
								{
									start_pos=(contigs_info->contigs_str[cur_ctg].size());
								}

								o_sc<<contigs_info->contigs_str[cur_ctg].substr(start_pos,200);//output all that can be extended
								//gap//
								string match1_s=contigs_info->contigs_str[cur_ctg].substr(codE,1000);
								int match1_sz=match1_s.size();
								int m_k=11;
								map<string,int> match_pos;	
								for(int i=0;i<=match1_sz - m_k;++i)
								{
									string substr=match1_s.substr(i,m_k);
									match_pos[substr]=i;
								}
								string nxt_ctg_s=contigs_info->contigs_str[abs(nxt_ctg)];
								if(nxt_ctg<0)
								{
									reverse(nxt_ctg_s.begin(),nxt_ctg_s.end());
									complement_str(nxt_ctg_s);
								}
								nxt_ctg_s=nxt_ctg_s.substr(0,200);
								int query_sz=nxt_ctg_s.size();
								map<int,int> gapoffset_cov;
								for(int i=0;i<query_sz-m_k;++i)
								{
									string substr=nxt_ctg_s.substr(i,m_k);
									if(match_pos.count(substr))
									{
										gap_closing_flag=1;
										gap_offset=i+match1_sz-match_pos[substr];	
										gapoffset_cov[gap_offset]++;
										if(gap_offset>=query_sz)
										{
											gap_closing_flag=0;
											gap_offset=0;

										}
										//break;
								
									}
								}
								map<int,int>::iterator offset_it;
								int best_offset=0,max_cov=0;
								for(offset_it=gapoffset_cov.begin();offset_it!=gapoffset_cov.end();++offset_it)
								{
									if((offset_it->second)>max_cov)
									{
										max_cov=offset_it->second;
										best_offset=offset_it->first;
									}

									
								}
								gap_offset=best_offset;
								if(gap_closing_flag==0)
								{
									
									//gap_closing_flag=1;
									gap_offset=0;
								}

							}
							
							


							if(gap_closing_flag==0)
							{
								for(int gg=0;gg<gap_sz;++gg)
								{o_sc<<'N';}
								if(gap_sz<0)
								{o_sc<<'N';}
							//	contigs_info->gaps_in_scaffolds[sc_cnt].push_back(gap_sz);
								
								o_sc_info<<"d: "<<gap_sz<<" ";
								scf_len+=gap_sz;
							}
							else
							{
								//contigs_info->gaps_in_scaffolds[sc_cnt].push_back(0);
								int zero=0;
								o_sc_info<<"d: "<<zero<<" ";
							}


						}
						else
						{
							//gap

						}
					}
				}
			}
			else
			{
			//	cout<<"-"<<endl;
				int cur_ctg=-ctgs_vt[j];
				int nxt_ctg=0;
				if(j!=ctgs_vt.size()-1)
				{
					nxt_ctg=ctgs_vt[j+1];
				}

				int codB;
				if(nxt_ctg==0)
				{
					codB=0;//contigs_info->contigs_hp_b[cur_ctg];
				}
				else
				{
					codB=contigs_info->contigs_hp_b[cur_ctg];
				}
				int codE;
				if(j==0)
				{
					codE=contigs_info->contig_sz_vt[cur_ctg];
				}
				else
				{
					codE=contigs_info->contigs_hp_e[cur_ctg];//-K_size+1
				}

				if(gap_closing_flag)
				{
					codE=contigs_info->contig_sz_vt[cur_ctg]-K_size-gap_offset;
					if(codE<codB)
					{
						codE=codB;
					}
					gap_offset=0;
					gap_closing_flag=0;
					//double check;
				}



				string cont_s;
				cont_s=contigs_info->contigs_str[cur_ctg].substr(codB,codE-codB);
				reverse(cont_s.begin(),cont_s.end());
				complement_str(cont_s);

				o_sc<<cont_s;
				if(nxt_ctg!=0)
				{
					if(contigs_info->contig_adjacency_left[cur_ctg].count(-nxt_ctg)>0)
					{
						string bridge=contigs_info->contig_adjacency_left[cur_ctg][-nxt_ctg].bridge;
						reverse(bridge.begin(),bridge.end());
						complement_str(bridge);


						
						o_sc<<bridge;
						int zero=0;
						o_sc_info<<"d: "<<zero<<" ";
						
					}
					else
					{
						if(contigs_info->scaffold_adjacency_left[cur_ctg].count(-nxt_ctg)>0)
						{
							//o_sc<<contigs_info->contigs_str[cur_ctg].substr(..,..);
							
							int gap_sz=contigs_info->scaffold_adjacency_left[cur_ctg][-nxt_ctg].dist_sum;
							int gap_cov=contigs_info->scaffold_adjacency_left[cur_ctg][-nxt_ctg].cov;
							gap_sz=gap_sz/gap_cov;

						//	cout<<gap_sz<<endl;
							if(gap_sz<150&&gap_sz>-150)
							{
						
								string tmp_str=contigs_info->contigs_str[cur_ctg].substr(0,codB);
								reverse(tmp_str.begin(),tmp_str.end());
								complement_str(tmp_str);
								o_sc<<tmp_str;//output all that can be extended
								//gap//
								string match1_s=contigs_info->contigs_str[cur_ctg].substr(0,200);
								reverse(match1_s.begin(),match1_s.end());
								complement_str(match1_s);
								int match1_sz=match1_s.size();
								int m_k=11;
								map<string,int> match_pos;	
								for(int i=0;i<=match1_sz - m_k;++i)
								{
									string substr=match1_s.substr(i,m_k);
									match_pos[substr]=i;
								}
								string nxt_ctg_s=contigs_info->contigs_str[abs(nxt_ctg)];
								if(nxt_ctg<0)
								{
									reverse(nxt_ctg_s.begin(),nxt_ctg_s.end());
									complement_str(nxt_ctg_s);
								}
								nxt_ctg_s=nxt_ctg_s.substr(0,200);
								int query_sz=nxt_ctg_s.size();
								map<int,int> gapoffset_cov;
								for(int i=0;i<query_sz-m_k;++i)
								{
									string substr=nxt_ctg_s.substr(i,m_k);
									if(match_pos.count(substr))
									{
										gap_closing_flag=1;
										gap_offset=i+match1_sz-match_pos[substr];
										if(gap_offset>=query_sz)
										{
											gap_closing_flag=0;
											gap_offset=0;
										}
										break;
								
									}
								}


								
								map<int,int>::iterator offset_it;
								int best_offset=0,max_cov=0;
								for(offset_it=gapoffset_cov.begin();offset_it!=gapoffset_cov.end();++offset_it)
								{
									if((offset_it->second)>max_cov)
									{
										max_cov=offset_it->second;
										best_offset=offset_it->first;
									}

									
								}
								gap_offset=best_offset;

								if(gap_closing_flag==0)
								{
									
									gap_offset=0;
								}


							}
							
							
							if(gap_closing_flag==0)
							{
								for(int gg=0;gg<gap_sz;++gg)
								{o_sc<<'N';}
								if(gap_sz<0)
								{o_sc<<'N';}
								
							//	contigs_info->gaps_in_scaffolds[sc_cnt].push_back(gap_sz);
								o_sc_info<<"d: "<<gap_sz<<" ";
								scf_len+=gap_sz;
							}
							else
							{
								//contigs_info->gaps_in_scaffolds[sc_cnt].push_back(0);
								int zero=0;
								o_sc_info<<"d: "<<zero<<" ";
							}
						}
						else
						{
						//gap

						}

					}


				}
			}

		}
		o_sc_info<<endl;
		o_sc<<endl;


	}



	for(int i=1;i<=tot_ctgs;++i )
	{
		if(contigs_info->c_info_vt[i].used==0)
		{
			sc_cnt++;
			o_sc_info<<">SuperContig_"<<sc_cnt<<endl;
			o_sc_info<<"1: ";
			o_sc_info<<i<<endl;
			o_sc<<">SuperContig_"<<sc_cnt<<endl;
			o_sc<<contigs_info->contigs_str[i]<<endl;

		}

	}



}


void ConstructContigGraph(struct hashtable *ht1,struct hashtable *merge_ht1, int K_size,contigs_info * contigs_info,string ContigFilename)
{

	time_t beg_time,read_time;
	string in_fname=ContigFilename;
	ifstream Contigs_in(in_fname.c_str());
	size_t num_Contigs=0;
	string tag,s,kmer,str,seq_s;
	uint64_t f_seq,hv;
	size_t hash_idx;
	bool found;
	bool flip_1,flip_2,flip_0;
	size_t ht_sz;
	size_t numReads=0;
	struct hashtable2 ht2,merge_ht2;

	int boundary=0,removed=0,bridge=0;

	ht_sz=ht1->ht_sz;

	cout<<"Contigs remapping."<<endl;
	if(ContigFilename=="Contigs.txt")
	{
		ContigsRemapping(ht1,&ht2, K_size, contigs_info,ContigFilename,0);
		
	}
	
	AppendMergeHT(ht1, merge_ht1);
	BuildContigAdjacency(ht1, &ht2, contigs_info, K_size,ContigFilename);	
		

}



void ConstructContigGraph0(struct hashtable0 *ht,struct hashtable0 *merge_ht, int K_size, contigs_info * contigs_info,string ContigFilename)
{

	time_t beg_time,read_time;
	string in_fname=ContigFilename;
	ifstream Contigs_in(in_fname.c_str());
	size_t num_Contigs=0;
	string tag,s,kmer,str,seq_s;
	uint64_t f_seq,hv;
	size_t hash_idx;
	bool found;
	bool flip_1,flip_2,flip_0;
	size_t ht_sz;
	size_t numReads=0;
	
	
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
	 AppendMergeHT0(ht, merge_ht,Kmer_arr_sz);
	//cout<<"Collecting informative reads."<<endl;
	BuildContigAdjacency0(ht,contigs_info, K_size,ContigFilename);
	

}



#endif
