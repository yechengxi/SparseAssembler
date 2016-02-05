#ifndef __GRAPH_CONSTRUCTION_H
#define __GRAPH_CONSTRUCTION_H


#include <iostream>
#include <string>
#include <string.h>
#include <stdint.h>
#include <vector>
#include <map>
#include <list>
#include <algorithm>
#include <fstream>
#include <sstream>
#include "time.h"
#include "BasicDataStructure.h"
using namespace std;


// initialize a hashtable
void Init_HT(struct hashtable* ht,size_t ht_sz)
{
	ht->ht_sz=ht_sz;
	ht->store_pos=(struct bucket**)calloc(ht_sz,sizeof(struct bucket*));

	for(size_t i=0;i<ht_sz;++i)
	{
		ht->store_pos[i]=NULL;
	}
}

void Init_HT2(struct hashtable2* ht,size_t ht_sz)
{
	ht->ht_sz=ht_sz;
	ht->store_pos=(struct bucket2**)calloc(ht_sz,sizeof(struct bucket2*));
	for(size_t i=0;i<ht_sz;++i)
	{
		ht->store_pos[i]=NULL;
	}
}

void Init_HT3(struct hashtable3* ht,size_t ht_sz)
{
	ht->ht_sz=ht_sz;
	ht->store_pos=(struct bucket3**)calloc(ht_sz,sizeof(struct bucket3*));
	for(size_t i=0;i<ht_sz;++i)
	{
		ht->store_pos[i]=NULL;
	}
}

void Init_HT4(struct hashtable4* ht,size_t ht_sz)
{
	ht->ht_sz=ht_sz;
	ht->store_pos=(struct bucket4**)calloc(ht_sz,sizeof(struct bucket4*));
	for(size_t i=0;i<ht_sz;++i)
	{
		ht->store_pos[i]=NULL;
	}
}

void Init_HT0(struct hashtable0* ht,size_t ht_sz)
{
	ht->ht_sz=ht_sz;
	ht->store_pos=(struct bucket0**)calloc(ht_sz,sizeof(struct bucket0*));
	for(size_t i=0;i<ht_sz;++i)
	{
		ht->store_pos[i]=NULL;
	}
}






//look up for a k-mer in a hashtable, if exists: 1, otherwise: 0. for round 2
bool look_up_in_a_list(uint64_t seq,struct bucket *** ptr)
{
	bool found=0;
	//struct bucket ** bktptr;


	while((**ptr)!=NULL)
	{
		if((**ptr)->kmer_t.kmer==seq)
		{
		//	bktptr=ptr;
			break;
		}

		(*ptr)=&((**ptr)->nxt_bucket);
	}
	if((**ptr)==NULL)
	{
		found=0;
		//bktptr=ptr;
	}
	else
	{
		found=1;
	}
	return found;
}

bool look_up_in_a_list2(struct kmer_t2 *seq,struct bucket2 *** ptr)
{
	bool found=0;
	//struct bucket ** bktptr;

	while((**ptr)!=NULL)
	{
		if(memcmp(&((**ptr)->kmer_t2.kmer),&(seq->kmer),sizeof(uint64_t)*2)==0)
		{
			break;
		}

		(*ptr)=&((**ptr)->nxt_bucket);
	}
	if((**ptr)==NULL)
	{
		found=0;
	}
	else
	{
		found=1;
	}
	return found;
}

bool look_up_in_a_list3(struct kmer_t3 *seq,struct bucket3 *** ptr)
{
	bool found=0;
	//struct bucket ** bktptr;

	while((**ptr)!=NULL)
	{
		if(memcmp(&((**ptr)->kmer_t3.kmer),&(seq->kmer),sizeof(uint64_t)*3)==0)
		{
			break;
		}

		(*ptr)=&((**ptr)->nxt_bucket);
	}
	if((**ptr)==NULL)
	{
		found=0;
	}
	else
	{
		found=1;
	}
	return found;
}

bool look_up_in_a_list4(struct kmer_t4 *seq,struct bucket4 *** ptr)
{
	bool found=0;
	//struct bucket ** bktptr;

	while((**ptr)!=NULL)
	{
		if(memcmp(&((**ptr)->kmer_t4.kmer),&(seq->kmer),sizeof(uint64_t)*4)==0)
		{
			break;
		}

		(*ptr)=&((**ptr)->nxt_bucket);
	}
	if((**ptr)==NULL)
	{
		found=0;
	}
	else
	{
		found=1;
	}
	return found;
}


bool look_up_in_a_list0( uint64_t *seq,struct bucket0 *** ptr,int arr_sz)
{
	bool found=0;
	//struct bucket ** bktptr;

	while((**ptr)!=NULL)
	{
		if(memcmp((**ptr)->kmer_t,seq,sizeof(uint64_t)*arr_sz)==0)
		{
			break;
		}

		(*ptr)=&((**ptr)->nxt_bucket);
	}
	if((**ptr)==NULL)
	{
		found=0;
	}
	else
	{
		found=1;
	}
	return found;
}



bool look_up_in_a_list_rm(uint64_t seq,struct bucket_rm *** ptr)
{
	bool found=0;
	//struct bucket ** bktptr;


	while((**ptr)!=NULL)
	{
		if((**ptr)->kmer_t.kmer==seq)
		{
		//	bktptr=ptr;
			break;
		}

		(*ptr)=&((**ptr)->nxt_bucket);
	}
	if((**ptr)==NULL)
	{
		found=0;
		//bktptr=ptr;
	}
	else
	{
		found=1;
	}
	return found;
}



bool look_up_in_a_list_rm2(struct kmer_t2 *seq,struct bucket_rm2 *** ptr)
{
	bool found=0;
	//struct bucket ** bktptr;

	while((**ptr)!=NULL)
	{
		if(memcmp(&((**ptr)->kmer_t2.kmer),&(seq->kmer),sizeof(uint64_t)*2)==0)
		{
			break;
		}

		(*ptr)=&((**ptr)->nxt_bucket);
	}
	if((**ptr)==NULL)
	{
		found=0;
	}
	else
	{
		found=1;
	}
	return found;
}

bool look_up_in_a_list_rm3(struct kmer_t3 *seq,struct bucket_rm3 *** ptr)
{
	bool found=0;
	//struct bucket ** bktptr;

	while((**ptr)!=NULL)
	{
		if(memcmp(&((**ptr)->kmer_t3.kmer),&(seq->kmer),sizeof(uint64_t)*3)==0)
		{
			break;
		}

		(*ptr)=&((**ptr)->nxt_bucket);
	}
	if((**ptr)==NULL)
	{
		found=0;
	}
	else
	{
		found=1;
	}
	return found;
}

bool look_up_in_a_list_rm4(struct kmer_t4 *seq,struct bucket_rm4 *** ptr)
{
	bool found=0;
	//struct bucket ** bktptr;

	while((**ptr)!=NULL)
	{
		if(memcmp(&((**ptr)->kmer_t4.kmer),&(seq->kmer),sizeof(uint64_t)*4)==0)
		{
			break;
		}

		(*ptr)=&((**ptr)->nxt_bucket);
	}
	if((**ptr)==NULL)
	{
		found=0;
	}
	else
	{
		found=1;
	}
	return found;
}


bool look_up_in_a_list_rm0(uint64_t *seq,struct bucket_rm0 *** ptr,int Kmer_arr_sz)
{
	bool found=0;
	//struct bucket ** bktptr;

	while((**ptr)!=NULL)
	{
		if(memcmp(&((**ptr)->kmer_t),seq,sizeof(uint64_t)*Kmer_arr_sz)==0)
		{
			break;
		}

		(*ptr)=&((**ptr)->nxt_bucket);
	}
	if((**ptr)==NULL)
	{
		found=0;
	}
	else
	{
		found=1;
	}
	return found;
}

// compare 2 arrays.
int uint64_t_cmp(uint64_t* A,uint64_t* B,int Kmer_arr_sz)
{
	int flag=0;
	for(int jj=0;jj<Kmer_arr_sz;++jj)
	{
		if(A[jj]>B[jj])
		{
			flag=1;
			break;
		}
		if(A[jj]<B[jj])
		{
			flag=-1;
			break;
		}
		if(A[jj]==B[jj])
		{
			continue;
		}
	}
	return flag;


}


//look up for a k-mer in a hashtable, if exists: 1, otherwise: 0. for round 1
bool look_up_in_a_list_r1(uint64_t seq,struct bucket_r1 *** ptr)
{	
	bool found=0;

	while((**ptr)!=NULL)
	{
		if((**ptr)->kmer_t.kmer==seq)
		{
		//	bktptr=ptr;	
			break;
		}

		(*ptr)=&((**ptr)->nxt_bucket);
	}
	if((**ptr)==NULL)
	{
		found=0;
		//bktptr=ptr;
	}
	else
	{
		found=1;
	}
	return found;
}

bool look_up_in_a_list2_r1(struct kmer_t2 *seq,struct bucket2_r1 *** ptr)
{	
	bool found=0;
	
	while((**ptr)!=NULL)
	{
		if(memcmp(&((**ptr)->kmer_t2.kmer),&(seq->kmer),sizeof(uint64_t)*2)==0)
		{
			break;
		}

		(*ptr)=&((**ptr)->nxt_bucket);
	}
	if((**ptr)==NULL)
	{
		found=0;
	}
	else
	{
		found=1;
	}
	return found;
}

bool look_up_in_a_list3_r1(struct kmer_t3 *seq,struct bucket3_r1 *** ptr)
{	
	bool found=0;
	
	while((**ptr)!=NULL)
	{
		if(memcmp(&((**ptr)->kmer_t3.kmer),&(seq->kmer),sizeof(uint64_t)*3)==0)
		{
			break;
		}

		(*ptr)=&((**ptr)->nxt_bucket);
	}
	if((**ptr)==NULL)
	{
		found=0;
	}
	else
	{
		found=1;
	}
	return found;
}

bool look_up_in_a_list4_r1(struct kmer_t4 *seq,struct bucket4_r1 *** ptr)
{	
	bool found=0;
	
	while((**ptr)!=NULL)
	{
		if(memcmp(&((**ptr)->kmer_t4.kmer),&(seq->kmer),sizeof(uint64_t)*4)==0)
		{
			break;
		}

		(*ptr)=&((**ptr)->nxt_bucket);
	}
	if((**ptr)==NULL)
	{
		found=0;
	}
	else
	{
		found=1;
	}
	return found;
}


bool look_up_in_a_list0_r1(uint64_t *seq,struct bucket0_r1 *** ptr,int arr_sz)
{	
	bool found=0;
	
	while((**ptr)!=NULL)
	{
		if(memcmp((**ptr)->kmer_t,seq,sizeof(uint64_t)*arr_sz)==0)
		{
			break;
		}

		(*ptr)=&((**ptr)->nxt_bucket);
	}
	if((**ptr)==NULL)
	{
		found=0;
	}
	else
	{
		found=1;
	}
	return found;
}



// sparse k-mer graph construction function, it contains 2 rounds. In the first round we select the sparse k-mers,
//in the second, we build the links between the k-mers.
void Sparse_Kmer_Graph_Construction(struct read_t *read,struct hashtable *ht,int64_t *bucket_count,int64_t *edge_cnt,int K_size,int gap,BF_info * BF_info ,int round)
{
	int readLen=read->readLen;
	int OverlappingKmers=readLen-K_size+1;
	int pre_search_len=100;
	
	if(OverlappingKmers<pre_search_len||round==2)
	{
		pre_search_len=OverlappingKmers;
	}
	
	
	if(gap>=OverlappingKmers)
	{	return;}



	int Read_arr_sz=readLen/32+1;
	int rem=readLen%32;
	if(rem==0)
	{Read_arr_sz--;}
	int tot_bits=Read_arr_sz*64;
	size_t ht_sz=ht->ht_sz;

	bool flip_beg,flip_end,found,found_bf;
	size_t hash_idx;
	uint64_t seq,f_seq,hv;
	int saved_kmers=0;
	
	//bucket ** bktptr;
	bucket ** bktptr_beg;
	bucket ** bktptr_end;
	
	int beg=0,end=0,end_bf=0;

	int first_found=-1,first_found_bf=-1;;



	//check the read to see if there is a saved kmer in the hashtable or bloom filter
	
	for (int j=0;j<pre_search_len;j++)
	{
		get_sub_arr(read->read_bits,read->readLen,j,K_size,&(seq));
	
		f_seq=get_rev_comp_seq(seq,K_size);
		if(seq>f_seq)
		{
			uint64_t t=seq;
			seq=f_seq;
			f_seq=t;
			//flip_beg=1;
		}
	
		hv=MurmurHash64A(&seq,sizeof(seq),0);

		hash_idx=(size_t) (hv%ht_sz);

		bktptr_end= &(ht->store_pos[hash_idx]);

		if(round==1)
		{
			found=look_up_in_a_list_r1(seq, ((bucket_r1 ***) &bktptr_end ));

			if(first_found_bf<0)
			{
				//check if there is a saved one in the BF table
				bool bf_inserted=0;
				if(BF_info->Bloom)
				{
				
					bf_inserted=1;
					if(BF_info->Bloom) 
					{
				
						for(int b=0;b<BF_info->d;++b)
						{
							size_t BF_hv=(size_t)(MurmurHash64A(&seq,sizeof(seq),b)%(BF_info->m));
							//size_t n_byte=BF_hv/8;
							//int n_bit=(int)(BF_hv%8);
							//int mask=1<<n_bit;
							//if((BF_info->BF_HT[n_byte]&mask)==0)
					
							if((BF_info->BF_HT[BF_hv >> 3] & (1U << (BF_hv & 0x7U)))==0)
							{
								bf_inserted=0;
							}
							//BF_info->BF_HT[BF_hv >> 3] |= (1U << (BF_hv & 0x7U));
						}
					}
					if(bf_inserted)
					{
						first_found_bf=j;
						break;
					}

				}
			}

		}
		else
		{
			found=look_up_in_a_list(seq,&bktptr_end);
		}


		if(found==1&&first_found<0)
		{
			first_found=j;
			break;
		}

	}


	

	///////updated for better sparseness
	if(round==1)
	{
		if(first_found>0)
		{beg=first_found%gap;}
		if(first_found<0&&first_found_bf>0)
		{beg=first_found_bf%gap;}
	}
	else
	{
		if(first_found<0)
		{return;}
		else
		{beg=first_found;}

	}



	saved_kmers=0;
	//look for the next sparse kmer

	while(beg<OverlappingKmers)
	{

		//cout<<saved_kmers<<endl;
		//check if there is a saved k-mer, if so resume from that, otherwise, save the last one.
		if(round==1)
		{

			for (int k=0;k<=gap;++k)
			{
				if(saved_kmers>0&&k==0)
				{	
					continue;
					
				}
				else
				{
					;
				}
				
				end=beg+k;
			
				if( end == OverlappingKmers )//overbound
				{
					return;
				}


				get_sub_arr(read->read_bits,read->readLen,end,K_size,&(seq));
				flip_end=0;
				f_seq=get_rev_comp_seq(seq,K_size);
				if(seq>f_seq)
				{
					uint64_t t=seq;
					seq=f_seq;
					f_seq=t;
					flip_end=1;
				}
	
				hv=MurmurHash64A(&seq,sizeof(seq),0);

				hash_idx=(size_t) (hv%ht_sz);

				bktptr_end= &(ht->store_pos[hash_idx]);

				if(round==1)
				{
					found=look_up_in_a_list_r1(seq, ((bucket_r1 ***) &bktptr_end ));

					//if found==0, check if something is in the bloom filter,here may be a tricky place, we shall treat bf and ht equally.
			
					if(found==0&&BF_info->Bloom)
					{
						bool bf_inserted=1;

						if(BF_info->Bloom) 
						{
				
							for(int b=0;b<BF_info->d;++b)
							{
								size_t BF_hv=(size_t)(MurmurHash64A(&seq,sizeof(seq),b)%(BF_info->m));
					
								//size_t n_byte=BF_hv/8;
								//int n_bit=(int)(BF_hv%8);
								//int mask=1<<n_bit;
								//if((BF_info->BF_HT[n_byte]&mask)==0)
						
								if((BF_info->BF_HT[BF_hv >> 3] & (1U << (BF_hv & 0x7U)))==0)
								{
									bf_inserted=0;
								}
								//BF_info->BF_HT[BF_hv >> 3] |= (1U << (BF_hv & 0x7U));
							}
						}
						if(bf_inserted)
						{
							//end_bf=beg+k;
							break;
						}
			
					}

				}
				else
				{
					found=look_up_in_a_list(seq,  &bktptr_end );
				
				}
			
			
				if(found==1)
				{

					break;
				}


				if(saved_kmers==0&&k==0)
				{	
					break;
					
				}
				


			}


			found_bf=1;
			if(found==0&&BF_info->Bloom)
			{
				for(int b=0;b<BF_info->d;++b)
				{
					size_t BF_hv=(size_t)(MurmurHash64A(&seq,sizeof(seq),b)%(BF_info->m));
					//size_t n_byte=BF_hv/8;
					//int n_bit=(int)(BF_hv%8);
					//int mask=1<<n_bit;
					//if((BF_info->BF_HT[n_byte]&mask)==0)
				
			
					if((BF_info->BF_HT[BF_hv >> 3] & (1U << (BF_hv & 0x7U)))==0)
					{
						found_bf=0;
						BF_info->BF_HT[BF_hv >> 3] |= (1U << (BF_hv & 0x7U));
						//BF_info->BF_HT[BF_hv >> 3] |= (1U << (BF_hv & 0x7U));
					}
					
				}
				//insertion finished, skip
				if(found_bf==0)
				{
					found_bf=1;
					beg=end;
					bktptr_beg=bktptr_end;
					saved_kmers++;
					continue;
				}
			}

		
		}
		else
		{

			for (int k=0;k<=OverlappingKmers;++k)
			{

				if(saved_kmers>0&&k==0)
				{	
					continue;
					
				}
				else
				{
					;
				}

				end=beg+k;
				
			
				if( end == OverlappingKmers )//overbound
				{
					return;
				}


				get_sub_arr(read->read_bits,read->readLen,end,K_size,&(seq));
	
				f_seq=get_rev_comp_seq(seq,K_size);
				flip_end=0;
				if(seq>f_seq)
				{
					uint64_t t=seq;
					seq=f_seq;
					f_seq=t;
					flip_end=1;
				}
	
				hv=MurmurHash64A(&seq,sizeof(seq),0);

				hash_idx=(size_t) (hv%ht_sz);

				bktptr_end= &(ht->store_pos[hash_idx]);

				found=look_up_in_a_list(seq,  &bktptr_end );

				if(found==1)
				{
					break;
				}
				if(saved_kmers==0)
				{break;}
			}


		}

		if(found==0)
		{
			if(round==1)
			{



				*bktptr_end=(struct bucket*)malloc(sizeof(struct bucket_r1));
				memset(*bktptr_end,0,sizeof(struct bucket_r1));
				//memset((ptr->kmer_t),0,sizeof(struct ptr->kmer_t));
				//ptr=ptr->nxt_bucket;
				((struct bucket_r1*) *(bktptr_end))->nxt_bucket=NULL;
				((struct bucket_r1*) *(bktptr_end))->kmer_t.kmer=seq;
				//(*(bktptr[j]))->kmer_info.cov1=0;
				((struct bucket_r1*) *(bktptr_end))->kmer_info.cov1++;
				//(*(bktptr[j]))->kmer_info.left=0;
				//(*(bktptr[j]))->kmer_info.right=0;
				//(*(bktptr[j]))->kmer_info.split_left=0;
				//(*(bktptr[j]))->kmer_info.split_right=0;
				//(*(bktptr[j]))->kmer_info.used=0;
				//(*(bktptr[j]))->kmer_info.cod=0;
				//(*(bktptr[j]))->kmer_info.contig_no=0;
				(*bucket_count)++;
				
			
				
				beg=end;
				saved_kmers++;
				bktptr_beg=bktptr_end;
				continue;
			}
			else
			{
				return;
			}
		

		}
		else
		{
			saved_kmers++;
			if(round==1)
			{
				if(((struct bucket_r1*) *(bktptr_end))->kmer_info.cov1<0xffff)
				{((struct bucket_r1*) *(bktptr_end))->kmer_info.cov1++;}
				beg=end;
				bktptr_beg=bktptr_end;
				continue;
			}
			else
			{
				if((*(bktptr_end))->kmer_info.cov1<0xffff)
				{(*(bktptr_end))->kmer_info.cov1++;}
			}

		}




		if(round==2)
		{
			
			if(end-beg<=gap&&saved_kmers>1)
			{
				
				uint64_t edge_bits,edge_bits_rc;
				size_t edge_len=end-beg;
				get_sub_arr(read->read_bits,read->readLen,beg,edge_len,&edge_bits);
				edge_bits_rc=get_rev_comp_seq(edge_bits,edge_len);

				
				if(flip_end==0)
				{
					struct edge_node **edge_node_p2p=&((*(bktptr_end))->kmer_info.left);
					while((*edge_node_p2p)!=NULL)
					{
						if((*edge_node_p2p)->edge==(uint64_t)edge_bits&&(((*edge_node_p2p)->len+1)==edge_len))
						{
							if((*edge_node_p2p)->edge_cov<0xf)
							{(*edge_node_p2p)->edge_cov++;}
							break;
						}
						edge_node_p2p=&((*edge_node_p2p)->nxt_edge);
					}
					if((*edge_node_p2p)==NULL)
					{
						(*edge_node_p2p) = (struct edge_node*)malloc(sizeof(struct edge_node));
						(*edge_cnt)++;
						memset(*edge_node_p2p,0,sizeof(struct edge_node));
						(*edge_node_p2p)->edge=(uint64_t)edge_bits;
						(*edge_node_p2p)->edge_cov=1;
						(*edge_node_p2p)->len=edge_len-1;
					}
				}
				else
				{
					//left_bits=get_rev_comp_seq(left_bits,h);
					struct edge_node **edge_node_p2p=&((*(bktptr_end))->kmer_info.right);
					while((*edge_node_p2p)!=NULL)
					{
						if((*edge_node_p2p)->edge==(uint64_t)edge_bits_rc&&(((*edge_node_p2p)->len+1)==edge_len))
						{
							if((*edge_node_p2p)->edge_cov<0xf)
							{(*edge_node_p2p)->edge_cov++;}
							break;
						}
						edge_node_p2p=&((*edge_node_p2p)->nxt_edge);

					}
					if((*edge_node_p2p)==NULL)
					{
						(*edge_node_p2p) = (struct edge_node*)malloc(sizeof(struct edge_node));
						(*edge_cnt)++;
						memset(*edge_node_p2p,0,sizeof(struct edge_node));
						(*edge_node_p2p)->edge=(uint64_t)edge_bits_rc;
						(*edge_node_p2p)->edge_cov=1;;
						(*edge_node_p2p)->len=edge_len-1;

					}
				}


				//uint64_t right_bits;
				//get_sub_arr(read->read_bits,read->readLen,j+K_size,h,&right_bits);

				get_sub_arr(read->read_bits,read->readLen,beg+K_size,edge_len,&edge_bits);
				edge_bits_rc=get_rev_comp_seq(edge_bits,edge_len);


				if(flip_beg==1)
				{

					//right_bits=get_rev_comp_seq(right_bits,h);

					struct edge_node **edge_node_p2p=&((*(bktptr_beg))->kmer_info.left);
					while((*edge_node_p2p)!=NULL)
					{
						if((*edge_node_p2p)->edge==(uint64_t)edge_bits_rc&&((*edge_node_p2p)->len+1)==edge_len)
						{
							if((*edge_node_p2p)->edge_cov<0xf)
							{(*edge_node_p2p)->edge_cov++;}
							break;
						}
						edge_node_p2p=&((*edge_node_p2p)->nxt_edge);
					}

					if((*edge_node_p2p)==NULL)
					{
						(*edge_node_p2p) = (struct edge_node*)malloc(sizeof(struct edge_node));
						(*edge_cnt)++;
						memset(*edge_node_p2p,0,sizeof(struct edge_node));
						(*edge_node_p2p)->edge=(uint64_t)edge_bits_rc;
						(*edge_node_p2p)->edge_cov=1;
						(*edge_node_p2p)->len=edge_len-1;
					}
				}
				else
				{

					struct edge_node **edge_node_p2p=&((*(bktptr_beg))->kmer_info.right);
					while((*edge_node_p2p)!=NULL)
					{
						if((*edge_node_p2p)->edge==(uint64_t)edge_bits&&((*edge_node_p2p)->len+1)==edge_len)
						{
							if((*edge_node_p2p)->edge_cov<0xf)
							{(*edge_node_p2p)->edge_cov++;}
							break;
						}
						edge_node_p2p=&((*edge_node_p2p)->nxt_edge);

					}
					if((*edge_node_p2p)==NULL)
					{
						(*edge_node_p2p) = (struct edge_node*)malloc(sizeof(struct edge_node));
						(*edge_cnt)++;
						memset(*edge_node_p2p,0,sizeof(struct edge_node));
						(*edge_node_p2p)->edge=(uint64_t)edge_bits;
						(*edge_node_p2p)->edge_cov=1;
						(*edge_node_p2p)->len=edge_len-1;
					}
				}


			}
			else
			{
				//cout<<"";
			}
		
			

		}
		beg=end;
		bktptr_beg=bktptr_end;
		flip_beg=flip_end;

	}


}




void Sparse_Kmer_Graph_Construction0(struct read_t *read,struct hashtable0 *ht, key_table *key_table, int64_t *bucket_count,int64_t *edge_cnt,int K_size,int gap,BF_info * BF_info ,int round)
{
	int readLen=read->readLen;
	int OverlappingKmers=readLen-K_size+1;
	int pre_search_len=100;
	
	if(OverlappingKmers<pre_search_len||round==2)
	{
		pre_search_len=OverlappingKmers;
	}
		
	if(gap>=OverlappingKmers)
	{	return;}

	int Read_arr_sz=readLen/32+1;
	int rem=readLen%32;
	if(rem==0)
	{Read_arr_sz--;}

	int Kmer_arr_sz=K_size/32+1;
	int rem1=K_size%32;
	if(rem1==0)
	{Kmer_arr_sz--;}

	int tot_bits=Read_arr_sz*64;
	size_t ht_sz=ht->ht_sz;

	bool flip_beg,flip_end,found,found_bf;
	size_t hash_idx;
	uint64_t seq[100],f_seq[100],temp_bits[100],hv;
	int saved_kmers=0;
	
	//bucket ** bktptr;
	bucket0 ** bktptr_beg;
	bucket0 ** bktptr_end;
	
	int beg=0,end=0,end_bf=0;

	int first_found=-1,first_found_bf=-1;;



	//check the read to see if there is a saved kmer in the hashtable or bloom filter
	
	for (int j=0;j<pre_search_len;j++)
	{
		get_sub_arr(read->read_bits,read->readLen,j,K_size,seq);
	
		memcpy(f_seq,seq,Kmer_arr_sz*sizeof(uint64_t));
		get_rev_comp_seq_arr(f_seq,K_size,Kmer_arr_sz);
		
		if(uint64_t_cmp(seq,f_seq,Kmer_arr_sz)>0)
		{
			memcpy(seq,f_seq,Kmer_arr_sz*sizeof(uint64_t));			
		}

		hv=MurmurHash64A(seq,sizeof(uint64_t)*Kmer_arr_sz,0);

		hash_idx=(size_t) (hv%ht_sz);

		bktptr_end= &(ht->store_pos[hash_idx]);

		if(round==1)
		{
			found=look_up_in_a_list0_r1(seq, (bucket0_r1 ***) &bktptr_end ,Kmer_arr_sz);

			if(first_found_bf<0)
			{
				//check if there is a saved one in the BF table
				bool bf_inserted=0;
				if(BF_info->Bloom)
				{
				
					bf_inserted=1;
					if(BF_info->Bloom) 
					{
				
						for(int b=0;b<BF_info->d;++b)
						{
							size_t BF_hv=(size_t)(MurmurHash64A(seq,sizeof(uint64_t)*Kmer_arr_sz,b)%(BF_info->m));
							//size_t n_byte=BF_hv/8;
							//int n_bit=(int)(BF_hv%8);
							//int mask=1<<n_bit;
							//if((BF_info->BF_HT[n_byte]&mask)==0)
					
							if((BF_info->BF_HT[BF_hv >> 3] & (1U << (BF_hv & 0x7U)))==0)
							{
								bf_inserted=0;
							}
							//BF_info->BF_HT[BF_hv >> 3] |= (1U << (BF_hv & 0x7U));
						}
					}
					if(bf_inserted)
					{
						first_found_bf=j;
						break;
					}

				}
			}

		}
		else
		{
			found=look_up_in_a_list0(seq, &bktptr_end ,Kmer_arr_sz);

		}


		if(found==1&&first_found<0)
		{
			first_found=j;
			break;
		}

	}


	

	///////updated for better sparseness
	if(round==1)
	{
		if(first_found>0)
		{beg=first_found%gap;}
		if(first_found<0&&first_found_bf>0)
		{beg=first_found_bf%gap;}
	}
	else
	{
		if(first_found<0)
		{return;}
		else
		{beg=first_found;}

	}



	saved_kmers=0;
	//look for the next sparse kmer

	while(beg<OverlappingKmers)
	{
		//cout<<saved_kmers<<endl;
		//check if there is a saved k-mer, if so resume from that, otherwise, save the last one.
		if(round==1)
		{

			for (int k=0;k<=gap;++k)
			{
				if(saved_kmers>0&&k==0)
				{	
					continue;
				}
				else
				{
					;
				}
				
				end=beg+k;
			
				if( end == OverlappingKmers )//overbound
				{
					return;
				}


				get_sub_arr(read->read_bits,read->readLen,end,K_size,seq);
				flip_end=0;
				
				memcpy(f_seq,seq,Kmer_arr_sz*sizeof(uint64_t));
				get_rev_comp_seq_arr(f_seq,K_size,Kmer_arr_sz);


				if(uint64_t_cmp(seq,f_seq,Kmer_arr_sz)>0)
				{
					memcpy(seq,f_seq,Kmer_arr_sz*sizeof(uint64_t));			
					flip_end=1;
				}

				hv=MurmurHash64A(seq,sizeof(uint64_t)*Kmer_arr_sz,0);

				hash_idx=(size_t) (hv%ht_sz);

				bktptr_end= &(ht->store_pos[hash_idx]);

				if(round==1)
				{

					
					found=look_up_in_a_list0_r1(seq, (bucket0_r1 ***) &bktptr_end ,Kmer_arr_sz);
					//if found==0, check if something is in the bloom filter,here may be a tricky place, we shall treat bf and ht equally.
			
					if(found==0&&BF_info->Bloom)
					{
						bool bf_inserted=1;

						if(BF_info->Bloom) 
						{
				
							for(int b=0;b<BF_info->d;++b)
							{
								size_t BF_hv=(size_t)(MurmurHash64A(seq,sizeof(uint64_t)*Kmer_arr_sz,b)%(BF_info->m));
					
								//size_t n_byte=BF_hv/8;
								//int n_bit=(int)(BF_hv%8);
								//int mask=1<<n_bit;
								//if((BF_info->BF_HT[n_byte]&mask)==0)
						
								if((BF_info->BF_HT[BF_hv >> 3] & (1U << (BF_hv & 0x7U)))==0)
								{
									bf_inserted=0;
								}
								//BF_info->BF_HT[BF_hv >> 3] |= (1U << (BF_hv & 0x7U));
							}
						}
						if(bf_inserted)
						{
							//end_bf=beg+k;
							break;
						}
			
					}

				}
				else
				{
					found=look_up_in_a_list0(seq,  &bktptr_end,Kmer_arr_sz);
				
				}
			
			
				if(found==1)
				{

					break;
				}


				if(saved_kmers==0&&k==0)
				{	
					break;
					
				}
				


			}


			found_bf=1;
			if(found==0&&BF_info->Bloom)
			{
				for(int b=0;b<BF_info->d;++b)
				{
					
					size_t BF_hv=(size_t)(MurmurHash64A(seq,sizeof(uint64_t)*Kmer_arr_sz,b)%(BF_info->m));
					//size_t n_byte=BF_hv/8;
					//int n_bit=(int)(BF_hv%8);
					//int mask=1<<n_bit;
					//if((BF_info->BF_HT[n_byte]&mask)==0)
				
			
					if((BF_info->BF_HT[BF_hv >> 3] & (1U << (BF_hv & 0x7U)))==0)
					{
						found_bf=0;
						BF_info->BF_HT[BF_hv >> 3] |= (1U << (BF_hv & 0x7U));
						//BF_info->BF_HT[BF_hv >> 3] |= (1U << (BF_hv & 0x7U));
					}
					
				}
				//insertion finished, skip
				if(found_bf==0)
				{
					found_bf=1;
					beg=end;
					bktptr_beg=bktptr_end;
					saved_kmers++;
					continue;
				}
			}

		
		}
		else
		{

			for (int k=0;k<=OverlappingKmers;++k)
			{

				if(saved_kmers>0&&k==0)
				{	
					continue;
					
				}
				else
				{
					;
				}

				end=beg+k;
				
			
				if( end == OverlappingKmers )//overbound
				{
					return;
				}


				get_sub_arr(read->read_bits,read->readLen,end,K_size,seq);
	
				memcpy(f_seq,seq,Kmer_arr_sz*sizeof(uint64_t));
				get_rev_comp_seq_arr(f_seq,K_size,Kmer_arr_sz);

				flip_end=0;
				if(uint64_t_cmp(seq,f_seq,Kmer_arr_sz)>0)
				{
					memcpy(seq,f_seq,Kmer_arr_sz*sizeof(uint64_t));			
					flip_end=1;
				}
	
				hv=MurmurHash64A(seq,sizeof(uint64_t)*Kmer_arr_sz,0);

				hash_idx=(size_t) (hv%ht_sz);

				bktptr_end= &(ht->store_pos[hash_idx]);

				found=look_up_in_a_list0(seq,  &bktptr_end,Kmer_arr_sz);

				if(found==1)
				{
					break;
				}
				if(saved_kmers==0)
				{break;}
			}


		}

		if(found==0)
		{
			if(round==1)
			{



				*bktptr_end=(struct bucket0*)malloc(sizeof(struct bucket_r1));
				memset(*bktptr_end,0,sizeof(struct bucket_r1));



				((struct bucket_r1*) *(bktptr_end))->nxt_bucket=NULL;

				//deal with the key_table
				if(key_table->current_index==key_table->KeysPerBlock)
				{
					uint64_t *new_block_ptr=(uint64_t *)malloc(sizeof(uint64_t)*Kmer_arr_sz*key_table->KeysPerBlock);	
					key_table->pblocks.push_back(new_block_ptr);
					key_table->current_index=0;
					key_table->current_block=key_table->pblocks.size();
						
				}

				
				uint64_t *block_ptr=(key_table->pblocks).back();
				uint64_t *key_ptr=&block_ptr[key_table->current_index*Kmer_arr_sz];
				memcpy(key_ptr,seq,sizeof(uint64_t)*Kmer_arr_sz);
				(key_table->current_index)++;
					

				((struct bucket0_r1*) *(bktptr_end))->kmer_t=key_ptr;
				((struct bucket0_r1*) *(bktptr_end))->kmer_info.cov1++;
				
				
				(*bucket_count)++;
				
			
				
				beg=end;
				saved_kmers++;
				bktptr_beg=bktptr_end;
				continue;
			}
			else
			{
				return;
			}
		

		}
		else
		{
			saved_kmers++;
			if(round==1)
			{
				if(((struct bucket0_r1*) *(bktptr_end))->kmer_info.cov1<0xffff)
				{((struct bucket0_r1*) *(bktptr_end))->kmer_info.cov1++;}
				beg=end;
				bktptr_beg=bktptr_end;
				continue;
			}
			else
			{
				if((*(bktptr_end))->kmer_info.cov1<0xffff)
				{(*(bktptr_end))->kmer_info.cov1++;}
			}

		}




		if(round==2)
		{
			
			if(end-beg<=gap&&saved_kmers>1)
			{
				
				uint64_t edge_bits,edge_bits_rc;
				size_t edge_len=end-beg;
				get_sub_arr(read->read_bits,read->readLen,beg,edge_len,&edge_bits);
				edge_bits_rc=get_rev_comp_seq(edge_bits,edge_len);

				
				if(flip_end==0)
				{
					struct edge_node **edge_node_p2p=&((*(bktptr_end))->kmer_info.left);
					while((*edge_node_p2p)!=NULL)
					{
						if((*edge_node_p2p)->edge==(uint64_t)edge_bits&&(((*edge_node_p2p)->len+1)==edge_len))
						{
							if((*edge_node_p2p)->edge_cov<0xf)
							{(*edge_node_p2p)->edge_cov++;}
							break;
						}
						edge_node_p2p=&((*edge_node_p2p)->nxt_edge);
					}
					if((*edge_node_p2p)==NULL)
					{
						(*edge_node_p2p) = (struct edge_node*)malloc(sizeof(struct edge_node));
						(*edge_cnt)++;
						memset(*edge_node_p2p,0,sizeof(struct edge_node));
						(*edge_node_p2p)->edge=(uint64_t)edge_bits;
						(*edge_node_p2p)->edge_cov=1;
						(*edge_node_p2p)->len=edge_len-1;
					}
				}
				else
				{
					//left_bits=get_rev_comp_seq(left_bits,h);
					struct edge_node **edge_node_p2p=&((*(bktptr_end))->kmer_info.right);
					while((*edge_node_p2p)!=NULL)
					{
						if((*edge_node_p2p)->edge==(uint64_t)edge_bits_rc&&(((*edge_node_p2p)->len+1)==edge_len))
						{
							if((*edge_node_p2p)->edge_cov<0xf)
							{(*edge_node_p2p)->edge_cov++;}
							break;
						}
						edge_node_p2p=&((*edge_node_p2p)->nxt_edge);

					}
					if((*edge_node_p2p)==NULL)
					{
						(*edge_node_p2p) = (struct edge_node*)malloc(sizeof(struct edge_node));
						(*edge_cnt)++;
						memset(*edge_node_p2p,0,sizeof(struct edge_node));
						(*edge_node_p2p)->edge=(uint64_t)edge_bits_rc;
						(*edge_node_p2p)->edge_cov=1;;
						(*edge_node_p2p)->len=edge_len-1;

					}
				}


				//uint64_t right_bits;
				//get_sub_arr(read->read_bits,read->readLen,j+K_size,h,&right_bits);

				get_sub_arr(read->read_bits,read->readLen,beg+K_size,edge_len,&edge_bits);
				edge_bits_rc=get_rev_comp_seq(edge_bits,edge_len);


				if(flip_beg==1)
				{

					//right_bits=get_rev_comp_seq(right_bits,h);

					struct edge_node **edge_node_p2p=&((*(bktptr_beg))->kmer_info.left);
					while((*edge_node_p2p)!=NULL)
					{
						if((*edge_node_p2p)->edge==(uint64_t)edge_bits_rc&&((*edge_node_p2p)->len+1)==edge_len)
						{
							if((*edge_node_p2p)->edge_cov<0xf)
							{(*edge_node_p2p)->edge_cov++;}
							break;
						}
						edge_node_p2p=&((*edge_node_p2p)->nxt_edge);
					}

					if((*edge_node_p2p)==NULL)
					{
						(*edge_node_p2p) = (struct edge_node*)malloc(sizeof(struct edge_node));
						(*edge_cnt)++;
						memset(*edge_node_p2p,0,sizeof(struct edge_node));
						(*edge_node_p2p)->edge=(uint64_t)edge_bits_rc;
						(*edge_node_p2p)->edge_cov=1;
						(*edge_node_p2p)->len=edge_len-1;
					}
				}
				else
				{

					struct edge_node **edge_node_p2p=&((*(bktptr_beg))->kmer_info.right);
					while((*edge_node_p2p)!=NULL)
					{
						if((*edge_node_p2p)->edge==(uint64_t)edge_bits&&((*edge_node_p2p)->len+1)==edge_len)
						{
							if((*edge_node_p2p)->edge_cov<0xf)
							{(*edge_node_p2p)->edge_cov++;}
							break;
						}
						edge_node_p2p=&((*edge_node_p2p)->nxt_edge);

					}
					if((*edge_node_p2p)==NULL)
					{
						(*edge_node_p2p) = (struct edge_node*)malloc(sizeof(struct edge_node));
						(*edge_cnt)++;
						memset(*edge_node_p2p,0,sizeof(struct edge_node));
						(*edge_node_p2p)->edge=(uint64_t)edge_bits;
						(*edge_node_p2p)->edge_cov=1;
						(*edge_node_p2p)->len=edge_len-1;
					}
				}


			}
			else
			{
				//cout<<"";
			}
		
			

		}
		beg=end;
		bktptr_beg=bktptr_end;
		flip_beg=flip_end;

	}


}





void Sparse_Kmer_Graph_Construction2(struct read_t *read,struct hashtable2 *ht,int64_t *bucket_count,int64_t *edge_cnt,int K_size,int gap,BF_info * BF_info,int round)
{
//	bool DISPLAY_SPLIT=0;

	int readLen=read->readLen;
	int OverlappingKmers=readLen-K_size+1;
	if(gap>=OverlappingKmers)
	{	return;}//}gap=OverlappingKmers-1;}

	int Read_arr_sz=readLen/32+1;
	int rem=readLen%32;
	if(rem==0)
	{Read_arr_sz--;}

	int Kmer_arr_sz=2;//K_size/32+1;
	int tot_bits=Read_arr_sz*64;
	size_t ht_sz=ht->ht_sz;
	bool flip[5000],found[5000];
	size_t hash_idx[5000];
	memset(flip,0,sizeof(flip));
	memset(found,0,sizeof(found));

	kmer_t2 seq[5000],f_seq[5000];
	uint64_t hv[5000],temp_bits[5000];

	bucket2 ** bktptr[5000];
	int first_found=-1;
	//char c_str[1000];
	//get_kmer_vec(read->read,K_size,(seq[0].kmer));
	//#pragma omp parallel for
	for (int j=0;j<OverlappingKmers;j++)
	{

//		cout<<j<<endl;

		get_sub_arr(read->read_bits,read->readLen,j,K_size,seq[j].kmer);



		memcpy(&f_seq[j],&seq[j],Kmer_arr_sz*sizeof(uint64_t));

		get_rev_comp_seq_arr((f_seq[j].kmer),K_size,Kmer_arr_sz);




		if(uint64_t_cmp(seq[j].kmer,f_seq[j].kmer,Kmer_arr_sz)>0)
		{
			flip[j]=1;
		}


		if(flip[j]==1)
		{
			memcpy(temp_bits,&(seq[j].kmer),Kmer_arr_sz*sizeof(uint64_t));
			memcpy(&(seq[j].kmer),&(f_seq[j].kmer),Kmer_arr_sz*sizeof(uint64_t));
			memcpy(&(f_seq[j].kmer),temp_bits,Kmer_arr_sz*sizeof(uint64_t));

		}



		hv[j]=MurmurHash64A((seq[j].kmer),sizeof(seq[j]),0);
	//	hv[j]=MurmurHash64A(&(seq[j].kmer[1]),sizeof(seq[j].kmer[1]),0);

		hash_idx[j]=(size_t) (hv[j]%ht_sz);

		bktptr[j]= &(ht->store_pos[hash_idx[j]]);

		if(round==1)
		{
			found[j]=look_up_in_a_list2_r1(&seq[j],(bucket2_r1***)&bktptr[j]);
			if(found[j]==1)
			{
				if(first_found<0)
				{first_found=j;}
			}
		}
		else
		{
			found[j]=look_up_in_a_list2(&seq[j],&bktptr[j]);
		}


	}



	int g,h;

	g=0;

	for (int k=0;k<gap;++k)
	{
		if(round==1)
		{
			found[k]=look_up_in_a_list2_r1(&seq[k],(struct bucket2_r1***) &bktptr[k]);
		}
		if(found[k]==1)
		{
			g=k;
			break;
		}

	}

	//check if there is a saved one in the BF table
	if(round==1&&found[g]==0&&BF_info->Bloom)
	{
		for (int k=0;k<gap;++k)
		{
			bool bf_inserted=1;
			if(BF_info->Bloom) 
			{
				
				for(int b=0;b<BF_info->d;++b)
				{
					size_t BF_hv=(size_t)(MurmurHash64A(&seq[k],sizeof(seq[k]),b)%(BF_info->m));
					//size_t n_byte=BF_hv/8;
					//int n_bit=(int)(BF_hv%8);
					//int mask=1<<n_bit;
					//if((BF_info->BF_HT[n_byte]&mask)==0)
					
					if((BF_info->BF_HT[BF_hv >> 3] & (1U << (BF_hv & 0x7U)))==0)
					{
						bf_inserted=0;
					}
					//BF_info->BF_HT[BF_hv >> 3] |= (1U << (BF_hv & 0x7U));
				}
			}
			if(bf_inserted)
			{
				g=k;break;
			}

		}
	}

	//check if there is a saved one in the BF table
	if(round==1&&found[g]==0&&BF_info->Bloom)
	{
		for (int k=0;k<gap;++k)
		{
			bool bf_inserted=1;
			if(BF_info->Bloom) 
			{
				
				for(int b=0;b<BF_info->d;++b)
				{
					size_t BF_hv=(size_t)(MurmurHash64A(&seq[k],sizeof(seq[k]),b)%(BF_info->m));
					//size_t n_byte=BF_hv/8;
					//int n_bit=(int)(BF_hv%8);
					//int mask=1<<n_bit;
					//if((BF_info->BF_HT[n_byte]&mask)==0)
					
					if((BF_info->BF_HT[BF_hv >> 3] & (1U << (BF_hv & 0x7U)))==0)
					{
						bf_inserted=0;
					}
					//BF_info->BF_HT[BF_hv >> 3] |= (1U << (BF_hv & 0x7U));
				}
			}
			if(bf_inserted)
			{
				g=k;break;
			}

		}
	}
	//look for the next sparse kmer

	//check if there is a saved one in the BF table
	if(round==1&&found[g]==0&&BF_info->Bloom)
	{
		for (int k=0;k<gap;++k)
		{
			bool bf_inserted=1;
			if(BF_info->Bloom) 
			{
				
				for(int b=0;b<BF_info->d;++b)
				{
					size_t BF_hv=(size_t)(MurmurHash64A(&seq[k],sizeof(seq[k]),b)%(BF_info->m));
					//size_t n_byte=BF_hv/8;
					//int n_bit=(int)(BF_hv%8);
					//int mask=1<<n_bit;
					//if((BF_info->BF_HT[n_byte]&mask)==0)
					
					if((BF_info->BF_HT[BF_hv >> 3] & (1U << (BF_hv & 0x7U)))==0)
					{
						bf_inserted=0;
					}
					//BF_info->BF_HT[BF_hv >> 3] |= (1U << (BF_hv & 0x7U));
				}
			}
			if(bf_inserted)
			{
				g=k;break;
			}

		}
	}
	if(round==1&&first_found>0&&found[g]==0)
	{g=first_found%gap;}
	//look for the next sparse kmer

	 for(int j=g;j<OverlappingKmers;)
	{

		h=gap;

		for (int k=1;k<=gap;++k)
		{
			if( (j+k)>=OverlappingKmers-1 &&(found[j+k]==0)&&k!=gap)
			{
				h=k+1;
				break;
			}
			if(round==1)
			{
				found[j+k]=look_up_in_a_list2_r1(&seq[j+k],(bucket2_r1***)&bktptr[j+k]);
			}

			if(k>0&&found[j+k]==1)
			{
				h=k;
				break;
			}

		}

		if(round==1&&found[j+h]==0&&BF_info->Bloom)
		{
			for (int k=1;k<=gap;++k)
			{

				bool bf_inserted=1;
				if(BF_info->Bloom) 
				{
				
					for(int b=0;b<BF_info->d;++b)
					{
						size_t BF_hv=(size_t)(MurmurHash64A(&seq[j+k],sizeof(seq[j+k]),b)%(BF_info->m));
						//size_t n_byte=BF_hv/8;
						//int n_bit=(int)(BF_hv%8);
						//int mask=1<<n_bit;
						//if((BF_info->BF_HT[n_byte]&mask)==0)
						
						if((BF_info->BF_HT[BF_hv >> 3] & (1U << (BF_hv & 0x7U)))==0)
						{
							bf_inserted=0;
						}
						//BF_info->BF_HT[BF_hv >> 3] |= (1U << (BF_hv & 0x7U));
					}
				}
				if(bf_inserted)
				{
					h=k;break;
				}
			}
		}


		if(round==1)
		{
			bool bf_inserted=1;
			if(BF_info->Bloom)
			{
				for(int b=0;b<BF_info->d;++b)
				{
					size_t BF_hv=(size_t)(MurmurHash64A(&seq[j],sizeof(seq[j]),b)%(BF_info->m));
					//size_t n_byte=BF_hv/8;
					//int n_bit=(int)(BF_hv%8);
					//int mask=1<<n_bit;
					//if((BF_info->BF_HT[n_byte]&mask)==0)
					
					if((BF_info->BF_HT[BF_hv >> 3] & (1U << (BF_hv & 0x7U)))==0)
					{
						bf_inserted=0;
					}
					BF_info->BF_HT[BF_hv >> 3] |= (1U << (BF_hv & 0x7U));
				}
				if(bf_inserted==0)
				{
					j=j+h;
					continue;
				}
			}
			//else, check if it has been inserted
			
			found[j]=look_up_in_a_list2_r1(&seq[j],(bucket2_r1***)&bktptr[j]);
		}

		if(found[j]==0)
		{
			if(round==1)
			{


				*(bktptr[j])=(struct bucket2*)malloc(sizeof(struct bucket2_r1));
				//memset((ptr->kmer_t),0,sizeof(struct ptr->kmer_t));
				//ptr=ptr->nxt_bucket;
				memset(*(bktptr[j]),0,sizeof(struct bucket2_r1));

				memcpy(&(((struct bucket2_r1*) *(bktptr[j]))->kmer_t2.kmer),&(seq[j].kmer),Kmer_arr_sz*sizeof(uint64_t));
	//			(*(bktptr[j]))->kmer_info.cov1=0;
				( (struct bucket2_r1*) *(bktptr[j]))->kmer_info.cov1++;
	//			(*(bktptr[j]))->kmer_info.left=0;
	//			(*(bktptr[j]))->kmer_info.right=0;
	//			(*(bktptr[j]))->kmer_info.split_left=0;
	//			(*(bktptr[j]))->kmer_info.split_right=0;
	//			(*(bktptr[j]))->kmer_info.used=0;
	//			(*(bktptr[j]))->kmer_info.cod=0;
	//			(*(bktptr[j]))->kmer_info.contig_no=0;
				(*bucket_count)++;
			}
		}
		else
		{
			if(round==1)
			{
				if(((struct bucket2_r1*) *(bktptr[j]))->kmer_info.cov1<0xffff)
				{
					((struct bucket2_r1*) *(bktptr[j]))->kmer_info.cov1++;
				}
			}
			else
			{
				if((*(bktptr[j]))->kmer_info.cov1<0xffff)
				{
					(*(bktptr[j]))->kmer_info.cov1++;
				}
			
			}

		}





///////////////////////////////////////////////////////////////////////////////////////////////////////////

		if(round==2)
		{
			if(h>0&&found[j]&&found[j+h])
			{
				uint64_t left_bits;
				get_sub_arr(read->read_bits,read->readLen,j,h,&left_bits);
				if(flip[j+h]==0)
				{
					struct edge_node **edge_node_p2p=&((*(bktptr[j+h]))->kmer_info.left);
					while((*edge_node_p2p)!=NULL)
					{
						if((*edge_node_p2p)->edge==(uint64_t)left_bits&&((*edge_node_p2p)->len+1)==h)
						{
							if((*edge_node_p2p)->edge_cov<0xf)
							{(*edge_node_p2p)->edge_cov++;}
							break;
						}
						edge_node_p2p=&((*edge_node_p2p)->nxt_edge);
					}
					if((*edge_node_p2p)==NULL)
					{
						(*edge_node_p2p) = (struct edge_node*)malloc(sizeof(struct edge_node));
						(*edge_cnt)++;
						memset(*edge_node_p2p,0,sizeof(struct edge_node));
						(*edge_node_p2p)->edge=(uint64_t)left_bits;
						(*edge_node_p2p)->edge_cov=1;
						(*edge_node_p2p)->len=h-1;
					}
				}
				else
				{
					left_bits=get_rev_comp_seq(left_bits,h);
					struct edge_node **edge_node_p2p=&((*(bktptr[j+h]))->kmer_info.right);
					while((*edge_node_p2p)!=NULL)
					{
						if((*edge_node_p2p)->edge==(uint64_t)left_bits&&((*edge_node_p2p)->len+1)==h)
						{
							if((*edge_node_p2p)->edge_cov<0xf)
							{(*edge_node_p2p)->edge_cov++;}
							break;
						}
						edge_node_p2p=&((*edge_node_p2p)->nxt_edge);

					}
					if((*edge_node_p2p)==NULL)
					{
						(*edge_node_p2p) = (struct edge_node*)malloc(sizeof(struct edge_node));
						(*edge_cnt)++;
						memset(*edge_node_p2p,0,sizeof(struct edge_node));
						(*edge_node_p2p)->edge=(uint64_t)left_bits;
						(*edge_node_p2p)->edge_cov=1;
						(*edge_node_p2p)->len=h-1;
					}
				}


			}


			if(h>0&&j+h<OverlappingKmers&&found[j]&&found[j+h])
			{
				uint64_t right_bits;
				get_sub_arr(read->read_bits,read->readLen,j+K_size,h,&right_bits);
				if(flip[j]==1)
				{

					right_bits=get_rev_comp_seq(right_bits,h);





					struct edge_node **edge_node_p2p=&((*(bktptr[j]))->kmer_info.left);
					while((*edge_node_p2p)!=NULL)
					{
						if((*edge_node_p2p)->edge==(uint64_t)right_bits&&((*edge_node_p2p)->len+1)==h)
						{
							if((*edge_node_p2p)->edge_cov<0xf)
							{(*edge_node_p2p)->edge_cov++;}
							break;
						}
						edge_node_p2p=&((*edge_node_p2p)->nxt_edge);
					}

					if((*edge_node_p2p)==NULL)
					{
						(*edge_node_p2p) = (struct edge_node*)malloc(sizeof(struct edge_node));
						(*edge_cnt)++;
						memset(*edge_node_p2p,0,sizeof(struct edge_node));
						(*edge_node_p2p)->edge=(uint64_t)right_bits;
						(*edge_node_p2p)->edge_cov=1;
						(*edge_node_p2p)->len=h-1;
					}
				}
				else
				{

					struct edge_node **edge_node_p2p=&((*(bktptr[j]))->kmer_info.right);
					while((*edge_node_p2p)!=NULL)
					{
						if((*edge_node_p2p)->edge==(uint64_t)right_bits&&((*edge_node_p2p)->len+1==h))
						{
							if((*edge_node_p2p)->edge_cov<0xf)
							{(*edge_node_p2p)->edge_cov++;}
							break;
						}
						edge_node_p2p=&((*edge_node_p2p)->nxt_edge);

					}
					if((*edge_node_p2p)==NULL)
					{
						(*edge_node_p2p) = (struct edge_node*)malloc(sizeof(struct edge_node));
						(*edge_cnt)++;
						memset(*edge_node_p2p,0,sizeof(struct edge_node));
						(*edge_node_p2p)->edge=(uint64_t)right_bits;
						(*edge_node_p2p)->edge_cov=1;
						(*edge_node_p2p)->len=h-1;
					}
				}


			}
		}


		j=j+h;
		//if(j>=OverlappingKmers-1)
		//{break;}
	}

}

void Sparse_Kmer_Graph_Construction3(struct read_t *read,struct hashtable3 *ht,int64_t *bucket_count,int64_t *edge_cnt,int K_size,int gap,BF_info * BF_info, int round)
{
//	bool DISPLAY_SPLIT=0;

	int readLen=read->readLen;
	int OverlappingKmers=readLen-K_size+1;
	if(gap>=OverlappingKmers)
	{	return;}//}gap=OverlappingKmers-1;}

	int Read_arr_sz=readLen/32+1;
	int rem=readLen%32;
	if(rem==0)
	{Read_arr_sz--;}

	int Kmer_arr_sz=3;

	int tot_bits=Read_arr_sz*64;
	size_t ht_sz=ht->ht_sz;
	bool flip[1000],found[1000];
	size_t hash_idx[1000];
	memset(flip,0,sizeof(flip));
	memset(found,0,sizeof(found));

	kmer_t3 seq[1000],f_seq[1000];
	memset(seq,0,sizeof(seq));
	memset(f_seq,0,sizeof(f_seq));

	uint64_t hv[1000],temp_bits[1000];

	bucket3 ** bktptr[1000];
	int first_found=-1;
	//char c_str[1000];

	//#pragma omp parallel for
	//get_kmer_vec(read->read,K_size,seq[0].kmer);
	for (int j=0;j<OverlappingKmers;j++)
	{

//		cout<<j<<endl;

		get_sub_arr(read->read_bits,read->readLen,j,K_size,seq[j].kmer);

		memcpy(&f_seq[j],&seq[j],Kmer_arr_sz*sizeof(uint64_t));

		get_rev_comp_seq_arr((f_seq[j].kmer),K_size,Kmer_arr_sz);


		if(uint64_t_cmp(seq[j].kmer,f_seq[j].kmer,Kmer_arr_sz)>0)
		{
			flip[j]=1;
		}


		if(flip[j]==1)
		{
			memcpy(temp_bits,&(seq[j].kmer),Kmer_arr_sz*sizeof(uint64_t));
			memcpy(&(seq[j].kmer),&(f_seq[j].kmer),Kmer_arr_sz*sizeof(uint64_t));
			memcpy(&(f_seq[j].kmer),temp_bits,Kmer_arr_sz*sizeof(uint64_t));

		}

		hv[j]=MurmurHash64A((seq[j].kmer),sizeof(seq[j]),0);
	//	hv[j]=MurmurHash64A(&(seq[j].kmer[1]),sizeof(seq[j].kmer[1]),0);

		hash_idx[j]=(size_t) (hv[j]%ht_sz);

		bktptr[j]= &(ht->store_pos[hash_idx[j]]);

		if(round==1)
		{

			found[j]=look_up_in_a_list3_r1(&seq[j],(bucket3_r1***)&bktptr[j]);
			if(found[j]==1)
			{
				if(first_found<0)
				{first_found=j;}
			}
		}
		else
		{
			found[j]=look_up_in_a_list3(&seq[j],&bktptr[j]);
		}


	}



	int g,h;

	g=0;

	for (int k=0;k<gap;++k)
	{
		if(round==1)
		{
			found[k]=look_up_in_a_list3_r1(&seq[k],(struct bucket3_r1***) &bktptr[k]);
		}
		if(found[k]==1)
		{
			g=k;
			break;
		}

	}


	if(round==1&&first_found>0&&found[g]==0)
	{g=first_found%gap;}
	 for(int j=g;j<OverlappingKmers;)
	{

		g=gap;

		for (int k=1;k<=gap;++k)
		{
			if( (j-k)<0 )
			{
				g=0;

				break;
			}
			//if(round==1)
			//{
			//	found[j-k]=look_up_in_a_list3_r1(&seq[j-k],(bucket3_r1***) &bktptr[j-k]);
			//}
			if(k>0&&(k<=j)&&found[j-k]==1)
			{
				g=k;
				break;
			}

		}

		h=gap;

		for (int k=1;k<=gap;++k)
		{
			if( (j+k)>=OverlappingKmers-1 &&(found[j+k]==0)&&k!=gap)
			{
				h=k+1;
				break;
			}
			if(round==1)
			{
				found[j+k]=look_up_in_a_list3_r1(&seq[j+k],(bucket3_r1***)&bktptr[j+k]);
			}

			if(k>0&&found[j+k]==1)
			{
				h=k;
				break;
			}

		}

		if(round==1&&found[j+h]==0&&BF_info->Bloom)
		{
			for (int k=1;k<=gap;++k)
			{

				bool bf_inserted=1;
				if(BF_info->Bloom) 
				{
				
					for(int b=0;b<BF_info->d;++b)
					{
						size_t BF_hv=(size_t)(MurmurHash64A(&seq[j+k],sizeof(seq[j+k]),b)%(BF_info->m));
						//size_t n_byte=BF_hv/8;
						//int n_bit=(int)(BF_hv%8);
						//int mask=1<<n_bit;
						//if((BF_info->BF_HT[n_byte]&mask)==0)
						
						if((BF_info->BF_HT[BF_hv >> 3] & (1U << (BF_hv & 0x7U)))==0)
						{
							bf_inserted=0;
						}
						//BF_info->BF_HT[BF_hv >> 3] |= (1U << (BF_hv & 0x7U));
					}
				}
				if(bf_inserted)
				{
					h=k;break;
				}
			}
		}


		if(round==1)
		{
			bool bf_inserted=1;
			if(BF_info->Bloom)
			{
				for(int b=0;b<BF_info->d;++b)
				{
					size_t BF_hv=(size_t)(MurmurHash64A(&seq[j],sizeof(seq[j]),b)%(BF_info->m));
					//size_t n_byte=BF_hv/8;
					//int n_bit=(int)(BF_hv%8);
					//int mask=1<<n_bit;
					//if((BF_info->BF_HT[n_byte]&mask)==0)
					
					if((BF_info->BF_HT[BF_hv >> 3] & (1U << (BF_hv & 0x7U)))==0)
					{
						bf_inserted=0;
					}
					BF_info->BF_HT[BF_hv >> 3] |= (1U << (BF_hv & 0x7U));
				}
				if(bf_inserted==0)
				{
					j=j+h;
					continue;
				}
			}
			//else, check if it has been inserted
			
			found[j]=look_up_in_a_list3_r1(&seq[j],(bucket3_r1***)&bktptr[j]);
		}

		if(found[j]==0)
		{
			if(round==1)
			{


				*(bktptr[j])=(struct bucket3*)malloc(sizeof(struct bucket3_r1));
				//memset((ptr->kmer_t),0,sizeof(struct ptr->kmer_t));
				//ptr=ptr->nxt_bucket;
				memset(*(bktptr[j]),0,sizeof(struct bucket3_r1));

				memcpy(&(((struct bucket3_r1*) *(bktptr[j]))->kmer_t3.kmer),&(seq[j].kmer),Kmer_arr_sz*sizeof(uint64_t));
	//			(*(bktptr[j]))->kmer_info.cov1=0;
				( (struct bucket3_r1*) *(bktptr[j]))->kmer_info.cov1++;
	//			(*(bktptr[j]))->kmer_info.left=0;
	//			(*(bktptr[j]))->kmer_info.right=0;
	//			(*(bktptr[j]))->kmer_info.split_left=0;
	//			(*(bktptr[j]))->kmer_info.split_right=0;
	//			(*(bktptr[j]))->kmer_info.used=0;
	//			(*(bktptr[j]))->kmer_info.cod=0;
	//			(*(bktptr[j]))->kmer_info.contig_no=0;
				(*bucket_count)++;
			}
		}
		else
		{
			if(round==1)
			{
				if(((struct bucket3_r1*) *(bktptr[j]))->kmer_info.cov1<0xffff)
				{
					((struct bucket3_r1*) *(bktptr[j]))->kmer_info.cov1++;
				}
			}
			else
			{
				if((*(bktptr[j]))->kmer_info.cov1<0xffff)
				{
					(*(bktptr[j]))->kmer_info.cov1++;
				}
			
			}

		}





///////////////////////////////////////////////////////////////////////////////////////////////////////////

		if(round==2)
		{
			if(h>0&&found[j]&&found[j+h])
			{

				uint64_t left_bits;
				get_sub_arr(read->read_bits,read->readLen,j,h,&left_bits);
				if(flip[j+h]==0)
				{

					struct edge_node **edge_node_p2p=&((*(bktptr[j+h]))->kmer_info.left);
					while((*edge_node_p2p)!=NULL)
					{
						if((*edge_node_p2p)->edge==(uint64_t)left_bits&&((*edge_node_p2p)->len+1)==h)
						{
							if((*edge_node_p2p)->edge_cov<0xf)
							{(*edge_node_p2p)->edge_cov++;}
							break;
						}
						edge_node_p2p=&((*edge_node_p2p)->nxt_edge);
					}
					if((*edge_node_p2p)==NULL)
					{
						(*edge_node_p2p) = (struct edge_node*)malloc(sizeof(struct edge_node));
						(*edge_cnt)++;
						memset(*edge_node_p2p,0,sizeof(struct edge_node));
						(*edge_node_p2p)->edge=(uint64_t)left_bits;
						(*edge_node_p2p)->edge_cov=1;
						(*edge_node_p2p)->len=h-1;
					}
				}
				else
				{
					left_bits=get_rev_comp_seq(left_bits,h);
					struct edge_node **edge_node_p2p=&((*(bktptr[j+h]))->kmer_info.right);
					while((*edge_node_p2p)!=NULL)
					{
						if((*edge_node_p2p)->edge==(uint64_t)left_bits&&((*edge_node_p2p)->len+1)==h)
						{
							if((*edge_node_p2p)->edge_cov<0xf)
							{(*edge_node_p2p)->edge_cov++;}
							break;
						}
						edge_node_p2p=&((*edge_node_p2p)->nxt_edge);

					}
					if((*edge_node_p2p)==NULL)
					{
						(*edge_node_p2p) = (struct edge_node*)malloc(sizeof(struct edge_node));
						(*edge_cnt)++;
						memset(*edge_node_p2p,0,sizeof(struct edge_node));
						(*edge_node_p2p)->edge=(uint64_t)left_bits;
						(*edge_node_p2p)->edge_cov=1;
						(*edge_node_p2p)->len=h-1;
					}
				}


			}


			if(h>0&&j+h<OverlappingKmers&&found[j]&&found[j+h])
			{
				uint64_t right_bits;
				get_sub_arr(read->read_bits,read->readLen,j+K_size,h,&right_bits);
				if(flip[j]==1)
				{

					right_bits=get_rev_comp_seq(right_bits,h);

					struct edge_node **edge_node_p2p=&((*(bktptr[j]))->kmer_info.left);
					while((*edge_node_p2p)!=NULL)
					{
						if((*edge_node_p2p)->edge==(uint64_t)right_bits&&((*edge_node_p2p)->len+1)==h)
						{
							if((*edge_node_p2p)->edge_cov<0xf)
							{(*edge_node_p2p)->edge_cov++;}
							break;
						}
						edge_node_p2p=&((*edge_node_p2p)->nxt_edge);
					}

					if((*edge_node_p2p)==NULL)
					{
						(*edge_node_p2p) = (struct edge_node*)malloc(sizeof(struct edge_node));
						(*edge_cnt)++;
						memset(*edge_node_p2p,0,sizeof(struct edge_node));
						(*edge_node_p2p)->edge=(uint64_t)right_bits;
						(*edge_node_p2p)->edge_cov=1;
						(*edge_node_p2p)->len=h-1;
					}
				}
				else
				{

					struct edge_node **edge_node_p2p=&((*(bktptr[j]))->kmer_info.right);
					while((*edge_node_p2p)!=NULL)
					{
						if((*edge_node_p2p)->edge==(uint64_t)right_bits&&((*edge_node_p2p)->len+1==h))
						{
							if((*edge_node_p2p)->edge_cov<0xf)
							{(*edge_node_p2p)->edge_cov++;}
							break;
						}
						edge_node_p2p=&((*edge_node_p2p)->nxt_edge);

					}
					if((*edge_node_p2p)==NULL)
					{
						(*edge_node_p2p) = (struct edge_node*)malloc(sizeof(struct edge_node));
						(*edge_cnt)++;
						memset(*edge_node_p2p,0,sizeof(struct edge_node));
						(*edge_node_p2p)->edge=(uint64_t)right_bits;
						(*edge_node_p2p)->edge_cov=1;
						(*edge_node_p2p)->len=h-1;
					}
				}


			}
		}


		j=j+h;
		//if(j>=OverlappingKmers-1)
		//{break;}
	}

}

void Sparse_Kmer_Graph_Construction4(struct read_t *read,struct hashtable4 *ht,int64_t *bucket_count,int64_t *edge_cnt,int K_size,int gap,BF_info * BF_info,int round)
{
//	bool DISPLAY_SPLIT=0;

	int readLen=read->readLen;
	int OverlappingKmers=readLen-K_size+1;
	if(gap>=OverlappingKmers)
	{	return;}//}gap=OverlappingKmers-1;}

	int Read_arr_sz=readLen/32+1;
	int rem=readLen%32;
	if(rem==0)
	{Read_arr_sz--;}

	int Kmer_arr_sz=4;

	int tot_bits=Read_arr_sz*64;
	size_t ht_sz=ht->ht_sz;
	bool flip[1000],found[1000];
	size_t hash_idx[1000];
	memset(flip,0,sizeof(flip));
	memset(found,0,sizeof(found));

	kmer_t4 seq[1000],f_seq[1000];
	memset(seq,0,sizeof(seq));
	memset(f_seq,0,sizeof(f_seq));

	uint64_t hv[1000],temp_bits[1000];

	bucket4 ** bktptr[1000];
	int first_found=-1;
	//char c_str[1000];
	//get_kmer_vec(read->read,K_size,seq[0].kmer);
	//#pragma omp parallel for
	for (int j=0;j<OverlappingKmers;j++)
	{

//		cout<<j<<endl;

		get_sub_arr(read->read_bits,read->readLen,j,K_size,seq[j].kmer);

		memcpy(&f_seq[j],&seq[j],Kmer_arr_sz*sizeof(uint64_t));

		get_rev_comp_seq_arr((f_seq[j].kmer),K_size,Kmer_arr_sz);


		if(uint64_t_cmp(seq[j].kmer,f_seq[j].kmer,Kmer_arr_sz)>0)
		{
			flip[j]=1;
		}


		if(flip[j]==1)
		{
			memcpy(temp_bits,&(seq[j].kmer),Kmer_arr_sz*sizeof(uint64_t));
			memcpy(&(seq[j].kmer),&(f_seq[j].kmer),Kmer_arr_sz*sizeof(uint64_t));
			memcpy(&(f_seq[j].kmer),temp_bits,Kmer_arr_sz*sizeof(uint64_t));

		}

		hv[j]=MurmurHash64A((seq[j].kmer),sizeof(seq[j]),0);
	//	hv[j]=MurmurHash64A(&(seq[j].kmer[1]),sizeof(seq[j].kmer[1]),0);

		hash_idx[j]=(size_t) (hv[j]%ht_sz);

		bktptr[j]= &(ht->store_pos[hash_idx[j]]);

		if(round==1)
		{
			found[j]=look_up_in_a_list4_r1(&seq[j],(bucket4_r1***)&bktptr[j]);
			if(found[j]==1)
			{
				if(first_found<0)
				{first_found=j;}
			}
		}
		else
		{
			found[j]=look_up_in_a_list4(&seq[j],&bktptr[j]);
		}


	}



	int g,h;

	g=0;

	for (int k=0;k<gap;++k)
	{
		if(round==1)
		{
			found[k]=look_up_in_a_list4_r1(&seq[k],(struct bucket4_r1***) &bktptr[k]);
		}
		if(found[k]==1)
		{
			g=k;
			break;
		}

	}

	if(round==1&&first_found>0&&found[g]==0)
	{g=first_found%gap;}
	 for(int j=g;j<OverlappingKmers;)
	{

		g=gap;

		for (int k=1;k<=gap;++k)
		{
			if( (j-k)<0 )
			{
				g=0;

				break;
			}
			//if(round==1)
			//{
			//	found[j-k]=look_up_in_a_list4_r1(&seq[j-k],(bucket4_r1***) &bktptr[j-k]);
			//}
			if(k>0&&(k<=j)&&found[j-k]==1)
			{
				g=k;
				break;
			}

		}

		h=gap;

		for (int k=1;k<=gap;++k)
		{
			if( (j+k)>=OverlappingKmers-1 &&(found[j+k]==0)&&k!=gap)
			{
				h=k+1;
				break;
			}
			if(round==1)
			{
				found[j+k]=look_up_in_a_list4_r1(&seq[j+k],(bucket4_r1***)&bktptr[j+k]);
			}

			if(k>0&&found[j+k]==1)
			{
				h=k;
				break;
			}

		}

		if(round==1&&found[j+h]==0&&BF_info->Bloom)
		{
			for (int k=1;k<=gap;++k)
			{

				bool bf_inserted=1;
				if(BF_info->Bloom) 
				{
				
					for(int b=0;b<BF_info->d;++b)
					{
						size_t BF_hv=(size_t)(MurmurHash64A(&seq[j+k],sizeof(seq[j+k]),b)%(BF_info->m));
						//size_t n_byte=BF_hv/8;
						//int n_bit=(int)(BF_hv%8);
						//int mask=1<<n_bit;
						//if((BF_info->BF_HT[n_byte]&mask)==0)
						
						if((BF_info->BF_HT[BF_hv >> 3] & (1U << (BF_hv & 0x7U)))==0)
						{
							bf_inserted=0;
						}
						//BF_info->BF_HT[BF_hv >> 3] |= (1U << (BF_hv & 0x7U));
					}
				}
				if(bf_inserted)
				{
					h=k;break;
				}
			}
		}

		if(round==1)
		{
			bool bf_inserted=1;
			if(BF_info->Bloom)
			{
				for(int b=0;b<BF_info->d;++b)
				{
					size_t BF_hv=(size_t)(MurmurHash64A(&seq[j],sizeof(seq[j]),b)%(BF_info->m));
					//size_t n_byte=BF_hv/8;
					//int n_bit=(int)(BF_hv%8);
					//int mask=1<<n_bit;
					//if((BF_info->BF_HT[n_byte]&mask)==0)
					
					if((BF_info->BF_HT[BF_hv >> 3] & (1U << (BF_hv & 0x7U)))==0)
					{
						bf_inserted=0;
					}
					BF_info->BF_HT[BF_hv >> 3] |= (1U << (BF_hv & 0x7U));
				}
				if(bf_inserted==0)
				{
					j=j+h;
					continue;
				}
			}
			//else, check if it has been inserted
			
			found[j]=look_up_in_a_list4_r1(&seq[j],(bucket4_r1***)&bktptr[j]);
		}

		if(found[j]==0)
		{
			if(round==1)
			{


				*(bktptr[j])=(struct bucket4*)malloc(sizeof(struct bucket4_r1));
				//memset((ptr->kmer_t),0,sizeof(struct ptr->kmer_t));
				//ptr=ptr->nxt_bucket;
				memset(*(bktptr[j]),0,sizeof(struct bucket4_r1));

				memcpy(&(((struct bucket4_r1*) *(bktptr[j]))->kmer_t4.kmer),&(seq[j].kmer),Kmer_arr_sz*sizeof(uint64_t));
	//			(*(bktptr[j]))->kmer_info.cov1=0;
				( (struct bucket4_r1*) *(bktptr[j]))->kmer_info.cov1++;
	//			(*(bktptr[j]))->kmer_info.left=0;
	//			(*(bktptr[j]))->kmer_info.right=0;
	//			(*(bktptr[j]))->kmer_info.split_left=0;
	//			(*(bktptr[j]))->kmer_info.split_right=0;
	//			(*(bktptr[j]))->kmer_info.used=0;
	//			(*(bktptr[j]))->kmer_info.cod=0;
	//			(*(bktptr[j]))->kmer_info.contig_no=0;
				(*bucket_count)++;
			}
		}
		else
		{
			if(round==1)
			{
				if(((struct bucket4_r1*) *(bktptr[j]))->kmer_info.cov1<0xffff)
				{
					((struct bucket4_r1*) *(bktptr[j]))->kmer_info.cov1++;
				}
			}
			else
			{
				if((*(bktptr[j]))->kmer_info.cov1<0xffff)
				{
					(*(bktptr[j]))->kmer_info.cov1++;
				}
			
			}

		}





///////////////////////////////////////////////////////////////////////////////////////////////////////////

		if(round==2)
		{
			if(h>0&&found[j]&&found[j+h])
			{

				uint64_t left_bits;
				get_sub_arr(read->read_bits,read->readLen,j,h,&left_bits);
				if(flip[j+h]==0)
				{

					struct edge_node **edge_node_p2p=&((*(bktptr[j+h]))->kmer_info.left);
					while((*edge_node_p2p)!=NULL)
					{
						if((*edge_node_p2p)->edge==(uint64_t)left_bits&&((*edge_node_p2p)->len+1)==h)
						{
							if((*edge_node_p2p)->edge_cov<0xf)
							{(*edge_node_p2p)->edge_cov++;}
							break;
						}
						edge_node_p2p=&((*edge_node_p2p)->nxt_edge);
					}
					if((*edge_node_p2p)==NULL)
					{
						(*edge_node_p2p) = (struct edge_node*)malloc(sizeof(struct edge_node));
						(*edge_cnt)++;
						memset(*edge_node_p2p,0,sizeof(struct edge_node));
						(*edge_node_p2p)->edge=(uint64_t)left_bits;
						(*edge_node_p2p)->edge_cov=1;
						(*edge_node_p2p)->len=h-1;
					}
				}
				else
				{
					left_bits=get_rev_comp_seq(left_bits,h);
					struct edge_node **edge_node_p2p=&((*(bktptr[j+h]))->kmer_info.right);
					while((*edge_node_p2p)!=NULL)
					{
						if((*edge_node_p2p)->edge==(uint64_t)left_bits&&((*edge_node_p2p)->len+1)==h)
						{
							if((*edge_node_p2p)->edge_cov<0xf)
							{(*edge_node_p2p)->edge_cov++;}
							break;
						}
						edge_node_p2p=&((*edge_node_p2p)->nxt_edge);

					}
					if((*edge_node_p2p)==NULL)
					{
						(*edge_node_p2p) = (struct edge_node*)malloc(sizeof(struct edge_node));
						(*edge_cnt)++;
						memset(*edge_node_p2p,0,sizeof(struct edge_node));
						(*edge_node_p2p)->edge=(uint64_t)left_bits;
						(*edge_node_p2p)->edge_cov=1;
						(*edge_node_p2p)->len=h-1;
					}
				}


			}


			if(h>0&&j+h<OverlappingKmers&&found[j]&&found[j+h])
			{
				uint64_t right_bits;
				get_sub_arr(read->read_bits,read->readLen,j+K_size,h,&right_bits);
				if(flip[j]==1)
				{

					right_bits=get_rev_comp_seq(right_bits,h);

					struct edge_node **edge_node_p2p=&((*(bktptr[j]))->kmer_info.left);
					while((*edge_node_p2p)!=NULL)
					{
						if((*edge_node_p2p)->edge==(uint64_t)right_bits&&((*edge_node_p2p)->len+1)==h)
						{
							if((*edge_node_p2p)->edge_cov<0xf)
							{(*edge_node_p2p)->edge_cov++;}
							break;
						}
						edge_node_p2p=&((*edge_node_p2p)->nxt_edge);
					}

					if((*edge_node_p2p)==NULL)
					{
						(*edge_node_p2p) = (struct edge_node*)malloc(sizeof(struct edge_node));
						(*edge_cnt)++;
						memset(*edge_node_p2p,0,sizeof(struct edge_node));
						(*edge_node_p2p)->edge=(uint64_t)right_bits;
						(*edge_node_p2p)->edge_cov=1;
						(*edge_node_p2p)->len=h-1;
					}
				}
				else
				{

					struct edge_node **edge_node_p2p=&((*(bktptr[j]))->kmer_info.right);
					while((*edge_node_p2p)!=NULL)
					{
						if((*edge_node_p2p)->edge==(uint64_t)right_bits&&((*edge_node_p2p)->len+1==h))
						{
							if((*edge_node_p2p)->edge_cov<0xf)
							{(*edge_node_p2p)->edge_cov++;}
							break;
						}
						edge_node_p2p=&((*edge_node_p2p)->nxt_edge);

					}
					if((*edge_node_p2p)==NULL)
					{
						(*edge_node_p2p) = (struct edge_node*)malloc(sizeof(struct edge_node));
						(*edge_cnt)++;
						memset(*edge_node_p2p,0,sizeof(struct edge_node));
						(*edge_node_p2p)->edge=(uint64_t)right_bits;
						(*edge_node_p2p)->edge_cov=1;
						(*edge_node_p2p)->len=h-1;
					}
				}


			}
		}


		j=j+h;
		//if(j>=OverlappingKmers-1)
		//{break;}
	}

}












void Sparse_Kmer_Index_Construction(struct read_t *read,struct hashtable *ht,int64_t *bucket_count,int K_size,int gap,BF_info * BF_info ,int round,int ref_pos,struct read_index *read_index)
{
	int readLen=read->readLen;
	int OverlappingKmers=readLen-K_size+1;
	int pre_search_len=100;
	
	if(OverlappingKmers<pre_search_len||round==2)
	{
		pre_search_len=OverlappingKmers;
	}
	
	
	if(gap>=OverlappingKmers)
	{	return;}



	int Read_arr_sz=readLen/32+1;
	int rem=readLen%32;
	if(rem==0)
	{Read_arr_sz--;}
	int tot_bits=Read_arr_sz*64;
	size_t ht_sz=ht->ht_sz;

	bool flip_beg,flip_end,found,found_bf;
	size_t hash_idx;
	uint64_t seq,f_seq,hv;
	int saved_kmers=0;
	
	//bucket ** bktptr;
	bucket ** bktptr_beg;
	bucket ** bktptr_end;
	
	int beg=0,end=0,end_bf=0;

	int first_found=-1,first_found_bf=-1;;



	//check the read to see if there is a saved kmer in the hashtable or bloom filter
	
	for (int j=0;j<pre_search_len;j++)
	{
		get_sub_arr(read->read_bits,read->readLen,j,K_size,&(seq));
	
		f_seq=get_rev_comp_seq(seq,K_size);
		if(seq>f_seq)
		{
			uint64_t t=seq;
			seq=f_seq;
			f_seq=t;
			//flip_beg=1;
		}
	
		hv=MurmurHash64A(&seq,sizeof(seq),0);

		hash_idx=(size_t) (hv%ht_sz);

		bktptr_end= &(ht->store_pos[hash_idx]);

		if(round==1)
		{
			found=look_up_in_a_list_r1(seq, ((bucket_r1 ***) &bktptr_end ));

			if(first_found_bf<0)
			{
				//check if there is a saved one in the BF table
				bool bf_inserted=0;
				if(BF_info->Bloom)
				{
				
					bf_inserted=1;
					if(BF_info->Bloom) 
					{
				
						for(int b=0;b<BF_info->d;++b)
						{
							size_t BF_hv=(size_t)(MurmurHash64A(&seq,sizeof(seq),b)%(BF_info->m));
							//size_t n_byte=BF_hv/8;
							//int n_bit=(int)(BF_hv%8);
							//int mask=1<<n_bit;
							//if((BF_info->BF_HT[n_byte]&mask)==0)
					
							if((BF_info->BF_HT[BF_hv >> 3] & (1U << (BF_hv & 0x7U)))==0)
							{
								bf_inserted=0;
							}
							//BF_info->BF_HT[BF_hv >> 3] |= (1U << (BF_hv & 0x7U));
						}
					}
					if(bf_inserted)
					{
						first_found_bf=j;
						break;
					}

				}
			}

		}
		else
		{
			found=look_up_in_a_list(seq,&bktptr_end);
		}


		if(found==1&&first_found<0)
		{
			first_found=j;
			break;
		}

	}


	

	///////updated for better sparseness
	if(round==1)
	{
		if(first_found>0)
		{beg=first_found%gap;}
		if(first_found<0&&first_found_bf>0)
		{beg=first_found_bf%gap;}
	}
	else
	{
		if(first_found<0)
		{return;}
		else
		{beg=first_found;}

	}



	saved_kmers=0;
	//look for the next sparse kmer

	while(beg<OverlappingKmers)
	{
		//cout<<saved_kmers<<endl;
		//check if there is a saved k-mer, if so resume from that, otherwise, save the last one.
		if(round==1)
		{

			for (int k=0;k<=gap;++k)
			{
				if(saved_kmers>0&&k==0)
				{	
					continue;
					
				}
				else
				{
					;
				}
				
				end=beg+k;
			
				if( end == OverlappingKmers )//overbound
				{
					return;
				}


				get_sub_arr(read->read_bits,read->readLen,end,K_size,&(seq));
				flip_end=0;
				f_seq=get_rev_comp_seq(seq,K_size);
				if(seq>f_seq)
				{
					uint64_t t=seq;
					seq=f_seq;
					f_seq=t;
					flip_end=1;
				}
	
				hv=MurmurHash64A(&seq,sizeof(seq),0);

				hash_idx=(size_t) (hv%ht_sz);

				bktptr_end= &(ht->store_pos[hash_idx]);

				if(round==1)
				{
					found=look_up_in_a_list_r1(seq, ((bucket_r1 ***) &bktptr_end ));

					//if found==0, check if something is in the bloom filter,here may be a tricky place, we shall treat bf and ht equally.
			
					if(found==0&&BF_info->Bloom)
					{
						bool bf_inserted=1;

						if(BF_info->Bloom) 
						{
				
							for(int b=0;b<BF_info->d;++b)
							{
								size_t BF_hv=(size_t)(MurmurHash64A(&seq,sizeof(seq),b)%(BF_info->m));
					
								//size_t n_byte=BF_hv/8;
								//int n_bit=(int)(BF_hv%8);
								//int mask=1<<n_bit;
								//if((BF_info->BF_HT[n_byte]&mask)==0)
						
								if((BF_info->BF_HT[BF_hv >> 3] & (1U << (BF_hv & 0x7U)))==0)
								{
									bf_inserted=0;
								}
								//BF_info->BF_HT[BF_hv >> 3] |= (1U << (BF_hv & 0x7U));
							}
						}
						if(bf_inserted)
						{
							//end_bf=beg+k;
							break;
						}
			
					}

				}
				else
				{
					found=look_up_in_a_list(seq,  &bktptr_end );
				
				}
			
			
				if(found==1)
				{

					break;
				}


				if(saved_kmers==0&&k==0)
				{	
					break;
					
				}
				


			}


			found_bf=1;
			if(found==0&&BF_info->Bloom)
			{
				for(int b=0;b<BF_info->d;++b)
				{
					size_t BF_hv=(size_t)(MurmurHash64A(&seq,sizeof(seq),b)%(BF_info->m));
					//size_t n_byte=BF_hv/8;
					//int n_bit=(int)(BF_hv%8);
					//int mask=1<<n_bit;
					//if((BF_info->BF_HT[n_byte]&mask)==0)
				
			
					if((BF_info->BF_HT[BF_hv >> 3] & (1U << (BF_hv & 0x7U)))==0)
					{
						found_bf=0;
						BF_info->BF_HT[BF_hv >> 3] |= (1U << (BF_hv & 0x7U));
						//BF_info->BF_HT[BF_hv >> 3] |= (1U << (BF_hv & 0x7U));
					}
					
				}
				//insertion finished, skip
				if(found_bf==0)
				{
					found_bf=1;
					beg=end;
					bktptr_beg=bktptr_end;
					saved_kmers++;
					continue;
				}
			}

		
		}
		else
		{

			for (int k=0;k<=OverlappingKmers;++k)
			{

				if(saved_kmers>0&&k==0)
				{	
					continue;
					
				}
				else
				{
					;
				}

				end=beg+k;
				
			
				if( end == OverlappingKmers )//overbound
				{
					return;
				}


				get_sub_arr(read->read_bits,read->readLen,end,K_size,&(seq));
	
				f_seq=get_rev_comp_seq(seq,K_size);
				flip_end=0;
				if(seq>f_seq)
				{
					uint64_t t=seq;
					seq=f_seq;
					f_seq=t;
					flip_end=1;
				}
	
				hv=MurmurHash64A(&seq,sizeof(seq),0);

				hash_idx=(size_t) (hv%ht_sz);

				bktptr_end= &(ht->store_pos[hash_idx]);

				found=look_up_in_a_list(seq,  &bktptr_end );

				if(found==1)
				{
					break;
				}
				if(saved_kmers==0)
				{break;}
			}


		}

		if(found==0)
		{
			if(round==1)
			{



				*bktptr_end=(struct bucket*)malloc(sizeof(struct bucket_r1));
				memset(*bktptr_end,0,sizeof(struct bucket_r1));
				//memset((ptr->kmer_t),0,sizeof(struct ptr->kmer_t));
				//ptr=ptr->nxt_bucket;
				((struct bucket_r1*) *(bktptr_end))->nxt_bucket=NULL;
				((struct bucket_r1*) *(bktptr_end))->kmer_t.kmer=seq;
				//(*(bktptr[j]))->kmer_info.cov1=0;
				((struct bucket_r1*) *(bktptr_end))->kmer_info.cov1++;
				//(*(bktptr[j]))->kmer_info.left=0;
				//(*(bktptr[j]))->kmer_info.right=0;
				//(*(bktptr[j]))->kmer_info.split_left=0;
				//(*(bktptr[j]))->kmer_info.split_right=0;
				//(*(bktptr[j]))->kmer_info.used=0;
				//(*(bktptr[j]))->kmer_info.cod=0;
				//(*(bktptr[j]))->kmer_info.contig_no=0;
				(*bucket_count)++;
				
			
				
				beg=end;
				saved_kmers++;
				bktptr_beg=bktptr_end;
				continue;
			}
			else
			{
				return;
			}
		

		}
		else
		{
			saved_kmers++;
			if(round==1)
			{
				if(((struct bucket_r1*) *(bktptr_end))->kmer_info.cov1<0xffff)
				{((struct bucket_r1*) *(bktptr_end))->kmer_info.cov1++;}
				beg=end;
				bktptr_beg=bktptr_end;
				continue;
			}
			else
			{
				if((*(bktptr_end))->kmer_info.cov1<0xffff)
				{(*(bktptr_end))->kmer_info.cov1++;}
			}

		}




		if(round==2)
		{


			if((*(bktptr_end))->kmer_info.cov1<0xffff)
			{
				(*(bktptr_end))->kmer_info.cov1++;
				if((*(bktptr_end))->kmer_info.cod==0)
				{
					(*(bktptr_end))->kmer_info.cod=read_index->repeat_cnt;
					map<uint64_t,bool> r_map;
					read_index->repeat_maps.push_back(r_map);
					(read_index->repeat_cnt)++;
				}
				size_t map_idx=((*(bktptr_end))->kmer_info.cod);
				uint64_t tmp_cod=end;//ref_pos+
				tmp_cod<<=32;
				tmp_cod|=read->read_idx;
				read_index->repeat_maps[map_idx][tmp_cod]=flip_end;
			}
		}
		beg=end;
		bktptr_beg=bktptr_end;
		flip_beg=flip_end;

	}


}







//convert the bucket type from round 1 to round 2. The buckets in round 2 are more expensive

void SwitchBuckets(hashtable *ht,hashtable2 *ht2,int K_size)
{
	size_t ht_sz;
	if(K_size<=32)
	{
		ht_sz=ht->ht_sz;
		bucket_r1 *store_pos_o,*store_pos_t;
		bucket *store_pos_n;
		bucket **bktp2p;
		for(size_t i=0;i<ht_sz;++i)
		{
			bktp2p=&(ht->store_pos[i]);
			store_pos_o=(bucket_r1*) ht->store_pos[i];
			while(store_pos_o!=NULL)
			{
				store_pos_n=(bucket*) malloc(sizeof(struct bucket));
				memset(store_pos_n,0,sizeof(struct bucket));
				store_pos_n->kmer_t=store_pos_o->kmer_t;
				store_pos_n->kmer_info.cov1=store_pos_o->kmer_info.cov1;
				store_pos_n->kmer_info.left=NULL;
				store_pos_n->kmer_info.right=NULL;
				*bktp2p=store_pos_n;
				bktp2p=&(store_pos_n->nxt_bucket);
				store_pos_t=store_pos_o;
				store_pos_o=store_pos_o->nxt_bucket;
				free(store_pos_t);
			}
		}
	}

	if(K_size>32&&K_size<64)
	{
		ht_sz=ht2->ht_sz;
		bucket2_r1 *store_pos_o,*store_pos_t;
		bucket2 *store_pos_n;
		bucket2 **bktp2p;
		for(size_t i=0;i<ht_sz;++i)
		{
			bktp2p=&(ht2->store_pos[i]);
			store_pos_o=(bucket2_r1*) ht2->store_pos[i];

			/*
			
			int n_buckets=0;
			while(store_pos_o!=NULL)
			{
				n_buckets++;
				store_pos_o=store_pos_o->nxt_bucket;
			
			}
			store_pos_o=(bucket2_r1*) ht2->store_pos[i];
			store_pos_n=(bucket2*) malloc(sizeof(struct bucket2)*n_buckets);
			bucket2 * init_bucket=store_pos_n;
			memset(store_pos_n,0,sizeof(bucket2)*n_buckets);
			n_buckets=0;	
			while(store_pos_o!=NULL)
			{
				store_pos_n->kmer_t2=store_pos_o->kmer_t2;
				store_pos_n->kmer_info.cov1=store_pos_o->kmer_info.cov1;
				store_pos_n->kmer_info.left=NULL;
				store_pos_n->kmer_info.right=NULL;
				*bktp2p=store_pos_n;
				bktp2p=&(store_pos_n->nxt_bucket);
				store_pos_t=store_pos_o;
				store_pos_o=store_pos_o->nxt_bucket;
				free(store_pos_t);
				//store_pos_n+=sizeof(bucket2);
				n_buckets++;
				store_pos_n=&(init_bucket[n_buckets]);
				
			}
			*/


			
			while(store_pos_o!=NULL)
			{
				store_pos_n=(bucket2*) malloc(sizeof(struct bucket2));
				memset(store_pos_n,0,sizeof(bucket2));
				store_pos_n->kmer_t2=store_pos_o->kmer_t2;
				store_pos_n->kmer_info.cov1=store_pos_o->kmer_info.cov1;
				store_pos_n->kmer_info.left=NULL;
				store_pos_n->kmer_info.right=NULL;
				*bktp2p=store_pos_n;
				bktp2p=&(store_pos_n->nxt_bucket);
				store_pos_t=store_pos_o;
				store_pos_o=store_pos_o->nxt_bucket;
				free(store_pos_t);
			}
			
		}
	}

}
void SwitchBuckets3(hashtable3 *ht3,int K_size)
{
	size_t ht_sz;
	
	if(K_size>64&&K_size<=96)
	{
		ht_sz=ht3->ht_sz;
		bucket3_r1 *store_pos_o,*store_pos_t;
		bucket3 *store_pos_n;
		bucket3 **bktp2p;
		for(size_t i=0;i<ht_sz;++i)
		{
			bktp2p=&(ht3->store_pos[i]);
			store_pos_o=(bucket3_r1*) ht3->store_pos[i];
			while(store_pos_o!=NULL)
			{
				store_pos_n=(bucket3*) malloc(sizeof(struct bucket3));
				memset(store_pos_n,0,sizeof(bucket3));
				store_pos_n->kmer_t3=store_pos_o->kmer_t3;
				store_pos_n->kmer_info.cov1=store_pos_o->kmer_info.cov1;
				*bktp2p=store_pos_n;
				bktp2p=&(store_pos_n->nxt_bucket);
				store_pos_t=store_pos_o;
				store_pos_o=store_pos_o->nxt_bucket;
				free(store_pos_t);
			}
		}
	}

}


void SwitchBuckets4(hashtable4 *ht4,int K_size)
{
	size_t ht_sz;
	
	if(K_size>96&&K_size<=128)
	{
		ht_sz=ht4->ht_sz;
		bucket4_r1 *store_pos_o,*store_pos_t;
		bucket4 *store_pos_n;
		bucket4 **bktp2p;
		for(size_t i=0;i<ht_sz;++i)
		{
			bktp2p=&(ht4->store_pos[i]);
			store_pos_o=(bucket4_r1*) ht4->store_pos[i];
			while(store_pos_o!=NULL)
			{
				store_pos_n=(bucket4*) malloc(sizeof(struct bucket4));
				memset(store_pos_n,0,sizeof(bucket4));
				store_pos_n->kmer_t4=store_pos_o->kmer_t4;
				store_pos_n->kmer_info.cov1=store_pos_o->kmer_info.cov1;
				*bktp2p=store_pos_n;
				bktp2p=&(store_pos_n->nxt_bucket);
				store_pos_t=store_pos_o;
				store_pos_o=store_pos_o->nxt_bucket;
				free(store_pos_t);
			}
		}
	}

}





void SwitchBuckets0(hashtable0 *ht0,int K_size)
{
	size_t ht_sz;
	
	ht_sz=ht0->ht_sz;
	bucket0_r1 *store_pos_o,*store_pos_t;
	bucket0 *store_pos_n;
	bucket0 **bktp2p;
	for(size_t i=0;i<ht_sz;++i)
	{
		bktp2p=&(ht0->store_pos[i]);
		store_pos_o=(bucket0_r1*) ht0->store_pos[i];
		while(store_pos_o!=NULL)
		{
			store_pos_n=(bucket0*) malloc(sizeof(struct bucket0));
			memset(store_pos_n,0,sizeof(bucket0));
			store_pos_n->kmer_t=store_pos_o->kmer_t;
			store_pos_n->kmer_info.cov1=store_pos_o->kmer_info.cov1;
			*bktp2p=store_pos_n;
			bktp2p=&(store_pos_n->nxt_bucket);
			store_pos_t=store_pos_o;
			store_pos_o=store_pos_o->nxt_bucket;
			free(store_pos_t);
		}
	}

}





void OutputSparseKmers(hashtable *ht,int K_size,bool Bloom)
{
	ofstream  o_kmers;
	//File *ht_content;
	map<int,int> cov_hist,edge_cov_hist;

	string o_kmers_name="KmerTable.txt";
	//ht_idx=fopen(ht_idx_name.c_str(),"w");
	o_kmers.open(o_kmers_name.c_str());
	
	bucket * bktptr=NULL;
	struct edge_node *edge_ptr;

	for(size_t i=0;i<ht->ht_sz;++i)
	{
		size_t list_sz=0;
		bktptr=ht->store_pos[i];
		while(bktptr!=NULL)
		{
			char kmer[1000];
			bitsarr2str(&(bktptr->kmer_t.kmer),K_size,kmer,1);
			if(!Bloom)
			{
				o_kmers<<kmer<<" "<<bktptr->kmer_info.cov1<<endl;			
			}
			else
			{
				o_kmers<<kmer<<" "<<(bktptr->kmer_info.cov1+1)<<endl;
			}
			bktptr=bktptr->nxt_bucket;
		}
	}



}

void OutputSparseKmers2(hashtable2 *ht,int K_size,bool Bloom)
{
	ofstream  o_kmers;
	//File *ht_content;

	string o_kmers_name="KmerTable.txt";
	//ht_idx=fopen(ht_idx_name.c_str(),"w");
	o_kmers.open(o_kmers_name.c_str());
	
	bucket2 * bktptr=NULL;
	
	for(size_t i=0;i<ht->ht_sz;++i)
	{
		size_t list_sz=0;
		bktptr=ht->store_pos[i];
		while(bktptr!=NULL)
		{
			char kmer[1000];
			bitsarr2str((bktptr->kmer_t2.kmer),K_size,kmer,2);
			if(!Bloom)
			{
				o_kmers<<kmer<<" "<<bktptr->kmer_info.cov1<<endl;			
			}
			else
			{
				o_kmers<<kmer<<" "<<(bktptr->kmer_info.cov1+1)<<endl;
			}
			bktptr=bktptr->nxt_bucket;
		}
	}



}

void OutputSparseKmers3(hashtable3 *ht,int K_size,bool Bloom)
{
	ofstream  o_kmers;
	//File *ht_content;

	string o_kmers_name="KmerTable.txt";
	//ht_idx=fopen(ht_idx_name.c_str(),"w");
	o_kmers.open(o_kmers_name.c_str());
	
	bucket3 * bktptr=NULL;
	
	for(size_t i=0;i<ht->ht_sz;++i)
	{
		size_t list_sz=0;
		bktptr=ht->store_pos[i];
		while(bktptr!=NULL)
		{
			char kmer[1000];
			bitsarr2str((bktptr->kmer_t3.kmer),K_size,kmer,3);
			if(!Bloom)
			{
				o_kmers<<kmer<<" "<<bktptr->kmer_info.cov1<<endl;			
			}
			else
			{
				o_kmers<<kmer<<" "<<(bktptr->kmer_info.cov1+1)<<endl;
			}
			bktptr=bktptr->nxt_bucket;
		}
	}



}
void OutputSparseKmers4(hashtable4 *ht,int K_size,bool Bloom)
{
	ofstream  o_kmers;
	//File *ht_content;

	string o_kmers_name="KmerTable.txt";
	//ht_idx=fopen(ht_idx_name.c_str(),"w");
	o_kmers.open(o_kmers_name.c_str());
	
	bucket4 * bktptr=NULL;
	
	for(size_t i=0;i<ht->ht_sz;++i)
	{
		size_t list_sz=0;
		bktptr=ht->store_pos[i];
		while(bktptr!=NULL)
		{
			char kmer[1000];
			bitsarr2str((bktptr->kmer_t4.kmer),K_size,kmer,4);
			if(!Bloom)
			{
				o_kmers<<kmer<<" "<<bktptr->kmer_info.cov1<<endl;			
			}
			else
			{
				o_kmers<<kmer<<" "<<(bktptr->kmer_info.cov1+1)<<endl;
			}
			bktptr=bktptr->nxt_bucket;
		}
	}



}

//save the sparse graph into disk
void SavingSparseKmerGraph(hashtable *ht,string &fname)
{
	ofstream  o_ht_idx,o_ht_content;
	//File *ht_content;
	map<int,int> cov_hist,edge_cov_hist;

	string ht_idx_name=fname+"HT_idx.txt",ht_content_name=fname+"HT_content";
	//ht_idx=fopen(ht_idx_name.c_str(),"w");
	o_ht_idx.open(ht_idx_name.c_str(),ios_base::out|ios_base::binary);
	o_ht_content.open(ht_content_name.c_str(),ios_base::out|ios_base::binary);
	o_ht_idx<<"Hashtable size: "<<endl<<ht->ht_sz<<endl;

	bucket * bktptr=NULL;
	struct edge_node *edge_ptr;

	for(size_t i=0;i<ht->ht_sz;++i)
	{
		size_t list_sz=0;
		bktptr=ht->store_pos[i];
		while(bktptr!=NULL)
		{
			o_ht_content.write((char*) bktptr,sizeof(struct bucket));
			edge_ptr=bktptr->kmer_info.left;
			while(edge_ptr!=NULL)
			{
				o_ht_content.write((char*) edge_ptr,sizeof(struct edge_node));
				edge_cov_hist[edge_ptr->edge_cov]++;
				edge_ptr=edge_ptr->nxt_edge;
			}
			edge_ptr=bktptr->kmer_info.right;
			while(edge_ptr!=NULL)
			{
				o_ht_content.write((char*) edge_ptr,sizeof(struct edge_node));
				edge_cov_hist[edge_ptr->edge_cov]++;
				edge_ptr=edge_ptr->nxt_edge;
			}
			int cov=bktptr->kmer_info.cov1;
			
			cov_hist[cov]++;
			bktptr=bktptr->nxt_bucket;
			list_sz++;
		}
		o_ht_idx<<list_sz<<endl;
	}
	ofstream o_cov("CovHist.txt");
	
	map<int,int >::iterator mit;
	for(mit=cov_hist.begin();mit!=cov_hist.end();++mit)
	{
		o_cov<<mit->first<<" "<<mit->second<<endl;
	}
	ofstream o_edge_cov("EdgeCovHist.txt");
	for(mit=edge_cov_hist.begin();mit!=edge_cov_hist.end();++mit)
	{
		o_edge_cov<<mit->first<<" "<<mit->second<<endl;
	}



}



//save the bubble removal information into disk

void SavingMergeHT(hashtable *ht)
{
	ofstream  o_ht_idx,o_ht_content;
	
	string ht_idx_name="MergeHT_idx.txt",ht_content_name="MergeHT_content";
	
	o_ht_idx.open(ht_idx_name.c_str(),ios_base::out|ios_base::binary);
	o_ht_content.open(ht_content_name.c_str(),ios_base::out|ios_base::binary);
	o_ht_idx<<"MergeHashTable size: "<<endl<<ht->ht_sz<<endl;

	bucket_rm * bktptr=NULL;
	
	for(size_t i=0;i<ht->ht_sz;++i)
	{
		size_t list_sz=0;
		bktptr=(bucket_rm *) ht->store_pos[i];
		while(bktptr!=NULL)
		{
			o_ht_content.write((char*) bktptr,sizeof(struct bucket_rm));
			bktptr=bktptr->nxt_bucket;
			list_sz++;
		}
		o_ht_idx<<list_sz<<endl;
	}
	

}




void SavingSparseKmerGraph_E2(hashtable *ht,string &fname)
{
	ofstream  o_ht_idx,o_ht_content;
	//File *ht_content;
	map<int,int> cov_hist,edge_cov_hist;

	string ht_idx_name=fname+"HT_idx.txt",ht_content_name=fname+"HT_content";
	//ht_idx=fopen(ht_idx_name.c_str(),"w");
	o_ht_idx.open(ht_idx_name.c_str(),ios_base::out|ios_base::binary);
	o_ht_content.open(ht_content_name.c_str(),ios_base::out|ios_base::binary);
	o_ht_idx<<"Hashtable size: "<<endl<<ht->ht_sz<<endl;

	bucket * bktptr=NULL;
	struct edge_node2 *edge_ptr;

	for(size_t i=0;i<ht->ht_sz;++i)
	{
		size_t list_sz=0;
		bktptr=ht->store_pos[i];
		while(bktptr!=NULL)
		{
			o_ht_content.write((char*) bktptr,sizeof(struct bucket));
			edge_ptr=(edge_node2 *)bktptr->kmer_info.left;
			while(edge_ptr!=NULL)
			{
				o_ht_content.write((char*) edge_ptr,sizeof(struct edge_node2));
				edge_cov_hist[edge_ptr->edge_cov]++;
				edge_ptr=edge_ptr->nxt_edge;
			}
			edge_ptr=(edge_node2 *)bktptr->kmer_info.right;
			while(edge_ptr!=NULL)
			{
				o_ht_content.write((char*) edge_ptr,sizeof(struct edge_node2));
				edge_cov_hist[edge_ptr->edge_cov]++;
				edge_ptr=edge_ptr->nxt_edge;
			}
			int cov=bktptr->kmer_info.cov1;
			
			cov_hist[cov]++;
			bktptr=bktptr->nxt_bucket;
			list_sz++;
		}
		o_ht_idx<<list_sz<<endl;
	}
	ofstream o_cov("CovHist.txt");
	
	map<int,int >::iterator mit;
	for(mit=cov_hist.begin();mit!=cov_hist.end();++mit)
	{
		o_cov<<mit->first<<" "<<mit->second<<endl;
	}
	ofstream o_edge_cov("EdgeCovHist.txt");
	for(mit=edge_cov_hist.begin();mit!=edge_cov_hist.end();++mit)
	{
		o_edge_cov<<mit->first<<" "<<mit->second<<endl;
	}



}


void SavingSparseKmerGraph2(hashtable2 *ht,string &fname)
{
	ofstream  o_ht_idx,o_ht_content;
	map<int,int> cov_hist,edge_cov_hist;

	//File *ht_content;
	string ht_idx_name=fname+"HT_idx.txt",ht_content_name=fname+"HT_content";
	//ht_idx=fopen(ht_idx_name.c_str(),"w");
	o_ht_idx.open(ht_idx_name.c_str(),ios_base::out|ios_base::binary);
	o_ht_content.open(ht_content_name.c_str(),ios_base::out|ios_base::binary);
	o_ht_idx<<"Hashtable size: "<<endl<<ht->ht_sz<<endl;

	bucket2 * bktptr=NULL;
	struct edge_node *edge_ptr;
	for(size_t i=0;i<ht->ht_sz;++i)
	{
		size_t list_sz=0;
		bktptr=ht->store_pos[i];
		while(bktptr!=NULL)
		{
			if(o_ht_content.write((char*) bktptr,sizeof(struct bucket2)))
			{
				edge_ptr=bktptr->kmer_info.left;
				while(edge_ptr!=NULL)
				{
					o_ht_content.write((char*) edge_ptr,sizeof(struct edge_node));
					edge_cov_hist[edge_ptr->edge_cov]++;
					edge_ptr=edge_ptr->nxt_edge;
				}
				edge_ptr=bktptr->kmer_info.right;
				while(edge_ptr!=NULL)
				{
					o_ht_content.write((char*) edge_ptr,sizeof(struct edge_node));
					edge_cov_hist[edge_ptr->edge_cov]++;
					edge_ptr=edge_ptr->nxt_edge;
				}

				int cov=bktptr->kmer_info.cov1;
				cov_hist[cov]++;
				bktptr=bktptr->nxt_bucket;
				list_sz++;
			}
			else
			{cout<<"Write error!"<<endl;}
		}
		o_ht_idx<<list_sz<<endl;
	}
	ofstream o_cov("CovHist.txt");
	map<int,int >::iterator mit;
	for(mit=cov_hist.begin();mit!=cov_hist.end();++mit)
	{
		o_cov<<mit->first<<" "<<mit->second<<endl;
	}
	ofstream o_edge_cov("EdgeCovHist.txt");
	for(mit=edge_cov_hist.begin();mit!=edge_cov_hist.end();++mit)
	{
		o_edge_cov<<mit->first<<" "<<mit->second<<endl;
	}

}



void SavingMergeHT2(hashtable2 *ht)
{
	ofstream  o_ht_idx,o_ht_content;
	string ht_idx_name="MergeHT_idx.txt",ht_content_name="MergeHT_content";
	//ht_idx=fopen(ht_idx_name.c_str(),"w");
	o_ht_idx.open(ht_idx_name.c_str(),ios_base::out|ios_base::binary);
	o_ht_content.open(ht_content_name.c_str(),ios_base::out|ios_base::binary);
	o_ht_idx<<"MergeHashTable size: "<<endl<<ht->ht_sz<<endl;

	bucket_rm2 * bktptr=NULL;
	for(size_t i=0;i<ht->ht_sz;++i)
	{
		size_t list_sz=0;
		bktptr=(bucket_rm2 *) ht->store_pos[i];
		while(bktptr!=NULL)
		{
			if(o_ht_content.write((char*) bktptr,sizeof(struct bucket_rm2)))
			{
				bktptr=bktptr->nxt_bucket;
				list_sz++;
			}
			else
			{cout<<"Write error!"<<endl;}
		}
		o_ht_idx<<list_sz<<endl;
	}
}





void SavingSparseKmerGraph3(hashtable3 *ht,string &fname)
{
	ofstream  o_ht_idx,o_ht_content;
	map<int,int> cov_hist,edge_cov_hist;
	//File *ht_content;
	string ht_idx_name=fname+"HT_idx.txt",ht_content_name=fname+"HT_content";
	//ht_idx=fopen(ht_idx_name.c_str(),"w");
	o_ht_idx.open(ht_idx_name.c_str(),ios_base::out|ios_base::binary);
	o_ht_content.open(ht_content_name.c_str(),ios_base::out|ios_base::binary);
	o_ht_idx<<"Hashtable size: "<<endl<<ht->ht_sz<<endl;

	bucket3 * bktptr=NULL;
	struct edge_node *edge_ptr;
	for(size_t i=0;i<ht->ht_sz;++i)
	{
		size_t list_sz=0;
		bktptr=ht->store_pos[i];
		while(bktptr!=NULL)
		{
			if(o_ht_content.write((char*) bktptr,sizeof(struct bucket3)))
			{
				edge_ptr=bktptr->kmer_info.left;
				while(edge_ptr!=NULL)
				{
					o_ht_content.write((char*) edge_ptr,sizeof(struct edge_node));
					edge_cov_hist[edge_ptr->edge_cov]++;
					edge_ptr=edge_ptr->nxt_edge;
				}
				edge_ptr=bktptr->kmer_info.right;
				while(edge_ptr!=NULL)
				{
					o_ht_content.write((char*) edge_ptr,sizeof(struct edge_node));
					edge_cov_hist[edge_ptr->edge_cov]++;
					edge_ptr=edge_ptr->nxt_edge;
				}

				int cov=bktptr->kmer_info.cov1;
				cov_hist[cov]++;
				bktptr=bktptr->nxt_bucket;
				list_sz++;
			}
			else
			{cout<<"Write error!"<<endl;}
		}
		o_ht_idx<<list_sz<<endl;
	}
	ofstream o_cov("CovHist.txt");
	map<int,int >::iterator mit;
	for(mit=cov_hist.begin();mit!=cov_hist.end();++mit)
	{
		o_cov<<mit->first<<" "<<mit->second<<endl;
	}
	ofstream o_edge_cov("EdgeCovHist.txt");
	for(mit=edge_cov_hist.begin();mit!=edge_cov_hist.end();++mit)
	{
		o_edge_cov<<mit->first<<" "<<mit->second<<endl;
	}

}


void SavingMergeHT3(hashtable3 *ht)
{
	ofstream  o_ht_idx,o_ht_content;
	string ht_idx_name="MergeHT_idx.txt",ht_content_name="MergeHT_content";
	//ht_idx=fopen(ht_idx_name.c_str(),"w");
	o_ht_idx.open(ht_idx_name.c_str(),ios_base::out|ios_base::binary);
	o_ht_content.open(ht_content_name.c_str(),ios_base::out|ios_base::binary);
	o_ht_idx<<"MergeHashTable size: "<<endl<<ht->ht_sz<<endl;

	bucket_rm3 * bktptr=NULL;
	for(size_t i=0;i<ht->ht_sz;++i)
	{
		size_t list_sz=0;
		bktptr=(bucket_rm3 *)ht->store_pos[i];
		while(bktptr!=NULL)
		{
			if(o_ht_content.write((char*) bktptr,sizeof(struct bucket_rm3)))
			{
				bktptr=bktptr->nxt_bucket;
				list_sz++;
			}
			else
			{cout<<"Write error!"<<endl;}
		}
		o_ht_idx<<list_sz<<endl;
	}
}



void SavingSparseKmerGraph4(hashtable4 *ht,string &fname)
{
	ofstream  o_ht_idx,o_ht_content;
	map<int,int> cov_hist,edge_cov_hist;
	//File *ht_content;
	string ht_idx_name=fname+"HT_idx.txt",ht_content_name=fname+"HT_content";
	//ht_idx=fopen(ht_idx_name.c_str(),"w");
	o_ht_idx.open(ht_idx_name.c_str(),ios_base::out|ios_base::binary);
	o_ht_content.open(ht_content_name.c_str(),ios_base::out|ios_base::binary);
	o_ht_idx<<"Hashtable size: "<<endl<<ht->ht_sz<<endl;

	bucket4 * bktptr=NULL;
	struct edge_node *edge_ptr;
	for(size_t i=0;i<ht->ht_sz;++i)
	{
		size_t list_sz=0;
		bktptr=ht->store_pos[i];
		while(bktptr!=NULL)
		{
			if(o_ht_content.write((char*) bktptr,sizeof(struct bucket4)))
			{
				edge_ptr=bktptr->kmer_info.left;
				while(edge_ptr!=NULL)
				{
					o_ht_content.write((char*) edge_ptr,sizeof(struct edge_node));
					edge_cov_hist[edge_ptr->edge_cov]++;
					edge_ptr=edge_ptr->nxt_edge;
				}
				edge_ptr=bktptr->kmer_info.right;
				while(edge_ptr!=NULL)
				{
					o_ht_content.write((char*) edge_ptr,sizeof(struct edge_node));
					edge_cov_hist[edge_ptr->edge_cov]++;
					edge_ptr=edge_ptr->nxt_edge;
				}

				int cov=bktptr->kmer_info.cov1;
				cov_hist[cov]++;
				bktptr=bktptr->nxt_bucket;
				list_sz++;
			}
			else
			{cout<<"Write error!"<<endl;}
		}
		o_ht_idx<<list_sz<<endl;
	}
	ofstream o_cov("CovHist.txt");
	map<int,int >::iterator mit;
	for(mit=cov_hist.begin();mit!=cov_hist.end();++mit)
	{
		o_cov<<mit->first<<" "<<mit->second<<endl;
	}
	ofstream o_edge_cov("EdgeCovHist.txt");
	for(mit=edge_cov_hist.begin();mit!=edge_cov_hist.end();++mit)
	{
		o_edge_cov<<mit->first<<" "<<mit->second<<endl;
	}

}



void SavingMergeHT4(hashtable4 *ht)
{
	ofstream  o_ht_idx,o_ht_content;
	string ht_idx_name="MergeHT_idx.txt",ht_content_name="MergeHT_content";
	//ht_idx=fopen(ht_idx_name.c_str(),"w");
	o_ht_idx.open(ht_idx_name.c_str(),ios_base::out|ios_base::binary);
	o_ht_content.open(ht_content_name.c_str(),ios_base::out|ios_base::binary);
	o_ht_idx<<"MergeHashTable size: "<<endl<<ht->ht_sz<<endl;

	bucket_rm4 * bktptr=NULL;
	for(size_t i=0;i<ht->ht_sz;++i)
	{
		size_t list_sz=0;
		bktptr=(bucket_rm4 *) ht->store_pos[i];
		while(bktptr!=NULL)
		{
			if(o_ht_content.write((char*) bktptr,sizeof(struct bucket_rm4)))
			{
				bktptr=bktptr->nxt_bucket;
				list_sz++;
			}
			else
			{cout<<"Write error!"<<endl;}
		}
		o_ht_idx<<list_sz<<endl;
	}
}










void SavingSparseKmerGraph0(hashtable0 *ht,string &fname,int Kmer_arr_sz)
{
	ofstream  o_ht_idx,o_ht_content;
	map<int,int> cov_hist,edge_cov_hist;
	//File *ht_content;
	string ht_idx_name=fname+"HT_idx.txt",ht_content_name=fname+"HT_content";
	//ht_idx=fopen(ht_idx_name.c_str(),"w");
	o_ht_idx.open(ht_idx_name.c_str(),ios_base::out|ios_base::binary);
	o_ht_content.open(ht_content_name.c_str(),ios_base::out|ios_base::binary);
	o_ht_idx<<"Hashtable size: "<<endl<<ht->ht_sz<<endl;

	bucket0 * bktptr=NULL;
	struct edge_node *edge_ptr;
	for(size_t i=0;i<ht->ht_sz;++i)
	{
		size_t list_sz=0;
		bktptr=ht->store_pos[i];
		while(bktptr!=NULL)
		{
			if(o_ht_content.write((char*) (bktptr->kmer_t),sizeof(uint64_t)*Kmer_arr_sz))
			{
				o_ht_content.write((char*) &(bktptr->kmer_info),sizeof(bktptr->kmer_info));
				edge_ptr=bktptr->kmer_info.left;
				while(edge_ptr!=NULL)
				{
					o_ht_content.write((char*) edge_ptr,sizeof(struct edge_node));
					edge_cov_hist[edge_ptr->edge_cov]++;
					edge_ptr=edge_ptr->nxt_edge;
				}
				edge_ptr=bktptr->kmer_info.right;
				while(edge_ptr!=NULL)
				{
					o_ht_content.write((char*) edge_ptr,sizeof(struct edge_node));
					edge_cov_hist[edge_ptr->edge_cov]++;
					edge_ptr=edge_ptr->nxt_edge;
				}

				int cov=bktptr->kmer_info.cov1;
				cov_hist[cov]++;
				bktptr=bktptr->nxt_bucket;
				list_sz++;
			}
			else
			{cout<<"Write error!"<<endl;}
		}
		o_ht_idx<<list_sz<<endl;
	}
	ofstream o_cov("CovHist.txt");
	map<int,int >::iterator mit;
	for(mit=cov_hist.begin();mit!=cov_hist.end();++mit)
	{
		o_cov<<mit->first<<" "<<mit->second<<endl;
	}
	ofstream o_edge_cov("EdgeCovHist.txt");
	for(mit=edge_cov_hist.begin();mit!=edge_cov_hist.end();++mit)
	{
		o_edge_cov<<mit->first<<" "<<mit->second<<endl;
	}

}



void SavingMergeHT0(hashtable0 *ht,int Kmer_arr_sz)
{
	ofstream  o_ht_idx,o_ht_content;
	string ht_idx_name="MergeHT_idx.txt",ht_content_name="MergeHT_content";
	//ht_idx=fopen(ht_idx_name.c_str(),"w");
	o_ht_idx.open(ht_idx_name.c_str(),ios_base::out|ios_base::binary);
	o_ht_content.open(ht_content_name.c_str(),ios_base::out|ios_base::binary);
	o_ht_idx<<"MergeHashTable size: "<<endl<<ht->ht_sz<<endl;

	bucket_rm0 * bktptr=NULL;
	for(size_t i=0;i<ht->ht_sz;++i)
	{
		size_t list_sz=0;
		bktptr=(bucket_rm0 *) ht->store_pos[i];
		while(bktptr!=NULL)
		{
			if(o_ht_content.write((char*) bktptr->kmer_t,sizeof(uint64_t)*Kmer_arr_sz))
			{
				o_ht_content.write((char*) bktptr->merged_kmer,sizeof(uint64_t)*Kmer_arr_sz);
				o_ht_content.write((char*) &(bktptr->flip),sizeof(bktptr->flip));

				bktptr=bktptr->nxt_bucket;
				list_sz++;
			}
			else
			{cout<<"Write error!"<<endl;}
		}
		o_ht_idx<<list_sz<<endl;
	}
}





void SavingSparseKmerIndex(hashtable *ht,read_index * read_index, string &fname)
{
	ofstream  o_ht_idx,o_ht_content,o_reads_idx,o_reads_content;
	//File *ht_content;
	map<int,int> cov_hist,edge_cov_hist;

	string ht_idx_name=fname+"HT_idx.txt",ht_content_name=fname+"HT_content";
	string read_idx_name=fname+"Reads_idx.txt",read_idx_name_content_name=fname+"Reads_idx_content";
	//ht_idx=fopen(ht_idx_name.c_str(),"w");
	o_ht_idx.open(ht_idx_name.c_str(),ios_base::out|ios_base::binary);
	o_ht_content.open(ht_content_name.c_str(),ios_base::out|ios_base::binary);
	o_reads_idx.open(read_idx_name.c_str(),ios_base::out|ios_base::binary);
	o_reads_content.open(read_idx_name_content_name.c_str(),ios_base::out|ios_base::binary);
	o_ht_idx<<"Hashtable size: "<<endl<<ht->ht_sz<<endl;

	bucket * bktptr=NULL;

	for(size_t i=0;i<ht->ht_sz;++i)
	{
		size_t list_sz=0;
		bktptr=ht->store_pos[i];
		while(bktptr!=NULL)
		{
			o_ht_content.write((char*) bktptr,sizeof(struct bucket));
			int cov=bktptr->kmer_info.cov1;
			cov_hist[cov]++;
			bktptr=bktptr->nxt_bucket;
			list_sz++;
		}
		o_ht_idx<<list_sz<<endl;
	}
	ofstream o_cov("CovHist.txt");
	
	o_reads_idx<<read_index->repeat_maps.size();
	size_t read_idx_sz=read_index->repeat_maps.size();
	for(size_t i=0;i<read_index->repeat_maps.size();++i)
	{
		size_t repeat_cnt=0;
		for(map<uint64_t,bool>::iterator repeat_map_it=read_index->repeat_maps[i].begin();repeat_map_it!=read_index->repeat_maps[i].end();++repeat_map_it)
		{
			o_reads_content<<repeat_map_it->first<<repeat_map_it->second;
			repeat_cnt++;
		}
		o_reads_idx<<repeat_cnt;
	}


	map<int,int >::iterator mit;
	for(mit=cov_hist.begin();mit!=cov_hist.end();++mit)
	{
		o_cov<<mit->first<<" "<<mit->second<<endl;
	}
	


}




void LoadingSparseKmerIndex(hashtable *ht,read_index * read_index, string &fname)
{
	ifstream  in_ht_idx,in_ht_content,in_reads_idx,in_reads_content;
	//File *ht_content;
	string ht_idx_name=fname+"HT_idx.txt",ht_content_name=fname+"HT_content";
	string read_idx_name=fname+"Reads_idx.txt",read_idx_name_content_name=fname+"Reads_idx_content";
	//ht_idx=fopen(ht_idx_name.c_str(),"w");
	in_ht_idx.open(ht_idx_name.c_str(),ios_base::in|ios_base::binary);
	in_ht_content.open(ht_content_name.c_str(),ios_base::in|ios_base::binary);
	in_reads_idx.open(read_idx_name.c_str(),ios_base::in|ios_base::binary);
	in_reads_content.open(read_idx_name_content_name.c_str(),ios_base::in|ios_base::binary);
	

	size_t ht_sz;
	string s;
	getline(in_ht_idx,s);
	getline(in_ht_idx,s);
	//ht_sz=atoi(s.c_str());
	stringstream strValue;
	strValue << s.c_str();
	strValue >> ht_sz;
	Init_HT(ht,ht_sz);
	struct edge_node **edge_p2p;
	for(size_t i=0;i<ht_sz;++i)
	{
		int list_sz;
		getline(in_ht_idx,s);
		if(s[s.size()-1]=='\r'||s[s.size()-1]=='\n')
		{s.resize(s.size()-1);}
		list_sz=atoi(s.c_str());
		struct bucket **bktp2p=&(ht->store_pos[i]);
		*bktp2p=NULL;
		for (int j=0;j<list_sz;++j)
		{
			*bktp2p=(struct bucket*)malloc(sizeof(struct bucket));
			in_ht_content.read((char*) (*bktp2p),sizeof(struct bucket));

			(*bktp2p)->kmer_info.used=0;
			(*bktp2p)->nxt_bucket=NULL;
			
			bktp2p=&((*bktp2p)->nxt_bucket);
		}
	}



	size_t read_idx_sz;
	in_reads_idx>>read_idx_sz;
	read_index->repeat_maps.resize(read_idx_sz);

	for(size_t i=0;i<read_idx_sz;++i)
	{
		size_t repeat_cnt;
		in_reads_idx>>repeat_cnt;

		for(size_t j=0;j<repeat_cnt;++j)//map<uint64_t,bool>::iterator repeat_map_it=read_index->repeat_maps[i].begin();repeat_map_it!=read_index->repeat_maps[i].end();++repeat_map_it)
		{
			uint64_t tmp1;
			bool tmp2;
			in_reads_content>>tmp1>>tmp2;
			read_index->repeat_maps[i][tmp1]=tmp2;
			
		}
	}


}




//load the saved information
void LoadingSparseKmerGraph(hashtable *ht,string &fname)
{
	string ht_idx_name=fname+"HT_idx.txt",ht_content_name=fname+"HT_content";
	ifstream in_ht_idx(ht_idx_name.c_str(),ios_base::in|ios_base::binary),in_ht_content(ht_content_name.c_str(),ios_base::in|ios_base::binary);
	
	size_t ht_sz;
	string s;
	getline(in_ht_idx,s);
	getline(in_ht_idx,s);
	//ht_sz=atoi(s.c_str());
	stringstream strValue;
	strValue << s.c_str();
	strValue >> ht_sz;
	Init_HT(ht,ht_sz);
	struct edge_node **edge_p2p;
	for(size_t i=0;i<ht_sz;++i)
	{
		int list_sz;
		getline(in_ht_idx,s);
		if(s[s.size()-1]=='\r'||s[s.size()-1]=='\n')
		{s.resize(s.size()-1);}
		list_sz=atoi(s.c_str());
		struct bucket **bktp2p=&(ht->store_pos[i]);
		*bktp2p=NULL;
		for (int j=0;j<list_sz;++j)
		{
			*bktp2p=(struct bucket*)malloc(sizeof(struct bucket));
			in_ht_content.read((char*) (*bktp2p),sizeof(struct bucket));

			(*bktp2p)->kmer_info.used=0;
			(*bktp2p)->nxt_bucket=NULL;
			edge_p2p=&((*bktp2p)->kmer_info.left);
			while((*edge_p2p)!=NULL)
			{
				(*edge_p2p)=(struct edge_node*)malloc(sizeof(struct edge_node));
				in_ht_content.read((char*) (*edge_p2p),sizeof(struct edge_node));
				edge_p2p=&((*edge_p2p)->nxt_edge);
			}
			edge_p2p=&((*bktp2p)->kmer_info.right);
			while((*edge_p2p)!=NULL)
			{
				(*edge_p2p)=(struct edge_node*)malloc(sizeof(struct edge_node));
				in_ht_content.read((char*) (*edge_p2p),sizeof(struct edge_node));
				edge_p2p=&((*edge_p2p)->nxt_edge);
			}
			bktp2p=&((*bktp2p)->nxt_bucket);
		}
	}

}



void LoadingMergeHT(hashtable *ht)
{
	ifstream in_ht_idx("MergeHT_idx.txt",ios_base::in|ios_base::binary),in_ht_content("MergeHT_content",ios_base::in|ios_base::binary);
	size_t ht_sz;
	string s;
	getline(in_ht_idx,s);
	getline(in_ht_idx,s);
	//ht_sz=atoi(s.c_str());
	
	if (s.size() == 0)
	{
		return;
	}
	stringstream strValue;
	strValue << s.c_str();
	strValue >> ht_sz;
	Init_HT(ht,ht_sz);
	for(size_t i=0;i<ht_sz;++i)
	{
		int list_sz;
		getline(in_ht_idx,s);
		
		if(s[s.size()-1]=='\r'||s[s.size()-1]=='\n')
		{s.resize(s.size()-1);}
		list_sz=atoi(s.c_str());
		struct bucket_rm **bktp2p=(bucket_rm **) &( ht->store_pos[i]);
		*bktp2p=NULL;
		for (int j=0;j<list_sz;++j)
		{
			*bktp2p=(struct bucket_rm*)malloc(sizeof(struct bucket_rm));
			in_ht_content.read((char*) (*bktp2p),sizeof(struct bucket_rm));

			(*bktp2p)->nxt_bucket=NULL;
			
			bktp2p=&((*bktp2p)->nxt_bucket);
		}
	}

}





void LoadingSparseKmerGraph2(hashtable2 *ht,string &fname)
{
	string ht_idx_name=fname+"HT_idx.txt",ht_content_name=fname+"HT_content";
	ifstream in_ht_idx(ht_idx_name.c_str(),ios_base::in|ios_base::binary),in_ht_content(ht_content_name.c_str(),ios_base::in|ios_base::binary);
	
	size_t ht_sz;
	string s;
	getline(in_ht_idx,s);
	getline(in_ht_idx,s);
	//ht_sz=atoi(s.c_str());
	stringstream strValue;
	strValue << s.c_str();
	strValue >> ht_sz;
	Init_HT2(ht,ht_sz);
//	cout<<ht_sz<<endl;
	struct edge_node **edge_p2p;
	for(size_t i=0;i<ht_sz;++i)
	{

		int list_sz;
		getline(in_ht_idx,s);
		if(s[s.size()-1]=='\r'||s[s.size()-1]=='\n')
		{s.resize(s.size()-1);}

		list_sz=atoi(s.c_str());
		struct bucket2 **bktp2p=&(ht->store_pos[i]);
		*bktp2p=NULL;
		for (int j=0;j<list_sz;++j)
		{
			*bktp2p=(struct bucket2*)malloc(sizeof(struct bucket2));

			if(in_ht_content.read((char*) (*bktp2p),sizeof(struct bucket2)))
			{

				(*bktp2p)->nxt_bucket=NULL;
				(*bktp2p)->kmer_info.used=0;
				edge_p2p=&((*bktp2p)->kmer_info.left);
				while((*edge_p2p)!=NULL)
				{
					(*edge_p2p)=(struct edge_node*)malloc(sizeof(struct edge_node));
					in_ht_content.read((char*) (*edge_p2p),sizeof(struct edge_node));
					edge_p2p=&((*edge_p2p)->nxt_edge);
				}
				edge_p2p=&((*bktp2p)->kmer_info.right);
				while((*edge_p2p)!=NULL)
				{
					(*edge_p2p)=(struct edge_node*)malloc(sizeof(struct edge_node));
					in_ht_content.read((char*) (*edge_p2p),sizeof(struct edge_node));
					edge_p2p=&((*edge_p2p)->nxt_edge);
				}


				bktp2p=&((*bktp2p)->nxt_bucket);
			}
			else
			{cout<<"Read error!"<<endl;}
		}
	}

}



void LoadingMergeHT2(hashtable2 *ht)
{
	ifstream in_ht_idx("MergeHT_idx.txt",ios_base::in|ios_base::binary),in_ht_content("MergeHT_content",ios_base::in|ios_base::binary);
	size_t ht_sz;
	string s;
	getline(in_ht_idx,s);
	getline(in_ht_idx,s);
	if (s.size() == 0)
	{
		return;
	}
	//ht_sz=atoi(s.c_str());
	stringstream strValue;
	strValue << s.c_str();
	strValue >> ht_sz;
	Init_HT2(ht,ht_sz);
	for(size_t i=0;i<ht_sz;++i)
	{
		int list_sz;
		getline(in_ht_idx,s);
		if(s[s.size()-1]=='\r'||s[s.size()-1]=='\n')
		{s.resize(s.size()-1);}
		list_sz=atoi(s.c_str());
		struct bucket_rm2 **bktp2p=(bucket_rm2 **)&(ht->store_pos[i]);
		*bktp2p=NULL;
		for (int j=0;j<list_sz;++j)
		{
			*bktp2p=(struct bucket_rm2*)malloc(sizeof(struct bucket_rm2));
			in_ht_content.read((char*) (*bktp2p),sizeof(struct bucket_rm2));

			(*bktp2p)->nxt_bucket=NULL;
			
			bktp2p=&((*bktp2p)->nxt_bucket);
		}
	}

}



void LoadingSparseKmerGraph3(hashtable3 *ht,string &fname)
{
	
	string ht_idx_name=fname+"HT_idx.txt",ht_content_name=fname+"HT_content";
	ifstream in_ht_idx(ht_idx_name.c_str(),ios_base::in|ios_base::binary),in_ht_content(ht_content_name.c_str(),ios_base::in|ios_base::binary);
	
	size_t ht_sz;
	string s;
	getline(in_ht_idx,s);
	getline(in_ht_idx,s);
//	ht_sz=atoi(s.c_str());
	stringstream strValue;
	strValue << s.c_str();
	strValue >> ht_sz;
	Init_HT3(ht,ht_sz);
//	cout<<ht_sz<<endl;
	struct edge_node **edge_p2p;
	for(size_t i=0;i<ht_sz;++i)
	{

		int list_sz;
		getline(in_ht_idx,s);
		if(s[s.size()-1]=='\r'||s[s.size()-1]=='\n')
		{s.resize(s.size()-1);}
		list_sz=atoi(s.c_str());
		struct bucket3 **bktp2p=&(ht->store_pos[i]);
		*bktp2p=NULL;
		for (int j=0;j<list_sz;++j)
		{
			*bktp2p=(struct bucket3*)malloc(sizeof(struct bucket3));

			if(in_ht_content.read((char*) (*bktp2p),sizeof(struct bucket3)))
			{

				(*bktp2p)->nxt_bucket=NULL;
				(*bktp2p)->kmer_info.used=0;
				edge_p2p=&((*bktp2p)->kmer_info.left);
				while((*edge_p2p)!=NULL)
				{
					(*edge_p2p)=(struct edge_node*)malloc(sizeof(struct edge_node));
					in_ht_content.read((char*) (*edge_p2p),sizeof(struct edge_node));
					edge_p2p=&((*edge_p2p)->nxt_edge);
				}
				edge_p2p=&((*bktp2p)->kmer_info.right);
				while((*edge_p2p)!=NULL)
				{
					(*edge_p2p)=(struct edge_node*)malloc(sizeof(struct edge_node));
					in_ht_content.read((char*) (*edge_p2p),sizeof(struct edge_node));
					edge_p2p=&((*edge_p2p)->nxt_edge);
				}


				bktp2p=&((*bktp2p)->nxt_bucket);
			}
			else
			{cout<<"Read error!"<<endl;}
		}
	}

}

void LoadingMergeHT3(hashtable3 *ht)
{
	ifstream in_ht_idx("MergeHT_idx.txt",ios_base::in|ios_base::binary),in_ht_content("MergeHT_content",ios_base::in|ios_base::binary);
	size_t ht_sz;
	string s;
	getline(in_ht_idx,s);
	getline(in_ht_idx,s);
	if (s.size() == 0)
	{
		return;
	}
	//ht_sz=atoi(s.c_str());
	stringstream strValue;
	strValue << s.c_str();
	strValue >> ht_sz;
	Init_HT3(ht,ht_sz);
	for(size_t i=0;i<ht_sz;++i)
	{
		int list_sz;
		getline(in_ht_idx,s);
		if(s[s.size()-1]=='\r'||s[s.size()-1]=='\n')
		{s.resize(s.size()-1);}
		list_sz=atoi(s.c_str());
		struct bucket_rm3 **bktp2p=(bucket_rm3 **)&(ht->store_pos[i]);
		*bktp2p=NULL;
		for (int j=0;j<list_sz;++j)
		{
			*bktp2p=(struct bucket_rm3*)malloc(sizeof(struct bucket_rm3));
			in_ht_content.read((char*) (*bktp2p),sizeof(struct bucket_rm3));

			(*bktp2p)->nxt_bucket=NULL;
			
			bktp2p=&((*bktp2p)->nxt_bucket);
		}
	}

}

void LoadingSparseKmerGraph4(hashtable4 *ht,string &fname)
{
	string ht_idx_name=fname+"HT_idx.txt",ht_content_name=fname+"HT_content";
	ifstream in_ht_idx(ht_idx_name.c_str(),ios_base::in|ios_base::binary),in_ht_content(ht_content_name.c_str(),ios_base::in|ios_base::binary);
	
	size_t ht_sz;
	string s;
	getline(in_ht_idx,s);
	getline(in_ht_idx,s);
	//ht_sz=atoi(s.c_str());
	stringstream strValue;
	strValue << s.c_str();
	strValue >> ht_sz;
	Init_HT4(ht,ht_sz);
//	cout<<ht_sz<<endl;
	struct edge_node **edge_p2p;
	for(size_t i=0;i<ht_sz;++i)
	{

		int list_sz;
		getline(in_ht_idx,s);
		if(s[s.size()-1]=='\r'||s[s.size()-1]=='\n')
		{s.resize(s.size()-1);}
		list_sz=atoi(s.c_str());
		struct bucket4 **bktp2p=&(ht->store_pos[i]);
		*bktp2p=NULL;
		for (int j=0;j<list_sz;++j)
		{
			*bktp2p=(struct bucket4*)malloc(sizeof(struct bucket4));

			if(in_ht_content.read((char*) (*bktp2p),sizeof(struct bucket4)))
			{

				(*bktp2p)->nxt_bucket=NULL;
				(*bktp2p)->kmer_info.used=0;
				edge_p2p=&((*bktp2p)->kmer_info.left);
				while((*edge_p2p)!=NULL)
				{
					(*edge_p2p)=(struct edge_node*)malloc(sizeof(struct edge_node));
					in_ht_content.read((char*) (*edge_p2p),sizeof(struct edge_node));
					edge_p2p=&((*edge_p2p)->nxt_edge);
				}
				edge_p2p=&((*bktp2p)->kmer_info.right);
				while((*edge_p2p)!=NULL)
				{
					(*edge_p2p)=(struct edge_node*)malloc(sizeof(struct edge_node));
					in_ht_content.read((char*) (*edge_p2p),sizeof(struct edge_node));
					edge_p2p=&((*edge_p2p)->nxt_edge);
				}


				bktp2p=&((*bktp2p)->nxt_bucket);
			}
			else
			{cout<<"Read error!"<<endl;}
		}
	}

}

void LoadingMergeHT4(hashtable4 *ht)
{
	ifstream in_ht_idx("MergeHT_idx.txt",ios_base::in|ios_base::binary),in_ht_content("MergeHT_content",ios_base::in|ios_base::binary);
	size_t ht_sz;
	string s;
	getline(in_ht_idx,s);
	getline(in_ht_idx,s);
	if (s.size() == 0)
	{
		return;
	}
	//ht_sz=atoi(s.c_str());
	stringstream strValue;
	strValue << s.c_str();
	strValue >> ht_sz;
	Init_HT4(ht,ht_sz);
	for(size_t i=0;i<ht_sz;++i)
	{
		int list_sz;
		getline(in_ht_idx,s);
		if(s[s.size()-1]=='\r'||s[s.size()-1]=='\n')
		{s.resize(s.size()-1);}
		list_sz=atoi(s.c_str());
		struct bucket_rm4 **bktp2p=(bucket_rm4 **)&(ht->store_pos[i]);
		*bktp2p=NULL;
		for (int j=0;j<list_sz;++j)
		{
			*bktp2p=(struct bucket_rm4*)malloc(sizeof(struct bucket_rm4));
			in_ht_content.read((char*) (*bktp2p),sizeof(struct bucket_rm4));

			(*bktp2p)->nxt_bucket=NULL;
			
			bktp2p=&((*bktp2p)->nxt_bucket);
		}
	}

}


void LoadingSparseKmerGraph0(hashtable0 *ht,key_table *key_table,string &fname,int Kmer_arr_sz)
{
	string ht_idx_name=fname+"HT_idx.txt",ht_content_name=fname+"HT_content";
	ifstream in_ht_idx(ht_idx_name.c_str(),ios_base::in|ios_base::binary),in_ht_content(ht_content_name.c_str(),ios_base::in|ios_base::binary);
	
	size_t ht_sz;
	string s;
	getline(in_ht_idx,s);
	getline(in_ht_idx,s);
	//ht_sz=atoi(s.c_str());
	stringstream strValue;
	strValue << s.c_str();
	strValue >> ht_sz;
	Init_HT0(ht,ht_sz);

	

	uint64_t * block_add=(uint64_t *) malloc(sizeof(uint64_t)*Kmer_arr_sz*key_table->KeysPerBlock);
	key_table->pblocks.push_back(block_add);

//	cout<<ht_sz<<endl;
	struct edge_node **edge_p2p;
	for(size_t i=0;i<ht_sz;++i)
	{

		int list_sz;
		getline(in_ht_idx,s);
		if(s[s.size()-1]=='\r'||s[s.size()-1]=='\n')
		{s.resize(s.size()-1);}
		list_sz=atoi(s.c_str());
		struct bucket0 **bktp2p=&(ht->store_pos[i]);
		*bktp2p=NULL;
		for (int j=0;j<list_sz;++j)
		{
			*bktp2p=(struct bucket0*)malloc(sizeof(struct bucket0));

			if(key_table->current_index==key_table->KeysPerBlock)
			{
				uint64_t *new_block_ptr=(uint64_t *)malloc(sizeof(uint64_t)*Kmer_arr_sz*key_table->KeysPerBlock);	
				key_table->pblocks.push_back(new_block_ptr);
				key_table->current_index=0;
				key_table->current_block=key_table->pblocks.size();						
			}

			uint64_t *block_ptr=(key_table->pblocks).back();
			uint64_t *key_ptr=&block_ptr[key_table->current_index*Kmer_arr_sz];


			//if(in_ht_content.read((char*) (*bktp2p),sizeof(struct bucket0)))
			if(in_ht_content.read((char*) key_ptr,sizeof(uint64_t)*Kmer_arr_sz))
			{
				(key_table->current_index)++;
				(*bktp2p)->kmer_t=key_ptr;
				in_ht_content.read((char*) &((*bktp2p)->kmer_info),sizeof((*bktp2p)->kmer_info));

				(*bktp2p)->nxt_bucket=NULL;
				(*bktp2p)->kmer_info.used=0;
				edge_p2p=&((*bktp2p)->kmer_info.left);
				while((*edge_p2p)!=NULL)
				{
					(*edge_p2p)=(struct edge_node*)malloc(sizeof(struct edge_node));
					in_ht_content.read((char*) (*edge_p2p),sizeof(struct edge_node));
					edge_p2p=&((*edge_p2p)->nxt_edge);
				}
				edge_p2p=&((*bktp2p)->kmer_info.right);
				while((*edge_p2p)!=NULL)
				{
					(*edge_p2p)=(struct edge_node*)malloc(sizeof(struct edge_node));
					in_ht_content.read((char*) (*edge_p2p),sizeof(struct edge_node));
					edge_p2p=&((*edge_p2p)->nxt_edge);
				}


				bktp2p=&((*bktp2p)->nxt_bucket);
			}
			else
			{cout<<"Read error!"<<endl;}
		}
	}

}

void LoadingMergeHT0(hashtable0 *ht,key_table *key_table,int Kmer_arr_sz)
{
	ifstream in_ht_idx("MergeHT_idx.txt",ios_base::in|ios_base::binary),in_ht_content("MergeHT_content",ios_base::in|ios_base::binary);
	size_t ht_sz;
	string s;
	getline(in_ht_idx,s);
	getline(in_ht_idx,s);
	if (s.size() == 0)
	{
		return;
	}
	//ht_sz=atoi(s.c_str());
	stringstream strValue;
	strValue << s.c_str();
	strValue >> ht_sz;
	Init_HT0(ht,ht_sz);
	for(size_t i=0;i<ht_sz;++i)
	{
		int list_sz;
		getline(in_ht_idx,s);
		if(s[s.size()-1]=='\r'||s[s.size()-1]=='\n')
		{s.resize(s.size()-1);}
		list_sz=atoi(s.c_str());
		struct bucket_rm0 **bktp2p=(bucket_rm0 **)&(ht->store_pos[i]);
		*bktp2p=NULL;
		for (int j=0;j<list_sz;++j)
		{

			
			*bktp2p=(struct bucket_rm0*)malloc(sizeof(struct bucket_rm0));

			if(key_table->current_index==key_table->KeysPerBlock)
			{
				uint64_t *new_block_ptr=(uint64_t *)malloc(sizeof(uint64_t)*Kmer_arr_sz*key_table->KeysPerBlock);	
				key_table->pblocks.push_back(new_block_ptr);
				key_table->current_index=0;
				key_table->current_block=key_table->pblocks.size();						
			}
			uint64_t *block_ptr=(key_table->pblocks).back();
			uint64_t *key_ptr=&block_ptr[key_table->current_index*Kmer_arr_sz];

			in_ht_content.read((char*) key_ptr,sizeof(uint64_t)*Kmer_arr_sz);
		
			(key_table->current_index)++;
			(*bktp2p)->kmer_t=key_ptr;
		



			if(key_table->current_index==key_table->KeysPerBlock)
			{
				uint64_t *new_block_ptr=(uint64_t *)malloc(sizeof(uint64_t)*Kmer_arr_sz*key_table->KeysPerBlock);	
				key_table->pblocks.push_back(new_block_ptr);
				key_table->current_index=0;
				key_table->current_block=key_table->pblocks.size();						
			}
			block_ptr=(key_table->pblocks).back();
			key_ptr=&block_ptr[key_table->current_index*Kmer_arr_sz];

			in_ht_content.read((char*) key_ptr,sizeof(uint64_t)*Kmer_arr_sz);
		
			(key_table->current_index)++;
			(*bktp2p)->merged_kmer=key_ptr;

			in_ht_content.read((char*) &((*bktp2p)->flip),sizeof((*bktp2p)->flip));

			(*bktp2p)->nxt_bucket=NULL;
			
			bktp2p=&((*bktp2p)->nxt_bucket);
		}
	}

}




/*
void LoadingSparseKmerGraph0(hashtable0 *ht,string &fname)
{
	string ht_idx_name=fname+"HT_idx.txt",ht_content_name=fname+"HT_content";
	ifstream in_ht_idx(ht_idx_name.c_str(),ios_base::in|ios_base::binary),in_ht_content(ht_content_name.c_str(),ios_base::in|ios_base::binary);
	
	size_t ht_sz;
	string s;
	getline(in_ht_idx,s);
	getline(in_ht_idx,s);
	//ht_sz=atoi(s.c_str());
	stringstream strValue;
	strValue << s.c_str();
	strValue >> ht_sz;
	Init_HT(ht,ht_sz);
//	cout<<ht_sz<<endl;
	struct edge_node **edge_p2p;
	for(size_t i=0;i<ht_sz;++i)
	{

		int list_sz;
		getline(in_ht_idx,s);
		if(s[s.size()-1]=='\r'||s[s.size()-1]=='\n')
		{s.resize(s.size()-1);}
		list_sz=atoi(s.c_str());
		struct bucket4 **bktp2p=&(ht->store_pos[i]);
		*bktp2p=NULL;
		for (int j=0;j<list_sz;++j)
		{
			*bktp2p=(struct bucket4*)malloc(sizeof(struct bucket4));

			if(in_ht_content.read((char*) (*bktp2p),sizeof(struct bucket4)))
			{

				(*bktp2p)->nxt_bucket=NULL;
				(*bktp2p)->kmer_info.used=0;
				edge_p2p=&((*bktp2p)->kmer_info.left);
				while((*edge_p2p)!=NULL)
				{
					(*edge_p2p)=(struct edge_node*)malloc(sizeof(struct edge_node));
					in_ht_content.read((char*) (*edge_p2p),sizeof(struct edge_node));
					edge_p2p=&((*edge_p2p)->nxt_edge);
				}
				edge_p2p=&((*bktp2p)->kmer_info.right);
				while((*edge_p2p)!=NULL)
				{
					(*edge_p2p)=(struct edge_node*)malloc(sizeof(struct edge_node));
					in_ht_content.read((char*) (*edge_p2p),sizeof(struct edge_node));
					edge_p2p=&((*edge_p2p)->nxt_edge);
				}


				bktp2p=&((*bktp2p)->nxt_bucket);
			}
			else
			{cout<<"Read error!"<<endl;}
		}
	}

}

void LoadingMergeHT0(hashtable0 *ht)
{
	ifstream in_ht_idx("MergeHT_idx.txt",ios_base::in|ios_base::binary),in_ht_content("MergeHT_content",ios_base::in|ios_base::binary);
	size_t ht_sz;
	string s;
	getline(in_ht_idx,s);
	getline(in_ht_idx,s);
	//ht_sz=atoi(s.c_str());
	stringstream strValue;
	strValue << s.c_str();
	strValue >> ht_sz;
	Init_HT4(ht,ht_sz);
	for(size_t i=0;i<ht_sz;++i)
	{
		int list_sz;
		getline(in_ht_idx,s);
		if(s[s.size()-1]=='\r'||s[s.size()-1]=='\n')
		{s.resize(s.size()-1);}
		list_sz=atoi(s.c_str());
		struct bucket_rm4 **bktp2p=(bucket_rm4 **)&(ht->store_pos[i]);
		*bktp2p=NULL;
		for (int j=0;j<list_sz;++j)
		{
			*bktp2p=(struct bucket_rm4*)malloc(sizeof(struct bucket_rm4));
			in_ht_content.read((char*) (*bktp2p),sizeof(struct bucket_rm4));

			(*bktp2p)->nxt_bucket=NULL;
			
			bktp2p=&((*bktp2p)->nxt_bucket);
		}
	}

}

*/


// free the graph
void FreeSparseKmerGraph(struct hashtable *ht)
{

	bucket * bktptr,*n_bktptr;//=NULL;
	edge_node *edgeptr,*n_edgeptr;
	for(size_t i=0;i<ht->ht_sz;++i)
	{
		size_t list_sz=0;
		bktptr=ht->store_pos[i];
		while(bktptr!=NULL)
		{
			n_bktptr=bktptr->nxt_bucket;
			edgeptr=bktptr->kmer_info.left;
			while(edgeptr!=NULL)
			{
				n_edgeptr=edgeptr->nxt_edge;
				free(edgeptr);
				edgeptr=n_edgeptr;
			}
			edgeptr=bktptr->kmer_info.right;
			while(edgeptr!=NULL)
			{
				n_edgeptr=edgeptr->nxt_edge;
				free(edgeptr);
				edgeptr=n_edgeptr;
			}

			free(bktptr);
			bktptr=n_bktptr;
		}

	}
	free(ht->store_pos);


}

void FreeSparseKmerGraph2(struct hashtable2 *ht)
{


	bucket2 * bktptr,*n_bktptr;//=NULL;

	for(size_t i=0;i<ht->ht_sz;++i)
	{
		size_t list_sz=0;
		bktptr=ht->store_pos[i];
		edge_node *edgeptr,*n_edgeptr;
		while(bktptr!=NULL)
		{
			n_bktptr=bktptr->nxt_bucket;

			edgeptr=bktptr->kmer_info.left;
			while(edgeptr!=NULL)
			{
				n_edgeptr=edgeptr->nxt_edge;
				free(edgeptr);
				edgeptr=n_edgeptr;
			}
			edgeptr=bktptr->kmer_info.right;
			while(edgeptr!=NULL)
			{
				n_edgeptr=edgeptr->nxt_edge;
				free(edgeptr);
				edgeptr=n_edgeptr;
			}
			free(bktptr);
			bktptr=n_bktptr;
		}

	}
	free(ht->store_pos);


}

void FreeSparseKmerGraph3(struct hashtable3 *ht)
{

	bucket3 * bktptr,*n_bktptr;//=NULL;

	for(size_t i=0;i<ht->ht_sz;++i)
	{
		size_t list_sz=0;
		bktptr=ht->store_pos[i];
		edge_node *edgeptr,*n_edgeptr;
		while(bktptr!=NULL)
		{
			n_bktptr=bktptr->nxt_bucket;

			edgeptr=bktptr->kmer_info.left;
			while(edgeptr!=NULL)
			{
				n_edgeptr=edgeptr->nxt_edge;
				free(edgeptr);
				edgeptr=n_edgeptr;
			}
			edgeptr=bktptr->kmer_info.right;
			while(edgeptr!=NULL)
			{
				n_edgeptr=edgeptr->nxt_edge;
				free(edgeptr);
				edgeptr=n_edgeptr;
			}
			free(bktptr);
			bktptr=n_bktptr;
		}

	}
	free(ht->store_pos);


}

void FreeSparseKmerGraph4(struct hashtable4 *ht)
{

	bucket4 * bktptr,*n_bktptr;//=NULL;

	for(size_t i=0;i<ht->ht_sz;++i)
	{
		size_t list_sz=0;
		bktptr=ht->store_pos[i];
		edge_node *edgeptr,*n_edgeptr;
		while(bktptr!=NULL)
		{
			n_bktptr=bktptr->nxt_bucket;

			edgeptr=bktptr->kmer_info.left;
			while(edgeptr!=NULL)
			{
				n_edgeptr=edgeptr->nxt_edge;
				free(edgeptr);
				edgeptr=n_edgeptr;
			}
			edgeptr=bktptr->kmer_info.right;
			while(edgeptr!=NULL)
			{
				n_edgeptr=edgeptr->nxt_edge;
				free(edgeptr);
				edgeptr=n_edgeptr;
			}
			free(bktptr);
			bktptr=n_bktptr;
		}

	}
	free(ht->store_pos);


}


void FreeSparseKmerGraph0(struct hashtable0 *ht,key_table *key_table)
{

	bucket0 * bktptr,*n_bktptr;//=NULL;

	for(size_t i=0;i<ht->ht_sz;++i)
	{
		size_t list_sz=0;
		bktptr=ht->store_pos[i];
		edge_node *edgeptr,*n_edgeptr;
		while(bktptr!=NULL)
		{
			n_bktptr=bktptr->nxt_bucket;

			edgeptr=bktptr->kmer_info.left;
			while(edgeptr!=NULL)
			{
				n_edgeptr=edgeptr->nxt_edge;
				free(edgeptr);
				edgeptr=n_edgeptr;
			}
			edgeptr=bktptr->kmer_info.right;
			while(edgeptr!=NULL)
			{
				n_edgeptr=edgeptr->nxt_edge;
				free(edgeptr);
				edgeptr=n_edgeptr;
			}
			free(bktptr);
			bktptr=n_bktptr;
		}

	}
	free(ht->store_pos);


	while(key_table->pblocks.size()>0)
	{
		uint64_t *block_add=key_table->pblocks.front();
		free(block_add);
		key_table->pblocks.pop_front();
	
	}

}


void ScanDataset(vector<string > & in_filenames_vt,uint64_t *tot_bases,uint64_t *numReads,uint64_t *totReads,int *maxlen)
{
	(*tot_bases)=0;
	(*numReads)=0;
	(*maxlen)=0;
	uint64_t tot_bases2=0,totReads2=(*totReads);
	size_t numReads2=0;
	cout<<"Scanning the dataset."<<endl;
	string seq_s,str,tag1,tag2;

	
	for(size_t jj=0;jj<in_filenames_vt.size();++jj)
	{
		uint64_t nLines=0;
		int seq_sz=0;

		ifstream infile(in_filenames_vt[jj].c_str());		
		cout<<jj+1<<"/"<<in_filenames_vt.size()<<" files."<<endl;
		cout<<"Processing file: "<<in_filenames_vt[jj]<<endl;
			
		seq_s.clear();

		bool fq_flag=0;
	

		getline(infile, str);
		if (fq_flag == 0 && str[0] == '@')
		{
			fq_flag = 1;
		}
		infile.close();

		infile.clear();
		infile.open(in_filenames_vt[jj].c_str());

		bool read_success = 0;

		read_success = 1;

		string tag, qs, n_tag;
		string QS_s;

		while (read_success)
		{
			if (fq_flag)
			{
				read_success = get_a_fastq_read(infile, tag, seq_s, QS_s);

			}
			else
			{
				read_success = get_a_fasta_read(infile, tag, seq_s, n_tag);

			}
			if (read_success == 0)
			{
				break;
			}

			seq_sz = seq_s.size();
			

			if (seq_s.size() == 0)
			{
				cout << "Empty sequence!" << endl;
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
					
			}
			if(bad_flag)
			{continue;}
						
				
			numReads2++;


			seq_s.clear();
		
				
			(tot_bases2)+=seq_sz;				

			if(seq_sz>*maxlen)
			{*maxlen=seq_sz;}

			if ((numReads2)%50000000==0)
			{
				cout<<"Reading: "<<(numReads2)<<endl;				
			}
				
			if ((totReads2) != 0 && (numReads2) >= (totReads2))
			{
				//(*numReads)=0;
				break;
			}

		}
		//fclose(infile);
		infile.close();
		infile.clear();
				 
				 

	}
	//cout<<"Scan finished."<<endl;
	
	(*tot_bases)=tot_bases2;
	(*numReads)=numReads2;
}








#endif
