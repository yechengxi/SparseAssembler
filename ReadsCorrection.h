#ifndef __ReadsCorrection_H
#define __ReadsCorrection_H

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


inline void replace_one_nc(uint64_t * bitsarr_in,int bitsarr_len,int begin_pos,uint64_t nt)
{
	
	int arr_sz_in=bitsarr_len/32+1;
	int rem=bitsarr_len%32;
	if(rem==0)
	{arr_sz_in--;}


	uint64_t temp_arr[10];
	//memset(temp_arr,0,sizeof(temp_arr));

	

	int rem2=(32-rem+begin_pos)%32;
	int block_beg=(32-rem+begin_pos)/32;
	if(rem==0)
	{block_beg--;}

	int rem3=(32-rem+begin_pos+1)%32;
	int block_end=(32-rem+begin_pos+1)/32;
	if(rem3!=0)
	{rem3=32-rem3;}
	else
	{
		block_end--;
	}
	if(rem==0)
	{block_end--;}

	int orig_sz=(block_end-block_beg+1);

	//memcpy(temp_arr,&bitsarr_in[block_beg],orig_sz*sizeof(uint64_t));
	

	L_shift_NB(temp_arr,rem2*2,orig_sz);
	bitsarr_in[block_beg]&=~(((uint64_t)0x3)<<(32-rem2-1)*2);
	bitsarr_in[block_beg]|=(((uint64_t)nt)<<(32-rem2-1)*2);

}

int BFsearch_denoising(struct read_t *read,struct hashtable *ht,struct hashtable2 *ht2,size_t *CovTh,size_t *CorrTh,int K_size,int gap,int *last_correct,bool SEARCH_RIGHT)
{

	bool CORRECTABLE=1;
	bool SELECT_OPTIMAL=0;
	bool CORRECT=1;
	bool ALLOW_WEAK_COV=1;
	bool BFS=1;
	bool Cov_Search=0;


//	char ch[100];

	size_t ht_sz;
	if(K_size<=32)
	{
		ht_sz=ht->ht_sz;
	}
	else
	{
		ht_sz=ht2->ht_sz;
	}
	
	int offset_nt[1000][2];
	//int dubious_nt[500][2];
	int correction_cnt=0;
	
	int most_cov=0;
	int first_edges_cnt=0;
	bool first_edge_found=0;
	edge_node *opt_edge=NULL,*first_edge=NULL;
	//uint64_t opt_edge_bits=0;
	//int opt_edge_len=0;
	int optimal_choice=-1;
	int opt_it1=-1,opt_it2=-1;

	int OverlapKmers=read->readLen-K_size+1;

	uint64_t ckey,f_ckey,temp_ckey;
	struct kmer_t2 ckey_t2,f_ckey_t2,temp_ckey_t2;
	struct bucket** bkt_last_correct;
	struct bucket2** bkt_last_correct_t2;
	bool flip_last_correct=0,found;
	uint64_t hv,hash_idx;
	if(K_size<=32)
	{
		get_sub_arr(read->read_bits,read->readLen,*last_correct,K_size,&ckey);
		f_ckey=get_rev_comp_seq(ckey,K_size);

		if(ckey>f_ckey)
		{
			flip_last_correct=1;
			temp_ckey=ckey;
			ckey=f_ckey;
			f_ckey=temp_ckey;
		}

		hv=MurmurHash64A(&ckey,sizeof(ckey),0);

		hash_idx=(size_t) (hv%ht_sz);

		bkt_last_correct= &(ht->store_pos[hash_idx]);
		found=look_up_in_a_list(ckey,&bkt_last_correct);
	}
	else
	{
	
		get_sub_arr(read->read_bits,read->readLen,*last_correct,K_size,ckey_t2.kmer);
		f_ckey_t2=ckey_t2;
		get_rev_comp_seq_arr(f_ckey_t2.kmer,K_size,2);
		
		if(uint64_t_cmp(ckey_t2.kmer,f_ckey_t2.kmer,2)>0)
		{
			flip_last_correct=1;
			temp_ckey_t2=ckey_t2;
			ckey_t2=f_ckey_t2;
			f_ckey_t2=temp_ckey_t2;
		}

		hv=MurmurHash64A(&ckey_t2,sizeof(ckey_t2),0);

		hash_idx=(size_t) (hv%ht_sz);

		bkt_last_correct_t2= &(ht2->store_pos[hash_idx]);
		found=look_up_in_a_list2(&ckey_t2,&bkt_last_correct_t2);
	}
	if(found==0)
	{return -1;}

	if(BFS)
	{
		if(SEARCH_RIGHT)
		{
			//search_right
		
			//int correction_cnt=0;

			//////////////////////////////BFS 
			
			for (int round=0;round<=2;++round)
			{
				edge_node* edge_ptr;
				int edge_len;
				if(flip_last_correct==0)
				{
					if(K_size<=32)
					{
						edge_ptr=(*bkt_last_correct)->kmer_info.right;
					}
					else
					{
						edge_ptr=(*bkt_last_correct_t2)->kmer_info.right;					
					}

				}
				else
				{
					if(K_size<=32)
					{
						edge_ptr=(*bkt_last_correct)->kmer_info.left;
					}
					else
					{
						edge_ptr=(*bkt_last_correct_t2)->kmer_info.left;					
					}
			
				}

				while(edge_ptr!=NULL)
				{
					edge_len=edge_ptr->len+1;
					if((*last_correct+edge_len)>=OverlapKmers-1)
					{
						edge_ptr=edge_ptr->nxt_edge;
						continue;
					}
					uint64_t read_bits,edge_bits;
				
					get_sub_arr(read->read_bits,read->readLen,*last_correct+K_size,edge_len,&read_bits);
					edge_bits=edge_ptr->edge;
					if(flip_last_correct)
					{
						read_bits=get_rev_comp_seq(read_bits,edge_len);
					}
					int mod_cnt=0;
					//int last_mod=-1;
					
					for(int l=0;l<edge_len;++l)
					{
						uint64_t mask=0x3;
						mask<<=2*l;
						if((read_bits&mask)!=(edge_bits&mask))
						{
							mod_cnt++;
						}
					}
				

					if(mod_cnt>round)
					{
						edge_ptr=edge_ptr->nxt_edge;
						continue;
					}
					else
					{
						if(mod_cnt==round&&(first_edge_found==0))
						{
							first_edge=edge_ptr;
							first_edges_cnt++;
						}
					
					}
				
					if(K_size<=32)
					{
						uint64_t t_kmer,t,f_kmer;
	
						t_kmer=(*bkt_last_correct)->kmer_t.kmer;
						if(flip_last_correct==0)
						{
							for(int j=edge_len-1;j>=0;--j)
							{
								uint64_t right_bits=(edge_bits>>(2*j));

								switch(right_bits&0x3)
								{
									case 0:
										t=3;
										t<<=(K_size-1)*2;
										t_kmer&=(~t);
										t_kmer<<=2;

										break;
									case 1:
										t=3;
										t<<=(K_size-1)*2;
										t_kmer&=(~t);
										t_kmer<<=2;
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

						}
						else
						{
						
							for(int j=0;j<edge_len;++j)
							{
								uint64_t left_bits=(edge_bits>>(2*j));
								switch(left_bits&0x3)
								{
									case 0:
										t_kmer>>=2;
										t=((uint64_t)0)<<((K_size-1)*2);
										t_kmer|=t;
									
										break;
									case 1:
										t_kmer>>=2;
										t=((uint64_t)1)<<((K_size-1)*2);
										t_kmer|=t;
								
										break;
									case 2:
										t_kmer>>=2;
										t=((uint64_t)2)<<(K_size-1)*2;
										t_kmer|=t;
								
										break;
									case 3:
										t_kmer>>=2;
										t=((uint64_t)3)<<((K_size-1)*2);
										t_kmer|=t;
								
										break;

								}
							}

						
						}
						f_kmer=get_rev_comp_seq(t_kmer,K_size);
						if(t_kmer>f_kmer)
						{
							t_kmer=f_kmer;
						}
						hv=MurmurHash64A(&t_kmer,sizeof(t_kmer),0);
						hash_idx=(size_t) (hv%ht_sz);

						struct bucket ** ptr;

						ptr= &(ht->store_pos[hash_idx]);
						
						found=look_up_in_a_list(t_kmer,&ptr);
						if(found==0)
						{
							edge_ptr=edge_ptr->nxt_edge;
							continue;
						}
						if((*ptr)->kmer_info.cov1>*CorrTh)
						{
							
							if((!SELECT_OPTIMAL)&&most_cov>=*CovTh)
							{
								CORRECTABLE=0;return -1;
							}
							most_cov=(*ptr)->kmer_info.cov1;
							optimal_choice=correction_cnt;
							opt_edge=edge_ptr;
							correction_cnt++;
							
						
						}
						
						
					}	
					else
					{
						if(K_size>32&&K_size<=64)
						{

							kmer_t2 t_kmer,f_kmer;
							t_kmer=(*bkt_last_correct_t2)->kmer_t2;
							uint64_t t;

							if(flip_last_correct==0)
							{
								for(int j=0;j<edge_len;++j)
								{
									uint64_t right_bits=(edge_bits>>(2*j));
									switch(right_bits&0x3)
									{
										case 0:
											t=3;
											t<<=(K_size-1-32)*2;
											t_kmer.kmer[0]&=(~t);
											L_shift_NB(t_kmer.kmer,2,2);

											break;
										case 1:
											t=3;
											t<<=(K_size-1-32)*2;
											t_kmer.kmer[0]&=(~t);
											L_shift_NB(t_kmer.kmer,2,2);
											t=1;

											t_kmer.kmer[1]|=t;
									
											break;
										case 2:
											t=3;
											t<<=(K_size-1-32)*2;
											t_kmer.kmer[0]&=(~t);
											L_shift_NB(t_kmer.kmer,2,2);
											t=2;

											t_kmer.kmer[1]|=t;
								
											break;
										case 3:
											t=3;
											t<<=(K_size-1-32)*2;
											t_kmer.kmer[0]&=(~t);
											L_shift_NB(t_kmer.kmer,2,2);
											t=3;

											t_kmer.kmer[1]|=t;
								
											break;



									}


								}
							}
							else
							{
								for(int j=0;j<edge_len;++j)
								{

									uint64_t left_bits=(edge_bits>>(2*j));

									switch(left_bits&0x3)
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
							}

							if(uint64_t_cmp(t_kmer.kmer,f_kmer.kmer,2)>0)
							{
								t_kmer=f_kmer;
							}

							hv=MurmurHash64A(&t_kmer,sizeof(t_kmer),0);

							hash_idx=(size_t) (hv%ht_sz);

							struct bucket2 ** ptr;

							ptr= &(ht2->store_pos[hash_idx]);
						
							found=look_up_in_a_list2(&t_kmer,&ptr);

							if(found==0)
							{
								edge_ptr=edge_ptr->nxt_edge;
								continue;
							}

							if((*ptr)->kmer_info.cov1>=*CorrTh)
							{

								if((!SELECT_OPTIMAL)&&most_cov>=*CovTh)
								{
									CORRECTABLE=0;return -1;
								}

								most_cov=(*ptr)->kmer_info.cov1;
								optimal_choice=correction_cnt;
								opt_edge=edge_ptr;
								correction_cnt++;
							}
	
						}
					
					
					}

					edge_ptr=edge_ptr->nxt_edge;
				}

				if(first_edges_cnt>0)
				{
					first_edge_found=1;
				}
					
				if(correction_cnt>=1)
				{
					break;
				}

			}

			
		}

		

		if(!SEARCH_RIGHT)
		{
			//search left
			uint64_t tt;
			struct kmer_t2 temp_t2;
	
			//int correction_cnt=0;				
	

			//////////////////////////////BFS 
			for (int round=0;round<=2;++round)
			{
				edge_node* edge_ptr;
				int edge_len;
				if(flip_last_correct==1)
				{
					if(K_size<=32)
					{
						edge_ptr=(*bkt_last_correct)->kmer_info.right;
					}
					else
					{
						edge_ptr=(*bkt_last_correct_t2)->kmer_info.right;					
					}

				}
				else
				{
					if(K_size<=32)
					{
						edge_ptr=(*bkt_last_correct)->kmer_info.left;
					}
					else
					{
						edge_ptr=(*bkt_last_correct_t2)->kmer_info.left;					
					}
			
				}

				

				while(edge_ptr!=NULL)
				{
					edge_len=edge_ptr->len+1;
					if(*last_correct-edge_len<0||(*last_correct==0))
					{
						edge_ptr=edge_ptr->nxt_edge;
						continue;
					}
					uint64_t read_bits,edge_bits;
				
					get_sub_arr(read->read_bits,read->readLen,*last_correct-edge_len,edge_len,&read_bits);
					edge_bits=edge_ptr->edge;
					if(flip_last_correct)
					{
						read_bits=get_rev_comp_seq(read_bits,edge_len);
					}
					int mod_cnt=0;
					//int last_mod=-1;
					for(int l=0;l<edge_len;++l)
					{
						uint64_t mask=0x3;
						mask<<=2*l;
					
						if((read_bits&mask)!=(edge_bits&mask))
						{
							mod_cnt++;
							
							//last_mod=l;
							
						}
					}

				
				
					if(mod_cnt>round)
					{
						edge_ptr=edge_ptr->nxt_edge;
						continue;
					}
					else
					{
						if(mod_cnt==round&&(first_edge_found==0))
						{
							first_edge=edge_ptr;
							first_edges_cnt++;
						}
					
					}

					

					if(K_size<=32)
					{
						uint64_t t_kmer,t,f_kmer;
						if(flip_last_correct==1)
						{
							for(int j=edge_len-1;j>=0;--j)
							{

								uint64_t right_bits=(edge_bits>>(2*j));

								switch(right_bits&0x3)
								{
									case 0:
										t=3;
										t<<=(K_size-1)*2;
										t_kmer&=(~t);
										t_kmer<<=2;

										break;
									case 1:
										t=3;
										t<<=(K_size-1)*2;
										t_kmer&=(~t);
										t_kmer<<=2;
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

						}
						else
						{
						
							t_kmer=(*bkt_last_correct)->kmer_t.kmer;
							for(int j=0;j<edge_len;++j)
							{
								uint64_t left_bits=(edge_bits>>(2*j));
								switch(left_bits&0x3)
								{
									case 0:
										t_kmer>>=2;
										t=((uint64_t)0)<<((K_size-1)*2);
										t_kmer|=t;
									
										break;
									case 1:
										t_kmer>>=2;
										t=((uint64_t)1)<<((K_size-1)*2);
										t_kmer|=t;
								
										break;
									case 2:
										t_kmer>>=2;
										t=((uint64_t)2)<<(K_size-1)*2;
										t_kmer|=t;
								
										break;
									case 3:
										t_kmer>>=2;
										t=((uint64_t)3)<<((K_size-1)*2);
										t_kmer|=t;
								
										break;

								}
							}

						
						}
						f_kmer=get_rev_comp_seq(t_kmer,K_size);
						if(t_kmer>f_kmer)
						{
							t_kmer=f_kmer;
							
						}
						hv=MurmurHash64A(&t_kmer,sizeof(t_kmer),0);
						hash_idx=(size_t) (hv%ht_sz);

						struct bucket ** ptr;

						ptr= &(ht->store_pos[hash_idx]);
						
						found=look_up_in_a_list(t_kmer,&ptr);
						if(found==0)
						{
							edge_ptr=edge_ptr->nxt_edge;
							continue;
						}
						if((*ptr)->kmer_info.cov1>*CorrTh)
						{
							
							if((!SELECT_OPTIMAL)&&most_cov>=*CovTh)
							{
								CORRECTABLE=0;return -1;
							}
							most_cov=(*ptr)->kmer_info.cov1;
							optimal_choice=correction_cnt;
							opt_edge=edge_ptr;
							correction_cnt++;
								
						
						}
						
							
						
					}	
					else
					{
						if(K_size>32&&K_size<=64)
						{

							kmer_t2 t_kmer,f_kmer;
							t_kmer=(*bkt_last_correct_t2)->kmer_t2;
							uint64_t t;

							if(flip_last_correct==1)
							{
								for(int j=0;j<edge_len;++j)
								{

									uint64_t right_bits=(edge_bits>>(2*j));
									switch(right_bits&0x3)
									{
										case 0:
											t=3;
											t<<=(K_size-1-32)*2;
											t_kmer.kmer[0]&=(~t);
											L_shift_NB(t_kmer.kmer,2,2);

											break;
										case 1:
											t=3;
											t<<=(K_size-1-32)*2;
											t_kmer.kmer[0]&=(~t);
											L_shift_NB(t_kmer.kmer,2,2);
											t=1;

											t_kmer.kmer[1]|=t;
									
											break;
										case 2:
											t=3;
											t<<=(K_size-1-32)*2;
											t_kmer.kmer[0]&=(~t);
											L_shift_NB(t_kmer.kmer,2,2);
											t=2;

											t_kmer.kmer[1]|=t;
								
											break;
										case 3:
											t=3;
											t<<=(K_size-1-32)*2;
											t_kmer.kmer[0]&=(~t);
											L_shift_NB(t_kmer.kmer,2,2);
											t=3;

											t_kmer.kmer[1]|=t;
								
											break;



									}


								}
							}
							else
							{
								for(int j=0;j<edge_len;++j)
								{

									uint64_t left_bits=(edge_bits>>(2*j));

									switch(left_bits&0x3)
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
							}

							if(uint64_t_cmp(t_kmer.kmer,f_kmer.kmer,2)>0)
							{
								t_kmer=f_kmer;
							}


							hv=MurmurHash64A(&t_kmer,sizeof(t_kmer),0);

							hash_idx=(size_t) (hv%ht_sz);

							struct bucket2 ** ptr;

							ptr= &(ht2->store_pos[hash_idx]);
						
							found=look_up_in_a_list2(&t_kmer,&ptr);


							if(found==0)
							{
								edge_ptr=edge_ptr->nxt_edge;
								continue;
							}



							if((*ptr)->kmer_info.cov1>=*CorrTh)
							{
							

								if((!SELECT_OPTIMAL)&&most_cov>=*CovTh)
								{
									CORRECTABLE=0;return -1;
								}
								most_cov=(*ptr)->kmer_info.cov1;
								optimal_choice=correction_cnt;
								opt_edge=edge_ptr;
								correction_cnt++;

							//	break;
								//
								
							}
	
						}
					
					
					}

					edge_ptr=edge_ptr->nxt_edge;
				}
///
				if(first_edges_cnt>0)
				{
					first_edge_found=1;
				}

				if(correction_cnt>=1)
				{
					break;
				}


			}

		}
	}
	if(correction_cnt>1)//correction_cnt==0||
	{
		CORRECTABLE=0;
		return -1;
	}
	if(correction_cnt==0&&first_edge_found&&(first_edges_cnt==1))
	{
		return -1;
		CORRECTABLE=1;
		opt_edge=first_edge;

		if(!SEARCH_RIGHT)
		{
			
			int edge_len=opt_edge->len+1;
			uint64_t edge_bits=opt_edge->edge;

			if(flip_last_correct)
			{
				edge_bits=get_rev_comp_seq(edge_bits,edge_len);
			}
			int mod_cnt=0;
			for(int l=0;l<edge_len;++l)
			{
						
				uint64_t read_bits;
				get_sub_arr(read->read_bits,read->readLen,*last_correct-edge_len,edge_len,&read_bits);
					
				uint64_t mask=0x3;
				mask<<=2*l;
					
				uint64_t nt;
				if((read_bits&mask)!=(edge_bits&mask))
				{
					nt=(edge_bits>>2*l)&(0x3);
					replace_one_nc(read->read_bits,read->readLen,(*last_correct)-1-l,nt);
					
					mod_cnt++;
					if(mod_cnt>3)
					{
						return -1;
					}
				}
					
					


			}
			memset((&read->error_nt[(*last_correct)-edge_len]),0,edge_len*sizeof(bool));

		
		}
		else
		{
			int edge_len=opt_edge->len+1;
			uint64_t edge_bits=opt_edge->edge;
				
			if(flip_last_correct)
			{
				edge_bits=get_rev_comp_seq(edge_bits,edge_len);
			}	
				

			int mod_cnt=0;
			uint64_t read_bits;
			get_sub_arr(read->read_bits,read->readLen,*last_correct+K_size,edge_len,&read_bits);
					
			for(int l=0;l<edge_len;++l)
			{
				
						
				uint64_t mask=0x3;
				mask<<=2*l;
					
				uint64_t nt;
				if((read_bits&mask)!=(edge_bits&mask))
				{
					nt=(edge_bits>>2*l)&(0x3);
					replace_one_nc(read->read_bits,read->readLen,(edge_len-1-l)+(*last_correct)+K_size,nt);
				
					mod_cnt++;
					if(mod_cnt>3)
					{
						return -1;

					}
				}
					
			}
			memset((&read->error_nt[(*last_correct)+K_size]),0,edge_len*sizeof(bool));


					
		
		}
		

	}
	
	if(opt_edge!=NULL)
	{
		if(SEARCH_RIGHT)
		{
		
			
			int edge_len=opt_edge->len+1;
			uint64_t edge_bits=opt_edge->edge;
			if(flip_last_correct==1)
			{
				edge_bits=get_rev_comp_seq(edge_bits,edge_len);
			}
			for(int l=0;l<edge_len;++l)
			{
				
				uint64_t read_bits;
				get_sub_arr(read->read_bits,read->readLen,*last_correct+K_size,edge_len,&read_bits);
					
				uint64_t mask=0x3;
				mask<<=2*l;
					
				uint64_t nt;
				if((read_bits&mask)!=(edge_bits&mask))
				{
					nt=(edge_bits>>2*l)&(0x3);
					replace_one_nc(read->read_bits,read->readLen,(edge_len-1-l)+(*last_correct)+K_size,nt);
				
				}
					



			}
			
	

		}
		else
		{

			int edge_len=opt_edge->len+1;
			uint64_t edge_bits=opt_edge->edge;
			if(flip_last_correct==1)
			{
				edge_bits=get_rev_comp_seq(edge_bits,edge_len);
			}
			for(int l=0;l<edge_len;++l)
			{
				
				uint64_t read_bits;
				get_sub_arr(read->read_bits,read->readLen,*last_correct-edge_len,edge_len,&read_bits);
					
				uint64_t mask=0x3;
				mask<<=2*l;
					
				uint64_t nt;
				if((read_bits&mask)!=(edge_bits&mask))
				{
					nt=(edge_bits>>2*l)&(0x3);
					replace_one_nc(read->read_bits,read->readLen,(*last_correct)-1-l,nt);
				}
					



			}
		
			
	
		}
	}
	return 1;
}

bool Sparse_Denoising(struct read_t *read,struct hashtable *ht,struct hashtable2 *ht2,size_t *CovTh,size_t *CorrTh,int K_size,int gap,uint64_t *correction_cnt,bool Hybrid)
{
//	bool DISPLAY_SPLIT=0;
	bool Seq_Pad=0;
	bool CORRECTABLE=1;
	bool DETAILED_CHECK=1;
	bool SELECT_OPTIMAL=0;
	bool CORRECT=1;
	bool ALLOW_WEAK_COV=1;
	bool SEARCH_LEFT=0,SEARCH_RIGHT=0;
					
	if(Hybrid)
	{
		Seq_Pad=1;
	}
	else
	{
		Seq_Pad=1;
	}
	int readLen=read->readLen;
	int OverlapKmers=readLen-K_size+1;

	if(gap>=OverlapKmers)
	{return 0;}

	int Read_arr_sz=readLen/32+1;
	int rem=readLen%32;
	if(rem==0)
	{Read_arr_sz--;}
	int tot_bits=Read_arr_sz*64;
	size_t ht_sz;
	if(K_size<=32)
	{
		ht_sz=ht->ht_sz;
	}
	else
	{
		ht_sz=ht2->ht_sz;
	}
	//int Kmer_arr_sz=K_size/32+1;


	bool flip[1000],found[1000];//error_nt[200];
	size_t hash_idx[1000];
	memset(flip,0,sizeof(flip));
	memset(found,0,sizeof(found));
	memset(read->error_nt,1,sizeof(read->error_nt));

	uint64_t seq[1000],f_seq[1000],hv[1000];
	uint64_t temp_bits[1000];
	struct kmer_t2 seq_t2[1000],f_seq_t2[1000];

	bucket ** bktptr[1000];
	bucket2 ** bktptr_t2[1000];
	char c_str[1000];
	
	int first_correct=-1,last_correct=-1;
	int first_error=-1,second_error=-1;
	//#pragma omp parallel for
	for (int j=0;j<OverlapKmers;j++)
	{
		//uint64_t temp_bitsarr[5];
	
		
		if(K_size<=32)
		{
			get_sub_arr(read->read_bits,read->readLen,j,K_size,&(seq[j]));	
			
			f_seq[j]=get_rev_comp_seq(seq[j],K_size);
			
			flip[j]=0;
			if(seq[j]>f_seq[j])
			{
				uint64_t t=seq[j];
				seq[j]=f_seq[j];
				f_seq[j]=t;
				flip[j]=1;
			}
		}
		else
		{
			
			get_sub_arr(read->read_bits,read->readLen,j,K_size,seq_t2[j].kmer);	
			f_seq_t2[j]=seq_t2[j];
			get_rev_comp_seq_arr(f_seq_t2[j].kmer,K_size,2);
			flip[j]=0;
			if(uint64_t_cmp(seq_t2[j].kmer,f_seq_t2[j].kmer,2)>0)
			{
				kmer_t2 t=seq_t2[j];
				seq_t2[j]=f_seq_t2[j];
				f_seq_t2[j]=t;
				flip[j]=1;
			}
			
		}
		
		if(K_size<=32)
		{
		hv[j]=MurmurHash64A(&seq[j],sizeof(seq[j]),0);
		
		hash_idx[j]=(size_t) (hv[j]%ht_sz);
		}
		else
		{
			if(K_size<=64)
			{
				hv[j]=MurmurHash64A(&seq_t2[j],sizeof(seq_t2[j]),0);
		
				hash_idx[j]=(size_t) (hv[j]%ht_sz);
			}
		}
		struct bucket ** ptr;
		struct bucket2 ** ptr_t2;
		if(K_size<=32)
		{
			bktptr[j]= &(ht->store_pos[hash_idx[j]]);

			found[j]=look_up_in_a_list(seq[j],&bktptr[j]);
		}
		else
		{	
			bktptr_t2[j]= &(ht2->store_pos[hash_idx[j]]);
			found[j]=look_up_in_a_list2(&(seq_t2[j]),&bktptr_t2[j]);
		}
		if(found[j])
		{
			if(K_size<=32)
			{
				if((*bktptr[j])->kmer_info.cov1>=*CovTh)
				{
					memset(&(read->error_nt[j+1]),0,(K_size-1)*sizeof(bool));
				}
		
			}
			else
			{
				if((*bktptr_t2[j])->kmer_info.cov1>=*CovTh)
				{
					//memset(&(error_nt[j]),0,K_size*sizeof(bool));
					memset(&(read->error_nt[j+1]),0,(K_size-1)*sizeof(bool));
				}
			}
		}

	}
	for (int i=0;i<OverlapKmers;++i)
	{
		if(read->error_nt[i]==0)
		{
			first_correct=i-1;
			read->error_nt[i-1]=0;
			break;
		}
	}
	for (int i=max(gap-1,1);i<OverlapKmers;++i)//-g
	{
		if(read->error_nt[i]==1&&read->error_nt[i-1]==0)
		{
			last_correct=i-K_size;
			SEARCH_RIGHT=1;
			break;
		}
	}
	if(first_correct>=gap)
	{SEARCH_LEFT=1;}

	
	int l_success=0,r_success=0;
	if(SEARCH_LEFT==1)
	{
		while(first_correct>=gap)
		{
			bool SPLIT=0;
			if(flip[first_correct]==0)
			{
				if(K_size<=32)
					SPLIT=(*bktptr[first_correct])->kmer_info.split_left;
				else
					SPLIT=(*bktptr_t2[first_correct])->kmer_info.split_left;
			

				if(!SPLIT)
				{
					struct edge_node *opt_edge;
					if(K_size<=32)
						opt_edge=(*bktptr[first_correct])->kmer_info.left;
					else
						opt_edge=(*bktptr_t2[first_correct])->kmer_info.left;

					if(opt_edge==NULL)
					{
						break;
					}	
					int edge_len=opt_edge->len+1;
					uint64_t edge_bits=opt_edge->edge;

					if(flip[first_correct])
					{
						edge_bits=get_rev_comp_seq(edge_bits,edge_len);
					}
				
					int mod_cnt=0;
					for(int l=0;l<edge_len;++l)
					{
				
						uint64_t read_bits;
						get_sub_arr(read->read_bits,read->readLen,first_correct-edge_len,edge_len,&read_bits);
					
						uint64_t mask=0x3;
						mask<<=2*l;
					
						uint64_t nt;
						if((read_bits&mask)!=(edge_bits&mask))
						{
							nt=(edge_bits>>2*l)&(0x3);
							replace_one_nc(read->read_bits,read->readLen,(first_correct)-1-l,nt);
							(*correction_cnt)++;
							mod_cnt++;
							if(mod_cnt>3)
							{
								(*correction_cnt)-=mod_cnt;
								return 0;

							}
						}
					
						
					}
					memset((&read->error_nt[(first_correct)-edge_len]),0,edge_len*sizeof(bool));

					
				}


			}
			else
			{
				if(K_size<=32)
				SPLIT=(*bktptr[first_correct])->kmer_info.split_right;
				else
					SPLIT=(*bktptr_t2[first_correct])->kmer_info.split_right;
		

				if(!SPLIT)
				{
					struct edge_node *opt_edge;
					if(K_size<=32)
					opt_edge=(*bktptr[first_correct])->kmer_info.right;
					else
					opt_edge=(*bktptr_t2[first_correct])->kmer_info.right;

					if(opt_edge==NULL)
					{
						break;
					}	
					
					int edge_len=opt_edge->len+1;
					uint64_t edge_bits=opt_edge->edge;

					if(flip[first_correct])
					{
						edge_bits=get_rev_comp_seq(edge_bits,edge_len);
					}
					int mod_cnt=0;
					for(int l=0;l<edge_len;++l)
					{
						
						uint64_t read_bits;
						get_sub_arr(read->read_bits,read->readLen,first_correct-edge_len,edge_len,&read_bits);
					
						uint64_t mask=0x3;
						mask<<=2*l;
					
						uint64_t nt;
						if((read_bits&mask)!=(edge_bits&mask))
						{
							nt=(edge_bits>>2*l)&(0x3);
							replace_one_nc(read->read_bits,read->readLen,(first_correct)-1-l,nt);
							(*correction_cnt)++;
							mod_cnt++;
							if(mod_cnt>3)
							{
								(*correction_cnt)-=mod_cnt;
								return 0;

							}
						}
					
					


					}
					memset((&read->error_nt[(first_correct)-edge_len]),0,edge_len*sizeof(bool));

					
				}


			}


			if(SPLIT)
			{
				l_success=BFsearch_denoising(read,ht,ht2,CovTh,CorrTh,K_size,gap,&first_correct,0);
				if(l_success>0)
				{(*correction_cnt)++;}
			}
			

			if((!ALLOW_WEAK_COV)&&(!SPLIT))
			{
				l_success=-1;
			}
			if(ALLOW_WEAK_COV&&(!SPLIT))
			{
				l_success=1;
			}
			//update left
			if(l_success<0)
			{break;}

			if(flip[first_correct]==1)
			{
				if(K_size<=32)
					SPLIT=(*bktptr[first_correct])->kmer_info.split_right;
				else
					SPLIT=(*bktptr_t2[first_correct])->kmer_info.split_right;

			}
			else
			{
				if(K_size<=32)
					SPLIT=(*bktptr[first_correct])->kmer_info.split_left;
				else
					SPLIT=(*bktptr_t2[first_correct])->kmer_info.split_left;
			}

			for (int j=first_correct-1;j>=0;j--)
			{
				
				if(K_size<=32)
				{
					
					get_sub_arr(read->read_bits,read->readLen,j,K_size,&(seq[j]));
					f_seq[j]=get_rev_comp_seq(seq[j],K_size);
					flip[j]=0;
					if(seq[j]>f_seq[j])
					{
						uint64_t t=seq[j];
						seq[j]=f_seq[j];
						f_seq[j]=t;
						flip[j]=1;
					}
				}
				else
				{
					get_sub_arr(read->read_bits,read->readLen,j,K_size,seq_t2[j].kmer);
					f_seq_t2[j]=seq_t2[j];
					get_rev_comp_seq_arr(f_seq_t2[j].kmer,K_size,2);
					flip[j]=0;
					if(uint64_t_cmp(seq_t2[j].kmer,f_seq_t2[j].kmer,2)>0)
					{
						kmer_t2 t=seq_t2[j];
						seq_t2[j]=f_seq_t2[j];
						f_seq_t2[j]=t;
						flip[j]=1;
					}
			
				}
				
				
				if(K_size<=32)
				{
				hv[j]=MurmurHash64A(&seq[j],sizeof(seq[j]),0);
		
				hash_idx[j]=(size_t) (hv[j]%ht_sz);
				}
				else
				{
					if(K_size<=64)
					{
						hv[j]=MurmurHash64A(&seq_t2[j],sizeof(seq_t2[j]),0);
		
						hash_idx[j]=(size_t) (hv[j]%ht_sz);
					}
				}
		
				

				struct bucket ** ptr;
				struct bucket2 ** ptr_t2;
				if(K_size<=32)
				{
				bktptr[j]= &(ht->store_pos[hash_idx[j]]);

				found[j]=look_up_in_a_list(seq[j],&bktptr[j]);
				}
				else
				{	
					bktptr_t2[j]= &(ht2->store_pos[hash_idx[j]]);

					found[j]=look_up_in_a_list2(&(seq_t2[j]),&bktptr_t2[j]);
			
				}


				if(found[j])
				{
					if(K_size<=32)
					{
						if((*bktptr[j])->kmer_info.cov1>=*CovTh)
						{
							memset(&(read->error_nt[j]),0,K_size*sizeof(bool));
		
						}
					}
					else
					{
					
						if((*bktptr_t2[j])->kmer_info.cov1>=*CovTh)
						{
						
							memset(&(read->error_nt[j]),0,K_size*sizeof(bool));
		
						}
					}
					if(flip[j]==1)
					{
						if(K_size<=32)
						SPLIT=(*bktptr[j])->kmer_info.split_right;
						else
							SPLIT=(*bktptr_t2[j])->kmer_info.split_right;
					}
					else
					{
						if(K_size<=32)
						SPLIT=(*bktptr[j])->kmer_info.split_left;
						else
							SPLIT=(*bktptr_t2[j])->kmer_info.split_left;
					}

				}
				
			}
			bool updated=0;
			for (int i=first_correct-1;i>=0;--i)
			{
				if(read->error_nt[i]==1)
				{
					if(first_correct==i+1&&(updated==0))
					{updated=0;break;}
					first_correct=i+1;
					if(found[first_correct]==0)
					{updated=0;break;}
					SEARCH_LEFT=1;
					updated=1;
					break;
				}
				else
				{
					if(found[i]==1)
					{first_correct=i;updated=1;
					//break;
					}
				}
			}
			if(updated==0)
			{break;}

		}
	}

	if(SEARCH_RIGHT==1)
	{
		while(last_correct<OverlapKmers-gap)//-g
		{
			//cout<<last_correct<<endl;

			bool SPLIT=0;

			if(flip[last_correct]==1)
			{
				if(K_size<=32)
					SPLIT=(*bktptr[last_correct])->kmer_info.split_left;
				else
					SPLIT=(*bktptr_t2[last_correct])->kmer_info.split_left;
				//bug here! Limit the # of modifications
				if(!SPLIT)
				{
					struct edge_node *opt_edge;

					if(K_size<=32)
					{
						opt_edge=(*bktptr[last_correct])->kmer_info.left;
					}
					else
					{
						opt_edge=(*bktptr_t2[last_correct])->kmer_info.left;
					}
					if(opt_edge==NULL)
					{
						break;
					}
					int edge_len=opt_edge->len+1;
					uint64_t edge_bits=opt_edge->edge;
				
					if(flip[last_correct])
					{
						edge_bits=get_rev_comp_seq(edge_bits,edge_len);
					}

					int mod_cnt=0;
					uint64_t read_bits;
					get_sub_arr(read->read_bits,read->readLen,last_correct+K_size,edge_len,&read_bits);
					
					for(int l=0;l<edge_len;++l)
					{
				
						
						uint64_t mask=0x3;
						mask<<=2*l;
					
						uint64_t nt;
						if((read_bits&mask)!=(edge_bits&mask))
						{
							nt=(edge_bits>>2*l)&(0x3);
							replace_one_nc(read->read_bits,read->readLen,(edge_len-1-l)+(last_correct)+K_size,nt);
							(*correction_cnt)++;
							mod_cnt++;
							if(mod_cnt>3)
							{
								(*correction_cnt)-=mod_cnt;
								return 0;

							}
						}
					
//						read->error_nt[(edge_len-1-l)+(last_correct)+K_size]=0;
					}
					memset((&read->error_nt[(last_correct)+K_size]),0,edge_len*sizeof(bool));

					
				}

			}
			else
			{
				if(K_size<=32)
				SPLIT=(*bktptr[last_correct])->kmer_info.split_right;
				else
					SPLIT=(*bktptr_t2[last_correct])->kmer_info.split_right;
			//bug here! Limit the # of modifications
				
				if(!SPLIT)
				{
					struct edge_node *opt_edge;
					if(K_size<=32)
						opt_edge=(*bktptr[last_correct])->kmer_info.right;
					else
						opt_edge=(*bktptr_t2[last_correct])->kmer_info.right;

					if(opt_edge==NULL)
					{
						break;
					}	

					int edge_len=opt_edge->len+1;
					uint64_t edge_bits=opt_edge->edge;
				
					if(flip[last_correct])
					{
						edge_bits=get_rev_comp_seq(edge_bits,edge_len);
					}
				

					int mod_cnt=0;
					uint64_t read_bits;
					get_sub_arr(read->read_bits,read->readLen,last_correct+K_size,edge_len,&read_bits);
					
					for(int l=0;l<edge_len;++l)
					{
				
						
						uint64_t mask=0x3;
						mask<<=2*l;
					
						uint64_t nt;
						if((read_bits&mask)!=(edge_bits&mask))
						{
							nt=(edge_bits>>2*l)&(0x3);
							replace_one_nc(read->read_bits,read->readLen,(edge_len-1-l)+(last_correct)+K_size,nt);
							(*correction_cnt)++;
							mod_cnt++;
							if(mod_cnt>3)
							{
								(*correction_cnt)-=mod_cnt;
								return 0;

							}
						}
					
					}
					memset((&read->error_nt[(last_correct)+K_size]),0,edge_len*sizeof(bool));


					
					
				}


			}

			


			if(SPLIT)
			{
				r_success=BFsearch_denoising(read,ht,ht2,CovTh,CorrTh,K_size,gap,&last_correct,1);
				if(r_success>0)
				{(*correction_cnt)++;}
			}
			if((!ALLOW_WEAK_COV)&&(!SPLIT))
			{
				r_success=-1;
			}
			if(ALLOW_WEAK_COV&&(!SPLIT))
			{
				r_success=1;
			}


			if(r_success<0)
			{break;}
			if(r_success==1)
			{
				if(K_size<=32)
				{
					if(flip[last_correct]==0)
					{
						SPLIT=(*bktptr[last_correct])->kmer_info.split_right;
					}
					else
					{
						SPLIT=(*bktptr[last_correct])->kmer_info.split_left;
					}
				}
				else
				{
					if(flip[last_correct]==0)
					{
						SPLIT=(*bktptr_t2[last_correct])->kmer_info.split_right;
					}
					else
					{
						SPLIT=(*bktptr_t2[last_correct])->kmer_info.split_left;
					}
				
				}
				for (int j=last_correct+1;j<OverlapKmers;j++)
				{

					
				
					
					if(K_size<=32)
					{
						get_sub_arr(read->read_bits,read->readLen,j,K_size,&(seq[j]));
					
					//	char c_str[300];
						//bitsarr2str(&seq[j],K_size,c_str,1);
						
						f_seq[j]=get_rev_comp_seq(seq[j],K_size);
						flip[j]=0;
						if(seq[j]>f_seq[j])
						{
							uint64_t t=seq[j];
							seq[j]=f_seq[j];
							f_seq[j]=t;
							flip[j]=1;
						}
					}
					else
					{
						get_sub_arr(read->read_bits,read->readLen,j,K_size,seq_t2[j].kmer);
						f_seq_t2[j]=seq_t2[j];
						get_rev_comp_seq_arr(f_seq_t2[j].kmer,K_size,2);
						flip[j]=0;
						if(uint64_t_cmp(seq_t2[j].kmer,f_seq_t2[j].kmer,2)>0)
						{
							kmer_t2 t=seq_t2[j];
							seq_t2[j]=f_seq_t2[j];
							f_seq_t2[j]=t;
							flip[j]=1;
						}
			
					}

		
					if(K_size<=32)
					{
					hv[j]=MurmurHash64A(&seq[j],sizeof(seq[j]),0);
		
					hash_idx[j]=(size_t) (hv[j]%ht_sz);
					}
					else
					{
						if(K_size<=64)
						{
							hv[j]=MurmurHash64A(&seq_t2[j],sizeof(seq_t2[j]),0);
		
							hash_idx[j]=(size_t) (hv[j]%ht_sz);
						}
					}

					struct bucket ** ptr;
					struct bucket2 ** ptr_t2;
					if(K_size<=32)
					{
					bktptr[j]= &(ht->store_pos[hash_idx[j]]);

					found[j]=look_up_in_a_list(seq[j],&bktptr[j]);
					}
					else
					{	
						bktptr_t2[j]= &(ht2->store_pos[hash_idx[j]]);

						found[j]=look_up_in_a_list2(&(seq_t2[j]),&bktptr_t2[j]);
			
					}
				
					
					if(found[j])
					{
						if(K_size<=32)
						{
							if((*bktptr[j])->kmer_info.cov1>=*CovTh)
							{
								
								//memset(&(error_nt[j]),0,K_size*sizeof(bool));
								memset(&(read->error_nt[j+1]),0,(K_size-1)*sizeof(bool));
		
							}
						}
						else
						{
					
							if((*bktptr_t2[j])->kmer_info.cov1>=*CovTh)
							{
							
//								memset(&(error_nt[j]),0,K_size*sizeof(bool));
								memset(&(read->error_nt[j+1]),0,(K_size-1)*sizeof(bool));
		
							}
						}
						if(flip[j]==0)
						{
							if(K_size<=32)
							SPLIT=(*bktptr[j])->kmer_info.split_right;
							else
							{SPLIT=(*bktptr_t2[j])->kmer_info.split_right;}
						}
						else
						{
							if(K_size<=32)
							SPLIT=(*bktptr[j])->kmer_info.split_left;
							else
							{SPLIT=(*bktptr_t2[j])->kmer_info.split_left;}
						}


						//break;//added on 0624
					}

				}
				bool updated=0;
				for (int i=last_correct+1;i<OverlapKmers;++i)
				{
					

					if(read->error_nt[i+K_size]==1&&read->error_nt[i+K_size-1]==0)
					{
						if(last_correct>=i-1&&updated==0)
						{updated=0;break;}
						last_correct=i;
						if(found[last_correct]==0)
						{
	//						cout<<"Last correct not found."<<endl;
							updated=0;break;
						}
						SEARCH_RIGHT=1;
						updated=1;
//						cout<<"Last Found."<<endl;
						break;
					}
					else
					{
						if(found[i]==1)
						{
							last_correct=i;updated=1;
						//break;
						}//cout<<"Normal Found"<<endl;}
					}
				}


				if(updated==0)
				{break;}


			}
		}
	}


	if(first_correct==-1)
	{
		CORRECTABLE=0;
		return 0;
	}
	bool start_print=0;

	
	if(Seq_Pad)
	{

		last_correct=-1;
		for(int i=0;i<read->readLen;++i)
		{
			if((read->error_nt[i]==0&&start_print==0&&found[i]==1))
			{
				start_print=1;
				
				first_correct=i;
			}
			if((start_print&&read->error_nt[i]))
			{
				for(int ii=i;ii>i-gap;--ii)
				{
					if(found[ii-K_size])
					{
						last_correct=ii-K_size;
						break;
					}
				}
				if(last_correct>=0)
				{break;}
			}
		}
		
		
		if(first_correct>0)
		{
			bool Padding=0;
			if(K_size<=32)
			{
				if((*bktptr[first_correct])->kmer_info.cov1>*CorrTh)
				{Padding=1;}
			}
			else
			{
				if((*bktptr_t2[first_correct])->kmer_info.cov1>*CorrTh)
				{Padding=1;}
			}
			
			

			bool SPLIT=0;
			
			if(Padding)
			{
				
				int correction_cnt=0;
				edge_node *opt_edge;
				int opt_pad_len=0;
				
				for (int round=1;round<=1;++round)
				{
					edge_node *edge_ptr;
					int edge_len,pad_len=0;
					if(flip[first_correct]==1)
					{
						if(K_size<=32)
						{
							edge_ptr=(*bktptr[first_correct])->kmer_info.right;
						}
						else
						{
							edge_ptr=(*bktptr_t2[first_correct])->kmer_info.right;					
						}

					}
					else
					{
						if(K_size<=32)
						{
							edge_ptr=(*bktptr[first_correct])->kmer_info.left;
						}
						else
						{
							edge_ptr=(*bktptr_t2[first_correct])->kmer_info.left;					
						}
			
					}

					int opt_edge_len=0;
					uint64_t opt_edge_bits=0;

					while(edge_ptr!=NULL)
					{
						if((edge_ptr->edge_cov<*CovTh)&&(edge_ptr->edge_cov<0xff))
						{
							edge_ptr=edge_ptr->nxt_edge;
							continue;
						}
						edge_len=edge_ptr->len+1;
						
						if(first_correct<edge_len)
						{
							pad_len=first_correct;
						}
						else
						{
							pad_len=edge_len;
						}

						uint64_t read_bits,edge_bits;
				
						get_sub_arr(read->read_bits,read->readLen,first_correct-pad_len,pad_len,&read_bits);
						edge_bits=edge_ptr->edge;
						if(flip[first_correct])
						{
							edge_bits=get_rev_comp_seq(edge_bits,edge_len);
						}
						int mod_cnt=0;
						int last_mod=-1;
						for(int l=0;l<pad_len;++l)
						{
							uint64_t mask=0x3;
							mask<<=2*l;
					
							if((read_bits&mask)!=(edge_bits&mask))
							{
								mod_cnt++;
								last_mod=l;
							}
							if((l==pad_len-1)&&(last_mod==-1))
							{
								last_mod=l;
							}
						}

					
				
						if(mod_cnt>round)
						{
							edge_ptr=edge_ptr->nxt_edge;
							continue;
						}


						opt_edge=edge_ptr;
						opt_pad_len=pad_len;

						if(last_mod>=0&&DETAILED_CHECK)
						{
							uint64_t edge_bits_t=edge_bits<<(64-2*(1+last_mod));
							if((opt_edge_len!=last_mod+1)||(edge_bits_t!=opt_edge_bits))
							{
								correction_cnt++;
								opt_edge_len=last_mod+1;
								opt_edge_bits=edge_bits_t;

							}
						}
						else
						{

						correction_cnt++;
						}
						edge_ptr=edge_ptr->nxt_edge;
					}
	///
					if(correction_cnt>=1)
					{
						break;
					}


				}
				if(correction_cnt==1)
				{
					int edge_len=opt_edge->len+1;
					uint64_t edge_bits=opt_edge->edge;
					if(flip[first_correct]==1)
					{
						edge_bits=get_rev_comp_seq(edge_bits,edge_len);
					}

					uint64_t read_bits;
					get_sub_arr(read->read_bits,read->readLen,first_correct-opt_pad_len,opt_pad_len,&read_bits);
					
					for(int l=0;l<opt_pad_len;++l)
					{
						uint64_t mask=0x3;
						mask<<=2*l;
					
						uint64_t nt;
						if((read_bits&mask)!=(edge_bits&mask))
						{
							nt=(edge_bits>>2*l)&(0x3);
							replace_one_nc(read->read_bits,read->readLen,(first_correct)-1-l,nt);

						}
						read->error_nt[first_correct-1-l]=0;
					



					}
					
				}

			}

		}







		if(last_correct<OverlapKmers-1&&last_correct>0)
		{
			bool Padding=0;
			if(K_size<=32)
			{
				if((*bktptr[last_correct])->kmer_info.cov1>*CorrTh)
				{Padding=1;}
			}
			else
			{
				if((*bktptr_t2[last_correct])->kmer_info.cov1>*CorrTh)
				{Padding=1;}
			}
			
			

			bool SPLIT=0;
			
			if(Padding)
			{
				
				int correction_cnt=0;
				edge_node *opt_edge;
				int opt_pad_len=0;
				
				for (int round=1;round<=1;++round)
				{
					edge_node *edge_ptr;
					int edge_len,pad_len=0;
					if(flip[last_correct]==0)
					{
						if(K_size<=32)
						{
							edge_ptr=(*bktptr[last_correct])->kmer_info.right;
						}
						else
						{
							edge_ptr=(*bktptr_t2[last_correct])->kmer_info.right;					
						}

					}
					else
					{
						if(K_size<=32)
						{
							edge_ptr=(*bktptr[last_correct])->kmer_info.left;
						}
						else
						{
							edge_ptr=(*bktptr_t2[last_correct])->kmer_info.left;					
						}
			
					}

					int opt_edge_len=0;
					uint64_t opt_edge_bits=0;

					while(edge_ptr!=NULL)
					{
						if((edge_ptr->edge_cov<*CovTh)&&(edge_ptr->edge_cov<0xff))
						{
							edge_ptr=edge_ptr->nxt_edge;
							continue;
						}
						edge_len=edge_ptr->len+1;
						
						if(edge_len<OverlapKmers-1-last_correct)
						{
							pad_len=edge_len;
						}
						else
						{
							pad_len=OverlapKmers-1-last_correct;
						}

						uint64_t read_bits,edge_bits;
				
						get_sub_arr(read->read_bits,read->readLen,last_correct+K_size,pad_len,&read_bits);
						edge_bits=edge_ptr->edge;
						if(flip[last_correct])
						{
							edge_bits=get_rev_comp_seq(edge_bits,edge_len);
						}
						edge_bits>>=2*(edge_len-pad_len);
						int mod_cnt=0;
						int last_mod=-1;
						for(int l=0;l<pad_len;++l)
						{
							uint64_t mask=0x3;
							mask<<=2*l;
					
							if((read_bits&mask)!=(edge_bits&mask))
							{
								mod_cnt++;
								if(last_mod<0)
								{last_mod=l;}
							}
							if(last_mod==-1&&l==pad_len-1)
							{
								last_mod=0;
							}
						}

						
						if(mod_cnt>round)
						{
							edge_ptr=edge_ptr->nxt_edge;
							continue;
						}


						opt_edge=edge_ptr;
						opt_pad_len=pad_len;

						if(last_mod>=0&&DETAILED_CHECK)
						{
							uint64_t edge_bits_t=edge_bits>>(2*(last_mod));
							if((opt_edge_len!=pad_len-last_mod)||(edge_bits_t!=opt_edge_bits))
							{
								correction_cnt++;
								opt_edge_len=pad_len-last_mod;
								opt_edge_bits=edge_bits_t;

							}
						}
						else
						{
							correction_cnt++;
						}
					

						edge_ptr=edge_ptr->nxt_edge;
					}
	///
					if(correction_cnt>=1)
					{	
						break;
					}


				}
				if(correction_cnt==1)
				{
					int edge_len=opt_edge->len+1;
					uint64_t edge_bits=opt_edge->edge;
					if(flip[last_correct]==1)
					{
						edge_bits=get_rev_comp_seq(edge_bits,edge_len);
					}
					edge_bits>>=2*(edge_len-opt_pad_len);

					uint64_t read_bits;
					get_sub_arr(read->read_bits,read->readLen,last_correct+K_size,opt_pad_len,&read_bits);
					
					for(int l=0;l<opt_pad_len;++l)
					{
						uint64_t mask=0x3;
						mask<<=2*l;
					
						uint64_t nt;
						if((read_bits&mask)!=(edge_bits&mask))
						{
							nt=(edge_bits>>2*l)&(0x3);
							replace_one_nc(read->read_bits,read->readLen,last_correct+K_size+(opt_pad_len-1-l),nt);

						}
						read->error_nt[last_correct+K_size+(opt_pad_len-1-l)]=0;
					



					}
					
				}

			}

		}
		
	

	}
	//memcpy(read->error_nt,error_nt,(read->readLen)*sizeof(bool));

	bitsarr2str(read->read_bits,read->readLen,read->c_seq,Read_arr_sz);
	
	/*
	for(int i=0;i<read->readLen;++i)
	{
		if((error_nt[i]==0&&start_print==0)||Hybrid)
		{
			start_print=1;
		}
		if(start_print&&error_nt[i]&&(!Hybrid))
		{break;}

		if(start_print)
		{
			o_dn<<read->c_seq[i];	
		
		}

		
	}
	o_dn<<endl;
	*/
	return 1;

}




#endif
