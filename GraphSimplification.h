#ifndef __GRAPH_SIMPLIFICATION_H
#define __GRAPH_SIMPLIFICATION_H
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


// remove edges and nodes that are low-coverage, round2


//break the bubble links in the bfs.

void BreakLinks( map<struct bucket* ,int > &stacked_nodes, struct bucket* bktptr1, struct bucket* bktptr2,int K_size,int edge_len)
{
	if(stacked_nodes[bktptr1]>1)
	{
		struct edge_node * edge_ptr=bktptr1->kmer_info.right;
		edge_node ** edge_p2p=&(bktptr1->kmer_info.right);
		while(edge_ptr!=NULL)
		{
			if(edge_ptr->len==edge_len-1)
			{
				uint64_t kmer;
				kmer=(bktptr1)->kmer_t.kmer;

				for(int g=edge_len-1;g>=0;--g)
				{
					uint64_t b=(edge_ptr->edge>>2*g)&0x3;
					kmer&=(~((uint64_t)0x3<<(2*(K_size-1))));
					kmer<<=2;
					kmer|=b;

				}
				uint64_t f_kmer=kmer;
				f_kmer=get_rev_comp_seq(kmer,K_size);
				if(kmer>f_kmer)
				{
					kmer=f_kmer;
				}

				if(kmer==bktptr2->kmer_t.kmer)
				{
					(*edge_p2p)=(*edge_p2p)->nxt_edge;
					struct edge_node* f_edge_ptr=edge_ptr;
					edge_ptr=edge_ptr->nxt_edge;
					free(f_edge_ptr);
					stacked_nodes[bktptr1]--;
					break;

				}
			}


			(edge_p2p)=&((*edge_p2p)->nxt_edge);
			edge_ptr=edge_ptr->nxt_edge;

		}


	}


	if(stacked_nodes[bktptr1]<-1)
	{
		struct edge_node * edge_ptr=bktptr1->kmer_info.left;
		edge_node ** edge_p2p=&(bktptr1->kmer_info.left);
		while(edge_ptr!=NULL)
		{
			if(edge_ptr->len==edge_len-1)
			{
				uint64_t kmer;
				kmer=(bktptr1)->kmer_t.kmer;

				for(int g=0;g<edge_len;++g)
				{
					uint64_t b=(edge_ptr->edge>>2*g)&0x3;
					kmer>>=2;
					kmer|=(b<<(2*(K_size-1)));
				}
				uint64_t f_kmer=kmer;
				f_kmer=get_rev_comp_seq(kmer,K_size);
				if(kmer>f_kmer)
				{
					kmer=f_kmer;
				}

				if(kmer==bktptr2->kmer_t.kmer)
				{
					(*edge_p2p)=(*edge_p2p)->nxt_edge;
					struct edge_node* f_edge_ptr=edge_ptr;
					edge_ptr=edge_ptr->nxt_edge;
					free(f_edge_ptr);
					stacked_nodes[bktptr1]++;
					break;

				}
			}
			(edge_p2p)=&((*edge_p2p)->nxt_edge);
			edge_ptr=edge_ptr->nxt_edge;

		}


	}




	if(stacked_nodes[bktptr2]<0)
	{
		struct edge_node * edge_ptr=bktptr2->kmer_info.right;
		edge_node ** edge_p2p=&(bktptr2->kmer_info.right);
		while(edge_ptr!=NULL)
		{
			if(edge_ptr->len==edge_len-1)
			{
				uint64_t kmer;
				kmer=(bktptr2)->kmer_t.kmer;

				for(int g=edge_len-1;g>=0;--g)
				{
					uint64_t b=(edge_ptr->edge>>2*g)&0x3;
					kmer&=(~((uint64_t)0x3<<(2*(K_size-1))));
					kmer<<=2;
					kmer|=b;

				}
				uint64_t f_kmer=kmer;
				f_kmer=get_rev_comp_seq(kmer,K_size);
				if(kmer>f_kmer)
				{
					kmer=f_kmer;
				}

				if(kmer==bktptr1->kmer_t.kmer)
				{
					(*edge_p2p)=(*edge_p2p)->nxt_edge;
					struct edge_node* f_edge_ptr=edge_ptr;
					edge_ptr=edge_ptr->nxt_edge;
					free(f_edge_ptr);

					break;

				}
			}


			(edge_p2p)=&((*edge_p2p)->nxt_edge);
			edge_ptr=edge_ptr->nxt_edge;

		}


	}


	if(stacked_nodes[bktptr2]>0)
	{
		struct edge_node * edge_ptr=bktptr2->kmer_info.left;
		edge_node ** edge_p2p=&(bktptr2->kmer_info.left);
		while(edge_ptr!=NULL)
		{
			if(edge_ptr->len==edge_len-1)
			{
				uint64_t kmer;
				kmer=(bktptr2)->kmer_t.kmer;

				for(int g=0;g<edge_len;++g)
				{
					uint64_t b=(edge_ptr->edge>>2*g)&0x3;
					kmer>>=2;
					kmer|=(b<<(2*(K_size-1)));
				}
				uint64_t f_kmer=kmer;
				f_kmer=get_rev_comp_seq(kmer,K_size);
				if(kmer>f_kmer)
				{
					kmer=f_kmer;
				}

				if(kmer==bktptr1->kmer_t.kmer)
				{
					(*edge_p2p)=(*edge_p2p)->nxt_edge;
					struct edge_node* f_edge_ptr=edge_ptr;
					edge_ptr=edge_ptr->nxt_edge;
					free(f_edge_ptr);

					break;

				}
			}
			(edge_p2p)=&((*edge_p2p)->nxt_edge);
			edge_ptr=edge_ptr->nxt_edge;

		}


	}

}

void BreakLinks2( map<struct bucket2* ,int > &stacked_nodes, struct bucket2* bktptr1, struct bucket2* bktptr2,int K_size,int edge_len)
{
	if(stacked_nodes[bktptr1]>1)
	{
		struct edge_node * edge_ptr=bktptr1->kmer_info.right;
		edge_node ** edge_p2p=&(bktptr1->kmer_info.right);
		while(edge_ptr!=NULL)
		{
			if(edge_ptr->len==edge_len-1)
			{
				kmer_t2 kmer;
				kmer=(bktptr1)->kmer_t2;

				for(int g=edge_len-1;g>=0;--g)
				{
					uint64_t b=(edge_ptr->edge>>2*g)&0x3;
					kmer.kmer[0]&=(~((uint64_t)0x3<<(2*(K_size-1-32))));
					L_shift_NB(kmer.kmer,2,2);
					kmer.kmer[1]|=b;

				}
				kmer_t2 f_kmer=kmer;
				get_rev_comp_seq_arr(f_kmer.kmer,K_size,2);
				if(uint64_t_cmp(kmer.kmer,f_kmer.kmer,2)>0)
				{
					kmer=f_kmer;
				}

				if(uint64_t_cmp(kmer.kmer,bktptr2->kmer_t2.kmer,2)==0)
				{
					(*edge_p2p)=(*edge_p2p)->nxt_edge;
					struct edge_node* f_edge_ptr=edge_ptr;
					edge_ptr=edge_ptr->nxt_edge;
					free(f_edge_ptr);
					stacked_nodes[bktptr1]--;
					break;

				}
			}


			(edge_p2p)=&((*edge_p2p)->nxt_edge);
			edge_ptr=edge_ptr->nxt_edge;

		}


	}


	if(stacked_nodes[bktptr1]<-1)
	{
		struct edge_node * edge_ptr=bktptr1->kmer_info.left;
		edge_node ** edge_p2p=&(bktptr1->kmer_info.left);
		while(edge_ptr!=NULL)
		{
			if(edge_ptr->len==edge_len-1)
			{
				kmer_t2 kmer;
				kmer=(bktptr1)->kmer_t2;

				for(int g=0;g<edge_len;++g)
				{
					uint64_t b=(edge_ptr->edge>>2*g)&0x3;
					R_shift_NB(kmer.kmer,2,2);
					kmer.kmer[0]|=(b<<(2*(K_size-1-32)));
				}
				kmer_t2 f_kmer=kmer;
				get_rev_comp_seq_arr(f_kmer.kmer,K_size,2);
				if(uint64_t_cmp(kmer.kmer,f_kmer.kmer,2)>0)
				{
					kmer=f_kmer;
				}

				if(uint64_t_cmp(kmer.kmer,bktptr2->kmer_t2.kmer,2)==0)
				{
					(*edge_p2p)=(*edge_p2p)->nxt_edge;
					struct edge_node* f_edge_ptr=edge_ptr;
					edge_ptr=edge_ptr->nxt_edge;
					free(f_edge_ptr);
					stacked_nodes[bktptr1]++;
					break;

				}
			}
			(edge_p2p)=&((*edge_p2p)->nxt_edge);
			edge_ptr=edge_ptr->nxt_edge;

		}


	}




	if(stacked_nodes[bktptr2]<0)
	{
		struct edge_node * edge_ptr=bktptr2->kmer_info.right;
		edge_node ** edge_p2p=&(bktptr2->kmer_info.right);
		while(edge_ptr!=NULL)
		{
			if(edge_ptr->len==edge_len-1)
			{
				kmer_t2 kmer;
				kmer=(bktptr2)->kmer_t2;

				for(int g=edge_len-1;g>=0;--g)
				{
					uint64_t b=(edge_ptr->edge>>2*g)&0x3;
					kmer.kmer[0]&=(~((uint64_t)0x3<<(2*(K_size-1-32))));
					L_shift_NB(kmer.kmer,2,2);
					kmer.kmer[1]|=b;

				}
				kmer_t2 f_kmer=kmer;
				get_rev_comp_seq_arr(f_kmer.kmer,K_size,2);
				if(uint64_t_cmp(kmer.kmer,f_kmer.kmer,2)>0)
				{
					kmer=f_kmer;
				}

				if(uint64_t_cmp(kmer.kmer,bktptr1->kmer_t2.kmer,2)==0)
				{
					(*edge_p2p)=(*edge_p2p)->nxt_edge;
					struct edge_node* f_edge_ptr=edge_ptr;
					edge_ptr=edge_ptr->nxt_edge;
					free(f_edge_ptr);

					break;

				}
			}


			(edge_p2p)=&((*edge_p2p)->nxt_edge);
			edge_ptr=edge_ptr->nxt_edge;

		}


	}


	if(stacked_nodes[bktptr2]>0)
	{
		struct edge_node * edge_ptr=bktptr2->kmer_info.left;
		edge_node ** edge_p2p=&(bktptr2->kmer_info.left);
		while(edge_ptr!=NULL)
		{
			if(edge_ptr->len==edge_len-1)
			{
				kmer_t2 kmer;
				kmer=(bktptr2)->kmer_t2;

				for(int g=0;g<edge_len;++g)
				{
					uint64_t b=(edge_ptr->edge>>2*g)&0x3;
					R_shift_NB(kmer.kmer,2,2);
					kmer.kmer[0]|=(b<<(2*(K_size-1-32)));
				}
				kmer_t2 f_kmer=kmer;
				get_rev_comp_seq_arr(f_kmer.kmer,K_size,2);
				if(uint64_t_cmp(kmer.kmer,f_kmer.kmer,2)>0)
				{
					kmer=f_kmer;
				}

				if(uint64_t_cmp(kmer.kmer,bktptr1->kmer_t2.kmer,2)==0)
				{
					(*edge_p2p)=(*edge_p2p)->nxt_edge;
					struct edge_node* f_edge_ptr=edge_ptr;
					edge_ptr=edge_ptr->nxt_edge;
					free(f_edge_ptr);

					break;

				}
			}
			(edge_p2p)=&((*edge_p2p)->nxt_edge);
			edge_ptr=edge_ptr->nxt_edge;

		}


	}

}

void BreakLinks3( map<struct bucket3* ,int > &stacked_nodes, struct bucket3* bktptr1, struct bucket3* bktptr2,int K_size,int edge_len)
{
	if(stacked_nodes[bktptr1]>1)
	{
		struct edge_node * edge_ptr=bktptr1->kmer_info.right;
		edge_node ** edge_p2p=&(bktptr1->kmer_info.right);
		while(edge_ptr!=NULL)
		{
			if(edge_ptr->len==edge_len-1)
			{
				kmer_t3 kmer;
				kmer=(bktptr1)->kmer_t3;

				for(int g=edge_len-1;g>=0;--g)
				{
					uint64_t b=(edge_ptr->edge>>2*g)&0x3;
					kmer.kmer[0]&=(~((uint64_t)0x3<<(2*(K_size-1-64))));
					L_shift_NB(kmer.kmer,2,3);
					kmer.kmer[2]|=b;

				}
				kmer_t3 f_kmer=kmer;
				get_rev_comp_seq_arr(f_kmer.kmer,K_size,3);
				if(uint64_t_cmp(kmer.kmer,f_kmer.kmer,3)>0)
				{
					kmer=f_kmer;
				}

				if(uint64_t_cmp(kmer.kmer,bktptr2->kmer_t3.kmer,3)==0)
				{
					(*edge_p2p)=(*edge_p2p)->nxt_edge;
					struct edge_node* f_edge_ptr=edge_ptr;
					edge_ptr=edge_ptr->nxt_edge;
					free(f_edge_ptr);
					stacked_nodes[bktptr1]--;
					break;

				}
			}


			(edge_p2p)=&((*edge_p2p)->nxt_edge);
			edge_ptr=edge_ptr->nxt_edge;

		}


	}


	if(stacked_nodes[bktptr1]<-1)
	{
		struct edge_node * edge_ptr=bktptr1->kmer_info.left;
		edge_node ** edge_p2p=&(bktptr1->kmer_info.left);
		while(edge_ptr!=NULL)
		{
			if(edge_ptr->len==edge_len-1)
			{
				kmer_t3 kmer;
				kmer=(bktptr1)->kmer_t3;

				for(int g=0;g<edge_len;++g)
				{
					uint64_t b=(edge_ptr->edge>>2*g)&0x3;
					R_shift_NB(kmer.kmer,2,3);
					kmer.kmer[0]|=(b<<(2*(K_size-1-64)));
				}
				kmer_t3 f_kmer=kmer;
				get_rev_comp_seq_arr(f_kmer.kmer,K_size,3);
				if(uint64_t_cmp(kmer.kmer,f_kmer.kmer,3)>0)
				{
					kmer=f_kmer;
				}

				if(uint64_t_cmp(kmer.kmer,bktptr2->kmer_t3.kmer,3)==0)
				{
					(*edge_p2p)=(*edge_p2p)->nxt_edge;
					struct edge_node* f_edge_ptr=edge_ptr;
					edge_ptr=edge_ptr->nxt_edge;
					free(f_edge_ptr);
					stacked_nodes[bktptr1]++;
					break;

				}
			}
			(edge_p2p)=&((*edge_p2p)->nxt_edge);
			edge_ptr=edge_ptr->nxt_edge;

		}


	}




	if(stacked_nodes[bktptr2]<0)
	{
		struct edge_node * edge_ptr=bktptr2->kmer_info.right;
		edge_node ** edge_p2p=&(bktptr2->kmer_info.right);
		while(edge_ptr!=NULL)
		{
			if(edge_ptr->len==edge_len-1)
			{
				kmer_t3 kmer;
				kmer=(bktptr2)->kmer_t3;

				for(int g=edge_len-1;g>=0;--g)
				{
					uint64_t b=(edge_ptr->edge>>2*g)&0x3;
					kmer.kmer[0]&=(~((uint64_t)0x3<<(2*(K_size-1-64))));
					L_shift_NB(kmer.kmer,2,3);
					kmer.kmer[2]|=b;

				}
				kmer_t3 f_kmer=kmer;
				get_rev_comp_seq_arr(f_kmer.kmer,K_size,3);
				if(uint64_t_cmp(kmer.kmer,f_kmer.kmer,3)>0)
				{
					kmer=f_kmer;
				}

				if(uint64_t_cmp(kmer.kmer,bktptr1->kmer_t3.kmer,3)==0)
				{
					(*edge_p2p)=(*edge_p2p)->nxt_edge;
					struct edge_node* f_edge_ptr=edge_ptr;
					edge_ptr=edge_ptr->nxt_edge;
					free(f_edge_ptr);

					break;

				}
			}


			(edge_p2p)=&((*edge_p2p)->nxt_edge);
			edge_ptr=edge_ptr->nxt_edge;

		}


	}


	if(stacked_nodes[bktptr2]>0)
	{
		struct edge_node * edge_ptr=bktptr2->kmer_info.left;
		edge_node ** edge_p2p=&(bktptr2->kmer_info.left);
		while(edge_ptr!=NULL)
		{
			if(edge_ptr->len==edge_len-1)
			{
				kmer_t3 kmer;
				kmer=(bktptr2)->kmer_t3;

				for(int g=0;g<edge_len;++g)
				{
					uint64_t b=(edge_ptr->edge>>2*g)&0x3;
					R_shift_NB(kmer.kmer,2,3);
					kmer.kmer[0]|=(b<<(2*(K_size-1-64)));
				}
				kmer_t3 f_kmer=kmer;
				get_rev_comp_seq_arr(f_kmer.kmer,K_size,3);
				if(uint64_t_cmp(kmer.kmer,f_kmer.kmer,3)>0)
				{
					kmer=f_kmer;
				}

				if(uint64_t_cmp(kmer.kmer,bktptr1->kmer_t3.kmer,3)==0)
				{
					(*edge_p2p)=(*edge_p2p)->nxt_edge;
					struct edge_node* f_edge_ptr=edge_ptr;
					edge_ptr=edge_ptr->nxt_edge;
					free(f_edge_ptr);

					break;

				}
			}
			(edge_p2p)=&((*edge_p2p)->nxt_edge);
			edge_ptr=edge_ptr->nxt_edge;

		}


	}

}


void BreakLinks4( map<struct bucket4* ,int > &stacked_nodes, struct bucket4* bktptr1, struct bucket4* bktptr2,int K_size,int edge_len)
{
	if(stacked_nodes[bktptr1]>1)
	{
		struct edge_node * edge_ptr=bktptr1->kmer_info.right;
		edge_node ** edge_p2p=&(bktptr1->kmer_info.right);
		while(edge_ptr!=NULL)
		{
			if(edge_ptr->len==edge_len-1)
			{
				kmer_t4 kmer;
				kmer=(bktptr1)->kmer_t4;

				for(int g=edge_len-1;g>=0;--g)
				{
					uint64_t b=(edge_ptr->edge>>2*g)&0x3;
					kmer.kmer[0]&=(~((uint64_t)0x3<<(2*(K_size-1-96))));
					L_shift_NB(kmer.kmer,2,4);
					kmer.kmer[3]|=b;

				}
				kmer_t4 f_kmer=kmer;
				get_rev_comp_seq_arr(f_kmer.kmer,K_size,4);
				if(uint64_t_cmp(kmer.kmer,f_kmer.kmer,4)>0)
				{
					kmer=f_kmer;
				}

				if(uint64_t_cmp(kmer.kmer,bktptr2->kmer_t4.kmer,4)==0)
				{
					(*edge_p2p)=(*edge_p2p)->nxt_edge;
					struct edge_node* f_edge_ptr=edge_ptr;
					edge_ptr=edge_ptr->nxt_edge;
					free(f_edge_ptr);
					stacked_nodes[bktptr1]--;
					break;

				}
			}


			(edge_p2p)=&((*edge_p2p)->nxt_edge);
			edge_ptr=edge_ptr->nxt_edge;

		}


	}


	if(stacked_nodes[bktptr1]<-1)
	{
		struct edge_node * edge_ptr=bktptr1->kmer_info.left;
		edge_node ** edge_p2p=&(bktptr1->kmer_info.left);
		while(edge_ptr!=NULL)
		{
			if(edge_ptr->len==edge_len-1)
			{
				kmer_t4 kmer;
				kmer=(bktptr1)->kmer_t4;

				for(int g=0;g<edge_len;++g)
				{
					uint64_t b=(edge_ptr->edge>>2*g)&0x3;
					R_shift_NB(kmer.kmer,2,4);
					kmer.kmer[0]|=(b<<(2*(K_size-1-96)));
				}
				kmer_t4 f_kmer=kmer;
				get_rev_comp_seq_arr(f_kmer.kmer,K_size,4);
				if(uint64_t_cmp(kmer.kmer,f_kmer.kmer,4)>0)
				{
					kmer=f_kmer;
				}

				if(uint64_t_cmp(kmer.kmer,bktptr2->kmer_t4.kmer,4)==0)
				{
					(*edge_p2p)=(*edge_p2p)->nxt_edge;
					struct edge_node* f_edge_ptr=edge_ptr;
					edge_ptr=edge_ptr->nxt_edge;
					free(f_edge_ptr);
					stacked_nodes[bktptr1]++;
					break;

				}
			}
			(edge_p2p)=&((*edge_p2p)->nxt_edge);
			edge_ptr=edge_ptr->nxt_edge;

		}


	}




	if(stacked_nodes[bktptr2]<0)
	{
		struct edge_node * edge_ptr=bktptr2->kmer_info.right;
		edge_node ** edge_p2p=&(bktptr2->kmer_info.right);
		while(edge_ptr!=NULL)
		{
			if(edge_ptr->len==edge_len-1)
			{
				kmer_t4 kmer;
				kmer=(bktptr2)->kmer_t4;

				for(int g=edge_len-1;g>=0;--g)
				{
					uint64_t b=(edge_ptr->edge>>2*g)&0x3;
					kmer.kmer[0]&=(~((uint64_t)0x3<<(2*(K_size-1-96))));
					L_shift_NB(kmer.kmer,2,4);
					kmer.kmer[3]|=b;

				}
				kmer_t4 f_kmer=kmer;
				get_rev_comp_seq_arr(f_kmer.kmer,K_size,4);
				if(uint64_t_cmp(kmer.kmer,f_kmer.kmer,4)>0)
				{
					kmer=f_kmer;
				}

				if(uint64_t_cmp(kmer.kmer,bktptr1->kmer_t4.kmer,4)==0)
				{
					(*edge_p2p)=(*edge_p2p)->nxt_edge;
					struct edge_node* f_edge_ptr=edge_ptr;
					edge_ptr=edge_ptr->nxt_edge;
					free(f_edge_ptr);

					break;

				}
			}


			(edge_p2p)=&((*edge_p2p)->nxt_edge);
			edge_ptr=edge_ptr->nxt_edge;

		}


	}


	if(stacked_nodes[bktptr2]>0)
	{
		struct edge_node * edge_ptr=bktptr2->kmer_info.left;
		edge_node ** edge_p2p=&(bktptr2->kmer_info.left);
		while(edge_ptr!=NULL)
		{
			if(edge_ptr->len==edge_len-1)
			{
				kmer_t4 kmer;
				kmer=(bktptr2)->kmer_t4;

				for(int g=0;g<edge_len;++g)
				{
					uint64_t b=(edge_ptr->edge>>2*g)&0x3;
					R_shift_NB(kmer.kmer,2,4);
					kmer.kmer[0]|=(b<<(2*(K_size-1-96)));
				}
				kmer_t4 f_kmer=kmer;
				get_rev_comp_seq_arr(f_kmer.kmer,K_size,4);
				if(uint64_t_cmp(kmer.kmer,f_kmer.kmer,4)>0)
				{
					kmer=f_kmer;
				}

				if(uint64_t_cmp(kmer.kmer,bktptr1->kmer_t4.kmer,4)==0)
				{
					(*edge_p2p)=(*edge_p2p)->nxt_edge;
					struct edge_node* f_edge_ptr=edge_ptr;
					edge_ptr=edge_ptr->nxt_edge;
					free(f_edge_ptr);

					break;

				}
			}
			(edge_p2p)=&((*edge_p2p)->nxt_edge);
			edge_ptr=edge_ptr->nxt_edge;

		}


	}

}

void BreakLinks0( map<struct bucket0* ,int > &stacked_nodes, struct bucket0* bktptr1, struct bucket0* bktptr2,int K_size,int edge_len)
{
	int Kmer_arr_sz=K_size/32+1;
	int rem1=K_size%32;
	if(rem1==0)
	{Kmer_arr_sz--;}

	if(stacked_nodes[bktptr1]>1)
	{
		struct edge_node * edge_ptr=bktptr1->kmer_info.right;
		edge_node ** edge_p2p=&(bktptr1->kmer_info.right);
		while(edge_ptr!=NULL)
		{
			if(edge_ptr->len==edge_len-1)
			{
				uint64_t kmer[100];
				memcpy(kmer,(bktptr1)->kmer_t,sizeof(uint64_t)*Kmer_arr_sz);

				for(int g=edge_len-1;g>=0;--g)
				{
					uint64_t b=(edge_ptr->edge>>2*g)&0x3;
					kmer[0]&=(~((uint64_t)0x3<<(2*((K_size-1)%32))));
					L_shift_NB(kmer,2,Kmer_arr_sz);
					kmer[Kmer_arr_sz-1]|=b;

				}
				uint64_t f_kmer[100];
				
				memcpy(f_kmer,kmer,sizeof(uint64_t)*Kmer_arr_sz);
				get_rev_comp_seq_arr(f_kmer,K_size,Kmer_arr_sz);
				if(uint64_t_cmp(kmer,f_kmer,Kmer_arr_sz)>0)
				{
					memcpy(kmer,f_kmer,sizeof(uint64_t)*Kmer_arr_sz);
				}

				if(uint64_t_cmp(kmer,bktptr2->kmer_t,Kmer_arr_sz)==0)
				{
					(*edge_p2p)=(*edge_p2p)->nxt_edge;
					struct edge_node* f_edge_ptr=edge_ptr;
					edge_ptr=edge_ptr->nxt_edge;
					free(f_edge_ptr);
					stacked_nodes[bktptr1]--;
					break;

				}
			}


			(edge_p2p)=&((*edge_p2p)->nxt_edge);
			edge_ptr=edge_ptr->nxt_edge;

		}


	}


	if(stacked_nodes[bktptr1]<-1)
	{
		struct edge_node * edge_ptr=bktptr1->kmer_info.left;
		edge_node ** edge_p2p=&(bktptr1->kmer_info.left);
		while(edge_ptr!=NULL)
		{
			if(edge_ptr->len==edge_len-1)
			{
				uint64_t kmer[100];
				
				memcpy(kmer,(bktptr1)->kmer_t,sizeof(uint64_t)*Kmer_arr_sz);


				for(int g=0;g<edge_len;++g)
				{
					uint64_t b=(edge_ptr->edge>>2*g)&0x3;
					R_shift_NB(kmer,2,Kmer_arr_sz);
					kmer[0]|=(b<<(2*((K_size-1)%32)));
				}



				uint64_t f_kmer[100];
				
				memcpy(f_kmer,kmer,sizeof(uint64_t)*Kmer_arr_sz);
				get_rev_comp_seq_arr(f_kmer,K_size,Kmer_arr_sz);
				if(uint64_t_cmp(kmer,f_kmer,Kmer_arr_sz)>0)
				{
					memcpy(kmer,f_kmer,sizeof(uint64_t)*Kmer_arr_sz);
				}

				if(uint64_t_cmp(kmer,bktptr2->kmer_t,Kmer_arr_sz)==0)
				{
					(*edge_p2p)=(*edge_p2p)->nxt_edge;
					struct edge_node* f_edge_ptr=edge_ptr;
					edge_ptr=edge_ptr->nxt_edge;
					free(f_edge_ptr);
					stacked_nodes[bktptr1]++;
					break;

				}

			}
			(edge_p2p)=&((*edge_p2p)->nxt_edge);
			edge_ptr=edge_ptr->nxt_edge;

		}


	}







	if(stacked_nodes[bktptr2]<0)
	{
		struct edge_node * edge_ptr=bktptr2->kmer_info.right;
		edge_node ** edge_p2p=&(bktptr2->kmer_info.right);
		while(edge_ptr!=NULL)
		{
			if(edge_ptr->len==edge_len-1)
			{
				uint64_t kmer[100];
				memcpy(kmer,(bktptr2)->kmer_t,sizeof(uint64_t)*Kmer_arr_sz);
				for(int g=edge_len-1;g>=0;--g)
				{
					uint64_t b=(edge_ptr->edge>>2*g)&0x3;
					kmer[0]&=(~((uint64_t)0x3<<(2*((K_size-1)%32))));
					L_shift_NB(kmer,2,Kmer_arr_sz);
					kmer[Kmer_arr_sz-1]|=b;

				}


				uint64_t f_kmer[100];
				
				memcpy(f_kmer,kmer,sizeof(uint64_t)*Kmer_arr_sz);
				get_rev_comp_seq_arr(f_kmer,K_size,Kmer_arr_sz);
				if(uint64_t_cmp(kmer,f_kmer,Kmer_arr_sz)>0)
				{
					memcpy(kmer,f_kmer,sizeof(uint64_t)*Kmer_arr_sz);
				}

				if(uint64_t_cmp(kmer,bktptr1->kmer_t,Kmer_arr_sz)==0)
				{
					(*edge_p2p)=(*edge_p2p)->nxt_edge;
					struct edge_node* f_edge_ptr=edge_ptr;
					edge_ptr=edge_ptr->nxt_edge;
					free(f_edge_ptr);

					break;

				}
			}


			(edge_p2p)=&((*edge_p2p)->nxt_edge);
			edge_ptr=edge_ptr->nxt_edge;

		}


	}


	if(stacked_nodes[bktptr2]>0)
	{
		struct edge_node * edge_ptr=bktptr2->kmer_info.left;
		edge_node ** edge_p2p=&(bktptr2->kmer_info.left);
		while(edge_ptr!=NULL)
		{
			if(edge_ptr->len==edge_len-1)
			{
				uint64_t kmer[100];
				memcpy(kmer,(bktptr2)->kmer_t,sizeof(uint64_t)*Kmer_arr_sz);



				for(int g=0;g<edge_len;++g)
				{
					uint64_t b=(edge_ptr->edge>>2*g)&0x3;
					R_shift_NB(kmer,2,Kmer_arr_sz);
					kmer[0]|=(b<<(2*((K_size-1)%32)));
				}

				uint64_t f_kmer[100];
				
				memcpy(f_kmer,kmer,sizeof(uint64_t)*Kmer_arr_sz);
				get_rev_comp_seq_arr(f_kmer,K_size,Kmer_arr_sz);
				if(uint64_t_cmp(kmer,f_kmer,Kmer_arr_sz)>0)
				{
					memcpy(kmer,f_kmer,sizeof(uint64_t)*Kmer_arr_sz);
				}



				if(uint64_t_cmp(kmer,bktptr1->kmer_t,Kmer_arr_sz)==0)
				{
					(*edge_p2p)=(*edge_p2p)->nxt_edge;
					struct edge_node* f_edge_ptr=edge_ptr;
					edge_ptr=edge_ptr->nxt_edge;
					free(f_edge_ptr);

					break;

				}
			}
			(edge_p2p)=&((*edge_p2p)->nxt_edge);
			edge_ptr=edge_ptr->nxt_edge;

		}


	}

}

void RemovingWeakNodesAndEdges(hashtable *ht,int K_size,int NodeCovTh, int EdgeCovTh,int64_t *bucket_cnt, int64_t *edge_cnt)
{
	bool ROBUST=1;
	int Removed_Nodes_cnt=0,Removed_Edges_cnt=0;
	bool DISPLAY=0;
	bucket * bktptr=NULL;

	bucket ** bktp2p=NULL;

//	struct edge_node *edge_ptr;

	if(DISPLAY)
	{
		cout<<"Removing weak nodes and edges..."<<endl;
	}
	for(size_t i=0;i<ht->ht_sz;++i)
	{
		
		size_t list_sz=0;
		bktptr=ht->store_pos[i];
		bktp2p=&(ht->store_pos[i]);
		while(bktptr!=NULL)
		{
			bucket * n_bktptr;
			struct edge_node *edge_ptr,*n_edge_ptr;
			struct edge_node **edge_p2p;
			
			if(bktptr->kmer_info.cov1<=NodeCovTh||bktptr->kmer_info.removed==1)
			{
				if(1)//!ROBUST||(bktptr->kmer_info.removed==1))
				{

					n_bktptr=bktptr->nxt_bucket;
					edge_ptr=bktptr->kmer_info.left;
					edge_p2p=&(bktptr->kmer_info.left);
					while(edge_ptr!=NULL)
					{
						n_edge_ptr=edge_ptr->nxt_edge;
						free(edge_ptr);
						edge_ptr=n_edge_ptr;
						Removed_Edges_cnt++;
					}
					*edge_p2p=NULL;
					edge_ptr=bktptr->kmer_info.right;
					edge_p2p=&(bktptr->kmer_info.right);
					while(edge_ptr!=NULL)
					{
						n_edge_ptr=edge_ptr->nxt_edge;
						free(edge_ptr);
						Removed_Edges_cnt++;
						edge_ptr=n_edge_ptr;
					}
					//*edge_p2p=NULL;
					if(1)//!bktptr->kmer_info.removed)///////////////////new....
					{
						free(bktptr);
						bktptr=n_bktptr;
						(*bktp2p)=n_bktptr;
						Removed_Nodes_cnt++;
					

					}
					else
					{
						bktp2p=&(bktptr->nxt_bucket);
						bktptr=bktptr->nxt_bucket;
				
					}
				}
				else
				{
					n_bktptr=bktptr->nxt_bucket;
					edge_ptr=bktptr->kmer_info.left;
					edge_p2p=&(bktptr->kmer_info.left);
					while(edge_ptr!=NULL)
					{
						n_edge_ptr=edge_ptr->nxt_edge;
						edge_ptr->masked=1;
						edge_ptr=n_edge_ptr;
						//Removed_Edges_cnt++;
					}
					*edge_p2p=NULL;
					edge_ptr=bktptr->kmer_info.right;
					edge_p2p=&(bktptr->kmer_info.right);
					while(edge_ptr!=NULL)
					{
						n_edge_ptr=edge_ptr->nxt_edge;
						//free(edge_ptr);
						edge_ptr->masked=1;
						//Removed_Edges_cnt++;
						edge_ptr=n_edge_ptr;
					}
					
					bktptr->kmer_info.masked=1;
					bktp2p=&(bktptr->nxt_bucket);
					bktptr=bktptr->nxt_bucket;
				
					
				}
			}
			else
			{

				edge_p2p=&(bktptr->kmer_info.left);
				edge_ptr=bktptr->kmer_info.left;
				int l_max_cov=0;
				int EdgeCovTh_t=EdgeCovTh;
				if(ROBUST)
				{
					while(edge_ptr!=NULL)
					{
						if(edge_ptr->edge_cov>l_max_cov)
						{
							l_max_cov=edge_ptr->edge_cov;
						}
						edge_ptr=edge_ptr->nxt_edge;
					}
					if(l_max_cov>20*EdgeCovTh)
					{
						EdgeCovTh_t=l_max_cov/20;
					}

				}
				
				edge_ptr=bktptr->kmer_info.left;
				while((*edge_p2p)!=NULL)
				{
					bool removed=0;
					if((*edge_p2p)->edge_cov<=EdgeCovTh_t)
					{
						removed=1;
						
					}
					if(removed==0||(ROBUST&&removed))
					{

						uint64_t t_kmer=bktptr->kmer_t.kmer,f_kmer;
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

						uint64_t hv=MurmurHash64A(&t_kmer,sizeof(t_kmer),0);

						size_t hash_idx=(size_t) (hv%ht->ht_sz);


						struct bucket ** nc_ptr;

						nc_ptr= &(ht->store_pos[hash_idx]);

						bool found=look_up_in_a_list(t_kmer,&nc_ptr);

						if(found==0||(*nc_ptr)->kmer_info.removed==1)
						{
							removed=1;
						}
						else
						{
							if(ROBUST&&removed)
							{
								//breaklinks
								if(EdgeCovTh!=EdgeCovTh_t)
								{
									map<struct bucket* ,int > stacked_nodes;
									stacked_nodes[bktptr]=-1;
									if(flip_nc)
									{
										stacked_nodes[(*nc_ptr)]=1;
									}
									else
									{
										stacked_nodes[(*nc_ptr)]=-1;
									}
									int edge_len=edge_ptr->len+1;
									BreakLinks( stacked_nodes, bktptr, (*nc_ptr),K_size,edge_len);
								}
								else
								{
									edge_ptr->masked=1;	
									removed=0;
								}

							}
						}
					}




					if(removed)
					{
						n_edge_ptr=(*edge_p2p)->nxt_edge;
						free(*edge_p2p);
						(*edge_p2p)=n_edge_ptr;
						edge_ptr=(*edge_p2p);
						Removed_Edges_cnt++;
					}
					else
					{
						edge_p2p=&((*edge_p2p)->nxt_edge);
						edge_ptr=(*edge_p2p);
					}
				}



				edge_ptr=bktptr->kmer_info.right;				
				edge_p2p=&(bktptr->kmer_info.right);

				int r_max_cov=0;
				EdgeCovTh_t=EdgeCovTh;
				if(ROBUST)
				{
					while(edge_ptr!=NULL)
					{
						if(edge_ptr->edge_cov>r_max_cov)
						{
							r_max_cov=edge_ptr->edge_cov;
						}
						edge_ptr=edge_ptr->nxt_edge;
					}
					if(r_max_cov>20*EdgeCovTh)
					{
						EdgeCovTh_t=r_max_cov/20;
					}

				}
				

				edge_ptr=bktptr->kmer_info.right;
				while((*edge_p2p)!=NULL)
				{
					bool removed=0;

					if((*edge_p2p)->edge_cov<=EdgeCovTh_t)
					{
						removed=1;
					}


					if(removed==0||(ROBUST&&removed))
					{

						uint64_t t_kmer=bktptr->kmer_t.kmer,f_kmer;
						uint64_t t=0;
						bool flip_nc=0;

						for(int j=edge_ptr->len;j>=0;--j)
						{
							int right=(int)((edge_ptr->edge)>>2*j);

							switch(right&0x3)
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

						uint64_t hv=MurmurHash64A(&t_kmer,sizeof(t_kmer),0);

						size_t hash_idx=(size_t) (hv%ht->ht_sz);


						struct bucket ** nc_ptr;

						nc_ptr= &(ht->store_pos[hash_idx]);

						bool found=look_up_in_a_list(t_kmer,&nc_ptr);


						if(found==0||(*nc_ptr)->kmer_info.removed==1)
						{
							removed=1;
						}
						if(ROBUST&&removed)
						{
							//breaklinks
							if(EdgeCovTh!=EdgeCovTh_t)
							{
								map<struct bucket* ,int > stacked_nodes;
								stacked_nodes[bktptr]=1;
								if(flip_nc)
								{
									stacked_nodes[(*nc_ptr)]=-1;
								}
								else
								{
									stacked_nodes[(*nc_ptr)]=1;
								}
								int edge_len=edge_ptr->len+1;
								BreakLinks( stacked_nodes, bktptr, (*nc_ptr),K_size,edge_len);
							}
							else
							{
								edge_ptr->masked=1;	
								removed=0;
							}

						}



					}




					if(removed)
					{
						n_edge_ptr=(*edge_p2p)->nxt_edge;
						free(*edge_p2p);
						Removed_Edges_cnt++;
						(*edge_p2p)=n_edge_ptr;
						edge_ptr=(*edge_p2p);
					}
					else
					{
						edge_p2p=&((*edge_p2p)->nxt_edge);
						edge_ptr=(*edge_p2p);
					}
				}



				if(EdgeCovTh>0&&(bktptr->kmer_info.left==NULL)&&(bktptr->kmer_info.right==NULL))
				{

					n_bktptr=bktptr->nxt_bucket;

					free(bktptr);
					bktptr=n_bktptr;
					(*bktp2p)=n_bktptr;
					Removed_Nodes_cnt++;

				}
				else
				{
					bktp2p=&(bktptr->nxt_bucket);
					bktptr=bktptr->nxt_bucket;

				}
			}

		}
	}
	if(DISPLAY)
	{
		cout<<Removed_Nodes_cnt<<" nodes removed."<<endl;
		cout<<Removed_Edges_cnt<<" edges removed."<<endl;
	}
	(*bucket_cnt)-=Removed_Nodes_cnt;
	(*edge_cnt)-=Removed_Edges_cnt;

}

void RemovingWeakNodesAndEdges2(hashtable2 *ht,int K_size,int NodeCovTh, int EdgeCovTh,int64_t *bucket_cnt, int64_t *edge_cnt)
{
	bool ROBUST=1;
	int Removed_Nodes_cnt=0,Removed_Edges_cnt=0;
	bool DISPLAY=0;
	bucket2 * bktptr=NULL;

	bucket2 ** bktp2p=NULL;

//	struct edge_node *edge_ptr;
	if(DISPLAY)
	{
	cout<<"Removing weak nodes and edges..."<<endl;
	}
	for(size_t i=0;i<ht->ht_sz;++i)
	{
		//if(i%1000000==0)
		//{cout<<i<<endl;}

		size_t list_sz=0;
		bktptr=ht->store_pos[i];
		bktp2p=&(ht->store_pos[i]);
		while(bktptr!=NULL)
		{
			bucket2 * n_bktptr;
			struct edge_node *edge_ptr,*n_edge_ptr;
			struct edge_node **edge_p2p;
			if(bktptr->kmer_info.cov1<=NodeCovTh||bktptr->kmer_info.removed==1)
			{
				if(1)//!ROBUST||(bktptr->kmer_info.removed==1))
				{

					n_bktptr=bktptr->nxt_bucket;
					edge_ptr=bktptr->kmer_info.left;
					while(edge_ptr!=NULL)
					{
						n_edge_ptr=edge_ptr->nxt_edge;
						free(edge_ptr);
						Removed_Edges_cnt++;
						edge_ptr=n_edge_ptr;
					}
					edge_ptr=bktptr->kmer_info.right;
					while(edge_ptr!=NULL)
					{
						n_edge_ptr=edge_ptr->nxt_edge;
						free(edge_ptr);
						Removed_Edges_cnt++;
						edge_ptr=n_edge_ptr;
					}
					//*edge_p2p=NULL;
					if(1)//!bktptr->kmer_info.removed)///////////////////new....
					{
						free(bktptr);
						bktptr=n_bktptr;
						(*bktp2p)=n_bktptr;
						Removed_Nodes_cnt++;
					

					}
					else
					{
						bktp2p=&(bktptr->nxt_bucket);
						bktptr=bktptr->nxt_bucket;
				
					}
				}
				else
				{
					n_bktptr=bktptr->nxt_bucket;
					edge_ptr=bktptr->kmer_info.left;
					edge_p2p=&(bktptr->kmer_info.left);
					while(edge_ptr!=NULL)
					{
						n_edge_ptr=edge_ptr->nxt_edge;
						edge_ptr->masked=1;
						edge_ptr=n_edge_ptr;
						//Removed_Edges_cnt++;
					}
					*edge_p2p=NULL;
					edge_ptr=bktptr->kmer_info.right;
					edge_p2p=&(bktptr->kmer_info.right);
					while(edge_ptr!=NULL)
					{
						n_edge_ptr=edge_ptr->nxt_edge;
						//free(edge_ptr);
						edge_ptr->masked=1;
						//Removed_Edges_cnt++;
						edge_ptr=n_edge_ptr;
					}
					
					bktptr->kmer_info.masked=1;
					bktp2p=&(bktptr->nxt_bucket);
					bktptr=bktptr->nxt_bucket;
				
					
				}

			}
			else
			{

				edge_p2p=&(bktptr->kmer_info.left);
				edge_node *edge_ptr=bktptr->kmer_info.left;

				int l_max_cov=0;
				int EdgeCovTh_t=EdgeCovTh;
				if(ROBUST)
				{
					while(edge_ptr!=NULL)
					{
						if(edge_ptr->edge_cov>l_max_cov)
						{
							l_max_cov=edge_ptr->edge_cov;
						}
						edge_ptr=edge_ptr->nxt_edge;
					}
					if(l_max_cov>20*EdgeCovTh)
					{
						EdgeCovTh_t=l_max_cov/20;
					}

				}
				edge_ptr=bktptr->kmer_info.left;
				while((*edge_p2p)!=NULL)
				{
					bool removed=0;

					if((*edge_p2p)->edge_cov<=EdgeCovTh_t)
					{
						removed=1;
					}
					if(removed==0||(ROBUST&&removed))
					{
						kmer_t2 t_kmer=bktptr->kmer_t2,f_kmer;
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

						uint64_t hv=MurmurHash64A(&t_kmer,sizeof(t_kmer),0);

						size_t hash_idx=(size_t) (hv%ht->ht_sz);


						struct bucket2 ** nc_ptr;

						nc_ptr= &(ht->store_pos[hash_idx]);

						bool found=look_up_in_a_list2(&t_kmer,&nc_ptr);


						if(found==0||(*nc_ptr)->kmer_info.removed==1)
						{
							removed=1;
						}
						else
						{
							if(ROBUST&&removed)
							{
								//breaklinks
								if(EdgeCovTh!=EdgeCovTh_t)
								{
									map<struct bucket2* ,int > stacked_nodes;
									stacked_nodes[bktptr]=-1;
									if(flip_nc)
									{
										stacked_nodes[(*nc_ptr)]=1;
									}
									else
									{
										stacked_nodes[(*nc_ptr)]=-1;
									}
									int edge_len=edge_ptr->len+1;
									BreakLinks2( stacked_nodes, bktptr, (*nc_ptr),K_size,edge_len);
								}
								else
								{
									edge_ptr->masked=1;	
									removed=0;
								}

							}
						}


					}




					if(removed)
					{

						n_edge_ptr=(*edge_p2p)->nxt_edge;
						free(*edge_p2p);
						Removed_Edges_cnt++;
						(*edge_p2p)=n_edge_ptr;
						edge_ptr=(*edge_p2p);
					}
					else
					{
						edge_p2p=&((*edge_p2p)->nxt_edge);
						edge_ptr=(*edge_p2p);
					}
				}



				edge_p2p=&(bktptr->kmer_info.right);
				edge_ptr=bktptr->kmer_info.right;

				int r_max_cov=0;
				EdgeCovTh_t=EdgeCovTh;
				if(ROBUST)
				{
					while(edge_ptr!=NULL)
					{
						if(edge_ptr->edge_cov>r_max_cov)
						{
							r_max_cov=edge_ptr->edge_cov;
						}
						edge_ptr=edge_ptr->nxt_edge;
					}
					if(r_max_cov>20*EdgeCovTh)
					{
						EdgeCovTh_t=r_max_cov/20;
					}

				}
				edge_ptr=bktptr->kmer_info.right;
				while((*edge_p2p)!=NULL)
				{
					bool removed=0;
					if((*edge_p2p)->edge_cov<=EdgeCovTh_t)
					{
						removed=1;
					}

					if(removed==0||(ROBUST&&removed))
					{
						kmer_t2 t_kmer=bktptr->kmer_t2,f_kmer;
						uint64_t t=0;
						bool flip_nc=0;

						for(int j=edge_ptr->len;j>=0;--j)
						{
							int right=(int)((edge_ptr->edge)>>2*j);

							//right=(edge_ptr->edge)>>2*j;

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

						uint64_t hv=MurmurHash64A(&t_kmer,sizeof(t_kmer),0);

						size_t hash_idx=(size_t) (hv%ht->ht_sz);


						struct bucket2 ** nc_ptr;

						nc_ptr= &(ht->store_pos[hash_idx]);

						bool found=look_up_in_a_list2(&t_kmer,&nc_ptr);


						if(found==0||(*nc_ptr)->kmer_info.removed==1)
						{

							removed=1;
						}
						if(ROBUST&&removed)
						{
							//breaklinks
							if(EdgeCovTh!=EdgeCovTh_t)
							{
								map<struct bucket2* ,int > stacked_nodes;
								stacked_nodes[bktptr]=1;
								if(flip_nc)
								{
									stacked_nodes[(*nc_ptr)]=-1;
								}
								else
								{
									stacked_nodes[(*nc_ptr)]=1;
								}
								int edge_len=edge_ptr->len+1;
								BreakLinks2( stacked_nodes, bktptr, (*nc_ptr),K_size,edge_len);
							}
							else
							{
								edge_ptr->masked=1;	
								removed=0;
							}

						}



					}



					if(removed)
					{
						n_edge_ptr=(*edge_p2p)->nxt_edge;
						free(*edge_p2p);
						Removed_Edges_cnt++;
						(*edge_p2p)=n_edge_ptr;
						edge_ptr=(*edge_p2p);
					}
					else
					{
						edge_p2p=&((*edge_p2p)->nxt_edge);
						edge_ptr=(*edge_p2p);
					}
				}

				if(EdgeCovTh>0&&(bktptr->kmer_info.left==NULL)&&(bktptr->kmer_info.right==NULL))
				{
					n_bktptr=bktptr->nxt_bucket;



					free(bktptr);

					Removed_Nodes_cnt++;
					bktptr=n_bktptr;
					(*bktp2p)=n_bktptr;

				}
				else
				{
					bktp2p=&(bktptr->nxt_bucket);
					bktptr=bktptr->nxt_bucket;

				}
			}

		}
	}
	if(DISPLAY)
	{
	cout<<Removed_Nodes_cnt<<" nodes removed."<<endl;
	cout<<Removed_Edges_cnt<<" edges removed."<<endl;
	}
	(*bucket_cnt)-=Removed_Nodes_cnt;
	(*edge_cnt)-=Removed_Edges_cnt;
}


void RemovingWeakNodesAndEdges3(hashtable3 *ht,int K_size,int NodeCovTh, int EdgeCovTh,int64_t *bucket_cnt, int64_t *edge_cnt)
{
	bool ROBUST=1;
	int Removed_Nodes_cnt=0,Removed_Edges_cnt=0;
	bool DISPLAY=0;
	bucket3 * bktptr=NULL;

	bucket3 ** bktp2p=NULL;

//	struct edge_node *edge_ptr;
	if(DISPLAY)
	{
	cout<<"Removing weak nodes and edges..."<<endl;
	}
	for(size_t i=0;i<ht->ht_sz;++i)
	{
		//if(i%1000000==0)
		//{cout<<i<<endl;}

		size_t list_sz=0;
		bktptr=ht->store_pos[i];
		bktp2p=&(ht->store_pos[i]);
		while(bktptr!=NULL)
		{
			bucket3 * n_bktptr;
			struct edge_node *edge_ptr,*n_edge_ptr;
			struct edge_node **edge_p2p;
			if(bktptr->kmer_info.cov1<=NodeCovTh||bktptr->kmer_info.removed==1)
			{
				if(1)//!ROBUST||(bktptr->kmer_info.removed==1))
				{
					n_bktptr=bktptr->nxt_bucket;
					edge_ptr=bktptr->kmer_info.left;
					while(edge_ptr!=NULL)
					{
						n_edge_ptr=edge_ptr->nxt_edge;
						free(edge_ptr);
						Removed_Edges_cnt++;
						edge_ptr=n_edge_ptr;
					}
					edge_ptr=bktptr->kmer_info.right;
					while(edge_ptr!=NULL)
					{
						n_edge_ptr=edge_ptr->nxt_edge;
						free(edge_ptr);
						Removed_Edges_cnt++;
						edge_ptr=n_edge_ptr;
					}	

					//*edge_p2p=NULL;
					if(1)//!bktptr->kmer_info.removed)///////////////////new....
					{
						free(bktptr);
						bktptr=n_bktptr;
						(*bktp2p)=n_bktptr;
						Removed_Nodes_cnt++;
					

					}
					else
					{
						bktp2p=&(bktptr->nxt_bucket);
						bktptr=bktptr->nxt_bucket;
				
					}
				}
				else
				{
					n_bktptr=bktptr->nxt_bucket;
					edge_ptr=bktptr->kmer_info.left;
					edge_p2p=&(bktptr->kmer_info.left);
					while(edge_ptr!=NULL)
					{
						n_edge_ptr=edge_ptr->nxt_edge;
						edge_ptr->masked=1;
						edge_ptr=n_edge_ptr;
						//Removed_Edges_cnt++;
					}
					*edge_p2p=NULL;
					edge_ptr=bktptr->kmer_info.right;
					edge_p2p=&(bktptr->kmer_info.right);
					while(edge_ptr!=NULL)
					{
						n_edge_ptr=edge_ptr->nxt_edge;
						//free(edge_ptr);
						edge_ptr->masked=1;
						//Removed_Edges_cnt++;
						edge_ptr=n_edge_ptr;
					}
					
					bktptr->kmer_info.masked=1;
					bktp2p=&(bktptr->nxt_bucket);
					bktptr=bktptr->nxt_bucket;
				
					
				}
			}
			else
			{

				edge_p2p=&(bktptr->kmer_info.left);
				edge_node *edge_ptr=bktptr->kmer_info.left;
				int l_max_cov=0;
				int EdgeCovTh_t=EdgeCovTh;
				if(ROBUST)
				{
					while(edge_ptr!=NULL)
					{
						if(edge_ptr->edge_cov>l_max_cov)
						{
							l_max_cov=edge_ptr->edge_cov;
						}
						edge_ptr=edge_ptr->nxt_edge;
					}
					if(l_max_cov>20*EdgeCovTh)
					{
						EdgeCovTh_t=l_max_cov/20;
					}

				}
				edge_ptr=bktptr->kmer_info.left;
				while((*edge_p2p)!=NULL)
				{
					bool removed=0;

					if((*edge_p2p)->edge_cov<=EdgeCovTh_t)
					{
						removed=1;
					}
					if(removed==0||(ROBUST&&removed))
					{
						kmer_t3 t_kmer=bktptr->kmer_t3,f_kmer;
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

						uint64_t hv=MurmurHash64A(&t_kmer,sizeof(t_kmer),0);

						size_t hash_idx=(size_t) (hv%ht->ht_sz);


						struct bucket3 ** nc_ptr;

						nc_ptr= &(ht->store_pos[hash_idx]);

						bool found=look_up_in_a_list3(&t_kmer,&nc_ptr);


						if(found==0||(*nc_ptr)->kmer_info.removed==1)
						{
							removed=1;
						}
						else
						{
							if(ROBUST&&removed)
							{
								//breaklinks
								if(EdgeCovTh!=EdgeCovTh_t)
								{
									map<struct bucket3* ,int > stacked_nodes;
									stacked_nodes[bktptr]=-1;
									if(flip_nc)
									{
										stacked_nodes[(*nc_ptr)]=1;
									}
									else
									{
										stacked_nodes[(*nc_ptr)]=-1;
									}
									int edge_len=edge_ptr->len+1;
									BreakLinks3( stacked_nodes, bktptr, (*nc_ptr),K_size,edge_len);
								}
								else
								{
									edge_ptr->masked=1;	
									removed=0;
								}

							}
						}


					}




					if(removed)
					{

						n_edge_ptr=(*edge_p2p)->nxt_edge;
						free(*edge_p2p);
						Removed_Edges_cnt++;
						(*edge_p2p)=n_edge_ptr;
						edge_ptr=(*edge_p2p);
					}
					else
					{
						edge_p2p=&((*edge_p2p)->nxt_edge);
						edge_ptr=(*edge_p2p);
					}
				}



				edge_p2p=&(bktptr->kmer_info.right);
				edge_ptr=bktptr->kmer_info.right;

				int r_max_cov=0;
				EdgeCovTh_t=EdgeCovTh;
				if(ROBUST)
				{
					while(edge_ptr!=NULL)
					{
						if(edge_ptr->edge_cov>r_max_cov)
						{
							r_max_cov=edge_ptr->edge_cov;
						}
						edge_ptr=edge_ptr->nxt_edge;
					}
					if(r_max_cov>20*EdgeCovTh)
					{
						EdgeCovTh_t=r_max_cov/20;
					}

				}

				edge_ptr=bktptr->kmer_info.right;
				while((*edge_p2p)!=NULL)
				{
					bool removed=0;
					if((*edge_p2p)->edge_cov<=EdgeCovTh_t)
					{
						removed=1;
					}

					if(removed==0||(ROBUST&&removed))
					{
						kmer_t3 t_kmer=bktptr->kmer_t3,f_kmer;
						uint64_t t=0;
						bool flip_nc=0;

						for(int j=edge_ptr->len;j>=0;--j)
						{
							int right=(int)((edge_ptr->edge)>>2*j);

							//right=(edge_ptr->edge)>>2*j;

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

						uint64_t hv=MurmurHash64A(&t_kmer,sizeof(t_kmer),0);

						size_t hash_idx=(size_t) (hv%ht->ht_sz);


						struct bucket3 ** nc_ptr;

						nc_ptr= &(ht->store_pos[hash_idx]);

						bool found=look_up_in_a_list3(&t_kmer,&nc_ptr);


						if(found==0||(*nc_ptr)->kmer_info.removed==1)
						{

							removed=1;
						}
						if(ROBUST&&removed)
						{
							//breaklinks
							if(EdgeCovTh!=EdgeCovTh_t)
							{
								map<struct bucket3* ,int > stacked_nodes;
								stacked_nodes[bktptr]=1;
								if(flip_nc)
								{
									stacked_nodes[(*nc_ptr)]=-1;
								}
								else
								{
									stacked_nodes[(*nc_ptr)]=1;
								}
								int edge_len=edge_ptr->len+1;
								BreakLinks3( stacked_nodes, bktptr, (*nc_ptr),K_size,edge_len);
							}
							else
							{
								edge_ptr->masked=1;	
								removed=0;
							}

						}



					}



					if(removed)
					{
						n_edge_ptr=(*edge_p2p)->nxt_edge;
						free(*edge_p2p);
						Removed_Edges_cnt++;
						(*edge_p2p)=n_edge_ptr;
						edge_ptr=(*edge_p2p);
					}
					else
					{
						edge_p2p=&((*edge_p2p)->nxt_edge);
						edge_ptr=(*edge_p2p);
					}
				}

				if(EdgeCovTh>0&&(bktptr->kmer_info.left==NULL)&&(bktptr->kmer_info.right==NULL))
				{
					n_bktptr=bktptr->nxt_bucket;



					free(bktptr);

					Removed_Nodes_cnt++;
					bktptr=n_bktptr;
					(*bktp2p)=n_bktptr;

				}
				else
				{
					bktp2p=&(bktptr->nxt_bucket);
					bktptr=bktptr->nxt_bucket;

				}
			}

		}
	}
	if(DISPLAY)
	{
	cout<<Removed_Nodes_cnt<<" nodes removed."<<endl;
	cout<<Removed_Edges_cnt<<" edges removed."<<endl;
	}
	(*bucket_cnt)-=Removed_Nodes_cnt;
	(*edge_cnt)-=Removed_Edges_cnt;
}


void RemovingWeakNodesAndEdges4(hashtable4 *ht,int K_size,int NodeCovTh, int EdgeCovTh,int64_t *bucket_cnt, int64_t *edge_cnt)
{
	bool ROBUST=1;
	int Removed_Nodes_cnt=0,Removed_Edges_cnt=0;
	bool DISPLAY=0;
	bucket4 * bktptr=NULL;

	bucket4 ** bktp2p=NULL;

//	struct edge_node *edge_ptr;
	if(DISPLAY)
	{
	cout<<"Removing weak nodes and edges..."<<endl;
	}
	for(size_t i=0;i<ht->ht_sz;++i)
	{
		//if(i%1000000==0)
		//{cout<<i<<endl;}

		size_t list_sz=0;
		bktptr=ht->store_pos[i];
		bktp2p=&(ht->store_pos[i]);
		while(bktptr!=NULL)
		{
			bucket4 * n_bktptr;
			struct edge_node *edge_ptr,*n_edge_ptr;
			struct edge_node **edge_p2p;
			if(bktptr->kmer_info.cov1<=NodeCovTh||bktptr->kmer_info.removed==1)
			{
				if(1)//!ROBUST||(bktptr->kmer_info.removed==1))
				{
					n_bktptr=bktptr->nxt_bucket;
					edge_ptr=bktptr->kmer_info.left;
					while(edge_ptr!=NULL)
					{
						n_edge_ptr=edge_ptr->nxt_edge;
						free(edge_ptr);
						Removed_Edges_cnt++;
						edge_ptr=n_edge_ptr;
					}
					edge_ptr=bktptr->kmer_info.right;
					while(edge_ptr!=NULL)
					{
						n_edge_ptr=edge_ptr->nxt_edge;
						free(edge_ptr);
						Removed_Edges_cnt++;
						edge_ptr=n_edge_ptr;
					}
					//*edge_p2p=NULL;
					if(1)//!bktptr->kmer_info.removed)///////////////////new....
					{
						free(bktptr);
						bktptr=n_bktptr;
						(*bktp2p)=n_bktptr;
						Removed_Nodes_cnt++;
					

					}
					else
					{
						bktp2p=&(bktptr->nxt_bucket);
						bktptr=bktptr->nxt_bucket;
				
					}
				}
				else
				{
					n_bktptr=bktptr->nxt_bucket;
					edge_ptr=bktptr->kmer_info.left;
					edge_p2p=&(bktptr->kmer_info.left);
					while(edge_ptr!=NULL)
					{
						n_edge_ptr=edge_ptr->nxt_edge;
						edge_ptr->masked=1;
						edge_ptr=n_edge_ptr;
						//Removed_Edges_cnt++;
					}
					*edge_p2p=NULL;
					edge_ptr=bktptr->kmer_info.right;
					edge_p2p=&(bktptr->kmer_info.right);
					while(edge_ptr!=NULL)
					{
						n_edge_ptr=edge_ptr->nxt_edge;
						//free(edge_ptr);
						edge_ptr->masked=1;
						//Removed_Edges_cnt++;
						edge_ptr=n_edge_ptr;
					}
					
					bktptr->kmer_info.masked=1;
					bktp2p=&(bktptr->nxt_bucket);
					bktptr=bktptr->nxt_bucket;
				
					
				}

			}
			else
			{

				edge_p2p=&(bktptr->kmer_info.left);
				edge_node *edge_ptr=bktptr->kmer_info.left;

				int l_max_cov=0;
				int EdgeCovTh_t=EdgeCovTh;
				if(ROBUST)
				{
					while(edge_ptr!=NULL)
					{
						if(edge_ptr->edge_cov>l_max_cov)
						{
							l_max_cov=edge_ptr->edge_cov;
						}
						edge_ptr=edge_ptr->nxt_edge;
					}
					if(l_max_cov>20*EdgeCovTh)
					{
						EdgeCovTh_t=l_max_cov/20;
					}

				}
				edge_ptr=bktptr->kmer_info.left;
				while((*edge_p2p)!=NULL)
				{
					bool removed=0;

					if((*edge_p2p)->edge_cov<=EdgeCovTh_t)
					{
						removed=1;
					}
					if(removed==0||(ROBUST&&removed))
					{
						kmer_t4 t_kmer=bktptr->kmer_t4,f_kmer;
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

						uint64_t hv=MurmurHash64A(&t_kmer,sizeof(t_kmer),0);

						size_t hash_idx=(size_t) (hv%ht->ht_sz);


						struct bucket4 ** nc_ptr;

						nc_ptr= &(ht->store_pos[hash_idx]);

						bool found=look_up_in_a_list4(&t_kmer,&nc_ptr);


						if(found==0||(*nc_ptr)->kmer_info.removed==1)
						{
							removed=1;
						}
						else
						{
							if(ROBUST&&removed)
							{
								//breaklinks
								if(EdgeCovTh!=EdgeCovTh_t)
								{
									map<struct bucket4* ,int > stacked_nodes;
									stacked_nodes[bktptr]=-1;
									if(flip_nc)
									{
										stacked_nodes[(*nc_ptr)]=1;
									}
									else
									{
										stacked_nodes[(*nc_ptr)]=-1;
									}
									int edge_len=edge_ptr->len+1;
									BreakLinks4( stacked_nodes, bktptr, (*nc_ptr),K_size,edge_len);
								}
								else
								{
									edge_ptr->masked=1;	
									removed=0;
								}

							}
						}


					}




					if(removed)
					{

						n_edge_ptr=(*edge_p2p)->nxt_edge;
						free(*edge_p2p);
						Removed_Edges_cnt++;
						(*edge_p2p)=n_edge_ptr;
						edge_ptr=(*edge_p2p);
					}
					else
					{
						edge_p2p=&((*edge_p2p)->nxt_edge);
						edge_ptr=(*edge_p2p);
					}
				}



				edge_p2p=&(bktptr->kmer_info.right);
				edge_ptr=bktptr->kmer_info.right;
				int r_max_cov=0;
				EdgeCovTh_t=EdgeCovTh;
				if(ROBUST)
				{
					while(edge_ptr!=NULL)
					{
						if(edge_ptr->edge_cov>r_max_cov)
						{
							r_max_cov=edge_ptr->edge_cov;
						}
						edge_ptr=edge_ptr->nxt_edge;
					}
					if(r_max_cov>20*EdgeCovTh)
					{
						EdgeCovTh_t=r_max_cov/20;
					}

				}
				edge_ptr=bktptr->kmer_info.right;
				while((*edge_p2p)!=NULL)
				{
					bool removed=0;
					if((*edge_p2p)->edge_cov<=EdgeCovTh_t)
					{
						removed=1;
					}

					if(removed==0||(ROBUST&&removed))
					{
						kmer_t4 t_kmer=bktptr->kmer_t4,f_kmer;
						uint64_t t=0;
						bool flip_nc=0;

						for(int j=edge_ptr->len;j>=0;--j)
						{
							int right=(int)((edge_ptr->edge)>>2*j);

							//right=(edge_ptr->edge)>>2*j;

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

						uint64_t hv=MurmurHash64A(&t_kmer,sizeof(t_kmer),0);

						size_t hash_idx=(size_t) (hv%ht->ht_sz);


						struct bucket4 ** nc_ptr;

						nc_ptr= &(ht->store_pos[hash_idx]);

						bool found=look_up_in_a_list4(&t_kmer,&nc_ptr);


						if(found==0||(*nc_ptr)->kmer_info.removed==1)
						{

							removed=1;
						}
						if(ROBUST&&removed)
						{
							//breaklinks
							if(EdgeCovTh!=EdgeCovTh_t)
							{
								map<struct bucket4* ,int > stacked_nodes;
								stacked_nodes[bktptr]=1;
								if(flip_nc)
								{
									stacked_nodes[(*nc_ptr)]=-1;
								}
								else
								{
									stacked_nodes[(*nc_ptr)]=1;
								}
								int edge_len=edge_ptr->len+1;
								BreakLinks4( stacked_nodes, bktptr, (*nc_ptr),K_size,edge_len);
							}
							else
							{
								edge_ptr->masked=1;	
								removed=0;
							}

						}



					}



					if(removed)
					{
						n_edge_ptr=(*edge_p2p)->nxt_edge;
						free(*edge_p2p);
						Removed_Edges_cnt++;
						(*edge_p2p)=n_edge_ptr;
						edge_ptr=(*edge_p2p);
					}
					else
					{
						edge_p2p=&((*edge_p2p)->nxt_edge);
						edge_ptr=(*edge_p2p);
					}
				}

				if(EdgeCovTh>0&&(bktptr->kmer_info.left==NULL)&&(bktptr->kmer_info.right==NULL))
				{
					n_bktptr=bktptr->nxt_bucket;



					free(bktptr);

					Removed_Nodes_cnt++;
					bktptr=n_bktptr;
					(*bktp2p)=n_bktptr;

				}
				else
				{
					bktp2p=&(bktptr->nxt_bucket);
					bktptr=bktptr->nxt_bucket;

				}
			}

		}
	}
	if(DISPLAY)
	{
	cout<<Removed_Nodes_cnt<<" nodes removed."<<endl;
	cout<<Removed_Edges_cnt<<" edges removed."<<endl;
	}
	(*bucket_cnt)-=Removed_Nodes_cnt;
	(*edge_cnt)-=Removed_Edges_cnt;
}



void RemovingWeakNodesAndEdges0(hashtable0 *ht,int K_size,int NodeCovTh, int EdgeCovTh,int64_t *bucket_cnt, int64_t *edge_cnt)
{
	bool ROBUST=1;
	int Removed_Nodes_cnt=0,Removed_Edges_cnt=0;
	bool DISPLAY=0;
	bucket0 * bktptr=NULL;

	bucket0 ** bktp2p=NULL;

	int Kmer_arr_sz=K_size/32+1;
	int rem1=K_size%32;
	if(rem1==0)
	{Kmer_arr_sz--;}

//	struct edge_node *edge_ptr;
	if(DISPLAY)
	{
	cout<<"Removing weak nodes and edges..."<<endl;
	}
	for(size_t i=0;i<ht->ht_sz;++i)
	{
		//if(i%1000000==0)
		//{cout<<i<<endl;}

		size_t list_sz=0;
		bktptr=ht->store_pos[i];
		bktp2p=&(ht->store_pos[i]);
		while(bktptr!=NULL)
		{
			bucket0 * n_bktptr;
			struct edge_node *edge_ptr,*n_edge_ptr;
			struct edge_node **edge_p2p;
			if(bktptr->kmer_info.cov1<=NodeCovTh||bktptr->kmer_info.removed==1)
			{
				if(1)//!ROBUST||(bktptr->kmer_info.removed==1))
				{
					n_bktptr=bktptr->nxt_bucket;
					edge_ptr=bktptr->kmer_info.left;
					while(edge_ptr!=NULL)
					{
						n_edge_ptr=edge_ptr->nxt_edge;
						free(edge_ptr);
						Removed_Edges_cnt++;
						edge_ptr=n_edge_ptr;
					}
					edge_ptr=bktptr->kmer_info.right;
					while(edge_ptr!=NULL)
					{
						n_edge_ptr=edge_ptr->nxt_edge;
						free(edge_ptr);
						Removed_Edges_cnt++;
						edge_ptr=n_edge_ptr;
					}
					//*edge_p2p=NULL;
					if(1)//!bktptr->kmer_info.removed)///////////////////new....
					{
						free(bktptr);
						bktptr=n_bktptr;
						(*bktp2p)=n_bktptr;
						Removed_Nodes_cnt++;
					

					}
					else
					{
						bktp2p=&(bktptr->nxt_bucket);
						bktptr=bktptr->nxt_bucket;
				
					}
				}
				else
				{
					n_bktptr=bktptr->nxt_bucket;
					edge_ptr=bktptr->kmer_info.left;
					edge_p2p=&(bktptr->kmer_info.left);
					while(edge_ptr!=NULL)
					{
						n_edge_ptr=edge_ptr->nxt_edge;
						edge_ptr->masked=1;
						edge_ptr=n_edge_ptr;
						//Removed_Edges_cnt++;
					}
					*edge_p2p=NULL;
					edge_ptr=bktptr->kmer_info.right;
					edge_p2p=&(bktptr->kmer_info.right);
					while(edge_ptr!=NULL)
					{
						n_edge_ptr=edge_ptr->nxt_edge;
						//free(edge_ptr);
						edge_ptr->masked=1;
						//Removed_Edges_cnt++;
						edge_ptr=n_edge_ptr;
					}
					
					bktptr->kmer_info.masked=1;
					bktp2p=&(bktptr->nxt_bucket);
					bktptr=bktptr->nxt_bucket;
				
					
				}

			}
			else
			{

				edge_p2p=&(bktptr->kmer_info.left);
				edge_node *edge_ptr=bktptr->kmer_info.left;

				int l_max_cov=0;
				int EdgeCovTh_t=EdgeCovTh;
				if(ROBUST)
				{
					while(edge_ptr!=NULL)
					{
						if(edge_ptr->edge_cov>l_max_cov)
						{
							l_max_cov=edge_ptr->edge_cov;
						}
						edge_ptr=edge_ptr->nxt_edge;
					}
					if(l_max_cov>20*EdgeCovTh)
					{
						EdgeCovTh_t=l_max_cov/20;
					}

				}
				edge_ptr=bktptr->kmer_info.left;
				while((*edge_p2p)!=NULL)
				{
					bool removed=0;

					if((*edge_p2p)->edge_cov<=EdgeCovTh_t)
					{
						removed=1;
					}
					if(removed==0||(ROBUST&&removed))
					{
						
							
						uint64_t t_kmer[100],f_kmer[100];
						memcpy(t_kmer,bktptr->kmer_t,sizeof(uint64_t)*Kmer_arr_sz);
						uint64_t t=0;
						bool flip_nc=0;

						for(int j=0;j<=edge_ptr->len;++j)
						{
							int left=(int)((edge_ptr->edge)>>2*j);
							uint64_t b=left&0x3;
							R_shift_NB(t_kmer,2,Kmer_arr_sz);
							t=b<<(((K_size-1)%32)*2);
							t_kmer[0]|=t;

						}

						flip_nc=0;
						memcpy(f_kmer,t_kmer,sizeof(uint64_t)*Kmer_arr_sz);

						get_rev_comp_seq_arr(f_kmer,K_size,Kmer_arr_sz);
						if(uint64_t_cmp(t_kmer,f_kmer,Kmer_arr_sz)>0)
						{
							memcpy(t_kmer,f_kmer,sizeof(uint64_t)*Kmer_arr_sz);
							
							flip_nc=!flip_nc;
						}

						uint64_t hv=MurmurHash64A(t_kmer,sizeof(uint64_t)*Kmer_arr_sz,0);

						size_t hash_idx=(size_t) (hv%ht->ht_sz);

						struct bucket0 ** nc_ptr;

						nc_ptr= &(ht->store_pos[hash_idx]);

						bool found=look_up_in_a_list0(t_kmer,&nc_ptr,Kmer_arr_sz);


						if(found==0||(*nc_ptr)->kmer_info.removed==1)
						{
							removed=1;
						}
						else
						{
							if(ROBUST&&removed)
							{
								//breaklinks
								if(EdgeCovTh!=EdgeCovTh_t)
								{
									map<struct bucket0* ,int > stacked_nodes;
									stacked_nodes[bktptr]=-1;
									if(flip_nc)
									{
										stacked_nodes[(*nc_ptr)]=1;
									}
									else
									{
										stacked_nodes[(*nc_ptr)]=-1;
									}
									int edge_len=edge_ptr->len+1;
									BreakLinks0( stacked_nodes, bktptr, (*nc_ptr),K_size,edge_len);
								}
								else
								{
									edge_ptr->masked=1;	
									removed=0;
								}

							}
						}


					}




					if(removed)
					{

						n_edge_ptr=(*edge_p2p)->nxt_edge;
						free(*edge_p2p);
						Removed_Edges_cnt++;
						(*edge_p2p)=n_edge_ptr;
						edge_ptr=(*edge_p2p);
					}
					else
					{
						edge_p2p=&((*edge_p2p)->nxt_edge);
						edge_ptr=(*edge_p2p);
					}
				}



				edge_p2p=&(bktptr->kmer_info.right);
				edge_ptr=bktptr->kmer_info.right;
				int r_max_cov=0;
				EdgeCovTh_t=EdgeCovTh;
				if(ROBUST)
				{
					while(edge_ptr!=NULL)
					{
						if(edge_ptr->edge_cov>r_max_cov)
						{
							r_max_cov=edge_ptr->edge_cov;
						}
						edge_ptr=edge_ptr->nxt_edge;
					}
					if(r_max_cov>20*EdgeCovTh)
					{
						EdgeCovTh_t=r_max_cov/20;
					}

				}
				edge_ptr=bktptr->kmer_info.right;
				while((*edge_p2p)!=NULL)
				{
					bool removed=0;
					if((*edge_p2p)->edge_cov<=EdgeCovTh_t)
					{
						removed=1;
					}

					if(removed==0||(ROBUST&&removed))
					{



						uint64_t t_kmer[100],f_kmer[100];
						memcpy(t_kmer,bktptr->kmer_t,sizeof(uint64_t)*Kmer_arr_sz);
						uint64_t t=0;
						bool flip_nc=0;

				

						

						for(int j=edge_ptr->len;j>=0;--j)
						{
							int right=(int)((edge_ptr->edge)>>2*j);

							//right=(edge_ptr->edge)>>2*j;
							uint64_t b=right&0x3;
							t=3;
							t<<=((K_size-1)%32)*2;
							t_kmer[0]&=(~t);
							L_shift_NB(t_kmer,2,Kmer_arr_sz);//
							t_kmer[Kmer_arr_sz-1]|=b;



						}

						flip_nc=0;

						memcpy(f_kmer,t_kmer,sizeof(uint64_t)*Kmer_arr_sz);

						get_rev_comp_seq_arr(f_kmer,K_size,Kmer_arr_sz);
						if(uint64_t_cmp(t_kmer,f_kmer,Kmer_arr_sz)>0)
						{
							memcpy(t_kmer,f_kmer,sizeof(uint64_t)*Kmer_arr_sz);
							
							flip_nc=!flip_nc;
						}

						uint64_t hv=MurmurHash64A(t_kmer,sizeof(uint64_t)*Kmer_arr_sz,0);


						size_t hash_idx=(size_t) (hv%ht->ht_sz);

						struct bucket0 ** nc_ptr;

						nc_ptr= &(ht->store_pos[hash_idx]);

						bool found=look_up_in_a_list0(t_kmer,&nc_ptr,Kmer_arr_sz);


						if(found==0||(*nc_ptr)->kmer_info.removed==1)
						{

							removed=1;
						}
						if(ROBUST&&removed)
						{
							//breaklinks
							if(EdgeCovTh!=EdgeCovTh_t)
							{
								map<struct bucket0* ,int > stacked_nodes;
								stacked_nodes[bktptr]=1;
								if(flip_nc)
								{
									stacked_nodes[(*nc_ptr)]=-1;
								}
								else
								{
									stacked_nodes[(*nc_ptr)]=1;
								}
								int edge_len=edge_ptr->len+1;
								BreakLinks0( stacked_nodes, bktptr, (*nc_ptr),K_size,edge_len);
							}
							else
							{
								edge_ptr->masked=1;	
								removed=0;
							}

						}



					}



					if(removed)
					{
						n_edge_ptr=(*edge_p2p)->nxt_edge;
						free(*edge_p2p);
						Removed_Edges_cnt++;
						(*edge_p2p)=n_edge_ptr;
						edge_ptr=(*edge_p2p);
					}
					else
					{
						edge_p2p=&((*edge_p2p)->nxt_edge);
						edge_ptr=(*edge_p2p);
					}
				}

				if(EdgeCovTh>0&&(bktptr->kmer_info.left==NULL)&&(bktptr->kmer_info.right==NULL))
				{
					n_bktptr=bktptr->nxt_bucket;



					free(bktptr);

					Removed_Nodes_cnt++;
					bktptr=n_bktptr;
					(*bktp2p)=n_bktptr;

				}
				else
				{
					bktp2p=&(bktptr->nxt_bucket);
					bktptr=bktptr->nxt_bucket;

				}
			}

		}
	}
	if(DISPLAY)
	{
	cout<<Removed_Nodes_cnt<<" nodes removed."<<endl;
	cout<<Removed_Edges_cnt<<" edges removed."<<endl;
	}
	(*bucket_cnt)-=Removed_Nodes_cnt;
	(*edge_cnt)-=Removed_Edges_cnt;
}



void RemovingWeakNodes_r1(hashtable *ht,hashtable2 *ht2,int K_size,int NodeCovTh, int64_t *bucket_cnt)
{
	int Removed_Nodes_cnt=0;
	if(K_size<=32)
	{

		for(size_t i=0;i<ht->ht_sz;++i)
		{
			//if(i%1000000==0)
			//{cout<<i<<endl;}

			size_t list_sz=0;
			bucket_r1 *bktptr=(bucket_r1 *) ht->store_pos[i];
			bucket_r1 **bktp2p=(bucket_r1 **) &(ht->store_pos[i]);
			while(bktptr!=NULL)
			{
				
				if( ((bucket_r1 *)bktptr)->kmer_info.cov1<=NodeCovTh)
				{
					bucket_r1 * n_bktptr;
			
					n_bktptr=((bucket_r1 *)bktptr)->nxt_bucket;
	
					free(bktptr);
					bktptr=n_bktptr;
					(*bktp2p)=n_bktptr;
					Removed_Nodes_cnt++;
				}
				else
				{
					bktp2p=&(bktptr->nxt_bucket);
					bktptr=bktptr->nxt_bucket;
				}
				
				
			}
		}	
	}
	else
	{
		for(size_t i=0;i<ht2->ht_sz;++i)
		{
			//if(i%1000000==0)
			//{cout<<i<<endl;}

			size_t list_sz=0;
			bucket2_r1 *bktptr=(bucket2_r1 *) ht2->store_pos[i];
			bucket2_r1 **bktp2p=(bucket2_r1 **) &(ht2->store_pos[i]);
			while(bktptr!=NULL)
			{
				
				if( ((bucket2_r1 *)bktptr)->kmer_info.cov1<=NodeCovTh)
				{
	
					bucket2_r1 * n_bktptr;
					n_bktptr=((bucket2_r1 *)bktptr)->nxt_bucket;
					free(bktptr);
					bktptr=n_bktptr;
					(*bktp2p)=n_bktptr;
					Removed_Nodes_cnt++;
				}
				else
				{
					bktp2p=&(bktptr->nxt_bucket);
					bktptr=bktptr->nxt_bucket;
				}
			}
		}
	}

	(*bucket_cnt)-=Removed_Nodes_cnt;

}


void RemovingWeakNodes3_r1(hashtable3 *ht3,int K_size,int NodeCovTh, int64_t *bucket_cnt)
{
	int Removed_Nodes_cnt=0;
	if(K_size<=96&&K_size>64)
	{
		for(size_t i=0;i<ht3->ht_sz;++i)
		{

			size_t list_sz=0;
			bucket3_r1 *bktptr=(bucket3_r1 *) ht3->store_pos[i];
			bucket3_r1 **bktp2p=(bucket3_r1 **) &(ht3->store_pos[i]);
			while(bktptr!=NULL)
			{
				
				if( ((bucket3_r1 *)bktptr)->kmer_info.cov1<=NodeCovTh)
				{
	
					bucket3_r1 * n_bktptr;
					n_bktptr=((bucket3_r1 *)bktptr)->nxt_bucket;
					free(bktptr);
					bktptr=n_bktptr;
					(*bktp2p)=n_bktptr;
					Removed_Nodes_cnt++;
				}
				else
				{
					bktp2p=&(bktptr->nxt_bucket);
					bktptr=bktptr->nxt_bucket;
				}
			}
		}
	}
	(*bucket_cnt)-=Removed_Nodes_cnt;

}


void RemovingWeakNodes4_r1(hashtable4 *ht4,int K_size,int NodeCovTh, int64_t *bucket_cnt)
{
	int Removed_Nodes_cnt=0;
	if(K_size<=128)
	{
		for(size_t i=0;i<ht4->ht_sz;++i)
		{
			//if(i%1000000==0)
			//{cout<<i<<endl;}

			size_t list_sz=0;
			bucket4_r1 *bktptr=(bucket4_r1 *) ht4->store_pos[i];
			bucket4_r1 **bktp2p=(bucket4_r1 **) &(ht4->store_pos[i]);
			while(bktptr!=NULL)
			{
				
				if( ((bucket4_r1 *)bktptr)->kmer_info.cov1<=NodeCovTh)
				{
	
					bucket4_r1 * n_bktptr;
					n_bktptr=((bucket4_r1 *)bktptr)->nxt_bucket;
					free(bktptr);
					bktptr=n_bktptr;
					(*bktp2p)=n_bktptr;
					Removed_Nodes_cnt++;
				}
				else
				{
					bktp2p=&(bktptr->nxt_bucket);
					bktptr=bktptr->nxt_bucket;
				}
			}
		}
	}
	(*bucket_cnt)-=Removed_Nodes_cnt;

}

void RemovingWeakNodes0_r1(hashtable0 *ht0,int K_size,int NodeCovTh, int64_t *bucket_cnt)
{
	int Removed_Nodes_cnt=0;

	for(size_t i=0;i<ht0->ht_sz;++i)
	{


		size_t list_sz=0;
		bucket0_r1 *bktptr=(bucket0_r1 *) ht0->store_pos[i];
		bucket0_r1 **bktp2p=(bucket0_r1 **) &(ht0->store_pos[i]);
		while(bktptr!=NULL)
		{
				
			if( ((bucket0_r1 *)bktptr)->kmer_info.cov1<=NodeCovTh)
			{
	
				bucket0_r1 * n_bktptr;
				n_bktptr=((bucket0_r1 *)bktptr)->nxt_bucket;
				free(bktptr);
				bktptr=n_bktptr;
				(*bktp2p)=n_bktptr;
				Removed_Nodes_cnt++;
			}
			else
			{
				bktp2p=&(bktptr->nxt_bucket);
				bktptr=bktptr->nxt_bucket;
			}
		}
	}
	
	(*bucket_cnt)-=Removed_Nodes_cnt;

}

void MergeNode(hashtable *merge_ht,uint64_t kmer,uint64_t merge_kmer,bool flip)
{
	
	uint64_t hv=MurmurHash64A(&kmer,sizeof(kmer),0);
	uint64_t hash_idx=(size_t) (hv%(merge_ht->ht_sz));

	struct bucket_rm ** ptr;
	ptr=(bucket_rm **) &(merge_ht->store_pos[hash_idx]);
	bool r_found=look_up_in_a_list_rm(kmer,&ptr);
	if(r_found==0)
	{
		(*ptr)=(struct bucket_rm*)malloc(sizeof(struct bucket_rm));
		memset(*ptr,0,sizeof(struct bucket_rm));
		((struct bucket_rm*) *(ptr))->kmer_t.kmer=kmer;
		((struct bucket_rm*) *(ptr))->merged_kmer.kmer=merge_kmer;
		(*(ptr))->flip=flip;
	}
}



void MergeNode2(hashtable2 *merge_ht,kmer_t2 kmer,kmer_t2 merge_kmer,bool flip)
{
	uint64_t hv=MurmurHash64A(&kmer,sizeof(kmer),0);
	uint64_t hash_idx=(size_t) (hv%(merge_ht->ht_sz));

	struct bucket_rm2 ** ptr;
	ptr=(bucket_rm2 **) &(merge_ht->store_pos[hash_idx]);

	bool r_found=look_up_in_a_list_rm2(&kmer,&ptr);

	if(r_found==0)
	{
	
		(*ptr)=(struct bucket_rm2*)malloc(sizeof(struct bucket_rm2));

		memset(*ptr,0,sizeof(struct bucket_rm2));
		(*ptr)->nxt_bucket=NULL;
	
		(*(ptr))->kmer_t2=kmer;
	
		(*(ptr))->merged_kmer=merge_kmer;
		
		(*(ptr))->flip=flip;
	
	}
}



void MergeNode3(hashtable3 *merge_ht,kmer_t3 kmer,kmer_t3 merge_kmer,bool flip)
{
	
	uint64_t hv=MurmurHash64A(&kmer,sizeof(kmer),0);
	uint64_t hash_idx=(size_t) (hv%(merge_ht->ht_sz));

	struct bucket_rm3 ** ptr;
	ptr=(bucket_rm3 **) &(merge_ht->store_pos[hash_idx]);
	bool r_found=look_up_in_a_list_rm3(&kmer,&ptr);
	if(r_found==0)
	{
		(*ptr)=(struct bucket_rm3*)malloc(sizeof(struct bucket_rm3));
		memset(*ptr,0,sizeof(struct bucket_rm3));
		((struct bucket_rm3*) *(ptr))->kmer_t3=kmer;
		((struct bucket_rm3*) *(ptr))->merged_kmer=merge_kmer;
		(*(ptr))->flip=flip;
	}
}


void MergeNode4(hashtable4 *merge_ht,kmer_t4 kmer,kmer_t4 merge_kmer,bool flip)
{
	
	uint64_t hv=MurmurHash64A(&kmer,sizeof(kmer),0);
	uint64_t hash_idx=(size_t) (hv%(merge_ht->ht_sz));

	struct bucket_rm4 ** ptr;
	ptr=(bucket_rm4 **) &(merge_ht->store_pos[hash_idx]);
	bool r_found=look_up_in_a_list_rm4(&kmer,&ptr);
	if(r_found==0)
	{
		(*ptr)=(struct bucket_rm4*)malloc(sizeof(struct bucket_rm4));
		memset(*ptr,0,sizeof(struct bucket_rm4));
		((struct bucket_rm4*) *(ptr))->kmer_t4=kmer;
		((struct bucket_rm4*) *(ptr))->merged_kmer=merge_kmer;
		(*(ptr))->flip=flip;
	}
}

void MergeNode0(hashtable0 *merge_ht,uint64_t *kmer,uint64_t *merge_kmer,bool flip,int Kmer_arr_sz)
{
	
	uint64_t hv=MurmurHash64A(kmer,sizeof(uint64_t)*Kmer_arr_sz,0);
	uint64_t hash_idx=(size_t) (hv%(merge_ht->ht_sz));

	struct bucket_rm0 ** ptr;
	ptr=(bucket_rm0 **) &(merge_ht->store_pos[hash_idx]);
	bool r_found=look_up_in_a_list_rm0(kmer,&ptr,Kmer_arr_sz);
	if(r_found==0)
	{
		(*ptr)=(struct bucket_rm0*)malloc(sizeof(struct bucket_rm0));
		memset(*ptr,0,sizeof( bucket_rm0));
		
		
		((struct bucket_rm0*) *(ptr))->kmer_t=kmer;
		((struct bucket_rm0*) *(ptr))->merged_kmer=merge_kmer;
		(*(ptr))->flip=flip;
	}
}


bool isSimplePath(hashtable *ht,struct bucket* bktptr,map<bucket*,struct BFS_path_info > & Visited_Path , map<struct bucket* ,int > &stacked_nodes)
{
	//return 1;
	struct bucket* end_bkt=bktptr;
	int dep=Visited_Path[bktptr].depth;

	for(int l=dep;l>2;--l)
	{
		struct bucket* p_bktptr=bktptr; //3406059097041372864
	//	if(bktptr->kmer_t.kmer==1027722121587175282)
		//{cout<<"";}

		bktptr=Visited_Path[bktptr].last_bkt;
		if(bktptr==NULL)
		{return 1;}
		
		if(abs(stacked_nodes[bktptr])>2)
		{
			return 0;
		}// complicated bubble.
//		{return 1;}// only backtrack to here so return 1.
		
			
		if(stacked_nodes[bktptr]>0&&(bktptr->kmer_info.left!=NULL&&bktptr->kmer_info.left->nxt_edge!=NULL))
		{return 0;}
		if(stacked_nodes[bktptr]<0&&(bktptr->kmer_info.right!=NULL&&bktptr->kmer_info.right->nxt_edge!=NULL))
		{return 0;}
				
		

	}
	return 1;

}

bool isSimplePath2(hashtable2 *ht,struct bucket2* bktptr,map<bucket2*,struct BFS_path_info2 > & Visited_Path , map<struct bucket2* ,int > &stacked_nodes)
{
	//return 1;
	struct bucket2* end_bkt=bktptr;
	int dep=Visited_Path[bktptr].depth;

	for(int l=dep;l>2;--l)
	{
		struct bucket2* p_bktptr=bktptr; //3406059097041372864
	//	if(bktptr->kmer_t.kmer==1027722121587175282)
		//{cout<<"";}

		bktptr=Visited_Path[bktptr].last_bkt;
		if(bktptr==NULL)
		{return 1;}
		
		if (abs(stacked_nodes[bktptr])>2)
		{
			return 0;
		}// complicated bubble.
		//		{return 1;}// only backtrack to here so return 1.
		
			
		if(stacked_nodes[bktptr]>0&&(bktptr->kmer_info.left!=NULL&&bktptr->kmer_info.left->nxt_edge!=NULL))
		{return 0;}
		if(stacked_nodes[bktptr]<0&&(bktptr->kmer_info.right!=NULL&&bktptr->kmer_info.right->nxt_edge!=NULL))
		{return 0;}
				
		

	}
	return 1;

}
bool isSimplePath3(hashtable3 *ht,struct bucket3* bktptr,map<bucket3*,struct BFS_path_info3 > & Visited_Path , map<struct bucket3* ,int > &stacked_nodes)
{
	//return 1;
	struct bucket3* end_bkt=bktptr;
	int dep=Visited_Path[bktptr].depth;

	for(int l=dep;l>2;--l)
	{
		struct bucket3* p_bktptr=bktptr; //3406059097041372864
	//	if(bktptr->kmer_t.kmer==1027722121587175282)
		//{cout<<"";}

		bktptr=Visited_Path[bktptr].last_bkt;
		if(bktptr==NULL)
		{return 1;}
		
		if (abs(stacked_nodes[bktptr])>2)
		{
			return 0;
		}// complicated bubble.
		//		{return 1;}// only backtrack to here so return 1.
		
			
		if(stacked_nodes[bktptr]>0&&(bktptr->kmer_info.left!=NULL&&bktptr->kmer_info.left->nxt_edge!=NULL))
		{return 0;}
		if(stacked_nodes[bktptr]<0&&(bktptr->kmer_info.right!=NULL&&bktptr->kmer_info.right->nxt_edge!=NULL))
		{return 0;}
				
		

	}
	return 1;

}
bool isSimplePath4(hashtable4 *ht,struct bucket4* bktptr,map<bucket4*,struct BFS_path_info4 > & Visited_Path , map<struct bucket4* ,int > &stacked_nodes)
{
	//return 1;
	struct bucket4* end_bkt=bktptr;
	int dep=Visited_Path[bktptr].depth;

	for(int l=dep;l>2;--l)
	{
		struct bucket4* p_bktptr=bktptr; //3406059097041372864
	//	if(bktptr->kmer_t.kmer==1027722121587175282)
		//{cout<<"";}

		bktptr=Visited_Path[bktptr].last_bkt;
		if(bktptr==NULL)
		{return 1;}
		
		if (abs(stacked_nodes[bktptr])>2)
		{
			return 0;
		}// complicated bubble.
		//		{return 1;}// only backtrack to here so return 1.
		
			
		if(stacked_nodes[bktptr]>0&&(bktptr->kmer_info.left!=NULL&&bktptr->kmer_info.left->nxt_edge!=NULL))
		{return 0;}
		if(stacked_nodes[bktptr]<0&&(bktptr->kmer_info.right!=NULL&&bktptr->kmer_info.right->nxt_edge!=NULL))
		{return 0;}
				
		

	}
	return 1;

}

bool isSimplePath0(hashtable0 *ht,struct bucket0* bktptr,map<bucket0*,struct BFS_path_info0 > & Visited_Path , map<struct bucket0* ,int > &stacked_nodes)
{
	//return 1;
	struct bucket0* end_bkt=bktptr;
	int dep=Visited_Path[bktptr].depth;

	for(int l=dep;l>2;--l)
	{
		struct bucket0* p_bktptr=bktptr; //3406059097041372864
	//	if(bktptr->kmer_t.kmer==1027722121587175282)
		//{cout<<"";}

		bktptr=Visited_Path[bktptr].last_bkt;
		if(bktptr==NULL)
		{return 1;}
		
		if (abs(stacked_nodes[bktptr])>2)
		{
			return 0;
		}// complicated bubble.
		//		{return 1;}// only backtrack to here so return 1.
		
			
		if(stacked_nodes[bktptr]>0&&(bktptr->kmer_info.left!=NULL&&bktptr->kmer_info.left->nxt_edge!=NULL))
		{return 0;}
		if(stacked_nodes[bktptr]<0&&(bktptr->kmer_info.right!=NULL&&bktptr->kmer_info.right->nxt_edge!=NULL))
		{return 0;}
				
		

	}
	return 1;

}

// backtrack operation in the bubble removal
void BacktrackBubbleRemoval(hashtable *ht,struct hashtable* merge_ht,struct bucket* bktptr,struct bucket* bktptr_merged,struct bucket* beg_bkt,map<bucket*,struct BFS_path_info > & Visited_Path , map<struct bucket* ,int > &stacked_nodes, int K_size)
{
	struct bucket* end_bkt=bktptr;
	int dep=Visited_Path[bktptr].depth;

	for(int l=dep;l>1;--l)
	{
		struct bucket* p_bktptr=bktptr; //3406059097041372864
	//	if(bktptr->kmer_t.kmer==1027722121587175282)
	//	{cout<<"";}

		bktptr=Visited_Path[bktptr].last_bkt;
		if(bktptr==NULL)
		{return;}
		
		
		
		if(stacked_nodes[bktptr]>=1||stacked_nodes[bktptr]<=-1)
		{
			uint64_t edge_bits=Visited_Path[p_bktptr].last_bkt_edge->edge;
			int edge_len=(int)(Visited_Path[p_bktptr].last_bkt_edge->len+1);

			if(abs(stacked_nodes[bktptr])>2)
			{
				BreakLinks(stacked_nodes,bktptr,p_bktptr,K_size,edge_len);
				break;
			}
			//else
			
			BreakLinks(stacked_nodes,bktptr,p_bktptr,K_size,edge_len);
			

			if(Visited_Path[bktptr].last_bkt==NULL)
			{
				break;
			}

			struct bucket *freebkt=bktptr;

			if(beg_bkt==freebkt)
			{
				break;
			}

			if(freebkt!=end_bkt)
			{
	//			if(freebkt->kmer_t.kmer==1027722121587175282)
		//		{cout<<"";}
			//	Free_A_Node( ht, freebkt);
		
				stacked_nodes[freebkt]=stacked_nodes[freebkt]/abs(stacked_nodes[freebkt]);
	
				
				freebkt->kmer_info.removed=1;
				
				bool flip_rm=0;
				if(stacked_nodes[freebkt]*stacked_nodes[bktptr_merged]<0)//double check
				{
					flip_rm=1;
				}
				MergeNode(merge_ht,freebkt->kmer_t.kmer, bktptr_merged->kmer_t.kmer,flip_rm);
				
			}


		}


	}

}


//check for self loops in the bubble removal

bool BackCheckLoop(struct bucket* bktptr,struct bucket* end_bkt,map<bucket*,struct BFS_path_info > & Visited_Path )
{

	int dep=Visited_Path[end_bkt].depth;
	bucket *cur_bkt=end_bkt;
	for(int l=dep;l>=1;--l)
	{
		if(bktptr==cur_bkt)
		{
			return 1;
		}
		
		cur_bkt=Visited_Path[cur_bkt].last_bkt;
		
		if(cur_bkt==NULL)
		{break;}
	}

	return 0;

}


void BFSearchBubbleRemoval(struct hashtable* ht,struct hashtable* merge_ht, struct bucket* bktptr,int K_size, int gap,list<struct stacked_bucket>& kmer_stack,int PathCovTh,int max_depth,int PathSim)
{
	map<bucket*,struct BFS_path_info > Visited_Path;
	map<struct bucket* ,int > stacked_nodes;
	struct bucket *beg_bkt= bktptr;
	int max_stack=300;
	int DepthTh=max_depth;//min(300/gap,20);
	int LenTh=300;
	bool RIGHT=0;
	struct stacked_bucket stacked_bkt=kmer_stack.front();

	map<int , list<stacked_bucket> > dist_ctgs;//neighborset
	dist_ctgs[0].push_back(kmer_stack.front());
	int NBs=1;

	int dist_searched=0;

	bucket* new_node=stacked_bkt.bktptr;
	
	uint64_t kmer,f_kmer;
	//char c_str[100],fc_str[100];
	if(stacked_bkt.RightSearch)
	{
		stacked_nodes[new_node]=1;
	}
	else
	{
		stacked_nodes[new_node]=-1;
	}

	Visited_Path[new_node].cov=0;
	Visited_Path[new_node].depth=1;
	Visited_Path[new_node].len=K_size;
	Visited_Path[new_node].last_bkt=NULL;
	Visited_Path[new_node].last_bkt_edge=NULL;

	map<int , list<stacked_bucket> >::iterator NB_it=dist_ctgs.begin();
	while(1)
	{
		NB_it=dist_ctgs.begin();
		if(NB_it==dist_ctgs.end())
		{break;}
		if(NB_it->second.size()==0)
		{dist_ctgs.erase(NB_it->first);continue;}
		//if(kmer_stack.size()>max_stack)
		if(NBs>max_stack)
		{
			break;
		}
		stacked_bkt=NB_it->second.front();

	/*	
		uint64_t kk=stacked_bkt.bktptr->kmer_t.kmer;
		bitsarr2str(&kk,K_size,c_str,1);
		uint64_t f_kk=get_rev_comp_seq(kk,K_size);
		bitsarr2str(&f_kk,K_size,fc_str,1);
		cout<<c_str<<endl;
		string str=c_str;
		cout<<fc_str<<endl;
		string fstr=fc_str;

	*/

		//kmer_stack.pop_front();
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
			int rb=0;
			while(edge_ptr!=NULL)
			{
				rb++;
				edge_ptr=edge_ptr->nxt_edge;
			}
			edge_ptr=(new_node)->kmer_info.right;

			if(stacked_nodes[new_node]==1&&rb>0)
			{
				stacked_nodes[new_node]=1+rb;
			}
			if(rb==0)
			{
				stacked_nodes[new_node]=2;
				struct bucket * freebkt=(new_node);
				if(freebkt==beg_bkt)
				{					
					continue;
				}
				//tip end reached so backtrack to the branching position.
				if(!isSimplePath(ht,new_node,Visited_Path,stacked_nodes))
				{
					continue;
				}
				BacktrackBubbleRemoval(ht,merge_ht,new_node,beg_bkt,beg_bkt,Visited_Path,stacked_nodes,K_size);
				
				edge_ptr=NULL;

				stacked_nodes[new_node]=1;
				freebkt->kmer_info.removed=1;
				//edge_ptr=NULL;
				continue;

			}

			while(edge_ptr!=NULL)
			{

				kmer=(new_node)->kmer_t.kmer;

				f_kmer=kmer;
				f_kmer=get_rev_comp_seq(kmer,K_size);
				

				int edge_len=(int)(edge_ptr->len+1);


				for(int g=edge_len-1;g>=0;--g)
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





				/*

				uint64_t kk=kmer;
				bitsarr2str(&kk,K_size,c_str,1);
				uint64_t f_kk=get_rev_comp_seq(kk,K_size);
				bitsarr2str(&f_kk,K_size,fc_str,1);
				cout<<c_str<<endl;
				string str=c_str;
				cout<<fc_str<<endl;
				string fstr=fc_str;

				*/






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
			//	if(kmer2==14437493651275778)
			//	{cout<<"";}
				r_found=look_up_in_a_list(kmer2,&ptr);


				if(r_found)
				{
					// not in stack
					if(stacked_nodes[*ptr]==0)
					{

						Visited_Path[*ptr].cov=(int)(Visited_Path[new_node].cov+edge_ptr->edge_cov);
						Visited_Path[*ptr].depth=(Visited_Path[new_node].depth+1);
						int cum_len=(int)(Visited_Path[new_node].len+edge_ptr->len+1);
						Visited_Path[*ptr].len=cum_len;//(Visited_Path[new_node].len+edge_ptr->len+1);

						Visited_Path[*ptr].last_bkt=new_node;
						Visited_Path[*ptr].last_bkt_edge=edge_ptr;
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
						//kmer_stack.push_back(stacked_bkt);
						dist_ctgs[cum_len].push_back(stacked_bkt);
						NBs++;

					}
					else
					{
						if((stacked_nodes[*ptr]>0&&r_flip==0)||(stacked_nodes[*ptr]<0&&r_flip==1))
						{
							//backtrack if the same direction is found
							if((Visited_Path[new_node].cov+edge_ptr->edge_cov<=Visited_Path[*ptr].cov)||( BackCheckLoop(*ptr,new_node,Visited_Path)==1))//loop
							{
								//if(0)
									
								if((((Visited_Path[new_node].cov+edge_ptr->edge_cov)/(Visited_Path[new_node].depth))>PathCovTh))
								{
									edge_ptr=edge_ptr->nxt_edge;
									continue;
								}
								if(abs(int((Visited_Path[new_node].len+edge_ptr->len+1)-Visited_Path[*ptr].len))>PathSim)
								{
									edge_ptr=edge_ptr->nxt_edge;
									continue;
								}
								if(!isSimplePath(ht,new_node,Visited_Path,stacked_nodes))
								{
									edge_ptr=edge_ptr->nxt_edge;
									continue;

								}

								//backtrack the current path, common operation in this search
								if(stacked_nodes[new_node]>2)
								{
									edge_ptr=edge_ptr->nxt_edge;
									BreakLinks(stacked_nodes,new_node,(*ptr),K_size,edge_len);
									continue;
								}
								else
								{
									if(stacked_nodes[new_node]<-2)
									{
										edge_ptr=edge_ptr->nxt_edge;
										BreakLinks(stacked_nodes,new_node,(*ptr),K_size,edge_len);
										continue;
									}
									else
									{
										struct bucket *freebkt=(new_node);

										if(freebkt==beg_bkt)
										{
											struct edge_node **edge_p2p=&(new_node->kmer_info.right);
											while(*edge_p2p!=edge_ptr)
											{
												edge_p2p=&((*edge_p2p)->nxt_edge);
											}
											if(*edge_p2p==edge_ptr)
											{
  												struct edge_node* f_edge_ptr=edge_ptr;
												edge_ptr=edge_ptr->nxt_edge;
												*edge_p2p=(*edge_p2p)->nxt_edge;
												free(f_edge_ptr);
											}
													
											continue;
										}
										edge_ptr=edge_ptr->nxt_edge;
										//free the node and edge.
										BreakLinks(stacked_nodes,new_node,(*ptr),K_size,edge_len);
										BacktrackBubbleRemoval(ht,merge_ht,new_node,beg_bkt,beg_bkt,Visited_Path , stacked_nodes,K_size);
										//Free_A_Node( ht, freebkt);

										//if(new_node->kmer_t.kmer==1027722121587175282)
										//{cout<<"";}

										freebkt->kmer_info.removed=1;
										//edge_ptr=NULL;
										continue;
				//
									}
								}
							}
							else
							{
								//backtrack the original path, rare operation in this search, can lead to errors
								//if(0)
								if(Visited_Path[*ptr].depth>1&&(((Visited_Path[*ptr].cov)/(Visited_Path[*ptr].depth-1))>PathCovTh))
								{
									edge_ptr=edge_ptr->nxt_edge;
									continue;
								}
								if(!isSimplePath(ht,*ptr,Visited_Path,stacked_nodes))
								{
									edge_ptr=edge_ptr->nxt_edge;
									continue;

								}
								BacktrackBubbleRemoval(ht,merge_ht,*ptr,beg_bkt,beg_bkt,Visited_Path , stacked_nodes,K_size);
								
								//marginal 
								Visited_Path[*ptr].cov=(int)(Visited_Path[new_node].cov+edge_ptr->edge_cov);
								Visited_Path[*ptr].depth=Visited_Path[new_node].depth+1;
								Visited_Path[*ptr].len=(int)(Visited_Path[new_node].len+edge_ptr->len+1);
								//marginal 

								Visited_Path[*ptr].last_bkt=new_node;
								Visited_Path[*ptr].last_bkt_edge=edge_ptr;
								edge_ptr=edge_ptr->nxt_edge;
								continue;
							}

						}
						else
						{

							//don't do anything,since both strands are visited.

						}


					}

				}





				if(r_found==0||(*ptr)->kmer_info.removed==1)
				{

					//cout<<"WarningR"<<endl;


					struct edge_node **edge_p2p=&(new_node->kmer_info.right);
					while(*edge_p2p!=edge_ptr)
					{
						edge_p2p=&((*edge_p2p)->nxt_edge);
					}
					if(*edge_p2p==edge_ptr)
					{
  						struct edge_node* f_edge_ptr=edge_ptr;
						edge_ptr=edge_ptr->nxt_edge;
						*edge_p2p=(*edge_p2p)->nxt_edge;
						free(f_edge_ptr);
					}

					if(stacked_nodes[new_node]>2)
					{
						stacked_nodes[new_node]--;
						continue;
					}
					if(stacked_nodes[new_node]<-2)
					{
						stacked_nodes[new_node]++;
						continue;
					}

					//else
					struct bucket * freebkt=(new_node);

					if(freebkt==beg_bkt)
					{
						

						continue;
					}
					if(!isSimplePath(ht,new_node,Visited_Path,stacked_nodes))
					{					
						continue;
					}
					BacktrackBubbleRemoval(ht,merge_ht,new_node,beg_bkt,beg_bkt,Visited_Path, stacked_nodes,K_size);

					//Free_A_Node( ht, freebkt);
					//if(new_node->kmer_t.kmer==1027722121587175282)
				//{cout<<"";}
					freebkt->kmer_info.removed=1;
					edge_ptr=NULL;
					continue;

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

			if(stacked_nodes[new_node]==-1&&lb>0)
			{
				stacked_nodes[new_node]=-1-lb;
			}
			if(lb==0)
			{
				stacked_nodes[new_node]=-2;

				struct bucket * freebkt=(new_node);

				if(freebkt==beg_bkt)
				{
							
					continue;
				}
				if(!isSimplePath(ht,new_node,Visited_Path,stacked_nodes))
				{
					continue;

				}
				BacktrackBubbleRemoval(ht,merge_ht,new_node,beg_bkt,beg_bkt,Visited_Path , stacked_nodes,K_size);
				stacked_nodes[new_node]=-1;
				//Free_A_Node( ht, freebkt);
			//	if(new_node->kmer_t.kmer==1027722121587175282)
			//	{cout<<"";}
				freebkt->kmer_info.removed=1;
				edge_ptr=NULL;
				continue;
			}

			while(edge_ptr!=NULL)
			{
				kmer=(new_node)->kmer_t.kmer;
				f_kmer=kmer;
				f_kmer=get_rev_comp_seq(kmer,K_size);
				int edge_len=(edge_ptr->len+1);

				for(int g=0;g<edge_len;++g)
				{
					uint64_t b=(edge_ptr->edge>>2*g)&0x3;

					//kmer.kmer[0]&=(~((uint64_t)0x3<<(K_size-1-32)));
					kmer>>=2;
					//R_shift_NB(kmer.kmer,2,2);
					kmer|=(b<<(2*(K_size-1)));
				}

				bool l_flip=0;
				uint64_t kmer2,f_kmer2;;

				kmer2=kmer;
				//f_kmer2=kmer;

				f_kmer2=get_rev_comp_seq(kmer2,K_size);






				/*

				uint64_t kk=kmer;
				bitsarr2str(&kk,K_size,c_str,1);
				uint64_t f_kk=get_rev_comp_seq(kk,K_size);
				bitsarr2str(&f_kk,K_size,fc_str,1);
				cout<<c_str<<endl;
		string str=c_str;
		cout<<fc_str<<endl;
		string fstr=fc_str;


		*/






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
			//	if(kmer2==14437493651275778)
			//	{cout<<"";}
				l_found=look_up_in_a_list(kmer2,&ptr);
				//626211953659051268
				if(l_found)
				{

					if(stacked_nodes[*ptr]==0)
					{
						Visited_Path[*ptr].cov=(int)(Visited_Path[new_node].cov+edge_ptr->edge_cov);
						Visited_Path[*ptr].depth=(Visited_Path[new_node].depth+1);
						int cum_len=(Visited_Path[new_node].len+edge_ptr->len+1);
						Visited_Path[*ptr].len=cum_len;

						Visited_Path[*ptr].last_bkt=new_node;
						Visited_Path[*ptr].last_bkt_edge=edge_ptr;
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
						if((stacked_nodes[*ptr]>0&&l_flip==1)||(stacked_nodes[*ptr]<0&&l_flip==0))
						{
							//backtrack if the same direction is found
							if((Visited_Path[new_node].cov+edge_ptr->edge_cov<=Visited_Path[*ptr].cov)||( BackCheckLoop(*ptr,new_node,Visited_Path)==1))
							{
								//if(0)
								if((((Visited_Path[new_node].cov+edge_ptr->edge_cov)/(Visited_Path[new_node].depth))>PathCovTh))
								{
									edge_ptr=edge_ptr->nxt_edge;
									continue;
								}
								if(abs(int((Visited_Path[new_node].len+edge_ptr->len+1)-Visited_Path[*ptr].len))>PathSim)
								{
									edge_ptr=edge_ptr->nxt_edge;
									continue;
								}
								if(!isSimplePath(ht,new_node,Visited_Path,stacked_nodes))
								{
									edge_ptr=edge_ptr->nxt_edge;
									continue;

								}

								//backtrack the current path
								if(stacked_nodes[new_node]>2)
								{
									edge_ptr=edge_ptr->nxt_edge;
									BreakLinks(stacked_nodes,new_node,(*ptr),K_size,edge_len);
									continue;
								}
								else
								{
									if(stacked_nodes[new_node]<-2)
									{
										edge_ptr=edge_ptr->nxt_edge;
										BreakLinks(stacked_nodes,new_node,(*ptr),K_size,edge_len);
										continue;
									}
									else
									{
										struct bucket * freebkt=(new_node);
										if(freebkt==beg_bkt)
										{
											struct edge_node **edge_p2p=&(new_node->kmer_info.left);
											while(*edge_p2p!=edge_ptr)
											{
												edge_p2p=&((*edge_p2p)->nxt_edge);
											}
											if(*edge_p2p==edge_ptr)
											{
  												struct edge_node* f_edge_ptr=edge_ptr;
												edge_ptr=edge_ptr->nxt_edge;
												*edge_p2p=(*edge_p2p)->nxt_edge;
												free(f_edge_ptr);
											}
														
											continue;
										}

										edge_ptr=edge_ptr->nxt_edge;
										BreakLinks(stacked_nodes,new_node,(*ptr),K_size,edge_len);
										BacktrackBubbleRemoval(ht,merge_ht,new_node,beg_bkt,beg_bkt,Visited_Path , stacked_nodes,K_size);
										//Free_A_Node( ht, freebkt);
									//	if(new_node->kmer_t.kmer==1027722121587175282)
				//{cout<<"";}
										freebkt->kmer_info.removed=1;
										//edge_ptr=NULL;
										continue;


									}
								}
							}
							else
							{
								//backtrack the original path
								//if(0)
								if(Visited_Path[*ptr].depth>1&&(((Visited_Path[*ptr].cov)/(Visited_Path[*ptr].depth-1))>PathCovTh))
								{
									edge_ptr=edge_ptr->nxt_edge;
									continue;
								}
								if(!isSimplePath(ht,*ptr,Visited_Path,stacked_nodes))
								{
									edge_ptr=edge_ptr->nxt_edge;
									continue;

								}
								BacktrackBubbleRemoval(ht, merge_ht,*ptr,beg_bkt,beg_bkt,Visited_Path , stacked_nodes,K_size);
								Visited_Path[*ptr].cov=(int)(Visited_Path[new_node].cov+edge_ptr->edge_cov);
								Visited_Path[*ptr].depth=Visited_Path[new_node].depth+1;
								Visited_Path[*ptr].len=(Visited_Path[new_node].len+edge_ptr->len+1);

								Visited_Path[*ptr].last_bkt=new_node;
								Visited_Path[*ptr].last_bkt_edge=edge_ptr;

								edge_ptr=edge_ptr->nxt_edge;
								continue;

							}

						}
						else
						{

							//don't do anything,since both strands are visited.

						}



					}

				}





				if(l_found==0||(*ptr)->kmer_info.removed==1)
				{
					//cout<<"WarningL"<<endl;


					struct edge_node **edge_p2p=&(new_node->kmer_info.left);
					while(*edge_p2p!=edge_ptr)
					{
						edge_p2p=&((*edge_p2p)->nxt_edge);
					}
					if(*edge_p2p==edge_ptr)
					{
  						struct edge_node* f_edge_ptr=edge_ptr;
						edge_ptr=edge_ptr->nxt_edge;
						*edge_p2p=(*edge_p2p)->nxt_edge;
						free(f_edge_ptr);
					}

					if(stacked_nodes[new_node]>2)
					{
						stacked_nodes[new_node]--;
						continue;
					}
					if(stacked_nodes[new_node]<-2)
					{
						stacked_nodes[new_node]++;
						continue;
					}
					//else
					struct bucket * freebkt=(new_node);
					if(freebkt==beg_bkt)
					{
						
						continue;
					}
					if(!isSimplePath(ht,new_node,Visited_Path,stacked_nodes))
					{
					
						continue;

					}
					BacktrackBubbleRemoval(ht,merge_ht,new_node,beg_bkt,beg_bkt,Visited_Path , stacked_nodes,K_size);

					//Free_A_Node(ht,freebkt);
				//	if(new_node->kmer_t.kmer==1027722121587175282)
				//{cout<<"";}
					freebkt->kmer_info.removed=1;
					edge_ptr=NULL;
					continue;

				}

				edge_ptr=edge_ptr->nxt_edge;
			}




		}

	}

}


void BacktrackBubbleRemoval2(hashtable2 *ht,struct hashtable2* merge_ht,struct bucket2* bktptr,struct bucket2* bktptr_merged,struct bucket2* beg_bkt,map<bucket2*,struct BFS_path_info2 > & Visited_Path , map<struct bucket2* ,int > &stacked_nodes, int K_size)
{
	struct bucket2* end_bkt=bktptr;
	int dep=Visited_Path[bktptr].depth;

	for(int l=dep;l>1;--l)
	{
		struct bucket2* p_bktptr=bktptr; //3406059097041372864
		//if(bktptr->kmer_t.kmer==766213858554745376)
		//{cout<<"";}

		bktptr=Visited_Path[bktptr].last_bkt;
		if(bktptr==NULL)
		{return;}
		if(stacked_nodes[bktptr]>=1||stacked_nodes[bktptr]<=-1)
		{
			uint64_t edge_bits=Visited_Path[p_bktptr].last_bkt_edge->edge;
			int edge_len=(Visited_Path[p_bktptr].last_bkt_edge->len+1);

			if(abs(stacked_nodes[bktptr])>2)
			{
				BreakLinks2(stacked_nodes,bktptr,p_bktptr,K_size,edge_len);
				break;
			}
			else
			{
				BreakLinks2(stacked_nodes,bktptr,p_bktptr,K_size,edge_len);
			}

			if(Visited_Path[bktptr].last_bkt==NULL)
			{
				break;
			}

			struct bucket2 *freebkt=bktptr;

			if(beg_bkt==freebkt)
			{
				break;
			}

			if(freebkt!=end_bkt)
			{

			//	Free_A_Node( ht, freebkt);
				freebkt->kmer_info.removed=1;
				stacked_nodes[freebkt]=stacked_nodes[freebkt]/abs(stacked_nodes[freebkt]);


				bool flip_rm=0;
				if(stacked_nodes[freebkt]*stacked_nodes[bktptr_merged]<0)//double check
				{
					flip_rm=1;
				}
				MergeNode2(merge_ht,freebkt->kmer_t2, bktptr_merged->kmer_t2,flip_rm);

			}


		}


	}

}


bool BackCheckLoop2(struct bucket2* bktptr,struct bucket2* end_bkt,map<bucket2*,struct BFS_path_info2 > & Visited_Path )
{

	int dep=Visited_Path[end_bkt].depth;
	bucket2 *cur_bkt=end_bkt;
	for(int l=dep;l>=1;--l)
	{
		if(bktptr==cur_bkt)
		{
			return 1;
		}
		
		cur_bkt=Visited_Path[cur_bkt].last_bkt;
		
		if(cur_bkt==NULL)
		{break;}
	}

	return 0;

}

// search operation in the bubble removal

void BFSearchBubbleRemoval2(struct hashtable2* ht,struct hashtable2* merge_ht,struct bucket2* bktptr,int K_size, int gap,list<struct stacked_bucket2>& kmer_stack,int PathCovTh,int max_depth,int PathSim)
{
	map<bucket2*,struct BFS_path_info2 > Visited_Path;
	map<struct bucket2* ,int > stacked_nodes;
	int max_stack=300;
	struct bucket2 *beg_bkt= bktptr;
	int DepthTh=max_depth;
	int LenTh=300;
	bool RIGHT=0;
	struct stacked_bucket2 stacked_bkt=kmer_stack.front();


	map<int , list<stacked_bucket2> > dist_ctgs;//neighborset
	dist_ctgs[0].push_back(kmer_stack.front());
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

	Visited_Path[new_node].cov=0;
	Visited_Path[new_node].depth=1;
	Visited_Path[new_node].len=K_size;
	Visited_Path[new_node].last_bkt=NULL;
	Visited_Path[new_node].last_bkt_edge=NULL;


	map<int , list<stacked_bucket2> >::iterator NB_it=dist_ctgs.begin();
	while(1)
	{
		NB_it=dist_ctgs.begin();
		if(NB_it==dist_ctgs.end())
		{break;}
		if(NB_it->second.size()==0)
		{dist_ctgs.erase(NB_it->first);continue;}
		//if(kmer_stack.size()>max_stack)
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
			int rb=0;
			while(edge_ptr!=NULL)
			{
				rb++;
				edge_ptr=edge_ptr->nxt_edge;
			}
			edge_ptr=(new_node)->kmer_info.right;

			if(stacked_nodes[new_node]==1&&rb>0)
			{
				stacked_nodes[new_node]=1+rb;
			}
			if(rb==0)
			{
				stacked_nodes[new_node]=2;
				struct bucket2 * freebkt=(new_node);
				if(freebkt==beg_bkt)
				{									
					continue;
			 	}
				//tip end reached so backtrack to the branching position.

				if(!isSimplePath2(ht,new_node,Visited_Path,stacked_nodes))
				{
					continue;
				}

				BacktrackBubbleRemoval2(ht,merge_ht,new_node,beg_bkt,beg_bkt,Visited_Path,stacked_nodes,K_size);

				//Free_A_Node(ht,freebkt);
				stacked_nodes[new_node]=1;
				freebkt->kmer_info.removed=1;
				edge_ptr=NULL;
				continue;

			}

			while(edge_ptr!=NULL)
			{


				kmer=(new_node)->kmer_t2;

				f_kmer=kmer;
				get_rev_comp_seq_arr(f_kmer.kmer,K_size,2);
				int edge_len=(edge_ptr->len+1);


				for(int g=edge_len-1;g>=0;--g)
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

				if(r_found)
				{
					// not in stack, put in stack 
					if(stacked_nodes[*ptr]==0)
					{
						Visited_Path[*ptr].cov=(int)(Visited_Path[new_node].cov+edge_ptr->edge_cov);
						Visited_Path[*ptr].depth=(Visited_Path[new_node].depth+1);
						
						int cum_len=(Visited_Path[new_node].len+edge_ptr->len+1);
						Visited_Path[*ptr].len=cum_len;//(Visited_Path[new_node].len+edge_ptr->len+1);


						Visited_Path[*ptr].last_bkt=new_node;
						Visited_Path[*ptr].last_bkt_edge=edge_ptr;
						if(r_flip)
						{
							stacked_nodes[*ptr]=-1;
							stacked_bkt.bktptr=*ptr;
							stacked_bkt.RightSearch=0;
						//	kmer_stack.push_back(stacked_bkt);
						}
						else
						{
							stacked_nodes[*ptr]=1;
							stacked_bkt.bktptr=*ptr;
							stacked_bkt.RightSearch=1;
							//kmer_stack.push_back(stacked_bkt);
						}
						dist_ctgs[cum_len].push_back(stacked_bkt);
						NBs++;

					}
					else
					{
						if((stacked_nodes[*ptr]>0&&r_flip==0)||(stacked_nodes[*ptr]<0&&r_flip==1))
						{
							//backtrack if the same direction is found
							if(((Visited_Path[new_node].cov+edge_ptr->edge_cov)<=Visited_Path[*ptr].cov)||( BackCheckLoop2(*ptr,new_node,Visited_Path)==1))//loop
							{
								//if(0)
								if((((Visited_Path[new_node].cov+edge_ptr->edge_cov)/(Visited_Path[new_node].depth))>PathCovTh))
								{
									edge_ptr=edge_ptr->nxt_edge;
									continue;
								}
								if(abs(int((Visited_Path[new_node].len+edge_ptr->len+1)-Visited_Path[*ptr].len))>PathSim)
								{
									edge_ptr=edge_ptr->nxt_edge;
									continue;
								}
								if(!isSimplePath2(ht,new_node,Visited_Path,stacked_nodes))
								{
									edge_ptr=edge_ptr->nxt_edge;
									continue;

								}

								//backtrack the current path
								if(stacked_nodes[new_node]>2)
								{
									edge_ptr=edge_ptr->nxt_edge;
									BreakLinks2(stacked_nodes,new_node,(*ptr),K_size,edge_len);
									continue;
								}
								else
								{
									if(stacked_nodes[new_node]<-2)
									{
										edge_ptr=edge_ptr->nxt_edge;
										BreakLinks2(stacked_nodes,new_node,(*ptr),K_size,edge_len);
							
										continue;
									}
									else
									{
										struct bucket2 *freebkt=(new_node);

										if(freebkt==beg_bkt)
										{

											struct edge_node **edge_p2p=&(new_node->kmer_info.right);
											while(*edge_p2p!=edge_ptr)
											{
												edge_p2p=&((*edge_p2p)->nxt_edge);
											}
											if(*edge_p2p==edge_ptr)
											{
  												struct edge_node* f_edge_ptr=edge_ptr;
												edge_ptr=edge_ptr->nxt_edge;
												*edge_p2p=(*edge_p2p)->nxt_edge;
												free(f_edge_ptr);
											}
														
											continue;
										}
										edge_ptr=edge_ptr->nxt_edge;
										//free the node and edge.
										
										BreakLinks2(stacked_nodes,new_node,(*ptr),K_size,edge_len);
									
										BacktrackBubbleRemoval2(ht,merge_ht,new_node,beg_bkt,beg_bkt,Visited_Path, stacked_nodes,K_size);
								
										//Free_A_Node( ht, freebkt);
										freebkt->kmer_info.removed=1;
										//edge_ptr=NULL;
										continue;
				//
									}
								}
							}
							else
							{
								//backtrack the original path
								
								if(Visited_Path[*ptr].depth>1&&(((Visited_Path[*ptr].cov)/(Visited_Path[*ptr].depth-1))>PathCovTh))
								{
									edge_ptr=edge_ptr->nxt_edge;
									continue;
								}
								if(!isSimplePath2(ht,*ptr,Visited_Path,stacked_nodes))
								{
									edge_ptr=edge_ptr->nxt_edge;
									continue;

								}
													
								BacktrackBubbleRemoval2(ht,merge_ht,*ptr,beg_bkt,beg_bkt,Visited_Path , stacked_nodes,K_size);
						
								Visited_Path[*ptr].cov=(int)(Visited_Path[new_node].cov+edge_ptr->edge_cov);
								Visited_Path[*ptr].depth=Visited_Path[new_node].depth+1;
								Visited_Path[*ptr].len=(int)(Visited_Path[new_node].len+edge_ptr->len+1);
								Visited_Path[*ptr].last_bkt=new_node;
								Visited_Path[*ptr].last_bkt_edge=edge_ptr;
								edge_ptr=edge_ptr->nxt_edge;
								continue;
							}

						}
						else
						{

							//don't do anything,since both strands are visited.

						}


					}

				}





				if(r_found==0||(*ptr)->kmer_info.removed==1)
				{
					//cout<<"WarningR"<<endl;


					struct edge_node **edge_p2p=&(new_node->kmer_info.right);
					while(*edge_p2p!=edge_ptr)
					{
						edge_p2p=&((*edge_p2p)->nxt_edge);
					}
					if(*edge_p2p==edge_ptr)
					{
  						struct edge_node* f_edge_ptr=edge_ptr;
						edge_ptr=edge_ptr->nxt_edge;
						*edge_p2p=(*edge_p2p)->nxt_edge;
						free(f_edge_ptr);
					}

					if(stacked_nodes[new_node]>2)
					{
						stacked_nodes[new_node]--;
						continue;
					}
					if(stacked_nodes[new_node]<-2)
					{
						stacked_nodes[new_node]++;
						continue;
					}

					//else
					struct bucket2 * freebkt=(new_node);

					if(freebkt==beg_bkt)
					{
			
						continue;
					}
					if(!isSimplePath2(ht,new_node,Visited_Path,stacked_nodes))
					{					
						continue;
					}
							
					BacktrackBubbleRemoval2(ht,merge_ht,new_node,beg_bkt,beg_bkt,Visited_Path, stacked_nodes,K_size);

					freebkt->kmer_info.removed=1;
					edge_ptr=NULL;
					continue;

				}

				edge_ptr=edge_ptr->nxt_edge;
			}


		}
		else
		{

			bool l_found=0;
			struct edge_node *edge_ptr;
//			struct edge_node **edge_p2p;
			edge_ptr=(new_node)->kmer_info.left;


			int lb=0;
			while(edge_ptr!=NULL)
			{
				lb++;
				edge_ptr=edge_ptr->nxt_edge;
			}
			edge_ptr=(new_node)->kmer_info.left;

			if(stacked_nodes[new_node]==-1&&lb>0)
			{
				stacked_nodes[new_node]=-1-lb;
			}
			if(lb==0)
			{
				stacked_nodes[new_node]=-2;

				struct bucket2 * freebkt=(new_node);

				if(freebkt==beg_bkt)
				{
					
							
					continue;
				}
				if(!isSimplePath2(ht,new_node,Visited_Path,stacked_nodes))
				{
					continue;

				}
	
				BacktrackBubbleRemoval2(ht,merge_ht,new_node,beg_bkt,beg_bkt,Visited_Path , stacked_nodes,K_size);

				stacked_nodes[new_node]=-1;
				//Free_A_Node( ht, freebkt);
				freebkt->kmer_info.removed=1;
				edge_ptr=NULL;
				continue;
			}

			while(edge_ptr!=NULL)
			{
				kmer=(new_node)->kmer_t2;
				f_kmer=kmer;
				get_rev_comp_seq_arr(f_kmer.kmer,K_size,2);
				int edge_len=(edge_ptr->len+1);

				for(int g=0;g<edge_len;++g)
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
				//626211953659051268
				if(l_found)
				{

					if(stacked_nodes[*ptr]==0)
					{
						Visited_Path[*ptr].cov=(int)(Visited_Path[new_node].cov+edge_ptr->edge_cov);
						Visited_Path[*ptr].depth=(Visited_Path[new_node].depth+1);
						int cum_len=(Visited_Path[new_node].len+edge_ptr->len+1);
						Visited_Path[*ptr].len=cum_len;
						Visited_Path[*ptr].last_bkt=new_node;
						Visited_Path[*ptr].last_bkt_edge=edge_ptr;
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
						if((stacked_nodes[*ptr]>0&&l_flip==1)||(stacked_nodes[*ptr]<0&&l_flip==0))
						{
							

							//backtrack if the same direction is found
							if((Visited_Path[new_node].cov+edge_ptr->edge_cov<=Visited_Path[*ptr].cov)||(BackCheckLoop2(*ptr,new_node,Visited_Path)==1))
							{
								//if(0)
								if((((Visited_Path[new_node].cov+edge_ptr->edge_cov)/(Visited_Path[new_node].depth))>PathCovTh))
								{
									edge_ptr=edge_ptr->nxt_edge;
									continue;
								}
								if(abs(int((Visited_Path[new_node].len+edge_ptr->len+1)-Visited_Path[*ptr].len))>PathSim)
								{
									edge_ptr=edge_ptr->nxt_edge;
									continue;
								}
								if(!isSimplePath2(ht,new_node,Visited_Path,stacked_nodes))
								{
									edge_ptr=edge_ptr->nxt_edge;
									continue;

								}
								//backtrack the current path
								if(stacked_nodes[new_node]>2)
								{
									edge_ptr=edge_ptr->nxt_edge;
			
									BreakLinks2(stacked_nodes,new_node,(*ptr),K_size,edge_len);
					
									continue;
								}
								else
								{
									if(stacked_nodes[new_node]<-2)
									{
										edge_ptr=edge_ptr->nxt_edge;
						
										BreakLinks2(stacked_nodes,new_node,(*ptr),K_size,edge_len);
						
										continue;
									}
									else
									{
										struct bucket2 * freebkt=(new_node);
										if(freebkt==beg_bkt)
										{
											struct edge_node **edge_p2p=&(new_node->kmer_info.left);
											while(*edge_p2p!=edge_ptr)
											{
												edge_p2p=&((*edge_p2p)->nxt_edge);
											}
											if(*edge_p2p==edge_ptr)
											{
  												struct edge_node* f_edge_ptr=edge_ptr;
												edge_ptr=edge_ptr->nxt_edge;
												*edge_p2p=(*edge_p2p)->nxt_edge;
												free(f_edge_ptr);
											}
											
											continue;
										}

										edge_ptr=edge_ptr->nxt_edge;
						
										BreakLinks2(stacked_nodes,new_node,(*ptr),K_size,edge_len);
			
										BacktrackBubbleRemoval2(ht,merge_ht,new_node,beg_bkt,beg_bkt,Visited_Path , stacked_nodes,K_size);
							
										//Free_A_Node( ht, freebkt);
										freebkt->kmer_info.removed=1;
										//edge_ptr=NULL;
										continue;


									}
								}
							}
							else
							{
								//backtrack the original path
								//if(0)
								if(Visited_Path[*ptr].depth>1&&(((Visited_Path[*ptr].cov)/(Visited_Path[*ptr].depth-1))>PathCovTh))
								{
									edge_ptr=edge_ptr->nxt_edge;
									continue;
								}
								if(!isSimplePath2(ht,*ptr,Visited_Path,stacked_nodes))
								{
									edge_ptr=edge_ptr->nxt_edge;
									continue;

								}
								
								BacktrackBubbleRemoval2(ht,merge_ht, *ptr,beg_bkt,beg_bkt,Visited_Path , stacked_nodes,K_size);
					
								Visited_Path[*ptr].cov=(int)(Visited_Path[new_node].cov+edge_ptr->edge_cov);
								Visited_Path[*ptr].depth=Visited_Path[new_node].depth+1;
								Visited_Path[*ptr].len=(Visited_Path[new_node].len+edge_ptr->len+1);
								Visited_Path[*ptr].last_bkt=new_node;
								Visited_Path[*ptr].last_bkt_edge=edge_ptr;
								edge_ptr=edge_ptr->nxt_edge;
								continue;

							}

						}
						else
						{

							//don't do anything,since both strands are visited.

						}



					}

				}





				if(l_found==0||(*ptr)->kmer_info.removed==1)
				{
					//cout<<"WarningL"<<endl;
					//char c_str[200];
					//bitsarr2str(kmer2.kmer,K_size,c_str,2);
					//cout<<c_str<<endl;
					//bitsarr2str(bktptr->kmer_t2.kmer,K_size,c_str,2);
					//cout<<c_str<<endl;



					struct edge_node **edge_p2p=&(new_node->kmer_info.left);
					while(*edge_p2p!=edge_ptr)
					{
						//uint64_t edge=edge_ptr->edge;
						//bitsarr2str(&edge,16,c_str,1);
						//cout<<c_str<<endl;

						edge_p2p=&((*edge_p2p)->nxt_edge);
					}
					if(*edge_p2p==edge_ptr)
					{
  						struct edge_node* f_edge_ptr=edge_ptr;
						edge_ptr=edge_ptr->nxt_edge;
						*edge_p2p=(*edge_p2p)->nxt_edge;
						free(f_edge_ptr);
					}

					if(stacked_nodes[new_node]>2)
					{
						stacked_nodes[new_node]--;
						continue;
					}
					if(stacked_nodes[new_node]<-2)
					{
						stacked_nodes[new_node]++;
						continue;
					}
					//else
					struct bucket2 * freebkt=(new_node);
					if(freebkt==beg_bkt)
					{
						
							
						continue;
					}

					if(!isSimplePath2(ht,new_node,Visited_Path,stacked_nodes))
					{
						//edge_ptr=edge_ptr->nxt_edge;
						continue;

					}
					BacktrackBubbleRemoval2(ht,merge_ht,new_node,beg_bkt,beg_bkt,Visited_Path , stacked_nodes,K_size);
			

					//Free_A_Node(ht,freebkt);
					freebkt->kmer_info.removed=1;
					edge_ptr=NULL;
					continue;

				}

				edge_ptr=edge_ptr->nxt_edge;
			}




		}

	}

}




void BacktrackBubbleRemoval3(hashtable3 *ht,struct hashtable3* merge_ht, struct bucket3* bktptr,struct bucket3* bktptr_merged,struct bucket3* beg_bkt,map<bucket3*,struct BFS_path_info3 > & Visited_Path , map<struct bucket3* ,int > &stacked_nodes, int K_size)
{
	struct bucket3* end_bkt=bktptr;
	int dep=Visited_Path[bktptr].depth;

	for(int l=dep;l>1;--l)
	{
		struct bucket3* p_bktptr=bktptr; //3406059097041372864
		//if(bktptr->kmer_t.kmer==766213858554745376)
		//{cout<<"";}

		bktptr=Visited_Path[bktptr].last_bkt;
		if(bktptr==NULL)
		{return;}
		if(stacked_nodes[bktptr]>=1||stacked_nodes[bktptr]<=-1)
		{
			uint64_t edge_bits=Visited_Path[p_bktptr].last_bkt_edge->edge;
			int edge_len=(Visited_Path[p_bktptr].last_bkt_edge->len+1);

			if(abs(stacked_nodes[bktptr])>2)
			{
				BreakLinks3(stacked_nodes,bktptr,p_bktptr,K_size,edge_len);
				break;
			}
			else
			{
				BreakLinks3(stacked_nodes,bktptr,p_bktptr,K_size,edge_len);
			}

			if(Visited_Path[bktptr].last_bkt==NULL)
			{
				break;
			}

			struct bucket3 *freebkt=bktptr;

			if(beg_bkt==freebkt)
			{
				break;
			}

			if(freebkt!=end_bkt)
			{

			//	Free_A_Node( ht, freebkt);
				freebkt->kmer_info.removed=1;
				stacked_nodes[freebkt]=stacked_nodes[freebkt]/abs(stacked_nodes[freebkt]);

				
				bool flip_rm=0;
				if(stacked_nodes[freebkt]*stacked_nodes[bktptr_merged]<0)//double check
				{
					flip_rm=1;
				}
				MergeNode3(merge_ht,freebkt->kmer_t3, bktptr_merged->kmer_t3,flip_rm);
			}


		}


	}

}


bool BackCheckLoop3(struct bucket3* bktptr,struct bucket3* end_bkt,map<bucket3*,struct BFS_path_info3 > & Visited_Path )
{

	int dep=Visited_Path[end_bkt].depth;
	bucket3 *cur_bkt=end_bkt;
	for(int l=dep;l>=1;--l)
	{
		if(bktptr==cur_bkt)
		{
			return 1;
		}
		
		cur_bkt=Visited_Path[cur_bkt].last_bkt;
		
		if(cur_bkt==NULL)
		{break;}
	}

	return 0;

}


void BFSearchBubbleRemoval3(struct hashtable3* ht,struct hashtable3* merge_ht,struct bucket3* bktptr,int K_size, int gap,list<struct stacked_bucket3>& kmer_stack,int PathCovTh,int max_depth,int PathSim)
{
	map<bucket3*,struct BFS_path_info3 > Visited_Path;
	map<struct bucket3* ,int > stacked_nodes;
	int max_stack=300;
	struct bucket3 *beg_bkt= bktptr;
	int DepthTh=max_depth;
	int LenTh=300;
	bool RIGHT=0;
	struct stacked_bucket3 stacked_bkt=kmer_stack.front();

	map<int , list<stacked_bucket3> > dist_ctgs;//neighborset
	dist_ctgs[0].push_back(kmer_stack.front());
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

	Visited_Path[new_node].cov=0;
	Visited_Path[new_node].depth=1;
	Visited_Path[new_node].len=K_size;
	Visited_Path[new_node].last_bkt=NULL;
	Visited_Path[new_node].last_bkt_edge=NULL;

	map<int , list<stacked_bucket3> >::iterator NB_it=dist_ctgs.begin();
	while(1)
	{
		NB_it=dist_ctgs.begin();
		if(NB_it==dist_ctgs.end())
		{break;}
		if(NB_it->second.size()==0)
		{dist_ctgs.erase(NB_it->first);continue;}

		if(NBs>max_stack)
		{
			break;
		}
		//stacked_bkt=kmer_stack.front();
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
			int rb=0;
			while(edge_ptr!=NULL)
			{
				rb++;
				edge_ptr=edge_ptr->nxt_edge;
			}
			edge_ptr=(new_node)->kmer_info.right;

			if(stacked_nodes[new_node]==1&&rb>0)
			{
				stacked_nodes[new_node]=1+rb;
			}
			if(rb==0)
			{
				stacked_nodes[new_node]=2;
				struct bucket3 * freebkt=(new_node);
				if(freebkt==beg_bkt)
				{									
					continue;
			 	}
				//tip end reached so backtrack to the branching position.
				if(!isSimplePath3(ht,new_node,Visited_Path,stacked_nodes))
				{
					continue;
				}
				BacktrackBubbleRemoval3(ht,merge_ht,new_node,beg_bkt,beg_bkt,Visited_Path,stacked_nodes,K_size);

				//Free_A_Node(ht,freebkt);
				stacked_nodes[new_node]=1;
				freebkt->kmer_info.removed=1;
				edge_ptr=NULL;
				continue;

			}

			while(edge_ptr!=NULL)
			{


				kmer=(new_node)->kmer_t3;

				f_kmer=kmer;
				get_rev_comp_seq_arr(f_kmer.kmer,K_size,3);
				int edge_len=(edge_ptr->len+1);


				for(int g=edge_len-1;g>=0;--g)
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

				if(r_found)
				{
					// not in stack, put in stack 
					if(stacked_nodes[*ptr]==0)
					{
						Visited_Path[*ptr].cov=(int)(Visited_Path[new_node].cov+edge_ptr->edge_cov);
						Visited_Path[*ptr].depth=(Visited_Path[new_node].depth+1);
						int cum_len=(Visited_Path[new_node].len+edge_ptr->len+1);
						Visited_Path[*ptr].len=cum_len;//(Visited_Path[new_node].len+edge_ptr->len+1);

						
						Visited_Path[*ptr].last_bkt=new_node;
						Visited_Path[*ptr].last_bkt_edge=edge_ptr;
						if(r_flip)
						{
							stacked_nodes[*ptr]=-1;
							stacked_bkt.bktptr=*ptr;
							stacked_bkt.RightSearch=0;
							//kmer_stack.push_back(stacked_bkt);
						}
						else
						{
							stacked_nodes[*ptr]=1;
							stacked_bkt.bktptr=*ptr;
							stacked_bkt.RightSearch=1;
							//kmer_stack.push_back(stacked_bkt);
						}
						dist_ctgs[cum_len].push_back(stacked_bkt);
						NBs++;
					}
					else
					{
						if((stacked_nodes[*ptr]>0&&r_flip==0)||(stacked_nodes[*ptr]<0&&r_flip==1))
						{
							//backtrack if the same direction is found
							if(((Visited_Path[new_node].cov+edge_ptr->edge_cov)<=Visited_Path[*ptr].cov)||( BackCheckLoop3(*ptr,new_node,Visited_Path)==1))//loop
							{
								//if(0)
								if((((Visited_Path[new_node].cov+edge_ptr->edge_cov)/(Visited_Path[new_node].depth))>PathCovTh))
								{
									edge_ptr=edge_ptr->nxt_edge;
									continue;
								}
								if(abs(int((Visited_Path[new_node].len+edge_ptr->len+1)-Visited_Path[*ptr].len))>PathSim)
								{
									edge_ptr=edge_ptr->nxt_edge;
									continue;
								}
								if(!isSimplePath3(ht,new_node,Visited_Path,stacked_nodes))
								{
									edge_ptr=edge_ptr->nxt_edge;
									continue;

								}
								//backtrack the current path
								if(stacked_nodes[new_node]>2)
								{
									edge_ptr=edge_ptr->nxt_edge;
									BreakLinks3(stacked_nodes,new_node,(*ptr),K_size,edge_len);
									continue;
								}
								else
								{
									if(stacked_nodes[new_node]<-2)
									{
										edge_ptr=edge_ptr->nxt_edge;
										BreakLinks3(stacked_nodes,new_node,(*ptr),K_size,edge_len);
							
										continue;
									}
									else
									{
										struct bucket3 *freebkt=(new_node);

										if(freebkt==beg_bkt)
										{

											struct edge_node **edge_p2p=&(new_node->kmer_info.right);
											while(*edge_p2p!=edge_ptr)
											{
												edge_p2p=&((*edge_p2p)->nxt_edge);
											}
											if(*edge_p2p==edge_ptr)
											{
  												struct edge_node* f_edge_ptr=edge_ptr;
												edge_ptr=edge_ptr->nxt_edge;
												*edge_p2p=(*edge_p2p)->nxt_edge;
												free(f_edge_ptr);
											}
														
											continue;
										}
										edge_ptr=edge_ptr->nxt_edge;
										//free the node and edge.
										
										BreakLinks3(stacked_nodes,new_node,(*ptr),K_size,edge_len);
									
										BacktrackBubbleRemoval3(ht,merge_ht,new_node,beg_bkt,beg_bkt,Visited_Path, stacked_nodes,K_size);
								
										//Free_A_Node( ht, freebkt);
										freebkt->kmer_info.removed=1;
										//edge_ptr=NULL;
										continue;
				//
									}
								}
							}
							else
							{
								//backtrack the original path
								
								if(Visited_Path[*ptr].depth>1&&(((Visited_Path[*ptr].cov)/(Visited_Path[*ptr].depth-1))>PathCovTh))
								{
									edge_ptr=edge_ptr->nxt_edge;
									continue;
								}
								if(!isSimplePath3(ht,*ptr,Visited_Path,stacked_nodes))
								{
									edge_ptr=edge_ptr->nxt_edge;
									continue;

								}
													
								BacktrackBubbleRemoval3(ht,merge_ht,*ptr,beg_bkt,beg_bkt,Visited_Path , stacked_nodes,K_size);
						
								Visited_Path[*ptr].cov=(int)(Visited_Path[new_node].cov+edge_ptr->edge_cov);
								Visited_Path[*ptr].depth=Visited_Path[new_node].depth+1;
								Visited_Path[*ptr].len=(Visited_Path[new_node].len+edge_ptr->len+1);
								Visited_Path[*ptr].last_bkt=new_node;
								Visited_Path[*ptr].last_bkt_edge=edge_ptr;
								edge_ptr=edge_ptr->nxt_edge;
								continue;
							}

						}
						else
						{

							//don't do anything,since both strands are visited.

						}


					}

				}





				if(r_found==0||(*ptr)->kmer_info.removed==1)
				{
					//cout<<"WarningR"<<endl;


					struct edge_node **edge_p2p=&(new_node->kmer_info.right);
					while(*edge_p2p!=edge_ptr)
					{
						edge_p2p=&((*edge_p2p)->nxt_edge);
					}
					if(*edge_p2p==edge_ptr)
					{
  						struct edge_node* f_edge_ptr=edge_ptr;
						edge_ptr=edge_ptr->nxt_edge;
						*edge_p2p=(*edge_p2p)->nxt_edge;
						free(f_edge_ptr);
					}

					if(stacked_nodes[new_node]>2)
					{
						stacked_nodes[new_node]--;
						continue;
					}
					if(stacked_nodes[new_node]<-2)
					{
						stacked_nodes[new_node]++;
						continue;
					}

					//else
					struct bucket3 * freebkt=(new_node);

					if(freebkt==beg_bkt)
					{
			
						continue;
					}
					if(!isSimplePath3(ht,new_node,Visited_Path,stacked_nodes))
					{					
						continue;
					}
							
					BacktrackBubbleRemoval3(ht,merge_ht,new_node,beg_bkt,beg_bkt,Visited_Path, stacked_nodes,K_size);

					freebkt->kmer_info.removed=1;
					edge_ptr=NULL;
					continue;

				}

				edge_ptr=edge_ptr->nxt_edge;
			}


		}
		else
		{

			bool l_found=0;
			struct edge_node *edge_ptr;
			//struct edge_node **edge_p2p;
			edge_ptr=(new_node)->kmer_info.left;


			int lb=0;
			while(edge_ptr!=NULL)
			{
				lb++;
				edge_ptr=edge_ptr->nxt_edge;
			}
			edge_ptr=(new_node)->kmer_info.left;

			if(stacked_nodes[new_node]==-1&&lb>0)
			{
				stacked_nodes[new_node]=-1-lb;
			}
			if(lb==0)
			{
				stacked_nodes[new_node]=-2;

				struct bucket3 * freebkt=(new_node);

				if(freebkt==beg_bkt)
				{
					
							
					continue;
				}
				if(!isSimplePath3(ht,new_node,Visited_Path,stacked_nodes))
				{					
					continue;
				}
	
				BacktrackBubbleRemoval3(ht,merge_ht,new_node,beg_bkt,beg_bkt,Visited_Path , stacked_nodes,K_size);

				stacked_nodes[new_node]=-1;
				//Free_A_Node( ht, freebkt);
				freebkt->kmer_info.removed=1;
				edge_ptr=NULL;
				continue;
			}

			while(edge_ptr!=NULL)
			{
				kmer=(new_node)->kmer_t3;
				f_kmer=kmer;
				get_rev_comp_seq_arr(f_kmer.kmer,K_size,3);
				int edge_len=(edge_ptr->len+1);

				for(int g=0;g<edge_len;++g)
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
				//626211953659051268
				if(l_found)
				{

					if(stacked_nodes[*ptr]==0)
					{
						Visited_Path[*ptr].cov=(int)(Visited_Path[new_node].cov+edge_ptr->edge_cov);
						Visited_Path[*ptr].depth=(Visited_Path[new_node].depth+1);
							int cum_len=(Visited_Path[new_node].len+edge_ptr->len+1);
						Visited_Path[*ptr].len=cum_len;//(Visited_Path[new_node].len+edge_ptr->len+1);

						
						Visited_Path[*ptr].last_bkt=new_node;
						Visited_Path[*ptr].last_bkt_edge=edge_ptr;
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
						if((stacked_nodes[*ptr]>0&&l_flip==1)||(stacked_nodes[*ptr]<0&&l_flip==0))
						{
							

							//backtrack if the same direction is found
							if((Visited_Path[new_node].cov+edge_ptr->edge_cov<=Visited_Path[*ptr].cov)||(BackCheckLoop3(*ptr,new_node,Visited_Path)==1))
							{
								//if(0)
								if((((Visited_Path[new_node].cov+edge_ptr->edge_cov)/(Visited_Path[new_node].depth))>PathCovTh))
								{
									edge_ptr=edge_ptr->nxt_edge;
									continue;
								}
								if(abs(int((Visited_Path[new_node].len+edge_ptr->len+1)-Visited_Path[*ptr].len))>PathSim)
								{
									edge_ptr=edge_ptr->nxt_edge;
									continue;
								}
								if(!isSimplePath3(ht,new_node,Visited_Path,stacked_nodes))
								{
									edge_ptr=edge_ptr->nxt_edge;
									continue;

								}
								//backtrack the current path
								if(stacked_nodes[new_node]>2)
								{
									edge_ptr=edge_ptr->nxt_edge;
			
									BreakLinks3(stacked_nodes,new_node,(*ptr),K_size,edge_len);
					
									continue;
								}
								else
								{
									if(stacked_nodes[new_node]<-2)
									{
										edge_ptr=edge_ptr->nxt_edge;
						
										BreakLinks3(stacked_nodes,new_node,(*ptr),K_size,edge_len);
						
										continue;
									}
									else
									{
										struct bucket3 * freebkt=(new_node);
										if(freebkt==beg_bkt)
										{
											struct edge_node **edge_p2p=&(new_node->kmer_info.left);
											while(*edge_p2p!=edge_ptr)
											{
												edge_p2p=&((*edge_p2p)->nxt_edge);
											}
											if(*edge_p2p==edge_ptr)
											{
  												struct edge_node* f_edge_ptr=edge_ptr;
												edge_ptr=edge_ptr->nxt_edge;
												*edge_p2p=(*edge_p2p)->nxt_edge;
												free(f_edge_ptr);
											}
											
											continue;
										}

										edge_ptr=edge_ptr->nxt_edge;
						
										BreakLinks3(stacked_nodes,new_node,(*ptr),K_size,edge_len);
			
										BacktrackBubbleRemoval3(ht,merge_ht,new_node,beg_bkt,beg_bkt,Visited_Path , stacked_nodes,K_size);
							
										//Free_A_Node( ht, freebkt);
										freebkt->kmer_info.removed=1;
										//edge_ptr=NULL;
										continue;


									}
								}
							}
							else
							{
								//backtrack the original path
								//if(0)
								if(Visited_Path[*ptr].depth>1&&(((Visited_Path[*ptr].cov)/(Visited_Path[*ptr].depth-1))>PathCovTh))
								{
									edge_ptr=edge_ptr->nxt_edge;
									continue;
								}
								if(!isSimplePath3(ht,*ptr,Visited_Path,stacked_nodes))
								{
									edge_ptr=edge_ptr->nxt_edge;
									continue;

								}
								BacktrackBubbleRemoval3(ht,merge_ht, *ptr,beg_bkt,beg_bkt,Visited_Path , stacked_nodes,K_size);
					
								Visited_Path[*ptr].cov=(int)(Visited_Path[new_node].cov+edge_ptr->edge_cov);
								Visited_Path[*ptr].depth=Visited_Path[new_node].depth+1;
								Visited_Path[*ptr].len=(Visited_Path[new_node].len+edge_ptr->len+1);
								Visited_Path[*ptr].last_bkt=new_node;
								Visited_Path[*ptr].last_bkt_edge=edge_ptr;
								edge_ptr=edge_ptr->nxt_edge;
								continue;

							}

						}
						else
						{

							//don't do anything,since both strands are visited.

						}



					}

				}





				if(l_found==0||(*ptr)->kmer_info.removed==1)
				{
					//cout<<"WarningL"<<endl;
					//char c_str[200];
					//bitsarr2str(kmer2.kmer,K_size,c_str,2);
					//cout<<c_str<<endl;
					//bitsarr2str(bktptr->kmer_t2.kmer,K_size,c_str,2);
					//cout<<c_str<<endl;



					struct edge_node **edge_p2p=&(new_node->kmer_info.left);
					while(*edge_p2p!=edge_ptr)
					{
						//uint64_t edge=edge_ptr->edge;
						//bitsarr2str(&edge,16,c_str,1);
						//cout<<c_str<<endl;

						edge_p2p=&((*edge_p2p)->nxt_edge);
					}
					if(*edge_p2p==edge_ptr)
					{
  						struct edge_node* f_edge_ptr=edge_ptr;
						edge_ptr=edge_ptr->nxt_edge;
						*edge_p2p=(*edge_p2p)->nxt_edge;
						free(f_edge_ptr);
					}

					if(stacked_nodes[new_node]>2)
					{
						stacked_nodes[new_node]--;
						continue;
					}
					if(stacked_nodes[new_node]<-2)
					{
						stacked_nodes[new_node]++;
						continue;
					}
					//else
					struct bucket3 * freebkt=(new_node);
					if(freebkt==beg_bkt)
					{
						
							
						continue;
					}
					if(!isSimplePath3(ht,new_node,Visited_Path,stacked_nodes))
					{					
						continue;
					}
					BacktrackBubbleRemoval3(ht,merge_ht,new_node,beg_bkt,beg_bkt,Visited_Path , stacked_nodes,K_size);
			

					//Free_A_Node(ht,freebkt);
					freebkt->kmer_info.removed=1;
					edge_ptr=NULL;
					continue;

				}

				edge_ptr=edge_ptr->nxt_edge;
			}




		}

	}

}


void BacktrackBubbleRemoval4(hashtable4 *ht,struct hashtable4* merge_ht,struct bucket4* bktptr,struct bucket4* bktptr_merged,struct bucket4* beg_bkt,map<bucket4*,struct BFS_path_info4 > & Visited_Path , map<struct bucket4* ,int > &stacked_nodes, int K_size)
{
	struct bucket4* end_bkt=bktptr;
	int dep=Visited_Path[bktptr].depth;

	for(int l=dep;l>1;--l)
	{
		struct bucket4* p_bktptr=bktptr; //3406059097041372864
		//if(bktptr->kmer_t.kmer==766213858554745376)
		//{cout<<"";}

		bktptr=Visited_Path[bktptr].last_bkt;
		if(bktptr==NULL)
		{return;}
		if(stacked_nodes[bktptr]>=1||stacked_nodes[bktptr]<=-1)
		{
			uint64_t edge_bits=Visited_Path[p_bktptr].last_bkt_edge->edge;
			int edge_len=(Visited_Path[p_bktptr].last_bkt_edge->len+1);

			if(abs(stacked_nodes[bktptr])>2)
			{
				BreakLinks4(stacked_nodes,bktptr,p_bktptr,K_size,edge_len);
				break;
			}
			else
			{
				BreakLinks4(stacked_nodes,bktptr,p_bktptr,K_size,edge_len);
			}

			if(Visited_Path[bktptr].last_bkt==NULL)
			{
				break;
			}

			struct bucket4 *freebkt=bktptr;

			if(beg_bkt==freebkt)
			{
				break;
			}

			if(freebkt!=end_bkt)
			{

			//	Free_A_Node( ht, freebkt);
				freebkt->kmer_info.removed=1;
				stacked_nodes[freebkt]=stacked_nodes[freebkt]/abs(stacked_nodes[freebkt]);
				
				bool flip_rm=0;
				if(stacked_nodes[freebkt]*stacked_nodes[bktptr_merged]<0)//double check
				{
					flip_rm=1;
				}
				MergeNode4(merge_ht,freebkt->kmer_t4, bktptr_merged->kmer_t4,flip_rm);
			}


		}


	}

}


bool BackCheckLoop4(struct bucket4* bktptr,struct bucket4* end_bkt,map<bucket4*,struct BFS_path_info4 > & Visited_Path )
{

	int dep=Visited_Path[end_bkt].depth;
	bucket4 *cur_bkt=end_bkt;
	for(int l=dep;l>=1;--l)
	{
		if(bktptr==cur_bkt)
		{
			return 1;
		}
		
		cur_bkt=Visited_Path[cur_bkt].last_bkt;
		
		if(cur_bkt==NULL)
		{break;}
	}

	return 0;

}


void BFSearchBubbleRemoval4(struct hashtable4* ht,struct hashtable4* merge_ht,struct bucket4* bktptr,int K_size, int gap,list<struct stacked_bucket4>& kmer_stack,int PathCovTh,int max_depth,int PathSim)
{
	map<bucket4*,struct BFS_path_info4 > Visited_Path;
	map<struct bucket4* ,int > stacked_nodes;
	int max_stack=300;
	struct bucket4 *beg_bkt= bktptr;
	int DepthTh=max_depth;
	int LenTh=300;
	bool RIGHT=0;
	struct stacked_bucket4 stacked_bkt=kmer_stack.front();
	
	map<int , list<stacked_bucket4> > dist_ctgs;//neighborset
	dist_ctgs[0].push_back(kmer_stack.front());
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

	Visited_Path[new_node].cov=0;
	Visited_Path[new_node].depth=1;
	Visited_Path[new_node].len=K_size;
	Visited_Path[new_node].last_bkt=NULL;
	Visited_Path[new_node].last_bkt_edge=NULL;

	map<int , list<stacked_bucket4> >::iterator NB_it=dist_ctgs.begin();
	while(1)
	{
		NB_it=dist_ctgs.begin();
		if(NB_it==dist_ctgs.end())
		{break;}
		if(NB_it->second.size()==0)
		{dist_ctgs.erase(NB_it->first);continue;}
		//if(kmer_stack.size()>max_stack)
		if(NBs>max_stack)
		{
			break;
		}
		//stacked_bkt=kmer_stack.front();
		stacked_bkt=NB_it->second.front();

		

		//kmer_stack.pop_front();
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
			int rb=0;
			while(edge_ptr!=NULL)
			{
				rb++;
				edge_ptr=edge_ptr->nxt_edge;
			}
			edge_ptr=(new_node)->kmer_info.right;

			if(stacked_nodes[new_node]==1&&rb>0)
			{
				stacked_nodes[new_node]=1+rb;
			}
			if(rb==0)
			{
				stacked_nodes[new_node]=2;
				struct bucket4 * freebkt=(new_node);
				if(freebkt==beg_bkt)
				{									
					continue;
			 	}
				//tip end reached so backtrack to the branching position.
				if(!isSimplePath4(ht,new_node,Visited_Path,stacked_nodes))
				{
					continue;
				}
				BacktrackBubbleRemoval4(ht,merge_ht,new_node,beg_bkt,beg_bkt,Visited_Path,stacked_nodes,K_size);

				//Free_A_Node(ht,freebkt);
				stacked_nodes[new_node]=1;
				freebkt->kmer_info.removed=1;
				edge_ptr=NULL;
				continue;

			}

			while(edge_ptr!=NULL)
			{


				kmer=(new_node)->kmer_t4;

				f_kmer=kmer;
				get_rev_comp_seq_arr(f_kmer.kmer,K_size,4);
				int edge_len=(edge_ptr->len+1);


				for(int g=edge_len-1;g>=0;--g)
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

				if(r_found)
				{
					// not in stack, put in stack 
					if(stacked_nodes[*ptr]==0)
					{
						Visited_Path[*ptr].cov=(int)(Visited_Path[new_node].cov+edge_ptr->edge_cov);
						Visited_Path[*ptr].depth=(Visited_Path[new_node].depth+1);
						int cum_len=(Visited_Path[new_node].len+edge_ptr->len+1);
						Visited_Path[*ptr].len=cum_len;//(Visited_Path[new_node].len+edge_ptr->len+1);

						Visited_Path[*ptr].last_bkt=new_node;
						Visited_Path[*ptr].last_bkt_edge=edge_ptr;
						if(r_flip)
						{
							stacked_nodes[*ptr]=-1;
							stacked_bkt.bktptr=*ptr;
							stacked_bkt.RightSearch=0;
							//kmer_stack.push_back(stacked_bkt);
						}
						else
						{
							stacked_nodes[*ptr]=1;
							stacked_bkt.bktptr=*ptr;
							stacked_bkt.RightSearch=1;
							//kmer_stack.push_back(stacked_bkt);
						}
						dist_ctgs[cum_len].push_back(stacked_bkt);
						NBs++;
					}
					else
					{
						if((stacked_nodes[*ptr]>0&&r_flip==0)||(stacked_nodes[*ptr]<0&&r_flip==1))
						{
							//backtrack if the same direction is found
							if(((Visited_Path[new_node].cov+edge_ptr->edge_cov)<=Visited_Path[*ptr].cov)||( BackCheckLoop4(*ptr,new_node,Visited_Path)==1))//loop
							{
								//if(0)
								if((((Visited_Path[new_node].cov+edge_ptr->edge_cov)/(Visited_Path[new_node].depth))>PathCovTh))
								{
									edge_ptr=edge_ptr->nxt_edge;
									continue;
								}
								if(abs(int((Visited_Path[new_node].len+edge_ptr->len+1)-Visited_Path[*ptr].len))>PathSim)
								{
									edge_ptr=edge_ptr->nxt_edge;
									continue;
								}
								if(!isSimplePath4(ht,new_node,Visited_Path,stacked_nodes))
								{
									edge_ptr=edge_ptr->nxt_edge;
									continue;

								}
								//backtrack the current path
								if(stacked_nodes[new_node]>2)
								{
									edge_ptr=edge_ptr->nxt_edge;
									BreakLinks4(stacked_nodes,new_node,(*ptr),K_size,edge_len);
									continue;
								}
								else
								{
									if(stacked_nodes[new_node]<-2)
									{
										edge_ptr=edge_ptr->nxt_edge;
										BreakLinks4(stacked_nodes,new_node,(*ptr),K_size,edge_len);
							
										continue;
									}
									else
									{
										struct bucket4 *freebkt=(new_node);

										if(freebkt==beg_bkt)
										{

											struct edge_node **edge_p2p=&(new_node->kmer_info.right);
											while(*edge_p2p!=edge_ptr)
											{
												edge_p2p=&((*edge_p2p)->nxt_edge);
											}
											if(*edge_p2p==edge_ptr)
											{
  												struct edge_node* f_edge_ptr=edge_ptr;
												edge_ptr=edge_ptr->nxt_edge;
												*edge_p2p=(*edge_p2p)->nxt_edge;
												free(f_edge_ptr);
											}
														
											continue;
										}
										edge_ptr=edge_ptr->nxt_edge;
										//free the node and edge.
										
										BreakLinks4(stacked_nodes,new_node,(*ptr),K_size,edge_len);
									
										BacktrackBubbleRemoval4(ht,merge_ht,new_node,beg_bkt,beg_bkt,Visited_Path, stacked_nodes,K_size);
								
										//Free_A_Node( ht, freebkt);
										freebkt->kmer_info.removed=1;
										//edge_ptr=NULL;
										continue;
				//
									}
								}
							}
							else
							{
								//backtrack the original path
								
								if(Visited_Path[*ptr].depth>1&&(((Visited_Path[*ptr].cov)/(Visited_Path[*ptr].depth-1))>PathCovTh))
								{
									edge_ptr=edge_ptr->nxt_edge;
									continue;
								}
								if(!isSimplePath4(ht,*ptr,Visited_Path,stacked_nodes))
								{
									edge_ptr=edge_ptr->nxt_edge;
									continue;

								}
													
								BacktrackBubbleRemoval4(ht,merge_ht,*ptr,beg_bkt,beg_bkt,Visited_Path , stacked_nodes,K_size);
						
								Visited_Path[*ptr].cov=(int)(Visited_Path[new_node].cov+edge_ptr->edge_cov);
								Visited_Path[*ptr].depth=Visited_Path[new_node].depth+1;
								Visited_Path[*ptr].len=(Visited_Path[new_node].len+edge_ptr->len+1);
								Visited_Path[*ptr].last_bkt=new_node;
								Visited_Path[*ptr].last_bkt_edge=edge_ptr;
								edge_ptr=edge_ptr->nxt_edge;
								continue;
							}

						}
						else
						{

							//don't do anything,since both strands are visited.

						}


					}

				}





				if(r_found==0||(*ptr)->kmer_info.removed==1)
				{
					//cout<<"WarningR"<<endl;


					struct edge_node **edge_p2p=&(new_node->kmer_info.right);
					while(*edge_p2p!=edge_ptr)
					{
						edge_p2p=&((*edge_p2p)->nxt_edge);
					}
					if(*edge_p2p==edge_ptr)
					{
  						struct edge_node* f_edge_ptr=edge_ptr;
						edge_ptr=edge_ptr->nxt_edge;
						*edge_p2p=(*edge_p2p)->nxt_edge;
						free(f_edge_ptr);
					}

					if(stacked_nodes[new_node]>2)
					{
						stacked_nodes[new_node]--;
						continue;
					}
					if(stacked_nodes[new_node]<-2)
					{
						stacked_nodes[new_node]++;
						continue;
					}

					//else
					struct bucket4 * freebkt=(new_node);

					if(freebkt==beg_bkt)
					{
			
						continue;
					}
					if(!isSimplePath4(ht,new_node,Visited_Path,stacked_nodes))
					{					
						continue;
					}
							
					BacktrackBubbleRemoval4(ht,merge_ht,new_node,beg_bkt,beg_bkt,Visited_Path, stacked_nodes,K_size);

					freebkt->kmer_info.removed=1;
					edge_ptr=NULL;
					continue;

				}

				edge_ptr=edge_ptr->nxt_edge;
			}


		}
		else
		{

			bool l_found=0;
			struct edge_node *edge_ptr;
//			struct edge_node **edge_p2p;
			edge_ptr=(new_node)->kmer_info.left;


			int lb=0;
			while(edge_ptr!=NULL)
			{
				lb++;
				edge_ptr=edge_ptr->nxt_edge;
			}
			edge_ptr=(new_node)->kmer_info.left;

			if(stacked_nodes[new_node]==-1&&lb>0)
			{
				stacked_nodes[new_node]=-1-lb;
			}
			if(lb==0)
			{
				stacked_nodes[new_node]=-2;

				struct bucket4 * freebkt=(new_node);

				if(freebkt==beg_bkt)
				{
					
							
					continue;
				}
	
				if(!isSimplePath4(ht,new_node,Visited_Path,stacked_nodes))
				{					
					continue;
				}
				BacktrackBubbleRemoval4(ht,merge_ht,new_node,beg_bkt,beg_bkt,Visited_Path , stacked_nodes,K_size);

				stacked_nodes[new_node]=-1;
				//Free_A_Node( ht, freebkt);
				freebkt->kmer_info.removed=1;
				edge_ptr=NULL;
				continue;
			}

			while(edge_ptr!=NULL)
			{
				kmer=(new_node)->kmer_t4;
				f_kmer=kmer;
				get_rev_comp_seq_arr(f_kmer.kmer,K_size,4);
				int edge_len=(edge_ptr->len+1);

				for(int g=0;g<edge_len;++g)
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
				//626211953659051268
				if(l_found)
				{

					if(stacked_nodes[*ptr]==0)
					{
						Visited_Path[*ptr].cov=(int)(Visited_Path[new_node].cov+edge_ptr->edge_cov);
						Visited_Path[*ptr].depth=(Visited_Path[new_node].depth+1);
						int cum_len=(int)(Visited_Path[new_node].len+edge_ptr->len+1);
						Visited_Path[*ptr].len=cum_len;//(Visited_Path[new_node].len+edge_ptr->len+1);

						Visited_Path[*ptr].last_bkt=new_node;
						Visited_Path[*ptr].last_bkt_edge=edge_ptr;
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
						if((stacked_nodes[*ptr]>0&&l_flip==1)||(stacked_nodes[*ptr]<0&&l_flip==0))
						{
							

							//backtrack if the same direction is found
							if((Visited_Path[new_node].cov+edge_ptr->edge_cov<=Visited_Path[*ptr].cov)||(BackCheckLoop4(*ptr,new_node,Visited_Path)==1))
							{
								//if(0)
								if((((Visited_Path[new_node].cov+edge_ptr->edge_cov)/(Visited_Path[new_node].depth))>PathCovTh))
								{
									edge_ptr=edge_ptr->nxt_edge;
									continue;
								}
								if(abs(int((Visited_Path[new_node].len+edge_ptr->len+1)-Visited_Path[*ptr].len))>PathSim)
								{
									edge_ptr=edge_ptr->nxt_edge;
									continue;
								}
								if(!isSimplePath4(ht,new_node,Visited_Path,stacked_nodes))
								{
									edge_ptr=edge_ptr->nxt_edge;
									continue;

								}
								//backtrack the current path
								if(stacked_nodes[new_node]>2)
								{
									edge_ptr=edge_ptr->nxt_edge;
			
									BreakLinks4(stacked_nodes,new_node,(*ptr),K_size,edge_len);
					
									continue;
								}
								else
								{
									if(stacked_nodes[new_node]<-2)
									{
										edge_ptr=edge_ptr->nxt_edge;
						
										BreakLinks4(stacked_nodes,new_node,(*ptr),K_size,edge_len);
						
										continue;
									}
									else
									{
										struct bucket4 * freebkt=(new_node);
										if(freebkt==beg_bkt)
										{
											struct edge_node **edge_p2p=&(new_node->kmer_info.left);
											while(*edge_p2p!=edge_ptr)
											{
												edge_p2p=&((*edge_p2p)->nxt_edge);
											}
											if(*edge_p2p==edge_ptr)
											{
  												struct edge_node* f_edge_ptr=edge_ptr;
												edge_ptr=edge_ptr->nxt_edge;
												*edge_p2p=(*edge_p2p)->nxt_edge;
												free(f_edge_ptr);
											}
											
											continue;
										}

										edge_ptr=edge_ptr->nxt_edge;
						
										BreakLinks4(stacked_nodes,new_node,(*ptr),K_size,edge_len);
			
										BacktrackBubbleRemoval4(ht,merge_ht,new_node,beg_bkt,beg_bkt,Visited_Path , stacked_nodes,K_size);
							
										//Free_A_Node( ht, freebkt);
										freebkt->kmer_info.removed=1;
										//edge_ptr=NULL;
										continue;


									}
								}
							}
							else
							{
								//backtrack the original path
								//if(0)
								if(Visited_Path[*ptr].depth>1&&(((Visited_Path[*ptr].cov)/(Visited_Path[*ptr].depth-1))>PathCovTh))
								{
									edge_ptr=edge_ptr->nxt_edge;
									continue;
								}
								if(!isSimplePath4(ht,*ptr,Visited_Path,stacked_nodes))
								{
									edge_ptr=edge_ptr->nxt_edge;
									continue;

								}
								
								BacktrackBubbleRemoval4(ht,merge_ht, *ptr,beg_bkt,beg_bkt,Visited_Path , stacked_nodes,K_size);
					
								Visited_Path[*ptr].cov=(int)(Visited_Path[new_node].cov+edge_ptr->edge_cov);
								Visited_Path[*ptr].depth=Visited_Path[new_node].depth+1;
								Visited_Path[*ptr].len=(int)(Visited_Path[new_node].len+edge_ptr->len+1);
								Visited_Path[*ptr].last_bkt=new_node;
								Visited_Path[*ptr].last_bkt_edge=edge_ptr;
								edge_ptr=edge_ptr->nxt_edge;
								continue;

							}

						}
						else
						{

							//don't do anything,since both strands are visited.

						}



					}

				}





				if(l_found==0||(*ptr)->kmer_info.removed==1)
				{
					//cout<<"WarningL"<<endl;
					//char c_str[200];
					//bitsarr2str(kmer2.kmer,K_size,c_str,2);
					//cout<<c_str<<endl;
					//bitsarr2str(bktptr->kmer_t2.kmer,K_size,c_str,2);
					//cout<<c_str<<endl;



					struct edge_node **edge_p2p=&(new_node->kmer_info.left);
					while(*edge_p2p!=edge_ptr)
					{
						//uint64_t edge=edge_ptr->edge;
						//bitsarr2str(&edge,16,c_str,1);
						//cout<<c_str<<endl;

						edge_p2p=&((*edge_p2p)->nxt_edge);
					}
					if(*edge_p2p==edge_ptr)
					{
  						struct edge_node* f_edge_ptr=edge_ptr;
						edge_ptr=edge_ptr->nxt_edge;
						*edge_p2p=(*edge_p2p)->nxt_edge;
						free(f_edge_ptr);
					}

					if(stacked_nodes[new_node]>2)
					{
						stacked_nodes[new_node]--;
						continue;
					}
					if(stacked_nodes[new_node]<-2)
					{
						stacked_nodes[new_node]++;
						continue;
					}
					//else
					struct bucket4 * freebkt=(new_node);
					if(freebkt==beg_bkt)
					{
						
							
						continue;
					}
					if(!isSimplePath4(ht,new_node,Visited_Path,stacked_nodes))
					{					
						continue;
					}

					BacktrackBubbleRemoval4(ht,merge_ht,new_node,beg_bkt,beg_bkt,Visited_Path , stacked_nodes,K_size);
			

					//Free_A_Node(ht,freebkt);
					freebkt->kmer_info.removed=1;
					edge_ptr=NULL;
					continue;

				}

				edge_ptr=edge_ptr->nxt_edge;
			}




		}

	}

}








void BacktrackBubbleRemoval0(hashtable0 *ht,struct hashtable0* merge_ht,struct bucket0* bktptr,struct bucket0* bktptr_merged,struct bucket0* beg_bkt,map<bucket0*,struct BFS_path_info0 > & Visited_Path , map<struct bucket0* ,int > &stacked_nodes, int K_size)
{
	struct bucket0* end_bkt=bktptr;
	int dep=Visited_Path[bktptr].depth;
	int Kmer_arr_sz=K_size/32+1;
	int rem1=K_size%32;
	if(rem1==0)
	{Kmer_arr_sz--;}
	for(int l=dep;l>1;--l)
	{
		struct bucket0* p_bktptr=bktptr; //3406059097041372864
		//if(bktptr->kmer_t.kmer==766213858554745376)
		//{cout<<"";}

		bktptr=Visited_Path[bktptr].last_bkt;
		if(bktptr==NULL)
		{return;}
		if(stacked_nodes[bktptr]>=1||stacked_nodes[bktptr]<=-1)
		{
			uint64_t edge_bits=Visited_Path[p_bktptr].last_bkt_edge->edge;
			int edge_len=(Visited_Path[p_bktptr].last_bkt_edge->len+1);

			if(abs(stacked_nodes[bktptr])>2)
			{
				BreakLinks0(stacked_nodes,bktptr,p_bktptr,K_size,edge_len);
				break;
			}
			else
			{
				BreakLinks0(stacked_nodes,bktptr,p_bktptr,K_size,edge_len);
			}

			if(Visited_Path[bktptr].last_bkt==NULL)
			{
				break;
			}

			struct bucket0 *freebkt=bktptr;

			if(beg_bkt==freebkt)
			{
				break;
			}

			if(freebkt!=end_bkt)
			{

			//	Free_A_Node( ht, freebkt);
				freebkt->kmer_info.removed=1;
				stacked_nodes[freebkt]=stacked_nodes[freebkt]/abs(stacked_nodes[freebkt]);
				
				bool flip_rm=0;
				if(stacked_nodes[freebkt]*stacked_nodes[bktptr_merged]<0)//double check
				{
					flip_rm=1;
				}
				MergeNode0(merge_ht,freebkt->kmer_t, bktptr_merged->kmer_t,flip_rm,Kmer_arr_sz);
			}


		}


	}

}


bool BackCheckLoop0(struct bucket0* bktptr,struct bucket0* end_bkt,map<bucket0*,struct BFS_path_info0 > & Visited_Path )
{

	int dep=Visited_Path[end_bkt].depth;
	bucket0 *cur_bkt=end_bkt;
	for(int l=dep;l>=1;--l)
	{
		if(bktptr==cur_bkt)
		{
			return 1;
		}
		
		cur_bkt=Visited_Path[cur_bkt].last_bkt;
		
		if(cur_bkt==NULL)
		{break;}
	}

	return 0;

}


void BFSearchBubbleRemoval0(struct hashtable0* ht,struct hashtable0* merge_ht,struct bucket0* bktptr,int K_size, int gap,list<struct stacked_bucket0>& kmer_stack,int PathCovTh,int max_depth,int PathSim)
{
	int Kmer_arr_sz=K_size/32+1;
	int rem1=K_size%32;
	if(rem1==0)
	{Kmer_arr_sz--;}
	map<bucket0*,struct BFS_path_info0 > Visited_Path;
	map<struct bucket0* ,int > stacked_nodes;
	int max_stack=300;
	struct bucket0 *beg_bkt= bktptr;
	int DepthTh=max_depth;
	int LenTh=300;
	bool RIGHT=0;
	struct stacked_bucket0 stacked_bkt=kmer_stack.front();
	
	map<int , list<stacked_bucket0> > dist_ctgs;//neighborset
	dist_ctgs[0].push_back(kmer_stack.front());
	int NBs=1;

	int dist_searched=0;


	bucket0* new_node=stacked_bkt.bktptr;
	uint64_t kmer[100],f_kmer[100];

	if(stacked_bkt.RightSearch)
	{
		stacked_nodes[new_node]=1;
	}
	else
	{
		stacked_nodes[new_node]=-1;
	}

	Visited_Path[new_node].cov=0;
	Visited_Path[new_node].depth=1;
	Visited_Path[new_node].len=K_size;
	Visited_Path[new_node].last_bkt=NULL;
	Visited_Path[new_node].last_bkt_edge=NULL;

	map<int , list<stacked_bucket0> >::iterator NB_it=dist_ctgs.begin();
	while(1)
	{
		NB_it=dist_ctgs.begin();
		if(NB_it==dist_ctgs.end())
		{break;}
		if(NB_it->second.size()==0)
		{dist_ctgs.erase(NB_it->first);continue;}
		//if(kmer_stack.size()>max_stack)
		if(NBs>max_stack)
		{
			break;
		}
		//stacked_bkt=kmer_stack.front();
		stacked_bkt=NB_it->second.front();

		

		//kmer_stack.pop_front();
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
			int rb=0;
			while(edge_ptr!=NULL)
			{
				rb++;
				edge_ptr=edge_ptr->nxt_edge;
			}
			edge_ptr=(new_node)->kmer_info.right;

			if(stacked_nodes[new_node]==1&&rb>0)
			{
				stacked_nodes[new_node]=1+rb;
			}
			if(rb==0)
			{
				stacked_nodes[new_node]=2;
				struct bucket0 * freebkt=(new_node);
				if(freebkt==beg_bkt)
				{									
					continue;
			 	}
				//tip end reached so backtrack to the branching position.
				if(!isSimplePath0(ht,new_node,Visited_Path,stacked_nodes))
				{
					continue;
				}
				BacktrackBubbleRemoval0(ht,merge_ht,new_node,beg_bkt,beg_bkt,Visited_Path,stacked_nodes,K_size);

				//Free_A_Node(ht,freebkt);
				stacked_nodes[new_node]=1;
				freebkt->kmer_info.removed=1;
				edge_ptr=NULL;
				continue;

			}

			while(edge_ptr!=NULL)
			{


				memcpy(kmer,(new_node)->kmer_t,sizeof(uint64_t)*Kmer_arr_sz);
				memcpy(f_kmer,kmer,sizeof(uint64_t)*Kmer_arr_sz);
				get_rev_comp_seq_arr(f_kmer,K_size,Kmer_arr_sz);

				int edge_len=(edge_ptr->len+1);


				for(int g=edge_len-1;g>=0;--g)
				{
					uint64_t b=(edge_ptr->edge>>2*g)&0x3;
					kmer[0]&=(~((uint64_t)0x3<<(2*((K_size-1)%32))));
					L_shift_NB(kmer,2,Kmer_arr_sz);
					kmer[Kmer_arr_sz-1]|=b;

				}
				bool r_flip=0;
				uint64_t kmer2[100],f_kmer2[100];

				memcpy(kmer2,kmer,sizeof(uint64_t)*Kmer_arr_sz);
				memcpy(f_kmer2,kmer2,sizeof(uint64_t)*Kmer_arr_sz);

				get_rev_comp_seq_arr(f_kmer2,K_size,Kmer_arr_sz);

				if(uint64_t_cmp(kmer2,f_kmer2,Kmer_arr_sz)>0)
				{
					memcpy(kmer2,f_kmer2,sizeof(uint64_t)*Kmer_arr_sz);
					r_flip=1;
				}
				else
				{r_flip=0;}

				uint64_t hv=MurmurHash64A(kmer2,sizeof(uint64_t)*Kmer_arr_sz,0);
				uint64_t hash_idx=(size_t) (hv%ht->ht_sz);

				struct bucket0** ptr;
				ptr= &(ht->store_pos[hash_idx]);
				r_found=look_up_in_a_list0(kmer2,&ptr,Kmer_arr_sz);

				if(r_found)
				{
					// not in stack, put in stack 
					if(stacked_nodes[*ptr]==0)
					{
						Visited_Path[*ptr].cov=(int)(Visited_Path[new_node].cov+edge_ptr->edge_cov);
						Visited_Path[*ptr].depth=(Visited_Path[new_node].depth+1);
						int cum_len=(Visited_Path[new_node].len+edge_ptr->len+1);
						Visited_Path[*ptr].len=cum_len;//(Visited_Path[new_node].len+edge_ptr->len+1);

						Visited_Path[*ptr].last_bkt=new_node;
						Visited_Path[*ptr].last_bkt_edge=edge_ptr;
						if(r_flip)
						{
							stacked_nodes[*ptr]=-1;
							stacked_bkt.bktptr=*ptr;
							stacked_bkt.RightSearch=0;
							//kmer_stack.push_back(stacked_bkt);
						}
						else
						{
							stacked_nodes[*ptr]=1;
							stacked_bkt.bktptr=*ptr;
							stacked_bkt.RightSearch=1;
							//kmer_stack.push_back(stacked_bkt);
						}
						dist_ctgs[cum_len].push_back(stacked_bkt);
						NBs++;
					}
					else
					{
						if((stacked_nodes[*ptr]>0&&r_flip==0)||(stacked_nodes[*ptr]<0&&r_flip==1))
						{
							//backtrack if the same direction is found
							if(((Visited_Path[new_node].cov+edge_ptr->edge_cov)<=Visited_Path[*ptr].cov)||( BackCheckLoop0(*ptr,new_node,Visited_Path)==1))//loop
							{
								//if(0)
								if((((Visited_Path[new_node].cov+edge_ptr->edge_cov)/(Visited_Path[new_node].depth))>PathCovTh))
								{
									edge_ptr=edge_ptr->nxt_edge;
									continue;
								}
								if(abs(int((Visited_Path[new_node].len+edge_ptr->len+1)-Visited_Path[*ptr].len))>PathSim)
								{
									edge_ptr=edge_ptr->nxt_edge;
									continue;
								}
								if(!isSimplePath0(ht,new_node,Visited_Path,stacked_nodes))
								{
									edge_ptr=edge_ptr->nxt_edge;
									continue;

								}
								//backtrack the current path
								if(stacked_nodes[new_node]>2)
								{
									edge_ptr=edge_ptr->nxt_edge;
									BreakLinks0(stacked_nodes,new_node,(*ptr),K_size,edge_len);
									continue;
								}
								else
								{
									if(stacked_nodes[new_node]<-2)
									{
										edge_ptr=edge_ptr->nxt_edge;
										BreakLinks0(stacked_nodes,new_node,(*ptr),K_size,edge_len);
							
										continue;
									}
									else
									{
										struct bucket0 *freebkt=(new_node);

										if(freebkt==beg_bkt)
										{

											struct edge_node **edge_p2p=&(new_node->kmer_info.right);
											while(*edge_p2p!=edge_ptr)
											{
												edge_p2p=&((*edge_p2p)->nxt_edge);
											}
											if(*edge_p2p==edge_ptr)
											{
  												struct edge_node* f_edge_ptr=edge_ptr;
												edge_ptr=edge_ptr->nxt_edge;
												*edge_p2p=(*edge_p2p)->nxt_edge;
												free(f_edge_ptr);
											}
														
											continue;
										}
										edge_ptr=edge_ptr->nxt_edge;
										//free the node and edge.
										
										BreakLinks0(stacked_nodes,new_node,(*ptr),K_size,edge_len);
									
										BacktrackBubbleRemoval0(ht,merge_ht,new_node,beg_bkt,beg_bkt,Visited_Path, stacked_nodes,K_size);
								
										//Free_A_Node( ht, freebkt);
										freebkt->kmer_info.removed=1;
										//edge_ptr=NULL;
										continue;
				//
									}
								}
							}
							else
							{
								//backtrack the original path
								
								if(Visited_Path[*ptr].depth>1&&(((Visited_Path[*ptr].cov)/(Visited_Path[*ptr].depth-1))>PathCovTh))
								{
									edge_ptr=edge_ptr->nxt_edge;
									continue;
								}
								if(!isSimplePath0(ht,*ptr,Visited_Path,stacked_nodes))
								{
									edge_ptr=edge_ptr->nxt_edge;
									continue;

								}
													
								BacktrackBubbleRemoval0(ht,merge_ht,*ptr,beg_bkt,beg_bkt,Visited_Path , stacked_nodes,K_size);
						
								Visited_Path[*ptr].cov=(int)(Visited_Path[new_node].cov+edge_ptr->edge_cov);
								Visited_Path[*ptr].depth=Visited_Path[new_node].depth+1;
								Visited_Path[*ptr].len=(Visited_Path[new_node].len+edge_ptr->len+1);
								Visited_Path[*ptr].last_bkt=new_node;
								Visited_Path[*ptr].last_bkt_edge=edge_ptr;
								edge_ptr=edge_ptr->nxt_edge;
								continue;
							}

						}
						else
						{

							//don't do anything,since both strands are visited.

						}


					}

				}





				if(r_found==0||(*ptr)->kmer_info.removed==1)
				{
					//cout<<"WarningR"<<endl;


					struct edge_node **edge_p2p=&(new_node->kmer_info.right);
					while(*edge_p2p!=edge_ptr)
					{
						edge_p2p=&((*edge_p2p)->nxt_edge);
					}
					if(*edge_p2p==edge_ptr)
					{
  						struct edge_node* f_edge_ptr=edge_ptr;
						edge_ptr=edge_ptr->nxt_edge;
						*edge_p2p=(*edge_p2p)->nxt_edge;
						free(f_edge_ptr);
					}

					if(stacked_nodes[new_node]>2)
					{
						stacked_nodes[new_node]--;
						continue;
					}
					if(stacked_nodes[new_node]<-2)
					{
						stacked_nodes[new_node]++;
						continue;
					}

					//else
					struct bucket0 * freebkt=(new_node);

					if(freebkt==beg_bkt)
					{
			
						continue;
					}
					if(!isSimplePath0(ht,new_node,Visited_Path,stacked_nodes))
					{					
						continue;
					}
							
					BacktrackBubbleRemoval0(ht,merge_ht,new_node,beg_bkt,beg_bkt,Visited_Path, stacked_nodes,K_size);

					freebkt->kmer_info.removed=1;
					edge_ptr=NULL;
					continue;

				}

				edge_ptr=edge_ptr->nxt_edge;
			}


		}
		else
		{

			bool l_found=0;
			struct edge_node *edge_ptr;
//			struct edge_node **edge_p2p;
			edge_ptr=(new_node)->kmer_info.left;


			int lb=0;
			while(edge_ptr!=NULL)
			{
				lb++;
				edge_ptr=edge_ptr->nxt_edge;
			}
			edge_ptr=(new_node)->kmer_info.left;

			if(stacked_nodes[new_node]==-1&&lb>0)
			{
				stacked_nodes[new_node]=-1-lb;
			}
			if(lb==0)
			{
				stacked_nodes[new_node]=-2;

				struct bucket0 * freebkt=(new_node);

				if(freebkt==beg_bkt)
				{
					
							
					continue;
				}
	
				if(!isSimplePath0(ht,new_node,Visited_Path,stacked_nodes))
				{					
					continue;
				}
				BacktrackBubbleRemoval0(ht,merge_ht,new_node,beg_bkt,beg_bkt,Visited_Path , stacked_nodes,K_size);

				stacked_nodes[new_node]=-1;
				//Free_A_Node( ht, freebkt);
				freebkt->kmer_info.removed=1;
				edge_ptr=NULL;
				continue;
			}

			while(edge_ptr!=NULL)
			{


				memcpy(kmer,(new_node)->kmer_t,sizeof(uint64_t)*Kmer_arr_sz);
				memcpy(f_kmer,kmer,sizeof(uint64_t)*Kmer_arr_sz);
				get_rev_comp_seq_arr(f_kmer,K_size,Kmer_arr_sz);

				int edge_len=(edge_ptr->len+1);
				for(int g=0;g<edge_len;++g)
				{
					uint64_t b=(edge_ptr->edge>>2*g)&0x3;

					//kmer.kmer[0]&=(~((uint64_t)0x3<<(K_size-1-32)));

					R_shift_NB(kmer,2,Kmer_arr_sz);
					//R_shift_NB(kmer.kmer,2,2);
					kmer[0]|=(b<<(2*((K_size-1)%32)));
				}

				bool l_flip=0;

				uint64_t kmer2[100],f_kmer2[100];
				memcpy(kmer2,kmer,sizeof(uint64_t)*Kmer_arr_sz);
				memcpy(f_kmer2,kmer2,sizeof(uint64_t)*Kmer_arr_sz);

				get_rev_comp_seq_arr(f_kmer2,K_size,Kmer_arr_sz);

				if(uint64_t_cmp(kmer2,f_kmer2,Kmer_arr_sz)>0)
				{
					memcpy(kmer2,f_kmer2,sizeof(uint64_t)*Kmer_arr_sz);
					l_flip=1;
				}
				else
				{l_flip=0;}

				uint64_t hv=MurmurHash64A(&kmer2,sizeof(uint64_t)*Kmer_arr_sz,0);

				uint64_t hash_idx=(size_t) (hv%ht->ht_sz);

				struct bucket0** ptr;

				ptr= &(ht->store_pos[hash_idx]);

				l_found=look_up_in_a_list0(kmer2,&ptr,Kmer_arr_sz);
				//626211953659051268
				if(l_found)
				{

					if(stacked_nodes[*ptr]==0)
					{
						Visited_Path[*ptr].cov=(int)(Visited_Path[new_node].cov+edge_ptr->edge_cov);
						Visited_Path[*ptr].depth=(Visited_Path[new_node].depth+1);
						int cum_len=(int)(Visited_Path[new_node].len+edge_ptr->len+1);
						Visited_Path[*ptr].len=cum_len;//(Visited_Path[new_node].len+edge_ptr->len+1);

						Visited_Path[*ptr].last_bkt=new_node;
						Visited_Path[*ptr].last_bkt_edge=edge_ptr;
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
						if((stacked_nodes[*ptr]>0&&l_flip==1)||(stacked_nodes[*ptr]<0&&l_flip==0))
						{
							

							//backtrack if the same direction is found
							if((Visited_Path[new_node].cov+edge_ptr->edge_cov<=Visited_Path[*ptr].cov)||(BackCheckLoop0(*ptr,new_node,Visited_Path)==1))
							{
								//if(0)
								if((((Visited_Path[new_node].cov+edge_ptr->edge_cov)/(Visited_Path[new_node].depth))>PathCovTh))
								{
									edge_ptr=edge_ptr->nxt_edge;
									continue;
								}
								if(abs(int((Visited_Path[new_node].len+edge_ptr->len+1)-Visited_Path[*ptr].len))>PathSim)
								{
									edge_ptr=edge_ptr->nxt_edge;
									continue;
								}
								if(!isSimplePath0(ht,new_node,Visited_Path,stacked_nodes))
								{
									edge_ptr=edge_ptr->nxt_edge;
									continue;

								}
								//backtrack the current path
								if(stacked_nodes[new_node]>2)
								{
									edge_ptr=edge_ptr->nxt_edge;
			
									BreakLinks0(stacked_nodes,new_node,(*ptr),K_size,edge_len);
					
									continue;
								}
								else
								{
									if(stacked_nodes[new_node]<-2)
									{
										edge_ptr=edge_ptr->nxt_edge;
						
										BreakLinks0(stacked_nodes,new_node,(*ptr),K_size,edge_len);
						
										continue;
									}
									else
									{
										struct bucket0 * freebkt=(new_node);
										if(freebkt==beg_bkt)
										{
											struct edge_node **edge_p2p=&(new_node->kmer_info.left);
											while(*edge_p2p!=edge_ptr)
											{
												edge_p2p=&((*edge_p2p)->nxt_edge);
											}
											if(*edge_p2p==edge_ptr)
											{
  												struct edge_node* f_edge_ptr=edge_ptr;
												edge_ptr=edge_ptr->nxt_edge;
												*edge_p2p=(*edge_p2p)->nxt_edge;
												free(f_edge_ptr);
											}
											
											continue;
										}

										edge_ptr=edge_ptr->nxt_edge;
						
										BreakLinks0(stacked_nodes,new_node,(*ptr),K_size,edge_len);
			
										BacktrackBubbleRemoval0(ht,merge_ht,new_node,beg_bkt,beg_bkt,Visited_Path , stacked_nodes,K_size);
							
										//Free_A_Node( ht, freebkt);
										freebkt->kmer_info.removed=1;
										//edge_ptr=NULL;
										continue;


									}
								}
							}
							else
							{
								//backtrack the original path
								//if(0)
								if(Visited_Path[*ptr].depth>1&&(((Visited_Path[*ptr].cov)/(Visited_Path[*ptr].depth-1))>PathCovTh))
								{
									edge_ptr=edge_ptr->nxt_edge;
									continue;
								}
								if(!isSimplePath0(ht,*ptr,Visited_Path,stacked_nodes))
								{
									edge_ptr=edge_ptr->nxt_edge;
									continue;

								}
								
								BacktrackBubbleRemoval0(ht,merge_ht, *ptr,beg_bkt,beg_bkt,Visited_Path , stacked_nodes,K_size);
					
								Visited_Path[*ptr].cov=(int)(Visited_Path[new_node].cov+edge_ptr->edge_cov);
								Visited_Path[*ptr].depth=Visited_Path[new_node].depth+1;
								Visited_Path[*ptr].len=(int)(Visited_Path[new_node].len+edge_ptr->len+1);
								Visited_Path[*ptr].last_bkt=new_node;
								Visited_Path[*ptr].last_bkt_edge=edge_ptr;
								edge_ptr=edge_ptr->nxt_edge;
								continue;

							}

						}
						else
						{

							//don't do anything,since both strands are visited.

						}



					}

				}





				if(l_found==0||(*ptr)->kmer_info.removed==1)
				{
					//cout<<"WarningL"<<endl;
					//char c_str[200];
					//bitsarr2str(kmer2.kmer,K_size,c_str,2);
					//cout<<c_str<<endl;
					//bitsarr2str(bktptr->kmer_t2.kmer,K_size,c_str,2);
					//cout<<c_str<<endl;



					struct edge_node **edge_p2p=&(new_node->kmer_info.left);
					while(*edge_p2p!=edge_ptr)
					{
						//uint64_t edge=edge_ptr->edge;
						//bitsarr2str(&edge,16,c_str,1);
						//cout<<c_str<<endl;

						edge_p2p=&((*edge_p2p)->nxt_edge);
					}
					if(*edge_p2p==edge_ptr)
					{
  						struct edge_node* f_edge_ptr=edge_ptr;
						edge_ptr=edge_ptr->nxt_edge;
						*edge_p2p=(*edge_p2p)->nxt_edge;
						free(f_edge_ptr);
					}

					if(stacked_nodes[new_node]>2)
					{
						stacked_nodes[new_node]--;
						continue;
					}
					if(stacked_nodes[new_node]<-2)
					{
						stacked_nodes[new_node]++;
						continue;
					}
					//else
					struct bucket0 * freebkt=(new_node);
					if(freebkt==beg_bkt)
					{
						
							
						continue;
					}
					if(!isSimplePath0(ht,new_node,Visited_Path,stacked_nodes))
					{					
						continue;
					}

					BacktrackBubbleRemoval0(ht,merge_ht,new_node,beg_bkt,beg_bkt,Visited_Path , stacked_nodes,K_size);
			

					//Free_A_Node(ht,freebkt);
					freebkt->kmer_info.removed=1;
					edge_ptr=NULL;
					continue;

				}

				edge_ptr=edge_ptr->nxt_edge;
			}




		}

	}

}




// call the bubble removal

void GraphSimplification(hashtable *ht,hashtable * merge_ht, hashtable2 *ht2, hashtable2 * merge_ht2,int K_size,int gap,int PathCovTh,int max_depth,int PathSim)
{
	if(K_size<=32)
	{
		for(size_t i=0;i<ht->ht_sz;++i)
		{
			

		
			if((i+1)%10000000==0)//||i>6000000)
			{cout<<i+1<<endl;}

			struct bucket *bktptr=ht->store_pos[i];

			while(bktptr!=NULL)
			{
		
				int lb=0,rb=0;
				struct edge_node* left_edges=bktptr->kmer_info.left;
				//345438006799503592
				while(left_edges!=NULL)
				{
					lb++;
					left_edges=left_edges->nxt_edge;
				}

				struct edge_node* right_edges=bktptr->kmer_info.right;
				while(right_edges!=NULL)
				{
					rb++;
					right_edges=right_edges->nxt_edge;
				}
				if(lb>1)
				{
					struct stacked_bucket stacked_bkt;
					stacked_bkt.bktptr=bktptr;
					stacked_bkt.RightSearch=0;
					list<struct stacked_bucket> kmer_stack;
					kmer_stack.push_back(stacked_bkt);

					BFSearchBubbleRemoval(ht,merge_ht,bktptr,K_size,gap,kmer_stack,PathCovTh,max_depth,PathSim);

				}

				if(rb>1)
				{
					struct stacked_bucket stacked_bkt;
					stacked_bkt.bktptr=bktptr;
					//if(bktptr->kmer_t.kmer==1120767112706309364)
					//{cout<<"";}
					stacked_bkt.RightSearch=1;
					list<struct stacked_bucket> kmer_stack;
					kmer_stack.push_back(stacked_bkt);
					BFSearchBubbleRemoval(ht,merge_ht,bktptr,K_size,gap,kmer_stack,PathCovTh,max_depth,PathSim);
				}

				bktptr=bktptr->nxt_bucket;
			}
		}
	}
	else
	{
		if(K_size<=64)
		{
			for(size_t i=0;i<ht2->ht_sz;++i)
			{

				//if(i==13174)//11390)//2483)
				//{cout<<"";}
					//if(1)
				if(i%10000000==0)//||i>6000000)
				{cout<<i<<endl;}

				struct bucket2 *bktptr=ht2->store_pos[i];

				while(bktptr!=NULL)
				{
					int lb=0,rb=0;
					struct edge_node* left_edges=bktptr->kmer_info.left;
					//345438006799503592
					while(left_edges!=NULL)
					{
						lb++;
						left_edges=left_edges->nxt_edge;
					}

					struct edge_node* right_edges=bktptr->kmer_info.right;
					while(right_edges!=NULL)
					{
						rb++;
						right_edges=right_edges->nxt_edge;
					}
					if(lb>1)
					{
						struct stacked_bucket2 stacked_bkt;
						stacked_bkt.bktptr=bktptr;
						stacked_bkt.RightSearch=0;
						list<struct stacked_bucket2> kmer_stack;
						kmer_stack.push_back(stacked_bkt);
						
						BFSearchBubbleRemoval2(ht2,merge_ht2,bktptr,K_size,gap,kmer_stack,PathCovTh,max_depth,PathSim);

					}

					if(rb>1)
					{
						struct stacked_bucket2 stacked_bkt;
						stacked_bkt.bktptr=bktptr;
						stacked_bkt.RightSearch=1;
						list<struct stacked_bucket2> kmer_stack;
						kmer_stack.push_back(stacked_bkt);
						BFSearchBubbleRemoval2(ht2,merge_ht2,bktptr,K_size,gap,kmer_stack,PathCovTh,max_depth,PathSim);
					}

					bktptr=bktptr->nxt_bucket;
				}
			}


		}
	}

}


void GraphSimplification3(hashtable3 *ht3, hashtable3 * merge_ht3,int K_size,int gap,int PathCovTh,int max_depth,int PathSim)
{
		
	for(size_t i=0;i<ht3->ht_sz;++i)
	{
		if(i%10000000==0)//||i>6000000)
		{cout<<i<<endl;}

		struct bucket3 *bktptr=ht3->store_pos[i];

		while(bktptr!=NULL)
		{
			int lb=0,rb=0;
			struct edge_node* left_edges=bktptr->kmer_info.left;
			//345438006799503592
			while(left_edges!=NULL)
			{
				lb++;
				left_edges=left_edges->nxt_edge;
			}

			struct edge_node* right_edges=bktptr->kmer_info.right;
			while(right_edges!=NULL)
			{
				rb++;
				right_edges=right_edges->nxt_edge;
			}
			if(lb>1)
			{
				struct stacked_bucket3 stacked_bkt;
				stacked_bkt.bktptr=bktptr;
				stacked_bkt.RightSearch=0;
				list<struct stacked_bucket3> kmer_stack;
				kmer_stack.push_back(stacked_bkt);
						
				BFSearchBubbleRemoval3(ht3,merge_ht3,bktptr,K_size,gap,kmer_stack,PathCovTh,max_depth,PathSim);

			}

			if(rb>1)
			{
				struct stacked_bucket3 stacked_bkt;
				stacked_bkt.bktptr=bktptr;
				stacked_bkt.RightSearch=1;
				list<struct stacked_bucket3> kmer_stack;
				kmer_stack.push_back(stacked_bkt);
				BFSearchBubbleRemoval3(ht3,merge_ht3,bktptr,K_size,gap,kmer_stack,PathCovTh,max_depth,PathSim);
			}

			bktptr=bktptr->nxt_bucket;
		}
	}


}

void GraphSimplification4(hashtable4 *ht4, hashtable4 * merge_ht4,int K_size,int gap,int PathCovTh,int max_depth,int PathSim)
{
	
	
	for(size_t i=0;i<ht4->ht_sz;++i)
	{
		if(i%10000000==0)//||i>6000000)
		{cout<<i<<endl;}

		struct bucket4 *bktptr=ht4->store_pos[i];

		while(bktptr!=NULL)
		{
			int lb=0,rb=0;
			struct edge_node* left_edges=bktptr->kmer_info.left;
			//345438006799503592
			while(left_edges!=NULL)
			{
				lb++;
				left_edges=left_edges->nxt_edge;
			}

			struct edge_node* right_edges=bktptr->kmer_info.right;
			while(right_edges!=NULL)
			{
				rb++;
				right_edges=right_edges->nxt_edge;
			}
			if(lb>1)
			{
				struct stacked_bucket4 stacked_bkt;
				stacked_bkt.bktptr=bktptr;
				stacked_bkt.RightSearch=0;
				list<struct stacked_bucket4> kmer_stack;
				kmer_stack.push_back(stacked_bkt);
						
				BFSearchBubbleRemoval4(ht4,merge_ht4,bktptr,K_size,gap,kmer_stack,PathCovTh,max_depth,PathSim);

			}

			if(rb>1)
			{
				struct stacked_bucket4 stacked_bkt;
				stacked_bkt.bktptr=bktptr;
				stacked_bkt.RightSearch=1;
				list<struct stacked_bucket4> kmer_stack;
				kmer_stack.push_back(stacked_bkt);
				BFSearchBubbleRemoval4(ht4,merge_ht4,bktptr,K_size,gap,kmer_stack,PathCovTh,max_depth,PathSim);
			}

			bktptr=bktptr->nxt_bucket;
		}
	}


}


void GraphSimplification0(hashtable0 *ht, hashtable0 * merge_ht,int K_size,int gap,int PathCovTh,int max_depth,int PathSim)
{
	
	
	for(size_t i=0;i<ht->ht_sz;++i)
	{
		if(i%10000000==0&&i>0)//||i>6000000)
		{cout<<i<<endl;}

		struct bucket0 *bktptr=ht->store_pos[i];

		while(bktptr!=NULL)
		{
			int lb=0,rb=0;
			struct edge_node* left_edges=bktptr->kmer_info.left;
			//345438006799503592
			while(left_edges!=NULL)
			{
				lb++;
				left_edges=left_edges->nxt_edge;
			}

			struct edge_node* right_edges=bktptr->kmer_info.right;
			while(right_edges!=NULL)
			{
				rb++;
				right_edges=right_edges->nxt_edge;
			}
			if(lb>1)
			{
				struct stacked_bucket0 stacked_bkt;
				stacked_bkt.bktptr=bktptr;
				stacked_bkt.RightSearch=0;
				list<struct stacked_bucket0> kmer_stack;
				kmer_stack.push_back(stacked_bkt);
						
				BFSearchBubbleRemoval0(ht,merge_ht,bktptr,K_size,gap,kmer_stack,PathCovTh,max_depth,PathSim);

			}

			if(rb>1)
			{
				struct stacked_bucket0 stacked_bkt;
				stacked_bkt.bktptr=bktptr;
				stacked_bkt.RightSearch=1;
				list<struct stacked_bucket0> kmer_stack;
				kmer_stack.push_back(stacked_bkt);
				BFSearchBubbleRemoval0(ht,merge_ht,bktptr,K_size,gap,kmer_stack,PathCovTh,max_depth,PathSim);
			}

			bktptr=bktptr->nxt_bucket;
		}
	}


}




#endif
