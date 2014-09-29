/*
Original Code From:
Copyright (C) 2006 Pedro Felzenszwalb
Modifications (may have been made) Copyright (C) 2011, 2012
  Chenliang Xu, Jason Corso.
Modifications Copyright (C) 2014
  Haw-Shiuan Chang
  
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
*/

/* Implements the node data structure. */

#ifndef DISJOINT_SET_H
#define DISJOINT_SET_H

#include <vector>
#include <algorithm>

#include "image.h"
#include "misc.h"
#include "histogram.h"
#include "edges.h"
#include "pnmfile.h"
#include "parameters.h"

using namespace std;

#define max(a,b) ((a)>(b)?(a):(b))
#define min(a,b) ((a)<(b)?(a):(b))

bool operator<(const edge &a, const edge &b) {
	return a.w < b.w;
}



/* Each pixel has such a node structure. 
 * Once go up a level, changes cannot be made.
 */

typedef struct {

	unsigned int p;

	// use to update MST mergeing condition
	float mst;
	int size;

	// for representative nodes
	Histogram<float> *His_L;
	Histogram<float> *His_a;
	Histogram<float> *His_b;
	
	// position in the 3D volume
} region_node;

typedef struct {

	// contain its representative node in each level
	unsigned int p;

	// use to update MST mergeing condition
	float mst;
	int size;

	// record the pixel color
	rgb color;
	float Lab_L;
	float Lab_a;
	float Lab_b;

	// position in the 3D volume
	int x;
	int y;
	int z;
} voxel_node;

/* Node Forest for the 3D video volume */
class universe {
public:
	// all nodes in the forest
	voxel_node *p_node;
	region_node **r_node;
	int *node_num; // total nodes in the forest

	// initialize 1st level forest
	universe(int num_frame, int width, int height, image<float> *r[],
			image<float> *g[], image<float> *b[], int hie_num);

	//void segment_pixel_graph(vector<edge>* edges_remain, edge *edges, int num_edges, float c,int min_size);
	void segment_pixel_graph(edge_index *edges_remain_index, edge *edges, int num_edges, float c,int min_size);
	
	void segment_graph_region(edge_index *edges_remain_index, edge *edges, float c_reg, int level,int min_size);
	
	// update an initial status of higher level
	void update(int level);

	//void fill_edge_weight(vector<edge>& edges_region_out,vector<edge>& edges_region_in, int level);
	unsigned int fill_edge_weight(const edge_index edges_remain_index,edge *edges_down,edge **edges_up, int level);
	
	void generate_output(char *path, int num_frame, int width, int height) ;

	// destroy
	~universe();

	//float get_mst(int a) { return p_node[a].mst; }
	//int get_size(int a) { return p_node[a].size; }
	//Histogram<float> *get_His_L(int a) { return p_node[a].His_L; }
	//Histogram<float> *get_His_a(int a) { return p_node[a].His_a; }
	//Histogram<float> *get_His_b(int a) { return p_node[a].His_b; }

private:

	int level_total;
	
	//void RGBtoLab(int a);
	// find parent in one level
	int find_in_level(int a, int level);
	int find_in_level_done(int a, int level);
	// merging two segments (join two representative nodes)
	int join(int a, int b, float mst, int level);
	// get information
	int get_size(int a, int level){ 
		if(level==0)
			return p_node[a].size; 
		else
			return r_node[level-1][a].size; 
	}
};

universe::universe(int num_frame, int width, int height, image<float> *r[],
		image<float> *g[], image<float> *b[], int hie_num) {
	
	level_total = hie_num + 1;
	int p_num = num_frame * width * height;
	node_num=new int[level_total];
	node_num[0]=p_num;
	p_node=new voxel_node[p_num];
	r_node=new region_node*[hie_num];

	// initialize each pixel nodes
	for (int i = 0; i < p_num; i++) {
		//p_node[i].p=new unsigned int[level_total];
		//p_node[i].p[0]=i;
		p_node[i].p=i;
		p_node[i].mst=0;
		p_node[i].size=1;
		p_node[i].z= i / (width * height);
		p_node[i].y = (i % (width * height)) / width;
		p_node[i].x = (i % (width * height)) % width;

		p_node[i].color.r =
				imRef(r[p_node[i].z], p_node[i].x, p_node[i].y);
		p_node[i].color.g =
				imRef(g[p_node[i].z], p_node[i].x, p_node[i].y);
		p_node[i].color.b =
				imRef(b[p_node[i].z], p_node[i].x, p_node[i].y);

		Lab color_out=RGBtoLab(p_node[i].color);
		p_node[i].Lab_L=color_out.L;
		p_node[i].Lab_a=color_out.a;
		p_node[i].Lab_b=color_out.b;

		/*p_node.push_back(node());
		p_node.back().p.resize(level_total);
		p_node.back().p.front() = i;
		p_node.back().mst = 0;
		p_node.back().size = 1;*/


		// local position
		/*p_node.back().z = i / (width * height);
		p_node.back().y = (i % (width * height)) / width;
		p_node.back().x = (i % (width * height)) % width;

		p_node.back().color.r =
				imRef(r[p_node.back().z], p_node.back().x, p_node.back().y);
		p_node.back().color.g =
				imRef(g[p_node.back().z], p_node.back().x, p_node.back().y);
		p_node.back().color.b =
				imRef(b[p_node.back().z], p_node.back().x, p_node.back().y);
		RGBtoXYZ(i);
		XYZtoLab(i);*/

		//p_node.back().His_L = new Histogram<float>(20, 0, 100);
		//p_node.back().His_a = new Histogram<float>(20, -50, 50);
		//p_node.back().His_b = new Histogram<float>(20, -50, 50);
	}
}

universe::~universe() {
	delete [] p_node;
	for(int i = 0; i < level_total-1 ; i++){
		//for (int j = 0; j < node_num[i+1]; j++) {
		//	delete r_node[i][j].His_L;
		//	delete r_node[i][j].His_a;
		//	delete r_node[i][j].His_b;
		//}
		delete [] r_node[i];
	}
	delete [] r_node;
	delete [] node_num;
}

/* pixel level minimum spanning tree merge */
//void universe::segment_pixel_graph(vector<edge>* edges_remain, edge *edges, int num_edges, float c,int min_size) {
void universe::segment_pixel_graph(edge_index *edges_remain_index, edge *edges, int num_edges, float c,int min_size) {
	// new vector containing remain edges
	//edges_remain->clear();
	unsigned int counter=0;
	// sort edges by weight
	sort(edges, edges + num_edges);
	// for each edge, in non-decreasing weight order...
	for (int i = 0; i < num_edges; i++) {
		//edge *pedge = &edges[i];
		// components conected by this edge
		//int a = find_in_level(pedge->a, 0);
		//int b = find_in_level(pedge->b, 0);
		int a = find_in_level(edges[i].a, 0);
		int b = find_in_level(edges[i].b, 0);
		if (a != b) {
			if((edges[i].w <= p_node[a].mst + (c / p_node[a].size))
					&& (edges[i].w <= p_node[b].mst + (c / p_node[b].size)) ) 
			// merging objective function
				join(a, b, edges[i].w, 0);	
			else
				++counter;
		}
		//if (a != b) {
		//	// merging objective function
		//	//if ((pedge->w <= p_node[a].mst + (c / p_node[a].size))
		//	//		&& (pedge->w <= p_node[b].mst + (c / p_node[b].size))) {
		//	//	//if ( join(a, b, pedge->w, 0) == 1)
		//	//	//	edges_remain->push_back(*pedge);
		//	//	join(a, b, pedge->w, 0);	
		//	//} else {
		//	//	edges_remain->push_back(*pedge);
		//	//}
		//	if ((edges[i].w <= p_node[a].mst + (c / p_node[a].size))
		//			&& (edges[i].w <= p_node[b].mst + (c / p_node[b].size))) {
		//		join(a, b, edges[i].w, 0);	
		//	} else {
		//		//edges_remain->push_back(*pedge);
		//		edges_remain_index->index[counter]=i;
		//		counter++;
		//	}
		//}

	}
	edges_remain_index->index=new unsigned int[counter];
	counter=0;
	// optional merging small components
	for (int i = 0; i < num_edges; i++) {
		int a = find_in_level(edges[i].a, 0);
		int b = find_in_level(edges[i].b, 0);
		//if ((a != b) && ((p_node[a].size < min_size) || (p_node[b].size < min_size)))
		//	join(a, b, 0, 0);
		if (a != b) {
			if ((p_node[a].size < min_size) || (p_node[b].size < min_size)){
				join(a, b, 0, 0);
			} else {
				//edges_remain->push_back(*pedge);
				edges_remain_index->index[counter]=i;
				counter++;
			}	
		}
	}
	edges_remain_index->num=counter;

}

/* region graph level minimum spanning tree merge */
void universe::segment_graph_region(edge_index *edges_remain_index, edge *edges, float c_reg, int level, int min_size) {
	unsigned int num_edges=edges_remain_index->num;
	unsigned int counter=0;
	sort(edges, edges + num_edges);
	//sort(edges_region->begin(), edges_region->end());
	region_node* temp_r_node=r_node[level-1];
	//for (int i=0; i < (int) edges_region->size(); i++) {
	//char file_name[100];
	//sprintf(file_name,"debug_log_%d.txt",level);
	//FILE *log_fid=fopen(file_name,"w");

	for (int i=0; i < (int) num_edges; i++) {
		/*int a = find_in_level(edges_region->at(i).a, level);
		int b = find_in_level(edges_region->at(i).b, level);*/
		int a = find_in_level(edges[i].a, level);
		int b = find_in_level(edges[i].b, level);
		if (a != b) {
			//fprintf(log_fid,"%d, %d, %f\n",a,b,edges[i].w);
			if((edges[i].w <= temp_r_node[a].mst + (c_reg / temp_r_node[a].size )) && 
				(edges[i].w <= temp_r_node[b].mst + (c_reg / temp_r_node[b].size ))) 
				/*if ( join(a, b, edges_region->at(i).w, level) == 1)
					edges_remain->push_back(edges_region->at(i));*/
				 join(a, b, edges[i].w, level);
			else
				++counter;
		}
	}
	//fclose(log_fid);

	edges_remain_index->index=new unsigned int[counter];
	counter=0;
	for (int it = 0; it < (int) num_edges; it++) {
		int a = find_in_level(edges[it].a, level);
		int b = find_in_level(edges[it].b, level);
		/*if ((a != b) && ((temp_r_node[a].size < min_size) || (temp_r_node[b].size < min_size)))
				join(a, b, 0, level);*/
		if (a != b) {
			if ((temp_r_node[a].size < min_size) || (temp_r_node[b].size < min_size))
				join(a, b, 0, level);
			else {
				edges_remain_index->index[counter]=it;
				counter++;
				//edges_remain->push_back(edges_region->at(i));
			}
		}
	}
	edges_remain_index->num=counter;
}

void universe::update(int level) {

	unsigned int over_segment_num=0;
	region_node* temp_r_node;
	int n_num=node_num[level];
	int *index_array=new int[n_num];
	for (int i = 0; i < n_num; ++i) {
		find_in_level(i, level);
		if(get_size(i, level)!=0)
			index_array[i]=over_segment_num++;
		else
			index_array[i]=-1;
	}
	node_num[level+1]=over_segment_num;
	temp_r_node=new region_node[over_segment_num];
	for(int i = 0; i < over_segment_num; ++i) {
		temp_r_node[i].His_L = new Histogram<float>(PARA_HIST_BIN_NUM, 0, 100);
		temp_r_node[i].His_a = new Histogram<float>(PARA_HIST_BIN_NUM, -50, 50);
		temp_r_node[i].His_b = new Histogram<float>(PARA_HIST_BIN_NUM, -50, 50);
	}
	if (level == 0) {
		
		/*float max_L=-10000000000;
		float min_L=10000000000;
		float max_a=-10000000000;
		float min_a=10000000000;
		float max_b=-10000000000;
		float min_b=10000000000;*/
		
		//float max_L_c=0;
		//float min_L_c=0;
		//float max_a_c=0;
		//float min_a_c=0;
		//float max_b_c=0;
		//float min_b_c=0;

		for (int i = 0; i < n_num; ++i) {
			//int p=index_array[p_node[i].p[level]];
			//p_node[i].p[level]=p;
			int p=index_array[p_node[i].p];
			p_node[i].p=p;
			temp_r_node[p].His_L->addSample(p_node[i].Lab_L);
			temp_r_node[p].His_a->addSample(p_node[i].Lab_a);
			temp_r_node[p].His_b->addSample(p_node[i].Lab_b);
			//if(p_node[i].Lab_L>max_L)
			//	max_L=p_node[i].Lab_L;
			//if(p_node[i].Lab_a>max_a)
			//	max_a=p_node[i].Lab_a;
			//if(p_node[i].Lab_b>max_b)
			//	max_b=p_node[i].Lab_b;
			//if(p_node[i].Lab_L<min_L)
			//	min_L=p_node[i].Lab_L;
			//if(p_node[i].Lab_a<min_a)
			//	min_a=p_node[i].Lab_a;
			//if(p_node[i].Lab_b<min_b)
			//	min_b=p_node[i].Lab_b;

			//if(p_node[i].Lab_L>100)
			//	max_L_c++;
			//if(p_node[i].Lab_a>50)
			//	max_a_c++;
			//if(p_node[i].Lab_b>50)
			//	max_b_c++;
			//if(p_node[i].Lab_L<0)
			//	min_L_c++;
			//if(p_node[i].Lab_a<-50)
			//	min_a_c++;
			//if(p_node[i].Lab_b<-50)
			//	min_b_c++;

			if(p_node[i].size!=0)
				temp_r_node[p].size=p_node[i].size;
		}
		//printf("%f %f %f\n",max_L_c/n_num,max_a_c/n_num,max_b_c/n_num);
		//printf("%f %f %f\n",min_L_c/n_num,min_a_c/n_num,min_b_c/n_num);
	
	/*	for (int i = 0; i < num; i++) {
			int p = find_in_level(i, level);
			p_node[i].p[level + 1] = p_node[i].p[level];
			p_node[i].mst = 0;
			p_node[p].His_L->addSample(p_node[i].Lab_L);
			p_node[p].His_a->addSample(p_node[i].Lab_a);
			p_node[p].His_b->addSample(p_node[i].Lab_b);
		}*/
	} else {
		for (int i = 0; i < n_num; ++i) {
			int p=index_array[r_node[level-1][i].p];
			r_node[level-1][i].p=p;
			temp_r_node[p].His_L->mergeHistogram(*r_node[level-1][i].His_L);
			temp_r_node[p].His_a->mergeHistogram(*r_node[level-1][i].His_a);
			temp_r_node[p].His_b->mergeHistogram(*r_node[level-1][i].His_b);
			if(r_node[level-1][i].size!=0)
				temp_r_node[p].size=r_node[level-1][i].size;
		}
		
		//for (int i = 0; i < num; i++) {
		//	int p = find_in_level(i, level);
		//	p_node[i].p[level + 1] = p_node[i].p[level];
		//	p_node[i].mst = 0;
		//	if (i != p) {
		//		p_node[p].His_L->mergeHistogram(*p_node[i].His_L);
		//		p_node[p].His_a->mergeHistogram(*p_node[i].His_a);
		//		p_node[p].His_b->mergeHistogram(*p_node[i].His_b);
		//	}
		//}
		for (int j = 0; j < node_num[level]; j++) {
			delete r_node[level-1][j].His_L;
			delete r_node[level-1][j].His_a;
			delete r_node[level-1][j].His_b;
		}
	}
	for (int i = 0; i < over_segment_num; i++) {
		temp_r_node[i].p=i;
		temp_r_node[i].mst=0;
	}
	r_node[level]=temp_r_node;
}

/* fill region graph edges */

unsigned int universe::fill_edge_weight(const edge_index edges_remain_index,edge *edges_down,edge **edges_up_ptr, int level) {
	region_node* temp_r_node=r_node[level];
	unsigned int node_num_level=node_num[level+1];
	unsigned int *nei_index=new unsigned int[node_num_level*MAX_NEI_CACHE_SIZE];
	unsigned int *nei_num=new unsigned int[node_num_level];
	//float *nei_w=new float[node_num_level*MAX_NEI_CACHE_SIZE];
	for(int i=0;i<node_num_level;++i)
		nei_num[i]=0;
	//unsigned int count=0;
	//for (int i = 0; i < ((int) edges_region_in.size()); i++) {
	unsigned int counter=0;
	for (int i = 0; i < ((int) edges_remain_index.num); i++) {
		unsigned int index=edges_remain_index.index[i];
		int a = edges_down[index].a;
		int b = edges_down[index].b;
		int a_p = find_in_level_done(a, level);
		int b_p = find_in_level_done(b, level);
		if (a_p != b_p) {
			//searching neighbors
			unsigned int nei_loc;
			for(nei_loc=0;nei_loc<nei_num[a_p];++nei_loc)
				if(nei_index[a_p*MAX_NEI_CACHE_SIZE+nei_loc]==b_p)
					break;
			
			if(nei_loc==nei_num[a_p])
			{
				//if(nei_num[a_p]>=MAX_NEI_CACHE_SIZE)
				//	nei_num[a_p]=0;
				if(nei_num[a_p]<MAX_NEI_CACHE_SIZE)
				{
					nei_index[a_p*MAX_NEI_CACHE_SIZE+nei_num[a_p]]=b_p;
					//nei_w[a_p*MAX_NEI_CACHE_SIZE+nei_num[a_p]]=w;
					++nei_num[a_p];
				}
				//if(nei_num[b_p]>=MAX_NEI_CACHE_SIZE)
				//	nei_num[b_p]=0;
				if(nei_num[b_p]<MAX_NEI_CACHE_SIZE)
				{
					nei_index[b_p*MAX_NEI_CACHE_SIZE+nei_num[b_p]]=a_p;
					//nei_w[b_p*MAX_NEI_CACHE_SIZE+nei_num[b_p]]=w;
					++nei_num[b_p];
				}
				//edges_up[counter].w =temp_r_node[a_p].His_L->chiSquared( *(temp_r_node[b_p].His_L))
				edges_down[counter].w =temp_r_node[a_p].His_L->chiSquared( *(temp_r_node[b_p].His_L))
					+ temp_r_node[a_p].His_a->chiSquared( *(temp_r_node[b_p].His_a))
					+ temp_r_node[a_p].His_b->chiSquared( *(temp_r_node[b_p].His_b));
				edges_down[counter].a=a_p;
				edges_down[counter].b=b_p;
				//edges_up[counter].a=a_p;
				//edges_up[counter].b=b_p;
				++counter;
			}
			/*else
			{
				w=nei_w[a_p*MAX_NEI_CACHE_SIZE+nei_loc];
			}*/
			
		}
	}
	delete [] nei_index;
	delete [] nei_num;
	*edges_up_ptr=new edge[counter];
	for(int i=0;i<counter;++i)
	{
		(*edges_up_ptr)[i].w=edges_down[i].w;
		(*edges_up_ptr)[i].a=edges_down[i].a;
		(*edges_up_ptr)[i].b=edges_down[i].b;
	}
	return counter;
}


///* fill region graph edges */
//void universe::fill_edge_weight(vector<edge>& edges_region_out,vector<edge>& edges_region_in, int level) {
//	region_node* temp_r_node=r_node[level];
//	//unsigned int count=0;
//	for (int i = 0; i < ((int) edges_region_in.size()); i++) {
//		int a = edges_region_in[i].a;
//		int b = edges_region_in[i].b;
//		/*int a_p = find_in_level(a, level);
//		int b_p = find_in_level(b, level);*/
//		int a_p = find_in_level_done(a, level);
//		int b_p = find_in_level_done(b, level);
//		if (a_p != b_p) {
//			edges_region_in[i].w = temp_r_node[a_p].His_L->chiSquared( *(temp_r_node[b_p].His_L))
//					+ temp_r_node[a_p].His_a->chiSquared( *(temp_r_node[b_p].His_a))
//					+ temp_r_node[a_p].His_b->chiSquared( *(temp_r_node[b_p].His_b));
//			edges_region_in[i].a=a_p;
//			edges_region_in[i].b=b_p;
//			edges_region_out.push_back(edges_region_in[i]);
//			/*edges_region_out[count].w = temp_r_node[a_p].His_L->chiSquared( *(temp_r_node[b_p].His_L))
//					+ temp_r_node[a_p].His_a->chiSquared( *(temp_r_node[b_p].His_a))
//					+ temp_r_node[a_p].His_b->chiSquared( *(temp_r_node[b_p].His_b));
//			edges_region_out[count].a=a;
//			edges_region_out[count].b=b;*/
//			//count++;
//			/*edges_region[i].w = mess->get_His_L(a_p)->chiSquared(*mess->get_His_L(b_p))
//					+ mess->get_His_a(a_p)->chiSquared(*mess->get_His_a(b_p))
//					+ mess->get_His_b(a_p)->chiSquared(*mess->get_His_b(b_p));*/
//		//} else {
//		//	edges_region[i].w = 0;
//		//	//acceleration: erase such edge
//		}
//	}
//}

int universe::find_in_level_done(int a, int level) {
	if(level==0)
	{
		return p_node[a].p;
	}else
	{
		return r_node[level-1][a].p;
	}
}

int universe::find_in_level(int a, int level) {
	int element = a;
	//acceleration: write a different function
	if(level==0)
	{
		//while (element != p_node[element].p[level])
		//	element = p_node[element].p[level];
		//p_node[a].p[level] = element;
		while (element != p_node[element].p)
			element = p_node[element].p;
		p_node[a].p = element;
	}else
	{
		while (element != r_node[level-1][element].p)
			element = r_node[level-1][element].p;
		r_node[level-1][a].p=element;
	}
	return element;
}

int universe::join(int a, int b, float mst, int level) {

	//acceleration: write a different function
	if(level==0)
	{
		mst = max(max(p_node[a].mst, p_node[b].mst), mst);
		if (a < b) {
			//p_node[b].p[level] = a;
			p_node[b].p = a;
			p_node[a].size += p_node[b].size;
			p_node[b].size=0;
			p_node[a].mst = mst;
		} else {
			//p_node[a].p[level] = b;
			p_node[a].p = b;
			p_node[b].size += p_node[a].size;
			p_node[a].size = 0;
			p_node[b].mst = mst;
		}
	}
	else
	{
		region_node *temp_r_node=r_node[level-1];
		mst = max(max(temp_r_node[a].mst, temp_r_node[b].mst), mst);
		if (a < b) {
			temp_r_node[b].p = a;
			temp_r_node[a].size += temp_r_node[b].size;
			temp_r_node[b].size=0;
			temp_r_node[a].mst = mst;
		} else {
			temp_r_node[a].p = b;
			temp_r_node[b].size += temp_r_node[a].size;
			temp_r_node[a].size = 0;
			temp_r_node[b].mst = mst;
		}
	}

	return 0;

	/*if (p_node[a].rank > p_node[b].rank) {
		p_node[b].p[level] = a;
		p_node[a].size += p_node[b].size;
		p_node[a].mst = mst;
	} else {
		p_node[a].p[level] = b;
		p_node[b].size += p_node[a].size;
		p_node[b].mst = mst;
		if (p_node[a].rank == p_node[b].rank)
			p_node[b].rank++;
	}*/


}

/* Save Output */
void universe::generate_output(char *path, int num_frame, int width, int height) {

	int num_vertices=node_num[0];

	// Prepare vector of random/unique values for access by counter
	//vector<int> randomNumbers;
	const int Large_const=16777215;
	//randomNumbers.resize(Large_const);
	unsigned int *randomNumbers=new unsigned int[Large_const];
	for(int i=0;i<Large_const;++i)
		randomNumbers[i]=i;
	unsigned int range=node_num[1];
	my_RandomPermuteRange(randomNumbers,Large_const,range);
	//random_shuffle(randomNumbers.begin(),randomNumbers.end());
	//unsigned int r = rand() % 10000000;
	//juRandomPermuteRange(16777215, randomNumbers, &r);

	char savepath[1024];
	image<rgb>** output = new image<rgb>*[num_frame];
	rgb* colors = new rgb[num_vertices];

	unsigned int *p_level=new unsigned int[num_vertices];
	for(int i=0;i<num_vertices;++i)
	{
		p_level[i]=p_node[i].p;
	}

	/*unsigned int **p_all=new unsigned int*[level_total];
	for(int i=0;i<level_total;++i)
	{
		p_all[i]=new unsigned int[num_vertices];
	}
	for(int i=0;i<num_vertices;++i)
	{
		p_all[0][i]=p_node[i].p;
	}
	for(int level=1;level<level_total;++level)
	{
		for(int j=0;j<num_vertices;++j)
		{
			p_all[level][j]=r_node[level-1][p_all[level-1][j]].p;
		}
	}*/

	// write out the ppm files.
	for (int k = 0; k < level_total; k++) {
		//set all voxel colors to white
		for (int i = 0; i < num_vertices; i++){	
			colors[i].r = colors[i].g = colors[i].b = 255; 
		}

		int counter = 0;
		for (int i = 0; i < num_frame; i++) {
			sprintf(savepath,  "%s/%02d/%05d.ppm", path, k, i + 1);
			output[i] = new image<rgb>(width, height);
			for (int y = 0; y < height; y++) {
				for (int x = 0; x < width; x++) {
					/*int comp = find_in_level(
							y * width + x + i * (width * height), k);*/
					int comp =p_level[y * width + x + i * (width * height)];

					// choose color from random permuted numbers
					while(isWhite(colors[comp])){
						setColor(colors[comp], randomNumbers[counter]);
						counter++;
						if(counter==Large_const){printf("**WARNING!** All "
						"colors have been used. Colors will no longer "
						"be unique.**\n");}
					}	
					imRef(output[i], x, y) = colors[comp];
				}
			}
			savePPM(output[i], savepath);
		}
		for (int i = 0; i < num_frame; i++)
			delete output[i];
		if(k<level_total-1)
		{
			for(int j=0;j<num_vertices;++j)
			{
				p_level[j]=r_node[k][p_level[j]].p;
			}
		}
	}
	delete[] colors;
	delete[] output;
	delete[] p_level;

}


/* RGB to XYZ color space */
//void universe::RGBtoXYZ(int a) {
//	float var_R = p_node[a].color.r / 255.0;
//	float var_G = p_node[a].color.g / 255.0;
//	float var_B = p_node[a].color.b / 255.0;
//
//	if (var_R > 0.04045)
//		var_R = pow(((var_R + 0.055) / 1.055), 2.4);
//	else
//		var_R = var_R / 12.92;
//
//	if (var_G > 0.04045)
//		var_G = pow(((var_G + 0.055) / 1.055), 2.4);
//	else
//		var_G = var_G / 12.92;
//
//	if (var_B > 0.04045)
//		var_B = pow(((var_B + 0.055) / 1.055), 2.4);
//	else
//		var_B = var_B / 12.92;
//
//	var_R = var_R * 100;
//	var_G = var_G * 100;
//	var_B = var_B * 100;
//
//	//Observer. = 2¢X, Illuminant = D65
//	p_node[a].XYZ_X = var_R * 0.412453f + var_G * 0.357580f + var_B * 0.180423f;
//	p_node[a].XYZ_Y = var_R * 0.212671f + var_G * 0.715160f + var_B * 0.072169f;
//	p_node[a].XYZ_Z = var_R * 0.019334f + var_G * 0.119193f + var_B * 0.950227f;
//}
//
///* XYZ to Lab color space */
//void universe::XYZtoLab(int a) {
//	float var_X;
//	float var_Y;
//	float var_Z;
//
//	// Observer= 2¢X, Illuminant= D65
//	var_X = p_node[a].XYZ_X / ref_X;
//	var_Y = p_node[a].XYZ_Y / ref_Y;
//	var_Z = p_node[a].XYZ_Z / ref_Z;
//	if (var_X > 0.008856)
//		var_X = pow((double) var_X, 0.333333);
//	else
//		var_X = (7.787 * var_X) + (16 / 116);
//	if (var_Y > 0.008856)
//		var_Y = pow((double) var_Y, 0.333333);
//	else
//		var_Y = (7.787 * var_Y) + (16 / 116);
//	if (var_Z > 0.008856)
//		var_Z = pow((double) var_Z, 0.333333);
//	else
//		var_Z = (7.787 * var_Z) + (16 / 116);
//
//	p_node[a].Lab_L = (116 * var_Y) - 16;
//	p_node[a].Lab_a = 500 * (var_X - var_Y);
//	p_node[a].Lab_b = 200 * (var_Y - var_Z);
//}

// level start from 0


#endif /* DISJOINT_SET_H */
