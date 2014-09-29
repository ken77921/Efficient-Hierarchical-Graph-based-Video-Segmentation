/*
Original Code From:
Copyright (C) 2006 Pedro Felzenszwalb
Modifications (may have been made) Copyright (C) 2011,2012 
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

/* This is the main implementation of the streaming graph-based
 * hierarchical segmentation algorithm.
*/

#ifndef SEGMENT_IMAGE_H
#define SEGMENT_IMAGE_H

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include "unistd.h"

#include "image.h"
#include "misc.h"
#include "filter.h"

#include "edges.h"
//#include "segment-graph.h"
//#include "disjoint-set.h"
#include "disjoint-set_new.h"

using namespace std;




/* Gaussian Smoothing */
void smooth_images(image<rgb> *im[], int num_frame, image<float> *smooth_r[],
		image<float> *smooth_g[], image<float> *smooth_b[], float sigma) {

	int width = im[0]->width();
	int height = im[0]->height();

	image<float>** r = new image<float>*[num_frame];
	image<float>** g = new image<float>*[num_frame];
	image<float>** b = new image<float>*[num_frame];
	for (int i = 0; i < num_frame; i++) {
		r[i] = new image<float>(width, height);
		g[i] = new image<float>(width, height);
		b[i] = new image<float>(width, height);
	}
	for (int i = 0; i < num_frame; i++) {
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				imRef(r[i], x, y) = imRef(im[i], x, y).r;
				imRef(g[i], x, y) = imRef(im[i], x, y).g;
				imRef(b[i], x, y) = imRef(im[i], x, y).b;
			}
		}
	}

	// smooth each color channel
	for (int i = 0; i < num_frame; i++) {
		smooth_r[i] = smooth(r[i], sigma);
		smooth_g[i] = smooth(g[i], sigma);
		smooth_b[i] = smooth(b[i], sigma);
	}
	for (int i = 0; i < num_frame; i++) {
		delete r[i];
		delete g[i];
		delete b[i];
	}
	delete[] r;
	delete[] g;
	delete[] b;
}



/* main operation steps */
void segment_image(char *path, image<rgb> *im[], int num_frame, float c,
		float c_reg, int min_size, float sigma, int hie_num) {

	// step 1 -- Get information
	int width = im[0]->width();
	int height = im[0]->height();

	// ----- node number
	int num_vertices = num_frame * width * height;
	// ----- edge number
	int num_edges_plane = (width - 1) * (height - 1) * 2 + width * (height - 1)
			+ (width - 1) * height;
	int num_edges_layer = (width - 2) * (height - 2) * 9 + (width - 2) * 2 * 6
			+ (height - 2) * 2 * 6 + 4 * 4;
	int num_edges = num_edges_plane * num_frame
			+ num_edges_layer * (num_frame - 1);

	// ----- hierarchy setup
	//vector<vector<edge>*> edges_region;
	edge **edges_region=new edge*[hie_num + 1];
	//edges_region.resize(hie_num + 1);

	// ------------------------------------------------------------------

	// step 2 -- smooth images
	image<float>** smooth_r = new image<float>*[num_frame];
	image<float>** smooth_g = new image<float>*[num_frame];
	image<float>** smooth_b = new image<float>*[num_frame];
	smooth_images(im, num_frame, smooth_r, smooth_g, smooth_b, sigma);
	// ------------------------------------------------------------------

	// step 3 -- build edges
	printf("start build edges\n");
	printf("edge number: %d\n",num_edges);
	edges_region[0]=new edge[num_edges];
	//edge* edges = new edge[num_edges];
	initialize_edges(edges_region[0], num_frame, width, height, smooth_r, smooth_g,
			smooth_b);
	printf("end build edges\n");
	// ------------------------------------------------------------------

	// step 4 -- build nodes
	printf("start build nodes\n");
	universe* mess = new universe(num_frame, width, height, smooth_r, smooth_g,
			smooth_b, hie_num);
	printf("node number: %d\n",mess->node_num[0]);
	printf("end build nodes\n");
	// ------------------------------------------------------------------

	// step 5 -- over-segmentation
	printf("start over-segmentation\n");
	//edges_region[0] = new vector<edge>();
	//segment_graph(mess, edges_region[0], edges, num_edges, c, 0);
	edge_index edges_remain_index;
	//edges_remain_index.index=new unsigned int[num_edges];
	mess->segment_pixel_graph(&edges_remain_index, edges_region[0], num_edges, c, min_size);

	printf("end over-segmentation\n");
	// ------------------------------------------------------------------

	// step 6 -- hierarchical segmentation
	for (int i = 0; i < hie_num; i++) {
		printf("level = %d\n", i+1);

		// incremental in each hierarchy
		min_size = min_size * 1.2;

		printf("start update\n");
		mess->update(i);
		printf("node number: %d\n",mess->node_num[i+1]);
		printf("end update\n");

		printf("start fill edge weight\n");
		/*vector<edge> *edges_region_new = new vector<edge>();
		mess->fill_edge_weight(*edges_region_new,*edges_region[i], i);
		delete edges_region[i];
		edges_region[i]=edges_region_new;*/
		//edges_region[i+1]=new edge[edges_remain_index.num];
		edges_remain_index.num=mess->fill_edge_weight(edges_remain_index,edges_region[i],&edges_region[i+1], i);
		printf("remaining edge number: %d\n",edges_remain_index.num);
		delete [] edges_region[i];
		delete [] edges_remain_index.index;
		//edges_remain_index.index=new unsigned int[edges_remain_index.num];
		
		printf("end fill edge weight\n");

		printf("start segment graph region and merging min_size\n");
		//edges_region[i + 1] = new vector<edge>();
		mess->segment_graph_region(&edges_remain_index, edges_region[i+1],  c_reg, i + 1, min_size);
		printf("end segment graph region and merging min_size\n");

		c_reg = c_reg * 1.4;
		//delete edges_region[i];
	}
	delete edges_region[hie_num];
	// ------------------------------------------------------------------

	// step 8 -- generate output
	printf("start output\n");
	mess->generate_output(path, num_frame, width, height);
	printf("end output\n");
	// ------------------------------------------------------------------

	// step 9 -- clear everything
	delete mess;
	//delete[] edges;
	for (int i = 0; i < num_frame; i++) {
		delete smooth_r[i];
		delete smooth_g[i];
		delete smooth_b[i];
	}
	delete[] smooth_r;
	delete[] smooth_g;
	delete[] smooth_b;

}

#endif /* SEGMENT_IMAGE_H */
