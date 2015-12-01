
/*AAC CODE starting from line 151*/

/*
    pbrt source code Copyright(c) 1998-2012 Matt Pharr and Greg Humphreys.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */

// accelerators/bvh.cpp*
#include "stdafx.h"
#include "accelerators/bvh.h"
#include "probes.h"
#include "paramset.h"
#include "timer.h"
#include <limits>
const int DELTA=4;
const float EPSILON=0.2;

struct BVHPrimitiveInfo {
	BVHPrimitiveInfo() { }
	BVHPrimitiveInfo(int pn, const BBox &b)
		: primitiveNumber(pn), bounds(b) {
		centroid = .5f * b.pMin + .5f * b.pMax;
		morton = 0;
	}
	uint32_t morton;
	int primitiveNumber;
	Point centroid;
	BBox bounds;
};

struct BVHBuildNode {
    // BVHBuildNode Public Methods
    BVHBuildNode() { children[0] = children[1] = NULL; }
    void InitLeaf(uint32_t first, uint32_t n, const BBox &b) {
        firstPrimOffset = first;
        nPrimitives = n;
        bounds = b;
    }
	
	void InitInterior(uint32_t axis, BVHBuildNode *c0, BVHBuildNode *c1) {
        children[0] = c0;
        children[1] = c1;
        bounds = Union(c0->bounds, c1->bounds);
		centroid = .5f * bounds.pMin + .5f * bounds.pMax;
        splitAxis = axis;
        nPrimitives = 0;
    }
    BBox bounds;
	Point centroid;
    BVHBuildNode *children[2];
    uint32_t splitAxis, firstPrimOffset, nPrimitives;
};

struct CompareToMid {
    CompareToMid(int d, float m) { dim = d; mid = m; }
    int dim;
    float mid;
    bool operator()(const BVHPrimitiveInfo &a) const {
        return a.centroid[dim] < mid;
    }
};

struct ComparePoints {
    ComparePoints(int d) { dim = d; }
    int dim;
    bool operator()(const BVHPrimitiveInfo &a,
                    const BVHPrimitiveInfo &b) const {
        return a.centroid[dim] < b.centroid[dim];
    }
};

struct CompareToBucket {
    CompareToBucket(int split, int num, int d, const BBox &b)
        : centroidBounds(b)
    { splitBucket = split; nBuckets = num; dim = d; }
    bool operator()(const BVHPrimitiveInfo &p) const;

    int splitBucket, nBuckets, dim;
    const BBox &centroidBounds;
};

bool CompareToBucket::operator()(const BVHPrimitiveInfo &p) const {
    int b = nBuckets * ((p.centroid[dim] - centroidBounds.pMin[dim]) /
            (centroidBounds.pMax[dim] - centroidBounds.pMin[dim]));
    if (b == nBuckets) b = nBuckets-1;
    Assert(b >= 0 && b < nBuckets);
    return b <= splitBucket;
}

struct LinearBVHNode {
    BBox bounds;
    union {
        uint32_t primitivesOffset;    // leaf
        uint32_t secondChildOffset;   // interior
    };

    uint8_t nPrimitives;  // 0 -> interior node
    uint8_t axis;         // interior node: xyz
    uint8_t pad[2];       // ensure 32 byte total size
};

static inline bool IntersectP(const BBox &bounds, const Ray &ray,
        const Vector &invDir, const uint32_t dirIsNeg[3]) {
    // Check for ray intersection against $x$ and $y$ slabs
    float tmin =  (bounds[  dirIsNeg[0]].x - ray.o.x) * invDir.x;
    float tmax =  (bounds[1-dirIsNeg[0]].x - ray.o.x) * invDir.x;
    float tymin = (bounds[  dirIsNeg[1]].y - ray.o.y) * invDir.y;
    float tymax = (bounds[1-dirIsNeg[1]].y - ray.o.y) * invDir.y;
    if ((tmin > tymax) || (tymin > tmax))
        return false;
    if (tymin > tmin) tmin = tymin;
    if (tymax < tmax) tmax = tymax;

    // Check for ray intersection against $z$ slab
    float tzmin = (bounds[  dirIsNeg[2]].z - ray.o.z) * invDir.z;
    float tzmax = (bounds[1-dirIsNeg[2]].z - ray.o.z) * invDir.z;
    if ((tmin > tzmax) || (tzmin > tmax))
        return false;
    if (tzmin > tmin)
        tmin = tzmin;
    if (tzmax < tmax)
        tmax = tzmax;
    return (tmin < ray.maxt) && (tmax > ray.mint);
}

//morton code implemented from stack overflow
void mortonInit(vector<BVHPrimitiveInfo> &data, BBox global){
	float minx = global.pMin.x;
	float miny = global.pMin.y;
	float minz = global.pMin.z;
	float x = 1024.0f / (global.pMax.x - minx);
	float y = 1024.0f / (global.pMax.y - miny);
	float z = 1024.0f / (global.pMax.z - minz);
	int size = data.size();
	for (int i = 0; i < size; i++){
		Point p = data[i].centroid;
		uint32_t ix = ((1 << 10) - 1)&(uint32_t)round(x*(p.x - minx));
		uint32_t iy = ((1 << 10) - 1)&(uint32_t)round(y*(p.y - miny));
		uint32_t iz = ((1 << 10) - 1)&(uint32_t)round(z*(p.z - minz));
        
        ix = (ix | (ix << 16)) & 0x030000FF;
        ix = (ix | (ix <<  8)) & 0x0300F00F;
        ix = (ix | (ix <<  4)) & 0x030C30C3;
        ix = (ix | (ix <<  2)) & 0x09249249;
        
        iy = (iy | (iy << 16)) & 0x030000FF;
        iy = (iy | (iy <<  8)) & 0x0300F00F;
        iy = (iy | (iy <<  4)) & 0x030C30C3;
        iy = (iy | (iy <<  2)) & 0x09249249;
		
        iz = (iz | (iz << 16)) & 0x030000FF;
        iz = (iz | (iz <<  8)) & 0x0300F00F;
        iz = (iz | (iz <<  4)) & 0x030C30C3;
        iz = (iz | (iz <<  2)) & 0x09249249;
		
        data[i].morton = (ix << 2 | iy << 1 | iz);
	}
}

void mortonSort(vector<BVHPrimitiveInfo> &data, int start, int end, int b){
	if (end-start < 2)return;
	int zeros = 0;
	
    //see if sort is necessary
	for (int i = start; i < end; i++){
		if (((data[i]).morton & (1 << b)) == 0) zeros++;
	}
	
    //if sort is unnecessary continue to next bucket
	if (zeros == 0 || zeros == end - start){
		if (b>0)mortonSort(data, start, end, b - 1);
		return;
	}
	
    //sort by bit b
	int zerobound = start; int onebound = end-1;
	while (onebound != zerobound){
		if ((data[zerobound]).morton & (1 << b)){
			swap(data[zerobound], data[onebound]);
			onebound--;
		}
		else{ zerobound++; }
	}
	mortonSort(data, start, start + zeros, b - 1);
	mortonSort(data, start + zeros, end, b - 1);
}

inline int func(int i){ 
    return (int)ceil(0.5*pow(DELTA, 0.5 + EPSILON)*pow(i, 0.5 - EPSILON)); 
}

vector<BVHBuildNode*> BVHAccel::aacBuild(MemoryArena &buildArena,
    vector<BVHPrimitiveInfo> &buildData, vector<Reference<Primitive> > &orderedPrims,
	uint32_t start, uint32_t end,
	uint32_t bit, uint32_t *totalNodes){
	vector<BVHBuildNode*> cluster;
    if (end - start < DELTA){
        cluster.reserve(end-start);	
        //initialize C
		(*totalNodes) += end - start;
		for (int i = start; i < end; i++){
			uint32_t primNum = buildData[i].primitiveNumber;
			BVHBuildNode *node = buildArena.Alloc<BVHBuildNode>();
			node->InitLeaf(orderedPrims.size(), 1, buildData[i].bounds);
            cluster.push_back(node);
			orderedPrims.push_back(primitives[primNum]);
		}
		combine(buildArena,cluster, func(DELTA),totalNodes);
		return cluster;
	}
	int changepoint = 0;
    if(((buildData[start].morton ^ buildData[end-1].morton) & (1 << bit))==0){
        changepoint = (start + end)/2 +1;
    }else{
        for (int i = start; i < end; i++){
            if ((buildData[i].morton & (1 << bit))!=0){
               changepoint = i;
                break;
            }
        }
    }
	cluster = 
        aacBuild(buildArena, buildData, orderedPrims, start, changepoint, bit-1, totalNodes);
	vector<BVHBuildNode*>rightCluster = 
        aacBuild(buildArena, buildData, orderedPrims, changepoint, end, bit-1, totalNodes);
    
    cluster.insert(cluster.end(), rightCluster.begin(), rightCluster.end());
    combine(buildArena, cluster, func(end-start),totalNodes);
    return cluster;
}

struct mdist{
    int index;
    float distance;
};

//finds best match in all clusters after index
void initMatch(vector<BVHBuildNode*> cluster, /*vector<vector<float> >&pairs,*/ vector<mdist>&min, int index){
    mdist result={};
    int size = cluster.size();
    //vector<float> vec;
    //vec.reserve(size-index-1);
    int mini = -1;
    float mind = std::numeric_limits<float>::max();
    BBox box = cluster[index]->bounds;
    for(int i = index +1; i <size;i++){
        BBox cmp = Union(box,cluster[i]->bounds);
        float dist = cmp.SurfaceArea();
    //    vec.push_back(dist);
        if(dist < mind){
            mini = i;
            mind = dist;
        }
    }
    //pairs.push_back(vec);
    result.index = mini;
    result.distance = mind;
    min.push_back(result);
}

void BVHAccel::findBestMatch(vector<BVHBuildNode*> &cluster, vector<mdist> &mins, int index){
    BBox box = cluster[index]->bounds;
    int mini = -1;
    float mind = std::numeric_limits<float>::max();
    for(int i = index+1; i < cluster.size();i++){
        BBox b = Union(box, cluster[i]->bounds);
        float d = b.SurfaceArea();
        if(d <mind){mind = d; mini = i;}
    }
    mins[index].index= mini;
    mins[index].distance= mind;
}

void BVHAccel::combine(MemoryArena &buildArena, vector<BVHBuildNode*> &cluster, int nNodes, uint32_t *totalNodes){
    int size = cluster.size();
	if(size <= 1 || size <= nNodes)return;
    vector<mdist> min;
    min.reserve(size-1);
    //vector<vector<float> > pairs;
    //pairs.reserve(size-1);

    for(int i = 0; i <size-1; i++){
        initMatch(cluster,/* pairs,*/ min, i);
    }//n-1 min values
    
    //loop and cluster
	while (size >nNodes){
        int d= size-2;//check min of all n-1 pairs
        for(int i = 0; i <size-2;i++){
            if(min[i].distance < min[d].distance) d = i;
        }
		//merge
        int left = d;
        int right = min[d].index;
        (*totalNodes)++;
		BVHBuildNode *node = buildArena.Alloc<BVHBuildNode>();
		node->InitInterior(0, cluster[left], cluster[right]);
        cluster[left]=node;
        cluster[right]=cluster.back();
    //    moveEnd(cluster, pairs, min, right);
    //   updateLeft(cluster, pairs, min, left);
    //   update(cluster, pairs, min, left, right);
        cluster.pop_back();
        min.pop_back();
        findBestMatch(cluster,min,left);
        if(right <size-1)findBestMatch(cluster,min,right);
        for(int i = 0; i <size-2;i++){
            if(min[i].index ==left ||min[i].index ==right)
                findBestMatch(cluster,min,i);
            else if( min[i].index == size-1 ){
                if(i < right) min[i].index = right;
                findBestMatch(cluster,min,i);
            }
        }
        size --;
    }
}

//optimization using double array that didn't end up speeding it up by much
/*
void updateLeft(vector<BVHBuildNode*> &cluster, vector<vector<float> > &pairs, vector<mdist> &min, int left){
    int size = cluster.size();//currently old size
    //original matrix[a][b] flattened [a][b-a-1]
    //update other nodes for combined node in left
    BBox box = cluster[left]->bounds;
    for(int i = 0; i < left;i++){
        BBox b = Union(cluster[i]->bounds,box);
        pairs[i][left-i-1]= b.SurfaceArea();
    }
    //update new left node
    int mini = -1;
    float mind = std::numeric_limits<float>::max();
    for(int i = left+1; i < size-1;i++){
        BBox b = Union(cluster[i]->bounds,box);
        pairs[left][i-left-1]= b.SurfaceArea();
        if(pairs[left][i-left-1] < mind){
            mind = pairs[left][i-left-1];
            mini = i;
        }
    }
    if(left < size -1){
        min[left].index = mini;
        min[left].distance = mind; 
    } 
}
void moveEnd(vector<BVHBuildNode*> &cluster, vector<vector<float> > &pairs, vector<mdist> &min, int right){
    int size = cluster.size();//currently old size
    //original matrix[a][b] flattened [a][b-a-1]
    //update other nodes for moving of end node to right node
    for(int i = 0; i < right;i++){
        pairs[i][right-i-1]=pairs[i][(size-1)-i-1];
        pairs[i].pop_back();
    }
    //update new right node
    if(right < size-1)pairs[right].pop_back();
    int mini = -1;
    int mind = std::numeric_limits<float>::max();
    for(int i = right+1; i < size-1;i++){
        pairs[right][i-right-1]=pairs[i][(size-1)-i-1];
        pairs[i].pop_back();
        if(pairs[right][i-right-1] < mind){
            mind = pairs[right][i-right-1];
            mini = i;
        }
    }
    if(right < size-2){
        min[right].index = mini;
        min[right].distance = mind; 
    }
    min.pop_back();
}
void update(vector<BVHBuildNode*> &cluster, vector<vector<float> > &pairs, vector<mdist> &min, int left, int right){
    int size = cluster.size();//currently old size
     //update any nodes with original left/right as mins
    for(int i = 0; i < size -1; i++){
        int prev = min[i].index;
        if(prev ==right || prev == left){
            int ind = -1;
            float dst = std::numeric_limits<float>::max();
            for(int j = i+1; j < size-1;j++){
                if(pairs[i][j-i-1] < dst){
                    dst = pairs[i][j-i-1];
                    ind= j;
                }
             }
            min[i].index = ind;
            min[i].distance = dst; 
            }
        else if(prev == size-1) min[i].index = right;
        }
}
*/


// BVHAccel Method Definitions
BVHAccel::BVHAccel(const vector<Reference<Primitive> > &p,
                   uint32_t mp, const string &sm) {
    maxPrimsInNode = min(255u, mp);
    for (uint32_t i = 0; i < p.size(); ++i)
        p[i]->FullyRefine(primitives);
    if (sm == "sah")         splitMethod = SPLIT_SAH;
    else if (sm == "middle") splitMethod = SPLIT_MIDDLE;
    else if (sm == "equal")  splitMethod = SPLIT_EQUAL_COUNTS;
    else if (sm == "aac")  splitMethod = SPLIT_AAC;
    else {
        Warning("BVH split method \"%s\" unknown.  Using \"sah\".",
                sm.c_str());
        splitMethod = SPLIT_SAH;
    }

    if (primitives.size() == 0) {
        nodes = NULL;
        return;
    }
    // Build BVH from _primitives_
    PBRT_BVH_STARTED_CONSTRUCTION(this, primitives.size());
    // Initialize _buildData_ array for primitives
    vector<BVHPrimitiveInfo> buildData;
    buildData.reserve(primitives.size());
	BBox global = primitives[0]->WorldBound();
    for (uint32_t i = 0; i < primitives.size(); ++i) {
        BBox bbox = primitives[i]->WorldBound();
		global = Union(global, bbox);
        buildData.push_back(BVHPrimitiveInfo(i, bbox));
    }

    // Recursively build BVH tree for primitives
    MemoryArena buildArena;
    uint32_t totalNodes = 0;
    vector<Reference<Primitive> > orderedPrims;
    orderedPrims.reserve(primitives.size());
	BVHBuildNode *root;
	
    if (splitMethod == SPLIT_AAC){
        mortonInit(buildData,global);
        mortonSort(buildData,0,primitives.size(),29);
        vector<BVHBuildNode*>cluster = 
            aacBuild(buildArena, buildData, orderedPrims, 0, primitives.size(), 29, &totalNodes);
		combine(buildArena,cluster,1,&totalNodes);
		root = cluster[0];
    }
	else{
		root = recursiveBuild(buildArena, buildData, 0,
			primitives.size(), &totalNodes,
			orderedPrims);
	}
    primitives.swap(orderedPrims);
        Info("BVH created with %d nodes for %d primitives (%.2f MB)", totalNodes,
             (int)primitives.size(), float(totalNodes * sizeof(LinearBVHNode))/(1024.f*1024.f));

    // Compute representation of depth-first traversal of BVH tree
    nodes = AllocAligned<LinearBVHNode>(totalNodes);
    for (uint32_t i = 0; i < totalNodes; ++i)
        new (&nodes[i]) LinearBVHNode;
    uint32_t offset = 0;
    flattenBVHTree(root, &offset);
    Assert(offset == totalNodes);
    PBRT_BVH_FINISHED_CONSTRUCTION(this);
}



BBox BVHAccel::WorldBound() const {
    return nodes ? nodes[0].bounds : BBox();
}


BVHBuildNode *BVHAccel::recursiveBuild(MemoryArena &buildArena,
        vector<BVHPrimitiveInfo> &buildData, uint32_t start,
        uint32_t end, uint32_t *totalNodes,
        vector<Reference<Primitive> > &orderedPrims) {
    Assert(start != end);
    (*totalNodes)++;
    BVHBuildNode *node = buildArena.Alloc<BVHBuildNode>();
    // Compute bounds of all primitives in BVH node
    BBox bbox;
    for (uint32_t i = start; i < end; ++i)
        bbox = Union(bbox, buildData[i].bounds);
    uint32_t nPrimitives = end - start;
    if (nPrimitives == 1) {
        // Create leaf _BVHBuildNode_
        uint32_t firstPrimOffset = orderedPrims.size();
        for (uint32_t i = start; i < end; ++i) {
            uint32_t primNum = buildData[i].primitiveNumber;
            orderedPrims.push_back(primitives[primNum]);
        }
        node->InitLeaf(firstPrimOffset, nPrimitives, bbox);
    }
    else {
        // Compute bound of primitive centroids, choose split dimension _dim_
        BBox centroidBounds;
        for (uint32_t i = start; i < end; ++i)
            centroidBounds = Union(centroidBounds, buildData[i].centroid);
        int dim = centroidBounds.MaximumExtent();

        // Partition primitives into two sets and build children
        uint32_t mid = (start + end) / 2;
        if (centroidBounds.pMax[dim] == centroidBounds.pMin[dim]) {
            // If nPrimitives is no greater than maxPrimsInNode,
            // then all the nodes can be stored in a compact bvh node.
            if (nPrimitives <= maxPrimsInNode) {
                // Create leaf _BVHBuildNode_
                uint32_t firstPrimOffset = orderedPrims.size();
                for (uint32_t i = start; i < end; ++i) {
                    uint32_t primNum = buildData[i].primitiveNumber;
                    orderedPrims.push_back(primitives[primNum]);
                }
                node->InitLeaf(firstPrimOffset, nPrimitives, bbox);
                return node;
            }
            else {
                // else if nPrimitives is greater than maxPrimsInNode, we
                // need to split it further to guarantee each node contains
                // no more than maxPrimsInNode primitives.
                node->InitInterior(dim,
                                   recursiveBuild(buildArena, buildData, start, mid,
                                                  totalNodes, orderedPrims),
                                   recursiveBuild(buildArena, buildData, mid, end,
                                                  totalNodes, orderedPrims));
                return node;
            }
        }

        // Partition primitives based on _splitMethod_
        switch (splitMethod) {
        case SPLIT_MIDDLE: {
            // Partition primitives through node's midpoint
            float pmid = .5f * (centroidBounds.pMin[dim] + centroidBounds.pMax[dim]);
            BVHPrimitiveInfo *midPtr = std::partition(&buildData[start],
                                                      &buildData[end-1]+1,
                                                      CompareToMid(dim, pmid));
            mid = midPtr - &buildData[0];
            if (mid != start && mid != end)
                // for lots of prims with large overlapping bounding boxes, this
                // may fail to partition; in that case don't break and fall through
                // to SPLIT_EQUAL_COUNTS
                break;
        }
        case SPLIT_EQUAL_COUNTS: {
            // Partition primitives into equally-sized subsets
            mid = (start + end) / 2;
            std::nth_element(&buildData[start], &buildData[mid],
                             &buildData[end-1]+1, ComparePoints(dim));
            break;
        }
        case SPLIT_SAH: default: {
            // Partition primitives using approximate SAH
            if (nPrimitives <= 4) {
                // Partition primitives into equally-sized subsets
                mid = (start + end) / 2;
                std::nth_element(&buildData[start], &buildData[mid],
                                 &buildData[end-1]+1, ComparePoints(dim));
            }
            else {
                // Allocate _BucketInfo_ for SAH partition buckets
                const int nBuckets = 12;
                struct BucketInfo {
                    BucketInfo() { count = 0; }
                    int count;
                    BBox bounds;
                };
                BucketInfo buckets[nBuckets];

                // Initialize _BucketInfo_ for SAH partition buckets
                for (uint32_t i = start; i < end; ++i) {
                    int b = nBuckets *
                        ((buildData[i].centroid[dim] - centroidBounds.pMin[dim]) /
                         (centroidBounds.pMax[dim] - centroidBounds.pMin[dim]));
                    if (b == nBuckets) b = nBuckets-1;
                    Assert(b >= 0 && b < nBuckets);
                    buckets[b].count++;
                    buckets[b].bounds = Union(buckets[b].bounds, buildData[i].bounds);
                }

                // Compute costs for splitting after each bucket
                float cost[nBuckets-1];
                for (int i = 0; i < nBuckets-1; ++i) {
                    BBox b0, b1;
                    int count0 = 0, count1 = 0;
                    for (int j = 0; j <= i; ++j) {
                        b0 = Union(b0, buckets[j].bounds);
                        count0 += buckets[j].count;
                    }
                    for (int j = i+1; j < nBuckets; ++j) {
                        b1 = Union(b1, buckets[j].bounds);
                        count1 += buckets[j].count;
                    }
                    cost[i] = .125f + (count0*b0.SurfaceArea() + count1*b1.SurfaceArea()) /
                              bbox.SurfaceArea();
                }

                // Find bucket to split at that minimizes SAH metric
                float minCost = cost[0];
                uint32_t minCostSplit = 0;
                for (int i = 1; i < nBuckets-1; ++i) {
                    if (cost[i] < minCost) {
                        minCost = cost[i];
                        minCostSplit = i;
                    }
                }

                // Either create leaf or split primitives at selected SAH bucket
                if (nPrimitives > maxPrimsInNode ||
                    minCost < nPrimitives) {
                    BVHPrimitiveInfo *pmid = std::partition(&buildData[start],
                        &buildData[end-1]+1,
                        CompareToBucket(minCostSplit, nBuckets, dim, centroidBounds));
                    mid = pmid - &buildData[0];
                }
                
                else {
                    // Create leaf _BVHBuildNode_
                    uint32_t firstPrimOffset = orderedPrims.size();
                    for (uint32_t i = start; i < end; ++i) {
                        uint32_t primNum = buildData[i].primitiveNumber;
                        orderedPrims.push_back(primitives[primNum]);
                    }
                    node->InitLeaf(firstPrimOffset, nPrimitives, bbox);
                    return node;
                }
            }
            break;
        }
        }
        node->InitInterior(dim,
                           recursiveBuild(buildArena, buildData, start, mid,
                                          totalNodes, orderedPrims),
                           recursiveBuild(buildArena, buildData, mid, end,
                                          totalNodes, orderedPrims));
    }
    return node;
}


uint32_t BVHAccel::flattenBVHTree(BVHBuildNode *node, uint32_t *offset) {
    LinearBVHNode *linearNode = &nodes[*offset];
    linearNode->bounds = node->bounds;
    uint32_t myOffset = (*offset)++;
    if (node->nPrimitives > 0) {
        Assert(!node->children[0] && !node->children[1]);
        linearNode->primitivesOffset = node->firstPrimOffset;
        linearNode->nPrimitives = node->nPrimitives;
    }
    else {
        // Creater interior flattened BVH node
        linearNode->axis = node->splitAxis;
        linearNode->nPrimitives = 0;
        flattenBVHTree(node->children[0], offset);
        linearNode->secondChildOffset = flattenBVHTree(node->children[1],
                                                       offset);
    }
    return myOffset;
}


BVHAccel::~BVHAccel() {
    FreeAligned(nodes);
}

bool BVHAccel::Intersect(const Ray &ray, Intersection *isect) const {
    if (!nodes) return false;
    PBRT_BVH_INTERSECTION_STARTED(const_cast<BVHAccel *>(this), const_cast<Ray *>(&ray));
    bool hit = false;
    Vector invDir(1.f / ray.d.x, 1.f / ray.d.y, 1.f / ray.d.z);
    uint32_t dirIsNeg[3] = { invDir.x < 0, invDir.y < 0, invDir.z < 0 };
    // Follow ray through BVH nodes to find primitive intersections
    uint32_t todoOffset = 0, nodeNum = 0;
    uint32_t todo[64];
    while (true) {
        const LinearBVHNode *node = &nodes[nodeNum];
        // Check ray against BVH node
        if (::IntersectP(node->bounds, ray, invDir, dirIsNeg)) {
            if (node->nPrimitives > 0) {
                // Intersect ray with primitives in leaf BVH node
                PBRT_BVH_INTERSECTION_TRAVERSED_LEAF_NODE(const_cast<LinearBVHNode *>(node));
                for (uint32_t i = 0; i < node->nPrimitives; ++i)
                {
                    PBRT_BVH_INTERSECTION_PRIMITIVE_TEST(const_cast<Primitive *>(primitives[node->primitivesOffset+i].GetPtr()));
                    if (primitives[node->primitivesOffset+i]->Intersect(ray, isect))
                    {
                        PBRT_BVH_INTERSECTION_PRIMITIVE_HIT(const_cast<Primitive *>(primitives[node->primitivesOffset+i].GetPtr()));
                        hit = true;
                    }
                    else {
                        PBRT_BVH_INTERSECTION_PRIMITIVE_MISSED(const_cast<Primitive *>(primitives[node->primitivesOffset+i].GetPtr()));
                   }
                }
                if (todoOffset == 0) break;
                nodeNum = todo[--todoOffset];
            }
            else {
                // Put far BVH node on _todo_ stack, advance to near node
                PBRT_BVH_INTERSECTION_TRAVERSED_INTERIOR_NODE(const_cast<LinearBVHNode *>(node));
                if (dirIsNeg[node->axis]) {
                   todo[todoOffset++] = nodeNum + 1;
                   nodeNum = node->secondChildOffset;
                }
                else {
                   todo[todoOffset++] = node->secondChildOffset;
                   nodeNum = nodeNum + 1;
                }
            }
        }
        else {
            if (todoOffset == 0) break;
            nodeNum = todo[--todoOffset];
        }
    }
    PBRT_BVH_INTERSECTION_FINISHED();
    return hit;
}

bool BVHAccel::IntersectP(const Ray &ray) const {
    if (!nodes) return false;
    PBRT_BVH_INTERSECTIONP_STARTED(const_cast<BVHAccel *>(this), const_cast<Ray *>(&ray));
    Vector invDir(1.f / ray.d.x, 1.f / ray.d.y, 1.f / ray.d.z);
    uint32_t dirIsNeg[3] = { invDir.x < 0, invDir.y < 0, invDir.z < 0 };
    uint32_t todo[64];
    uint32_t todoOffset = 0, nodeNum = 0;
    while (true) {
        const LinearBVHNode *node = &nodes[nodeNum];
        if (::IntersectP(node->bounds, ray, invDir, dirIsNeg)) {
            // Process BVH node _node_ for traversal
            if (node->nPrimitives > 0) {
                PBRT_BVH_INTERSECTIONP_TRAVERSED_LEAF_NODE(const_cast<LinearBVHNode *>(node));
                  for (uint32_t i = 0; i < node->nPrimitives; ++i) {
                    PBRT_BVH_INTERSECTIONP_PRIMITIVE_TEST(const_cast<Primitive *>(primitives[node->primitivesOffset + i].GetPtr()));
                    if (primitives[node->primitivesOffset+i]->IntersectP(ray)) {
                        PBRT_BVH_INTERSECTIONP_PRIMITIVE_HIT(const_cast<Primitive *>(primitives[node->primitivesOffset+i].GetPtr()));
                        return true;
                    }
                else {
                        PBRT_BVH_INTERSECTIONP_PRIMITIVE_MISSED(const_cast<Primitive *>(primitives[node->primitivesOffset + i].GetPtr()));
                    }
                }
                if (todoOffset == 0) break;
                nodeNum = todo[--todoOffset];
            }
            else {
                PBRT_BVH_INTERSECTIONP_TRAVERSED_INTERIOR_NODE(const_cast<LinearBVHNode *>(node));
                if (dirIsNeg[node->axis]) {
                   /// second child first
                   todo[todoOffset++] = nodeNum + 1;
                   nodeNum = node->secondChildOffset;
                }
                else {
                   todo[todoOffset++] = node->secondChildOffset;
                   nodeNum = nodeNum + 1;
                }
            }
        }
        else {
            if (todoOffset == 0) break;
            nodeNum = todo[--todoOffset];
        }
    }
    PBRT_BVH_INTERSECTIONP_FINISHED();
    return false;
}


BVHAccel *CreateBVHAccelerator(const vector<Reference<Primitive> > &prims,
        const ParamSet &ps) {
    string splitMethod = ps.FindOneString("splitmethod", "sah");
    uint32_t maxPrimsInNode = ps.FindOneInt("maxnodeprims", 4);
    return new BVHAccel(prims, maxPrimsInNode, splitMethod);
}


