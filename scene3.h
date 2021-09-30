/***************************************************************************
 *
 * krt - Kens Raytracer - Coursework Edition. (C) Copyright 1997-2018.
 *
 * Do what you like with this code as long as you retain this comment.
 */


#include "photon.h"
#include "colour.h"
#include "ray.h"
#include "object.h"
#include "light.h"
#include "hit.h"
#include "vector.h"
#include "vertex.h"
#include "nanoflann-master/include/nanoflann.hpp"

using namespace nanoflann;
struct PointCloud{
	struct Point{            
		float  x,y,z;
		Vector incomingray;
		Colour intensity; 
	};

	std::vector<Point>  nodes;
	inline size_t kdtree_get_point_count() const { return nodes.size(); }
	inline float kdtree_get_pt(const size_t idx, const size_t dim) const{
		if (dim == 0) return nodes[idx].x;
		else if (dim == 1) return nodes[idx].y;
		else return nodes[idx].z;
	}
	template <class BBOX>
	bool kdtree_get_bbox(BBOX&) const { return false; }

};

class Scene3 {
public:
	typedef KDTreeSingleIndexAdaptor< L2_Simple_Adaptor<float, PointCloud > , PointCloud,	3 > kdTree;
	
	kdTree *treeN;
	kdTree *treeC;
	Object *objects;
	PointCloud nMap;
	PointCloud cMap;
	Scene3();
	void trace(Ray ray, Object *objects, Hit &hit);
    void photonTracer(Photon photon, int level, int depthD, bool depthS, bool caustic);
	void raytrace(Ray ray, Vertex camera, Colour &colour, Vertex l, int level, bool caustic);
	Vector reflection(Vector direction, Vector normal);
	double clamp(double x, double upper, double lower);
	Vector transparent(Vector direction, Vector N, float ior);
	float fresnel(Vector &direction, Vector &normal, float &ior,float &kr);
};
