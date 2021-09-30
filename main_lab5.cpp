/***************************************************************************
 *
 * krt - Kens Raytracer - Coursework Edition. (C) Copyright 1997-2019.
 *
 * Do what you like with this code as long as you retain this comment.
 */

/* This is the entry point function for the program you need to create for lab two.
 * You should not need to modify this code.
 * It creates a framebuffer, loads an triangle mesh object, calls the drawing function to render the object and then outputs the framebuffer as a ppm file.
 *
 * On linux.bath.ac.uk:
 *
 * Compile the code using g++ -o lab2executable main_lab2.cpp framebuffer.cpp linedrawer.cpp polymesh.cpp -lm
 *
 * Execute the code using ./lab2executable
 *
 * This will produce an image file called test.ppm. You can convert this a png file for viewing using
 *
 * pbmropng test.ppm > test.png
 *
 * You are expected to fill in the missing code in polymesh.cpp.
 */



/*
to run:

cd path to file
g++ -o lab main_lab5.cpp framebuffer.cpp directional_light.cpp phong.cpp polymesh.cpp scene3.cpp sphere.cpp -lm
./lab

*/


#include "framebuffer.h"
#include "ray.h"
#include "hit.h"
#include "polymesh.h"
#include "sphere.h"
#include "light.h"
#include "directional_light.h"
#include "material.h"
#include "phong.h"
#include "math.h"
#include "scene3.h"
#include "polymesh2.h"


#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include "nanoflann-master/include/nanoflann.hpp"
using namespace nanoflann;
using namespace std;

void object_test(Ray ray, Object *objects, Hit &best_hit){
  Object *obj = objects;
  best_hit.flag = false;
  while(obj != 0){
    Hit obj_hit;
    obj_hit.flag=false;
    obj->intersection(ray, obj_hit);
    if (obj_hit.flag){
      if (obj_hit.t > 0.0f){
        if (best_hit.flag == false){
	        best_hit = obj_hit;
	      } else if (obj_hit.t < best_hit.t){
	        best_hit = obj_hit;
      	}
      }
    }
    obj = obj->next;
  }
  return;
}


int main(int argc, char *argv[])
{
  int width = 1012;
  int height = 1012;
  // Create a framebuffer
  FrameBuffer *fb = new FrameBuffer(width,height);

  // The following transform allows 4D homogeneous coordinates to be transformed. It moves the supplied teapot model to somewhere visible.
  Transform *transform = new Transform(1.0f, 0.0f, 0.0f, 0.0f,
				                               0.0f, 0.0f, 1.0f,-2.7f,
				                               0.0f, 1.0f, 0.0f, 5.0f,
				                               0.0f, 0.0f, 0.0f, 1.0f);
  
  Transform *transform3 = new Transform(1.0f, 0.0f, 0.0f, 1.3f,
				                                0.0f, 0.0f, 1.0f, 1.4f,
				                                0.0f, 1.0f, 0.0f, 5.0f,
				                                0.0f, 0.0f, 0.0f, 1.0f);
  
  Transform *transform22 = new Transform(1.0f, 0.0f, 0.0f, 1.0f,
                                        0.0f, 0.0f, 1.0f,-5.5f,
                                        0.0f, 1.0f, 0.0f, 5.5f,
                                        0.0f, 0.0f, 0.0f, 1.0f);

  Transform *transform2 = new Transform(1.2f, 0.0f, 0.0f, -0.5f,//-1.5
				                                0.0f, 0.0f, 1.2f,-6.0f,
				                                0.0f, 1.2f, 0.0f, 5.0f,
				                                0.0f, 0.0f, 0.0f, 1.0f);

  Transform *beeth_transform = new Transform(-0.3f, 0.0f, 0.0f, 2.8f,
                                              0.0f, 0.3f, 0.0f,-2.5f,
                                              0.0f, 0.0f,-0.3f,-3.5f,//0.0 o pune invers si deasupra// 1.0 o pune in spate
                                              0.0f, 0.0f, 0.0f, 1.0f); 

  //  Read in the teapot model.
  // PolyMesh2 *pm2 = new PolyMesh2((char *)"ply/teapot.ply", transform22);
  PolyMesh *pm = new PolyMesh((char *)"ply/teapot_smaller.ply", transform2);
  PolyMesh *floor = new PolyMesh((char *)"ply/floor.ply", transform);
  PolyMesh *wallright = new PolyMesh((char *)"ply/wallright.ply", transform);
  PolyMesh *wallleft = new PolyMesh((char *)"ply/wallleft.ply", transform);
  PolyMesh *ceilingbackfloor = new PolyMesh((char *)"ply/ceilingbackfloor.ply", transform);
  // PolyMesh *beeth = new PolyMesh((char *)"ply/beethoven.ply", beeth_transform);
 
  Vertex v1;
  v1.x = 3.6f;  
  v1.y = -0.7f;  
  v1.z = 7.4f; //big reflective
  Sphere *sphere1 = new Sphere(v1,1.0f);

  Vertex v2;
  v2.x = 4.3f;  
  v2.y = 2.0f;  
  v2.z = 6.5f; //middle underneathreflective
  Sphere *sphere2 = new Sphere(v2,1.5f);
  
  Vertex v3;
  v3.x = 2.3f;  
  v3.y = 1.0f;  
  v3.z = 5.5f; //small reflective
  Sphere *sphere3 = new Sphere(v3,0.7f);

  Vertex v4;
  v4.x = 4.3f;  
  v4.y = -4.9f;  
  v4.z = 3.5f; // big transparent
  Sphere *sphere4 = new Sphere(v4,1.7f);//was 1.5

  Vertex v5;
  v5.x = 2.2f;  
  v5.y = -3.2f;  
  v5.z = 1.2f; //small transarent
  Sphere *sphere5 = new Sphere(v5,0.5f);
  
  Vertex v6;
  v6.x = 2.0f;  
  v6.y = -5.5f;  
  v6.z = 1.9f; //middle transparent
  Sphere *sphere6 = new Sphere(v6,1.0f);
 

  Phong red;//red
  red.ambient = red.a; 
	red.d.r = 1.0f;
  red.d.g = 0.0f;
  red.d.b = 0.0f;
  red.diffuse = red.d;
  red.s.r = 1.0f;
  red.s.g = 1.0f;
  red.s.b = 1.0f;
  red.specular = red.s;
  red.power = 40.0f;
  red.ior = 1.0f;
  
  

  Phong greenwall;//green
  greenwall.ambient = greenwall.a; 
	greenwall.d.r = 0.0f;
  greenwall.d.g = 1.0f;
  greenwall.d.b = 0.0f;
  greenwall.diffuse = greenwall.d;
  greenwall.s.r = 1.0f;
  greenwall.s.g = 1.0f;
  greenwall.s.b = 1.0f;
  greenwall.specular = greenwall.s;
  greenwall.power = 40.0f;
  greenwall.ior = 1.0f;
  

  Phong whitewall;
  whitewall.ambient = whitewall.a; 
	whitewall.d.r = 1.0f;
  whitewall.d.g = 1.0f;
  whitewall.d.b = 1.0f;
  whitewall.diffuse = whitewall.d;
  whitewall.s.r = 1.0f;
  whitewall.s.g = 1.0f;
  whitewall.s.b = 1.0f;
  whitewall.specular = whitewall.s;
  whitewall.power = 40.0f;
  whitewall.ior = 1.0f;
  
    
  Phong metal;
  metal.a.r = 0.33f;
  metal.a.g = 0.33f;
  metal.a.b = 0.33f;
  metal.ambient = metal.a; 
	metal.d.r = 0.66f;
  metal.d.g = 0.66f;
  metal.d.b = 0.66f;
  metal.diffuse = metal.d;
  metal.s.r = 1.0f;
  metal.s.g = 1.0f;
  metal.s.b = 1.0f;
  metal.specular = metal.s;
  metal.power = 40.0f;
  metal.ior = 2.5f;
  metal.reflect = true;


  Phong water;
  water.ambient = water.a; 
	water.d.r = 0.0f;
  water.d.g = 0.0f;
  water.d.b = 0.0f;
  water.diffuse = water.d;
  water.s.r = 1.0f;
  water.s.g = 1.0f;
  water.s.b = 1.0f;
  water.specular = water.s;
  water.power = 40.0f;
  water.ior = 1.33f;
  water.reflect = true; 
  water.refract = true;

  Phong glass;
  // glass.a.r = 1.0f;
  // glass.a.g = 0.5;
  // glass.a.b = 0.9;
  // glass.d.r = 1.0f;
  // glass.d.g = 0.5f;
  // glass.d.b = 0.9f;
  glass.ambient = glass.a; 
	glass.d.r = 0.0f;
  glass.d.g = 0.0f;
  glass.d.b = 0.0f;
  glass.diffuse = water.d;
  glass.s.r = 1.0f;
  glass.s.g = 1.0f;
  glass.s.b = 1.0f;
  glass.specular = glass.s;
  glass.power = 40.0f;
  glass.ior = 1.5f;
  glass.reflect = true; 
  glass.refract = true; 
  
  sphere1 -> material = &metal;
  sphere2 -> material = &metal;
  sphere3 -> material = &metal;
  sphere4 -> material = &water; //was glass
  sphere5 -> material = &water;
  sphere6 -> material = &water;

 
  pm -> material = &water;


  floor->material = &whitewall;
  wallright ->material = &red;
  wallleft ->material = &greenwall;
  ceilingbackfloor ->material = &whitewall;
  
  floor -> next = wallright;
  wallright -> next = wallleft;
  wallleft -> next = ceilingbackfloor;
  ceilingbackfloor -> next = sphere1;

  sphere1 -> next = sphere2;
  sphere2 -> next = sphere3;
  sphere3 -> next = sphere4;
  sphere4 -> next = sphere5;
  sphere5 -> next = sphere6;

  sphere6 -> next = pm;
  // ceilingbackfloor -> next = pm;



  Vertex l; 
  l.x = 1.0f; l.y = 4.5f; l.z = 1.0f;//3
  DirectionalLight *dl = new DirectionalLight(Vector(1.01f, 5.0f, -1.0f),Colour(1.0f, 1.0f, 1.0f, 0.0f));
  Scene3 *scene = new Scene3();
  scene->objects = floor;


  int photonsc = 300000;
  int photons = 450000;
  int level = 7;

  float x,y,z;
  for(int i=0;i<photonsc;i++){ //fire #photonsc in the scene not up	
    x = ((float)rand() / (RAND_MAX))*2 - 1;
    y = -abs(((float)rand() / (RAND_MAX))*2 - 1);
    z = ((float)rand() / (RAND_MAX))*2 - 1;
    Photon cp;
    cp.ray.position = l; // set attributes for photon cp
    cp.ray.direction = Vector(x,y,z);
    cp.ray.direction.normalise();
    float intensity =1.0f/(float)photonsc; 
    int depthD = 1;
    bool depthS = true;
    bool isCaustic = true;
    scene->photonTracer(cp,level,depthD,depthS,isCaustic); 
  }

  for(int i=0;i<photons;i++){	//fire #photons in the scene not up	
    
    x = ((float)rand() / (RAND_MAX))*2 - 1; 
    y = -abs(((float)rand() / (RAND_MAX))*2 - 1);
    z = ((float)rand() / (RAND_MAX))*2 - 1;
    Photon p;
    p.ray.position = l; // set attributes for photon p
    p.ray.direction = Vector(x,y,z);
    p.ray.direction.normalise();
    float intensity = 1.0f/365.0f;
    p.intensity = Colour(intensity,intensity,intensity,0.0f); 
    int depthD = 670; 
    bool isCaustic = false;
    int depthS = false; 
    scene->photonTracer(p,level,depthD,depthS, isCaustic); 
    
  }

  
	typedef KDTreeSingleIndexAdaptor<	L2_Simple_Adaptor<float, PointCloud > ,	PointCloud, 3 /* dim */	> kdTree;
    
    // normal tree
    kdTree  median(3, scene->nMap, KDTreeSingleIndexAdaptorParams(15) );
	  median.buildIndex();
    scene->treeN = &median;

    //caustic tree
    kdTree cmedian(3, scene->cMap, KDTreeSingleIndexAdaptorParams(15));
    cmedian.buildIndex();
    scene->treeC = &cmedian;
  
    Vector camera = Vector(1.01,-5,-7);
    Ray ray;
    ray.position.x = 0.0001f;
    ray.position.y = 0.0f;
    ray.position.z = 0.0f;


  for (int y = 0; y < height; y += 1){
    for (int x = 0; x < width; x += 1){
      float fx = (float)x/(float)width;
      float fy = (float)y/(float)height;

      Vector direction;
      ray.direction.x = (fx-0.5f);
      ray.direction.y = (0.5f-fy);
      ray.direction.z = 0.5f;
      ray.direction.normalise();

      Colour colour;
      float depth;
      bool isCaustic = true;

      scene->raytrace(ray,Vertex(camera.x,camera.y,camera.z),colour,l,level,isCaustic);

      fb->plotPixel(x, y, colour.r, colour.g, colour.b);
      fb->plotDepth(x,y, depth);
    }
    cerr << "*" << flush;
  }
  
  // Output the framebuffer.
  fb->writeRGBFile((char *)"final_light_up.ppm");
  //  fb->writeDepthFile((char *)"depth.ppm");
  return 0;
  
}
