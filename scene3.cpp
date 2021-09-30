/***************************************************************************
 *
 * krt - Kens Raytracer - Coursework Edition. (C) Copyright 1997-2018.
 *
 * Do what you like with this code as long as you retain this comment.
 */
// #define _USE_MATH_DEFINES
#include "vector.h"
#include "vertex.h"
#include "scene3.h"
#include "photon.h"
#include <math.h>

#include "nanoflann-master/include/nanoflann.hpp" 

using namespace std;

//methods reflection, transparent and fresnel are not written by me
//ray.direction made just a few modifications to match my requirements
// they are originating from https://www.scratchapixel.com/lessons/3d-basic-rendering/introduction-to-shading/reflection-refraction-fresnel 

Scene3::Scene3(){
	objects = 0;
}

void Scene3::trace(Ray ray, Object *objects, Hit &hit){
	Hit current_hit;
	hit.flag = false;
	hit.t = 0.0f;
	hit.what = 0;
	while (objects != 0){
		Hit hit_current;
		objects->intersection(ray, hit_current);
		if (hit_current.flag == true){
			if (hit.flag == false){
				hit = hit_current;
			} else if (hit_current.t < hit.t){
				hit = hit_current;
			}
		}
		objects = objects->next;
	}
}

Vector Scene3::reflection(Vector direction, Vector normal){
  Vector rDirection;
  rDirection = direction - 2 * normal.dot(direction)*normal;//astea erau invers pe
  rDirection.normalise();
  return rDirection;
}

double Scene3::clamp(double x, double upper, double lower){
    return min(upper, max(x, lower));
}

Vector Scene3::transparent(Vector direction, Vector normal, float ior){
  direction.normalise();
  normal.normalise();
  float cosi = clamp(-1.0, 1.0, normal.dot(direction)); 
  float etai = 1, etat = ior, eta; 
  if (cosi < 0) { 
	  cosi = -cosi; 
  } else { 
	std::swap(etai, etat); normal.negate(); 
  } 
  eta = etai / etat; 
  float k = 1 - eta * eta * (1 - cosi * cosi); 
  Vector nullVect;
  nullVect.x = 0.0; nullVect.y = 0.0; nullVect.z = 0.0;
  Vector tDirection = (eta * direction + (eta * cosi - sqrtf(k)) * normal);
  tDirection.normalise();
  return k < 0 ? nullVect : tDirection; 
}

float Scene3::fresnel(Vector &direction, Vector &normal, float &ior,float &kr){
    // float kr;
    float cosi = clamp(-1, 1, normal.dot(direction));//dotProduct(direction, normal)); 
    float etai = 1, etat = ior; 
    if (cosi > 0) { 
		std::swap(etai, etat); 
	} 
    // Compute sini using Snell's law
    float sint = etai / etat * sqrtf(std::max(0.f, 1 - cosi * cosi)); 
    // Total internal reflection
    if (sint >= 1) { 
        kr = 1; 
    } 
    else { 
        float cost = sqrtf(std::max(0.f, 1 - sint * sint)); 
        cosi = fabsf(cosi); 
        float Rs = ((etat * cosi) - (etai * cost)) / ((etat * cosi) + (etai * cost)); 
        float Rp = ((etai * cosi) - (etat * cost)) / ((etai * cosi) + (etat * cost)); 
        kr = (Rs * Rs + Rp * Rp) / 2; 
    } 
    return kr;
} 

void Scene3::photonTracer(Photon p, int level, int depthD, bool depthS, bool isCaustic){

    Ray ray = p.ray;
	Hit hit;
	float small_e = 0.0001f;
    level = level - 1;
    if(level < 0){return; }
	
	if(depthD <= 0 && isCaustic) {return;} //daca nu ai asta apar puncte albe pe pereti/
	
    // Need to trace ray through scene
    // First need to find first primative
    trace(ray,objects,hit); 
	if(hit.flag == false){return;}

	p.ray.position = hit.position;// photon position
    
    Colour diffuse = hit.what->material->d;
    Colour specular = hit.what->material->s;
	Colour ambient = hit.what->material->a;

	// Russian Roulette
    float Pr = max({diffuse.r + specular.r, diffuse.g + specular.g,diffuse.b + specular.b});
    float D = diffuse.r + diffuse.g + diffuse.b;
    float S = specular.r + specular.g + specular.b;
    float Pd = (D/(D + S))*Pr;
    float Ps = Pr - Pd;    
    float random = (float) (rand()) / (float) (RAND_MAX);

	
    if(random>=0 && random<=Pd) {
		if(isCaustic && depthS == false){ return;}///////////////////////// need to have at least one specular bounce.///daca nu ai asta apar puncte albe pe pereti
		depthD = depthD - 1; // diffuse reflection////////////////////////////////////////////////////
        
		// photon becomes node in kd tree
		PointCloud::Point node;
		node.x = p.ray.position.x;
		node.y = p.ray.position.y;
		node.z = p.ray.position.z;
		node.incomingray = p.ray.direction;
		node.intensity = p.intensity;
		if(isCaustic){ cMap.nodes.push_back(node);}else{ nMap.nodes.push_back(node);}

        Colour rintensity; // reflected intensity
        Photon rp;  // reflected photon
		Vector light = ray.direction; light.negate();
        rintensity.r = p.intensity.r*diffuse.r*max(light.dot(hit.normal),0.0f);
        rintensity.g = p.intensity.g*diffuse.g*max(light.dot(hit.normal),0.0f);
        rintensity.b = p.intensity.b*diffuse.b*max(light.dot(hit.normal),0.0f);
		

        Vector newDirection;//photon bounces in this random direction:
		float x,y,z;
		x = ((float)rand()/(RAND_MAX))*2 - 1;
		y = abs(((float)rand()/(RAND_MAX))*2 - 1); // avoids going IN object
		z = ((float)rand()/(RAND_MAX))*2 - 1; 
		newDirection.x = x; newDirection.y = y; newDirection.z = z;
        rp.ray.direction = newDirection;
		rp.ray.direction.normalise();

        rp.ray.position.x = hit.position.x + small_e*rp.ray.direction.x;
        rp.ray.position.y = hit.position.y + small_e*rp.ray.direction.y;
        rp.ray.position.z = hit.position.z + small_e*rp.ray.direction.z;
        rp.intensity = rintensity; 
        
        photonTracer(rp,level,depthD,depthS,isCaustic);


    }else if(random>Pd && random <= Ps + Pd) {
		depthS = true;
		

		float kr, kt;
		kr = fresnel(ray.direction, hit.normal, hit.what->material->ior, kr);
		kt = 1-kr;
		float e = (float) (rand()) / (float) (RAND_MAX);
		
		Colour i; //intensity
		i.r = (p.intensity.r*specular.r)/Ps;
		i.g = (p.intensity.g*specular.g)/Ps;
		i.b = (p.intensity.b*specular.b)/Ps;

		
		//if (kr<1){////////////////////////////////////////////////////////
		if(e<kr && hit.what->material->reflect){ // reflect

			Photon rp;//reflected photon
			rp.ray.direction = reflection(ray.direction,hit.normal);
			rp.ray.position.x = hit.position.x + small_e*rp.ray.direction.x;
			rp.ray.position.y = hit.position.y + small_e*rp.ray.direction.y;
			rp.ray.position.z = hit.position.z + small_e*rp.ray.direction.z;
			rp.intensity = p.intensity;
			photonTracer(rp,level,depthD,depthS,isCaustic);
	
		}
		// else{// refract
		if (e>kr && hit.what->material->refract){

			Photon tp;//transarent/refracted photon
			tp.intensity = p.intensity;
			tp.ray.direction = transparent(ray.direction,hit.normal,hit.what->material->ior);
			tp.ray.position.x = hit.position.x + small_e*tp.ray.direction.x;
			tp.ray.position.y = hit.position.y + small_e*tp.ray.direction.y;
			tp.ray.position.z = hit.position.z + small_e*tp.ray.direction.z;
			photonTracer(tp,level,depthD,depthS,isCaustic);

		}

    } else {
		
        // light doesent escape -> no bounce
		PointCloud::Point node;
		node.x = p.ray.position.x; node.y = p.ray.position.y; node.z = p.ray.position.z;
		node.incomingray = p.ray.direction;
		node.intensity = p.intensity;
		if(isCaustic){
			cMap.nodes.push_back(node);	
		}else{
			nMap.nodes.push_back(node);
		}
		
    }
    Photon shadow;
    Hit shadow_hit;
    shadow.ray.position = hit.position;
    shadow.ray.direction = ray.direction;
    shadow.intensity.r = 0.0f; 
	shadow.intensity.b = 0.0f; 
	shadow.intensity.g = 0.0f;
	trace(shadow.ray,objects,shadow_hit);
    if(shadow_hit.flag == true){ /////////////////////////////////////////////////////////// if our shadow ray intersects we store the photon at the hit position. No shadow rays for refractive	
		shadow.ray.position = shadow_hit.position;
		PointCloud::Point n;
		n.x = shadow.ray.position.x;
		n.y = shadow.ray.position.y;
		n.z = shadow.ray.position.z;
		n.incomingray = shadow.ray.direction;
		n.intensity = shadow.intensity;
		if(isCaustic){
			cMap.nodes.push_back(n);	
		}else{
			nMap.nodes.push_back(n);
		}
	}
	
}





void Scene3::raytrace(Ray ray, Vertex camera, Colour &colour, Vertex l, int level, bool isCaustic){
	
	level = level - 1; // recursed too far
	if(level<0){return;}
	float small_e = 0.0001f; //multiplier
	Vector viewer,reflected,incom; // vectors for brdf

	Hit hit;
	trace(ray,objects,hit);
	if(hit.flag==false){return;}

	Ray tolight; // compute ray backwords, from camera to light
	tolight.direction.x = l.x - hit.position.x;
	tolight.direction.y = l.y - hit.position.y;
	tolight.direction.z = l.z - hit.position.z;
	
	tolight.direction.normalise();
	tolight.position.x = hit.position.x + 0.0001*tolight.direction.x;
	tolight.position.y = hit.position.y + 0.0001*tolight.direction.y;
	tolight.position.z = hit.position.z + 0.0001*tolight.direction.z;
	
	float kr,kt; 
	kr = fresnel(ray.direction, hit.normal, hit.what->material->ior, kr);
	kt = 1-kr;	
	bool bool_refract = true;
	

	viewer.x = - hit.position.x;
	viewer.y = - hit.position.y;
	viewer.z = - hit.position.z;
	viewer.normalise();
	

	Colour diffuse = hit.what->material->d;
	Colour specular = hit.what->material->s;
	Colour ambient = hit.what->material->a;
	
	Colour direct;
	direct.r = 0.0f;
	direct.g = 0.0f;
	direct.b = 0.0f;

	Colour room;
	room.r = 0.0f;
	room.g = 0.0f;
	room.b = 0.0f;

	Colour indirect;
	indirect.r = 0.0f;
	indirect.g = 0.0f;
	indirect.b = 0.0f;

	Colour caustic;
	caustic.r = 0.0f;
	caustic.g = 0.0f;
	caustic.b = 0.0f;

	Hit first;
	trace(tolight,objects,first);
	if(first.t > tolight.direction.length()){
		first.flag = false;
	}

	if(first.flag == false){ //compute direct colour
	
		tolight.direction.negate();
		hit.what->material->compute_light_colour(viewer,hit.normal,tolight.direction,direct);
		direct.r = min(direct.r,1.0f);
		direct.g = min(direct.g,1.0f);
		direct.b = min(direct.b,1.0f);

	}
	

	if(kr!=0 && hit.what->material->reflect){ //reflection
		Colour rcolour;
		Ray rray;
		rray.direction = reflection(ray.direction,hit.normal); 
		rray.position.x = hit.position.x + small_e*rray.direction.x;
		rray.position.y = hit.position.y + small_e*rray.direction.y;
		rray.position.z = hit.position.z + small_e*rray.direction.z;
		raytrace(rray,camera,rcolour,l,level,isCaustic);
		rcolour.r = kr*rcolour.r;
		rcolour.g = kr*rcolour.g;
		rcolour.b = kr*rcolour.b;
		room.add(rcolour); 
	}

	if(kt!=0 && hit.what->material->refract){ //transparent
	
		Colour tcolour;
		Ray tray;
		tray.direction = transparent(ray.direction,hit.normal,hit.what->material->ior);

		tray.position.x = hit.position.x + small_e*tray.direction.x;
		tray.position.y = hit.position.y + small_e*tray.direction.y;
		tray.position.z = hit.position.z + small_e*tray.direction.z; 

		raytrace(tray,camera,tcolour,l,level,isCaustic);

		tcolour.r = kt*tcolour.r;
		tcolour.g = kt*tcolour.g;
		tcolour.b = kt*tcolour.b;

		room.add(tcolour);
	}

	room.r = min(room.r,1.0f);
	room.g = min(room.g,1.0f);
	room.b = min(room.b,1.0f);
	
	const float intersection[3] = {hit.position.x,hit.position.y,hit.position.z};
	size_t k = 1000; 
	std::vector<size_t> boundary(k);
	std::vector<float> radius(k);
	
	k = treeN->knnSearch(&intersection[0],k,&boundary[0],&radius[0]); // perform knn search
	float r = *max_element(radius.begin(),radius.end());
	float area  = M_PI*pow(r,2);
	
	for(size_t i=0;i<k;i++) {	
		size_t median = boundary[i];
		float x = nMap.nodes[median].x;
		float y = nMap.nodes[median].y;
		float z = nMap.nodes[median].z;
		Vector incomingray = nMap.nodes[median].incomingray;
		Colour intensity = nMap.nodes[median].intensity;
		int photon_map_size = nMap.nodes.size();

		// Bidirectional reflectance distribution function
		
		incomingray.negate();
		incomingray.normalise();
		
		viewer.x = camera.x - x;
		viewer.y = camera.y - y;
		viewer.z = camera.z - z;
		viewer.normalise();
		Vector reflected = reflection(incomingray, hit.normal);


		indirect.r += min(0.5f*(diffuse.r*(max(incomingray.dot(hit.normal),0.0f))*intensity.r + 
						specular.r*pow(max(reflected.dot(viewer),0.0f),20.0f)*intensity.r + 
						ambient.r*intensity.r),1.0f);
		indirect.g += min(0.5f*(diffuse.g*(max(incomingray.dot(hit.normal),0.0f))*intensity.g + 
						specular.g*pow(max(reflected.dot(viewer),0.0f),20.0f)*intensity.g + 
						ambient.r*intensity.r),1.0f);
		indirect.b += min(0.5f*(diffuse.b*(max(incomingray.dot(hit.normal),0.0f))*intensity.b + 
						specular.b*pow(max(reflected.dot(viewer),0.0f),20.0f)*intensity.b + 
						ambient.r*intensity.r),1.0f);


	}

	indirect.r = min(indirect.r/area,1.0f);
	indirect.g = min(indirect.g/area,1.0f);
	indirect.b = min(indirect.b/area,1.0f);
	
	if(isCaustic){ 
		
		size_t caustick = 35;
		std::vector<size_t> boundaryC(caustick);
		std::vector<float> radiusC(caustick);

		caustick = treeC->knnSearch(&intersection[0],caustick,&boundaryC[0],&radiusC[0]);
		float rC = *max_element(radiusC.begin(),radiusC.end());
		float areaC  = M_PI*pow(rC,2);
		float k = 1;
		for(size_t i = 0; i<caustick;i++){
			
			size_t cmedian = boundaryC[i];
			Colour intensityC = cMap.nodes[cmedian].intensity;
			
			float distance = sqrt(pow(abs(cMap.nodes[cmedian].x - hit.position.x),2) + pow(abs(cMap.nodes[cmedian].y - hit.position.y),2) + pow(abs(cMap.nodes[cmedian].z - hit.position.z),2));
			float direction = 1 - distance/(sqrt(rC)); 
			incom = cMap.nodes[cmedian].incomingray;
			incom.negate();
			incom.normalise();

			viewer.x = - cMap.nodes[cmedian].x;
			viewer.y = - cMap.nodes[cmedian].y;
			viewer.z = - cMap.nodes[cmedian].z;
			viewer.normalise();
			reflected = reflection(cMap.nodes[cmedian].incomingray, hit.normal);
			reflected.normalise();

			caustic.r += min(0.5f*(direction*diffuse.r*(max(incom.dot(hit.normal),0.0f))*intensityC.r + specular.r*pow(max(reflected.dot(viewer),0.0f),20.0f)*intensityC.r),1.0f); 
			caustic.g += min(0.5f*(direction*diffuse.g*(max(incom.dot(hit.normal),0.0f))*intensityC.g + specular.g*pow(max(reflected.dot(viewer),0.0f),20.0f)*intensityC.g),1.0f);
			caustic.b += min(0.5f*(direction*diffuse.b*(max(incom.dot(hit.normal),0.0f))*intensityC.b + specular.b*pow(max(reflected.dot(viewer),0.0f),20.0f)*intensityC.b),1.0f);

		}

		caustic.r = min(caustic.r/areaC,1.0f);
		caustic.g = min(caustic.g/areaC,1.0f);
		caustic.b = min(caustic.b/areaC,1.0f);

	}

	colour.add(direct);
	colour.add(indirect);
	colour.add(room);
	colour.add(caustic);

}

