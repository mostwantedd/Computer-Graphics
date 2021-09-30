#pragma once
#include "ray.h"
#include "colour.h"

class Photon {
public:
    Ray ray;
    Colour intensity;
};