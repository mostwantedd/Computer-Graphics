# Computer-Graphics
• The aim of the project was to create a rendering program in C++. The program is able to generate a rasterised wireframe image of any model given only a .ply file, compute the depth map of any model, perform shading, illuminate any scene using photon mapping with integrated BRDF that also supports reflection and refraction effects.

• The main priority of this project was obtaining photo-realistic images (minimum 720x720p) in under 40 minutes. This requirement was achieved by optimising the path tracing algorithm through k-d trees and the Russian Roulette technique, as well as using a large enough number of photons (600k normal photons, 300k caustic photons) and the BRDF function.


![alt text](https://github.com/mostwantedd/Computer-Graphics/blob/main/farak9.jpg?raw=true)


![alt text](https://github.com/mostwantedd/Computer-Graphics/blob/main/a12_1big.jpg?raw=true)


![alt text](https://github.com/mostwantedd/Computer-Graphics/blob/main/farak7.jpg?raw=true)
