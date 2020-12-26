# Computer Graphics 
<a href="https://www.dis.uniroma1.it/"><img src="http://www.dis.uniroma1.it/sites/default/files/marchio%20logo%20eng%20jpg.jpg" width="500"></a>

AIRO course of Fundamentals for Computer Graphics:

Three homeworks (Raytracing, to do, to do) and final project (to do).

Tools we used: C++, Visual Studio Code, CMake, gcc and ninja.

## Raytracing
In 3D computer graphics, ray tracing is a rendering technique for generating an image by tracing the path of light as pixels in an image plane and simulating the effects of its encounters with virtual objects. The technique is capable of producing a high degree of visual realism, more so than typical scanline rendering methods, but at a greater computational cost. This makes ray tracing best suited for applications where taking a relatively long time to render can be tolerated, such as in still computer-generated images, and film and television visual effects (VFX), but generally more poorly suited to real-time applications such as video games, where speed is critical in rendering each frame. Ray tracing is capable of simulating a variety of optical effects, such as reflection and refraction, scattering, and dispersion phenomena (such as chromatic aberration).


First part: implement images with own colors, textures... as we can see here:
![](Raytrace/out/lowres/03_texture_720_256.jpg)

Second part: implementn of raytracing shader, with refraction for non-thin materials.
![](Raytrace/out/Refraction/glass_(notThin).jpg)

Third part: implement previous images with a cartoon shader.
![](Raytrace/out/shade_cartoon/materialsb.png)
.



## Pathtracing
First part: implement naive pathtracing and refraction

![](PathTrace/out/naive/02_matte_720_256.jpg)
![](PathTrace/out/Refraction/naive.jpg)



Second part: implement the pathtrace with lights

![](PathTrace/out/path/01_cornellbox_512_256.jpg)
![](PathTrace/out/path/15_classroom_720_256.jpg)

Third part: creation of personal scenes to render 

![](PathTrace/out/OwnScenes/1own.png) 
![](PathTrace/out/OwnScenes/2own.png)
![](PathTrace/out/OwnScenes/3own.png)


Fourt part: test the pathtrace on large scenes to see the effects of the lights

![](PathTrace/out/LargeScenesHD/bisEXT_hd.jpg)




## Volumetric pathtracing 

## Hairshading
