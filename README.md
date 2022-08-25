# GzLibrary: A C++ 3D Graphics Render Library

## Features

An OpenGL-like lightweight library that supports the complete model-to-image rendering pipeline.

Given an `.asc` file which defines triangle meshes of a 3D model, this library produces a 2D image of that model in `.ppm` format and display on the screen by flushing the image into framebuffer as well.

- **Rasterization**: Performed Linear Expression Evaluation (LEE) with vertex sorting for rasterization; solved hidden surface problem with Z-buffer; implemented Backface Culling algorithm to remove invisible faces thereby speedup the process.
- **Transformation**: Achieved translation, rotation, and scale by leveraging homogeneous coordinate; supported transformation between different spaces (_i.e._, Model Space, World Space, Image (Camera) Space, Perspective (NDC) Space, Screen Space).

   |                      Image w/ translation                      |                 Image w/ rotation                 |                 Image w/ scaling                 |
   :----------------------------------------------------:|:--------------------------------------------------------:|:------------------------------------------------:|
   ![image with translation](/docs/translation_sample1.gif) | ![image with rotation](/docs/rotation_sample2.gif) | ![image with scaling](/docs/scaling_sample1.gif) |


- **Shading**: Implemented Phong shading, Gouraud shading, and Flat shading; supported changing material appearance by altering specular/diffuse/ambient terms in shading equation.

   |          Image w/ Phong shading                |                 Image w/ Gouraud shading                 |               Image w/ Flat shading                |
   :-----------------------------------------------:|:--------------------------------------------------------:|:--------------------------------------------------:|
   ![image with phong shading](/docs/phong_shading.png) | ![image with gouraud shading](/docs/gouraud_shading.png) | ![image with flat shading](/docs/flat_shading.png) | 

   | Image w/ metal appearance |                 Image w/ plastic appearance                  | Image w/ matte appearance |
   :--------------------------:|:------------------------------------------------------------:|:-------------------------:|
   ![image with metal appearance](/docs/metal_material.png) | ![image with plastic appearance](/docs/plastic_material.png) | ![image with matte appearance](/docs/matte_material.png)


- **Texture**: Provided functionalities on applying image texture or procedure texture; performed perspective correction in affine space to correct distortion.

  Image w/ perspective correction             |  Image w/o perspective correction
  :-------------------------:|:-------------------------:
  ![image with perspective correction](/docs/w_perspective_correction.png) | ![image without perspective correction](/docs/wo_perspective_correction.png)
  **Image w/ image texture** | **Image w/ procedure texture**
  ![image with image texture](/docs/image_texture.png) | ![image with procedure texture](/docs/procedure_texture.png)

- **Anti-aliasing**: Used regular-jittered super-sampling for anti-aliasing.

## Quick Start

This project requires Windows platform and Visual Studio 2019 to build and run.

### How to build

1. Download and install the toolkit for Visual Studio 2019 by checking `Desktop development with C++` (also `C++ MFC for latest v142 build tools` in the installation details) in `Tools -> Get tools and features`.
2. Open the solution file, click `Build`.

### How to run

To visualize `.ppm` file, you can download the free software [Irfanview](https://www.irfanview.com/).

1. Change or put your `.asc` file and texture file in `src/`.
2. Open solution file, in Visual Studio, click `Debug -> Start Debugging` to launch the program.
3. Click `Renderer -> RunRender` to visualize image.
4. To perform transformation, click `Edit` and select corresponding transformations, then click `Renderer -> RunRender` to view the result.

## Acknowledge

The skeleton code (_i.e._, API and window) comes from course [CSCI 580 3-D Graphics and Rendering, taught by Ulrich Neumann](https://www.cs.usc.edu/directory/faculty/profile/?lname=Neumann&fname=Ulrich).


