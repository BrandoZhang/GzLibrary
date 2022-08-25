#include	"stdafx.h"
#include	"stdio.h"
#include	"math.h"
#include	"Gz.h"
#include	"rend.h"
/*   CS580 HW   */
#include    "stdafx.h"  
#include	"Gz.h"

/* Convention for HW1 */
#define PIXEL_DEPTH 3  // 3 for RGB, 4 for RGBA

GzRender::GzRender(int xRes, int yRes)
{
/* HW1.1 create a framebuffer for MS Windows display:
 -- set display resolution
 -- allocate memory for framebuffer : 3 bytes(b, g, r) x width x height
 -- allocate memory for pixel buffer
 */
    // Boundary check on image size of xres and yres
    xres = min(max(xRes, 0), MAXXRES);  // TODO: cast?
    yres = min(max(yRes, 0), MAXYRES);  // TODO: cast?
    framebuffer = new char[PIXEL_DEPTH * xres * yres];
    pixelbuffer = new GzPixel[xres * yres];
}

GzRender::~GzRender()
{
/* HW1.2 clean up, free buffer memory */

	delete[] framebuffer;
	delete[] pixelbuffer;
}

int GzRender::GzDefault()
{
/* HW1.3 set pixel buffer to some default values - start a new frame */

    // Define default pixel value (R, G, B, A, Z)
    GzPixel defaultPixel = {0, 0, 0, 1, 0};  // Black background
    for (int i = 0; i < xres * yres; ++i) {
        // Set default value to each pixel of pixelbuffer
        pixelbuffer[i] = defaultPixel;
        // Set default value to each pixel of framebuffer (BGR order)
        framebuffer[PIXEL_DEPTH * i + 0] = (unsigned char) defaultPixel.blue;
        framebuffer[PIXEL_DEPTH * i + 1] = (unsigned char) defaultPixel.green;
        framebuffer[PIXEL_DEPTH * i + 2] = (unsigned char) defaultPixel.red;
    }
	return GZ_SUCCESS;
}


int GzRender::GzPut(int i, int j, GzIntensity r, GzIntensity g, GzIntensity b, GzIntensity a, GzDepth z)
{
/* HW1.4 write pixel values into the buffer */

    if (i >= 0 && i < xres && j >= 0 && j < yres) {  // boundary check
        int idx = ARRAY(i, j);
        GzPixel pixelVal = {r, g, b, a, z};
        pixelbuffer[idx] = pixelVal;
    }
	return GZ_SUCCESS;
}


int GzRender::GzGet(int i, int j, GzIntensity *r, GzIntensity *g, GzIntensity *b, GzIntensity *a, GzDepth *z)
{
/* HW1.5 retrieve a pixel information from the pixel buffer */

	if (i >= 0 && i < xres && j >= 0 && j < yres) {  // boundary check
		int idx = ARRAY(i, j);
		GzPixel pixelVal = pixelbuffer[idx];
		*r = pixelVal.red;
		*g = pixelVal.green;
		*b = pixelVal.blue;
		*a = pixelVal.alpha;
		*z = pixelVal.z;
	}
	return GZ_SUCCESS;
}


int GzRender::GzFlushDisplay2File(FILE* outfile)
{
/* HW1.6 write image to ppm file -- "P6 %d %d 255\r" */

    fprintf(outfile, "P6 %d %d 255\n", xres, yres);
    // Ignore off-screen coordinate
    for (int i = 0; i < xres * yres; ++i) {
        GzPixel pixelVal = pixelbuffer[i];

        // Clamp to 0-4095 (valid range of 12-bits)
        GzIntensity red = min(max(pixelVal.red, 0), 4095);
        GzIntensity green = min(max(pixelVal.green, 0), 4095);
        GzIntensity blue = min(max(pixelVal.blue, 0), 4095);

        // Convert to 8-bits
        red = red >> 4;
        green = green >> 4;
        blue = blue >> 4;
        unsigned char redVal = (unsigned char) red;
        unsigned char greenVal = (unsigned char) green;
        unsigned char blueVal = (unsigned char) blue;

        // Write ppm file in RGB order
        unsigned char color[PIXEL_DEPTH] = {redVal, greenVal, blueVal};
        fwrite(color, 1, PIXEL_DEPTH, outfile);
    }
	return GZ_SUCCESS;
}

int GzRender::GzFlushDisplay2FrameBuffer()
{
/* HW1.7 write pixels to framebuffer: 
	- put the pixels into the frame buffer
	- CAUTION: when storing the pixels into the frame buffer, the order is blue, green, and red 
	- NOT red, green, and blue !!!
*/

    // ignore off-screen coordinate
    for (int i = 0; i < xres * yres; ++i) {
        GzPixel pixelVal = pixelbuffer[i];

        // Clamp to 0-4095 (valid range of 12-bits)
        GzIntensity red = min(max(pixelVal.red, 0), 4095);
        GzIntensity green = min(max(pixelVal.green, 0), 4095);
        GzIntensity blue = min(max(pixelVal.blue, 0), 4095);

        // Convert to 8-bits
        red = red >> 4;
        green = green >> 4;
        blue = blue >> 4;
        unsigned char redVal = (unsigned char) red;
        unsigned char greenVal = (unsigned char) green;
        unsigned char blueVal = (unsigned char) blue;

        // Framebuffer in BGR order
        framebuffer[PIXEL_DEPTH * i + 0] = blueVal;
        framebuffer[PIXEL_DEPTH * i + 1] = greenVal;
        framebuffer[PIXEL_DEPTH * i + 2] = redVal;
    }
	return GZ_SUCCESS;
}