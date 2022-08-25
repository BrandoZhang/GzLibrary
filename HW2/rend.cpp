#include	"stdafx.h"
#include	"stdio.h"
#include	"math.h"
#include	"Gz.h"
#include	"rend.h"
#include    <algorithm>

/***********************************************/
/* HW1 methods: copy here the methods from HW1 */
/* Conventions */
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
    GzPixel defaultPixel = {0, 0, 0, 1, INT_MAX};  // Black background
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
        if (z < pixelbuffer[idx].z) {  // z-buffer checking
            GzPixel pixelVal = {r, g, b, a, z};
            pixelbuffer[idx] = pixelVal;
        }
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


/***********************************************/
/* HW2 methods: implement from here */

int GzRender::GzPutAttribute(int numAttributes, GzToken	*nameList, GzPointer *valueList) 
{
/* HW 2.1
-- Set renderer attribute states (e.g.: GZ_RGB_COLOR default color)
-- In later homeworks set shaders, interpolaters, texture maps, and lights
*/

    // Loop attributes
    for (int i = 0; i < numAttributes; ++i) {
        // Check the token name (GZ_RGB_COLOR in HW2)
        if (nameList[i] == GZ_RGB_COLOR) {
            // Cast the corresponding value (GzColor* in HW2)
            GzColor* colorVal = (GzColor*) valueList[i];
            // Sent to destination (flatcolor in HW2)
            flatcolor[RED] = (*colorVal)[RED];
            flatcolor[GREEN] = (*colorVal)[GREEN];
            flatcolor[BLUE] = (*colorVal)[BLUE];
        }
    }
	return GZ_SUCCESS;
}

float lineEquation(float A, float B, float C, float i, float j) {
    /* Test the relationship between `(i, j)` and a 2D reference line of  `Ax + By + C = 0`.
     * A, B, and C are coefficients of the reference line.
     * if returns 0: `(i, j)` is on the line.
     * if returns <0: `(i, j)` is in right/below half-plane of the line.
     * if returns >0: `(i, j)` is in left/above half-plane of the line.
     */
    return A * i + B * j + C;
}

int GzRender::GzPutTriangle(int	numParts, GzToken *nameList, GzPointer *valueList) 
/* numParts - how many names and values */
{
/* HW 2.2
-- Pass in a triangle description with tokens and values corresponding to
      GZ_NULL_TOKEN:		do nothing - no values
      GZ_POSITION:		3 vert positions in model space
-- Invoke the rastrizer/scanline framework
-- Return error code
*/

    // Approach 1: LEE
    for (int p = 0; p < numParts; ++p) {
        if (nameList[p] == GZ_POSITION) {
            // Rasterization
            GzCoord *vertices = (GzCoord *) valueList[0];

            // Sort vertices to be clockwise (CW), i.e., 1-2, 2-3, 3-1,
            //  where V1 is supposed to be the top vertex and V3 the bottom vertex.
            // Sort by Y-coords.
            if (vertices[0][1] > vertices[1][1]) {  // V1 is below V2
                for (int i = 0; i < 3; ++i) {
                    std::swap(vertices[0][i], vertices[1][i]);  // Swap V1 and V2
                }
            }  // Now, V1 should be above or horizontally equal to V2
            if (vertices[0][1] > vertices[2][1]) {  // V1 is below V3
                for (int i = 0; i < 3; ++i) {
                    std::swap(vertices[0][i], vertices[2][i]);  // Swap V1 and V3
                }
            }  // Now, V1 should be above both V2 and V3, or horizontally equal to one (or two) of them
            if (vertices[1][1] > vertices[2][1]) {  // V2 is below V3
                for (int i = 0; i < 3; ++i) {
                    std::swap(vertices[1][i], vertices[2][i]);  // Swap V2 and V3
                }
            }  // Now, V2 should be above or horizontally equal to V3
            // Determine L/R relationship.
            // The edge that most close to vertical is 3-1, which is treated as a reference edge.
            float dYEdge31 = vertices[0][1] - vertices[2][1];  // dY of Edge31
            float dXEdge31 = vertices[0][0] - vertices[2][0];  // dX of Edge31
            float CEdge31 = dXEdge31 * vertices[2][1] - dYEdge31 * vertices[2][0];  // coefficient C of Edge31: C = dX * Y - dY * X
            // Plug V2 into the line equation of Edge31
            float halfPlaneV2 = lineEquation(dYEdge31, -1.0f * dXEdge31, CEdge31, vertices[1][0], vertices[1][1]);
            if (halfPlaneV2 > 0) {  // V2 is in left half-plane of edge31.
                for (int i = 0; i < 3; ++i) {
                    std::swap(vertices[1][i], vertices[2][i]);  // Swap V2 and V3
                }
            }
            // TODO: Assert that halfPlaneV2 = 0 if dYEdge31 = 0
            // Now, the vertices are sorted clockwise as 1-2, 2-3, 3-1.

            // Compute A, B, C (the line equation coefficients) for each edge
            //  General 2D line equation: Ax + By + C = 0
            //  A = dY, B = -dX, C = dX * Y - dY * X
            //  where (X, Y) is the tail of that edge
            //  Edge 1-2, Vec = (dX, dY, dZ)
            float dXEdge12 = vertices[1][0] - vertices[0][0];  // dX of Edge12
            float dYEdge12 = vertices[1][1] - vertices[0][1];  // dY of Edge12
            float dZEdge12 = vertices[1][2] - vertices[0][2];  // dZ of Edge12
            float CEdge12 = dXEdge12 * vertices[0][1] - dYEdge12 * vertices[0][0];  // coefficient C of Edge12
            //  Edge 2-3, Vec = (dX, dY, dZ)
            float dXEdge23 = vertices[2][0] - vertices[1][0];  // dX of Edge23
            float dYEdge23 = vertices[2][1] - vertices[1][1];  // dY of Edge23
            float dZEdge23 = vertices[2][2] - vertices[1][2];  // dZ of Edge23
            float CEdge23 = dXEdge23 * vertices[1][1] - dYEdge23 * vertices[1][0];  // coefficient C of Edge23
            //  Edge 3-1
            dXEdge31 = vertices[0][0] - vertices[2][0];  // dX of Edge31
            dYEdge31 = vertices[0][1] - vertices[2][1];  // dY of Edge31
            CEdge31 = dXEdge31 * vertices[2][1] - dYEdge31 * vertices[2][0];  // coefficient C of Edge31

            // Compute plane coefficients
            //  General 3D plane equation: Ax + By + Cz + D = 0
            //  vector (A, B, C) is normal to the plane,
            //  and hence can be calculated by cross product of two edges.
            //  Reference: https://en.wikipedia.org/wiki/Cross_product
            float APlane = dYEdge12 * dZEdge23 - dZEdge12 * dYEdge23;
            float BPlane = dZEdge12 * dXEdge23 - dXEdge12 * dZEdge23;
            float CPlane = dXEdge12 * dYEdge23 - dYEdge12 * dXEdge23;
            float DPlane = -1.0f * (APlane * vertices[2][0] + BPlane * vertices[2][1] + CPlane * vertices[2][2]);

            // Find bounding box of the triangle
            float minX = min(min(vertices[0][0], vertices[1][0]), vertices[2][0]);
            float maxX = max(max(vertices[0][0], vertices[1][0]), vertices[2][0]);
            float minY = min(min(vertices[0][1], vertices[1][1]), vertices[2][1]);
            float maxY = max(max(vertices[0][1], vertices[1][1]), vertices[2][1]);
            // Round the bounding box
            int minXPix = (int) ceil(minX);
            int maxXPix = (int) ceil(maxX);
            int minYPix = (int) ceil(minY);
            int maxYPix = (int) ceil(maxY);

            // Iterate inside bounding box for rasterization
            for (int i = minXPix; i <= maxXPix; ++i) {
                for (int j = minYPix; j <= maxYPix; ++j) {

                    float sign12 = lineEquation(dYEdge12, -1.0f * dXEdge12, CEdge12, i, j);
                    float sign23 = lineEquation(dYEdge23, -1.0f * dXEdge23, CEdge23, i, j);
                    float sign31 = lineEquation(dYEdge31, -1.0f * dXEdge31, CEdge31, i, j);

                    if ((sign12 > 0 && sign23 > 0 && sign31 > 0 && CPlane != 0)  // inside and not perpendicular to view point
                    || (sign12 < 0 && sign23 < 0 && sign31 <0 && CPlane != 0)
                        || sign12 == 0 || sign23 == 0 || sign31 == 0) {  // edges
                        // Z value interpolation
                        float tempZ = -1.0f * (APlane * (float) i + BPlane * (float) j + DPlane) / CPlane;
                        int z = (int) ceil(tempZ);
                        // Color the pixel
                        GzIntensity red = ctoi(flatcolor[0]);
                        GzIntensity green = ctoi(flatcolor[1]);
                        GzIntensity blue = ctoi(flatcolor[2]);
                        GzPut(i, j, red, green, blue, 1, z);
                    }
                }
            }

        }
    }
	return GZ_SUCCESS;
}

