/* CS580 Homework 3 */

#include	"stdafx.h"
#include	"stdio.h"
#include	"math.h"
#include	"Gz.h"
#include	"rend.h"
#include    <vector>
#include    <algorithm>

#define PI (float) 3.14159265358979323846
/* Conventions */
#define PIXEL_DEPTH 3  // 3 for RGB, 4 for RGBA
using std::vector;
// TODO: Implement most inner logics with vector, array is nightmare.

int GzRender::GzRotXMat(float degree, GzMatrix mat)
{
/* HW 3.1
// Create rotate matrix : rotate along x axis
// Pass back the matrix using mat value
*/
//    mat = {
//            1.0,    0.0,            0.0,            0.0,
//            0.0,    cos_radian,     -sin_radian,    0.0,
//            0.0,    sin_radian,     cos_radian,     0.0,
//            0.0,    0.0,            0.0,            1.0
//    };
    // Initialize with all 0 homogenous matrix
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            mat[i][j] = 0;
        }
    }

    float radian = degree * PI / 180.0f;
    float cos_radian = cos(radian);
    float sin_radian = sin(radian);

    // Set to RotateX matrix
    mat[0][0] = 1.0;
    mat[1][1] = cos_radian;
    mat[1][2] = - sin_radian;
    mat[2][1] = sin_radian;
    mat[2][2] = cos_radian;
    mat[3][3] = 1.0;
	return GZ_SUCCESS;
}

int GzRender::GzRotYMat(float degree, GzMatrix mat)
{
/* HW 3.2
// Create rotate matrix : rotate along y axis
// Pass back the matrix using mat value
*/
    // Initialize with all 0 homogenous matrix
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            mat[i][j] = 0;
        }
    }

    float radian = degree * PI / 180.0f;
    float cos_radian = cos(radian);
    float sin_radian = sin(radian);

    // Set to RotateY matrix
    mat[0][0] = cos_radian;
    mat[0][2] = sin_radian;
    mat[1][1] = 1.0;
    mat[2][0] = - sin_radian;
    mat[2][2] = cos_radian;
    mat[3][3] = 1.0;
	return GZ_SUCCESS;
}

int GzRender::GzRotZMat(float degree, GzMatrix mat)
{
/* HW 3.3
// Create rotate matrix : rotate along z axis
// Pass back the matrix using mat value
*/
    // Initialize with all 0 homogenous matrix
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            mat[i][j] = 0;
        }
    }

    float radian = degree * PI / 180.0f;
    float cos_radian = cos(radian);
    float sin_radian = sin(radian);

    // Set to RotateZ matrix
    mat[0][0] = cos_radian;
    mat[0][1] = - sin_radian;
    mat[1][0] = sin_radian;
    mat[1][1] = cos_radian;
    mat[2][2] = 1.0;
    mat[3][3] = 1.0;
	return GZ_SUCCESS;
}

int GzRender::GzTrxMat(GzCoord translate, GzMatrix mat)
{
/* HW 3.4
// Create translation matrix
// Pass back the matrix using mat value
*/
    // Initialize with all 0 homogenous matrix
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            mat[i][j] = 0;
        }
    }

    // Set to Translation matrix
    mat[0][0] = 1.0;
    mat[1][1] = 1.0;
    mat[2][2] = 1.0;
    mat[3][3] = 1.0;
    mat[0][3] = translate[0];  // Tx
    mat[1][3] = translate[1];  // Ty
    mat[2][3] = translate[2];  // Tz
	return GZ_SUCCESS;
}


int GzRender::GzScaleMat(GzCoord scale, GzMatrix mat)
{
/* HW 3.5
// Create scaling matrix
// Pass back the matrix using mat value
*/
    // Initialize with all 0 homogenous matrix
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            mat[i][j] = 0;
        }
    }

    // Set to Scale matrix
    mat[0][0] = scale[0];  // Sx
    mat[1][1] = scale[1];  // Sy
    mat[2][2] = scale[2];  // Sz
    mat[3][3] = 1.0;
    return GZ_SUCCESS;
}


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
//	framebuffer = (char*) malloc (3 * sizeof(char) * xRes * yRes);
    framebuffer = (char*) malloc (3 * sizeof(char) * xres * yres);
    pixelbuffer = new GzPixel[xres * yres];

/* HW 3.6
- setup Xsp and anything only done once 
- init default camera 
*/
    // Set up Xsp
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            Xsp[i][j] = 0;
        }
    }
    Xsp[0][0] = (float) xres / 2.0f;
    Xsp[0][3] = (float) xres / 2.0f;
    Xsp[1][1] = -1.0f * (float) yres / 2.0f;
    Xsp[1][3] = (float) yres / 2.0f;
    Xsp[2][2] = (float) INT_MAX;
    Xsp[3][3] = 1.0f;
    // Initialize default camera
    m_camera.position[X] = DEFAULT_IM_X;
    m_camera.position[Y] = DEFAULT_IM_Y;
    m_camera.position[Z] = DEFAULT_IM_Z;
    m_camera.FOV = DEFAULT_FOV;
    for (int i = 0; i < 3; i++) {
        m_camera.lookat[i] = 0;  // default look-at point = (0, 0, 0)
        m_camera.worldup[i] = 0;
    }
    m_camera.worldup[1] = 1.0f;  // default world-up = (0, 1, 0)
    matlevel = -1;  // index of Ximage (Xsm)
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

int GzRender::GzBeginRender()
{
/* HW 3.7 
- setup for start of each frame - init frame buffer color,alpha,z
- compute Xiw and projection xform Xpi from camera definition 
- init Ximage - put Xsp at base of stack, push on Xpi and Xiw 
- now stack contains Xsw and app can push model Xforms when needed 
*/ 
    // Compute Xpi
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            m_camera.Xpi[i][j] = 0;
        }
    }
    m_camera.Xpi[0][0] = 1.0f;
    m_camera.Xpi[1][1] = 1.0f;
    m_camera.Xpi[2][2] = tan(m_camera.FOV * PI / 360.0f);  // (FOV * PI / 180) / 2
    m_camera.Xpi[3][2] = tan(m_camera.FOV * PI / 360.0f);
    m_camera.Xpi[3][3] = 1.0f;
    // TODO: Compute Xiw
    // Construct vector cl (align with camera Z axis)
    GzCoord cl;
    for (int i = 0; i < 3; ++i) {
        cl[i] = m_camera.lookat[i] - m_camera.position[i];
    }
    // Calculate vector Z by normalizing cl
    // TODO: write function for normalizing vectors
    GzCoord cameraZ;
    float clNorm = sqrt(cl[0] * cl[0] + cl[1] * cl[1] + cl[2] * cl[2]);
    for (int i = 0; i < 3; ++i) {
        cameraZ[i] = cl[i] / clNorm;
    }
    // TODO: write function for dot product and cross product
    // Calculate vector up' = up - (up * Z) Z
    GzCoord cameraUp;
    float upDotZ = m_camera.worldup[0] * cameraZ[0] + m_camera.worldup[1] * cameraZ[1] + m_camera.worldup[2] * cameraZ[2];
    for (int i = 0; i < 3; ++i) {
        cameraUp[i] = m_camera.worldup[i] - upDotZ * cameraZ[i];
    }
    // Calculate Y axis in camera space
    GzCoord cameraY;
    float cameraUpNorm = sqrt(cameraUp[0] * cameraUp[0] + cameraUp[1] * cameraUp[1] + cameraUp[2] * cameraUp[2]);
    for (int i = 0; i < 3; ++i) {
        cameraY[i] = cameraUp[i] / cameraUpNorm;
    }
    // Calculate X axis in camera space
    GzCoord cameraX;
    cameraX[0] = cameraY[1] * cameraZ[2] - cameraY[2] * cameraZ[1];
    cameraX[1] = cameraY[2] * cameraZ[0] - cameraY[0] * cameraZ[2];
    cameraX[2] = cameraY[0] * cameraZ[1] - cameraY[1] * cameraZ[0];
    // Fill into Xiw
    m_camera.Xiw[0][0] = cameraX[0];
    m_camera.Xiw[0][1] = cameraX[1];
    m_camera.Xiw[0][2] = cameraX[2];
    m_camera.Xiw[0][3] = -1.0f * (cameraX[0] * m_camera.position[0] + cameraX[1] * m_camera.position[1] + cameraX[2] * m_camera.position[2]);
    m_camera.Xiw[1][0] = cameraY[0];
    m_camera.Xiw[1][1] = cameraY[1];
    m_camera.Xiw[1][2] = cameraY[2];
    m_camera.Xiw[1][3] = -1.0f * (cameraY[0] * m_camera.position[0] + cameraY[1] * m_camera.position[1] + cameraY[2] * m_camera.position[2]);
    m_camera.Xiw[2][0] = cameraZ[0];
    m_camera.Xiw[2][1] = cameraZ[1];
    m_camera.Xiw[2][2] = cameraZ[2];
    m_camera.Xiw[2][3] = -1.0f * (cameraZ[0] * m_camera.position[0] + cameraZ[1] * m_camera.position[1] + cameraZ[2] * m_camera.position[2]);
    m_camera.Xiw[3][0] = 0;
    m_camera.Xiw[3][1] = 0;
    m_camera.Xiw[3][2] = 0;
    m_camera.Xiw[3][3] = 1;

    // Init Ximage
    int status = GZ_SUCCESS;
    // Push Xsp
    status |= GzPushMatrix(Xsp);
    // Push Xpi and Xiw
    status |= GzPushMatrix(m_camera.Xpi);
    status |= GzPushMatrix(m_camera.Xiw);

    return status;
}

int GzRender::GzPutCamera(GzCamera camera)
{
/* HW 3.8 
/*- overwrite renderer camera structure with new camera definition
*/
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            m_camera.Xiw[i][j] = camera.Xiw[i][j];
            m_camera.Xpi[i][j] = camera.Xpi[i][j];
        }
    }

    for (int i = 0; i < 3; ++i) {
        m_camera.position[i] = camera.position[i];
        m_camera.lookat[i] = camera.lookat[i];
        m_camera.worldup[i] = camera.worldup[i];
    }
    m_camera.FOV = camera.FOV;

	return GZ_SUCCESS;	
}

int GzRender::GzPushMatrix(GzMatrix	matrix)
{
/* HW 3.9 
- push a matrix onto the Ximage stack
- check for stack overflow
*/
    // Check for stack overflow
    if (matlevel > MATLEVELS) return GZ_FAILURE;
	// If the stack is empty
	if (matlevel == -1) {
	    for (int i = 0; i < 4; ++i) {
	        for (int j = 0; j < 4; ++j) {
	            Ximage[0][i][j] = matrix[i][j];  // push into the first stack
	        }
	    }
	} else {
	// Multiply the new `matrix` by the top of the stack
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                float element = 0;
                for (int k = 0; k < 4; ++k) {
                    element += Ximage[matlevel][i][k] * matrix[k][j];
                }
                Ximage[matlevel + 1][i][j] = element;
            }
        }
	}
	// Increase `matlevel`
	matlevel++;
//	// and push it into the stack
//	Ximage[matlevel] = newTosMatrix
	return GZ_SUCCESS;
}

int GzRender::GzPopMatrix()
{
/* HW 3.10
- pop a matrix off the Ximage stack
- check for stack underflow
*/
    // Check for stack underflow
    if (matlevel <= -1) return GZ_FAILURE;
    // Otherwise, decrease `matlevel`
    matlevel--;
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

bool discardTriangle(const vector<vector<float> > & triangle) {
    for (int i = 0; i < 3; ++i) {  // for 3 vertices in the given triangle
        if (triangle[i][2] < 0.0f) return true;  // only check z value
    }
    return false;
}

int GzRender::GzPutTriangle(int numParts, GzToken *nameList, GzPointer *valueList)
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
            GzCoord *vertices = (GzCoord *) valueList[0];
            // Change vertices into homogeneous to easy calculation
            float verticesHomogeneous[3][4];
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    verticesHomogeneous[i][j] = vertices[i][j];
                }
                verticesHomogeneous[i][3] = 1;
            }
            // Transform vertex coordinates using current TOS matrix in stack
            float resultHomogeneous[3][4];
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 4; ++j) {
                    float element = 0;
                    for (int k = 0; k < 4; ++k) {
                        element += Ximage[matlevel][j][k] * verticesHomogeneous[i][k];
                    }
                    resultHomogeneous[i][j] = element;  // Transposed of original matrix multiplication
                }
            }
            // Write back to vertices, change from homogeneous format (4D) back to Cartesian format (3D)
            vector<vector<float>> verticesVector(3, vector<float>(3));
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    if (resultHomogeneous[i][3] == 0) {  // avoid divided by 0
                        resultHomogeneous[i][3] += 0.00001;  // add a tiny number
                    }
                    verticesVector[i][j] = resultHomogeneous[i][j] / resultHomogeneous[i][3];
                    vertices[i][j] = verticesVector[i][j];
                }
            }

            // Discard any triangle with verts behind view plane (z < 0)
            if (discardTriangle(verticesVector)) continue;

            // Rasterization
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

