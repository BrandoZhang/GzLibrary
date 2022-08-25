/* Texture functions for cs580 GzLib	*/
#include    "stdafx.h" 
#include	"stdio.h"
#include	"Gz.h"

GzColor	*image=NULL;
int xs, ys;
int reset = 1;

/* Image texture function */
int tex_fun(float u, float v, GzColor color)
{
  unsigned char		pixel[3];
  unsigned char     dummy;
  char  		foo[8];
  int   		i, j;
  FILE			*fd;

  if (reset) {          /* open and load texture file */
    fd = fopen("texture", "rb");
    if (fd == NULL) {
      fprintf(stderr, "texture file not found\n");
      exit(-1);
    }
    fscanf(fd, "%s %d %d %c", foo, &xs, &ys, &dummy);
    image = (GzColor*)malloc(sizeof(GzColor)*(xs+1)*(ys+1));
    if (image == NULL) {
      fprintf(stderr, "malloc for texture image failed\n");
      exit(-1);
    }

    for (i = 0; i < xs*ys; i++) {	/* create array of GzColor values */
      fread(pixel, sizeof(pixel), 1, fd);
      image[i][RED] = (float)((int)pixel[RED]) * (1.0 / 255.0);
      image[i][GREEN] = (float)((int)pixel[GREEN]) * (1.0 / 255.0);
      image[i][BLUE] = (float)((int)pixel[BLUE]) * (1.0 / 255.0);
      }

    reset = 0;          /* init is done */
	fclose(fd);
  }

/* bounds-test u,v to make sure nothing will overflow image array bounds */
  u = min(max(u, 0.0f), 1.0f);
  v = min(max(v, 0.0f), 1.0f);
  u *= (xs - 1);
  v *= (ys - 1);

/* determine texture cell corner values and perform bi-linear interpolation */
/*
 *  A(x1, y1) -------------- B(x2, y1)
 *            |         |    |
 *            |         | t  |
 *            |--- s ---|    |
 *            |          p   |
 *  D(x1, y2) -------------- C(x2, y2)
 *
 */
  int x1 = (int) floor(u);  // left-most x
  int x2 = (int) ceil(u);  // right-most x
  int y1 = (int) floor(v);  // top y
  int y2 = (int) ceil(v);  // bottom y

  int idxA = y1 * xs + x1;
  int idxB = y1 * xs + x2;
  int idxC = y2 * xs + x2;
  int idxD = y2 * xs + x1;
  GzColor colorA = {image[idxA][RED], image[idxA][GREEN], image[idxA][BLUE]};
  GzColor colorB = {image[idxB][RED], image[idxB][GREEN], image[idxB][BLUE]};
  GzColor colorC = {image[idxC][RED], image[idxC][GREEN], image[idxC][BLUE]};
  GzColor colorD = {image[idxD][RED], image[idxD][GREEN], image[idxD][BLUE]};

  float s = u - x1;
  float t = v - y1;

/* set color to interpolated GzColor value and return */
  color[RED] = s * t * colorC[RED] + (1.0f - s) * t * colorD[RED] + s * (1.0f - t) * colorB[RED] + (1.0f - s) * (1.0f - t) * colorA[RED];
  color[GREEN] = s * t * colorC[GREEN] + (1.0f - s) * t * colorD[GREEN] + s * (1.0f - t) * colorB[GREEN] + (1.0f - s) * (1.0f - t) * colorA[GREEN];
  color[BLUE] = s * t * colorC[BLUE] + (1.0f - s) * t * colorD[BLUE] + s * (1.0f - t) * colorB[BLUE] + (1.0f - s) * (1.0f - t) * colorA[BLUE];

	return GZ_SUCCESS;
}

/* Procedural texture function */
int ptex_fun(float u, float v, GzColor color)
{
    u = min(max(u, 0.0f), 1.0f);
    v = min(max(v, 0.0f), 1.0f);

    /* Checkerboard of red and white */
    int checkerboardSize = 6;  // 6x6 checkerboard
    u *= checkerboardSize;
    v *= checkerboardSize;
    int intervalU = (int) floor(u);
    int intervalV = (int) floor(v);
    if ((intervalU % 2) == (intervalV % 2)) {
        // Both u and v fall into the same (odd/odd, even/even) intervals, set red
        color[RED] = 1.0f;
        color[GREEN] = 0.0f;
        color[BLUE] = 0.0f;
    } else {
        // Both u and v fall into different (odd/even, even/odd) intervals, set white
        color[RED] = 1.0f;
        color[GREEN] = 1.0f;
        color[BLUE] = 1.0f;
    }
	return GZ_SUCCESS;
}

/* Free texture memory */
int GzFreeTexture()
{
	if(image!=NULL)
		free(image);
	return GZ_SUCCESS;
}

