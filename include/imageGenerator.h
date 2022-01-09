#ifndef IMAGEGENERATOR_H
#define IMAGEGENERATOR_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <math.h>
#include "star.h"


int createImg(star* s, int nb_star, double size);
void generateBitmapImage(unsigned char* image, int height, int width, const char* imageFileName);
unsigned char* createBitmapFileHeader(int height, int stride);
unsigned char* createBitmapInfoHeader(int height, int width);



#endif