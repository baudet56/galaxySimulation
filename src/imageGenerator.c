#include "imageGenerator.h"
#include <stdio.h>
#include <stdlib.h>
int cpt = 0;


const int BYTES_PER_PIXEL = 3;
const int FILE_HEADER_SIZE = 14;
const int INFO_HEADER_SIZE = 40;



int createImg(star* s, int nb_star, double size)
{
    int height = 1000;
    int width = 1000;
    unsigned char* image;
    char str[90];
    char nb[6];
    int i;
    double x, y;

    sprintf(nb, "%d", cpt);
    strcpy(str, "./data/img:");
    strcat(str, nb);
    cpt +=1;
    strcat(str, ".bmp");



    image = (unsigned char*)malloc(3*height*width);/*[width][BYTES_PER_PIXEL];*/
    memset(image, 0, 3*width * height);

    for(i = 0; i<nb_star; i++){
        x = s[i].x;
        y = s[i].y;
        if(fabs(x) < size && fabs(y) < size){
            x = (x/size * width/2)+ width/2;
            y = (y/size * height/2) + height/2;
            image[3*( (int)y*width + (int)x) + 0] = 255;
            image[3*( (int)y*width + (int)x) + 1] = 255;
            image[3*( (int)y*width + (int)x) + 2] = 255;
        }
    }

    generateBitmapImage((unsigned char*) image, height, width, str);

    free(image);

    return 0;
}


void generateBitmapImage (unsigned char* image, int height, int width,const char* imageFileName)
{
    int widthInBytes = width * BYTES_PER_PIXEL;

    int i;
    unsigned char padding[3] = {0, 0, 0};
    int paddingSize = (4 - (widthInBytes) % 4) % 4;

    int stride = (widthInBytes) + paddingSize;

    FILE* imageFile = fopen(imageFileName, "wb");

    unsigned char* fileHeader = createBitmapFileHeader(height, stride);
    unsigned char* infoHeader = createBitmapInfoHeader(height, width);

    fwrite(fileHeader, 1, FILE_HEADER_SIZE, imageFile);
    fwrite(infoHeader, 1, INFO_HEADER_SIZE, imageFile);

    for (i = 0; i < height; i++) {
        fwrite(image + (width * (height - i - 1)*3 ), BYTES_PER_PIXEL, width, imageFile);
        fwrite(padding, 1, paddingSize, imageFile);
    }

    fclose(imageFile);
}

unsigned char* createBitmapFileHeader (int height, int stride)
{
    int fileSize = FILE_HEADER_SIZE + INFO_HEADER_SIZE + (stride * height);

    static unsigned char fileHeader[] = {
        0,0,
        0,0,0,0,
        0,0,0,0,
        0,0,0,0,
    };

    fileHeader[ 0] = (unsigned char)('B');
    fileHeader[ 1] = (unsigned char)('M');
    fileHeader[ 2] = (unsigned char)(fileSize      );
    fileHeader[ 3] = (unsigned char)(fileSize >>  8);
    fileHeader[ 4] = (unsigned char)(fileSize >> 16);
    fileHeader[ 5] = (unsigned char)(fileSize >> 24);
    fileHeader[10] = (unsigned char)(FILE_HEADER_SIZE + INFO_HEADER_SIZE);

    return fileHeader;
}

unsigned char* createBitmapInfoHeader (int height, int width)
{
    static unsigned char infoHeader[] = {
        0,0,0,0,
        0,0,0,0,
        0,0,0,0,
        0,0,
        0,0,
        0,0,0,0,
        0,0,0,0,
        0,0,0,0,
        0,0,0,0,
        0,0,0,0,
        0,0,0,0,
    };

    infoHeader[ 0] = (unsigned char)(INFO_HEADER_SIZE);
    infoHeader[ 4] = (unsigned char)(width      );
    infoHeader[ 5] = (unsigned char)(width >>  8);
    infoHeader[ 6] = (unsigned char)(width >> 16);
    infoHeader[ 7] = (unsigned char)(width >> 24);
    infoHeader[ 8] = (unsigned char)(height      );
    infoHeader[ 9] = (unsigned char)(height >>  8);
    infoHeader[10] = (unsigned char)(height >> 16);
    infoHeader[11] = (unsigned char)(height >> 24);
    infoHeader[12] = (unsigned char)(1);
    infoHeader[14] = (unsigned char)(BYTES_PER_PIXEL*8);

    return infoHeader;
}