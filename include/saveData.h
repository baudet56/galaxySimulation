#ifndef SAVEDATA_H
#define SAVEDATA_H

#include <stdio.h>
#include <stdlib.h>

void saveSetting(double Tmax, double dt, int nb_star, double softening, int cptCapt, int nb_thread);

void saveStar(star* stars, int nb_star);

#endif