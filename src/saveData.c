#include "star.h"
#include "saveData.h"

void saveSetting(double Tmax, double dt, int nb_star, double softening, int cptCapt, int nb_thread) {

    FILE* f = fopen("./save/setting.txt", "w+");

    fprintf(f, "simulationDuration: %.1lf\n", Tmax);
    fprintf(f, "dt: %.4lf\n", dt);
    fprintf(f, "numberOfStar: %i\n", nb_star);
    fprintf(f, "softenningCoef: %.3lf\n", softening);
    fprintf(f, "counterCapture: %i\n", cptCapt);
    fprintf(f, "numberOfThread: %i\n", nb_thread);
    fprintf(f, "reloadSimulation(0=no/1=yes): %i\n", 0);
    fclose(f);
    return;
}


void saveStar(star* stars, int nb_star) {
    int i;
    FILE* f = fopen("./save/star", "w+");
    star s;
    fprintf(f, "x;y;Vx;Vy;Ax;Ay;m\n");
    for(i = 0; i<nb_star; i++) {
        s = stars[i];
        fprintf(f, "%.5f;%.5f;%.5f;%.5f;%.5f;%.5f;%.5f\n", 
                    s.x ,s.y ,s.vx,s.vy,s.ax,s.ay,s.mass);    
    }



    fclose(f);
    return;
}