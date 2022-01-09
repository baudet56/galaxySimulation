#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <pthread.h>
#include <stdint.h>

#include "imageGenerator.h"
#include "star.h"
#include "thread_struct.h"
#include "saveData.h"

#define PI 3.1415

int cptCapt = 6;
double softening = 1e-1;
double pourcentageStar = 0.7;
double T = 20;
double dt = 0.003;
int nb_star = 100000;
int nb_thread = 4;
int reload = 0;

void *thread_compute(void* p);
void initStar(star *tab_star);
void loadSetting();
void load(star* tab_star);

int main()
{
    double t;
    star* tab_star;
    int cpt = 0;
    pthread_t *id_list;
    thread_struct *data;
    int thread_index;
    clock_t begin, end;
    double minute = 0;
    loadSetting();

    tab_star = (star* )malloc(sizeof(star)*nb_star);

    if(reload){
        load(tab_star);
    } else {
        initStar(tab_star);
        saveStar(tab_star, nb_star);
    }

    id_list = (pthread_t *)malloc( nb_thread*sizeof(pthread_t) );
    data = (thread_struct *)malloc( nb_thread*sizeof(thread_struct) );

    for(thread_index = 0; thread_index<nb_thread; thread_index++) {
        data[thread_index].tab_star = tab_star;
        data[thread_index].thread_id = thread_index;
    }

    begin = time(NULL);
    printf("start\n");
    for(t = 0; t<T; t+=dt) {

        if(cpt%cptCapt == 0) {
            end = time(NULL);
            if(t) minute = (double)(end-begin)*(T/t-1);
            printf("%.2f%% : time left: %i h %i m %i sec\n", 
            t/T*100, (int)(minute/60/60), (int)(minute/60)%60, (int)(minute)%60);
            createImg(tab_star, nb_star, 3);
        }

        if( cpt%100 == 0 ) saveStar(tab_star, nb_star);

        for(thread_index = 0; thread_index < nb_thread; ++thread_index){
            pthread_create(&id_list[thread_index], NULL, thread_compute,(void *)&data[thread_index]);
        }

        for(thread_index = 0; thread_index < nb_thread; ++thread_index){
            pthread_join(id_list[thread_index], NULL);
        }
        ++cpt;
    }

    saveStar(tab_star, nb_star);

    free(tab_star);
    free(data);
    free(id_list);
    return 0;
}

double invsqrtQuake( double number )
  {
      double y = number;
      double x2 = y * 0.5;
      int64_t i = *(int64_t *) &y;
      // The magic number is for doubles is from https://cs.uwaterloo.ca/~m32rober/rsqrt.pdf
      i = 0x5fe6eb50c7b537a9 - (i >> 1);
      y = *(double *) &i;
      y = y * (1.5 - (x2 * y * y));   // 1st iteration
      //      y  = y * ( 1.5 - ( x2 * y * y ) );   // 2nd iteration, this can be removed
      return y;
  }

void *thread_compute(void* p){
    int j, i;
    double dx, dy, ax, ay;
    star s, ns;
    double F;
    thread_struct* t = (thread_struct*)p;
    double dt2 = dt/2.;
    double softening2 = softening*softening;
    int start = (t->thread_id*nb_star)/nb_thread;
    int end = start + nb_star/nb_thread;
    int nb_star_reduce = (int)( pourcentageStar * nb_star );
    star * tab_star = t->tab_star;
    for(i = start; i < end; ++i){
        s = tab_star[i];
        s.vx += s.ax*dt2;
        s.vy += s.ay*dt2;

        s.x += s.vx*dt;
        s.y += s.vy*dt;

        ax = 0;
        ay = 0;
        for(j = 0; j < nb_star_reduce; ++j){
            if(j != i){
                ns = tab_star[j];
                dx = ns.x - s.x;
                dy = ns.y - s.y;
                F = dx*dx + dy*dy + softening2;
                F = invsqrtQuake(F*F*F);
                
                ax += ns.mass*(F*dx);
                ay += ns.mass*(F*dy);
            }
        }
        s.ax = ax;
        s.ay = ay;
        
        s.vx += ax*dt2;
        s.vy += ay*dt2;
        tab_star[i] = s;
    }
    pthread_exit(0);
}


double getRandomValue(){
    return (double)rand()/RAND_MAX;
}


void initStar(star* tab_star) {
    double u, phi;
    int i;

    srand(time(NULL));
    for(i = 0; i<nb_star; ++i){
        u = getRandomValue();
        phi = 2*PI*u;
        tab_star[i].x = pow(2*getRandomValue(), (double)1./3) * cos(phi);
        tab_star[i].y = pow(2*getRandomValue(), 1./3) * sin(phi);
        tab_star[i].mass = (double)1./nb_star*1./pourcentageStar;
        tab_star[i].vx = pow(2*getRandomValue(), (double)1./3) * ( -sin(phi) ) * 2*PI*dt;
        tab_star[i].vy = pow(2*getRandomValue(), (double)1./3) * cos(phi) * 2*PI*dt;
        tab_star[i].ax = 0;
        tab_star[i].ay = 0;

        if(i%2){
            tab_star[i].x += 1.5;
            tab_star[i].y += 1.5;
        } else {
            tab_star[i].x -= 1.5;
            tab_star[i].y -= 1.5;
        }
    }
    /* Black hole ! */
    tab_star[0].x = 1.5;
    tab_star[0].y = 1.5;
    tab_star[0].mass = 0.1;

    tab_star[1].x = -1.5;
    tab_star[1].y = -1.5;
    tab_star[1].mass = 0.1;
}

void loadSetting() {
    FILE* f = fopen("./save/setting.txt", "r");
    
    fscanf(f, "simulationDuration: %lf\n", &T);
    fscanf(f, "dt: %lf\n", &dt);
    printf("%lf\n", dt);
    fscanf(f, "numberOfStar: %i\n", &nb_star);
    fscanf(f, "softenningCoef: %lf\n", &softening);
    fscanf(f, "counterCapture: %i\n", &cptCapt);
    fscanf(f, "numberOfThread: %i\n", &nb_thread);
    fscanf(f, "reloadSimulation(0=no/1=yes): %i\n", &reload);
    fclose(f);
    return;
}

void load(star* tab_star) {
    FILE* f = fopen("./save/star", "r");
    int i;
    star *s;
    fscanf(f, "x;y;Vx;Vy;Ax;Ay;m\n");

    for(i = 0; i<nb_star; ++i){
        s = &tab_star[i];
        fscanf(f, "%lf;%lf;%lf;%lf;%lf;%lf;%lf\n", 
        &s->x, &s->y, &s->vx, &s->vy, &s->ax, &s->ay, &s->mass);
    }
    fclose(f);
}
