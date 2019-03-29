#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

struct Wrapper{
    float energysum;
    float phireturn[100][100];
  };


float f(float * in_phi, float * out_phi, float * out_en, int ite);
int modulo(int x,int N);



int main(void){


  return 0;
}

float f(float * in_phi, float * out_phi, float * out_en, int ite){
  float a,k,dx,dt,m;
  int iterations, dimensions,i,j,it;
  a = 0.1;
  k = 0.1;
  m = 0.1;
  dx = 1.0;
  dt = 2.0;
  iterations = ite;
  dimensions = 100;
  float phi[100][100] = {{0}};
  memcpy(phi,in_phi,sizeof(phi));
  float mu[100][100] = {{0}};
  float lastphi[100][100] = {{0}};
  float en[100][100] = {{0}};
  float ensum = 0;


  /*FILE *f = fopen("anim.txt", "a+");
  for(i=0 ; i<dimensions ; i++){
    for(j=0 ; j<dimensions ; j++){
      fprintf(f, "%lf ", phi[i][j]);
    }
    fprintf(f, "\n");
  }
  fprintf(f, "- \n");*/
  for(it=0 ; it<iterations ; it++){
    /*memcpy(lastphi,phi,sizeof(lastphi));*/
    for(i=0 ; i<dimensions ; i++){
      for(j=0 ; j<dimensions ; j++){
        lastphi[i][j] = phi[i][j];
      }
    }
    for(i=0 ; i<dimensions ; i++){
      for(j=0 ; j<dimensions ; j++){
        mu[i][j] = -a*lastphi[i][j] + a*lastphi[i][j]*lastphi[i][j]*lastphi[i][j] -(k/(dx*dx))* \
                  (lastphi[modulo((i+1),dimensions)][j] \
                  +lastphi[modulo((i-1),dimensions)][j] \
                  +lastphi[i][modulo((j+1),dimensions)] \
                  +lastphi[i][modulo((j-1),dimensions)] \
                  -4*lastphi[i][j]);
      }
    }
    for(i=0 ; i<dimensions ; i++){
      for(j=0 ; j<dimensions ; j++){
        phi[i][j] = lastphi[i][j] + ((m*dt)/(dx*dx))* \
                    (mu[modulo((i+1),dimensions)][j] \
                    +mu[modulo((i-1),dimensions)][j] \
                    +mu[i][modulo((j+1),dimensions)] \
                    +mu[i][modulo((j-1),dimensions)] \
                    -4*mu[i][j]);
      }
    }
    for(i=0 ; i<dimensions ; i++){
      for(j=0 ; j<dimensions ; j++){
        en[i][j] = -(a/2)*pow(phi[i][j],2) + (a/4)*pow(phi[i][j],4) + (k/2)* \
                    (pow((phi[modulo((i+1),dimensions)][j] \
                    -phi[modulo((i-1),dimensions)][j])/(2*dx), 2) \
                    +pow((phi[i][modulo((j+1),dimensions)] \
                    -phi[i][modulo((j-1),dimensions)])/(2*dx),2));
      }
    }
  }
  ensum = 0;
  for(i=0 ; i<dimensions ; i++){
        for(j=0 ; j<dimensions ; j++){
          ensum=ensum+en[i][j];
        }
      }

  memcpy(out_phi,phi,sizeof(phi));
  memcpy(out_en,en,sizeof(en));
  return ensum;
}

int modulo(int x,int N){
    return (x % N + N) %N;
}
