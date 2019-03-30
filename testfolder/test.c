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
int j(float * in_rhoj,float * out_phij, float * out_conv,float * out_Exj, float * out_Eyj, float treshold, int field);
int g(float * in_rhog, float * out_phig, float * out_convg, float * out_Exg, float * out_Eyg, float tresholdg, int field);
int r(float * in_rhog, float * out_phig, float * out_convg, float * out_Exg, float * out_Eyg, float tresholdg, int field, float dphi);


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

int j(float * in_rhoj,float * out_phij, float * out_conv, float * out_Exj, float * out_Eyj, float treshold, int field){

  int i,j,k, it;
  int dimensions = 50;
  int iterations = 10000;
  float phi[50][50][50] ={{{0}}};
  float lastphi[50][50][50] ={{{0}}};
  float rho[50][50][50] = {{{0}}};
  float Ex[50][50][50] = {{{0}}};
  float Ey[50][50][50] = {{{0}}};
  float Ez[50][50][50] = {{{0}}};
  float Mx[50][50][50] = {{{0}}};
  float My[50][50][50] = {{{0}}};
  double sum = 0;
  double diff = 1;
  float conv[10000] = {0};
  it = 0;

  memcpy(rho,in_rhoj,sizeof(rho));


  while(diff>treshold){
    memcpy(lastphi,phi,sizeof(phi));
    sum =0;


    for(i=1 ; i<dimensions-1 ; i++){
      for(j=1 ; j<dimensions-1 ; j++){
        for(k=1 ; k<dimensions-1 ; k++){
          phi[i][j][k] = (1.0/6.0)*(lastphi[modulo(i+1,dimensions)][j][k]+lastphi[modulo(i-1,dimensions)][j][k] \
                                  +lastphi[i][modulo(j+1,dimensions)][k]+lastphi[i][modulo(j-1,dimensions)][k] \
                                  +lastphi[i][j][modulo(k+1,dimensions)]+lastphi[i][j][modulo(k-1,dimensions)] \
                                  +rho[i][j][k]);
        }
      }
    }
    for(i=0 ; i<dimensions ; i++){
      for(j=0 ; j<dimensions ; j++){
        for(k=0 ; k<dimensions ; k++){
          sum+= fabs(lastphi[i][j][k]-phi[i][j][k]);
        }
      }
    }
    diff = sum;
    conv[it] = diff;
    it++;
  }
  if( field == 0){

  for(i=0 ; i<dimensions ; i++){
      for(j=0 ; j<dimensions ; j++){
        for(k=0 ; k<dimensions ; k++){
          Ex[i][j][k] = -1*((phi[modulo(i+1,dimensions)][j][k]-phi[modulo(i-1,dimensions)][j][k])/2);
          Ey[i][j][k] = -1*((phi[i][modulo(j+1,dimensions)][k]-phi[i][modulo(j-1,dimensions)][k])/2);
          Ez[i][j][k] = -1*((phi[i][j][modulo(k+1,dimensions)]-phi[i][j][modulo(k-1,dimensions)])/2);
        }
      }
    }
    memcpy(out_Exj,Ex,sizeof(Ex));
    memcpy(out_Eyj,Ey,sizeof(Ey));
  }
  else{
  for(i=0 ; i<dimensions ; i++){
      for(j=0 ; j<dimensions ; j++){
        for(k=0 ; k<dimensions ; k++){
          My[i][j][k] = -1*((phi[modulo(i+1,dimensions)][j][k]-phi[modulo(i-1,dimensions)][j][k])/2);
          Mx[i][j][k] = 1*((phi[i][modulo(j+1,dimensions)][k]-phi[i][modulo(j-1,dimensions)][k])/2);
        }
      }
    }
    memcpy(out_Exj,Mx,sizeof(Mx));
    memcpy(out_Eyj,My,sizeof(My));
  }

    memcpy(out_conv,conv,sizeof(conv));
    memcpy(out_phij,phi,sizeof(phi));

  return it;

}

int g(float * in_rhog,float * out_phig, float * out_convg, float * out_Exg, float * out_Eyg, float tresholdg, int field){

  int i,j,k, it;
  int dimensions = 50;
  int iterations = 10000;
  float phi[50][50][50] ={{{0}}};
  float lastphi[50][50][50] ={{{0}}};
  float rho[50][50][50] = {{{0}}};
  float Ex[50][50][50] = {{{0}}};
  float Ey[50][50][50] = {{{0}}};
  float Ez[50][50][50] = {{{0}}};
  float Mx[50][50][50] = {{{0}}};
  float My[50][50][50] = {{{0}}};
  double sum = 0;
  double diff = 1;
  float conv[10000] = {0};
  it = 0;

  memcpy(rho,in_rhog,sizeof(rho));


  while(diff>tresholdg){
    memcpy(lastphi,phi,sizeof(phi));
    sum =0;


    for(i=1 ; i<dimensions-1 ; i++){
      for(j=1 ; j<dimensions-1 ; j++){
        for(k=1 ; k<dimensions-1 ; k++){
          phi[i][j][k] = (1.0/6.0)*(phi[modulo(i+1,dimensions)][j][k]+phi[modulo(i-1,dimensions)][j][k] \
                                  +phi[i][modulo(j+1,dimensions)][k]+phi[i][modulo(j-1,dimensions)][k] \
                                  +phi[i][j][modulo(k+1,dimensions)]+phi[i][j][modulo(k-1,dimensions)] \
                                  +rho[i][j][k]);
        }
      }
    }
    for(i=0 ; i<dimensions ; i++){
      for(j=0 ; j<dimensions ; j++){
        for(k=0 ; k<dimensions ; k++){
          sum+= fabs(lastphi[i][j][k]-phi[i][j][k]);
        }
      }
    }

    diff = sum;
    conv[it] = diff;
    it++;
  }

  if( field == 0){

  for(i=0 ; i<dimensions ; i++){
      for(j=0 ; j<dimensions ; j++){
        for(k=0 ; k<dimensions ; k++){
          Ex[i][j][k] = -1*((phi[modulo(i+1,dimensions)][j][k]-phi[modulo(i-1,dimensions)][j][k])/2);
          Ey[i][j][k] = -1*((phi[i][modulo(j+1,dimensions)][k]-phi[i][modulo(j-1,dimensions)][k])/2);
          Ez[i][j][k] = -1*((phi[i][j][modulo(k+1,dimensions)]-phi[i][j][modulo(k-1,dimensions)])/2);
        }
      }
    }
    memcpy(out_Exg,Ex,sizeof(Ex));
    memcpy(out_Eyg,Ey,sizeof(Ey));
  }
  else{
  for(i=0 ; i<dimensions ; i++){
      for(j=0 ; j<dimensions ; j++){
        for(k=0 ; k<dimensions ; k++){
          My[i][j][k] = -1*((phi[modulo(i+1,dimensions)][j][k]-phi[modulo(i-1,dimensions)][j][k])/2);
          Mx[i][j][k] = 1*((phi[i][modulo(j+1,dimensions)][k]-phi[i][modulo(j-1,dimensions)][k])/2);
        }
      }
    }
    memcpy(out_Exg,Mx,sizeof(Ex));
    memcpy(out_Eyg,My,sizeof(Ey));
  }
    memcpy(out_convg,conv,sizeof(conv));
    memcpy(out_phig,phi,sizeof(phi));

  return it;

}

int r(float * in_rhog,float * out_phig, float * out_convg, float * out_Exg, float * out_Eyg, float tresholdg, int field, float dphi){

  int i,j,k, it;
  int dimensions = 50;
  int iterations = 100000;
  float phi[50][50][50] ={{{0}}};
  float lastphi[50][50][50] ={{{0}}};
  float rho[50][50][50] = {{{0}}};
  float Ex[50][50][50] = {{{0}}};
  float Ey[50][50][50] = {{{0}}};
  float Ez[50][50][50] = {{{0}}};
  float Mx[50][50][50] = {{{0}}};
  float My[50][50][50] = {{{0}}};
  double sum = 0;
  double diff = 1;
  float conv[10000] = {0};
  it = 0;

  memcpy(rho,in_rhog,sizeof(rho));


  while(diff>tresholdg){
    memcpy(lastphi,phi,sizeof(phi));
    sum =0;


    for(i=1 ; i<dimensions-1 ; i++){
      for(j=1 ; j<dimensions-1 ; j++){
        for(k=1 ; k<dimensions-1 ; k++){
          phi[i][j][k] = dphi*(1.0/6.0)*(phi[modulo(i+1,dimensions)][j][k]+phi[modulo(i-1,dimensions)][j][k] \
                                  +phi[i][modulo(j+1,dimensions)][k]+phi[i][modulo(j-1,dimensions)][k] \
                                  +phi[i][j][modulo(k+1,dimensions)]+phi[i][j][modulo(k-1,dimensions)] \
                                  +rho[i][j][k]) \
                                  +(1-dphi)*phi[i][j][k];
        }
      }
    }
    for(i=0 ; i<dimensions ; i++){
      for(j=0 ; j<dimensions ; j++){
        for(k=0 ; k<dimensions ; k++){
          sum+= fabs(lastphi[i][j][k]-phi[i][j][k]);
        }
      }
    }

    diff = sum;
    conv[it] = diff;
    it++;
  }

  if( field == 0){

  for(i=0 ; i<dimensions ; i++){
      for(j=0 ; j<dimensions ; j++){
        for(k=0 ; k<dimensions ; k++){
          Ex[i][j][k] = -1*((phi[modulo(i+1,dimensions)][j][k]-phi[modulo(i-1,dimensions)][j][k])/2);
          Ey[i][j][k] = -1*((phi[i][modulo(j+1,dimensions)][k]-phi[i][modulo(j-1,dimensions)][k])/2);
          Ez[i][j][k] = -1*((phi[i][j][modulo(k+1,dimensions)]-phi[i][j][modulo(k-1,dimensions)])/2);
        }
      }
    }
    memcpy(out_Exg,Ex,sizeof(Ex));
    memcpy(out_Eyg,Ey,sizeof(Ey));
  }
  else{
  for(i=0 ; i<dimensions ; i++){
      for(j=0 ; j<dimensions ; j++){
        for(k=0 ; k<dimensions ; k++){
          My[i][j][k] = -1*((phi[modulo(i+1,dimensions)][j][k]-phi[modulo(i-1,dimensions)][j][k])/2);
          Mx[i][j][k] = 1*((phi[i][modulo(j+1,dimensions)][k]-phi[i][modulo(j-1,dimensions)][k])/2);
        }
      }
    }
    memcpy(out_Exg,Mx,sizeof(Ex));
    memcpy(out_Eyg,My,sizeof(Ey));
  }
    memcpy(out_convg,conv,sizeof(conv));
    memcpy(out_phig,phi,sizeof(phi));

  return it;

}


int modulo(int x,int N){
    return (x % N + N) %N;
}
