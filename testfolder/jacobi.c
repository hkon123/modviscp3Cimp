#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


int modulo(int x,int N);

int main(void){

  int i,j,k, it;
  int dimensions = 50;
  int iterations = 10000;
  double phi[50][50][50] ={{{0}}};
  double lastphi[50][50][50] ={{{0}}};
  double rho[50][50][50] = {{{0}}};
  double E[50][50][50] = {{{0}}};
  double sum = 0;


  rho[25][25][25] = 1;


  for(it=0 ; it<iterations ; it++){
    //memcpy(lastphi,phi,sizeof(phi));
    for(i=0 ; i<dimensions ; i++){
      for(j=0 ; j<dimensions ; j++){
        for(k=0 ; k<dimensions ; k++){
          lastphi[i][j][k] = phi[i][j][k];
        }
      }
    }

    sum = 0;
    /*for(i=0 ; i<dimensions ; i++){
      for(j=0 ; j<dimensions ; j++){
        lastphi[0][i][j] = 0;
        lastphi[4][i][j] = 0;
        lastphi[i][0][j] = 0;
        lastphi[i][4][j] = 0;
        lastphi[i][j][0] = 0;
        lastphi[i][j][4] = 0;
      }
    }*/

    for(i=1 ; i<dimensions-1 ; i++){
      for(j=1 ; j<dimensions-1 ; j++){
        for(k=1 ; k<dimensions-1 ; k++){
          phi[i][j][k] = (1.0/6.0)*(phi[modulo(i+1,dimensions)][j][k]+phi[modulo(i-1,dimensions)][j][k] \
                                  +phi[i][modulo(j+1,dimensions)][k]+phi[i][modulo(j-1,dimensions)][k] \
                                  +phi[i][j][modulo(k+1,dimensions)]+phi[i][j][modulo(k-1,dimensions)] \
                                  +rho[i][j][k]);
          //printf("%lf    %lf\n",phi[i][j][k], lastphi[i][j][k]);
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

    if(modulo(it,50)==0){
      printf("%lf    %d\n", sum, it);
    }
  }

  for(i=0 ; i<dimensions ; i++){
      for(j=0 ; j<dimensions ; j++){
        for(k=0 ; k<dimensions ; k++){
          E[i][j][k] = -1*((phi[modulo(i+1,dimensions)][j][k]-phi[modulo(i-1,dimensions)][j][k])/2 \
                          +(phi[i][modulo(j+1,dimensions)][k]-phi[i][modulo(j-1,dimensions)][k])/2 \
                          +(phi[i][j][modulo(k+1,dimensions)]-phi[i][j][modulo(k-1,dimensions)])/2 );
        }
      }
    }



  return 0;
}




int modulo(int x,int N){
    return (x % N + N) %N;
}
