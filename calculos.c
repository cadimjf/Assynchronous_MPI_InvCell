#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "calculos.h"


double LeDados_col1(char *filename, double *vetor, int n) {

   FILE *f1;
   double *lixo;
   int i;

   f1 = fopen(filename, "r");
   if(!f1) {
      printf("falha na abertura do arquivo %s\n", filename);
      exit(1);
   }

   lixo = (double*) malloc(n*sizeof(double));

   for(i = 0; i < n; i++) {
      fscanf(f1, "%lf %lf", &vetor[i], &lixo[i]);
   }

   fclose(f1);
   free(lixo);

   return 0.0;   
}

double LeDados_col1_1d(char *filename, double *vetor, int n) {

   FILE *f1;
   int i;

   f1 = fopen(filename, "r");
   if(!f1) {
      printf("falha na abertura do arquivo %s\n", filename);
      exit(1);
   }

   for(i = 0; i < n; i++) {
      fscanf(f1, "%lf", &vetor[i]);
   }

   fclose(f1);

   return 0.0;   
}

double LeDados_col2(char *filename, double *vetor, int n) {

   FILE *f1;
   double *lixo;
   int i;

   f1 = fopen(filename, "r");
   if(!f1) {
      printf("falha na abertura do arquivo %s\n", filename);
      exit(1);
   }

   lixo = (double*) malloc(n*sizeof(double));

   for(i = 0; i < n; i++) {
      fscanf(f1, "%lf %lf", &lixo[i], &vetor[i]);
   }

   fclose(f1);
   free(lixo);

   return 0.0;   
}

double LeDados_col2_lim(char *filename, double *vetor, int n, int lim) {

   FILE *f1;
   double *lixo1,
          *lixo2;
   int i;

   f1 = fopen(filename, "r");
   if(!f1) {
      printf("falha na abertura do arquivo %s\n", filename);
      exit(1);
   }

   lixo1 = (double*) malloc(n*sizeof(double));
   lixo2 = (double*) malloc(lim*sizeof(double));

   for(i = 0; i < n; i++) {
      if(i < lim)
         fscanf(f1, "%lf %lf", &lixo1[i], &lixo2[i]);
      else
         fscanf(f1, "%lf %lf", &lixo1[i], &vetor[i-lim]);
   }

   fclose(f1);
   free(lixo1);
   free(lixo2);

   return 0.0;   
}

int ind_max(int n, double *dados, double *max) {
   
   int indice, i;

   *max = dados[0];
   indice = 0;
   
   for(i = 1; i < n; i++) {
      if(dados[i] > *max) {
         *max = dados[i];
         indice = i;
      }
   }

   return indice;
}

int ind_min(int n, double *dados, double *min) {
   
   int indice, i;
   
   *min = dados[0];
   indice = 0;
   
   for(i = 1; i < n; i++) {
      if(dados[i] < *min) {
         *min = dados[i];
         indice = i;
      }
   }

   return indice;
}

double* derivadaCentral(int n, double dt, double *vm) {

   double *derivada;
   int m;

   derivada = (double*) malloc(n*sizeof(double));

   derivada[0] = (vm[1] - vm[0])/dt;
   derivada[n-1] = (vm[n-1] - vm[n-2])/dt;
 
   for(m = 1; m < n - 1; m++) {
      derivada[m] = (vm[m+1] - vm[m-1])/(2*dt);
   }

   return derivada;
}

double* derivadaPraFrente(int n, double dt, double *vm) {

   double *derivada;
   int m;

   derivada = (double*) malloc(n*sizeof(double));

   for(m = 1; m < n; m++) {
      derivada[m] = (vm[m] - vm[m-1])/(dt);
   }

   return derivada;
}

double* derivadaPraTras(int n, double dt, double *vm) {

   double *derivada;
   int m;

   derivada = (double*) malloc(n*sizeof(double));

   for(m = 1; m < n; m++) {
      derivada[m] = (vm[m-1] - vm[m])/dt;
   }

   return derivada;
}

double APD(int n, double dt, double *vt, double *t, double *APD, int tipoAPD) {
   
   double max,
          ti,
          tf,
          tmin,
          tmax,
          vt_min,
          vt_max,
          PAmax,
          x,
          PAx,
          APDX;

   double *derivadas;

   int i,
       itmax,
       ind;

   // Calculando ti
   
   derivadas = derivadaCentral(n, dt, vt);
   ind = ind_max(n, derivadas, &max);
   ti = t[ind];

   // CALCULO DE tf

   // Encontrando vm_min e vm_max
   ind   = ind_min(n, vt, &vt_min);
   tmin = t[ind];
   itmax = ind_max(n, vt, &vt_max);
   tmax = t[itmax];

   // Calculo de PAmax
   PAmax = vt_max - vt_min;

   // Calculo de PAx
   x = tipoAPD*0.01;
   PAx = (1.0 - x)*PAmax;

   // Encontrando tf
   tf = tmax;
   for(i = itmax; i < n; i++)
      if(vt[i] <= (vt_min + PAx)) {
         tf = t[i];
         break;
      }

   // Calculo do APDX
   APDX = tf - ti;
   
   *APD = APDX;

   return 0.0;
}

double Max_dVdt(int n, double dt, double *vm) {

   double *derivadas = derivadaCentral(n, dt, vm);

   double max;
   ind_max(n, derivadas, &max);

   return max;
}

double Amplitude(int n, double *vm, double *pico, double *rest) {

   double ampl;

   ind_max(n, vm, pico);
   ind_min(n, vm, rest);
 
   ampl = *pico - *rest;

   return ampl;
}

double MediaDesvio(int n, double *vm, double *media, double *desvio) {

   int i;
   double dp, md;

   md = 0.0;
   for(i = 0; i < n; i++) {
      md += vm[i];
   }
   *media = md/n;

   dp = 0.0;
   for(i = 0; i < n; i++) {
      dp += pow((vm[i] - *media), 2.0);
   }
   *desvio = sqrt(dp/n);

   return 0.0;
}

double erroNorma1(double *sd, double *dd, int size) {

   double erro;
   double n1;
   double n2;
   int i;

   for(i = 0; i < size; i++) {
      n1 = fabs(dd[i]-sd[i]);
      n2 = fabs(dd[i]);
   }

   erro = n1/n2;

   return erro;
}

double erroNorma2(double *sd, double *dd, int size) {

   double erro;
   double n1 = 0.0;
   double n2 = 0.0;
   int i;

   for(i = 0; i < size; i++) {
      n1 += pow((dd[i]-sd[i]),2.0);
      n2 += pow(dd[i],2.0);
   }

   erro = sqrt(n1/n2);

   return erro;
}

double erroNormaMax(double *sd, double *dd, int size) {

   double erro, n2, max;
   double n1[size];
   int i, ind;

   n2 = 0.0;
   for(i = 0; i < size; i++) {
      n1[i] = fabs(dd[i]-sd[i]);
      n2 += n1[i];
   }

   max = n1[0];
   ind = 0;
   for(i = 1; i < size; size++) {
      if(n1[i] > max) {
         max = n1[i];
         ind = i;
      }
   }
 
   erro = max/n2;

   free(n1);

   return erro;
}

double* razaoX_Y(double *cd, double *id, int size) {

   double *razao;
   int i;

   razao = (double*) malloc(size*sizeof(double));

   for(i = 0; i < size; i++) {
      razao[i] = (id[i]-cd[i])/cd[i];
   }

   return razao;
}

double alinhaPAs(double *novo, double *PA1, double *PA2, int size1, int size2, double dt1, double dt2) {

   double *diff1,
          *diff2,
          max;
   int i,
       ind1,
       ind2,
       ini,
       fim;


   diff1 = (double*) malloc(size1*sizeof(double));
   diff2 = (double*) malloc(size2*sizeof(double));

   diff1 = derivadaCentral(size1, dt1, PA1);
   ind1  = ind_max(size1, diff1, &max); 
   diff2 = derivadaCentral(size2, dt2, PA2);
   ind2  = ind_max(size2, diff2, &max);

   ini = abs(ind1-ind2);

   fim = ini + size1;

   for(i = ini; i < fim; i++) {
      novo[i-ini] = PA2[i];
   }

   free(diff1);
   free(diff2);

   return 0.0;
}