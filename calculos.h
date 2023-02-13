#ifndef CALCULOS_H
#define CALCULOS_H

double LeDados_col1(char *filename, double *vetor, int n);

double LeDados_col1_1d(char *filename, double *vetor, int n);

double LeDados_col2(char *filename, double *vetor, int n);

double LeDados_col2_lim(char *filename, double *vetor, int n, int lim);

int ind_max(int n, double *dados, double *max);

int ind_min(int n, double *dados, double *min);

double* derivadaCentral(int n, double dt, double *vm);

double* derivadaPraFrente(int n, double dt, double *vm);

double* derivadaPraTras(int n, double dt, double *vm);

double APD(int n, double dt, double *vt, double *t, double *APD, int tipoAPD);

double Max_dVdt(int n, double dt, double *vm);

double Amplitude(int n, double *vm, double *pico, double *rest);

double MediaDesvio(int n, double *vm, double *media, double *desvio);

double erroNorma1(double *sd, double *dd, int size);

double erroNorma2(double *sd, double *dd, int size);

double erroNormaMax(double *sd, double *dd, int size);

double* razaoX_Y(double *cd, double *id, int size);

double alinhaPAs(double *novo, double *PA1, double *PA2, int size1, int size2, double dt1, double dt2);

#endif