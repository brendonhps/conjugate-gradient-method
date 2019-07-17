/**
 * @file utils.h
 * @author Brendon Henrique de Paula da Silva GRR20170203
 * @date 16 Setembro 2018
 * @brief Funções auxiliares
 * 
 */
#ifndef __UTILS_H__
#define __UTILS_H__

#include <stdlib.h>
#include <sys/time.h>

double timestamp(void);

void entradaDados(int opcao, int argc, char **argv, int *n, int *k, float *w, int *i, double *EPS, char **arquivoSaida);

void inicializaVetores(int n, int i, double *X, double *A, double *TxA,double *TxB, int k, double *residuo);

void generateMatrizes(int n, int k, double *A, double *B);

void calcTransposta(int n, double *TxB, double *A, double *B, double *TxA,int k,int numDiagonais);

void calcMatrizPrecondicionadora(float w, double *M, int n, double *TxA, int numDiagonais);

void calculaResiduo(double *A, double *B, double *X, int n, double *residuo, double *valorResiduo, double *tempoResiduo,int k);

#endif // __UTILS_H__
