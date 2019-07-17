/**
 * @file gradienteconj.h
 * @author Brendon Henrique de Paula da Silva GRR20170203
 * @date 16 Setembro 2018
 * @brief Funções utilizadas para o método de gradientes conjulgados com pré condicionadores
 * 
 */
#ifndef __GRADIENTECONJ__H
#define __GRADIENTECONJ__H

double generateRandomA(unsigned int i, unsigned int j, unsigned int k);

double generateRandomB(unsigned int k);

double *gradienteConjPreCond(double *A, double *B, int n, double *X,
							 double *Minvertida, int it, double *tIteracao, int *totalIt, double EPS,int numDiagonais);

#endif