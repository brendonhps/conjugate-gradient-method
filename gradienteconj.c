/**
 * @file gradienteconj.c
 * @author Brendon Henrique de Paula da Silva GRR20170203
 * @date 16 Setembro 2018
 * @brief Funções utilizadas para o método de gradientes conjulgados com pré condicionadores
 * 
 */

#include <stdio.h>
#include <math.h>
#include "utils.h"
#include "gradienteconj.h"

/**
 * @brief Função que gera os coeficientes de um sistema linear k-diagonal
 *
 * @param i,j: coordenadas do elemento a ser calculado (0<=i,j<n)
 * @param k: numero de diagonais da matriz A
 */
inline double generateRandomA(unsigned int i, unsigned int j, unsigned int k)
{
	static double invRandMax = 1.0 / (double)RAND_MAX;
	return ((i == j) ? ((double)(k << 1)) : (1.0)) * ((double)rand() * invRandMax);
}

/**
 * @brief Função que gera os termos independentes de um sistema linear k-diagonal
 *
 * @param k: numero de diagonais da matriz A
 */
inline double generateRandomB(unsigned int k)
{
	static double invRandMax = 1.0 / (double)RAND_MAX;
	return ((double)(k << 2)) * ((double)rand() * invRandMax);
}

/**
 * @brief Implementação do Método Gradiente Conjulgado com/sem Pré-Condicionador
 * @param A Matriz A .
 * @param B Matriz B.
 * @param n Dimensão da matriz A.
 * @param X Vetor de Solução X (Iniciado com zeros).
 * @param Minvertida Vetor com a diagonal da matriz pré-condicionadora.
 * @param it Número de iterações passado por linha de comando.
 * @param tIteracao Tempo de todas as iterações do método.
 * @param totalIt Total de iterações realizadas pelo método.
 * @param EPS Valor do EPS que pode ser passado por linha de comando ou não.
 * @return Vetor de doubles com o erro aproximado após cada iteração.
 */

double *gradienteConjPreCond(double *A, double *B, int n, double *X,
							 double *M, int it, double *tIteracao, int *totalIt, double EPS, int numDiagonais)
{

	double *z, *v, *y, *XNew, *iterK;
	double s = 0.0, m = 0.0, dif, aux2 = 0.0;
	double aux = 0.0;
	double aux1 = 0.0;
	double erro = 100.0;
	int pos = numDiagonais - 1;
	int deslocamento;

	v = (double *)malloc(sizeof(double) * n);
	y = (double *)malloc(sizeof(double) * n);
	z = (double *)malloc(sizeof(double) * n);
	XNew = (double *)malloc(sizeof(double) * n);
	iterK = (double *)malloc(sizeof(double) * it);
	for (int i = 0; i < n; i++)
	{
		z[i] = 0.0;
	}
#pragma GCC vector aligned
#pragma GCC ivdep
	for (int i = 0; i < n; i += 4)
	{
		y[i] = M[i] * B[i];
		y[i + 1] = M[i + 1] * B[i + 1];
		y[i + 2] = M[i + 2] * B[i + 2];
		y[i + 3] = M[i + 3] * B[i + 3];
		v[i] = y[i];
		v[i + 1] = y[i+1];
		v[i + 2] = y[i+2];
		v[i + 3] = y[i+3];
		aux += y[i] * B[i];
		aux += y[i + 1] * B[i + 1];
		aux += y[i + 2] * B[i + 2];
		aux += y[i + 3] * B[i + 3];
	}

	*tIteracao = timestamp();
	int k;

	for (k = 0; k < it; k++)
	{

		aux2 = 0.0;
		aux1 = 0.0;
		dif = 0.0;

		for (int linha = 0; linha < numDiagonais; linha++)
		{
			for (int i = linha; i > 0; i--)
			{
				deslocamento = linha - i;
				z[linha] += A[deslocamento * numDiagonais + i] * v[deslocamento];
			}

			deslocamento = 0;
			int tmp = linha * numDiagonais;
			for (int indice = linha; indice < linha + numDiagonais && indice < n; indice++)
			{

				z[linha] += A[tmp + deslocamento] * v[indice];
				deslocamento++;
			}

			aux2 += v[linha] * z[linha];
		}

		for (int linha = numDiagonais; linha < n; linha++)
		{

			for (int i = pos; i > 0; i--)
			{
				deslocamento = linha - i;
				z[linha] += A[deslocamento * numDiagonais + i] * v[deslocamento];
			}

			deslocamento = 0;

			int tmp = linha * numDiagonais;
			for (int indice = linha; indice < linha + numDiagonais && indice < n; indice++)
			{

				z[linha] += A[tmp + deslocamento] * v[indice];
				deslocamento++;
			}

			aux2 += v[linha] * z[linha];
		}

		s = aux / aux2;

		for (int i = 0; i < n; i++)
		{

			XNew[i] = X[i] + s * v[i];

			if (fabs(XNew[i] - X[i]) > dif)
			{
				dif = fabs(XNew[i] - X[i]);
			}
			X[i] = XNew[i];
			B[i] = B[i] - (s * z[i]);

			y[i] = M[i] * B[i];

			aux1 += y[i] * B[i];
			z[i] = 0.0;
		}

		m = aux1 / aux;
		aux = aux1;

#pragma GCC vector aligned
#pragma GCC ivdep
		for (int i = 0; i < n; i += 4)
		{
			v[i] = y[i] + m * v[i];
			v[i + 1] = y[i + 1] + m * v[i + 1];
			v[i + 2] = y[i + 2] + m * v[i + 2];
			v[i + 3] = y[i + 3] + m * v[i + 3];
		}
		iterK[k] = dif;
	}
	// }
	// else
	// {

	// 	while ((k < it) && (erro > EPS))
	// 	{

	// 		aux2 = 0.0;
	// 		aux1 = 0.0;
	// 		dif = 0.0;
	// 		for (int linha = 0; linha < n; linha++)
	// 		{
	// 			for (int i = 0; i < pos; i++)
	// 			{
	// 				z[linha] += A[(linha - i) * numDiagonais + i] * v[linha - i];
	// 			}

	// 			deslocamento = 0;
	// 			for (int indice = linha; indice < linha + numDiagonais; indice++)
	// 			{

	// 				if (A[(linha * numDiagonais) + deslocamento] != 0.0)
	// 				{
	// 					z[linha] += A[(linha * numDiagonais) + deslocamento] * v[indice];
	// 				}
	// 				deslocamento++;
	// 			}

	// 			aux2 += v[linha] * z[linha];
	// 		}
	// 		s = aux / aux2;

	// 		for (int i = 0; i < n; i++)
	// 		{

	// 			XNew[i] = X[i] + s * v[i];
	// 			if (fabs(XNew[i] - X[i]) > dif)
	// 			{
	// 				dif = fabs(XNew[i] - X[i]);
	// 				erro = dif * (1.0 / fabs(XNew[i]));
	// 			}
	// 			X[i] = XNew[i];
	// 			r[i] = r[i] - (s * z[i]);
	// 			y[i] = Minvertida[i] * r[i];
	// 			aux1 += y[i] * r[i];
	// 		}

	// 		m = aux1 / aux;
	// 		aux = aux1;

	// 		for (int i = 0; i < n; i++)
	// 		{
	// 			v[i] = y[i] + m * v[i];
	// 			z[i] = 0.0;
	// 		}
	// 		iterK[k] = dif;
	// 		k++;
	// 	}
	// }

	*tIteracao = timestamp() - *tIteracao;
	*totalIt = k - 1;

	free(z);
	free(v);
	free(y);
	free(XNew);
	return iterK;
}
