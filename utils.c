/**
 * @file utils.c
 * @author Brendon Henrique de Paula da Silva GRR20170203
 * @date 16 Setembro 2018
 * @brief Funções auxiliares
 * 
 */

#include "utils.h"
#include <math.h>
#include "gradienteconj.h"
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <unistd.h>

/**
 * @brief Função utilizada para calcular o tempo em milisegundos
 *  
 *@return Tempo em milisegundos
 */
double timestamp(void)
{
	struct timeval tp;
	gettimeofday(&tp, NULL);
	return ((double)(tp.tv_sec * 1000.0 + tp.tv_usec / 1000.0));
}

/**
 * @brief Função utilizada para tratar os dados passados pelo usuário por linha de comando
 * 
 * @param opcao Opção passada pelo usuário na linha de comando.
 * @param argc Número de argumentos passado para o programa.
 * @param argv Vetor de string contendo os argumentos.
 * @param n Dimensão da matriz A[n][n].
 * @param k Número de diagonais não nulas.
 * @param w Valor indicando o Pré-Condicionador utilizado.
 * @param i Número de iterações max que o método pode realizar.
 * @param EPS Erro aproximado absoluto máximo.
 * @param arquivoSaida Arquivo de saida passado pelo usuário.
 * @return Sem retorno.
 */
void entradaDados(int opcao, int argc, char **argv, int *n, int *k, float *w, int *i, double *EPS, char **arquivoSaida)
{

	while ((opcao = getopt(argc, argv, "n:k:p:i:e:o:")) != -1)
	{
		switch (opcao)
		{
		case 'n':
			*n = atoi(optarg);
			if (*n < 10)
			{
				fprintf(stderr, "%s\n", "Valor de n tem que ser maior que 10");
				exit(1);
			}
			break;
		case 'k':
			*k = atoi(optarg);
			if ((*k % 2 == 0) || (*k < 1))
			{
				fprintf(stderr, "%s\n", " Valor de k tem que ser ímpar e maior que 1");
				exit(1);
			}
			break;
		case 'p':
			*w = atof(optarg);
			if ((*w < 0.0) || (*w >= 1.0))
			{
				fprintf(stderr, "%s\n", "Pré-Condicionador não reconhecido");
				exit(1);
			}
			break;
		case 'i':
			*i = atoi(optarg);
			break;
		case 'e':
			*EPS = atof(optarg);
			break;
		case 'o':
			*arquivoSaida = optarg;
			break;
		default:
			fprintf(stderr, "%s\n", "Opção invalida ou faltando argumentos na linha de comando");
			exit(1);
		}
	}

	if ((*arquivoSaida == NULL))
	{
		fprintf(stderr, "%s\n", "Falta valores para determinados argumentos na linha de comando");
		exit(1);
	}
}

/**
 * @brief Função utilizada para inicialização dos vetores alocados dinamicamente.
 * 
 * @param n Dimensão da matriz A.
 * @param i Número de iterações max.
 * @param A Vetor para armazenar Matriz A .
 * @param X Vetor de Solução X (Iniciado com zeros).
 * @param M Vetor com a diagonal da matriz pré-condicionadora.
 * @param B Vetor para armazenar B.
 * @param TxA Vetor para armazenar Matriz transposta de A x A.
 * @param TxB Vetor para armazenar Matriz transposta de A x B.
 * @param vetorDeErros Vetor que irá armazenar vetor o erro de cada iteração.
 * @return Sem retorno
 */
void inicializaVetores(int n, int i, double *X, double *A, double *TxA, double *TxB, int k, double *residuo)
{

	for (int it = 0; it < k * n; ++it)
	{
		A[it] = 0.0;
	}
	for (int linha = 0; linha < n; ++linha)
	{
		residuo[linha] = 0.0;
		TxB[linha] = 0.0;
		X[linha] = 0.0;
	}

	for (int i = 0; i < (k + ((k - 1) * 0.5)) * n; ++i)
	{
		TxA[i] = 0.0;
	}
}

/**
 * @brief Função utilizada para gerar as matrizes A e B
 * 
 * @param n Dimensão da matriz A.
 * @param k Número de diagonais não nulas de A.
 * @param A Vetor para armazenar Matriz A .
 * @param B Vetor para armazenar B.
 * @return Sem retorno
 */
void generateMatrizes(int n, int k, double *A, double *B)
{
	int auxCordCima;  /** auxCordCima < Variavel utilizada para auxilio no preenchimento das k diagonais da matriz A  */
	int auxCordBaixo; /** auxCordBaixo < Variavel utilizada para auxilio no preenchimento das k diagonais da matriz A  */
	int aux;
	int var;
	for (int col = 0; col < n; col++)
	{

		B[col] = generateRandomB(k);

		//printf ("%d de B = %f\n", col,B[col]);
		//printf("%f bb\n",B[col]);

		var = col;
		if (col > (k - 1) / 2)
		{
			var = (k - 1) / 2;
		}
		aux = col * k + var;
		A[aux] = generateRandomA(col, col, k);
		// printf("pos %d %d linha %d coluna\n",aux,col ,col);
		auxCordCima = 1;
		auxCordBaixo = 1;

		while ((col - auxCordCima >= 0) && (auxCordCima < (k + 1) * 0.5))
		{

			A[aux - auxCordCima] = generateRandomA(col - auxCordCima, col, k);
			auxCordCima++;
		}

		while ((col + auxCordBaixo < n) && (auxCordBaixo < +(k + 1) * 0.5))
		{
			A[aux + auxCordBaixo] = generateRandomA(col + auxCordBaixo, col, k);
			auxCordBaixo++;
		}
	}
}

/**
 * @brief Função utilizada para calcular Transposta de A x A e Transposta de A x B.
 * 
 * @param n Dimensão da matriz A.
 * @param TxB Vetor para armazenar Transposta de A x B.
 * @param A Vetor para armazenar Matriz A .
 * @param B Vetor para armazenar B.
 * @param TxA Vetor para armazenar Transposta de A x A.
 */
void calcTransposta(int n, double *TxB, double *A, double *B, double *TxA, int k, int numDiagonais)
{

	//Transposta de A * B , Resultando no vetor TxB //TA FAZENDO CERTO AEHOO
	int aux = 0;
	for (int ite = 0; ite < (k + 1) * 0.5; ite++)
	{
		for (int pos = 0; pos < k; ++pos)
		{
			TxB[ite] += A[aux + pos] * B[pos];
		}
		aux += k;
	}

	int inc = 1;
	for (int ite = (k + 1) * 0.5; ite < n; ite++)
	{
		for (int pos = 0; pos < k; ++pos)
		{
			TxB[ite] += A[aux + pos] * B[pos + inc];
		}
		aux += k;
		inc++;
	}

	int pos = 0;
	int tempCol = 0, tempLin = 0;
	int aux1 = 0;
	int aux2;
	int desloca_lin = 0, desloca_col = 0;
	int j, saveDesloca_lin = 0, saveDesloca_col = 0;

	for (int linha = 0; linha < n; ++linha)
	{ //laco das colunas de TxA
		if (linha >= (k + 1) * 0.5)
		{
			desloca_lin++;
		}
		for (int coluna = linha; coluna < linha + k; ++coluna)
		{ //laco das linhas de TxA
			if (coluna >= (k + 1) * 0.5)
			{
				desloca_col++;
			}

			j = 0;
			if ((desloca_col - desloca_lin == 0) && (desloca_col != 0 || desloca_lin != 0))
			{
				saveDesloca_lin = desloca_lin;
				saveDesloca_col = desloca_col;
				desloca_col = 0;
				desloca_lin = 0;
			}
			if (linha == coluna)
			{
				while (j < k)
				{
					TxA[aux1] += A[linha * k + j] * A[coluna * k + j];
					j++;
				}
			}
			else
			{

				while (j < k - (desloca_col - desloca_lin) && linha < n && coluna < n)
				{
					TxA[aux1] += A[linha * k + j + desloca_col] * A[coluna * k + desloca_lin + j];
					j++;
				}
			}
			aux1++;
		}
		aux1 = (linha + 1) * numDiagonais;
		desloca_col = 0;
	}
	desloca_lin = saveDesloca_lin;
	desloca_col = saveDesloca_col;
}

/**
 * @brief Função utilizada para calcular Matriz Pré-Condicionadora.
 * 
 * @param w Valor indicando qual o pré condicionamento que deve ser utilizado.
 * @param M Vetor com a diagonal da matriz pré-condicionadora.
 * @param n Dimensão da matriz A.
 * @param TxA Vetor para armazenar Matriz transposta de A x A.
 */
void calcMatrizPrecondicionadora(float w, double *M, int n, double *TxA, int numDiagonais)
{

	for (int linha = 0; linha < n; linha++)
	{
		if (w == 0.0)
			M[linha] = 1.0;
		else if (w < 1.0)
			M[linha] = 1.0 / TxA[linha * numDiagonais];
	}
}

/**
 * @brief Função utilizada para inicialização dos vetores alocados dinamicamente.
 * 
 * @param A Vetor para armazenar Matriz A.
 * @param B Vetor para armazenar B.
 * @param X Vetor de Solução X (Iniciado com zeros).
 * @param n Dimensão da matriz A.
 * @param residuo Vetor utilizado para armazenar os residuos.
 * @param valorResiduo Variável que armazena a norma do vetor residuo.
 * @param tempoResiduo Variável utilizada para armazenar o tempo gasto no calculo do residuo.
 */
void calculaResiduo(double *A, double *B, double *X, int n, double *residuo,
					double *valorResiduo, double *tempoResiduo, int k)
{

	*tempoResiduo = timestamp();
	int pos_safe=(k+1) *0.5;

	for (int coluna = 0; coluna < pos_safe; coluna++)
	{
		int tmp = coluna*k;
		#pragma GCC ivdep
		#pragma GCC vector aligned
		#pragma GCC vector always
		for (int p = 0; p < k - k % 4; p+=4)
		{

			residuo[p] += A[tmp + p] * X[coluna];
			residuo[p+1] += A[tmp + p+1] * X[coluna];
			residuo[p+2] += A[tmp + p+2] * X[coluna];
			residuo[p+3] += A[tmp + p+3] * X[coluna];
		}

		for (int p = k- k % 4; p < k ; p++) {
			residuo[p] += A[tmp + p] * X[coluna];
		}
	}

	int auxiliar = 1;

	for (int coluna = pos_safe; coluna < n; coluna++)
	{
		int p = 0;
		int tmp1 = coluna*k;
		while ((p < k) && (p + auxiliar < n))
		{
			residuo[p + auxiliar] += A[tmp1 + p] * X[coluna];

			p++;
		}
		auxiliar++;
	}


	#pragma GCC ivdep
	#pragma vector aligned
	#pragma vector always
	for (int d = 0; d < n; d += 4)
	{
		residuo[d] = B[d] - residuo[d];
		residuo[d + 1] = B[d + 1] - residuo[d + 1];
		residuo[d + 2] = B[d + 2] - residuo[d + 2];
		residuo[d + 3] = B[d + 3] - residuo[d + 3];
	}
	double tmp3 = 0.0;
	for (int d = 0; d < n; d++) {
		tmp3 += residuo[d] *residuo[d];
	}
	*valorResiduo = sqrt(tmp3);

	*tempoResiduo = timestamp() - *tempoResiduo;
}
