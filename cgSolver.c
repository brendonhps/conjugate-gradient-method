/**
 * @file cgSolver.c
 * @author Brendon Henrique de Paula da Silva GRR20170203
 * @date 16 Setembro 2018
 * @brief Função Main
 * 
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include "gradienteconj.h"
#include "utils.h"
#include <getopt.h>
#include <likwid.h>



/**
 * @brief Função Main.
 * 		  Resposável por receber as opções de linha de comando.
 * 		  Abrir arquivo de saída passado por linha de comando.
 * 		  Chamar as funções auxiliares. 
 * 		  Escrever no arquivo os dados calculados.	
 * @param argc Número de argumentos passado para o programa.
 * @param argv Vetor de string contendo os argumentos.
 * @return Retorna 0 caso o programa funcione e 1 caso se tenha algum erro.
 */

int main(int argc, char *argv[])
{
 LIKWID_MARKER_INIT; 

	int opcao=0, i = 0, totalIt;
	/** opcao < Opção passada pelo usuário do programa*/
	/** totalIt < Total de Iterações realizadas pelo método.*/
	int n = 0, k = 0;
	/** n < Dimensão da matriz A */
	/** k < Número de diagonais não nulas de A*/
	double EPS = 0.0, *A, *TxA, *TxB, *B, *X, *M,*vetorDeErro,*residuo,tempoPC,tempoResiduo,tempoTotalIt,valorResiduo=0.0;
	/** EPS < Erro aproximado absoluto máximo passado pelo usuário.*/
	/** A < Vetor de tamanho n*n.*/
	/** TxA < Vetor que armazena a transposta de A*A.*/
	/** TxB < Vetor que armazena a transposta de A*B.*/
	/** B < Vetor utilizado para representar a coluna B.*/
	/** X < Vetor de solução do sistema Ax=B , inicializado com 0.*/
	/** M < Vetor utilizado para guardar a inversa da matriz pré-condicionante.*/
	/** vetorDeErro < Vetor utilizado para armazenar o erro aproximado em x após cada iteração.*/
	/** residuo < Vetor com os residuos. */
	/** tempoPC < Tempo para calcular a inversa da matriz pré-condicionante M e preparar o SL para o uso do pré-condicionante. */
	/** tempoResiduo < Tempo para calcular o residuo. */
	/** tempoTotalIt < Tempo total de todas as iterações.*/
	/** valorResiduo < Utilizado para armazenar a norma do vetor de residuos.*/
	float w = -1.0;
	/** w < Variável que indica qual matriz pré-condicionadora deve ser utilizada. */

	char *arquivoSaida = NULL;
	/** arquivoSaida < Armazena a string com o nome do arquivo de saida passado pelo usuário*/
	entradaDados(opcao, argc, argv, &n, &k, &w, &i, &EPS, &arquivoSaida);


	FILE *output;

	output = fopen(arquivoSaida, "w");

	if (!output)
	{
		fprintf(stderr, "%s\n", "Erro na abertura do arquivo de saida\n");
		
		exit(1);
	}

	srand(20182);

	B = (double *)malloc(sizeof(double) * n);
	A = (double *)malloc(sizeof(double)* k*n);
	
	TxB = (double *)malloc(sizeof(double) *n);
	TxA = (double *)malloc(sizeof(double) *(k+((k-1)/2))*n);

	M = (double *)malloc(sizeof(double) * n);
	X = (double *)malloc(sizeof(double) * n);

	residuo = (double *)malloc(sizeof(double) * n);
	
	inicializaVetores(n,i,X,A,TxA,TxB,k,residuo);  
	generateMatrizes(n, k, A, B); 
    


	tempoPC = timestamp();

	int numDiagonais = k+((k-1)/2);
	calcTransposta(n, TxB, A, B, TxA, k,numDiagonais);


	calcMatrizPrecondicionadora(w, M, n, TxA,numDiagonais); 

	tempoPC = timestamp() - tempoPC;

	//Método do gradiente conjulgado A.k.a Operação 1
	LIKWID_MARKER_START("CG");
	
	vetorDeErro = gradienteConjPreCond(TxA, TxB, n, X, M, i, &tempoTotalIt, &totalIt, EPS,numDiagonais); 
	
	LIKWID_MARKER_STOP("CG");
	free(TxB);
	free(TxA);
	free(M);

	//Calculo do Residuo A.k.a Operação 2
	
	LIKWID_MARKER_START("RESIDUO");

	calculaResiduo(A, B, X, n, residuo, &valorResiduo, &tempoResiduo,k); 
	
	LIKWID_MARKER_STOP("RESIDUO");
	free(residuo);
	free(A);
	free(X);
	free(B);

	//Imprime no arquivo as respostas
	fprintf(output, "# bhps17 Brendon Henrique\n#\n");
	for (int it = 1; it <= totalIt; it++)
	{
		fprintf(output, "# iter %d: <%.15g>\n", it, fabs(vetorDeErro[it]));
	}
	free(vetorDeErro);
	fprintf(output, "# residuo: <%.15g>\n", valorResiduo);
	fprintf(output, "# Tempo PC: <%.15g> \n", tempoPC/1000);
	fprintf(output, "# Tempo gasto no GC: <%.15g> \n", tempoTotalIt/1000);
	tempoTotalIt = tempoTotalIt/totalIt;
	fprintf(output, "# Tempo iter: <%.15g> \n", tempoTotalIt/1000);
	fprintf(output, "# Tempo residuo: <%.15g> \n#\n", tempoResiduo/1000);
	fprintf(output, "%d\n", n);

	for (int raiz = 0; raiz < n; raiz++)
	{
		fprintf(output, "%.15g\t  ",  X[raiz]);
	}
	fclose(output);
	
	LIKWID_MARKER_CLOSE; 
	return (0);
}
