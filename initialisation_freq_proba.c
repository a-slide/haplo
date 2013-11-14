#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "EM_main.h"

/***** initialisation_freq_proba ******************************************/

void initialisation_freq_proba (T_info* pvar, char mode)
{
	if ((mode == 'A') || (mode == 'a'))
		haplo_random_freq (pvar);
	
	else
		init_haplo_equi_freq (pvar);
	
	print_tab_haplo (pvar);
	
	
	///	init_geno_proba (pvar);
	
	print_tab_geno (pvar);
	
	return;
}

void init_haplo_equi_freq (T_info* pvar)
{
	int i;
	
	printf ("\nInitialisation equi-probable des frequences des haplotypes\n\n");
	for (i = 0; i < pvar->nb_haplo; i++)
		pvar->tab_haplo[i].frequence = 1.0/ pvar->nb_haplo;

	return;
}

void haplo_random_freq (T_info* pvar)
{
	int i;
	int sum = 0;
	srand (time(NULL)); // seed value pour rand
	
	printf ("\nInitialisation aleatoire des frequences des haplotypes\n\n");
	for (i = 0; i < pvar->nb_haplo; i++)
	{
		pvar->tab_haplo[i].frequence = rand() % 100 + 1;
		sum = sum + pvar->tab_haplo[i].frequence;
	}

	for (i = 0; i < pvar->nb_haplo; i++)
		pvar->tab_haplo[i].frequence = pvar->tab_haplo[i].frequence / sum;

	return;
}


 
/**** print_tab_haplo *************************************************/
void print_tab_haplo (T_info* pvar)
{
	int i;
	T_geno_expl* ptrj = NULL;
	
	printf("\n\n###############################################\n");
	printf("\nListe Des haplotypes et des genotypes expliqués\n");
	printf("\n###############################################\n\n");
	
	for (i = 0; i < pvar->nb_haplo; i++ )
	{
		ptrj = pvar->tab_haplo[i].tete;
		printf("Haplotype #%d \t Séquence : %s \t Fréquence %f \t Nombre de genotype(s) expliqué(s) : %d \n",
			i,
			pvar->tab_haplo[i].sequence,
			pvar->tab_haplo[i].frequence,
			pvar->tab_haplo[i].nb_geno_expl);
			
		printf("Liste des génotypes expliqués\n");
		while (ptrj)
		{
			printf("Genotype # %d (%s) avec Haplotype # %d (%s)\n", 
				ptrj -> num_geno_expl,
				pvar->tab_geno[ptrj->num_geno_expl].sequence,
				ptrj -> num_haplo_compl,
				pvar->tab_haplo[ptrj->num_haplo_compl].sequence );
						
			ptrj = ptrj -> suivant;
		}
		printf("\n");
	}
}

/**** print_tab_geno *************************************************/
void print_tab_geno (T_info* pvar)
{
	int i;
	T_diplo_expl* ptrj = NULL;
	
	printf("\n\n#################################################\n");
	printf("\nListe des genotypes et des diplotypes explicatifs\n");
	printf("\n#################################################\n\n");
	
	for (i = 0; i < pvar->nb_geno; i++ )
	{
		ptrj = pvar->tab_geno[i].tete;
		printf("Genotype #%d\t Séquence : %s\t Probabilité %f\t Nombre d'individu(s) concerné(s) : %d\t Nombre de diplotype(s) explicatif(s) : %d\n",
			i,
			pvar->tab_geno[i].sequence,
			pvar->tab_geno[i].proba,
			pvar->tab_geno[i].nb_ind,
			pvar->tab_geno[i].nb_diplo_expl);
				
		printf("Liste des diplotypes explicatifs\n");
		while (ptrj)
		{
			printf("Haplotype # %d (%s) avec Haplotype # %d (%s)\n", 
				ptrj->num_haplo_A,
				pvar->tab_haplo[ptrj->num_haplo_A].sequence,
				ptrj->num_haplo_B,
				pvar->tab_haplo[ptrj->num_haplo_A].sequence );
						
			ptrj = ptrj -> suivant;
		}
		printf("\n");
	}
}
