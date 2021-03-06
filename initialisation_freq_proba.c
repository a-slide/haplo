#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "inference_haplotype.h"

/***********************************************************************
 * initialisation_freq_proba
 **********************************************************************/

void initialisation_freq_proba (T_info* pvar, int mode)
{
	printf ("\nINITIALISATION DES FREQUENCES DES HAPLOTYPES\n");
	// Initialisation équiprobable ou aléatoire des fréquences des haplotypes
	(mode == 1) ? haplo_random_freq (pvar) : init_haplo_equi_freq (pvar);
	print_tab_haplo (pvar);
	
	printf ("\nINITIALISATION DES PROBABILITES DES GENOTYPES\n");
	// Calcul des probabilités des génotypes
	init_geno_proba (pvar);
	print_tab_geno (pvar);
	
	return;
}

/***********************************************************************
 * init_haplo_equi_freq
 **********************************************************************/
 
void init_haplo_equi_freq (T_info* pvar)
{
	int i;
	
	for (i = 0; i < pvar->nb_haplo; i++)
		pvar->tab_haplo[i].frequence_prec = 1.0/ pvar->nb_haplo;

	return;
}

/***********************************************************************
 * haplo_random_freq
 **********************************************************************/
 
void haplo_random_freq (T_info* pvar)
{
	int i;
	int sum = 0;
	srand (time(NULL)); // seed value pour rand
	
	for (i = 0; i < pvar->nb_haplo; i++)
	{
		pvar->tab_haplo[i].frequence_prec = rand() % 100 + 1; // initialisation aléatoire entre 1 et 100
		sum = sum + pvar->tab_haplo[i].frequence_prec;
	}

	for (i = 0; i < pvar->nb_haplo; i++)
		pvar->tab_haplo[i].frequence_prec = pvar->tab_haplo[i].frequence_prec / sum; // transformation du pourcentage en fréquence

	return;
}

/***********************************************************************
 * init_geno_proba
 **********************************************************************/
 
void init_geno_proba (T_info* pvar)
{
	int i;
	double f1, f2;
	T_diplo_expl* ptrj = NULL;
	
	for (i = 0; i < pvar->nb_geno; i++ ) // pour chaque génotype de tab_geno
	{
		ptrj = pvar->tab_geno[i].tete;
		
		while (ptrj != NULL) // parcours la liste de diplo explicatifs jusqu'à la fin
		{
			if (ptrj->num_haplo_A == ptrj->num_haplo_B) // Si Haplo A et Haplo B sont les mêmes
			{
				f1 = pvar->tab_haplo [ptrj->num_haplo_A].frequence_prec; // frequence de l'haplotype A
				pvar->tab_geno[i].proba_prec += (f1*f1); // mise à jour de la probabilité du genome i
			}
			else
			{	
				f1 = pvar->tab_haplo [ptrj->num_haplo_A].frequence_prec; // frequence de l'haplotype A
				f2 = pvar->tab_haplo [ptrj->num_haplo_B].frequence_prec; // frequence de l'haplotype A
				pvar->tab_geno[i].proba_prec += (2*f1*f2); //  mise à jour de la probabilité du genome i
			}
			ptrj = ptrj -> suivant;
		}
		//printf ("\nProbabilite initiale de geno %d = %e\n", i, pvar->tab_geno[i].proba_prec);
	}
	return;
}


/***********************************************************************
 * print_tab_haplo
 **********************************************************************/

void print_tab_haplo (T_info* pvar)
{
	int i;
	T_geno_expl* ptrj = NULL;
	
	printf("\n############### LISTE DES HAPLOTYPES ###############\n\n");
	
	for (i = 0; i < pvar->nb_haplo; i++ )
	{
		ptrj = pvar->tab_haplo[i].tete;
		printf("Haplotype #%d \t Séquence : %s \t Fréquence %.2e \t Nombre de genotype(s) expliqué(s) : %d \n",
			i,
			pvar->tab_haplo[i].sequence,
			pvar->tab_haplo[i].frequence_prec,
			pvar->tab_haplo[i].nb_geno_expl);
			
		printf("Liste des génotypes expliqués\n");
		while (ptrj != NULL)
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

/***********************************************************************
 * print_tab_geno
 **********************************************************************/
 
void print_tab_geno (T_info* pvar)
{
	int i;
	T_diplo_expl* ptrj = NULL;
	
	printf("\n############### LISTE DES GENOTYPES ################\n\n");
	
	for (i = 0; i < pvar->nb_geno; i++ )
	{
		ptrj = pvar->tab_geno[i].tete;
		printf("Genotype #%d\t Séquence : %s\t Probabilité %.2e\t Nombre d'individu(s) concerné(s) : %d\t Nombre de diplotype(s) explicatif(s) : %d\n",
			i,
			pvar->tab_geno[i].sequence,
			pvar->tab_geno[i].proba_prec,
			pvar->tab_geno[i].nb_ind,
			pvar->tab_geno[i].nb_diplo_expl);
				
		printf("Liste des diplotypes explicatifs\n");
		while (ptrj != NULL)
		{
			printf("Haplotype # %d (%s) avec Haplotype # %d (%s)\n", 
				ptrj->num_haplo_A,
				pvar->tab_haplo[ptrj->num_haplo_A].sequence,
				ptrj->num_haplo_B,
				pvar->tab_haplo[ptrj->num_haplo_B].sequence );
						
			ptrj = ptrj -> suivant;
		}
		printf("\n");
	}
}
