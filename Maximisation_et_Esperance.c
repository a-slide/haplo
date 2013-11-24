#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "EM_main.h"
	

/***********************************************************************
 * Maximisation et estimation esperance
 **********************************************************************/

void maximisation(T_info* pvar)
{
/* PRECONDITIONS
 * Le calcul qui va avoir lieu concerne l'étape it
 * Pour l'étape it -1 les fréquences des haplotypes sont connues (freq-prec)
 * Pour l'étape it -1 les probabilité des genotypes sont connues (proba-prec)
 * pvar->nb_ind  est le nbr total d'individus
 * pvar->tab_geno[i].nb_ind est le nbr d'individu possédant le genotype i
 * POSTCONDITIONS
 * Pour l'étape it les fréquences des haplotypes sont connues
 */
	
	int i;
	double freq;
	T_geno_expl* ptrj = NULL; // Pointeur sur la structure contenant les genotypes explicatifs avec l'haplotype complémentaire
	double contrib;

	printf("\n################### MAXIMISATION ###################\n\n");

	for (i = 0; i < pvar->nb_haplo; i++ ) // Parcours du tableau contenant la liste des haplotypes non redondant
	{
		freq = 0; //initialisation fréquence 
		ptrj = pvar->tab_haplo[i].tete;
		printf ("haplo numéro: %d \t sequence: (%s) \t fréquence: %.2e \n",
			i,	pvar->tab_haplo[i].sequence, pvar->tab_haplo[i].frequence_prec);
		printf ("Initialisation de la fréquence:\n");
		printf ("freq_initiale=%f\n", freq);
		while (ptrj != NULL)
		{
			printf("\nGenotype expliqué # %d (%s), possédé par %d individus\n",
				ptrj->num_geno_expl,
				pvar->tab_geno[ptrj->num_geno_expl].sequence,
				pvar->tab_geno[ptrj->num_geno_expl].nb_ind);
			printf("Haplotype_complémentaire numéro# %d dont la séquence est (%s)\n", 
				ptrj -> num_haplo_compl,
				pvar->tab_haplo[ptrj->num_haplo_compl].sequence);
			printf("Nombre d'indiv total:%d\n", pvar->nb_ind);
			
			if (i==ptrj -> num_haplo_compl) //cas homozygote
			{
				// Calcul de la contribution
				contrib=((2 * pvar->tab_haplo[i].frequence_prec * pvar->tab_haplo[i].frequence_prec) / 
				(pvar->tab_geno[ptrj->num_geno_expl].proba_prec)) * ((pvar->tab_geno[ptrj->num_geno_expl].nb_ind / (double)pvar->nb_ind));
				printf("contribution=%.2e\n", contrib);
			}
			else //cas hétérozygote
			{
				// Calcul de la contribution
				contrib=((2 * pvar->tab_haplo[i].frequence_prec * pvar->tab_haplo[ptrj->num_haplo_compl].frequence_prec) / 
				(pvar->tab_geno[ptrj->num_geno_expl].proba_prec)) * ((pvar->tab_geno[ptrj->num_geno_expl].nb_ind / (double)pvar->nb_ind));

				//printf("freq = %.2e\n", freq);
				printf("Contribution = %.2e\n", contrib);
			}
			freq=freq+contrib; //
			printf("freq_accumulée = %.2e\n", freq);
			ptrj = ptrj -> suivant;
		}
		//freq = freq/2.0;

		printf("\nfreq_totale=%.2e\n", freq);
		pvar->tab_haplo[i].frequence = freq/2.0;
		
		printf("nouvelle fréquence de l'haplotype %d:\t%.2e\n\n", i, pvar->tab_haplo[i].frequence);
	}
return;
}


double estimation_esperance(T_info* pvar)
{
/* PRECONDITIONS
 * Le calcul qui va avoir lieu concerne l'itération courante
 * Les fréquences des haplotypes sont connues pour l'itération courante
 * nb_ind est le nombre d'individus pour le génotype geno
 * POSTCONDITIONS
 * Les probabilités des génotypes sont connues à l'itération courante 
 * La vraisemblance des données a été calculée.
 */
	
	int i;
	double f1, f2, loglikelihood;
	T_diplo_expl* ptrj = NULL;

	loglikelihood = 0;

	printf("\n############### ESTIMATION ESPERANCE ###############\n\n");

	// Pour chaque genotypes
	for (i = 0; i < pvar->nb_geno; i++ )
	{
		pvar->tab_geno[i].proba = 0;
		ptrj = pvar->tab_geno[i].tete;
		
		printf ("\nProbabilites partielles de geno # %d", i);
		// Pour chaque diplotypes explicatifs
		while (ptrj != NULL) 
		{
			if (ptrj->num_haplo_A == ptrj->num_haplo_B) // Si Homozygote
			{
				f1 = pvar->tab_haplo [ptrj->num_haplo_A].frequence; // frequence de l'haplotype A
				pvar->tab_geno[i].proba += (f1*f1); // mise à jour de la probabilité du genome i
			}
			else // Si Heterozygote
			{	
				f1 = pvar->tab_haplo [ptrj->num_haplo_A].frequence; // frequence de l'haplotype A
				f2 = pvar->tab_haplo [ptrj->num_haplo_B].frequence; // frequence de l'haplotype A
				pvar->tab_geno[i].proba += (2*f1*f2); //  mise à jour de la probabilité du genome i
			}
			printf (" -> %.2e", pvar->tab_geno[i].proba);
			ptrj = ptrj -> suivant;
		}
		printf ("\nProbabilite finale de geno %d = %.2e\n", i, pvar->tab_geno[i].proba);
		
		// #####################################################
		// problème ici...
		
		// Calcul de la loglikelyhood par acumulation logarithmique de vraisemblance
		loglikelihood = loglikelihood + (pvar->tab_geno[i].nb_ind) * log (pvar->tab_geno[i].proba); // A revoir car pas sûr
		printf("LogLikelihood calculé lors de ce génotype: %.2e\n\n", loglikelihood);
	}
	return loglikelihood;
}


void update_proba_freq_vraisemblance (T_info* pvar)
{
	int i;

	printf("\n################ UPDATE FREQ ET PROBA ##############\n\n");

	// Mise à jour des probabilités des génotypes proba.prec qui prennent la valeur des proba courantes
	for (i = 0; i < pvar->nb_geno; i++ )
	{
		printf("Genotype #%d\tAncienne proba_prec : %.2e", i, pvar->tab_geno[i].proba_prec);
		pvar->tab_geno[i].proba_prec = pvar->tab_geno[i].proba;
		printf("\tNouvelle proba_prec : %.2e\n", pvar->tab_geno[i].proba_prec);
	}
	printf("\n");
	
	// Mise à jour des fréquences des haplotypes frequence.prec qui prennent la valeur des fréquences courantes
	for (i = 0; i < pvar->nb_haplo; i++ )
	{
		printf("Haplotype #%d\tAncienne freq_prec : %.2e", i, pvar->tab_haplo[i].frequence_prec);
		pvar->tab_haplo[i].frequence_prec = pvar->tab_haplo[i].frequence;
		printf("\tNouvelle freq_prec : %.2e\n", pvar->tab_haplo[i].frequence_prec);
	}
	printf("\n");
	
	// Mise à jour de la vraisemblance précédente qui prend la valeur de la vraisemblance courante
	printf("Ancienne vraisemblance : %.2e\t ", pvar->vraisemblance_prec);
	pvar->vraisemblance_prec = pvar->vraisemblance;
	printf("\tNouvelle vraisemblance : %.2e\n\n ", pvar->vraisemblance_prec);

	return;
}




