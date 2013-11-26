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

	for (i = 0; i < pvar->nb_haplo; i++ ) // Pour chaque haplotype de tab_haplo
	{
		freq = 0; //initialisation fréquence 
		ptrj = pvar->tab_haplo[i].tete;
		printf ("Haplotype #%d\t Sequence: %s\tFrequence Initiale: %e", i, pvar->tab_haplo[i].sequence, pvar->tab_haplo[i].frequence_prec);

		while (ptrj != NULL) // pour chaque genotype dont l'haplotype courant rentre dans l'explication
		{
			if (i==ptrj -> num_haplo_compl) //cas homozygote
			{
				// Calcul de la contribution
				contrib = ((2 * pvar->tab_haplo[i].frequence_prec * pvar->tab_haplo[i].frequence_prec) / (pvar->tab_geno[ptrj->num_geno_expl].proba_prec)) * ((pvar->tab_geno[ptrj->num_geno_expl].nb_ind / (double)pvar->nb_ind));
				/// printf("contribution=%.2e\n", contrib);
			}
			else //cas hétérozygote
			{
				// Calcul de la contribution
				contrib = ((2 * pvar->tab_haplo[i].frequence_prec * pvar->tab_haplo[ptrj->num_haplo_compl].frequence_prec) / (pvar->tab_geno[ptrj->num_geno_expl].proba_prec)) * ((pvar->tab_geno[ptrj->num_geno_expl].nb_ind / (double)pvar->nb_ind));

				//printf("freq = %.2e\n", freq);
				/// printf("Contribution = %.2e\n", contrib);
			}
			// Accumulation de contribution pour chaque genotype expliqués
			freq += contrib; 
			// Progression dans la liste
			ptrj = ptrj -> suivant;
		}
		pvar->tab_haplo[i].frequence = freq/2.0;
		printf("\tFrequence Finale: %e\n", pvar->tab_haplo[i].frequence);
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

	for (i = 0; i < pvar->nb_geno; i++ ) // Pour chaque genotypes de tab_geno
	{
		pvar->tab_geno[i].proba = 0;
		ptrj = pvar->tab_geno[i].tete;

		printf ("Genotype #%d\t Sequence: %s\tProbabilité Initiale: %e", i, pvar->tab_geno[i].sequence, pvar->tab_geno[i].proba_prec);

		while (ptrj != NULL) // Pour chaque diplotypes explicatifs
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
			ptrj = ptrj -> suivant;
		}
		printf("\tProbabilite Finale: %e\n", pvar->tab_geno[i].proba);
			
		// Calcul de la loglikelyhood par acumulation logarithmique de vraisemblance
		loglikelihood += ((pvar->tab_geno[i].nb_ind) * log (pvar->tab_geno[i].proba));
	}
	return loglikelihood;
}


void update_proba_freq (T_info* pvar)
{
	int i;

	printf("\n################ UPDATE FREQ ET PROBA ##############\n\n");

	// Mise à jour des probabilités des génotypes proba.prec qui prennent la valeur des proba courantes
	for (i = 0; i < pvar->nb_geno; i++ )
		pvar->tab_geno[i].proba_prec = pvar->tab_geno[i].proba;
	
	
	// Mise à jour des fréquences des haplotypes frequence.prec qui prennent la valeur des fréquences courantes
	for (i = 0; i < pvar->nb_haplo; i++ )
		pvar->tab_haplo[i].frequence_prec = pvar->tab_haplo[i].frequence;
	
	return;
}




