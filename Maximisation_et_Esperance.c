#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "EM_main.h"
	

/***********************************************************************
 * Maximisation et estimation esperance
 **********************************************************************/
void Maximisation_et_Esperance(T_info* pvar)
{

Maximisation (pvar);
pvar->vraisemblance=Estimation_Esperance (pvar);


return;
}

void Maximisation(T_info* pvar)
{
int i;
double freq;
T_geno_expl* ptrj = NULL; // Pointeur sur la structure contenant les genotypes explicatifs avec l'haplotype complémentaire
double contrib;

	for (i = 0; i < pvar->nb_haplo; i++ ) // Parcours du tableau contenant la liste des haplotypes non redondant
	{
		freq=0; //initialisation fréquence 
		ptrj = pvar->tab_haplo[i].tete;
		printf ("haplo numéro: %d \t sequence: (%s) \t fréquence: %.2e \n", i, 	pvar->tab_haplo[i].sequence, 
											pvar->tab_haplo[i].frequence_prec);
		printf ("Initialisation de la fréquence:\n");
		printf ("freq_initiale=%f\n", freq);
		while (ptrj != NULL)
		{
			printf("\nGenotype expliqué # %d (%s), possédé par %d individus \nHaplotype_complémentaire numéro#%i dont la séquence est (%s)\n", 
				ptrj -> num_geno_expl,
				pvar->tab_geno[ptrj->num_geno_expl].sequence,
				pvar->tab_geno[ptrj->num_geno_expl].nb_ind,
				ptrj -> num_haplo_compl,
				pvar->tab_haplo[ptrj->num_haplo_compl].sequence);
				
				printf("nombre d'indiv total:%d\n", pvar->nb_ind);
			
			if (i==ptrj -> num_haplo_compl) //cas homozygothie
			{
			// Calcul de la contribution
			contrib=((2*pvar->tab_haplo[i].frequence_prec*pvar->tab_haplo[i].frequence_prec)/
			(pvar->tab_geno[ptrj->num_geno_expl].proba_prec))*(pvar->tab_geno[ptrj->num_geno_expl].nb_ind/(float)pvar->nb_ind);

			
			printf("contribution=%.2e\n", contrib);

			}
			else //cas hétérozygothie
			{
			// Calcul de la contribution
			contrib=((2*pvar->tab_haplo[i].frequence_prec*pvar->tab_haplo[ptrj -> num_haplo_compl].frequence_prec)/
			(pvar->tab_geno[ptrj->num_geno_expl].proba_prec))*(pvar->tab_geno[ptrj->num_geno_expl].nb_ind/(float)pvar->nb_ind);

			//printf("freq=%.2e\n", freq);
			printf("contribution=%.2e\n", contrib);
			}
			freq=freq+contrib; //
			printf("freq_accumulée=%.2e\n", freq);
			ptrj = ptrj -> suivant;
		}
		//freq=freq/2.0;

		printf("\nfreq_totale=%.2e\n", freq);
		pvar->tab_haplo[i].frequence=freq/2.0;
		
		printf("nouvelle fréquence de l'haplotype %d:\t%.2e\n", i, pvar->tab_haplo[i].frequence);
		printf("########################################################################################################");
		printf("\n");
	}
return;
}


// Le calcul qui va avoir lieu concerne l'étape it
// Pour l'étape it -1 les fréquences des haplotypes sont connues (freq-prec)
// Pour l'étape it -1 les probabilité des genotypes sont connues (proba-prec)
// N est le nbr total d'individus
// Ngeno est le nbr d'individu possédant le genotype geno
// POSTCONDITIONS
// Pour l'étape it les fréquences des haplotypes sont connues
//
//	pour chaque haplotype h1
//	//	récupérer freq_prec (h1)
//	//	freq = 0
//	//	pour chaque genotype geno susceptible d'être expliqué par h1
//	//	//	si geno (h1,h1)
//	//	//	//	alors	contribution = ( 2 freq_prec (h1)2 / proba_prec ( geno )) . ( Ngeno / N ))
//	//	//	//	sinon	contribution = ( 2 freq_prec (h1) . freq_prec (h2)   / proba_prec ( geno ) . ( Ngeno / N ))
//	//	//	finsi
//	//	//	freq = freq + contribution
//	//	finpour
//	//	freq = ½ freq
//	finpour


double Estimation_Esperance(T_info* pvar)
{

int i;
double f1, f2, loglikelihood;
T_diplo_expl* ptrj = NULL;
loglikelihood = 0;

for (i = 0; i < pvar->nb_geno; i++ )
	{
		pvar->tab_geno[i].proba = 0;
		ptrj = pvar->tab_geno[i].tete;
		
		printf ("\nProbabilites partielles de geno %d", i);
		while (ptrj != NULL) // parcourt la liste de diplo explicatifs jusqu'à la fin
		{
			if (ptrj->num_haplo_A == ptrj->num_haplo_B) // Si Haplo A et Haplo B sont les mêmes
			{
				f1 = pvar->tab_haplo [ptrj->num_haplo_A].frequence; // frequence de l'haplotype A
				pvar->tab_geno[i].proba += (f1*f1); // mise à jour de la probabilité du genome i
			}
			else
			{	
				f1 = pvar->tab_haplo [ptrj->num_haplo_A].frequence; // frequence de l'haplotype A
				f2 = pvar->tab_haplo [ptrj->num_haplo_B].frequence; // frequence de l'haplotype A
				pvar->tab_geno[i].proba += (2*f1*f2); //  mise à jour de la probabilité du genome i A VERIFIER
			}
			
			printf (" -> %.2e", pvar->tab_geno[i].proba);
			ptrj = ptrj -> suivant;
		}
		printf ("\nProbabilite finale de geno %d = %.2e\n", i, pvar->tab_geno[i].proba);
		loglikelihood = loglikelihood + (pvar->tab_geno[i].nb_ind) * log(pvar->tab_geno[i].proba); // A revoir car pas sûr
		printf("loglikelihood calculé lors de ce génotype: %.2e\n\n", loglikelihood);
	}
	return loglikelihood;
}

//Précondition
//Le calcul qui va avoir lieu concerne l'étape it
//Les fréquences des haplotypes sont connues pour l'étape it-1
//Ngeno est le nombre d'individus pour le génotype geno
//Postcondition
//Les probabilités des génotypes sont connues à l'étape it 
//La vraisemblance des données a été calculée.

//	loglikelyhood = 0 // accumulateur logarithmique de vraisemblance
//	pour chaque genotype geno
//	//	proba = 0
//	//	pour chaque paires explicative (h1, h2) de geno
//	//	//	p1 = freq (h1)
//	//	//	p2 = freq (h2)
//	//	//	si (h1 == h2)
//	//	//	//	alors	ppart = p12
//	//	//	//	sinon	ppart = 2.p1.p2
//	//	//	finsi
//	//	//	p = p+ppart	
//	//	finpour
		//postcondition: La probabilités du génotype geno vient d'être calculée	
//	//	loglikelyhood =  loglikelyhood + Ngeno log (p) // accumule 
//	finpour


void Update_Hfreqpreq_Gprobaprec_vraisemblance_preq (T_info* pvar)
{
	int i;
	printf("########################################################################################\n");
	printf("######################## Mise à jour de freq_prec et proba_prec ########################\n");
	printf("########################################################################################\n\n");


	for (i = 0; i < pvar->nb_geno; i++ )
	{
		pvar->tab_geno[i].proba_prec = pvar->tab_geno[i].proba;
		printf("nouvelle proba_prec du génotype #%d: %.2e\n", i, pvar->tab_geno[i].proba_prec);
	}
	
	printf("\n\n");

	for (i = 0; i < pvar->nb_haplo; i++ )
	{
		pvar->tab_haplo[i].frequence_prec = pvar->tab_haplo[i].frequence;
		printf("nouvelle freq_prec de l'haplotype #%d: %.2e\n", i, pvar->tab_haplo[i].frequence_prec);
	}
	pvar->vraisemblance_prec=pvar->vraisemblance; // Mise à jour de la vraisemblance précédente qui prend la valeur de la vraisemblance 								courante
	return;
}




