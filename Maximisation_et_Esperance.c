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
			(pvar->tab_geno[ptrj->num_geno_expl].proba))*(pvar->tab_geno[ptrj->num_geno_expl].nb_ind/(float)pvar->nb_ind);

			
			printf("contribution=%.2e\n", contrib);

			}
			else //cas hétérozygothie
			{
			//Grosse question: la fréquence de la séquence complémentaire, doit elle être la même que la séquence courante ?(vu qu'avec 				l'option -E toutes les fréquences sont les mêmes) ou bien on autorise la fréquence complémentaire a avoir été mise à 				jour ?	(vu que la fréquence de l'haplotype courant peut se trouver être une séquence complémentaire d'un autre haplotype 				plus tard dans la boucle for.
			//Pour l'instant je considère la première proposition mais dans ce cas je suis obligé de faire l'égalité suivante:
			

			// Calcul de la contribution
			contrib=((2*pvar->tab_haplo[i].frequence_prec*pvar->tab_haplo[ptrj -> num_haplo_compl].frequence_prec)/
			(pvar->tab_geno[ptrj->num_geno_expl].proba))*(pvar->tab_geno[ptrj->num_geno_expl].nb_ind/(float)pvar->nb_ind);

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


