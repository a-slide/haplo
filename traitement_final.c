#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "inference_haplotype.h"
#include "ptr_allocation.h"

/***********************************************************************
 * diplotype_plus_probable
 **********************************************************************/

void diplotype_plus_probable (T_info* pvar)
{
/* PRECONDITIONS
 * Les fréquences des haplotypes ont été calculées jusqu'à l'itération de sortie
 * POSTCONDITIONS
 * Le diplotype explicatif le plus probable a été associé à chaque genotype
 */
	
	int i;
	double proba_diplo;
	T_diplo_expl* ptrj = NULL;

	for (i = 0; i < pvar->nb_geno; i++ ) // Pour chaque génotype de tab_geno
	{
		ptrj = pvar->tab_geno[i].tete;
		pvar->tab_geno[i].num_haplo_A_max = 0; // initialisation
		pvar->tab_geno[i].num_haplo_B_max = 0;
		pvar->tab_geno[i].proba_diplo_max = 0;
		printf ("GENOTYPE %d\n", i);
		printf ("Affichage uniquement en cas de probabilité plus haute\n");
		
		while (ptrj != NULL) // Pour chaque diplotype explicatif
		{
			// Si on a affaire à un homozygote
			if (ptrj->num_haplo_A == ptrj->num_haplo_B)
 				proba_diplo = pvar->tab_haplo[ptrj->num_haplo_A].frequence * pvar->tab_haplo[ptrj->num_haplo_A].frequence;
 				
			// Si on a affaire à un hétérozygote
			else
				proba_diplo = 2*(pvar->tab_haplo[ptrj->num_haplo_A].frequence * pvar->tab_haplo[ptrj->num_haplo_B].frequence);

			// Si le diplotype est plus fréquent que les diplotypes précédents.
			if (proba_diplo > pvar->tab_geno[i].proba_diplo_max)
			{
				printf ("Diplotype : %d / %d\tProba : %e\n",ptrj->num_haplo_A, ptrj->num_haplo_B, proba_diplo);
				pvar->tab_geno[i].proba_diplo_max = proba_diplo; // le diplotype les plus explicatif devient le diplotype courant
				pvar->tab_geno[i].num_haplo_A_max = ptrj->num_haplo_A;// Les numéros des 2 haplotypes correspondant sont mis à jour
				pvar->tab_geno[i].num_haplo_B_max = ptrj->num_haplo_B;	
			}
			ptrj = ptrj -> suivant;
		}
	}
	return;
}


/***********************************************************************
 * export_geno_diplo
 **********************************************************************/

void export_geno_diplo (T_info* pvar)
{
	int i;
	int num_geno, num_haplo_A, num_haplo_B;
	FILE* fp = init_file_ptr ("Liste_Diplo_Expl.txt", "w");
	
	for (i = 0; i < pvar->nb_ind; i++ ) // Pour chaque génotype de tab_geno
	{
		// Accès à tab_individus et mémorisation du numéro de génotype associé 
		fprintf(fp, "Individu # %d\n",i);
		num_geno = pvar->tab_individus[i].num_geno;
		
		// Accès au génotype associé dans tab_geno et mémorisation des numéros des haplotypes explicatifs 
		fprintf(fp, "Genotype : %d\t%s\n", num_geno, pvar->tab_individus[i].sequence);
		num_haplo_A = pvar->tab_geno[num_geno].num_haplo_A_max;
		num_haplo_B = pvar->tab_geno[num_geno].num_haplo_B_max;
		
		// Accès aux séquences des haplotypes explicatifs dans tab_haplo 
		fprintf(fp, "HaploA %d :\t%s\n", num_haplo_A, pvar->tab_haplo[num_haplo_A].sequence);
		fprintf(fp, "HaploB %d :\t%s\n\n", num_haplo_B, pvar->tab_haplo[num_haplo_B].sequence);
	}
	fclose(fp);
	return;
}


/***********************************************************************
 * comparaison_frequence = fonction utilisateur de comparaison pour qsort
 **********************************************************************/

int comparaison_frequence (void const* a, void const* b)
{
	// typage des pointeurs
	T_haplo const* pa = a;
	T_haplo const* pb = b;
	double difference = 0;
	
	// Calcul de la différence pour tri decroissant
	difference = (pb->frequence - pa->frequence); 
	
	// test multiple car la notion d'égalité n'existe pas en réel
	if (difference > 0)
		return 1;
	
	else if (difference < 0)
		return -1;
	
	else
		return 0;
}

/***********************************************************************
 * export_haplo
 **********************************************************************/

void export_haplo (T_info* pvar)
{
	int i;
	FILE* fp = init_file_ptr ("Liste_Haplo_Freq.txt", "w"); // Ouverture d'un fichier en mode écriture

	for (i = 0; i < pvar->nb_haplo; i++ ) // Ecriture des séquences des haplotypes suivis de leur fréquence
		fprintf(fp, "%s \t Fréquence %e\n", pvar->tab_haplo[i].sequence, pvar->tab_haplo[i].frequence);
	
	fclose(fp);
	return;
}
