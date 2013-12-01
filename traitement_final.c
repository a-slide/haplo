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
 * Les fréquences des haplotypes ont été calculés jusqu'à l'itération de sortie
 * POSTCONDITIONS
 * Un diplotype explicatif le plus probable a été associé à chaque genotype
 */
	
	int i;
	double proba_diplo;
	T_diplo_expl* ptrj = NULL;

	for (i = 0; i < pvar->nb_geno; i++ ) // Pour chaque genotypes de tab_geno
	{
		ptrj = pvar->tab_geno[i].tete;
		pvar->tab_geno[i].num_haplo_A_max = 0;
		pvar->tab_geno[i].num_haplo_B_max = 0;
		pvar->tab_geno[i].proba_diplo_max = 0;
		printf ("GENOTYPE %d\n", i);
		
		while (ptrj != NULL) // Pour chaque diplotypes explicatifs
		{
			// Si Homozygote
			if (ptrj->num_haplo_A == ptrj->num_haplo_B)
			{
 				proba_diplo = pvar->tab_haplo[ptrj->num_haplo_A].frequence * pvar->tab_haplo[ptrj->num_haplo_A].frequence;
 				printf ("Diplotype : %d / %d\tProba : %e\n",ptrj->num_haplo_A, ptrj->num_haplo_B, proba_diplo);
 						
			}
			// Si Heterozygote
			else
			{	
				proba_diplo = 2*(pvar->tab_haplo[ptrj->num_haplo_A].frequence * pvar->tab_haplo[ptrj->num_haplo_B].frequence);
 				printf ("\nDiplotype : %d / %d\tProba : %e",ptrj->num_haplo_A, ptrj->num_haplo_B, proba_diplo);
			}
			// si le diplotype est plus fréquent que les diplotypes precedants.
			if (proba_diplo > pvar->tab_geno[i].proba_diplo_max)
			{
				printf ("\tProbabilité plus haute,remplacement du couple max courant");
				pvar->tab_geno[i].proba_diplo_max = proba_diplo;
				pvar->tab_geno[i].num_haplo_A_max = ptrj->num_haplo_A;	
				pvar->tab_geno[i].num_haplo_B_max = ptrj->num_haplo_B;	
			}
			ptrj = ptrj -> suivant;
		}
		printf ("\nDiplotype final : %d / %d\tProba : %e\n\n",\
		pvar->tab_geno[i].num_haplo_A_max,\
		pvar->tab_geno[i].num_haplo_B_max,\
		pvar->tab_geno[i].proba_diplo_max);
	}
	return;
}


/***********************************************************************
 * export_geno_diplo
 **********************************************************************/

void export_geno_diplo (T_info* pvar, char* output)
{
	int i;
	int num_geno, num_haplo_A, num_haplo_B;
	
	FILE* fp = init_file_ptr ((strcat(output, "_diplotypes_inférés.txt")), "w");

	for (i = 0; i < pvar->nb_ind; i++ ) // Pour chaque genotypes de tab_geno
	{
		// Accès à tab_individus et memorisation du numéro de genotype associé 
		//fprintf(fp, "Individu # %d\n",i);
		num_geno = pvar->tab_individus[i].num_geno;
		
		// Accès au génotype associé dans tab_geno et memorisation des numéros de d'haplotypes explicatifs 
		fprintf(fp, "Genotype : %d\t%s\n", num_geno, pvar->tab_individus[i].sequence);
		num_haplo_A = pvar->tab_geno[num_geno].num_haplo_A_max;
		num_haplo_B = pvar->tab_geno[num_geno].num_haplo_B_max;
		
		// Accès aux sequences des haplotypes explicatifs dans tab_haplo 
		fprintf(fp, "HaploA %d :\t%s\n", num_haplo_A, pvar->tab_haplo[num_haplo_A].sequence);
		fprintf(fp, "HaploB %d :\t%s\n\n", num_haplo_B, pvar->tab_haplo[num_haplo_B].sequence);
	}
	fclose(fp);
	return;
}














