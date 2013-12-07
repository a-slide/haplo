#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "inference_haplotype.h"
#include "ptr_allocation.h"

/***********************************************************************
 * importation_genotypes
 **********************************************************************/
//PRE CONDITIONS
//POST CONDITIONS

//Importation des génotypes des individus depuis le fichier vers un tableau dynamique 
void importation_genotypes(char* genotype_file, T_info* pvar)
{
	int i;
	int taille; 
	int nb_ind; 
	FILE* file = NULL;
	
	file = init_file_ptr(genotype_file, "r");
	taille = nb_char(file);
	rewind(file); // Retour au début du fichier
	nb_ind = nb_ligne(file);
	rewind(file);
	
	printf("\nTaille des génotypes = %d, Nombre d'individus = %d\n\n", taille, nb_ind);
	
	pvar->tab_individus = malloc (sizeof (T_individu) * nb_ind); // Allocation mémoire pour le tableau dynamique
	if (pvar->tab_individus == NULL) error_and_exit(); // Si l'allocation mémoire n'a pas fonctionné, on affiche une erreur
	
	for (i = 0; i < nb_ind; i++) // Pour chaque individu
	{
		pvar->tab_individus[i].sequence = malloc_char_string(taille + 1);
		fgets(pvar->tab_individus[i].sequence, taille + 2, file); // Remplissage du tableau par fgets
		pvar->tab_individus[i].sequence[taille] = '\0'; // Remplacement de \n par \0
		///printf("Individu #%d\t Sequence: %s\t\n", i, pvar->tab_individus[i].sequence);
	}
	
	fclose(file);
	pvar->taille = taille;
	pvar->nb_ind = nb_ind;
	
	return;
}

/***********************************************************************
 * nb_lignes
 **********************************************************************/
// Compte le nombre de lignes du fichier
int nb_ligne (FILE *fp)
{
	int n = 0, c;
	while ((c = fgetc(fp)) != EOF)
		if (c == '\n')
			n++;
	return n;
}

/***********************************************************************
 * nb_char_per_line
 **********************************************************************/
// Compte le nombre de char de la première ligne du fichier
int nb_char (FILE *fp)
{
	int n = 0, c;
	while ((c = fgetc(fp)) != '\n')
		n++;
	return n;
}
