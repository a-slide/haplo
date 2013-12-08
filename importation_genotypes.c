#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "inference_haplotype.h"
#include "ptr_allocation.h"

/***********************************************************************
 * importation_genotypes
 **********************************************************************/
/*PRE CONDITIONS
Les génotypes se trouvent dans le fichier créé par le programme "generateur_genotypes"
POST CONDITIONS
Les génotypes sont importés depuis le fichier dans un tableau alloué dynamiquement nommé "tab_individus"
*/
//Importation des génotypes des individus depuis le fichier vers un tableau dynamique 
void importation_genotypes(char* genotype_file, T_info* pvar)
{
	int i;
	int taille; // Taille des génotypes
	int nb_ind; // Nombre d'individus dans le fichier
	FILE* file = NULL;
	
	file = init_file_ptr(genotype_file, "r"); // Ouverture du fichier en mode écriture
	taille = nb_char(file); // La taille des génotypes vaut le nombre de caractère de la première ligne du fichier
	rewind(file); // Retour au début du fichier
	nb_ind = nb_ligne(file); // Le nombre d'individus vaut le nombre de lignes dans le fichier
	rewind(file);
	
	printf("\nTaille des génotypes = %d, Nombre d'individus = %d\n\n", taille, nb_ind);
	
	pvar->tab_individus = malloc (sizeof (T_individu) * nb_ind); // Allocation mémoire pour le tableau dynamique
	if (pvar->tab_individus == NULL) error_and_exit(); // Si l'allocation mémoire n'a pas fonctionné, on affiche une erreur
	
	for (i = 0; i < nb_ind; i++) // Pour chaque individu
	{
		pvar->tab_individus[i].sequence = malloc_char_string(taille + 1); // Allocation mémoire pour la séquence dans le tableau d'individus
		fgets(pvar->tab_individus[i].sequence, taille + 2, file); // Remplissage du tableau par fgets
		pvar->tab_individus[i].sequence[taille] = '\0'; // Remplacement de \n par \0
		//printf("Individu #%d\t Sequence: %s\t\n", i, pvar->tab_individus[i].sequence);
	}
	
	fclose(file); // Fermeture du fichier
	pvar->taille = taille;
	pvar->nb_ind = nb_ind;
	
	return;
}

/***********************************************************************
 * nb_lignes
 **********************************************************************/
// Comptage du nombre de lignes dans le fichier
int nb_ligne (FILE *fp)
{
	int n = 0, c;
	while ((c = fgetc(fp)) != EOF) // Tant que la fin du fichier n'est pas atteinte
		if (c == '\n') // Si on recontre un retour à la ligne, on incrémente le nombre de lignes comptées
			n++;
	return n; // On retourne le nombre de lignes total
}

/***********************************************************************
 * nb_char_per_line
 **********************************************************************/
// Comptage du nombre de caractère de la première ligne du fichier
int nb_char (FILE *fp)
{
	int n = 0, c;
	while ((c = fgetc(fp)) != '\n') // Tant que la fin de la première ligne du fichier n'est pas atteinte
		n++; // On incrémente le nombre de caractère
	return n; // On retourne le nombre de caractère total de la première ligne
}
