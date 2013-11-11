#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include "EM_main.h"

////////////////////////////////////////////////////////////////////////
// importation_genotypes.c
////////////////////////////////////////////////////////////////////////

void importation_genotypes(char* genotype_file, T_info* pvar)
{
	int i = 0;
	FILE* file = NULL;
	
	file = init_file_ptr(genotype_file, "r");
	pvar->taille = nb_char(file);
	rewind(file); // Retour au début du fichier
	pvar->nb_ind = nb_ligne(file);
	rewind(file);
	
	printf("\nTaille des génotypes = %d, Nombre d'individus = %d\n\n", pvar->taille, pvar->nb_ind);
	
	pvar->tab_individus = malloc (sizeof (T_individu) * pvar->nb_ind);
	
	for (i = 0; i < pvar->nb_ind; i++) 
	{
		pvar->tab_individus[i].sequence = malloc (1); // creation d'un espace d'1 bit 
		fgets(pvar->tab_individus[i].sequence, pvar->taille + 2, file); // remplissage et réalocation automatique par fgets...
		
		if (pvar->tab_individus[i].sequence[pvar->taille] == '\n')
			pvar->tab_individus[i].sequence[pvar->taille] = '\0';
		
		printf("Individu #%d\t Sequence: %s\t\n", i, pvar->tab_individus[i].sequence);
	}
	
	fclose(file);
	return;
}

/***** init_file_ptr **************************************************/
// name = nom du fichier à ouvrir
// mode = mode d'ouverture d'un fichier r = lecture, w = ecriture a = append
FILE* init_file_ptr (char* name, char* mode)
{
	FILE* file = NULL;
	
	file = fopen (name, mode);
	if (file == NULL)
	{
		fprintf (stderr, "Impossible d'ouvrir ou de creer le fichier %s\n\n", name );
		exit (EXIT_FAILURE);
	}
	return file;
}
 
/***** nb_lignes ******************************************************/
// Compte le nombre de lignes du fichier
int nb_ligne (FILE *fp)
{
	int n = 0, c;
	while ((c = fgetc(fp)) != EOF)
		if (c == '\n')
			n++;
	return n;
}
 
/***** nb_char_per_line ***********************************************/
// Compte le nombre de char de la première ligne du fichier
int nb_char (FILE *fp)
{
	int n = 0, c;
	while ((c = fgetc(fp)) != '\n')
		n++;
	return n;
}
