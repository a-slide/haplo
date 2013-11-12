#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "EM_main.h"
 
int main(int argc, char** argv)
{
	/* Arguments demandés en entrée :
	- Identificateur du fichier contenant les genotypes (string)
	- Nombre d'itérations de EM maximum (int)
	- Seuil de convergence (float) */
	
	/******** Declaration des variables *********/
	//int nb_iterations;
	//int seuil;
	char* genotype_file;
	
	// Declaration de la structure T_info var et initialisation des pointeurs à NULL
	T_info var = { .tab_individus = NULL, .tab_geno = NULL, .tab_haplo = NULL };
		
	/******** Test d'usage *********/
	if (argc != 2) usage(argv[0]); // Affiche usage et sort si nombre de paramètres incorect
	
	/******** Initialisation des variables *********/
	genotype_file = argv[1];
	//nb_iterations = atoi (argv[2]);
	//seuil = atoi (argv[3]);
	
	/******** Importation et initialisation des données *********/
	importation_genotypes (genotype_file, &var);
	preparer_liste_geno_haplo (&var);
	//initialiser_frequence_haplotype
	//calculer_frequence_genotype
	
	/******** Expectation Maximisation *********/
	
	return 0;
}
	/*
	calculer_frequence_génotype
	Convergence = faux
	nbr_etapes = 0
	Vraisemblance_prec = VALEUR_MIN // valeur min est la plus petite valeur gérer par l'ordinateur
	Tantque (Convergence == Faux && nb_etapes <= nb_etapes_max)
	incr nb_etapes
	maximisation ()
	Vraisemblance = esperance ()
	Convergence = ((| Vraisemblance – Vraisemblance_prec | / Vraisemblance_prec) <= seuil) // A REVOIR ICI.......
	si (Convergence == Faux)
	alors
	Vraisemblance_prec = Vraisemblance
	Freq_prec = freq
	proba_prec = proba
	finsi
	finTantQue
	}
	*/
 
void usage (char* prog_name)
{
	fprintf (stderr, "\nUsage:\t%s [fichier_contenant_genotypes] [nombre_iterations] [seuil_de_convergence]\n\n", prog_name);
	exit (EXIT_FAILURE);
}

