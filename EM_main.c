#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "EM_main.h"
 
int main(int argc, char** argv)
{
	/* Arguments demandés en entrée :
	 * 1 = Identificateur du fichier contenant les genotypes (string)
	 * 2 = Mode d'initialisation des fréquences d'haplotypes (A)léatoire ou (E)qui-probable
	 * 3 = Nombre d'itérations de EM maximum (int)
	 * 4 = Seuil de convergence (double)
	 */
	
	/******** Declaration des variables *********/
	///int nb_iterations;
	///int seuil
	T_info var;
	
	/******** Test d'usage *********/
	// Affiche usage et sort si nombre de paramètres incorect
	if (argc != 3) usage(argv[0]);
	
	/******** Initialisation des variables *********/
	///nb_iterations = atoi (argv[3]);
	///seuil = atoi (argv[4]);
	var.tab_individus = NULL;
	var.tab_geno = NULL;
	var.tab_haplo = NULL;
	
	/******** Importation et initialisation des données *********/
	importation_genotypes (argv[1], &var);
	preparation_liste_geno_haplo (&var);
	initialisation_freq_proba (&var, argv[2][0]);
	
	/******** Expectation Maximisation *********/
	
	return 0;
}
	/*
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
	fprintf (stderr, "\nUsage:\t%s [fichier__genotypes] [mode_init] [nb_iterations] [seuil_convergence]\n\n", prog_name);
	fprintf (stderr, "\t[fichier__genotypes] (string)\n");
	fprintf (stderr, "\tListe de genotype pour une serie d'individus. Valeur max conseillée = ???\n\n");
	fprintf (stderr, "\t[mode_init] (char)\n");
	fprintf (stderr, "\tMode d'initialisation des fréquences d'haplotypes (A)léatoire ou (E)qui-probable\n\n");
	fprintf (stderr, "\t[nb_iterations] (int)\n");
	fprintf (stderr, "\tNombre d'iteration de la boucle EM. Valeur max conseillée = ???\n\n");
	fprintf (stderr, "\t[seuil_convergence] (double)\n");
	fprintf (stderr, "\tCondition de sortie de boucle EM si convergence atteinte. Valeur max conseillée = ???\n\n");
	exit (EXIT_FAILURE);
}
