#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "EM_main.h"
 
 /***********************************************************************
 * MAIN
 **********************************************************************/
int main(int argc, char** argv)
	/* Arguments demandés en entrée :
	 * 1 = Identificateur du fichier contenant les genotypes (string)
	 * 2 = Mode d'initialisation des fréquences d'haplotypes (A)léatoire ou (E)qui-probable
	 * 3 = Nombre d'itérations de EM maximum (int)
	 * 4 = Seuil de convergence (double) */
	 
	{
	/******** Declaration des variables *********/
	int nb_iterations = 1, nb_iterations_max;
	double seuil;
	T_info var;
	double convergence = 0, conv = 0;
	
	/******** Test d'usage *********/
	// Affiche usage et sort si nombre de paramètres incorect
	if (argc != 5) usage(argv[0]);
	
	/******** Initialisation des variables *********/
	nb_iterations_max = atoi (argv[3]);
	seuil = atof (argv[4]);
	var.tab_individus = NULL;
	var.tab_geno = NULL;
	var.tab_haplo = NULL;
	var.vraisemblance = 0;
	var.vraisemblance_prec = 0;
	
	/******** Importation et initialisation des données *********/
	printf("\n##############################################################################\n");
	printf("# IMPORTATION ET PREPARATION DES DONNES \n");
	printf("##############################################################################\n");
	importation_genotypes (argv[1], &var);
	preparation_liste_geno_haplo (&var);
	printf("\n##############################################################################\n");
	printf("# INITIALIZATION\n");
	printf("##############################################################################\n");
	initialisation_freq_proba (&var, argv[2][0]);	
	
	/******** Boucle Expectation/Maximisation *********/
	// Tant que le nombre d'itération entré en paramètre n'est pas atteind et tant qu'il n'y a pas de convergence
	while (nb_iterations <= nb_iterations_max && !convergence)
	{
		printf("\n##############################################################################\n");
		printf("# ITERATION EM %d\n", nb_iterations);
		printf("##############################################################################\n");
		
		// Calcul de la fréquence de chaque haplotype par maximisation
		maximisation (&var);
		// Calcul de la proba de chaque génotypes pour l'itération en cours et 
		var.vraisemblance = estimation_esperance (&var);

		printf("vraisemblance = %.9e\nvraisemblance_prec = %.9e\n", var.vraisemblance, var.vraisemblance_prec);
		// Calcule de la valeur de la convergence
		conv = (fabs(var.vraisemblance-var.vraisemblance_prec)/(var.vraisemblance_prec));
		// Test pour voir si la valeur de la convergence est inférieur ou égale au seuil entré en paramètre. Renvoie 1 ou 0
		convergence = (conv <= seuil);
		
		printf("seuil = %.2e\n", seuil);
		printf("valeur_convergence = %.2e\n", conv);
		printf("valeur_convergence = %.2e\n\n", convergence);
	
		/* Si il n'y a pas convergence, on met à jour les fréquences à l'itération
		 * précedente des haplotypes et les probabilités à l'itération précédente
		 * des génotypes pour qu'elles prennent les nouvelles valeurs de fréquences
		 * et de probabilités calculées dans la fonction "Maximisation_et_Esperance" */
		
		if (!convergence)
			// Mise à jour des valeurs prec qui prennent les valeurs courantes
			update_proba_freq_vraisemblance (&var);
		
		nb_iterations ++;
	}
	/******** Creation et export des fichiers de resultats *********/
	
	return 1;
}

void usage (char* prog_name)
{
	fprintf (stderr, "\nUsage:\t%s [fichier__genotypes] [mode_init] [nb_iterations] [seuil_convergence]\n\n", prog_name);
	fprintf (stderr, "\t[fichier__genotypes] (string)\n");
	fprintf (stderr, "\tListe de genotype pour une serie d'individus. Valeur max conseillée = ???\n\n");
	fprintf (stderr, "\t[mode_init] (char)\n");
	fprintf (stderr, "\tMode d'initialisation des fréquences d'haplotypes (A)léatoire ou (E)qui-probable\n\n");
	fprintf (stderr, "\t[nb_iterations] (int)\n");
	fprintf (stderr, "\tNombre d'iteration de la boucle EM. Valeur max conseillée = ???\n\n");
	fprintf (stderr, "\t[seuil_convergence] (double). Valeur conseillée < (-0.177)\n");
	fprintf (stderr, "\tCondition de sortie de boucle EM si convergence atteinte. Valeur max conseillée = ???\n\n");
	exit (EXIT_FAILURE);
}
