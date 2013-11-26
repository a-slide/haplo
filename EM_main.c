#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <unistd.h>
#include "EM_main.h"
#include "ptr_allocation.h"
 
/***********************************************************************
 *  MAIN
 **********************************************************************/

int main(int argc, char** argv)
	{
	/******** Declaration et initialisation des variables *********/
	int nb_iterations = 1; 					// Compteur d'iteration courante pour EM
	int nb_iterations_max = 10; 			// Nombre maximal de cycles EM pouvant être modifie par l'utilisateur
	int mode_init = 0; 						// Mode d'initialisation des frequences d'haplotype 0 = equiprobable / 1 = aleatoire 
	double seuil = -0.001;					// Seuil de convergence de EM pouvant être modifie par l'utilisateur
	double convergence = 0;					// Valeur de la convergence calcule à chaque iteration
	char* filename = NULL;					// Nom du fichier contenant la liste de genotypes des individus
	double vraisemblance = 0;				// Valeur de vraisenblance qui sera utilisé comme contrôle de sortie de boucle EM
	double vraisemblance_prec = -DBL_MAX;	// Valeur de vraisenblance de l'itération precedente initialisée à -infini
	T_info var; 							// Structure var contenant les variables et pointeurs de structures importants
		var.tab_individus = NULL;
		var.tab_geno = NULL;
		var.tab_haplo = NULL;
	
	/******** Parsing des options avec Getopt *********/
	int optch;
    extern int opterr;
    char format[] = "aehf:i:s:";
	opterr = 1;
	
    while ((optch = getopt(argc, argv, format)) != -1)
    switch (optch)
    {
        case 'a':
            mode_init = 1;
            break;
        case 'e':
            mode_init = 0;
            break;
		case 'i':
            nb_iterations_max = atoi(optarg);
            if (nb_iterations_max < 0) usage(argv[0]);
            break;
		case 'f':
            filename = optarg;
            break;
        case 's':
            seuil = atof (optarg);
            if ((seuil < -10)||(seuil > 0)) usage(argv[0]);
			break;
        case 'h':
			usage(argv[0]);
            break;
    }
	
	if (filename == NULL) usage(argv[0]);
	
	encadre ("PARAMETRES D'ENTREE ET REGLAGES");
	printf ("\tFichier contenant les genotypes = %s\n",filename);
	printf ("\tMode d'initialisation des haplotypes = %s\n", ((mode_init == 1) ? "Aleatoire" : "Equiprobable"));
	printf ("\tNombre d'iterations maximum = %d\n",nb_iterations_max );
	printf ("\tSeuil de convergence = %e\n",seuil);
	
	/******** Importation et initialisation des donnees *********/
	
	encadre ("IMPORTATION ET PREPARATION DES DONNEES");
	// Importation des genotypes des individus
	importation_genotypes (filename, &var);
	// Creation des structures de donnees contenant les listes non redondantes d'haplotypes et de genotypes
	preparation_liste_geno_haplo (&var);
	
	// Initialisation des frequences des haplotypes et des probabilites des genotypes
	encadre ("INITIALISATION");
	initialisation_freq_proba (&var, mode_init);
	
	/******** Boucle Expectation/Maximisation *********/
	
	// Tant que le nombre d'iteration entre en parametre n'est pas atteint et tant qu'il n'y a pas de convergence
	while (nb_iterations != nb_iterations_max) // Pour permettre l'exception nb_iterations_max = 0
	{
		encadre ("DEBUT NOUVELLE ITERATION EM");
		// Calcul de la frequence de chaque haplotype par maximisation
		maximisation (&var);

		// Calcul de la proba de chaque genotypes pour l'iteration en cours et de la vraisemblance
		vraisemblance = estimation_esperance (&var);
		// printf("\nvraisemblance = %e\nvraisemblance_prec = %e\n", var.vraisemblance, var.vraisemblance_prec);
		
		// Calcule de la valeur de la convergence
		convergence = (fabs(vraisemblance-vraisemblance_prec)/(vraisemblance_prec));
		printf("\nValeur_convergence = %e\n\n", convergence);
		
		// Si le seuil de convergence est atteint = sortie de EM
		if (convergence >= seuil)
			break;
		
		// Mise à jour des valeurs prec qui prennent les valeurs courantes
		update_proba_freq (&var);
		vraisemblance_prec = vraisemblance;
		nb_iterations ++;
	}	
	
	encadre ("SORTIE DE EM");
	printf("Iteration de sortie %d / convergence de sortie : %e (seuil : %e) \n", nb_iterations, convergence, seuil);

	
	/******** Creation et export des fichiers de resultats *********/
	
	// definir paire plus probable pour chaque genotypes
	// creer fichier listantnt les genotypes et la paire la plus probable pour chaque individus
	// creer fichier listant les haplotypes par ordre croissant
	// merci bonsoir!
	
	return 0;
}

/***********************************************************************
 *  USAGE
***********************************************************************/

void usage (char* prog_name)
{
	fprintf (stderr, "\nUSAGE:\t%s -f <chemin> [-a -e -i <entier> -s <decimal> -h]\n\n", prog_name);
	fprintf (stderr, "Ce programme permet d'inferer les haplotypes constituant une liste de génotypes\n\n");
	
	fprintf (stderr, "DETAILS DES OPTIONS\n\n");
	
	fprintf (stderr, "\t-f\t fichier texte contenant une liste de genotypes pour une serie d'individus\n");
	fprintf (stderr, "\t\tIl s'agit de la seule option obligatoire du programme\n");
	fprintf (stderr, "\t\tDans ce programme l'encodage des positions ambigues adopté est 1\n");
	fprintf (stderr, "\t\tPar exemple si le genotype G est constitué des haplotypes Ha = aBCd et Hb = abCd\n");
	fprintf (stderr, "\t\t ou les minuscules sont les alleles minoritaires et les majuscules les allèles\n");
	fprintf (stderr, "\t\t majoritaires, alors les Ha, Hb et G seront respectivement codés O110, 0010 et 0120\n\n");
	
	fprintf (stderr, "\t-e\tMode d'initialisation des frequences d'haplotypes equi-probable (option par defaut)\n");
	fprintf (stderr, "\t-a\tMode d'initialisation des frequences d'haplotypes aléatoire\n");
	fprintf (stderr, "\t\tLes options -a et -e sont mutuellement exclusives (en cas de conflit -e sera appliquée)\n\n");
	
	fprintf (stderr, "\t-i\tNombre d'iteration de la boucle EM. Valeur obligatoirement positive ou nulle\n");
	fprintf (stderr, "\t\tPar defaut = 10. la valeur peut être réglée à 0 pour boucler jusqu'à convergence\n\n");

	fprintf (stderr, "\t-s\tValeur seuil de convergence declenchant l'arret de la boucle EM. Min = -10 / Max = 0\n");
	fprintf (stderr, "\t\tPar defaut = -0.001. la valeur peut être réglée à 0 pour boucler jusqu'à l'iteration maximale\n\n");
	
	fprintf (stderr, "\t-h\tAffiche cet ecran d'aide\n\n");
	
	fprintf (stderr, "EXEMPLES\n\n");
	fprintf (stderr, "\t%s -f mon_dossier/genotypes.txt -a -i 20 -s -0.00001 \n", prog_name);
	fprintf (stderr, "\t%s -f mon_dossier/genotypes.txt -e -i 0 -s -0.000000001 \n", prog_name);
	
	exit (EXIT_FAILURE);
}

/***********************************************************************
 * ENCADRE = encadre un chaine de charactere par des #
 **********************************************************************/

void encadre (char* name)
{
	int i;
	
	printf("\n\n");
	
	for (i = 0; i < 75; i++)
		printf("#");
		
	printf("\n# %s\n", name);

	for (i = 0; i < 75; i++)
		printf("#");
	
	printf("\n\n");
	return;
}
