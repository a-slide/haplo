#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <unistd.h>
#include "inference_haplotype.h"
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
	
	double vraisemblance = 0;				// Valeur de vraisemblance qui sera utilisé comme contrôle de sortie de boucle EM
	double vraisemblance_prec = -DBL_MAX;	// Valeur de vraisemblance de l'itération precedente initialisée à -infini
	T_info var; 							// Structure var contenant les variables et pointeurs de structures importants
		var.tab_individus = NULL;
		var.tab_geno = NULL;
		var.tab_haplo = NULL;
	
	/******** Parsing des options avec Getopt *********/
	int optch;
    extern int opterr;
    char format[] = "aeh:f:i:s:";
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
	
	if ((filename == NULL))
	{
		printf ("\nles options -f (nom de fichier d'entrée) \n\n");
		usage(argv[0]);
	}
			
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
	
	// Tant que le nombre d'itération entré en parametre n'est pas atteint et tant qu'il n'y a pas de convergence
	while (nb_iterations != nb_iterations_max) // Pour permettre l'exception nb_iterations_max = 0
	{
		encadre ("DEBUT NOUVELLE ITERATION EM");
		// Calcul de la fréquence de chaque haplotype par maximisation
		maximisation (&var);

		// Calcul de la proba de chaque genotypes pour l'iteration en cours et de la vraisemblance
		vraisemblance = estimation_esperance (&var);
		// printf("\nvraisemblance = %e\nvraisemblance_prec = %e\n", var.vraisemblance, var.vraisemblance_prec);
		
		// Calcule de la valeur de la convergence
		convergence = (fabs(vraisemblance-vraisemblance_prec)/(vraisemblance_prec));
		printf("\nValeur_convergence = %e\n\n", convergence);
		
		// Si le seuil de convergence est atteint = sortie de EM
		if (convergence >= seuil) break;
		
		// Mise à jour des valeurs prec qui prennent les valeurs courantes
		update_proba_freq (&var);
		vraisemblance_prec = vraisemblance;
		nb_iterations ++;
	}	
	
	encadre ("SORTIE DE EM");
	printf("Iteration de sortie %d / convergence de sortie : %e (seuil : %e) \n", nb_iterations, convergence, seuil);

	/******** Creation et export des fichiers de resultats *********/
	
	encadre ("CREATION DES FICHIERS DE SORTIE ");
	
	// Definir paire plus probable pour chaques genotypes
	diplotype_plus_probable (&var);				

	// Creer fichier listant les genotypes et la paire la plus probable pour chaque individus
	export_geno_diplo (&var);
	
	qsort (var.tab_haplo, var.nb_haplo, (sizeof (T_haplo)), comparaison_frequence);
	
	// Creer fichier listant les haplotypes par ordre croissant
	export_haplo (&var);
	
	encadre ("FIN DE L'ALGORITHME");
	
	return (EXIT_SUCCESS);
}

/***********************************************************************
 *  USAGE
***********************************************************************/

void usage (char* prog_name)
{
	fprintf (stderr, "\nUSAGE:\t%s -f <chemin> -o <prefixe> [-a -e -i <entier> -s <decimal> -h]\n\n", prog_name);
	
	fprintf (stderr, "DESCRIPTION\n\n\
	Ce programme permet d'inférer les haplotypes constituant une liste de génotypes\n\n\
	La convention d'encodage des positions doit être la suivante:\n\
		0 si les 2 haplotypes possèdent l'allèle mineure\n\
		2 si les 2 haplotypes possèdent l'allèle majeure\n\
		1 si 1 haplotype possède l'allèle majeure et l'autre l'allèle mineure\n\
	Par exemple si le genotype G est constitué des haplotypes Ha = aBCd et Hb = abCd ou les minuscules\n\
	sont les allèles minneures et les majuscules les allèles majeures, alors les Ha, Hb et G\n\
	seront respectivement codés O110, 0010 et 0120.\n\n");
	
	fprintf (stderr, "DETAILS DES OPTIONS\n\n\
	-f	Fichier texte contenant une liste de génotypes pour une serie d'individus (Obligatoire)\n\
	-o	Prefixe du nom du fichier de sortie (Obligatoire)\n\n\
	-e	Mode d'initialisation des fréquences d'haplotypes equi-probable (option par defaut)\n\
	-a	Mode d'initialisation des fréquences d'haplotypes aléatoire\n\
		Les options -a et -e sont mutuellement exclusives (en cas de conflit -e sera appliquée)\n\n\
	-i	Nombre d'itération de la boucle EM. Valeur obligatoirement positive ou nulle\n\
		Par defaut = 10. la valeur peut être réglée à 0 pour boucler jusqu'à convergence\n\n\
	-s	Valeur seuil de convergence declenchant l'arret de la boucle EM. Min = -10 / Max = 0\n\
		Par defaut = -0.001. la valeur peut être réglée à 0 pour boucler jusqu'à l'iteration maximale\n\n\
	-h	Affiche cet ecran d'aide\n\n");
	
	fprintf (stderr, "EXEMPLES\n\n\
	%s -f mon_dossier/genotypes.txt -a -i 20 -s -0.00001\n\
	%s -f mon_dossier/genotypes.txt -e -i 0 -s -0.000000001\n",prog_name, prog_name);
	
	exit (EXIT_FAILURE);
}

