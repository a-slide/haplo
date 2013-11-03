#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
 
typedef struct diplo_expl T_diplo_expl;
struct diplo_expl
{
	int num_haplo_A; // paire d'haplotype explicatifs = diplotype
	int num_haplo_B;
	T_diplo_expl* suivant;
};

typedef struct geno T_geno;// tête de liste de type T_diplo_expl
struct geno
{
	char* sequence;
	float proba;
	int nb_ind; // compteur d'individus avec ce genotype
	int nb_diplo_expl; // compteur de diplotype explicatifs
	T_diplo_expl* tete;
	T_diplo_expl* queue;
	//T_individu* tete;
	//T_individu* queue;
};

typedef struct geno_expl T_geno_expl;
struct geno_expl
{
	int num_geno_expl; // genotype expliqué
	int num_haplo_compl; // haplotype complémentaire pour expliquer le genotype concerné
	T_geno_expl* suivant;
};
 
typedef struct haplo T_haplo;// tête de liste de type T_geno_expl
struct haplo
{
	char* sequence;
	float frequence;
	int nb_geno_expl;
	T_geno_expl* tete;
	T_geno_expl* queue;
};
 
// Fonctions partagées par plusieurs fichiers sources
char** create_char_mat (int, int);
void print_string_table (char**, int);
 
// Fonctions specifiques de main
void usage (char*);
 
// Fonctions specifiques de importation_genotypes
void importation_genotypes(char***, char*, int*, int*);
FILE* init_file_ptr (char* name, char* mode);
int nb_ligne (FILE*);
int nb_char (FILE*);
 
// Fonctions specifiques de initialiser_frequence_haplotype
void preparer_liste_geno_haplo (char**, int, int, T_geno**, T_haplo**);

void ajouter_tab_geno (char*, int*, T_geno**);
void ajouter_tab_haplo (char*, int*, int*, T_haplo**);

void init_tab_geno (char*, T_geno**);
void init_tab_haplo (char*, T_haplo**);

void extend_tab_geno (int , char* , T_geno**);
void extend_tab_haplo (int, char*, T_haplo**);

int compte_ambiguites (char*, int);
void haplotypes_possibles (char*, int, char***, int*);
void ajouter_diplo_a_geno (int, int, int, T_geno**);
void ajouter_geno_a_haplo (int, int, int, T_haplo**);
void print_tab_haplo (int, T_haplo**, char***);
 
 
////////////////////////////////////////////////////////////////////////
// main.c
////////////////////////////////////////////////////////////////////////
 
int main(int argc, char** argv)
{
	/* Arguments demandés en entrée :
	- Identificateur du fichier contenant les genotypes (string)
	- Nombre d'itérations de EM maximum (int)
	- Seuil de convergence (float) */
	
	/******** Declaration des variables *********/
	//int nb_iterations;
	//int seuil;
	int taille_geno;
	int nb_ind;
	char* genotype_file;
	char** tab_individus; 	// Tableau de string contenant la liste de tt les génotypes de tt les individus
	T_geno* tab_geno;		// Tableau de structure T_Geno contenant une liste non redondante de génotypes
	T_haplo* tab_haplo;		// Tableau de structure T_Haplo contenant une liste non redondante d'haplotypes
	
	/******** Test d'usage *********/
	if (argc != 2) usage(argv[0]); // Affiche usage et sort si nombre de paramètres incorect
	
	/******** Initialisation des variables *********/
	genotype_file = argv[1];
	//nb_iterations = atoi (argv[2]);
	//seuil = atoi (argv[3]);
	
	/******** Importation et initialisation des données *********/
	importation_genotypes (&tab_individus, genotype_file, &taille_geno, &nb_ind);
	preparer_liste_geno_haplo (tab_individus, taille_geno, nb_ind, &tab_geno, &tab_haplo);
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


////////////////////////////////////////////////////////////////////////
// importation_genotypes.c
////////////////////////////////////////////////////////////////////////
 
/**** importation_genotypes *******************************************/
 
void importation_genotypes(char*** p_tab_individus, char* genotype_file, int* p_taille_geno, int* p_nb_individus)
{
	int i = 0;
	int taille_geno, nb_individus;
	FILE* file = NULL;
	char** tab_individus;
	
	file = init_file_ptr(genotype_file, "r");
	taille_geno = nb_char(file);
	rewind(file); // Retour au début du fichier
	nb_individus = nb_ligne(file);
	tab_individus = create_char_mat (nb_individus, taille_geno+1); // Création tableau vide destiné à contenir les génotypes
	rewind(file);

	while( i < nb_individus &&(fgets(tab_individus[i], taille_geno + 2, file) != NULL) )
	{
		if (tab_individus[i][taille_geno] == '\n')
			tab_individus[i][taille_geno] = '\0';
		i++;
	}
	fclose(file);
	printf("\nTaille des génotypes = %d, Nombre de génotypes = %d\n\n", taille_geno, nb_individus);
	print_string_table (tab_individus, nb_individus); // Affichage du Tableau rempli
	
	*p_taille_geno = taille_geno; // Retour de la taille et du nombre de génotype au main
	*p_nb_individus = nb_individus;
	*p_tab_individus = tab_individus;
	return;
}
 
/**** create_char_mat *************************************************/
 
char** create_char_mat (int col, int line)
{
	int i;
	char** matrice;
	matrice = malloc (sizeof (char*) * col); // creation colonnes
	
	if (matrice == NULL)	// Verification de l'allocation
	{
		fprintf (stderr, "Allocation impossible\n\n");
		exit (EXIT_FAILURE);
	}
	for ( i = 0; i < col; i++ )
	{
		matrice[i] = malloc (sizeof (char) * line); // creation lignes
		if (matrice[i] == NULL)
		{
			fprintf (stderr, "Allocation impossible\n\n");
			exit (EXIT_FAILURE);
		}
	}
	return matrice;
}

/**** print_string_table **********************************************/
 
void print_string_table (char** tab, int line)
{
	int i;
	 
	for (i = 0 ; i < line ; i++)
		printf("#%i : %s\n", i, tab[i]);
	printf("\n");
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


////////////////////////////////////////////////////////////////////////
// preparer_liste_geno_haplo.c
////////////////////////////////////////////////////////////////////////
 
/**** preparer_liste_geno_haplo *********************************/
 
void preparer_liste_geno_haplo (char** tab_individus, int taille_geno, int nb_individus, T_geno** p_tab_geno, T_haplo** p_tab_haplo)
{
	int i, j;
	// Variables pour la liste de genotype non redondante (tab_geno)
	int nb_geno; // Nombre de genotypes totalno 
	///int num_geno; // Indice du genotype
	// Variables pour la liste d'haplotypes explicatifs générée par genotype (tab_haplo_expl)
	char** tab_haplo_expl; // Pour stocker les haplotypes explicatifs
	int nb_haplo_expl; // Nombre d'haplotypes explicatifs
	// Variables pour la liste d'haplotype non redondante (tab_haplo)
	int nb_haplo; // Nombre d'haplotype total
	int num_haplo_A, num_haplo_B; // Indices d'une paire d'haplotypes explicatifs
	
	*p_tab_geno = NULL;
	*p_tab_haplo = NULL;
	nb_geno = 0;
	nb_haplo = 0;
	
	// Création d'une liste de genotypes non redondants à partir de tab_individus
	for (i = 0 ; i < nb_individus ; i++)
		ajouter_tab_geno (tab_individus[i], &nb_geno, p_tab_geno);
	
	printf("\nListe non redondante de genotypes\n");
	for (i = 0; i < nb_geno ; i ++)
			printf("Genotype #%d\t Sequence %s\t Individus concernés %d\n", i, p_tab_geno[0][i].sequence, p_tab_geno[0][i].nb_ind);
	
	// Parcours du tableau de structure T_geno pour analyser les genotypes 1 par 1
	for (i = 0; i < nb_geno ; i ++)
	{
		// Genere une liste des haplotypes explicatifs potentiels pour le genotype i
		haplotypes_possibles (p_tab_geno[0][i].sequence, taille_geno, &tab_haplo_expl, &nb_haplo_expl);
		///print_string_table (tab_haplo_expl, nb_haplo_expl); // verification des haplotypes générés
	
		for (j = 0; j < nb_haplo_expl; j +=2 ) // Boucle parcourant les haplotypes explicatifs par paire complémentaire (diplotypes)
		{
			// Ajout des haplotypes A et B d'une paire explicative à la liste des haplotypes si ils n'existent pas et retour de leur indice
			ajouter_tab_haplo (tab_haplo_expl[j], &nb_haplo, &num_haplo_A, p_tab_haplo);
			ajouter_tab_haplo (tab_haplo_expl[j+1], &nb_haplo, &num_haplo_B, p_tab_haplo);
			///printf("Diplotype #%d: haplotype #%d haplotype #%d\n\n ", j/2, num_haplo_A, num_haplo_B);
			
			// Mettre les genotypes explicatifs avec l'haplo complementaire dans la liste chainée associée à chacun des 2 haplotypes
			ajouter_geno_a_haplo (num_haplo_A, num_haplo_B, i, p_tab_haplo);
			ajouter_geno_a_haplo (num_haplo_B, num_haplo_A, i, p_tab_haplo);
			
			//Mettre les paires d'haplotypes dans la liste chainée associée au génotype courant
			ajouter_diplo_a_geno (num_haplo_A, num_haplo_B, i, p_tab_geno);
		}
		// Desalouer l'espace mémoire de tab_haplo_expl
	}
	
	printf("\nListe non redondante d'haplotypes\n");
	for (i = 0; i < nb_haplo ; i ++)
			printf("Haplotype #%d\t Sequence %s\t Genotypes concernés %d\n", i, p_tab_haplo[0][i].sequence, p_tab_haplo[0][i].nb_geno_expl);
		
	///print_tab_haplo (nb_haplo, &tab_haplo, &tab_geno);
	return;
}
 
///** ajouter un genotype de tableau tab_geno *************************/
void ajouter_tab_geno (char* geno_seq, int* p_nb_geno, T_geno** p_tab_geno)
{
	int i;
	// Si la matrice est vide = ajout du premier genotype
	if (*p_tab_geno==NULL) 
	{
		///printf("Ajout du premier genotype en position 0\n");
		*p_nb_geno = 1;
		init_tab_geno (geno_seq, p_tab_geno);
		return;
	}
	// Sinon recherche si le genotype existe déjà dans la liste
	for (i = 0; i < *p_nb_geno; i++)
	{ 
		if ( strcmp(p_tab_geno[0][i].sequence, geno_seq) == 0 )
		{
			///printf("Le genotype existe déjà en position %d\n", i);
			p_tab_geno[0][i].nb_ind++;	
			return;	
		}
	}
	// Sinon stockage du genotype dans une nouvelle case à la suite du tableau
	///printf("Ajout d'un nouveau génotype en position %d\n", i);
	*p_nb_geno = i+1;
	extend_tab_geno (*p_nb_geno, geno_seq, p_tab_geno);
	return;
}
 
///** initialisation tableau genotype *********************************/
void init_tab_geno (char* geno_seq, T_geno** p_tab_geno)
{
	*p_tab_geno = malloc (sizeof (T_geno)); // taille = 1 pour le premier element
	
	if (*p_tab_geno == NULL)
	{
		fprintf (stderr, "Allocation impossible\n\n");
		exit (EXIT_FAILURE);
	}
	p_tab_geno[0][0].sequence = geno_seq;
	p_tab_geno[0][0].proba = 0;
	p_tab_geno[0][0].nb_ind = 1;
	p_tab_geno[0][0].nb_diplo_expl = 0;
	p_tab_geno[0][0].tete = NULL;
	p_tab_geno[0][0].queue = NULL;
	return;	
}
 
 /// Ajout d'une case au tableau genotype ///////////////////////////////////
void extend_tab_geno (int nb_geno, char* geno_seq, T_geno** p_tab_geno)
 {
	*p_tab_geno = realloc (*p_tab_geno, sizeof (T_geno) * nb_geno);
	
	if (*p_tab_geno == NULL)
	{
		fprintf (stderr, "Allocation impossible\n\n");
		exit (EXIT_FAILURE);
	}
	p_tab_geno[0][nb_geno-1].sequence = geno_seq;
	p_tab_geno[0][nb_geno-1].nb_ind = 1;
	p_tab_geno[0][nb_geno-1].nb_diplo_expl = 0;
	p_tab_geno[0][nb_geno-1].tete = NULL;
	p_tab_geno[0][nb_geno-1].queue = NULL;
	return;
} 
 
/**** haplotypes_possibles ********************************************/
void haplotypes_possibles (char* genotype, int taille_geno, char*** p_tab_haplo, int* p_nb_haplo)
{
	int nb_amb, amb; // nombre total d'ambiguité dans la séquence et compteur d'ambiguité courante
	int change_prog; // variable booléenne utilisé pour remplir le tableau des haplotype possible en suivant alternativement 1/0 ou 0/1
	int val_basc, basc; //valeur à attendre pour basculer change_prog et compteur de basculement
	int nb_haplo; // calculé à partir du nombre d'ambiguités
	int j, k; // variables de contrôle de boucle
	char** tab_haplo; // tableau temporaire permettant de stocker les haplotypes explicatifs générés.
	
	nb_amb = compte_ambiguites (genotype, taille_geno); // calcul du nombre d'ambiguités dans le genotype courant
	///printf("Nombre d'ambiguités : %d\n", nb_amb);
	
	if (nb_amb != 0) nb_haplo = exp2(nb_amb);
	else nb_haplo = 2; // si 0 ambiguité, il faut quand même générer 2 haplotypes
	 
	tab_haplo = create_char_mat (nb_haplo, taille_geno+1);
	amb = 0;
	
	for (j = 0; j < taille_geno ; j ++)
	{
		if (genotype[j] == '0') // Garnissage de tous les champs de la table haplo à cette même position avec 0
			for (k = 0; k < nb_haplo ; k ++)
				tab_haplo[k][j] = '0';

		else if (genotype[j] == '2') // Garnissage tous les champs de la table haplo à cette même position avec 0
			for (k = 0; k < nb_haplo ; k ++)
				tab_haplo[k][j] = '1';

		else  // cas le plus complexe de position ambigue
		{
			change_prog = 0;
			val_basc = nb_haplo/exp2(amb++); // Calcul le nbr d'haplotype aprés lequel il faut changer le sens de remplissage
			for (k = 0, basc = 0 ; k < nb_haplo ; k +=2, basc +=2)  // k = compteur de boucle
			{
				if (basc == val_basc) // Si le compteur de basculement atteint la valeur de bascule
				{ 
					change_prog = !change_prog ; // bascule du booléen
					basc = 0 ; // reinitialisation du compteur de basculement
				}
				if (change_prog == 0) // remplissage sens 0/1
				{
					tab_haplo[k][j] = '0';
					tab_haplo[k+1][j] = '1';
				}
				else
				{ // remplissage sens 1/0
					tab_haplo[k][j] = '1';
					tab_haplo[k+1][j] = '0';
				}
			}
		}
	}
	*p_tab_haplo = tab_haplo;
	*p_nb_haplo = nb_haplo;
	return;
}
 
/**** compte_ambiguites ***********************************************/
int compte_ambiguites (char* genotype, int taille_geno)
{
	int i;
	int nb_amb = 0;
	
	for (i = 0 ; i < taille_geno; i++)
		if (genotype[i] == '1')
			nb_amb ++;
	return nb_amb;
}
 
/**** ajouter_tab_haplo ***********************************************/
void ajouter_tab_haplo (char* haplo_seq, int* p_nb_haplo, int* p_haplo_num, T_haplo** p_tab_haplo)
{
	int i;
	// Si la matrice est vide = ajout du premier haplotype
	if (*p_tab_haplo == NULL) { 
		///printf("Ajout du premier haplotype en position 0\n");
		*p_haplo_num = 0;
		*p_nb_haplo = 1;
		init_tab_haplo (haplo_seq, p_tab_haplo);
		return;
	}
	// Sinon recherche si l'haplotype existe déjà dans la liste
	for (i = 0; i < *p_nb_haplo; i++)
	{
		if ( strcmp(p_tab_haplo[0][i].sequence, haplo_seq) == 0 ) // si chaines identiques
		{
			///printf("L'haplotype existe déjà en position %d\n", i);
			*p_haplo_num = i;
			p_tab_haplo[0][i].nb_geno_expl++;
			return;
		}
	}
	// Sinon stockage de l'haplotype dans une nouvelle case à la suite du tableau
	///printf("Ajout d'un nouvel haplotype en position %d\n", i);
	*p_haplo_num = i;
	*p_nb_haplo = i+1;
	extend_tab_haplo (*p_nb_haplo, haplo_seq, p_tab_haplo);
	return;
}
 
/**** init_tab_haplo ************************************************/
void init_tab_haplo (char* haplo_seq, T_haplo** p_tab_haplo)
{
	*p_tab_haplo = malloc (sizeof (T_haplo)); // taille = 1 pour le premier element
	if (*p_tab_haplo == NULL)
	{
		fprintf (stderr, "Allocation impossible\n\n");
		exit (EXIT_FAILURE);
	}
	p_tab_haplo[0][0].sequence = haplo_seq;
	p_tab_haplo[0][0].frequence = 0;
	p_tab_haplo[0][0].nb_geno_expl = 1;
	p_tab_haplo[0][0].tete = NULL;
	p_tab_haplo[0][0].queue = NULL;
	return;
}
 
/**** extend_tab_haplo **********************************************/
void extend_tab_haplo (int nb_haplo, char* haplo_seq, T_haplo** p_tab_haplo)
{
	*p_tab_haplo = realloc (*p_tab_haplo , sizeof (T_haplo) * nb_haplo);
	if (*p_tab_haplo == NULL)
	{
		fprintf (stderr, "Allocation impossible\n\n");
		exit (EXIT_FAILURE);
	}
	p_tab_haplo[0][nb_haplo-1].sequence = haplo_seq;
	p_tab_haplo[0][nb_haplo-1].frequence = 0;
	p_tab_haplo[0][nb_haplo-1].nb_geno_expl = 1;
	p_tab_haplo[0][nb_haplo-1].tete = NULL;
	p_tab_haplo[0][nb_haplo-1].queue = NULL;
	return;
}
 
/**** ajouter_diplo_a_geno ********************************************/
void ajouter_diplo_a_geno (int num_haplo_A, int num_haplo_B, int num_geno, T_geno** p_tab_geno)
{
	T_diplo_expl* cellule = malloc( sizeof (T_diplo_expl));
	if (cellule == NULL)
	{
		fprintf (stderr, "Allocation impossible\n\n");
		exit (EXIT_FAILURE);
	}
	cellule -> num_haplo_A = num_haplo_A;
	cellule -> num_haplo_B = num_haplo_B;
	cellule -> suivant = NULL;

	if (p_tab_geno[0][num_geno].tete == NULL) // ajout de la première cellule à la liste
	{
		p_tab_geno[0][num_geno].tete = cellule;
		p_tab_geno[0][num_geno].queue = cellule;
	}
	else // ajout en queue de liste = chainage simplifié
	{
		p_tab_geno[0][num_geno].queue -> suivant = cellule; // chainage avec cellule precedente
		p_tab_geno[0][num_geno].queue = cellule; // avancée du pointeur
	}
	return;
}
 
/**** ajouter_geno_a_haplo ********************************************/
void ajouter_geno_a_haplo (int num_haplo_principal, int num_haplo_compl, int num_geno_expl, T_haplo** p_tab_haplo)
{
	T_geno_expl* cellule = malloc( sizeof (T_geno_expl));
	if (cellule == NULL)
	{
		fprintf (stderr, "Allocation impossible\n\n");
		exit (EXIT_FAILURE);
	}
	cellule -> num_geno_expl = num_geno_expl;
	cellule -> num_haplo_compl = num_haplo_compl;
	cellule -> suivant = NULL;

	if (p_tab_haplo[0][num_haplo_principal].tete == NULL) // Ajout de la première cellule à la liste
	{
		p_tab_haplo[0][num_haplo_principal].tete = cellule;
		p_tab_haplo[0][num_haplo_principal].queue = cellule;
	}
	else // Ajout en queue de liste = chainage simplifié
	{
		p_tab_haplo[0][num_haplo_principal].queue -> suivant = cellule; // Chainage avec cellule precedente
		p_tab_haplo[0][num_haplo_principal].queue = cellule; // Avancée du pointeur
	}
	return;
}
 
/**** print_tab_haplo ********************************************/
void print_tab_haplo (int nb_haplo, T_haplo** p_tab_haplo, char*** tab_geno)
{
	int i;
	T_geno_expl* ptrj = p_tab_haplo[0][0].tete;
	printf("Liste non redondante des haplotypes et des genotype expliqué\n\n");
	for (i = 0; i < nb_haplo; i++ )
	{
		printf("Haplotype #%d \t Séquence : %s \t Fréquence %f \t Nombre de genotype(s) expliqué(s) : %d \n",i, p_tab_haplo[0][i].sequence, p_tab_haplo[0][i].frequence, p_tab_haplo[0][i].nb_geno_expl);
		printf("Liste des génotypes expliqués\n");
		while (ptrj != NULL)
		{
			//printf("Genotype # %d (%s) avec Haplotype # %d (%s)\n", ptrj -> num_geno_expl, *tab_geno[ptrj->num_geno_expl], ptrj -> num_haplo_compl, p_tab_haplo[0][ptrj->num_haplo_compl].sequence );
			printf("#");
			ptrj = ptrj -> suivant;
		}
		printf("\n");
	}
}
