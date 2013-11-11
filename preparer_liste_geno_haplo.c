#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include "EM_main.h"


////////////////////////////////////////////////////////////////////////
// preparer_liste_geno_haplo.c
////////////////////////////////////////////////////////////////////////
 
/**** preparer_liste_geno_haplo ***************************************/
// T_individu* tab_individus, int taille, int nb_individus, T_geno** p_tab_geno, T_haplo** p_tab_haplo
void preparer_liste_geno_haplo (T_info* pvar)
{
	int i, j;
	int num_geno, num_haplo_A, num_haplo_B; // Indices d'un genotype et d'une paire d'haplotypes explicatifs

	// Création d'une liste de genotypes non redondants à partir de tab_individus
	
	printf("\nListe des individus et numero de genotype correspondant\n");
	for (i = 0 ; i < pvar->nb_ind; i++)
	{
		num_geno = ajouter_tab_geno (pvar->tab_individus[i].sequence, pvar);
		printf("Individu #%d\t Genotype %d\n",i, num_geno);
		//ajouter_geno_a_individu (i, pvar);
	}	
	printf("\nListe non redondante de genotypes\n");
	for (i = 0; i < pvar->nb_geno ; i ++)
			printf("Genotype #%d\t Sequence %s\t Individus concernés %d\n",\
			i, pvar->tab_geno[i].sequence, pvar->tab_geno[i].nb_ind);
	
	// Parcours du tableau de structure T_geno pour analyser les genotypes 1 par 1
	for (i = 0; i < pvar->nb_geno ; i ++)
	{
		// Genere une liste des haplotypes explicatifs potentiels pour le genotype i
		haplotypes_possibles (pvar->tab_geno[i].sequence, pvar);
		//print_string_table (pvar->tab_haplo_expl, pvar->nb_haplo_expl); //  Verification des haplotypes générés
		pvar->tab_geno[i].nb_diplo_expl = pvar->nb_haplo_expl/2; // ajout du nbr de couple explivatifs à tab_geno 
		
		for (j = 0; j < pvar->nb_haplo_expl; j +=2 ) // Parcours des haplotypes explicatifs par paire complémentaire
		{
			// Ajout des haplotypes A et B d'une paire explicative à la liste des haplotypes si ils n'existent pas et retour de leur indice
			num_haplo_A = ajouter_tab_haplo (pvar->tab_haplo_expl[j], pvar);
			num_haplo_B = ajouter_tab_haplo (pvar->tab_haplo_expl[j+1], pvar);
			///printf("Diplotype #%d: haplotype #%d haplotype #%d\n\n ", j/2, num_haplo_A, num_haplo_B);
			
			// Mettre les genotypes explicatifs avec l'haplo complementaire dans la liste chainée associée à chacun des 2 haplotypes
			ajouter_geno_a_haplo (num_haplo_A, num_haplo_B, i, pvar);
			ajouter_geno_a_haplo (num_haplo_B, num_haplo_A, i, pvar);
			
			//Mettre les paires d'haplotypes dans la liste chainée associée au génotype courant
			ajouter_diplo_a_geno (num_haplo_A, num_haplo_B, i, pvar);
		}
		//liberer_char_mat (pvar->tab_haplo_expl, pvar->nb_haplo_expl);
	}
	
	printf("\nListe non redondante d'haplotypes\n");
	for (i = 0; i < pvar->nb_haplo ; i ++)
			printf("Haplotype #%d\t Sequence %s\t Genotypes concernés %d\n",\
				i, pvar->tab_haplo[i].sequence, pvar->tab_haplo[i].nb_geno_expl);
		
	print_tab_haplo (pvar);
	print_tab_geno (pvar);
	return;
}
 
///** ajouter_tab_geno  ***********************************************/
int ajouter_tab_geno (char* geno_seq, T_info* pvar)
//char* geno_seq, int* p_nb_geno, T_geno** p_tab_geno
{
	int i;
	// Si la matrice est vide = ajout du premier genotype
	if (pvar->tab_geno == NULL) 
	{
		///printf("Ajout du premier genotype en position 0\n");
		pvar->nb_geno = 1;
		init_tab_geno (geno_seq, pvar);
		return 0; // retour de la position 0
	}
	// Sinon recherche si le genotype existe déjà dans la liste
	for (i = 0; i < pvar->nb_geno; i++)
	{ 
		if ( strcmp(pvar->tab_geno[i].sequence, geno_seq) == 0 ) // si chaines identiques
		{
			///printf("Le genotype existe déjà en position %d\n", i);
			pvar->tab_geno[i].nb_ind++;	
			return i; // retour de la position ou le genotype a été trouvé
		}
	}
	// Sinon stockage du genotype dans une nouvelle case à la suite du tableau
	///printf("Ajout d'un nouveau génotype en position %d\n", i);
	pvar->nb_geno = i+1;
	extend_tab_geno (geno_seq, pvar);
	return i; // retour de la position ou le genotype a été ajouté
}

/**** ajouter_tab_haplo ***********************************************/
int ajouter_tab_haplo (char* haplo_seq, T_info* pvar)
//char* haplo_seq, int* p_nb_haplo, int* p_haplo_num, T_haplo** p_tab_haplo
{
	int i;
	// Si la matrice est vide = ajout du premier haplotype
	if (pvar->tab_haplo == NULL)
	{ 
		///printf("Ajout du premier haplotype en position 0\n");
		pvar->nb_haplo = 1;
		init_tab_haplo (haplo_seq, pvar);
		return 0; // retour de la position 0
	}
	// Sinon recherche si l'haplotype existe déjà dans la liste
	for (i = 0; i < pvar->nb_haplo; i++)
	{
		if ( strcmp(pvar->tab_haplo[i].sequence, haplo_seq) == 0 ) // si chaines identiques
		{
			///printf("L'haplotype existe déjà en position %d\n", i);
			pvar->tab_haplo[i].nb_geno_expl++;
			return i; // retour de la position ou l'haplotype a été trouvé
		}
	}
	// Sinon stockage de l'haplotype dans une nouvelle case à la suite du tableau
	///printf("Ajout d'un nouvel haplotype en position %d\n", i);
	pvar->nb_haplo = i+1;
	extend_tab_haplo (haplo_seq, pvar);
	return i;
}

///** init_tab_geno ***************************************************/
void init_tab_geno (char* geno_seq, T_info* pvar)
//char* geno_seq, T_geno** p_tab_geno
{
	pvar->tab_geno = malloc (sizeof (T_geno)); // taille = 1 pour le premier element
	
	if (pvar->tab_geno == NULL)
	{
		fprintf (stderr, "Allocation impossible\n\n");
		exit (EXIT_FAILURE);
	}
	
	pvar->tab_geno[0].sequence = geno_seq;
	pvar->tab_geno[0].proba = 0;
	pvar->tab_geno[0].nb_ind = 1;
	pvar->tab_geno[0].nb_diplo_expl = 0;
	pvar->tab_geno[0].tete = NULL;
	pvar->tab_geno[0].queue = NULL;
	return;	
}

/**** init_tab_haplo **************************************************/
void init_tab_haplo (char* haplo_seq, T_info* pvar)
//char* haplo_seq, T_haplo** p_tab_haplo
{
	pvar->tab_haplo = malloc (sizeof (T_haplo)); // taille = 1 pour le premier element
	
	if (pvar->tab_haplo == NULL)
	{
		fprintf (stderr, "Allocation impossible\n\n");
		exit (EXIT_FAILURE);
	}
	
	pvar->tab_haplo[0].sequence = haplo_seq;
	pvar->tab_haplo[0].frequence = 0;
	pvar->tab_haplo[0].nb_geno_expl = 1;
	pvar->tab_haplo[0].tete = NULL;
	pvar->tab_haplo[0].queue = NULL;
	return;
}
 
///** extend_tab_geno *************************************************/
void extend_tab_geno (char* geno_seq, T_info* pvar)
//int nb_geno, char* geno_seq, T_geno** p_tab_geno
 {
	int n = pvar->nb_geno; // Nombre de génotypes dans la tab_geno
	
	pvar->tab_geno = realloc (pvar->tab_geno, sizeof (T_geno) * n);
	
	if (pvar->tab_geno == NULL)
	{
		fprintf (stderr, "Allocation impossible\n\n");
		exit (EXIT_FAILURE);
	}
	
	pvar->tab_geno[n-1].sequence = geno_seq;
	pvar->tab_geno[n-1].proba = 0;
	pvar->tab_geno[n-1].nb_ind = 1;
	pvar->tab_geno[n-1].nb_diplo_expl = 0;
	pvar->tab_geno[n-1].tete = NULL;
	pvar->tab_geno[n-1].queue = NULL;
	return;
} 

/**** extend_tab_haplo ************************************************/
void extend_tab_haplo (char* haplo_seq, T_info* pvar)
//int nb_haplo, char* haplo_seq, T_haplo** p_tab_haplo
{
	int n = pvar->nb_haplo; // Nombre de génotypes dans la tab_geno
	
	pvar->tab_haplo = realloc (pvar->tab_haplo , sizeof (T_haplo) * n);
	
	if (pvar->tab_haplo == NULL)
	{
		fprintf (stderr, "Allocation impossible\n\n");
		exit (EXIT_FAILURE);
	}
	
	pvar->tab_haplo[n-1].sequence = haplo_seq;
	pvar->tab_haplo[n-1].frequence = 0;
	pvar->tab_haplo[n-1].nb_geno_expl = 1;
	pvar->tab_haplo[n-1].tete = NULL;
	pvar->tab_haplo[n-1].queue = NULL;
	return;
}

/**** haplotypes_possibles ********************************************/
void haplotypes_possibles (char* geno_seq, T_info* pvar)
//char* genotype, int taille, char*** p_tab_haplo, int* p_nb_haplo
{
	int nb_amb, amb; // nombre total d'ambiguité dans la séquence et compteur d'ambiguité courante
	int change_prog; // variable booléenne pour remplir le tableau des haplotype possible en suivant alternativement 1/0 ou 0/1
	int val_basc, basc; //valeur à attendre pour basculer change_prog et compteur de basculement
	int j, k; // variables de contrôle de boucle
	int n; // compteur d'haplotypes générable pour le génotype donné (a retourner par p_nb_haplo_expl)
	char** tab; // tableau permettant de stocker temporairement les haplo générés
	
	nb_amb = compte_ambiguites(geno_seq, pvar->taille); // calcul du nombre d'ambiguités dans le genotype courant
	///printf("Nombre d'ambiguités : %d\n", nb_amb);
	
	if (nb_amb != 0) n = exp2(nb_amb); // il existe 2^amb haplotypes possibles
	else n = 2; // si 0 ambiguité, il faut quand même générer 2 haplotypes
	 
	tab = create_char_mat (n, pvar->taille+1);
	amb = 0;
	
	for (j = 0; j < pvar->taille ; j ++)
	{
		if (geno_seq[j] == '0') // Garnissage de tous les champs de la table haplo à cette même position avec 0
			for (k = 0; k < n ; k ++)
				tab[k][j] = '0';

		else if (geno_seq[j] == '2') // Garnissage tous les champs de la table haplo à cette même position avec 0
			for (k = 0; k < n ; k ++)
				tab[k][j] = '1';

		else  // cas le plus complexe de position ambigue
		{
			change_prog = 0; // initialisation du booléen
			val_basc = exp2(nb_amb)/exp2(amb++); // Nb d'haplotype aprés lequel il faut changer le sens de remplissage
			for (k = 0, basc = 0 ; k < n ; k +=2, basc +=2)  // k = compteur de boucle
			{
				if (basc == val_basc) { // Si le compteur de basculement atteint la valeur de bascule
					change_prog = !change_prog ; // bascule du booléen
					basc = 0 ; // reinitialisation du compteur de basculement
				}
				if (change_prog == 0) { // remplissage sens 0/1
					tab[k][j] = '0';
					tab[k+1][j] = '1';
				}
				else { // remplissage sens 1/0
					tab[k][j] = '1';
					tab[k+1][j] = '0';
				}
			}
		}
	}
	// Pour terminer proprement les chaines de charactères
	for (k = 0; k < n ; k ++)
		tab[k][j] = '\0';
	// Retour par remplissage des variable de T_info var
	pvar->tab_haplo_expl = tab;
	pvar->nb_haplo_expl = n;
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

/**** liberer_char_mat *************************************************/
 void liberer_char_mat (char** tab, int line)
{
	int i;

	for (i = 0; i < line; i++)
			free(tab[i]);
			
	free(tab);
		return;
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

/**** compte_ambiguites ***********************************************/
int compte_ambiguites (char* geno_seq, int taille)
//char* geno_seq, int taille
{
	int i;
	int nb_amb = 0;
	
	for (i = 0 ; i < taille; i++)
		if (geno_seq[i] == '1')
			nb_amb ++;
	return nb_amb;
}
 
 
///** ajouter_geno_a_individu *****************************************/
///void ajouter_geno_a_individu (int num_geno, T_info* pvar)
// int num_geno, T_individu* tab_individus
// A CODER

/**** ajouter_diplo_a_geno ********************************************/
void ajouter_diplo_a_geno (int num_haplo_A, int num_haplo_B, int num_geno, T_info* pvar)
//int num_haplo_A, int num_haplo_B, int num_geno, T_geno** p_tab_geno
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

	if (pvar->tab_geno[num_geno].tete == NULL) // ajout de la première cellule à la liste
	{
		pvar->tab_geno[num_geno].tete = cellule;
		pvar->tab_geno[num_geno].queue = cellule;
	}
	else // ajout en queue de liste = chainage simplifié
	{
		pvar->tab_geno[num_geno].queue -> suivant = cellule; // chainage avec cellule precedente
		pvar->tab_geno[num_geno].queue = cellule; // avancée du pointeur
	}
	return;
}
 
/**** ajouter_geno_a_haplo ********************************************/
void ajouter_geno_a_haplo (int num_haplo_principal, int num_haplo_compl, int num_geno_expl, T_info* pvar)
//int num_haplo_principal, int num_haplo_compl, int num_geno_expl, T_haplo** p_tab_haplo
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

	if (pvar->tab_haplo[num_haplo_principal].tete == NULL) // Ajout de la première cellule à la liste
	{
		pvar->tab_haplo[num_haplo_principal].tete = cellule;
		pvar->tab_haplo[num_haplo_principal].queue = cellule;
	}
	else // Ajout en queue de liste = chainage simplifié
	{
		pvar->tab_haplo[num_haplo_principal].queue -> suivant = cellule; // Chainage avec cellule precedente
		pvar->tab_haplo[num_haplo_principal].queue = cellule; // Avancée du pointeur
	}
	return;
}
 
/**** print_tab_haplo *************************************************/
void print_tab_haplo (T_info* pvar)
{
	int i;
	T_geno_expl* ptrj = NULL;
	
	printf("\n\n###############################################\n");
	printf("\nListe Des haplotypes et des genotypes expliqués\n");
	printf("\n###############################################\n\n");
	
	for (i = 0; i < pvar->nb_haplo; i++ )
	{
		ptrj = pvar->tab_haplo[i].tete;
		printf("Haplotype #%d \t Séquence : %s \t Fréquence %f \t Nombre de genotype(s) expliqué(s) : %d \n",
			i,
			pvar->tab_haplo[i].sequence,
			pvar->tab_haplo[i].frequence,
			pvar->tab_haplo[i].nb_geno_expl);
			
		printf("Liste des génotypes expliqués\n");
		while (ptrj)
		{
			printf("Genotype # %d (%s) avec Haplotype # %d (%s)\n", 
				ptrj -> num_geno_expl,
				pvar->tab_geno[ptrj->num_geno_expl].sequence,
				ptrj -> num_haplo_compl,
				pvar->tab_haplo[ptrj->num_haplo_compl].sequence );
						
			ptrj = ptrj -> suivant;
		}
		printf("\n");
	}
}

/**** print_tab_geno *************************************************/
void print_tab_geno (T_info* pvar)
{
	int i;
	T_diplo_expl* ptrj = NULL;
	
	printf("\n\n#################################################\n");
	printf("\nListe des genotypes et des diplotypes explicatifs\n");
	printf("\n#################################################\n\n");
	
	for (i = 0; i < pvar->nb_geno; i++ )
	{
		ptrj = pvar->tab_geno[i].tete;
		printf("Haplotype #%d\t Séquence : %s\t Probabilité %f\t Nombre d'individu(s) concerné(s) : %d\t Nombre de diplotype(s) explicatif(s) : %d\n",
			i,
			pvar->tab_geno[i].sequence,
			pvar->tab_geno[i].proba,
			pvar->tab_geno[i].nb_ind,
			pvar->tab_geno[i].nb_diplo_expl);
				
		printf("Liste des diplotypes explicatifs\n");
		while (ptrj)
		{
			printf("Haplotype # %d (%s) avec Haplotype # %d (%s)\n", 
				ptrj->num_haplo_A,
				pvar->tab_haplo[ptrj->num_haplo_A].sequence,
				ptrj->num_haplo_B,
				pvar->tab_haplo[ptrj->num_haplo_A].sequence );
						
			ptrj = ptrj -> suivant;
		}
		printf("\n");
	}
}
