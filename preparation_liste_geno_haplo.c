#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "inference_haplotype.h"
#include "ptr_allocation.h"
 
/**** preparation_liste_geno_haplo ***************************************/
// PRE CONDITIONS
// POST CONDITIONS
void preparation_liste_geno_haplo (T_info* pvar)
{
	int i, j;
	int num_haplo_A, num_haplo_B; // Indices d'un genotype et d'une paire d'haplotypes explicatifs

	// Création d'une liste de genotypes non redondants à partir de tab_individus
	
	printf("\nListe des individus et numero de genotype correspondant\n");
	for (i = 0 ; i < pvar->nb_ind; i++)
	{
		pvar->tab_individus[i].num_geno = ajouter_tab_geno (pvar->tab_individus[i].sequence, pvar);
		printf("Individu #%d\tSequence: %s\tGenotype %d\n",i, pvar->tab_individus[i].sequence, pvar->tab_individus[i].num_geno);
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
			
			// Mettre les genotypes explicatifs avec l'haplo complementaire dans la liste chainée associée à chacun des 2 haplotypes
			ajouter_geno_a_haplo (num_haplo_A, num_haplo_B, i, pvar);
			ajouter_geno_a_haplo (num_haplo_B, num_haplo_A, i, pvar);
			
			//Mettre les paires d'haplotypes dans la liste chainée associée au génotype courant
			ajouter_diplo_a_geno (num_haplo_A, num_haplo_B, i, pvar);
		}
	//	free_char_mat (pvar->tab_haplo_expl, pvar->nb_haplo_expl);
	}
	
	return;
}
 
///** ajouter_tab_geno  ***********************************************/
int ajouter_tab_geno (char* geno_seq, T_info* pvar)
// PRE CONDITIONS
// POST CONDITIONS
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
// PRE CONDITIONS
// POST CONDITIONS
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
{
	pvar->tab_geno = malloc (sizeof (T_geno)); // taille = 1 pour le premier element
	if (pvar->tab_geno == NULL) error_and_exit();
	
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
{
	pvar->tab_haplo = malloc (sizeof (T_haplo)); // taille = 1 pour le premier element
	if (pvar->tab_haplo == NULL) error_and_exit();
	
	pvar->tab_haplo[0].sequence = haplo_seq;
	pvar->tab_haplo[0].frequence = 0;
	pvar->tab_haplo[0].nb_geno_expl = 1;
	pvar->tab_haplo[0].tete = NULL;
	pvar->tab_haplo[0].queue = NULL;
	return;
}
 
///** extend_tab_geno *************************************************/
void extend_tab_geno (char* geno_seq, T_info* pvar)
 {
	int n = pvar->nb_geno; // Nombre de génotypes dans la tab_geno
	
	pvar->tab_geno = realloc (pvar->tab_geno, sizeof (T_geno) * n);
	if (pvar->tab_geno == NULL) error_and_exit();
	
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
{
	int n = pvar->nb_haplo; // Nombre de génotypes dans la tab_geno
	
	pvar->tab_haplo = realloc (pvar->tab_haplo , sizeof (T_haplo) * n);
	if (pvar->tab_haplo == NULL) error_and_exit();
	
	pvar->tab_haplo[n-1].sequence = haplo_seq;
	pvar->tab_haplo[n-1].frequence = 0;
	pvar->tab_haplo[n-1].nb_geno_expl = 1;
	pvar->tab_haplo[n-1].tete = NULL;
	pvar->tab_haplo[n-1].queue = NULL;
	return;
}

/**** haplotypes_possibles ********************************************/
void haplotypes_possibles (char* geno_seq, T_info* pvar)
// PRE CONDITIONS
// POST CONDITIONS
{
	int nb_amb, amb; // nombre total d'ambiguité dans la séquence et compteur d'ambiguité courante
	int change_prog; // variable booléenne pour remplir le tableau des haplotype possible en suivant alternativement 1/0 ou 0/1
	int val_basc, basc; //valeur à attendre pour basculer change_prog et compteur de basculement
	int j, k; // variables de contrôle de boucle
	int n; // compteur d'haplotypes générable pour le génotype donné (a retourner par p_nb_haplo_expl)
	char** tab; // tableau permettant de stocker temporairement les haplo générés
	
	// calcul du nombre d'ambiguités dans le genotype courant
	nb_amb = compte_ambiguites(geno_seq, pvar->taille); 
	
	// il existe 2^nb_amb haplotypes possibles mais si nb_amb = 0 il faut quand même generer 2 haplotypes
	n = ((nb_amb == 0) ? 2 : exp2(nb_amb)); 

	tab = malloc_char_mat(n, pvar->taille+1);
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
	// Remplissage des variables de T_info var
	pvar->tab_haplo_expl = tab;
	pvar->nb_haplo_expl = n;
	
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
{
	int i;
	int nb_amb = 0;
	
	for (i = 0 ; i < taille; i++)
		if (geno_seq[i] == '1')
			nb_amb ++;
	return nb_amb;
}
 
/**** ajouter_diplo_a_geno ********************************************/
void ajouter_diplo_a_geno (int num_haplo_A, int num_haplo_B, int num_geno, T_info* pvar)
// PRE CONDITIONS
// POST CONDITIONS
{
	T_diplo_expl* cellule = malloc( sizeof (T_diplo_expl));
	if (cellule == NULL) error_and_exit();
	
	cellule -> num_haplo_A = num_haplo_A;
	cellule -> num_haplo_B = num_haplo_B;
	cellule -> suivant = NULL;

	if (pvar->tab_geno[num_geno].tete == NULL) // Ajout de la première cellule à la liste
	{
		pvar->tab_geno[num_geno].tete = cellule;
		pvar->tab_geno[num_geno].queue = cellule;
	}
	else // Ajout en queue de liste = chainage simplifié
	{
		pvar->tab_geno[num_geno].queue -> suivant = cellule; // Chainage avec la cellule precedente
		pvar->tab_geno[num_geno].queue = cellule; // Avancée du pointeur
	}
	return;
}
 
/**** ajouter_geno_a_haplo ********************************************/
void ajouter_geno_a_haplo (int num_haplo_principal, int num_haplo_compl, int num_geno_expl, T_info* pvar)
// PRE CONDITIONS
// POST CONDITIONS
{
	T_geno_expl* cellule = malloc( sizeof (T_geno_expl));
	if (cellule == NULL) error_and_exit();
	
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
		pvar->tab_haplo[num_haplo_principal].queue -> suivant = cellule; // Chainage avec la cellule precedente
		pvar->tab_haplo[num_haplo_principal].queue = cellule; // Avancée du pointeur
	}
	return;
}

