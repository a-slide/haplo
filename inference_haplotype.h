#ifndef INFERENCE_HAPLOTYPE_H
#define INFERENCE_HAPLOTYPE_H


	/***********************************************************************
	 * STRUCTURES
	 **********************************************************************/
	
	// Structure permettant de stocker le génotype d'un individu et le numéro de génotype de référence
	typedef struct individu T_individu;
	struct individu
	{
		char* sequence; // Séquence de l'individu
		int num_geno; // Numéro du génotype
	};
	
	// Structure permettant de stocker un diplotype explicatif d'un génotype 
	typedef struct diplo_expl T_diplo_expl;
	struct diplo_expl
	{
		int num_haplo_A; // paire d'haplotype explicatifs = diplotype
		int num_haplo_B;
		T_diplo_expl* suivant; // Pointeur vers le diplotype explicatif suivant
	};
	
	// Structure permettant de stocker un génotype, sa probabilité, son nombre d'individus expliqués...
	typedef struct geno T_geno;// tête de liste de type T_diplo_expl
	struct geno
	{
		char* sequence; // Séquence du génotype
		double proba; // Probabilité du génotype à l'itération i de la boucle EM
		double proba_prec; // Probabilité du génotype à l'itération i-1 de la boucle EM
		int nb_ind; // compteur d'individus avec ce génotype
		int nb_diplo_expl; // compteur de diplotype explicatifs
		int num_haplo_A_max; // numéro de l'haplotype A de la paire la plus explicative
		int num_haplo_B_max; // numéro de l'haplotype B de la paire la plus explicative
		double proba_diplo_max; // probabilité du diplotype le plus explicatif
		T_diplo_expl* tete; // Pointeur vers la tête de la liste chaînée de diplotypes explicatifs
		T_diplo_expl* queue; // Pointeur vers la queue de la liste chaînée de diplotypes explicatifs
	};
	
	// Structure permettant de stocker un génotype expliqué par un haplotype avec l'haplotype complémentaire 
	typedef struct geno_expl T_geno_expl;
	struct geno_expl
	{
		int num_geno_expl; // génotype expliqué
		int num_haplo_compl; // haplotype complémentaire pour expliquer le génotype concerné
		T_geno_expl* suivant;
	};
	
	// Structure permettant de stocker un haplotype, sa fréquence, le nombre de génotypes expliqués par cet haplotype....
	typedef struct haplo T_haplo;// tête de liste de type T_geno_expl
	struct haplo
	{
		char* sequence; // Séquence de l'haplotype
		double frequence; // Fréquence de l'haplotype
		double frequence_prec; // Fréquence de l'haplotype à l'itération i-1 de la boucle EM
		int nb_geno_expl; // Nombre de génotypes expliqués par l'haplotype
		int nb_haplo; 
		T_geno_expl* tete; // Pointeur vers la tête de la liste chaînée de génotypes explicatifs
		T_geno_expl* queue; // Pointeur vers la tête de la liste chaînée de génotypes explicatifs
	};
	
	// Structure contenant les variables et tableaux importants facilitant le passage de paramètres
	typedef struct info T_info;
	struct info
	{
		int taille; // Taille des génotypes et haplotypes
		int nb_ind; // Nombre d'individus
		int nb_geno; // Nombre de génotypes
		int nb_haplo; // Nombre d'haplotypes
		int nb_haplo_expl; // Nombre d'haplotypes explicatifs
		T_individu* tab_individus; // Tableau d'individus
		T_geno* tab_geno; // Tableau de génotypes
		T_haplo* tab_haplo; // Tableau d'haplotypes
		char** tab_haplo_expl; // Tableau d'haplotypes explicatifs
	};

	/***********************************************************************
	 * PROTOTYPES
	 **********************************************************************/

	// Dans main
	void usage (char*);
	void importation_genotypes(char*, T_info*);
	void preparation_liste_geno_haplo (T_info*);
	void initialisation_freq_proba (T_info*, int);

	//Dans importation_genotypes
	int nb_ligne (FILE*);
	int nb_char (FILE*);

	// Dans preparation_liste_geno_haplo
	int ajouter_tab_geno (char*, T_info*);
	void init_tab_geno (char*, T_info*);
	void extend_tab_geno (char*, T_info*);
	void haplotypes_possibles (char*, T_info*);
	void print_string_table (char**, int);
	int compte_ambiguites (char*, int);
	int ajouter_tab_haplo (char*, T_info*);
	void init_tab_haplo (char*, T_info*);
	void extend_tab_haplo (char*, T_info*);
	void ajouter_diplo_a_geno (int, int, int, T_info*);
	void ajouter_geno_a_haplo (int, int, int, T_info*);

	// Dans initialisation_freq_proba
	void init_haplo_equi_freq (T_info*);
	void haplo_random_freq (T_info*);
	void init_geno_proba (T_info*);
	void print_tab_haplo (T_info*);
	void print_tab_geno (T_info*);
	
	// Dans maximisation_estimation
	void maximisation (T_info*);
	double estimation_esperance(T_info*);
	void update_proba_freq (T_info*);
	
	// Dans traitement_final
	void diplotype_plus_probable (T_info*);
	void export_geno_diplo (T_info*);
	int comparaison_frequence (void const*, void const*);
	void export_haplo (T_info*);

#endif /* INFERENCE_HAPLOTYPE_H */
