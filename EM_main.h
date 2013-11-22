/***********************************************************************
 * STRUCTURES
 **********************************************************************/

typedef struct individu T_individu;
struct individu
{
	char* sequence;
	int num_geno;
};

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
	double proba;
	double proba_prec;
	int nb_ind; // compteur d'individus avec ce genotype
	int nb_diplo_expl; // compteur de diplotype explicatifs
	T_diplo_expl* tete;
	T_diplo_expl* queue;
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
	double frequence;
	double frequence_prec;
	int nb_geno_expl;
	int nb_haplo;
	T_geno_expl* tete;
	T_geno_expl* queue;
};

typedef struct info T_info;// structure contenant les variables et tableaux importants
struct info
{
	int taille;
	int nb_ind;
	int nb_geno;
	int nb_haplo;
	int nb_haplo_expl;  //Plutot nombre géno explicatifs ?
	T_individu* tab_individus;
	T_geno* tab_geno;
	T_haplo* tab_haplo;
	char** tab_haplo_expl;
	double vraisemblance_prec;
	double vraisemblance;
};

/***********************************************************************
 * PROTOTYPES
 **********************************************************************/

// Dans main
void usage (char*);
void importation_genotypes(char*, T_info*);
void preparation_liste_geno_haplo (T_info*);
void initialisation_freq_proba (T_info*, char);
void Maximisation_et_Esperance(T_info*);

///expectation;
///maximisation;

//Dans importation_genotypes 
FILE* init_file_ptr (char*, char*);
int nb_ligne (FILE*);
int nb_char (FILE*);

// Dans preparation_liste_geno_haplo
int ajouter_tab_geno (char*, T_info*);
void init_tab_geno (char*, T_info*);
void extend_tab_geno (char*, T_info*);
void haplotypes_possibles (char*, T_info*);
char** create_char_mat (int, int);
void liberer_char_mat (char**, int);
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

// Dans Maximisation
void Maximisation (T_info*);
double Estimation_Esperance(T_info*);
void Update_Hfreqpreq_Gprobaprec_vraisemblance_preq (T_info*);
