############################# VARIABLES ################################

CC = gcc
	# Compilateur
	
CFLAGS = -g -Wall -W -Wno-unused-parameter -lm
	# Options de compilation
	# -Wall -W et -Wno-unused-parameter pour apporter des informations

LDFLAGS = 
	# Options d'edition de lien

BIN = EM_main
	# Nom des executables a creer
	
OBJ = $(SRC:.c=.o)
	# Liste automatique des fichiers objects


#################### INSTRUCTIONS DE COMPILATION #######################


all: $(BIN)
	# Ensemble des executables à produire

EM_main: EM_main.o importation_genotypes.o preparer_liste_geno_haplo.o
	$(CC) -o EM_main EM_main.o importation_genotypes.o preparer_liste_geno_haplo.o $(CFLAGS)
	# Edition de lien a partir des fichiers objets

EM_main.o: EM_main.c
	$(CC) -o EM_main.o -c EM_main.c $(CFLAGS)
	# Création du fichier binaire EM_main.o 

importation_genotypes.o: importation_genotypes.c EM_main.h
	$(CC) -o importation_genotypes.o -c importation_genotypes.c $(CFLAGS)
	# Création du fichier binaire importation_genotypes.o 

preparer_liste_geno_haplo.o: preparer_liste_geno_haplo.c EM_main.h
	$(CC) -o preparer_liste_geno_haplo.o -c preparer_liste_geno_haplo.c $(CFLAGS)
	# Création du fichier binaire preparer_liste_geno_haplo.o 

##################### INSTRUCTIONS DE NETTOYAGE ########################

.PHONY: clean mrproper
	# PHONY = Dependances systematiquement reconstruites

clean:
	rm -rf $(OBJ)
	# Supprimer tous les fichiers intermédiaires

mrproper: clean
	rm -rf $(BIN) $(OBJ)
	# Supprimer tout ce qui peut être régénéré et reconstruit complètement
