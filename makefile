############################# VARIABLES ################################

CC = gcc
	# Compilateur
	
CFLAGS =  -c -g -Wall -W -Wno-unused-parameter
	# Options de compilation
	# -Wall -W et -Wno-unused-parameter pour apporter des informations

LDFLAGS = -lm
	# Options d'edition de lien

BIN = EM_main
	# Nom des executables a creer

#################### INSTRUCTIONS DE COMPILATION #######################
# $@ =  Cible # $^ = liste des dépendances # $< Première dépendance #


all: $(BIN)
	
EM_main: EM_main.o importation_genotypes.o preparer_liste_geno_haplo.o
	$(CC) $^ $(LDFLAGS) -o $@
	# Edition de lien a partir des fichiers objets

EM_main.o: EM_main.c EM_main.h
	$(CC) $< $(CFLAGS) -o $@ 
	# Compilation de EM_main.c

importation_genotypes.o: importation_genotypes.c EM_main.h
	$(CC) $< $(CFLAGS) -o $@
	# Compilation de importation_genotypes.c

preparer_liste_geno_haplo.o: preparer_liste_geno_haplo.c EM_main.h
	$(CC) $< $(CFLAGS) -o $@
	# Compilation de preparer_liste_geno_haplo.c

##################### INSTRUCTIONS DE NETTOYAGE ########################

.PHONY: clean mrproper
	# PHONY = Dependances systematiquement reconstruites

clean:
	rm -rf EM_main.o importation_genotypes.o preparer_liste_geno_haplo.o
	# Supprimer tous les fichiers intermédiaires

mrproper: clean
	rm -rf $(BIN)
	# Supprimer tout ce qui peut être régénéré et reconstruit complètement
