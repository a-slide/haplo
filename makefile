############################# VARIABLES ################################

CC = gcc
	# Compilateur
	
CFLAGS =  -c -g -Wall -W -Wno-unused-parameter
	# Options de compilation
	# -Wall -W et -Wno-unused-parameter pour apporter des informations

LDFLAGS = -lm
	# Options d'edition de lien

BIN = inference_haplotype
	# Nom des executables a creer

#################### INSTRUCTIONS DE COMPILATION #######################
# $@ =  Cible # $^ = liste des dépendances # $< Première dépendance #

all: $(BIN)
	# Compilation terminée
	
inference_haplotype: inference_haplotype.o importation_genotypes.o preparation_liste_geno_haplo.o initialisation_freq_proba.o maximisation_esperance.o ptr_allocation.o traitement_final.o
	$(CC) $^ $(LDFLAGS) -o $@
	# Edition de lien a partir des fichiers objets

inference_haplotype.o: inference_haplotype.c inference_haplotype.h ptr_allocation.h
	$(CC) $< $(CFLAGS) -o $@ 
	# Compilation de inference_haplotype.c

importation_genotypes.o: importation_genotypes.c inference_haplotype.h ptr_allocation.h
	$(CC) $< $(CFLAGS) -o $@
	# Compilation de importation_genotypes.c

preparation_liste_geno_haplo.o: preparation_liste_geno_haplo.c inference_haplotype.h ptr_allocation.h
	$(CC) $< $(CFLAGS) -o $@
	# Compilation de preparation_liste_geno_haplo.c

initialisation_freq_proba.o: initialisation_freq_proba.c inference_haplotype.h
	$(CC) $< $(CFLAGS) -o $@
	# Compilation de initialisation_freq_proba.c

maximisation_esperance.o: maximisation_esperance.c inference_haplotype.h
	$(CC) $< $(CFLAGS) -o $@
	# Compilation de maximisation_esperance.c
	
traitement_final.o: traitement_final.c inference_haplotype.h
	$(CC) $< $(CFLAGS) -o $@
	# Compilation de traitement_final.c

ptr_allocation.o: ptr_allocation.c ptr_allocation.h
	$(CC) $< $(CFLAGS) -o $@
	# Compilation de ptr_allocation.c

##################### INSTRUCTIONS DE NETTOYAGE ########################

.PHONY: clean mrproper
	# PHONY = Dependances systematiquement reconstruites

clean:
	rm -rf inference_haplotype.o importation_genotypes.o preparation_liste_geno_haplo.o initialisation_freq_proba.o maximisation_esperance.o ptr_allocation.o
	# Supprimer tous les fichiers intermédiaires

mrproper: clean
	rm -rf $(BIN)
	# Supprimer tout ce qui peut être régénéré et reconstruit complètement
