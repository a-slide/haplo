############################# VARIABLES ################################

CC = gcc
	# Compilateur
	
CFLAGS = -g -Wall -W -Wno-unused-parameter -lm
	# Options de compilation
	# -Wall -W et -Wno-unused-parameter pour apporter des informations

LDFLAGS = 
	# Options d'edition de lien

BIN = EM_Haplotype_inference
	# Nom des executables a creer

SRC = $(wildcard *.c)
	# Liste automatique des fichiers sources
	
OBJ = $(SRC:.c=.o)
	# Liste automatique des fichiers objects

HEAD = $(SRC:.c=.h)
	# Liste automatique des fichiers headers

#################### INSTRUCTIONS DE COMPILATION #######################
# $@ =  Cible # $^ = liste des dépendances # $< Première dépendance #

all: $(BIN)
	# Ensemble des executables à produire

EM_Haplotype_inference: $(OBJ)
	$(CC) $^ $(LDFLAGS) -o $@
	# Edition de lien a partir des fichiers objets

EM_main.o: $(HEAD) EM_structures.h
	# Si un des headers a changé alors recompiler EM_main.o
	
%.o: %.c %.h
	$(CC) -c $< $(CFLAGS) -o $@ 
	# Règle générique de compilation de tt les sources en objets	

##################### INSTRUCTIONS DE NETTOYAGE ########################

.PHONY: clean mrproper
	# PHONY = Dependances systematiquement reconstruites

clean:
	rm -rf $(OBJ)
	# Supprimer tous les fichiers intermédiaires

mrproper: clean
	rm -rf $(BIN) $(OBJ)
	# Supprimer tout ce qui peut être régénéré et reconstruit complètement
