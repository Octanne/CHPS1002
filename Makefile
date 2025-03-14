CC=gfortran
#CC=ifort

SRC_DIR = sources
OBJ_DIR = obj
BIN_DIR = bin

help:
	@echo "(C) JC.Boisson"
	@echo "Sous-commandes :"
	@echo "NomFichier       : compile et fait l'édition de lien du fichier \"NomFichier.f90\" correspondant pour générer un exécutable"
	@echo "make cleanSource : supprime les fichiers parasites (*~, *.old,#*,*.bak)"
	@echo "make clean       : supprime *tous* les fichiers reproductibles ici les.o et aussi les fichiers parasites"
	@echo "make clean_all   : supprime *tous* les fichiers reproductibles, les fichiers parasites et aussi les exécutables"

$(BIN_DIR)/lecteur_mol2: $(OBJ_DIR)/lecteur_mol2.o
	-@echo ""
	-@echo "Linking    $(@)"
	-@echo ""
	-@$(CC) -o $@ $+

$(BIN_DIR)/chargeur_covalence: $(OBJ_DIR)/chargeur_covalence.o
	-@echo ""
	-@echo "Linking    $(@)"
	-@echo ""
	-@$(CC) -o $@ $+

$(BIN_DIR)/affiche_topologie: $(OBJ_DIR)/affiche_topologie.o
	-@echo ""
	-@echo "Linking    $(@)"
	-@echo ""
	-@$(CC) -o $@ $+

$(BIN_DIR)/complete_mol2: $(OBJ_DIR)/complete_mol2.o
	-@echo ""
	-@echo "Linking    $(@)"
	-@echo ""
	-@$(CC) -o $@ $+

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90 
	-@echo ""
	-@echo "Generating $@"
	-@mkdir -p $(OBJ_DIR)
	-@$(CC) -c $< -o $@

EXEC = $(BIN_DIR)/lecteur_mol2
EXEC+= $(BIN_DIR)/chargeur_covalence
EXEC+= $(BIN_DIR)/affiche_topologie
EXEC+= $(BIN_DIR)/complete_mol2

###------------------------------
### Cleaning
###------------------------------------------------------------

clean:
	-@rm -rf $(OBJ_DIR)/*.o

clean_all: clean cleanSource
	-@rm -rf $(EXEC)

cleanSource:
	-@find . \( -name "*~" -o -name "*.old" -o -name "#*" \) -print -exec rm \{\} \;

.PHONY:  $(EXEC) clean clean_all cleanSource
