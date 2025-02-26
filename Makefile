CFLAGS = -I/usr/local/include -Wall -O0 -Wextra
LDFLAGS = -L/usr/local/lib -lflint -lgmp -lmpfr
CC = gcc
SRC_DIR = src
OBJ_DIR = obj

TARGET = model

SRC_FILES = $(shell find $(SRC_DIR) -name "*.c")
OBJ_FILES = $(SRC_FILES:$(SRC_DIR)/%.c=$(OBJ_DIR)/%.o)

all: $(TARGET)

$(TARGET): $(OBJ_FILES)
	$(CC) $(CFLAGS) $(OBJ_FILES) $(LDFLAGS) -o $(TARGET)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	@mkdir -p $(OBJ_DIR)  # Créer le dossier obj s'il n'existe pas
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -rf $(TARGET) $(OBJ_DIR)
