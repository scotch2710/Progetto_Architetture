#include <stdio.h>
#include <stdlib.h>

void readDs2File(const char* filename) {
    FILE* file = fopen(filename, "rb");
    if (!file) {
        perror("Errore nell'apertura del file");
        return;
    }

    int firstInt, secondInt;
    double value;
    size_t count = 0;

    // Legge i primi due interi
    if (fread(&firstInt, sizeof(int), 1, file) != 1) {
        perror("Errore nella lettura del primo intero");
        fclose(file);
        return;
    }
    if (fread(&secondInt, sizeof(int), 1, file) != 1) {
        perror("Errore nella lettura del secondo intero");
        fclose(file);
        return;
    }

    printf("Primo intero: %d\n", firstInt);
    printf("Secondo intero: %d\n", secondInt);

    // Legge i double successivi
    printf("Vettori di double:\n");
    while (fread(&value, sizeof(double), 1, file) == 1) {
        printf("Elemento %zu: %lf\n", count, value);
        count++;
    }

    fclose(file);
    printf("Totale elementi double letti: %zu\n", count);
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        printf("Uso: %s <nome_file>\n", argv[0]);
        return 1;
    }

    readDs2File(argv[1]);
    return 0;
}
