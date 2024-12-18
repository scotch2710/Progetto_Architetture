/**************************************************************************************
* 
* CdL Magistrale in Ingegneria Informatica
* Corso di Architetture e Programmazione dei Sistemi di Elaborazione - a.a. 2024/25
* 
* Progetto dell'algoritmo Predizione struttura terziaria proteine 221 231 a
* in linguaggio assembly x86-32 + SSE
* 
* F. Angiulli F. Fassetti S. Nisticò, novembre 2024
* 
**************************************************************************************/

/*
* 
* Software necessario per l'esecuzione:
* 
*    NASM (www.nasm.us)
*    GCC (gcc.gnu.org)
* 
* entrambi sono disponibili come pacchetti software 
* installabili mediante il packaging tool del sistema 
* operativo; per esempio, su Ubuntu, mediante i comandi:
* 
*    sudo apt-get install nasm
*    sudo apt-get install gcc
* 
* potrebbe essere necessario installare le seguenti librerie:
* 
*    sudo apt-get install lib32gcc-4.8-dev (o altra versione)
*    sudo apt-get install libc6-dev-i386
* 
* Per generare il file eseguibile:
* 
* nasm -f elf32 pst32.nasm && gcc -m32 -msse -O0 -no-pie sseutils32.o pst32.o pst32c.c -o pst32c -lm && ./pst32c $pars
* 
* oppure
* 
* ./runpst32
* 
*/

#define _USE_MATH_DEFINES
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <libgen.h>
#include <xmmintrin.h>

#define	type		float
#define	MATRIX		type*
#define	VECTOR		type*

#define random() (((type) rand())/RAND_MAX)

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

type hydrophobicity[] = {1.8, -1, 2.5, -3.5, -3.5, 2.8, -0.4, -3.2, 4.5, -1, -3.9, 3.8, 1.9, -3.5, -1, -1.6, -3.5, -4.5, -0.8, -0.7, -1, 4.2, -0.9, -1, -1.3, -1};		// hydrophobicity
type volume[] = {88.6, -1, 108.5, 111.1, 138.4, 189.9, 60.1, 153.2, 166.7, -1, 168.6, 166.7, 162.9, 114.1, -1, 112.7, 143.8, 173.4, 89.0, 116.1, -1, 140.0, 227.8, -1, 193.6, -1};		// volume
type charge[] = {0, -1, 0, -1, -1, 0, 0, 0.5, 0, -1, 1, 0, 0, 0, -1, 0, 0, 1, 0, 0, -1, 0, 0, -1, 0, -1};		// charge

typedef struct {
	char* seq;		// sequenza di amminoacidi
	int N;			// lunghezza sequenza
	unsigned int sd; 	// seed per la generazione casuale
	type to;		// temperatura INIZIALE
	type alpha;		// tasso di raffredamento
	type k;		// costante
	VECTOR hydrophobicity; // hydrophobicity
	VECTOR volume;		// volume
	VECTOR charge;		// charge
	VECTOR phi;		// vettore angoli phi
	VECTOR psi;		// vettore angoli psi
	type e;		// energy
	int display;
	int silent;

} params;


/*
* 
*	Le funzioni sono state scritte assumento che le matrici siano memorizzate 
* 	mediante un array (float*), in modo da occupare un unico blocco
* 	di memoria, ma a scelta del candidato possono essere 
* 	memorizzate mediante array di array (float**).
* 
* 	In entrambi i casi il candidato dovr� inoltre scegliere se memorizzare le
* 	matrici per righe (row-major order) o per colonne (column major-order).
*
* 	L'assunzione corrente � che le matrici siano in row-major order.
* 
*/

void* get_block(int size, int elements) { 
	return _mm_malloc(elements*size,16); 
}

void free_block(void* p) { 
	_mm_free(p);
}

MATRIX alloc_matrix(int rows, int cols) {
	return (MATRIX) get_block(sizeof(type),rows*cols);
}

int* alloc_int_matrix(int rows, int cols) {
	return (int*) get_block(sizeof(int),rows*cols);
}

char* alloc_char_matrix(int rows, int cols) {
	return (char*) get_block(sizeof(char),rows*cols);
}

void dealloc_matrix(void* mat) {
	free_block(mat);
}

/*
* 
* 	load_data
* 	=========
* 
*	Legge da file una matrice di N righe
* 	e M colonne e la memorizza in un array lineare in row-major order
* 
* 	Codifica del file:
* 	primi 4 byte: numero di righe (N) --> numero intero
* 	successivi 4 byte: numero di colonne (M) --> numero intero
* 	successivi N*M*4 byte: matrix data in row-major order --> numeri floating-point a precisione singola
* 
*****************************************************************************
*	Se lo si ritiene opportuno, � possibile cambiare la codifica in memoria
* 	della matrice. 
*****************************************************************************
* 
*/
MATRIX load_data(char* filename, int *n, int *k) {
	FILE* fp;
	int rows, cols, status, i;
	
	fp = fopen(filename, "rb");
	
	if (fp == NULL){
		printf("'%s': bad data file name!\n", filename);
		exit(0);
	}
	
	status = fread(&cols, sizeof(int), 1, fp);
	status = fread(&rows, sizeof(int), 1, fp);
	
	MATRIX data = alloc_matrix(rows,cols);
	status = fread(data, sizeof(type), rows*cols, fp);
	fclose(fp);
	
	*n = rows;
	*k = cols;
	
	return data;
}

/*
* 
* 	load_seq
* 	=========
* 
*	Legge da file una matrice di N righe
* 	e M colonne e la memorizza in un array lineare in row-major order
* 
* 	Codifica del file:
* 	primi 4 byte: numero di righe (N) --> numero intero
* 	successivi 4 byte: numero di colonne (M) --> numero intero
* 	successivi N*M*1 byte: matrix data in row-major order --> charatteri che compongono la stringa
* 
*****************************************************************************
*	Se lo si ritiene opportuno, � possibile cambiare la codifica in memoria
* 	della matrice. 
*****************************************************************************
* 
*/
char* load_seq(char* filename, int *n, int *k) {
	FILE* fp;
	int rows, cols, status, i;
	
	fp = fopen(filename, "rb");
	
	if (fp == NULL){
		printf("'%s': bad data file name!\n", filename);
		exit(0);
	}
	
	status = fread(&cols, sizeof(int), 1, fp);
	status = fread(&rows, sizeof(int), 1, fp);

	
	char* data = alloc_char_matrix(rows,cols);
	status = fread(data, sizeof(char), rows*cols, fp);
	fclose(fp);
	
	*n = rows;
	*k = cols;
	
	return data;
}

/*
* 	save_data
* 	=========
* 
*	Salva su file un array lineare in row-major order
*	come matrice di N righe e M colonne
* 
* 	Codifica del file:
* 	primi 4 byte: numero di righe (N) --> numero intero a 32 bit
* 	successivi 4 byte: numero di colonne (M) --> numero intero a 32 bit
* 	successivi N*M*4 byte: matrix data in row-major order --> numeri interi o floating-point a precisione singola
*/
void save_data(char* filename, void* X, int n, int k) {
	FILE* fp;
	int i;
	fp = fopen(filename, "wb");
	if(X != NULL){
		fwrite(&k, 4, 1, fp);
		fwrite(&n, 4, 1, fp);
		for (i = 0; i < n; i++) {
			fwrite(X, sizeof(type), k, fp);
			//printf("%i %i\n", ((int*)X)[0], ((int*)X)[1]);
			X += sizeof(type)*k;
		}
	}
	else{
		int x = 0;
		fwrite(&x, 4, 1, fp);
		fwrite(&x, 4, 1, fp);
	}
	fclose(fp);
}

/*
* 	save_out
* 	=========
* 
*	Salva su file un array lineare composto da k elementi.
* 
* 	Codifica del file:
* 	primi 4 byte: contenenti l'intero 1 		--> numero intero a 32 bit
* 	successivi 4 byte: numero di elementi k     --> numero intero a 32 bit
* 	successivi byte: elementi del vettore 		--> k numero floating-point a precisione singola
*/
void save_out(char* filename, MATRIX X, int k) {
	FILE* fp;
	int i;
	int n = 1;
	fp = fopen(filename, "wb");
	if(X != NULL){
		fwrite(&n, 4, 1, fp);
		fwrite(&k, 4, 1, fp);
		fwrite(X, sizeof(type), k, fp);
	}
	fclose(fp);
}

/*
* 	gen_rnd_mat
* 	=========
* 
*	Genera in maniera casuale numeri reali tra -pi e pi
*	per riempire una struttura dati di dimensione Nx1
* 
*/
void gen_rnd_mat(VECTOR v, int N){
	int i;

	for(i=0; i<N; i++){
		// Campionamento del valore + scalatura
		v[i] = (random() * 2 * M_PI) - M_PI;
	}
}

// PROCEDURE ASSEMBLY
extern void prova(params* input);

void vector_matrix_product(VECTOR v, MATRIX R, VECTOR result) {
    // Calcola il prodotto v * R
    for (int i = 0; i < 3; i++) {
        result[i] = v[0] * R[i] +
                    v[1] * R[i + 3] +
                    v[2] * R[i + 6];
    }
}

type approx_cos(type theta) {
    type x2 = theta * theta;
    return 1 - (x2 / 2.0) + (x2 * x2 / 24.0) - (x2 * x2 * x2 / 720.0);
}

type approx_sin(type theta) {
    type x2 = theta * theta;
    return theta - (x2 * theta / 6.0) + (x2 * x2 * theta / 120.0) - (x2 * x2 * x2 / 5040.0);

}

extern MATRIX rotation(VECTOR axis, type theta){
	type prod_scal= ((axis[0]*axis[0])+(axis[1]*axis[1])+(axis[2]*axis[2]));

	axis[0] = axis[0] / prod_scal;
	axis[1] = axis[1] / prod_scal;
	axis[2] = axis[2] / prod_scal;

	type a= approx_cos(theta/2.0);
	type s = -1.0 * approx_sin(theta / 2.0);
    type b = s * axis[0];
    type c = s * axis[1];
    type d = s * axis[2];
    
	MATRIX result = alloc_matrix(3, 3);

    
    result[0] = a * a + b * b - c * c - d * d;
    result[1] = 2 * (b * c + a * d);
    result[2] = 2 * (b * d - a * c);

    result[3] = 2 * (b * c - a * d);
    result[4] = a * a + c * c - b * b - d * d;
    result[5] = 2 * (c * d + a * b);

    result[6] = 2 * (b * d + a * c);
    result[7] = 2 * (c * d - a * b);
    result[8] = a * a + d * d - b * b - c * c;

    return result;
	}

MATRIX backbone(char* seq, VECTOR phi, VECTOR psi){
	const int n = 256;
	printf("%d", n);
	type r_CaN = 1.46;
	type r_CaC = 1.52;
	type r_CN = 1.33;

	type theta_CaCN = 2.028;
	type theta_CNCa = 2.124;
	type theta_NCaC = 1.940;
	printf("dentro");
	MATRIX coords= alloc_matrix(n*3,3);
	coords[0]=0;
	coords[1]=0;
	coords[2]=0;
	coords[3]=r_CaN;
	coords[4]=0;
	coords[5]=0;
	printf("dentro");
	for(int i=0; i<n; i++){
		int idx = i*9;
		if(i>0){

        	// Posizionamento di N
        	VECTOR v3 = alloc_matrix(1, 3);  //forse si potrebbe allocare un solo vettore fuori dal for e riutilizzarlo
        	v3[0] = coords[idx - 3] - coords[idx - 6]; // C(i-1) - Cα(i-1)
        	v3[1] = coords[idx - 2] - coords[idx - 5];
        	v3[2] = coords[idx - 1] - coords[idx - 4];
        	type norm_v1 = sqrt(v3[0] * v3[0] + v3[1] * v3[1] + v3[2] * v3[2]);
			v3[0]/=norm_v1;
			v3[1]/=norm_v1;
			v3[2]/=norm_v1;
			MATRIX rot= rotation(v3, theta_CNCa);
			VECTOR new_v = alloc_matrix(1, 3);
        	new_v[0] = 0;
       		new_v[1] = r_CN;
        	new_v[2] = 0;
			VECTOR res = alloc_matrix(1, 3);
			vector_matrix_product(new_v, rot, res);
			coords[idx] = coords[idx-3] + res[0];
			coords[idx+1] = coords[idx-2] + res[1];
			coords[idx+2] = coords[idx-1] + res[2];

			//Posizionamento Ca
			VECTOR v2 = alloc_matrix(1, 3);
			v2[0] = coords[idx] - coords[idx - 3]; // N-C
        	v2[1] = coords[idx + 1] - coords[idx - 2];
        	v2[2] = coords[idx + 2] - coords[idx - 1];
			type norm_v2 = sqrt(v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2]);
			v2[0]/=norm_v2;
			v2[1]/=norm_v2;
			v2[2]/=norm_v2;
			rot= rotation(v2, phi[i]);
			new_v[1] = r_CaN;
			vector_matrix_product(new_v, rot, res);
			coords[idx+3] = coords[idx] + res[0];
			coords[idx+4] = coords[idx+1] + res[1];
			coords[idx+5] = coords[idx+2] + res[2];
		}

		VECTOR v3 = alloc_matrix(1, 3);  //forse si potrebbe allocare un solo vettore fuori dal for e riutilizzarlo
        v3[0] = coords[idx + 3] - coords[idx]; // C(i-1) - Cα(i-1)
    	v3[1] = coords[idx + 4] - coords[idx + 1];
    	v3[2] = coords[idx + 5] - coords[idx + 2];
		type norm_v3 = sqrt(v3[0] * v3[0] + v3[1] * v3[1] + v3[2] * v3[2]);
		v3[0]/=norm_v3;
		v3[1]/=norm_v3;
		v3[2]/=norm_v3;
		MATRIX rot= rotation(v3, psi[i]);
		VECTOR new_v = alloc_matrix(1, 3);
        new_v[0] = 0;
		new_v[1] = r_CaC;
		new_v[2] = 0;
		VECTOR res = alloc_matrix(1, 3);
		vector_matrix_product(new_v, rot, res);
		coords[idx + 6] = coords[idx + 3] + res[0];
		coords[idx + 7] = coords[idx + 4] + res[1];
		coords[idx + 8] = coords[idx + 5] + res[2];


	}

	return coords;
}


type rama_energy(VECTOR phi, VECTOR psi) {
    // Costanti di Ramachandran
    
    const type alpha_phi = -57.8;
    const type alpha_psi = -47.0;
    const type beta_phi = -119.0;
    const type beta_psi = 113.0;
	const int n = 256;
    
	
	type energy = 0.0;

    // Itera su tutti gli elementi
    for (int i = 0; i < n; i++) {
        // Calcola la distanza alpha
        type alpha_dist = sqrt(pow(phi[i] - alpha_phi, 2) + pow(psi[i] - alpha_psi, 2));

        // Calcola la distanza beta
        type beta_dist = sqrt(pow(phi[i] - beta_phi, 2) + pow(psi[i] - beta_psi, 2));

        // Somma il contributo minimo all'energia con confronto esplicito
        if (alpha_dist < beta_dist) {
            energy += 0.5 * alpha_dist;
        } else {
            energy += 0.5 * beta_dist;
        }
    }

    return energy;
}

extern MATRIX coordsca(MATRIX coords) {
    const int n = 256;
	MATRIX Cacoords = alloc_matrix(n, 3);
    
	for (int i = 0; i < n; i++) {
        Cacoords[i * 3] = coords[i * 9 + 3]; //X
        Cacoords[i* 3 +1] = coords[i * 9 + 4]; //Y
        Cacoords[i * 3 + 2] = coords[i * 9 + 5]; //Z
    }
    return Cacoords; 
}


type distanza (MATRIX coordinateCa, int i, int j){
		int x_df = coordinateCa[3*i] - coordinateCa[3*j];
		int y_df = coordinateCa[3*i+1] - coordinateCa[3*j+1];
		int z_df = coordinateCa[3*i+2] - coordinateCa[3*j+2];
		return sqrt(pow(x_df,2) + pow(y_df,2 ) +pow(z_df,2));
}

type hydrofobic_energy (char* sequenza, MATRIX coordinate){
	type energy = 0;
	MATRIX coordinateCa = coordsca(coordinate);
	printf("coordinateCA");
	const int n = 256;

	for(int i=0; i< n; i++){
		for(int j= i+1; j<n; j++){
			type dist = distanza(coordinate, i, j);
			if(dist < 10.0){
				energy += (hydrophobicity[(int)sequenza[i]] * hydrophobicity[(int)sequenza[j]] )/ dist;
			}
		}
	}
	return energy;
}

extern type electrostatic_energy(char* s, MATRIX coords){
	MATRIX coordinateCa= coordsca(coords);
	printf("elecCoor");
	type energy= 0.0;
	const int n = 256;
	for(int i=0; i < n; i++){
		for(int j= i+1; j < n; j++){
			if(i!= j){
				type dist= distanza(coordinateCa, i, j);
				if(dist < 10.0 && charge[(int)s[i]] !=0 && charge[(int)s[j] != 0] ){
					energy += (charge[(int)s[i]]*charge[(int)s[j]])/(dist*4.0);
				}
			}
		}
	}
	return energy; 
}

extern type packing_energy(char*s,MATRIX coords) {
    const int n = 256; 
    MATRIX cacoords = coordsca(coords);
    type energy = 0.0;
    for (int i = 0; i < n; i++) {
		type  density = 0.0;
		for (int j = 0; j < n; j++) {
			if(i != j){
				type dist = distanza(cacoords, i, j);
				if (dist < 10.0) {
					density = density + volume[(int)s[j]] / (pow(dist, 3)); 
				}
			}
			energy = energy + pow((volume[(int)s[i]] - density), 2);
		}
    }
	return energy;
}



extern type energy(char* seq, VECTOR phi, VECTOR psi){
	printf("energy");
	MATRIX coords= backbone(seq, phi, psi);
	
	printf("backbone");
	type rama= rama_energy(phi, psi);
	printf("rama");
	type hydro = hydrofobic_energy(seq, coords);
	printf("hydro");
	type elec = electrostatic_energy(seq, coords);
	printf("elec");
	type pack = packing_energy(seq, coords);
	printf("pack");
	type w_rama= 1.0;
	type w_hydro= 0.5;
	type w_elec= 0.2;
	type w_pack= 0.3;

	type tot= (w_rama*rama) + (w_elec*elec)+(w_hydro*hydro)+(w_pack*pack);

	return tot;
}



void pst(params* input){
	
	int n = input->N;
	
	type T = input->to;
	printf("%d %f",n, T );
	VECTOR phi= input->phi;
	
	VECTOR psi= input->psi;
	
	type E= energy(input->seq, phi, psi);
	
	type t=0.0;
	
	while(T>0){
		srand((int)time(NULL));

    	// Genera un numero casuale tra 0 e n
    	int i = rand() % (n + 1);
		type theta_phi = 2.2;
		phi[i] = phi[i] + theta_phi;
		type theta_psi = 1.4;
		psi[i] = psi[i] + theta_psi;
		type deltaE= energy(input->seq, phi, psi) - E;

		if(deltaE<=0){
			E= energy(input->seq, phi, psi);
		}
		else{
			type P = pow(M_E, (-deltaE/(input->k*T)));
			type r = random();

			if (r<=P)
				E= energy(input->seq, phi, psi);
			else{
				phi[i] = phi[i] - theta_phi;
				psi[i] = psi[i] + theta_psi;
			}
		}
		t=t+1;
		T= input->to - sqrt(input->alpha*t);
	}
}

int main(int argc, char** argv) {
	char fname_phi[256];
	char fname_psi[256];
	char* seqfilename = NULL;
	clock_t t;
	float time;
	int d;
	
	//
	// Imposta i valori di default dei parametri
	//
	params* input = malloc(sizeof(params));
	input->seq = NULL;	
	input->N = -1;			
	input->to = -1;
	input->alpha = -1;
	input->k = -1;		
	input->sd = -1;		
	input->phi = NULL;		
	input->psi = NULL;
	input->silent = 0;
	input->display = 0;
	input->e = -1;
	input->hydrophobicity = hydrophobicity;
	input->volume = volume;
	input->charge = charge;


	//
	// Visualizza la sintassi del passaggio dei parametri da riga comandi
	//
	if(argc <= 1){
		printf("%s -seq <SEQ> -to <to> -alpha <alpha> -k <k> -sd <sd> [-s] [-d]\n", argv[0]);
		printf("\nParameters:\n");
		printf("\tSEQ: il nome del file ds2 contenente la sequenza amminoacidica\n");
		printf("\tto: parametro di temperatura\n");
		printf("\talpha: tasso di raffredamento\n");
		printf("\tk: costante\n");
		printf("\tsd: seed per la generazione casuale\n");
		printf("\nOptions:\n");
		printf("\t-s: modo silenzioso, nessuna stampa, default 0 - false\n");
		printf("\t-d: stampa a video i risultati, default 0 - false\n");
		exit(0);
	}

	//
	// Legge i valori dei parametri da riga comandi
	//

	int par = 1;
	while (par < argc) {
		if (strcmp(argv[par],"-s") == 0) {
			input->silent = 1;
			par++;
		} else if (strcmp(argv[par],"-d") == 0) {
			input->display = 1;
			par++;
		} else if (strcmp(argv[par],"-seq") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing dataset file name!\n");
				exit(1);
			}
			seqfilename = argv[par];
			par++;
		} else if (strcmp(argv[par],"-to") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing to value!\n");
				exit(1);
			}
			input->to = atof(argv[par]);
			par++;
		} else if (strcmp(argv[par],"-alpha") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing alpha value!\n");
				exit(1);
			}
			input->alpha = atof(argv[par]);
			par++;
		} else if (strcmp(argv[par],"-k") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing k value!\n");
				exit(1);
			}
			input->k = atof(argv[par]);
			par++;
		} else if (strcmp(argv[par],"-sd") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing seed value!\n");
				exit(1);
			}
			input->sd = atoi(argv[par]);
			par++;
		}else{
			printf("WARNING: unrecognized parameter '%s'!\n",argv[par]);
			par++;
		}
	}

	//
	// Legge i dati e verifica la correttezza dei parametri
	//
	if(seqfilename == NULL || strlen(seqfilename) == 0){
		printf("Missing ds file name!\n");
		exit(1);
	}

	input->seq = load_seq(seqfilename, &input->N, &d);

	
	if(d != 1){
		printf("Invalid size of sequence file, should be %ix1!\n", input->N);
		exit(1);
	} 

	if(input->to <= 0){
		printf("Invalid value of to parameter!\n");
		exit(1);
	}

	if(input->k <= 0){
		printf("Invalid value of k parameter!\n");
		exit(1);
	}

	if(input->alpha <= 0){
		printf("Invalid value of alpha parameter!\n");
		exit(1);
	}

	input->phi = alloc_matrix(input->N, 1);
	input->psi = alloc_matrix(input->N, 1);
	// Impostazione seed 
	srand(input->sd);
	// Inizializzazione dei valori
	gen_rnd_mat(input->phi, input->N);
	gen_rnd_mat(input->psi, input->N);

	//
	// Visualizza il valore dei parametri
	//

	if(!input->silent){
		printf("Dataset file name: '%s'\n", seqfilename);
		printf("Sequence lenght: %d\n", input->N);
	}

	// COMMENTARE QUESTA RIGA!
	//prova(input);
	//

	//
	// Predizione struttura terziaria
	//
	
	t = clock();
	printf("Ciao Angiulli");
	pst(input);
	t = clock() - t;
	time = ((float)t)/CLOCKS_PER_SEC;
	printf("Ciao Fassetti");



	if(!input->silent)
		printf("PST time = %.3f secs\n", time);
	else
		printf("%.3f\n", time);

	//
	// Salva il risultato
	//
	sprintf(fname_phi, "out32_%d_%d_phi.ds2", input->N, input->sd);
	save_out(fname_phi, input->phi, input->N);
	sprintf(fname_psi, "out32_%d_%d_psi.ds2", input->N, input->sd);
	save_out(fname_psi, input->psi, input->N);
	if(input->display){
		if(input->phi == NULL || input->psi == NULL)
			printf("out: NULL\n");
		else{
			int i,j;
			printf("energy: %f, phi: [", input->e);
			for(i=0; i<input->N; i++){
				printf("%f,", input->phi[i]);
			}
			printf("]\n");
			printf("psi: [");
			for(i=0; i<input->N; i++){
				printf("%f,", input->psi[i]);
			}
			printf("]\n");
		}
	}

	if(!input->silent)
		printf("\nDone.\n");

	dealloc_matrix(input->phi);
	dealloc_matrix(input->psi);
	free(input);

	return 0;
}
