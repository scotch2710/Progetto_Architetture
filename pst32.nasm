; ---------------------------------------------------------
; Predizione struttura terziaria con istruzioni SSE a 32 bit
; ---------------------------------------------------------
; F. Angiulli
; 23/11/2017
;

;
; Software necessario per l'esecuzione:
;
;     NASM (www.nasm.us)
;     GCC (gcc.gnu.org)
;
; entrambi sono disponibili come pacchetti software 
; installabili mediante il packaging tool del sistema 
; operativo; per esempio, su Ubuntu, mediante i comandi:
;
;     sudo apt-get install nasm
;     sudo apt-get install gcc
;
; potrebbe essere necessario installare le seguenti librerie:
;
;     sudo apt-get install lib32gcc-4.8-dev (o altra versione)
;     sudo apt-get install libc6-dev-i386
;
; Per generare file oggetto:
;
;     nasm -f elf32 fss32.nasm 
;
%include "sseutils32.nasm"

section .data			; Sezione contenente dati inizializzati
	unroll equ 4
	n 	   equ 256
	dim    equ 4

section .bss			; Sezione contenente dati non inizializzati
	alignb 16
	e		resd		1

section .text			; Sezione contenente il codice macchina


; ----------------------------------------------------------
; macro per l'allocazione dinamica della memoria
;
;	getmem	<size>,<elements>
;
; alloca un'area di memoria di <size>*<elements> bytes
; (allineata a 16 bytes) e restituisce in EAX
; l'indirizzo del primo bytes del blocco allocato
; (funziona mediante chiamata a funzione C, per cui
; altri registri potrebbero essere modificati)
;
;	fremem	<address>
;
; dealloca l'area di memoria che ha inizio dall'indirizzo
; <address> precedentemente allocata con getmem
; (funziona mediante chiamata a funzione C, per cui
; altri registri potrebbero essere modificati)

extern get_block
extern free_block

%macro	getmem	2
	mov	eax, %1
	push	eax
	mov	eax, %2
	push	eax
	call	get_block
	add	esp, 8
%endmacro

%macro	fremem	1
	push	%1
	call	free_block
	add	esp, 4
%endmacro

; ------------------------------------------------------------
; Funzioni
; ------------------------------------------------------------

global distanza1
global coordsca

N 		equ 	8
seq		equ		12
sd		equ		16
t0		equ		20
alpha	equ		24
k		equ		28
hydro	equ		32
volume	equ		36
charge 	equ		40
energy 	equ		44

;char* seq;		// sequenza di amminoacidi
;	int N;			// lunghezza sequenza
;	unsigned int sd; 	// seed per la generazione casuale
;	type to;		// temperatura INIZIALE
;	type alpha;		// tasso di raffredamento
;	type k;		// costante
;	VECTOR hydrophobicity; // hydrophobicity
;	VECTOR volume;		// volume
;	VECTOR charge;		// charge
;	VECTOR phi;		// vettore angoli phi
;	VECTOR psi;		// vettore angoli psi
;	type e;		// energy
;	int display;
;	int silent;


; ------------------------------------------------------------
; Funzioni
; ------------------------------------------------------------

; ------------------------------------------------------------
; Funzione distanza1
; ------------------------------------------------------------
distanza1:    
	push	ebp			; salva il Base Pointer
	mov		ebp, esp	; il Base Pointer punta al Record di Attivazione corrente
	push	ebx			; salva i registri da preservare
	push	esi
	push	edi                 
	
	; Carica i parametri    
	mov     esi, [ebp+8]               ; coordinateCa    
	mov     eax, [ebp+12]              ; i    
	mov     ebx, [ebp+16]              ; j    
	
	; Calcola 3*i e 3*j    
	lea     eax, [eax + eax * 2]       ; eax = 3*i    
	lea     ebx, [ebx + ebx * 2]       ; ebx = 3*j    
	
	; X_df    
	movss   xmm0, dword [esi + eax * 4] ; xmm0 = coordinateCa[3*i]    
	movss   xmm1, dword [esi + ebx * 4] ; xmm1 = coordinateCa[3*j]    
	subss   xmm0, xmm1                  ; xmm0 = x_df    
	
	; Y_df    
	add     eax, 1    
	add     ebx, 1    
	movss   edx, dword [esi + eax * 4]    
	movss   xmm3, dword [esi + ebx * 4]    
	subss   edx, xmm3                  ; edx = y_df    
	
	; Z_df    
	add     eax, 1    
	add     ebx, 1    
	movss   xmm4, dword [esi + eax * 4]    
	movss   xmm5, dword [esi + ebx * 4]    
	subss   xmm4, xmm5                  ; xmm4 = z_df    
	
	; Calcola x_df^2, y_df^2, z_df^2 e somma    
	mulss   xmm0, xmm0                  ; xmm0 = x_df^2    
	mulss   edx, edx                  ; edx = y_df^2    
	mulss   xmm4, xmm4                  ; xmm4 = z_df^2    
	addss   xmm0, edx                  ; xmm0 += y_df^2    
	addss   xmm0, xmm4                  ; xmm0 += z_df^2    
	; Radice quadrata    
	sqrtss  xmm0, xmm0                  ; xmm0 = sqrt(x_df^2 + y_df^2 + z_df^2)    
	mov 	eax, [ebp+20]				;[ebp+20] contiene l'indirizzo della variabile dist passata come parametro (&dist)
	movss 	[eax], xmm0					;[eax] inserisce nell'indirizzo passato in eax il valore risultante in xmm0

	
	pop edi
	pop	esi
	pop	ebx
	mov	esp, ebp	; ripristina lo Stack Pointer 
	pop ebp    
	ret

; ------------------------------------------------------------
; Funzione coordsca
; ------------------------------------------------------------
; void coordsca(MATRIX coords, MATRIX cacoords) {
;     const int n = 256;
; 	for (int i = 0; i < n; i++) {
;         cacoords[i * 3] = coords[i * 9 + 3]; //X
;         cacoords[i* 3 + 1] = coords[i * 9 + 4]; //Y
;         cacoords[i * 3 + 2] = coords[i * 9 + 5]; //Z
;     }
; }
coordsca:
	push	ebp			; salva il Base Pointer
	mov		ebp, esp	; il Base Pointer punta al Record di Attivazione corrente
	push	ebx			; salva i registri da preservare
	push	esi
	push	edi  

	mov ebx, [ebp+8]    ;coords
	mov eax, [ebp+12]		;cacoords

	xor esi, esi 		;ESI: i=0
	
    forCacoords:
		cmp esi, 256
		je fineCacoords
		mov ecx, esi ;ecx contatore di coords
		imul ecx, 9
		

		mov edx, esi ; edx contatore di cacords
		shl edx, 1
		add edx, esi 
		
		movss xmm0, [ebx + ecx*dim + 3*dim] ; mette in xmm0 coords[i*9+3] salvando x
		movss [eax + edx*dim], xmm0

        inc edx
		movss xmm0, [ebx + ecx*dim + 4*dim] ; mette in xmm0 coords[i*9+4] salvando y
		movss [eax + edx*dim], xmm0

		inc edx
		movss xmm0, [ebx + ecx*dim + 5*dim] ; mette in xmm0 coords[i*9+5] salvando z
		movss [eax + edx*dim], xmm0
		inc esi
		jmp forCacoords
		;--------fine ciclo for--------

	fineCacoords:
	;--- fine logica ---
	pop edi
	pop	esi
	pop	ebx
	mov	esp, ebp	; ripristina lo Stack Pointer 
	pop ebp    
	ret