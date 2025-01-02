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
	; Costanti di Ramachandran
	alpha_phi	 dd -57.8
	alpha_psi	 dd -47.0
	beta_phi	 dd -119.0
	beta_psi	 dd 113.0
	un_mezzo     dd 0.5
	dieci        dd 10.0
	uno 	     dd 1.0
	24s			 dd 24.0
	720s         dd 720.0	
	2s           dd 2.0
	6s 			 dd 6.0
	120s 	     dd 120.0
	5040s		 dd 5040.0
	meno1		 dd -1.0
	; Hydrophobicity
	alignb 16
	hydrophobicity1 dd 1.8, -1, 2.5, -3.5, -3.5, 2.8, -0.4, -3.2, 4.5, -1, -3.9, 3.8, 1.9, -3.5, -1, -1.6, -3.5, -4.5, -0.8, -0.7, -1, 4.2, -0.9, -1, -1.3, -1

	alignb 16 
	; Volume
	volume1 dd 88.6, -1, 108.5, 111.1, 138.4, 189.9, 60.1, 153.2, 166.7, -1, 168.6, 166.7, 162.9, 114.1, -1, 112.7, 143.8, 173.4, 89.0, 116.1, -1, 140.0, 227.8, -1, 193.6, -1
	
	alignb 16
	; Charge
	charge1 dd 0, -1, 0, -1, -1, 0, 0, 0.5, 0, -1, 1, 0, 0, 0, -1, 0, 0, 1, 0, 0, -1, 0, 0, -1, 0, -1


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
global rama_energy
global rotation
;global hydrofobic_energy

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
	movss   xmm2, dword [esi + eax * 4]    
	movss   xmm3, dword [esi + ebx * 4]    
	subss   xmm2, xmm3                  ; edx = y_df    
	
	; Z_df    
	add     eax, 1    
	add     ebx, 1    
	movss   xmm4, dword [esi + eax * 4]    
	movss   xmm5, dword [esi + ebx * 4]    
	subss   xmm4, xmm5                  ; xmm4 = z_df    
	
	; Calcola x_df^2, y_df^2, z_df^2 e somma    
	mulss   xmm0, xmm0                  ; xmm0 = x_df^2    
	mulss   xmm2, xmm2                  ; edx = y_df^2    
	mulss   xmm4, xmm4                  ; xmm4 = z_df^2    
	addss   xmm0, xmm2                  ; xmm0 += y_df^2    
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
		jge fineCacoords
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

; ------------------------------------------------------------
; Funzione rama_energy
; ------------------------------------------------------------
	;-----------implematation of the rama_energy method------------------------
	; 	extern type rama_energy(VECTOR phi, VECTOR psi) {
	; 	// Costanti di Ramachandran
		
	; 	const type alpha_phi = -57.8;
	; 	const type alpha_psi = -47.0;
	; 	const type beta_phi = -119.0;
	; 	const type beta_psi = 113.0;
	; 	const int n = 256;
		
		
	; 	type energy = 0.0;

	; 	// Itera su tutti gli elementi
	; 	for (int i = 0; i < n; i++) {
	; 		// Calcola la distanza alpha
	; 		type alpha_dist = sqrt((phi[i] - alpha_phi) * (phi[i] - alpha_phi) + (psi[i] - alpha_psi) * (psi[i] - alpha_psi));

	; 		// Calcola la distanza beta
	; 		type beta_dist = sqrt((phi[i] - beta_phi) * (phi[i] - beta_phi) + (psi[i] - beta_psi) * (psi[i] - beta_psi));

	; 		// Somma il contributo minimo all'energia con confronto esplicito
	; 		if (alpha_dist < beta_dist) {
	; 			energy += 0.5 * alpha_dist;
	; 		} else {
	; 			energy += 0.5 * beta_dist;
	; 		}
	; 	}

	; 	return energy;
	; }
	;----------------------------------------------------------------------------
		; IMPLEMENTAZIONE DI ROTATION
	;extern MATRIX rotation(VECTOR axis, type theta){
; 	//  {
; // 	type prod_scal= (axis[0]*axis[0])+(axis[1]*axis[1])+(axis[2]*axis[2]);
	
; // 	axis[0] = axis[0] / prod_scal;
; // 	axis[1] = axis[1] / prod_scal;
; // 	axis[2] = axis[2] / prod_scal;

; // 	type a = approx_cos(theta/2.0);
; // 	type s = -1.0 * approx_sin(theta / 2.0);
; //     type b = s * axis[0];
; //     type c = s * axis[1];
; //     type d = s * axis[2];

    
; //     result[0] = a * a + b * b - c * c - d * d; 
; //     result[1] = 2 * (b * c + a * d);
; //     result[2] = 2 * (b * d - a * c);

; //     result[3] = 2 * (b * c - a * d);
; //     result[4] = a * a + c * c - b * b - d * d;
; //     result[5] = 2 * (c * d + a * b);

; //     result[6] = 2 * (b * d + a * c);
; //     result[7] = 2 * (c * d - a * b);
; //     result[8] = a * a + d * d - b * b - c * c;
; //     return;
; // 	} a*a, b*b, c*c, d*d, b*c, a*d, b*d, a*c, c*d, a*b
; 0, 4, 8 --- 1, 3 -- 2, 6 -- 5, 7

rotation:
	push	ebp			; salva il Base Pointer
	mov		ebp, esp	; il Base Pointer punta al Record di Attivazione corrente
	push	ebx			; salva i registri da preservare
	push	esi
	push	edi

	mov ebx, [ebp+8]    ;axis
	mov ecx, [ebp+12]		;theta
	mov eax, [ebp+16]		;result

	; inizio logica

	; Calcola il prodotto scalare
	movss xmm0, [ebx] ; axis[0]
	mulss xmm0, xmm0

	movss xmm1, [ebx+dim] ; axis[1]
	mulss xmm1, xmm1
	addss xmm0, xmm1
	;pxor xmm1, xmm1

	movss xmm1, [ebx+2*dim] ; axis[2]
	mulss xmm1, xmm1
	addss xmm0, xmm1
	;pxor xmm1, xmm1
	; xmm0 ha il prodotto scalare

	; Calcola 1/prodotto scalare
	movss xmm1, [ebx] ; axis[0]
	divss xmm1, xmm0
	movss [ebx], xmm1

	;pxor xmm1, xmm1
	movss xmm1, [ebx+dim] ; axis[1]
	divss xmm1, xmm0
	movss [ebx+dim], xmm1

	;pxor xmm1, xmm1
	movss xmm1, [ebx+2*dim] ; axis[2]
	divss xmm1, xmm0
	movss [ebx+2*dim], xmm1

	; Calcola a con approx del coseno
	; return 1 - (x2 / 2.0) + (x2 * x2 / 24.0) - (x2 * x2 * x2 / 720.0);

	movss xmm1, [ecx]
	divss xmm1, [2s] ; in rotation theta/2.0
	movss xmm7, xmm1 ; xmm7 =theta
	mulss xmm1, xmm1 ; xmm1 = theta^2
	movss xmm2, xmm1
	divss xmm2, [2s] ; xmm2 = theta^2 / 2.0

	mulss xmm3, xmm1, xmm1
	movss xmm4, xmm3
	movss xmm6, xmm3 ; xmm6 = theta^2 * theta^2
	divss xmm3, [24s] ; xmm3 = theta^2 * theta^2 / 24.0

	mulss xmm4, xmm1
	divss xmm4, [720s] ; xmm4 = theta^2 * theta^2 * theta^2 / 720.0

	movss xmm5, [uno]
	subss xmm5, xmm2
	addss xmm5, xmm3
	subss xmm5, xmm4 ; xmm5 = 1 - (theta^2 / 2.0) + (theta^2 * theta^2 / 24.0) - (theta^2 * theta^2 * theta^2 / 720.0)

	; xmm5 ha l'approssimazione del coseno

	; Calcola s con approx del seno
	;return theta - (x2 * theta / 6.0) + (x2 * x2 * theta / 120.0) - (x2 * x2 * x2 * theta/ 5040.0);

	mulss xmm2, xmm7, xmm1 ; xmm2 = theta^3
	; free register , xmm4
	movss xmm3, xmm2
	divss xmm3, [6s] ; xmm3 = theta^3 / 6.0
	
	mulss  xmm6, xmm7 ; xmm6 = theta^5
	movss xmm4, xmm6
	divss xmm6, [120s] ; xmm6 = theta^5 / 120.0

	mulss xmm4, xmm1 ; xmm4 = theta^7
	divss xmm4, [5040s] ; xmm4 = theta^7 / 5040.0

	subss xmm7, xmm3
	addss xmm7, xmm6
	subbss xmm7, xmm4 ; xmm7 = theta - (theta^3 / 6.0) + (theta^5 / 120.0) - (theta^7 / 5040.0)

	; xmm7 ha l'approssimazione del seno

	imull xmm7, [meno1] ; xmm7 = -1.0 * approx_sin(theta / 2.0)

	; a = xmm5
	; s = xmm7

	; Calcol0 b, c, d

	movss xmm0, [ebx] ; axis[0]
	mulss xmm0, xmm7 ; b = s * axis[0]

	movss xmm1, [ebx+dim] ; axis[1]
	mulss xmm1, xmm7 ; c = s * axis[1]

	movss xmm2, [ebx+2*dim] ; axis[2]
	mulss xmm2, xmm7 ; d = s * axis[2]


	; calcolo delle righe del risultato; 1, 3
	; //     result[3] = 2 * (b * c - a * d);
	; //     result[1] = 2 * (b * c + a * d);

	movss xmm3, xmm0 ; xmm3 = b
	mulss xmm3, xmm1 ; xmm3 = b * c
	movss xmm4, xmm5; xmm4 = a
	mulss xmm4, xmm2 ; xmm4 = a * d
	movss xmm6, xmm3 ; xmm6 = b * c
	movss xmm7, xmm3 ; xmm7 = b * c
	subss xmm6, xmm4 ; xmm3 = b * c - a * d
	adds xmm7,  xmm4 ; xmm7 = b * c + a * d
	mulss xmm6, [2s] ; xmm6 = 2 * (b * c - a * d)
	mulss xmm7, [2s] ; xmm7 = 2 * (b * c + a * d)
	movss [eax+dim*3], xmm6 ; result[3] = 2 * (b * c - a * d)
	movss [eax+dim], xmm7   ; result[1] = 2 * (b * c + a * d

	; calcolo delle righe del risultato; 2, 6
	; //     result[6] = 2 * (b * d + a * c);
	; //     result[2] = 2 * (b * d - a * c);
	; a= xmm5, b= xmm0, c= xmm1, d= xmm2  ---> registri liberi xmm3, xmm4, xmm6, xmm7

	movss xmm3, xmm0 ; xmm3 = b
	mulss xmm3, xmm2 ; xmm3 = b * d
	movss xmm4, xmm5; xmm4 = a
	mulss xmm4, xmm1 ; xmm4 = a * c
	movss xmm6, xmm3 ; xmm6 = b * d
	movss xmm7, xmm3 ; xmm7 = b * d
	subss xmm6, xmm4 ; xmm3 = b * d - a * c
	addss xmm7,  xmm4 ; xmm7 = b * d + a * c
	mulss xmm6, [2s] ; xmm6 = 2 * (b * d - a * c)
	mulss xmm7, [2s] ; xmm7 = 2 * (b * d + a * c)
	movss [eax+dim*6], xmm6   ; result[6] = 2 * (b * d - a * c)
	movss [eax+dim*2], xmm7   ; result[2] = 2 * (b * d + a * c)

	; calcolo delle righe del risultato; 5, 7
	; //     result[7] = 2 * (c * d - a * b);
	; //     result[5] = 2 * (c * d + a * b);
	; a= xmm5, b= xmm0, c= xmm1, d= xmm2  ---> registri liberi xmm3, xmm4, xmm6, xmm7

	movss xmm3, xmm1 ; xmm3 = c
	mulss xmm3, xmm2 ; xmm3 = c * d
	movss xmm4, xmm5; xmm4 = a
	mulss xmm4, xmm0 ; xmm4 = a * b
	movss xmm6, xmm3 ; xmm6 = c * d
	movss xmm7, xmm3 ; xmm7 = c * d
	subss xmm6, xmm4 ; xmm3 = c * d - a * b
	addss xmm7,  xmm4 ; xmm7 = c * d + a * b
	mulss xmm6, [2s] ; xmm6 = 2 * (c * d - a * b)
	mulss xmm7, [2s] ; xmm7 = 2 * (c * d + a * b)
	movss [eax+dim*7], xmm6   ; result[7] = 2 * (c * d - a * b)
	movss [eax+dim*5], xmm7   ; result[5] = 2 * (c * d + a * b)


	; implementazioni delle seguenti righe del result: 0, 4, 8
	; a= xmm5, b= xmm0, c= xmm1, d= xmm2  ---> registri liberi xmm3, xmm4, xmm6, xmm7
	; calcolare a*a, b*b, c*c, d*d
	mulss xmm5, xmm5 ; xmm5 = a*a
	mulss xmm0, xmm0 ; xmm0 = b*b
	mulss xmm1, xmm1 ; xmm1 = c*c
	mulss xmm2, xmm2 ; xmm2 = d*d

	; //     result[0] = a * a + b * b - c * c - d * d; 
	movss xmm3, xmm5 ; xmm3 = a*a
	addss xmm3, xmm0 ; xmm3 = a*a + b*b
	subss xmm3, xmm1 ; xmm3 = a*a + b*b - c*c
	subss xmm3, xmm2 ; xmm3 = a*a + b*b - c*c - d*d
	movss [eax], xmm3 ; result[0] = a*a + b*b - c*c - d*d

	; //     result[4] = a * a + c * c - b * b - d * d;
	movss xmm4, xmm5 ; xmm4 = a*a
	adds xmm4, xmm1 ; xmm4 = a*a + c*c
	subss xmm4, xmm0 ; xmm4 = a*a + c*c - b*b
	subss xmm4, xmm2 ; xmm4 = a*a + c*c - b*b - d*d
	movss [eax+dim*4], xmm4 ; result[4] = a*a + c*c - b*b - d*d

	; //     result[8] = a * a + d * d - b * b - c * c;
	movss xmm6, xmm5 ; xmm6 = a*a
	adds xmm6, xmm2 ; xmm6 = a*a + d*d
	subss xmm6, xmm0 ; xmm6 = a*a + d*d - b*b
	subss xmm6, xmm1 ; xmm6 = a*a + d*d - b*b - c*c
	movss [eax+dim*8], xmm6 ; result[8] = a*a + d*d - b*b - c*c

	;--- fine logica ---

	pop edi
	pop	esi
	pop	ebx
	mov	esp, ebp	; ripristina lo Stack Pointer 
	pop ebp    
	ret

	




















	;--------------------------------------------------------------------
rama_energy: 
	push	ebp			; salva il Base Pointer
	mov		ebp, esp	; il Base Pointer punta al Record di Attivazione corrente
	push	ebx			; salva i registri da preservare
	push	esi
	push	edi  

    mov ebx, [ebp+8]    	;phi
	mov ecx, [ebp+12]		;psi

	
	xor esi, esi 			;esi: i=0
	xorps xmm0, xmm0        ;init energy = 0.0
	
    forRamaEnergy:
		cmp esi, 256
		jge fineRamaEnergy
		; uso xmm4 per salvare phi[i]
		; uso xmm5 per salvare psi[i]

		;--------calcolo alpha_dist--------
		movss xmm4, [ebx+esi*dim]   ; recupero phi[i]
		movss xmm6, xmm4 			; salvo phi[i] in xmm6
		subss xmm6, [alpha_phi] 			; phi[i] - alpha_phi
		mulss xmm6, xmm6			;(phi[i] - alpha_phi)^2

		movss xmm5, [ecx+esi*dim]   ; recupero psi[i]
		movss xmm7, xmm5 			; salvo psi[i] in xmm7
		subss xmm7, [alpha_psi]		; psi[i] - alpha_psi
		mulss xmm7, xmm7  	        ; (psi[i] - alpha_psi)^2

		addss xmm6, xmm7 			; (phi[i] - alpha_phi)^2 + (psi[i] - alpha_psi)^2
		
		; in xmm6 ho alpha_dist
   
		;--------calcolo beta_dist--------
		subss xmm4, [beta_phi]		; phi[i] - beta_phi
		mulss xmm4, xmm4			; (phi[i] - beta_phi)^2

		subss xmm5, [beta_psi]		; psi[i] - beta_psi
		mulss xmm5, xmm5			; (psi[i] - beta_psi)^2
		addss xmm4, xmm5 			; (phi[i] - beta_phi)^2 + (psi[i] - beta_psi)^2
		
		; in xmm4 ho beta_dist

		;--------confronto esplicito--------
		comiss xmm6, xmm4 			; xmm6 < xmm4
		jl alpha_dist_minore
		; se distanza alpha maggiore
		sqrtss xmm6, xmm6 			; sqrt((phi[i] - alpha_phi)^2 + (psi[i] - alpha_psi)^2)
		movss xmm5, [un_mezzo]
		mulss xmm6, xmm5	    	; 0.5 * alpha_dist
		addss xmm0, xmm6 			; energy += 0.5 * alpha_dist
		jmp fine_if
		

		;--------alpha_dist < beta_dist--------
		alpha_dist_minore:
		sqrtss xmm4, xmm4 			; sqrt((phi[i] - beta_phi)^2 + (psi[i] - beta_psi)^2)
		mulss xmm4, xmm5			; 0.5 * beta_dist
		addss xmm0, xmm4 			; energy += 0.5 * beta_dist
	
	fine_if:
		inc esi
		jmp forRamaEnergy

    fineRamaEnergy:
		mov eax, [ebp+16]
		movss [eax], xmm0
	;--- fine logica ---

	pop edi
	pop	esi
	pop	ebx
	mov	esp, ebp	; ripristina lo Stack Pointer 
	pop ebp    
	ret

; extern void hydrofobic_energy (char* sequenza, MATRIX coordinate, MATRIX cacoords, type* energy){
; 	//type energy = 0.0;
; 	const int n = 256;
	
; 	for(int i=0; i< n; i++){
; 		for(int j= i+1; j<n; j++){
; 			//type dist = distanza(cacoords, i, j);
; 			type dist = 0.0;
; 			distanza1(cacoords, i, j, &dist);
; 			//printf("distanza: %f\n", dist);
; 			if(dist < 10.0){
; 				*energy += (hydrophobicity[(int)sequenza[i]-65] * hydrophobicity[(int)sequenza[j]-65] )/ dist;
; 			}
; 		}
; 	}
; 	//printf("energy hhydro: %f\n", energy);
; 	return;
; ; }

hydrofobic_energy:
	push	ebp			; salva il Base Pointer
	mov		ebp, esp	; il Base Pointer punta al Record di Attivazione corrente
	push	ebx			; salva i registri da preservare
	push	esi
	push	edi  

	mov ebx, [ebp+8]    	;sequenza
	mov ecx, [ebp+12]		;coordinate
	mov edx, [ebp+16]		;cacoords

	xor esi, esi 			;esi: i=0
	xorps xmm2, xmm2        ;init energy = 0.0

	externalLoop:
		cmp esi, 256
		jge fineHydrofobicEnergy

		xor edi, edi 			;edi: j=0
		mov edi, esi			; edi = j = i
		inc edi					; edi = j= i+1
						
		internaloop:
			cmp edi, 256
			jge fine_internal_loop

			;--------calcolo distanza--------

			
			;uso xmm4 per salvare dist
			;xor eax, eax 
			push eax
			
			push edi
			push esi
			push edx
			
			
			call distanza1
			add esp, 16
			xorps xmm4, xmm4
			movss xmm4, [eax] ; xmm4 = dist
			;movss xmm4, [ebp+20] ; xmm4 = dist
			
			comiss xmm4,  [dieci]
			jge distanza_maggiore

			;--------calcolo energia--------
			;push eax

			movzx eax, byte [ebx+esi] ; sequenza[i]
			sub eax, 65				  ; sequenza[i] - 65
			movss xmm0, [hydrophobicity1 + eax*4] ; hydrophobicity[sequenza[i]-65]

			movzx eax, byte [ebx+edi] ; sequenza[j]
			sub eax, 65 ; sequenza[j] - 65
			movss xmm1, [hydrophobicity1 + eax*4] ; hydrophobicity[sequenza[j]-65]
			;pop eax
			mulss xmm0, xmm1 ; hydrophobicity[sequenza[i]-65] * hydrophobicity[sequenza[j]-65]
			divss xmm0, xmm4 
			addss xmm2, xmm0 
		distanza_maggiore:
			inc edi
			jmp internaloop
		
		fine_internal_loop:
			inc esi
			jmp externalLoop
	fineHydrofobicEnergy:
	;xor eax, eax 
	mov eax, [ebp+20]
	movss [eax], xmm2
	
	pop edi
	pop	esi
	pop	ebx
	mov	esp, ebp	; ripristina lo Stack Pointer 
	pop ebp    
	ret








