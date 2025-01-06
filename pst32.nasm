; ---------------------------------------------------------
; Predizione struttura terziaria con istruzioni SSE a 32 bit
; ---------------------------------------------------------
; F. Angiulli
; 23/11/2017
;

;
; Software[size]ecessario per l'esecuzione:
;
;    [size]ASM (www.nasm.us)
;     GCC (gcc.gnu.org)
;
; entrambi sono disponibili come pacchetti software 
; installabili mediante il packaging tool del sistema 
; operativo; per esempio, su Ubuntu, mediante i comandi:
;
;     sudo apt-get install[size]asm
;     sudo apt-get install gcc
;
; potrebbe essere[size]ecessario installare le seguenti librerie:
;
;     sudo apt-get install lib32gcc-4.8-dev (o altra versione)
;     sudo apt-get install libc6-dev-i386
;
; Per generare file oggetto:
;
;    [size]asm -f elf32 fss32.nasm 
;
%include "sseutils32.nasm"

section .data			; Sezione contenente dati inizializzati
	unroll equ 4
	dim    equ 4
	sessanta_cinque_int equ 65
	; Costanti di Ramachandran
	alpha_phi	dd -57.8
	alpha_psi	dd -47.0
	beta_phi	dd -119.0
	beta_psi	dd 113.0
	un_mezzo    dd 0.5
	align 4
	dieci       dd 10.0
	sessanta_cinque dd 65.0
	zero 		dd 0.0
	due         dd 2.0
	venti_quattro dd 24.0
	settecento_venti dd 720.0
	uno         dd 1.0
	sei         dd 6.0
	cento_venti dd 120.0
	cinquemila_quaranta dd 5040.0
	meno_uno    dd -1.0
	tmpStamp 	dd 0.0
	extern size 
	msg db "pack: "
	

	; Hydrophobicity
	alignb 16
	hydrophobicity1 dd 1.8, -1, 2.5, -3.5, -3.5, 2.8, -0.4, -3.2, 4.5, -1, -3.9, 3.8, 1.9, -3.5, -1, -1.6, -3.5, -4.5, -0.8, -0.7, -1, 4.2, -0.9, -1, -1.3, -1

	alignb 16 
	; Volume
	volume1 dd 88.6, -1, 108.5, 111.1, 138.4, 189.9, 60.1, 153.2, 166.7, -1, 168.6, 166.7, 162.9, 114.1, -1, 112.7, 143.8, 173.4, 89.0, 116.1, -1, 140.0, 227.8, -1, 193.6, -1
	
	alignb 16
	; Charge
	charge1 dd 0, -1, 0, -1, -1, 0, 0, 0.5, 0, -1, 1, 0, 0, 0, -1, 0, 0, 1, 0, 0, -1, 0, 0, -1, 0, -1

	align 16
	intArray dd 19, 21, 2, 15, 24, 5, 4, 15, 12, 7, 0, 11, 13, 24, 5, 13, 3, 5, 2, 7, 0, 15, 8, 15, 18, 17, 3, 3, 22, 5, 17, 18, 18, 17, 3, 22, 5, 13, 13, 22, 3, 17, 21, 13, 15, 2, 11, 3, 17, 19, 15, 10, 19, 8, 24, 0, 11, 0, 6, 24, 17, 0, 0, 10, 7, 11, 13, 7, 24, 6, 24, 7, 18, 12, 11, 22, 2, 11, 21, 12, 3, 7, 12, 18, 0, 18, 11, 22, 17, 11, 2, 5, 22, 21, 7, 19, 8, 4, 13, 11, 0, 3, 11, 12, 17, 12, 24, 10, 17, 4, 16, 5, 5, 21, 15, 22, 5, 24, 13, 13, 21, 0, 13, 24, 18, 22, 24, 2, 15, 16, 5, 22, 24, 7, 16, 13, 7, 22, 6, 0, 17, 15, 24, 16, 10, 0, 22, 6, 3, 2, 0, 16, 12, 19, 7, 24, 12, 22, 19, 21, 2, 6, 21, 7, 5, 16, 12, 13, 21, 11, 4, 6, 10, 2, 12, 5, 4, 5, 0, 11, 24, 19, 21, 24, 6, 6, 22, 19, 18, 19, 3, 11, 11, 0, 7, 19, 15, 24, 4, 19, 18, 4, 22, 4, 4, 3, 17, 7, 4, 19, 21, 3, 2, 21, 7, 11, 2, 22, 11, 7, 4, 13, 24, 12, 7, 2, 2, 8, 5, 16, 12, 0, 3, 12, 4, 22, 3, 0, 4, 18, 6, 4, 17, 10, 19, 3, 24, 6, 10, 15, 11, 10, 21, 11, 16, 16

section .bss			; Sezione contenente dati[size]on inizializzati
	alignb 16
	e		resd		1

	alignb 16
	dist   resd        1

section .text			; Sezione contenente il codice macchina


; ----------------------------------------------------------
; macro per l'allocazione dinamica della memoria
;
;	getmem	<[size]>,<elements>
;
; alloca un'area di memoria di <[size]>*<elements> bytes
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
global approx_cos
global prodottoScalare
;global approx_sin
;global hydrofobic_energy
;global electrostatic_energy
global packing_energy


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
;	int[size];			// lunghezza sequenza
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
	push 	ecx     
	push 	eax           
	
	; Carica i parametri    
	mov     esi, [ebp+8]               ; coordinateCa    
	mov     eax, [ebp+12]              ; i    
	mov     ebx, [ebp+16]              ; j    
	mov     ecx, [ebp+20]              ; dist
	

	; xor edi, edi
	; xor edx, edx
	; stampa:
	; 	mov edx, edi
	; 	imul edx, 3
	; 	cmp edx, size
	; 	jge nostampa
	; 	mov edx, esi
	; 	add edx, edi
	; 	movss xmm0, dword [edx] ; mette in xmm0 coords[i*9+5] salvando z
	; 	movss [tmpStamp], xmm0
	; 	printss tmpStamp
	; 	inc edi
	; 	jmp stampa
	; 	nostampa:
	; 	xor edi, edi
	; 	xor edx, edx
	


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
	
	; eax = i, ebx = j, esi = cacoords
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
	;mov 	eax, [ebp+20]				;[ebp+20] contiene l'indirizzo della variabile dist passata come parametro (&dist)
	movss 	[ecx], xmm0					;[eax] inserisce[size]ell'indirizzo passato in eax il valore risultante in xmm0

	pop eax
	pop ecx
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
;     const int[size] =[size];
; 	for (int i = 0; i <[size]; i++) {
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
	push 	eax
	push	edx
	push 	ecx

	mov ebx, [ebp+8]    ;coords
	mov eax, [ebp+12]	;cacoords

		 
		; ;Stampa del valore intero size convertito in valore float
		; mov edx, [size]    ; Carica il valore di 'size' in eax
		; cvtsi2ss xmm1, edx ; Converte 'size' (in eax) in float in xmm0
		; movss [tmpStamp], xmm0
		; printss tmpStamp
		; xor edx, edx


		 
	xor esi, esi 		;ESI: i=0

	
    forCacoords:
		cmp esi,[size]
		jge fineCacoords
		
		mov ecx, esi ;ecx contatore di coords
		imul ecx, 9	;i*9
		

		mov edx, esi ; edx contatore di cacords
		shl edx, 1
		add edx, esi 
		;ebx = coords, ecx= i
		movss xmm0, [ebx + ecx*dim + 3*dim] ; mette in xmm0 coords[i*9+3] salvando x
		movss [eax + edx*dim], xmm0
		;cvtsi2ss xmm6, esi ; Converte 'size' (in eax) in float in xmm0
		

        inc edx
		movss xmm0, [ebx + ecx*dim + 4*dim] ; mette in xmm0 coords[i*9+4] salvando y
		movss [eax + edx*dim], xmm0
		

		inc edx
		movss xmm0, [ebx + ecx*dim + 5*dim] ; mette in xmm0 coords[i*9+5] salvando z
		movss [eax + edx*dim], xmm0
		; movss [tmpStamp], xmm0
		; printss tmpStamp
		inc esi
		jmp forCacoords
		;--------fine ciclo for--------

	fineCacoords:
	;--- fine logica ---
	pop ecx
	pop edx
	pop eax
	pop edi
	pop	esi
	pop	ebx
	mov	esp, ebp	; ripristina lo Stack Pointer 
	pop ebp    
	ret

; ------------------------------------------------------------
; Funzione rama_energy
; ------------------------------------------------------------
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
		cmp esi,[size]
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

; ------------------------------------------------------------
; Funzione approx_cos
; ------------------------------------------------------------
approx_cos:    
	push	ebp			; salva il Base Pointer
	mov		ebp, esp	; il Base Pointer punta al Record di Attivazione corrente
	push	ebx			; salva i registri da preservare
	push 	edx
	push	esi
	push	edi
	
	; Recupero parametro theta    
	movss xmm1, [ebp+8] ; theta     
	mulss xmm1, xmm1    ; x2 = theta * theta    
	movss xmm2, xmm1    ; copia di x2 

	; Calcolo x2 / 2.0    
	divss xmm2, [due]   ; xmm2 = x2 / 2.0 

	; Calcolo x2 * x2 / 24.0    
	movss xmm3, xmm1    
	mulss xmm3, xmm3    ; xmm3 = x2 * x2    
	movss xmm4, xmm3    
	divss xmm4, [venti_quattro] ; xmm4 = x2 * x2 / 24.0    
	
	; Calcolo x2 * x2 * x2 / 720.0    
	mulss xmm3, xmm1    ; xmm3 = x2 * x2 * x2    
	movss xmm5, xmm3    
	divss xmm5, [settecento_venti] ; xmm5 = x2 * x2 * x2 / 720.0    
	
	; Combinazione dei risultati per il coseno approssimato    
	movss xmm6, [uno]   ; xmm6 = 1.0    
	subss xmm6, xmm2    ; xmm6 = 1 - (x2 / 2.0)    
	addss xmm6, xmm4    ; xmm6 = 1 - (x2 / 2.0) + (x2 * x2 / 24.0)    
	subss xmm6, xmm5    ; xmm6 = 1 - (x2 / 2.0) + (x2 * x2 / 24.0) - (x2 * x2 * x2 / 720.0)    
	
	; Salva il risultato in [ebp+12] (result) 
	mov eax, [ebp+12] ; result 
	movss [eax], xmm6    
	;movss [ebp+12], xmm6    
	
	pop edi
	pop	esi
	pop edx
	pop	ebx
	mov	esp, ebp	; ripristina lo Stack Pointer 
	pop ebp    
	ret


; ------------------------------------------------------------
; Funzione prod_scal
; ------------------------------------------------------------
prodottoScalare:
	push	ebp			; salva il Base Pointer
	mov		ebp, esp	; il Base Pointer punta al Record di Attivazione corrente
	push	ebx			; salva i registri da preservare
	push 	edx
	push	esi
	push	edi
 
 
	mov ebx, [ebp+8]    ;axis
	mov eax, [ebp+12]		;axis[size]ormalizzato
 

 
		; Calcola il prodotto scalare
	
	movss xmm0, [ebx] ; axis[0]
	mulss xmm0, xmm0
 
	movss xmm1, [ebx+4] ; axis[1]
	mulss xmm1, xmm1
	addss xmm0, xmm1
 
 
	movss xmm2, [ebx+2*4] ; axis[2]
	mulss xmm2, xmm2
	addss xmm0, xmm2

	;sqrtss xmm0, xmm0
 
	; xmm0 ha il prodotto scalare
 
	; Calcola 1/prodotto scalare
	movss xmm1, [ebx] ; axis[0]
	divss xmm1, xmm0
	movss xmm2, [ebx+dim] ; axis[1]
	divss xmm2, xmm0
	movss xmm3, [ebx+2*dim] ; axis[2]
	divss xmm3, xmm0
 
	movss [eax], xmm0

	;movss [eax], xmm1 ;[size]ew axis[0]
	;movss [eax+dim], xmm2  ;[size]ew axis[1]
	;movss [eax+2*dim], xmm3  ;[size]ew axis[2]
 
	pop edi
	pop	esi
	pop edx
	pop	ebx
	mov	esp, ebp	; ripristina lo Stack Pointer 
	pop ebp    
	ret

; ------------------------------------------------------------
; Funzione approx_sin
; ------------------------------------------------------------
;approx_sin:
	; push	ebp			; salva il Base Pointer
	; mov		ebp, esp	; il Base Pointer punta al Record di Attivazione corrente
	; push	ebx			; salva i registri da preservare
	; push 	edx
	; push	esi
	; push	edi

	; ; Calcolo del seno
	; movss xmm7, [ebp+8] ; theta
	; movss xmm1, xmm7
	; movss xmm5, xmm7 ; theta
	; mulss xmm1, xmm1 ; theta^2
	; movss xmm2, xmm7
	; mulss xmm2, xmm1 ; theta^3
	; movss xmm3, xmm2
	; divss xmm3, [sei]; theta^3 / 6.0
	; subss xmm7, xmm3 ; risultato parziale theta - theta^3 / 6.0

	; movss xmm6, xmm1 ; theta^2
	; mulss xmm6, xmm1 ; theta^4
	; mulss xmm6, xmm5 ; theta^5
	; movss xmm5, xmm6 ; theta^5
	; divss xmm6, [cento_venti] ; theta^5 / 120.0
	; addss xmm7, xmm6 ; risultato parziale theta - theta^3 / 6.0 + theta^5 / 120.0

	; mulss xmm5, xmm1 ; theta^7
	; divss xmm5, [cinquemila_quaranta] ; theta^7 / 5040.0 
	; subss xmm7, xmm5 ; risultato finale theta - theta^3 / 6.0 + theta^5 / 120.0 - theta^7 / 5040.0

	; mulss xmm7, [meno_uno]

	; mov eax, [ebp+12]
	; movss [eax], xmm7

	; pop edi
	; pop	esi
	; pop edx
	; pop	ebx
	; mov	esp, ebp	; ripristina lo Stack Pointer 
	; pop ebp    
	; ret

; ------------------------------------------------------------
; Funzione hydro_energy
; ------------------------------------------------------------

hydrofobic_energy:
	push	ebp			; salva il Base Pointer
	mov		ebp, esp	; il Base Pointer punta al Record di Attivazione corrente
	push	ebx			; salva i registri da preservare
	push 	edx
	push	esi
	push	edi

		;INPUT
	mov ebx, [ebp + 8] ;ebx = s
	mov ecx, [ebp + 12] ;ecx = cacoords

;Stampa del valore intero size convertito in valore float

	; xor esi, esi
	; mov esi, 21
	; imul esi, 4
	; movss xmm6, [charge1]
	; cvtsi2ss xmm7, esi
	; addss xmm7, xmm6
	; 	;movss xmm7, [edx]
	; 	movss [tmpStamp], xmm7
	; 	printss tmpStamp
	; stampa:
	; 	cmp esi, size
	; 	jge nostampa
	; 	movzx edx, byte [charge1 + esi*4]    ; Carica il valore di 'size' in eax
	; 	cvtsi2ss xmm7, edx ; Converte 'size' (in eax) in float in xmm0
	; 	;movss xmm7, [edx]
	; 	movss [tmpStamp], xmm7
	; 	printss tmpStamp
	; 	nostampa:
	; 		xor edx, edx
			; xor esi, esi

	; cvtsi2ss xmm7, ebx 
	; movups [tmpStamp], xmm7
	; printss tmpStamp

	xor esi, esi 	; esi = i = 0
	pxor xmm3, xmm3 ;energy = 0

	loopEsterno: 
		cmp esi,[size]
		jge fineloopEsterno
		xor edi, edi ;edi = j = 0
		mov edi, esi ;edi = i
		inc edi		 ;edi = i+1

		loopInterno: 
			cmp edi, [size]
			jge fineloopInterno

			;Chiamata alla funzione distanza
			push dist ; &dist
			push edi ;j
			push esi ;i
			push ecx ;cacoords

			call distanza1
			;fase di svuotamento dello stack
			add esp, 16
			;pxor xmm0, xmm0
			; ecx = distanza1
			movss xmm0, [dist] ;xmm0 = dist
			
			
			
			;if (dist <10.0)
			comiss xmm0, [dieci]
			ja incj
			;dist<10
			;movss [tmpStamp], xmm0
			;printss tmpStamp

			mov edx, [intArray + esi  * dim]
			;sub edx, [sessanta_cinque]
			; cvtsi2ss xmm6, [hydrophobicity1+4]
			; movups [tmpStamp], xmm6
			; printss tmpStamp
			movss xmm1, [hydrophobicity1 + edx * dim]
			; movups [tmpStamp], xmm1
			; printss tmpStamp


			

			 mov edx, [intArray + edi  * dim]
			; sub edx, [sessanta_cinque]
			 movss xmm2, [hydrophobicity1 + edx *dim]
			

			;energy += (hydrophobicity[(int)s[i]-65]*hydrophobicity[(int)s[j]-65])/(dist);
			mulss xmm2, xmm1
			divss xmm2, xmm0

			addss xmm3, xmm2

			incj: 
				inc edi
				jmp loopInterno
				fineloopInterno:
					inc esi
					jmp loopEsterno
					fineloopEsterno:
						mov eax, [ebp+16]
						movss [eax], xmm3

	pop edi
	pop	esi
	pop edx
	pop	ebx
	mov	esp, ebp	; ripristina lo Stack Pointer 
	pop ebp    
	ret


; ------------------------------------------------------------
; Funzione elec_energy
; ------------------------------------------------------------

electrostatic_energy:
	push	ebp			; salva il Base Pointer
	mov		ebp, esp	; il Base Pointer punta al Record di Attivazione corrente
	push	ebx			; salva i registri da preservare
	push 	edx
	push 	ecx
	push	eax
	push	esi
	push	edi

	;INPUT
	mov ebx, [ebp + 8] ;ebx = s
	mov ecx, [ebp + 12] ;ecx = cacoords

	; cvtsi2ss xmm7, ebx 
	; movups [tmpStamp], xmm7
	; printss tmpStamp

	xor esi, esi 	; esi = i = 0
	pxor xmm3, xmm3 ;energy = 0

	; cvtsi2ss xmm6, esi ; Converte 'size' (in eax) in float in xmm0
	; 	movss [tmpStamp], xmm6
	; 	printss tmpStamp

	fori: 
		cmp esi,[size]
		jge finefori
		xor edi, edi ;edi = j = 0
		mov edi, esi ;edi = i
		inc edi		 ;edi = i+1

		; ;Stampa del valore intero size convertito in valore float
		; mov edx, [size]    ; Carica il valore di 'size' in eax
		; cvtsi2ss xmm7, edi ; Converte 'size' (in eax) in float in xmm0
		; movss [tmpStamp], xmm7
		; printss tmpStamp
		; xor edx, edx

		forj: 
			cmp edi,[size]
			jge fineforj

			;Chiamata alla funzione distanza
			push dist ; &dist
			push edi ;j
			push esi ;i
			push ecx ;cacoords

			call distanza1
			;fase di svuotamento dello stack
			add esp, 16
			;pxor xmm0, xmm0
			movss xmm0, [dist] ;xmm0 = dist
			; stampa:
			; cmp edi, 20
			; jge nostampa
			; prints msg
			; movss [tmpStamp], xmm0
			; printss tmpStamp
			; nostampa:
			
		; cvtsi2ss xmm7, edi ; Converte 'size' (in eax) in float in xmm0
			
			
		
			;if (dist <10.0)
			comiss xmm0, [dieci]
			jae incrementoj

			
			
			; cvtsi2ss xmm7, edi 
			; movss [tmpStamp], xmm7
			; printss tmpStamp
				

			;charge[s[i]-65]			
			mov edx, [ebx + esi  * dim]

			; cvtsi2ss xmm7, edx 
			; movss [tmpStamp], xmm7
			; printss tmpStamp
			;sub edx, 65
			; cvtsi2ss xmm7, edx 
			; movss [tmpStamp], xmm7
			; printss tmpStamp
		
			
			; movss xmm4, [ebx + esi  * dim]
			; subss xmm4, [sessanta_cinque]
			; cvtss2si edx, xmm4
			; prints msg1
			; cvtsi2ss xmm5, edx 
			; movss [tmpStamp], xmm5
			; printss tmpStamp
			movss xmm1, [charge1 + edx * dim]
			; movss [tmpStamp], xmm0
			; printss tmpStamp
			
			

			;if charge[s[i]-65] != 0
			comiss xmm1, [zero]
			je incrementoj 

			;charge[s[j]-65]
			mov edx, [ebx + edi  * dim]
			;sub edx, [sessanta_cinque]
			; cvtss2si edx, xmm6
			; cvtsi2ss xmm5, edx 
			; movss [tmpStamp], xmm5
			; printss tmpStamp
			movss xmm2, [charge1 + edx *dim]
			
			;if charge[s[j]-65] !=0
			comiss xmm2, [zero]
			je incrementoj 

			;energy += (charge[(int)s[i]-65]*charge[(int)s[j]-65])/(dist*4.0);
			mulss xmm2, xmm1
			mulss xmm0, [dim]
			divss xmm2, xmm0

			addss xmm3, xmm2

			incrementoj: 
				inc edi
				jmp forj
				fineforj:
					inc esi
					jmp fori
					finefori:
						mov eax, [ebp+16]
						movss [eax], xmm3

	pop edi
	pop	esi
	pop eax
	pop ecx
	pop edx
	pop	ebx
	mov	esp, ebp	; ripristina lo Stack Pointer 
	pop ebp    
	ret

; ------------------------------------------------------------
; Funzione packing_energy
; ------------------------------------------------------------
; extern void packing_energy(char*s, MATRIX cacoords, type *pack); 
; /*{
;     const int[size] =[size]; 
;     type energy = 0.0;
;     for (int i = 0; i <[size]; i++) {
; 		type  density = 0.0;
; 		for (int j = 0; j <[size]; j++) {
; 			if(i != j){
; 				//type dist = distanza(cacoords, i, j);
; 				type dist = 0.0;
; 				distanza1(cacoords, i, j, &dist);
; 				if (dist < 10.0) {
; 					density  += volume[(int)s[j]-65] / (dist * dist * dist); 
; 				}
; 			}
; 		}
; 		energy  += ((volume[(int)s[i]-65] - density) * (volume[(int)s[i]-65] - density));
;     }
; 	//printf("energy pack %f\n", energy);
; 	*pack = energy;
; 	return ;
; }*/

packing_energy:
	push	ebp			; salva il Base Pointer
	mov		ebp, esp	; il Base Pointer punta al Record di Attivazione corrente
	push	ebx			; salva i registri da preservare
	push 	edx
	push	esi
	push	edi

	;INPUT
	mov ebx, [ebp + 8] ;ebx = s
	mov ecx, [ebp + 12] ;ecx = cacoords

		; movss  xmm4, dword [ecx]    
		; movss [tmpStamp], xmm4
		; printss tmpStamp

	;stampa il primo valore di ebx, cioé della sequenza
	; movzx eax, byte [ebx]
	; sub eax, 65
	; cvtsi2ss xmm5, eax 
	; movss [tmpStamp], xmm5
	; printss tmpStamp
	; xor eax, eax

	xor esi, esi 	; esi = i = 0
	pxor xmm3, xmm3 ;energy = 0

	for_i: 
		cmp esi,[size]
		jge fine_fori
		pxor xmm4, xmm4 ; densità
		xor edi, edi ;edi = j = 0
		

		for_j: 
			cmp edi,[size]
			jge fine_forj
			cmp edi, esi ; if j==i
			je incremento_j
			;Chiamata alla funzione distanza
			push eax ; &dist
			push edi ;j
			push esi ;i
			push ecx ;cacoords

			call distanza1
			;fase di svuotamento dello stack
			add esp, 16
			;pxor xmm0, xmm0
			movss xmm0, [eax] ;xmm0 = dist
			
			
			;if (dist <10.0)
			comiss xmm0, [dieci]
			jae incremento_j
			;dist<10
			
			;Stampa distanze minori di dieci (corrette)
			; movss [tmpStamp], xmm0
			; printss tmpStamp


			;volume[s[j]]
			xor edx, edx
			movzx edx, byte [ebx + edi]
			;mov edx, [ebx + esi  * dim]
			sub edx, 65
			 cvtsi2ss xmm5, edx 			;xmm5 = s[j]
			; movss [tmpStamp], xmm5
			; printss tmpStamp
			movss xmm1, [volume1 + edx * 4]  ;xmm1 = volume[s[j]]
			; movss [tmpStamp], xmm1
			; printss tmpStamp
			movss xmm7, xmm0
			mulss xmm7, xmm0

			
			mulss xmm7, xmm0 ; xmm0= dist^3
			; movss [tmpStamp], xmm7
			; printss tmpStamp
			movss xmm7, xmm7
			divss xmm1, xmm7 
			addss xmm4, xmm1 ; densità +=
			; movss [tmpStamp], xmm4
			; printss tmpStamp
			
			incremento_j: 
				inc edi
				jmp for_j
				fine_forj:
					xor edx, edx
					movzx edx, byte [ebx + esi]
					; ;mov edx, [ebx + esi  * dim]
					 sub edx, 65
					; cvtsi2ss xmm5, edx 
					; movss [tmpStamp], xmm5
					; printss tmpStamp
					
					movss xmm1, [volume1 + edx * 4]
					subss xmm1, xmm4; volume-densità
					mulss xmm1, xmm1
					addss xmm3, xmm1; energia+=
					inc esi

					jmp for_i
					fine_fori:
						mov eax, [ebp+16]
						movss [eax], xmm3
						;prints msg
						

	pop edi
	pop	esi
	pop edx
	pop	ebx
	mov	esp, ebp	; ripristina lo Stack Pointer 
	pop ebp    
	ret