; ---------------------------------------------------------
; Regression con istruzioni AVX a 64 bit
; ---------------------------------------------------------
; F. Angiulli, F. Fassetti, S. Nisticò
; 12/11/2024
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
;     nasm -f elf64 pst64.nasm
;

%include "sseutils64.nasm"




section .bss			; Sezione contenente dati non inizializzati

alignb 32
e		resq		1

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
	mov	rdi, %1
	mov	rsi, %2
	call	get_block
%endmacro

%macro	fremem	1
	mov	rdi, %1
	call	free_block
%endmacro

section .data			; Sezione contenente dati inizializzati
	unroll equ 4
	n 	   equ 256
	dim    equ 8
	; Costanti di Ramachandran
	alpha_phi	dq -57.8
	alpha_psi	dq -47.0
	beta_phi	dq -119.0
	beta_psi	dq 113.0
	un_mezzo    dq 0.5
	dieci       dq 10.0
	sessanta_cinque dq 65.0
	zero 		dq 0.0
	due         dq 2.0
	venti_quattro dq 24.0
	settecento_venti dq 720.0
	uno         dq 1.0
	sei         dq 6.0
	cento_venti dq 120.0
	cinquemila_quaranta dq 5040.0
	meno_uno    dq -1.0
	tmp         dd  0.0
	; Hydrophobicity
	alignb 32
	hydrophobicity1 dq 1.8, -1, 2.5, -3.5, -3.5, 2.8, -0.4, -3.2, 4.5, -1, -3.9, 3.8, 1.9, -3.5, -1, -1.6, -3.5, -4.5, -0.8, -0.7, -1, 4.2, -0.9, -1, -1.3, -1

	alignb 32 
	; Volume
	volume1 dq 88.6, -1, 108.5, 111.1, 138.4, 189.9, 60.1, 153.2, 166.7, -1, 168.6, 166.7, 162.9, 114.1, -1, 112.7, 143.8, 173.4, 89.0, 116.1, -1, 140.0, 227.8, -1, 193.6, -1
	
	alignb 32
	; Charge
	charge1 dq 0, -1, 0, -1, -1, 0, 0, 0.5, 0, -1, 1, 0, 0, 0, -1, 0, 0, 1, 0, 0, -1, 0, 0, -1, 0, -1



; ------------------------------------------------------------
; Funzioni
; ------------------------------------------------------------

global distanza1
global coordsca
global rama_energy
global approx_cos
global prodotto_scal
global approx_sin
; global hydrofobic_energy
; global electrostatic_energy
; global packing_energy





; ORDINE INPUT: RDI - RSI - RDX - RCX - R8 - R9

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
section .text
extern size
; La funzione accetta i seguenti argomenti:
; RDI: pointer a coordinateCa (array di float)
; RSI: i
; RDX: j
; R8:  pointer a dist (variabile di output)

distanza1:
    ; Carica i parametri    
    ; coordinateCa è passato in RDI
    ; i è passato in RSI
    ; j è passato in RDX

    push	rbp			; salva il Base Pointer
	mov		rbp, rsp	; il Base Pointer punta al Record di Attivazione corrente
	pushaq  
    
	mov     rax, rsi                ; rax = i
    lea     rax, [rax * 3]          ; rax = 3*i
    mov     rbx, rdx                ; rbx = j
    lea     rbx, [rbx * 3]          ; rbx = 3*j

    ; X_df    
    movsd   xmm0, qword [rdi + rax * dim]  ; xmm0 = coordinateCa[3*i]    
    movsd   xmm1, qword [rdi + rbx * dim]  ; xmm1 = coordinateCa[3*j]    
    subsd   xmm0, xmm1                   ; xmm0 = x_df    

    ; Y_df    
    inc     rax                          ; incrementa 3*i
    inc     rbx                          ; incrementa 3*j
    movsd   xmm2, qword [rdi + rax * dim]  ; xmm2 = coordinateCa[3*i + 1]    
    movsd   xmm3, qword [rdi + rbx * dim]  ; xmm3 = coordinateCa[3*j + 1]    
    subsd   xmm2, xmm3                   ; xmm2 = y_df    

    ; Z_df    
    inc     rax                          ; incrementa 3*i
    inc     rbx                          ; incrementa 3*j
    movsd   xmm4, qword [rdi + rax * dim]  ; xmm4 = coordinateCa[3*i + 2]    
    movsd   xmm5, qword [rdi + rbx * dim]  ; xmm5 = coordinateCa[3*j + 2]    
    subsd   xmm4, xmm5                   ; xmm4 = z_df    

    ; Calcola x_df^2, y_df^2, z_df^2 e somma    
    mulsd   xmm0, xmm0                   ; xmm0 = x_df^2    
    mulsd   xmm2, xmm2                   ; xmm2 = y_df^2    
    mulsd   xmm4, xmm4                   ; xmm4 = z_df^2    
    addsd   xmm0, xmm2                   ; xmm0 += y_df^2    
    addsd   xmm0, xmm4                   ; xmm0 += z_df^2    

    ; Radice quadrata    
    sqrtsd  xmm0, xmm0                   ; xmm0 = sqrt(x_df^2 + y_df^2 + z_df^2)    

    ; Salva il risultato
    movsd   qword [rcx], xmm0             ; salva il risultato in dist

    popaq
	mov	rsp, rbp	; ripristina lo Stack Pointer 
 	pop rbp
	ret                                   ; ritorna



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
	push	rbp			; salva il Base Pointer
	mov		rbp, rsp	; il Base Pointer punta al Record di Attivazione corrente
	pushaq 

	;coords RDI
	;cacoords RSI  
	; rcx= ecx  rbx= esi   
	xor rbx, rbx 		;ESI: i=0
	
	
    forCacoords:
		
		cmp rbx, [size]
		jge fineCacoords
		mov rcx, rbx ;ecx contatore di coords
		imul rcx, 9
		

		mov rdx, rbx ; rdx contatore di cacords
		shl rdx, 1
		add rdx, rbx 
		
		movsd xmm0, qword[rdi + rcx*dim + 3*dim] ; mette in xmm0 coords[i*9+3] salvando x
		movsd [rsi + rdx*dim], xmm0

        inc rdx
		movsd xmm0, qword[rdi + rcx*dim + 4*dim] ; mette in xmm0 coords[i*9+4] salvando y
		movsd [rsi + rdx*dim], xmm0

		inc rdx
		movsd xmm0, qword[rdi + rcx*dim + 5*dim] ; mette in xmm0 coords[i*9+5] salvando z
		movsd [rsi + rdx*dim], xmm0
		inc rbx
		jmp forCacoords
		;--------fine ciclo for--------

	fineCacoords:
	;--- fine logica ---
	popaq
	mov	rsp, rbp	; ripristina lo Stack Pointer 
 	pop rbp
	ret

; ; ------------------------------------------------------------
; ; Funzione rama_energy
; ; ------------------------------------------------------------
rama_energy: 
	push	rbp			; salva il Base Pointer
	mov		rbp, rsp	; il Base Pointer punta al Record di Attivazione corrente
	pushaq  

    ;mov ebx, [ebp+8]    	;phi  RDI
	;mov ecx, [ebp+12]		;psi  RSI
	; output                      RDX
	
	xor rcx, rcx 			;esi: i=0
	xorps xmm0, xmm0        ;init energy = 0.0
	
    forRamaEnergy:
		cmp rcx, [size]
		jge fineRamaEnergy
		; uso xmm4 per salvare phi[i]
		; uso xmm5 per salvare psi[i]

		;--------calcolo alpha_dist--------
		vmovsd xmm4, [rdi+rcx*dim]   ; recupero phi[i]
		vmovsd xmm6, xmm4 			; salvo phi[i] in xmm6
		vsubsd xmm6, [alpha_phi] 			; phi[i] - alpha_phi
		vmulsd xmm6, xmm6			;(phi[i] - alpha_phi)^2

		vmovsd xmm5, [rsi+rcx*dim]   ; recupero psi[i]
		vmovsd xmm7, xmm5 			; salvo psi[i] in xmm7
		vsubsd xmm7, [alpha_psi]		; psi[i] - alpha_psi
		vmulsd xmm7, xmm7  	        ; (psi[i] - alpha_psi)^2

		vaddsd xmm6, xmm7 			; (phi[i] - alpha_phi)^2 + (psi[i] - alpha_psi)^2
		
		; in xmm6 ho alpha_dist
   
		;--------calcolo beta_dist--------
		vsubsd xmm4, [beta_phi]		; phi[i] - beta_phi
		vmulsd xmm4, xmm4			; (phi[i] - beta_phi)^2

		vsubsd xmm5, [beta_psi]		; psi[i] - beta_psi
		vmulsd xmm5, xmm5			; (psi[i] - beta_psi)^2
		vaddsd xmm4, xmm5 			; (phi[i] - beta_phi)^2 + (psi[i] - beta_psi)^2
		
		; in xmm4 ho beta_dist

		;--------confronto esplicito--------
		comisd xmm6, xmm4 			; xmm6 < xmm4
		jl alpha_dist_minore
		; se distanza alpha maggiore
		sqrtsd xmm6, xmm6 			; sqrt((phi[i] - alpha_phi)^2 + (psi[i] - alpha_psi)^2)
		vmovsd xmm5, [un_mezzo]
		vmulsd xmm6, xmm5	    	; 0.5 * alpha_dist
		vaddsd xmm0, xmm6 			; energy += 0.5 * alpha_dist
		jmp fine_if
		

		;--------alpha_dist < beta_dist--------
		alpha_dist_minore:
		sqrtsd xmm4, xmm4 			; sqrt((phi[i] - beta_phi)^2 + (psi[i] - beta_psi)^2)
		vmulsd xmm4, xmm5			; 0.5 * beta_dist
		vaddsd xmm0, xmm4 			; energy += 0.5 * beta_dist
	
	fine_if:
		inc rcx
		jmp forRamaEnergy

    fineRamaEnergy:
		vmovsd   qword [rdx], xmm0
	;--- fine logica ---

	popaq
	mov	rsp, rbp	; ripristina lo Stack Pointer 
 	pop rbp
	ret

	ret

; ; ------------------------------------------------------------
; ; Funzione approx_cos
; ; ------------------------------------------------------------
approx_cos:    
	push	rbp			; salva il Base Pointer
	mov		rbp, rsp	; il Base Pointer punta al Record di Attivazione corrente
	pushaq  
	
	; theta   XMM0
	; Recupero parametro theta    
	movsd xmm7, xmm0    ; theta     
	mulsd xmm7, xmm7    ; x2 = theta * theta    
	movsd xmm2, xmm7    ; copia di x2 

	; Calcolo x2 / 2.0    
	divsd xmm2, [due]   ; xmm2 = x2 / 2.0 

	; Calcolo x2 * x2 / 24.0    
	movsd xmm3, xmm7    
	mulsd xmm3, xmm3    ; xmm3 = x2 * x2    
	movsd xmm4, xmm3    
	divsd xmm4, [venti_quattro] ; xmm4 = x2 * x2 / 24.0    
	
	; Calcolo x2 * x2 * x2 / 720.0    
	mulsd xmm3, xmm7    ; xmm3 = x2 * x2 * x2    
	movsd xmm5, xmm3    
	divsd xmm5, [settecento_venti] ; xmm5 = x2 * x2 * x2 / 720.0    
	
	; Combinazione dei risultati per il coseno approssimato    
	movsd xmm6, [uno]   ; xmm6 = 1.0    
	subsd xmm6, xmm2    ; xmm6 = 1 - (x2 / 2.0)    
	addsd xmm6, xmm4    ; xmm6 = 1 - (x2 / 2.0) + (x2 * x2 / 24.0)    
	subsd xmm6, xmm5    ; xmm6 = 1 - (x2 / 2.0) + (x2 * x2 / 24.0) - (x2 * x2 * x2 / 720.0)    
	
	; Salva il risultato in [ebp+12] (result) 
	; result 
	movsd qword[rdi], xmm6    
	    
	popaq
	mov	rsp, rbp	; ripristina lo Stack Pointer 
 	pop rbp
	   
	ret


; ; ------------------------------------------------------------
; ; Funzione approx_sin
; ; ------------------------------------------------------------
approx_sin:
	push	rbp			; salva il Base Pointer
	mov		rbp, rsp	; il Base Pointer punta al Record di Attivazione corrente
	pushaq  

	; Calcolo del seno
	movsd xmm7, xmm0 ; theta
	movsd xmm1, xmm7
	movsd xmm5, xmm7 ; theta
	mulsd xmm1, xmm1 ; theta^2
	movsd xmm2, xmm7
	mulsd xmm2, xmm1 ; theta^3
	movsd xmm3, xmm2
	divsd xmm3, [sei]; theta^3 / 6.0
	subsd xmm7, xmm3 ; risultato parziale theta - theta^3 / 6.0

	movsd xmm6, xmm1 ; theta^2
	mulsd xmm6, xmm1 ; theta^4
	mulsd xmm6, xmm5 ; theta^5
	movsd xmm5, xmm6 ; theta^5
	divsd xmm6, [cento_venti] ; theta^5 / 120.0
	addsd xmm7, xmm6 ; risultato parziale theta - theta^3 / 6.0 + theta^5 / 120.0

	mulsd xmm5, xmm1 ; theta^7
	divsd xmm5, [cinquemila_quaranta] ; theta^7 / 5040.0 
	subsd xmm7, xmm5 ; risultato finale theta - theta^3 / 6.0 + theta^5 / 120.0 - theta^7 / 5040.0

	mulsd xmm7, [meno_uno]

	
	movsd qword[rdi], xmm7

	popaq
	mov	rsp, rbp	; ripristina lo Stack Pointer 
 	pop rbp
	ret


; ; ------------------------------------------------------------
; ; Funzione prod_scal
; ; ------------------------------------------------------------
prodotto_scal:

	push	rbp			; salva il Base Pointer
	mov		rbp, rsp	; il Base Pointer punta al Record di Attivazione corrente
	pushaq 
 
	movsd xmm0, [rdi]    ;axis 0                 ;axis rdi
	movsd xmm1, [rdi+8]	 ;axis 1
	movsd xmm2, [rdi+16] ;axis 2
			             ;res RDI 
 
	; Calcola il prodotto scalare
	mulsd xmm0, xmm0
 
 	mulsd xmm1, xmm1
	addsd xmm0, xmm1
 
 	mulsd xmm2, xmm2
	addsd xmm0, xmm2
	; xmm0 ha il prodotto scalare
 	movsd [rsi], xmm0
	
	popaq
	mov	rsp, rbp	; ripristina lo Stack Pointer 
 	pop rbp
	ret



; ; ------------------------------------------------------------
; ; Funzione hydro_energy
; ; ------------------------------------------------------------

; hydrofobic_energy:
; 	push	ebp			; salva il Base Pointer
; 	mov		ebp, esp	; il Base Pointer punta al Record di Attivazione corrente
; 	push	ebx			; salva i registri da preservare
; 	push 	edx
; 	push	esi
; 	push	edi

; 	;INPUT
; 	mov ebx, [ebp + 8] ;ebx = s
; 	mov ecx, [ebp + 12] ;ecx = cacoords

; 	xor esi, esi 	; esi = i = 0
; 	pxor xmm3, xmm3 ;energy = 0

; 	loopEsterno: 
; 		cmp esi, 256
; 		jge fineloopEsterno
; 		xor edi, edi ;edi = j = 0
; 		mov edi, esi ;edi = i
; 		inc edi		 ;edi = i+1

; 		loopInterno: 
; 			cmp edi, 255
; 			jge fineloopInterno

; 			;Chiamata alla funzione distanza
; 			push eax ; &dist
; 			push edi ;j
; 			push esi ;i
; 			push ecx ;cacoords

; 			call distanza1
; 			;fase di svuotamento dello stack
; 			add esp, 16
; 			;pxor xmm0, xmm0
; 			movsd xmm0, [eax] ;xmm0 = dist
			
			
; 			;if (dist <10.0)
; 			comiss xmm0, [dieci]
; 			jge incj

; 			mov edx, [ebx + esi  * dim]
; 			sub edx, [sessanta_cinque]
; 			movsd xmm1, [hydrophobicity1 + edx * dim]

; 			mov edx, [ebx + edi  * dim]
; 			sub edx, [sessanta_cinque]
; 			movsd xmm2, [hydrophobicity1 + edx *dim]
			

; 			;energy += (hydrophobicity[(int)s[i]-65]*hydrophobicity[(int)s[j]-65])/(dist);
; 			mulss xmm2, xmm1
; 			divss xmm2, xmm0

; 			addss xmm3, xmm2

; 			incj: 
; 				inc edi
; 				jmp loopInterno
; 				fineloopInterno:
; 					inc esi
; 					jmp loopEsterno
; 					fineloopEsterno:
; 						mov eax, [ebp+16]
; 						movsd [eax], xmm3

; 	pop edi
; 	pop	esi
; 	pop edx
; 	pop	ebx
; 	mov	esp, ebp	; ripristina lo Stack Pointer 
; 	pop ebp    
; 	ret


; ; ------------------------------------------------------------
; ; Funzione elec_energy
; ; ------------------------------------------------------------

; electrostatic_energy:
; 	push	ebp			; salva il Base Pointer
; 	mov		ebp, esp	; il Base Pointer punta al Record di Attivazione corrente
; 	push	ebx			; salva i registri da preservare
; 	push 	edx
; 	push	esi
; 	push	edi

; 	;INPUT
; 	mov ebx, [ebp + 8] ;ebx = s
; 	mov ecx, [ebp + 12] ;ecx = cacoords

; 	xor esi, esi 	; esi = i = 0
; 	pxor xmm3, xmm3 ;energy = 0

; 	fori: 
; 		cmp esi, n
; 		jge finefori
; 		xor edi, edi ;edi = j = 0
; 		mov edi, esi ;edi = i
; 		inc edi		 ;edi = i+1

; 		forj: 
; 			cmp edi, n
; 			jge fineforj

; 			;Chiamata alla funzione distanza
; 			push eax ; &dist
; 			push edi ;j
; 			push esi ;i
; 			push ecx ;cacoords

; 			call distanza1
; 			;fase di svuotamento dello stack
; 			add esp, 16
; 			;pxor xmm0, xmm0
; 			movsd xmm0, [eax] ;xmm0 = dist
			
			
; 			;if (dist <10.0)
; 			comiss xmm0, [dieci]
; 			jge incrementoj

; 			;if charge[s[i]-65] !=0			
; 			mov edx, [ebx + esi  * dim]
; 			sub edx, [sessanta_cinque]
; 			movsd xmm1, [charge1 + edx * dim]
			

; 			comiss xmm1, [zero]
; 			je incrementoj 

; 			;if charge[s[j]-65] !=0
; 			mov edx, [ebx + edi  * dim]
; 			sub edx, [sessanta_cinque]
; 			movsd xmm2, [charge1 + edx *dim]
			

; 			comiss xmm2, [zero]
; 			je incrementoj 

; 			;energy += (charge[(int)s[i]-65]*charge[(int)s[j]-65])/(dist*4.0);
; 			mulss xmm2, xmm1
; 			mulss xmm0, [dim]
; 			divss xmm2, xmm0

; 			addss xmm3, xmm2

; 			incrementoj: 
; 				inc edi
; 				jmp forj
; 				fineforj:
; 					inc esi
; 					jmp fori
; 					finefori:
; 						mov eax, [ebp+16]
; 						movsd [eax], xmm3

; 	pop edi
; 	pop	esi
; 	pop edx
; 	pop	ebx
; 	mov	esp, ebp	; ripristina lo Stack Pointer 
; 	pop ebp    
; 	ret

; ; ------------------------------------------------------------
; ; Funzione packing_energy
; ; ------------------------------------------------------------
; ; extern void packing_energy(char*s, MATRIX cacoords, type *pack); 
; ; /*{
; ;     const int n = 256; 
; ;     type energy = 0.0;
; ;     for (int i = 0; i < n; i++) {
; ; 		type  density = 0.0;
; ; 		for (int j = 0; j < n; j++) {
; ; 			if(i != j){
; ; 				//type dist = distanza(cacoords, i, j);
; ; 				type dist = 0.0;
; ; 				distanza1(cacoords, i, j, &dist);
; ; 				if (dist < 10.0) {
; ; 					density  += volume[(int)s[j]-65] / (dist * dist * dist); 
; ; 				}
; ; 			}
; ; 		}
; ; 		energy  += ((volume[(int)s[i]-65] - density) * (volume[(int)s[i]-65] - density));
; ;     }
; ; 	//printf("energy pack %f\n", energy);
; ; 	*pack = energy;
; ; 	return ;
; ; }*/

; packing_energy:
; 	push	ebp			; salva il Base Pointer
; 	mov		ebp, esp	; il Base Pointer punta al Record di Attivazione corrente
; 	push	ebx			; salva i registri da preservare
; 	push 	edx
; 	push	esi
; 	push	edi

; 	;INPUT
; 	mov ebx, [ebp + 8] ;ebx = s
; 	mov ecx, [ebp + 12] ;ecx = cacoords

; 	xor esi, esi 	; esi = i = 0
; 	pxor xmm3, xmm3 ;energy = 0

; 	for_i: 
; 		cmp esi, n
; 		jge fine_fori
; 		pxor xmm4, xmm4 ; densità
; 		xor edi, edi ;edi = j = 0
		

; 		for_j: 
; 			cmp edi, n
; 			jge fine_forj
; 			cmp edi, esi ; if j==i
; 			je incremento_j
; 			;Chiamata alla funzione distanza
; 			push eax ; &dist
; 			push edi ;j
; 			push esi ;i
; 			push ecx ;cacoords

; 			call distanza1
; 			;fase di svuotamento dello stack
; 			add esp, 16
; 			;pxor xmm0, xmm0
; 			movsd xmm0, [eax] ;xmm0 = dist
			
			
; 			;if (dist <10.0)
; 			comiss xmm0, [dieci]
; 			jge incremento_j

; 			;if volume[s[i]-65] !=0			
; 			xor edx, edx
; 			mov edx, [ebx + edi  * dim]
; 			sub edx, [sessanta_cinque]
; 			movsd xmm1, [volume1 + edx * dim]
; 			movsd xmm7, xmm0 
; 			mulss xmm7, xmm0
			
; 			mulss xmm7, xmm0 ; xmm0= dist^3
; 			divss xmm1, xmm7 
; 			addss xmm4, xmm1 ; densità +=
			
; 			incremento_j: 
; 				inc edi
; 				jmp for_j
; 				fine_forj:
; 					xor edx, edx
; 					mov edx, [ebx + esi  * dim]
; 					sub edx, [sessanta_cinque]
; 					movsd xmm1, [volume1 + edx * dim]
; 					subss xmm1, xmm4; volume-densità
; 					mulss xmm1, xmm1
; 					addss xmm3, xmm1; energia+=
					
; 					inc esi

; 					jmp for_i
; 					fine_fori:
; 						mov eax, [ebp+16]
; 						movsd [eax], xmm3

; 	pop edi
; 	pop	esi
; 	pop edx
; 	pop	ebx
; 	mov	esp, ebp	; ripristina lo Stack Pointer 
; 	pop ebp    
; 	ret