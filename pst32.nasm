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

input		equ		8

msg	db	'e:',32,0
nl	db	10,0


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
	subss   xmm2, xmm3                  ; xmm2 = y_df    
	
	; Z_df    
	add     eax, 1    
	add     ebx, 1    
	movss   xmm4, dword [esi + eax * 4]    
	movss   xmm5, dword [esi + ebx * 4]    
	subss   xmm4, xmm5                  ; xmm4 = z_df    
	
	; Calcola x_df^2, y_df^2, z_df^2 e somma    
	mulss   xmm0, xmm0                  ; xmm0 = x_df^2    
	mulss   xmm2, xmm2                  ; xmm2 = y_df^2    
	mulss   xmm4, xmm4                  ; xmm4 = z_df^2    
	addss   xmm0, xmm2                  ; xmm0 += y_df^2    
	addss   xmm0, xmm4                  ; xmm0 += z_df^2    
	; Radice quadrata    
	sqrtss  xmm0, xmm0                  ; xmm0 = sqrt(x_df^2 + y_df^2 + z_df^2)    
	movss   dword [ebp-4], xmm0         ; Salva risultato    
	mov     eax, [ebp-4]   
	
	pop edi
	pop	esi
	pop	ebx
	mov	esp, ebp	; ripristina lo Stack Pointer 
	pop ebp    
	ret