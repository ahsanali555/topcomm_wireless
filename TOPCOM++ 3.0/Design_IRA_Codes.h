#pragma once
#include <math.h>
#define min(u,v) ((u)<(v)?(u):(v))
typedef struct LISTA_CONTENUTO
{
	int intero; //contiene il grado nella sentinella ed il numero di edge
	int riga;
	int colonna;
}lista_base;
typedef struct LISTA_ELEM
{
	lista_base elemento;
	struct LISTA_ELEM* next_orizzontale;
	struct LISTA_ELEM* next_verticale;
}lista_unita;
typedef struct CHECKNODE
{
	int vn_padre;
	int cn_selfid;
	int tipo_ramo; /*puo' essere SYS oppure PAR*/
}check_lista;
typedef struct CHECKNODE_ELEM
{
	check_lista elemento;
	struct CHECKNODE_ELEM* next;
	struct CHECKNODE_ELEM* prev;
}check_nodo_lista;
typedef struct DUEINTELEM
{
	int posizione;
	int grado;
}dueint_elem;
typedef dueint_elem KEY_S;

/* Design an optimal IRA code*/
int Optimize_IRA(
			int KCOD, 
			const int NCOD, 
			const int MM, 
			const int ndeg, 
			int * deg, 
			double * P, 
			const char * prefix, 
			int & g_syst, 
			int & g_parity, 
			int & d_w1, 
			int & d_w2, 
			int & MNDR, 
			int & MRD, 
			const int maxtr=100, 
			const bool tailbiting=false,
			const bool optimize=false
			);

void init_genrand(unsigned long s);
static void next_state_random(void);
unsigned long genrand_int32(void);
long genrand_int31(void);
void init_by_array(unsigned long* init_key, unsigned long key_length);
float canalehard(double ebn, double rate, int numbits, float* vetts, float* vettd);
//fixed canalesoft(double ebn,double rate,int numbits,fixed* vett);
void lfsr(int periodo, int length, unsigned short* vett);
double Gaussian2();
//fixed AWGN_func(double ebn,double rate,int numbits,complex_s* vett,fixed limite,int dimmod);

void Initialize(const int KCODi, const int NCODi, const int MMi);
void Release();
void inizializza_vettore_pos(lista_unita* pos, int n);
void inserisci_pos(lista_unita* vett_orizz, lista_unita* vett_vert, int riga, int colonna, int numero);
void carica_matrice_H(char* sorgente, lista_unita* variables, lista_unita* checks);
void stampa_pos(lista_unita* variables, lista_unita* checks);
int spanning_tree(lista_unita* variables, lista_unita* checks, int colonna, int g_s, int g_p, int* g);
void spanning_tree_ultimostep(lista_unita* variables, lista_unita* checks, int colonna, int g_s, int g_p, dueint_elem* elenco);
int scandisci_modo1(int* vett, lista_unita* checks);
int scandisci_modo2(int* vett, lista_unita* checks);
int calcola_posizione(int riga, int numero);
void inizializza_lista_cn(check_nodo_lista* testa);
void inserisci_lista_cn(check_nodo_lista* testa, int vn_padre, int cn_selfid, int tipo_ramo);
void cancella_lista_cn(check_nodo_lista* pos);
void imposta_doppia_diagonale(lista_unita* variables,lista_unita* checks);
int imposta_vnode(lista_unita* variables,lista_unita* checks,int vn_numero,int deg,int g_syst,int g_parity,int d_w1,int d_w2,int mndr,int mrd, const int maxtr);
void sistema_sottomatrice(lista_unita* variables,lista_unita* checks,int col,int riga);
void annulla_sottomatrice(lista_unita* variables,lista_unita* checks,int col);
int scegli_riga(int* g,lista_unita* checks);
void rimuovi_ultimo_sottoblocco(lista_unita* variables,lista_unita* checks,int col,int riga_parametro);
void aggiungi_nodo_fittizio(lista_unita* variables,lista_unita* checks);
void rimuovi_nodo_fittizio(lista_unita* variables,lista_unita* checks);



/*
* These are the COMPARISON macros
* Replace these macros by YOUR comparison operations.
* e.g. if you are sorting an array of pointers to strings
* you should define:
*
*	GT(x, y)  as   (strcmp((x),(y)) > 0) 	Greater than
*	LT(x, y)  as   (strcmp((x),(y)) < 0) 	Less than
*	GE(x, y)  as   (strcmp((x),(y)) >= 0) 	Greater or equal
*	LE(x, y)  as   (strcmp((x),(y)) <= 0) 	Less or equal
*	EQ(x, y)  as   (strcmp((x),(y)) == 0) 	Equal
*	NE(x, y)  as   (strcmp((x),(y)) != 0) 	Not Equal
*/
#define LT2(x, y) ((x.grado) < (y.grado))
#define SWAP2(x, y) temp = (x); (x) = (y); (y) = temp
void quicksort2(KEY_S *array, int lo, int hi);
void stampa_matrice(lista_unita * variables, lista_unita * checks);
void scrivi_matrice_H(FILE*f, lista_unita * checks);
void scrivi_matrice_H2(FILE*f,lista_unita * checks);


