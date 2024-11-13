
#include "syspara.h"

typedef double Number;

void initial_mem()
{

	// ina_fast
	TkC2O=(Number *)calloc(VNMAX,sizeof(Number));
	TkOC=(Number *)calloc(VNMAX,sizeof(Number));
	TkOI2=(Number *)calloc(VNMAX,sizeof(Number));
	TkIsb=(Number *)calloc(VNMAX,sizeof(Number));
	TkIsf=(Number *)calloc(VNMAX,sizeof(Number));
	TfC=(Number *)calloc(VNMAX,sizeof(Number));
	if( TkC2O == NULL || TkOC == NULL || TkOI2 == NULL ||
		TkIsb == NULL || TkIsf == NULL || TfC == NULL ) exit(1);
	
	// iKur
	Tass=(Number *)calloc(VNMAX,sizeof(Number));
	Tiss=(Number *)calloc(VNMAX,sizeof(Number));
	Ttau_aur=(Number *)calloc(VNMAX,sizeof(Number));
	Ttau_iur=(Number *)calloc(VNMAX,sizeof(Number));
	if( Tass ==NULL || Ttau_aur==NULL 
		|| Tiss==NULL || Ttau_iur==NULL ) exit(1);

	// ito
	Tpr_inf=(Number *)calloc(VNMAX,sizeof(Number));
	Ts_inf=(Number *)calloc(VNMAX,sizeof(Number));
	Tss_inf=(Number *)calloc(VNMAX,sizeof(Number));
	Ttau_pr=(Number *)calloc(VNMAX,sizeof(Number));
	Ttau_s=(Number *)calloc(VNMAX,sizeof(Number));
	Ttau_ss=(Number *)calloc(VNMAX,sizeof(Number));
	if( Tpr_inf ==NULL || Ttau_pr==NULL 
		|| Ts_inf==NULL || Ttau_s==NULL 
		|| Tss_inf==NULL || Ttau_ss==NULL ) exit(1);
	
	// ical
	TfVmACT=(Number *)calloc(VNMAX,sizeof(Number));
	TexpF=(Number *)calloc(VNMAX,sizeof(Number));
	Tkooco_iz=(Number *)calloc(VNMAX,sizeof(Number));
	Tkccco_iz=(Number *)calloc(VNMAX,sizeof(Number));
	Tkooco_blk=(Number *)calloc(VNMAX,sizeof(Number));
	Tkccco_blk=(Number *)calloc(VNMAX,sizeof(Number));
	if( TfVmACT == NULL || TexpF == NULL || Tkooco_iz == NULL || 
		Tkccco_iz == NULL || Tkooco_blk == NULL || Tkccco_blk == NULL ) exit(1);

	// ikr
	Tay1=(Number *)calloc(VNMAX,sizeof(Number));
	Tay2=(Number *)calloc(VNMAX,sizeof(Number));
	Tay3=(Number *)calloc(VNMAX,sizeof(Number));
	Tby1=(Number *)calloc(VNMAX,sizeof(Number));
	Tby2=(Number *)calloc(VNMAX,sizeof(Number));
	Tby3=(Number *)calloc(VNMAX,sizeof(Number));
	if( Tay1 == NULL || Tay2 == NULL || Tay3 == NULL || Tby1 == NULL || Tby2 == NULL || Tby3 == NULL) exit(1);

	// iClh 
	Ta1_Clh=(Number *)calloc(VNMAX,sizeof(Number));
	Ta2_Clh=(Number *)calloc(VNMAX,sizeof(Number));
	Tb1_Clh=(Number *)calloc(VNMAX,sizeof(Number));
	Tb2_Clh=(Number *)calloc(VNMAX,sizeof(Number));
	if( Ta1_Clh == NULL || Ta2_Clh == NULL || Tb1_Clh == NULL || Tb2_Clh == NULL ) exit(1);

	// INCX
	Tk1_NCX = (Number *)calloc(VNMAX,sizeof(Number));
	Tk2_NCX = (Number *)calloc(VNMAX,sizeof(Number));
	if( Tk1_NCX == NULL || Tk2_NCX == NULL ) exit(1);

	// RyR 
	Tkvoc  = (Number *)calloc(VNMAX,sizeof(Number));
	Tkcaco = (Number *)calloc(VNMAX,sizeof(Number));
	if( Tkvoc == NULL || Tkcaco == NULL ) exit(1);
	// ikb	
	// icab
	// inab
}


void closed_mem()
{

	free(TkC2O);free(TkOC);free(TkOI2);free(TkIsb);free(TkIsf);free(TfC);

	free(Tass); free(Tiss); free(Ttau_aur); free(Ttau_iur);

	free(Tpr_inf); free(Ts_inf); free(Tss_inf); free(Ttau_pr); free(Ttau_s);free(Ttau_ss); 

	free(TfVmACT);free(TexpF);free(Tkooco_iz);free(Tkccco_iz);free(Tkooco_blk);free(Tkccco_blk);

	free(Tay1); free(Tay2); free(Tay3); free(Tby1); free(Tby2); free(Tby3); 

	free(Ta1_Clh); free(Ta2_Clh); free(Tb1_Clh); free(Tb2_Clh);

	free(Tk1_NCX); free(Tk2_NCX);
	free(Tkvoc); free(Tkcaco);

}

// メモリ領域が連続な2次元配列
void *malloc2d(size_t size, int row, int col)
{
	char **a, *b;
	int  t = size * col;
	int i;

	// インデックスと要素を一気に確保
	a = (char**)malloc((sizeof(*a) + t) * row);

	if (a) {
		// [インデックス, インデックス, ..., 要素, 要素, 要素, ...]と整列させるため要素の開始位置をずらす
		b = (char*)(a + row);

		// 各行の先頭アドレスを与える
		for (i = 0; i < row; i++) {
			a[i] = b;
			b += t; // 要素のサイズ×列の長さの分だけずらす
		}
		return a;
	}

	return NULL;
}

// メモリ領域が連続な3次元配列
void *malloc3d(size_t size, int i, int j, int k)
{
	char ***a, **b, *c;
	int  t = size * k;
	int idx1,idx2;

	// インデックスと要素を一気に確保
	a = (char***)malloc((sizeof(*a) + sizeof(**a) * j + t * j) * i);

	if (a) {
		b = (char**)(a + i);
		c = (char *)(b + i * j);

		for (idx1 = 0; idx1 < i; idx1++) {
			a[idx1] = b;
			for (idx2 = 0; idx2 < j; idx2++) {
				b[idx2] = c;
				c += t;
			}
			b += j;
		}

		return a;
	}

	return NULL;
}
