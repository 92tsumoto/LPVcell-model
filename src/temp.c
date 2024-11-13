#include "syspara.h"

typedef double Number;
typedef long long Lint;

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
