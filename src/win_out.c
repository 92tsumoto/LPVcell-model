#include <string.h>
#include "syspara.h"

void orbit(int *mode, double m[], double x2)
{
	if(*mode==1)   xhplot(PSET,m[13],m[0],CYAN);
	if(*mode==2)   xhplot(PSET,m[13],m[0],MAGENTA);
} 

void draw_p(int *mode, int P, double x[], double x2)
{
	switch(*mode){
		case 1:
			xhplot(P, x[13],x[0],YELLOW);
			xhplot(POINT,x[13],x[0],WHITE);
			break;
		case 2:	
			xhplot(P, x[13],x[0],RED);
			xhplot(POINT,x[13],x[0],WHITE);
			break;
	}
}

void mouse(int *mode, double x[], double x2)
{
	if (!xddpret.key){
		switch(*mode){
			case 1:	x[13] = xddpret.x; x[0] = xddpret.y; break;
			case 2:	x[13] = xddpret.x; x[0] = xddpret.y; break;
		}
	}
}

