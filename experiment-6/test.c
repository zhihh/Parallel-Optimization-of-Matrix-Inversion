#include <stdio.h>
#include <stdlib.h>

int main()
{
	int(*p)[4] = (int(*)[4])malloc(3 * 4 * sizeof(int));
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			printf("%p\n", &p[i][j]);
		}
	}
    p[0][0] = 3;
    int** m;
    m[0] = p[0];
    m[1] = p[1];
    m[2] = p[2];
    printf("%d", m[0][0]);
	free (p);
	return 0;
}
