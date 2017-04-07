/* arithm.c
 * --------
 *
 * some easy functions.
 */


/* Max
 * ---
 *
 * computes the maximum of a, b and c.
 */

int Max (int a, int b, int c)
{
   int res;
   res=a;

   if (a<b)
      res=b;
   if (res<c)
      res=c;
   return res;
}


/* Max2
 * ----
 *
 * computes the maximum of a and b.
 */

int Max2 (int a, int b)
{
   if (b>a)
      return b;
   else
      return a;
}


/* Min
 * ---
 *
 * computes the minimum of a, b and c.
 */

int Min (int a, int b, int c)
{
   int res;
   res=a;

   if (b<a)
      res=b;
   if (c<res)
      res=c;
   return res;
}


/* Min2
 * ----
 *
 * computes the minimum of a and b.
 */

int Min2 (int a, int b)
{
   if (b<a)
      return b;
   else
      return a;
}


/* Exp_next_power2
 * ---------------
 *
 * computes for a the next bigger number with b = 2^m for one m and returns m.
 */

int Exp_next_power2 (int a)
{
   int b=1;
   int m=0;

   while (b<a) {
      m++;
      b*=2;
   }
   return m;
}


/* Power_to_2
 * ----------
 */

int Power_to_2 (int n)
{
   int res, i;

   if (n==0)
      return 1;
   else
   {
      res=1;
      for (i=1; i<=n; i++)
         res*=2;
   }  
   return res;
}


/* Power_to_8
 * ----------
 */

int Power_to_8 (int n)
{
   int res, i;

   if (n==0)
      return 1;
   else
   {
      res=1;
      for (i=1; i<=n; i++)
         res*=8;
   }
   return res;
}
