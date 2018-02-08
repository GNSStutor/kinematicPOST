//////////////////////////////////////////////////////
//
// ‹ts—ñ‚ğŒvZ‚·‚éŠÖ”
//
//////////////////////////////////////////////////////

#include <stdio.h>
#include <math.h>

//---matrix inversion & determinant---//
//‹ts—ñ‚ÌŒvZ//
int minv(double xx[2034],double eps, int nsv)
{
	int work[4000],i,j,k,l,m,r,iw,s,t,u,v;
	double w,wmax,pivot,api,w1,det;
//	extern double fabs();
//	double a[]={xx[]};
//	l = nsv*2+3;
//	m = nsv*2+3;

	l=nsv;
	m=nsv;

	if(m<2 || m>500 || eps<=0.0)return(999);
	w1 = 1.0;
	for(i=0;i<m;i++)
		work[i]=i;
	for(k=0;k<m;k++)
	{
		wmax=0.0;
		for(i=k;i<m;i++)
		{
			w=fabs(xx[i*l+k]);
			if(w>wmax)
			{
				wmax=w;
				r=i;
			}
		}
		pivot=xx[r*l+k];
		api=fabs(pivot);
		if(api<=eps)
		{
			det=w1;
			return(1);
		}
		w1*=pivot;
		u=k*l;
		v=r*l;
		if(r!=k)
		{
			w1=-w1;
			iw=work[k];
			work[k]=work[r];
			work[r]=iw;
			for(j=0;j<m;j++)
			{
				s=u+j;
				t=v+j;
				w=xx[s];
				xx[s]=xx[t];
				xx[t]=w;
			}
		}
		for(i=0;i<m;i++)
			xx[u+i]/=pivot;
		for(i=0;i<m;i++)
		{
			if(i!=k)
			{
				v=i*l;
				s=v+k;
				w=xx[s];
				if(w!=0.0)
				{
					for(j=0;j<m;j++)
						if(j!=k)xx[v+j]-=w*xx[u+j];
						xx[s]=-w/pivot;
				}
			}
		}
		xx[u+k]=1.0/pivot;
	}
	for(i=0;i<m;i++)
	{
		while(1)
		{
			k=work[i];
			if(k==i)break;
			iw=work[k];
			work[k]=work[i];
			work[i]=iw;
			for(j=0;j<m;j++)
			{
				u=j*l;
				s=u+i;
				t=u+k;
				w=xx[s];
				xx[s]=xx[t];
				xx[t]=w;
			}
		}
	}
	det=w1;

	return(0);
}
