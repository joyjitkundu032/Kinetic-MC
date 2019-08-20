#include<mpi.h>
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include</home/jkundu/mt19937ar.h>
#include</home/jkundu/mt19937ar.c>

/* Lattice structure: every alternate site is a binding site; the lattice is periodic only along Y */

#define Lx 10
#define Ly 100
#define R 30000000000 /* # of steps */
#define N (Lx*Ly)
#define Rin 5 /* sets deposition probability */
#define Rd 30 /* sets diffusion probability */
#define Rout 2 /* sets evaporation probability */
#define dir 4
#define COM 3 /* three species */
#define T 388
#define KB 0.000086173324
#define beta (1.0/KB/T)
#define gap (R/10000)
#define INC 0
#define NS 1 /* # of copies */
#define DA (R/gap)

int NN[N][dir],Lat[N],np_type[COM],ntotal,bind[N],nbinds,last,diff[N],taskid,numtasks;
double Eb[COM],Eint[COM][COM],Rtot,Pin,Pd,Pout;
double tot_time,rho,rhob,rho0,rho1,rho2;
char outfile1[100];
double Rh[DA+1],Rhb[DA+1],Rh0[DA+1],Rh1[DA+1],Rh2[DA+1];
long seedval;

/* Create different threads of random numbers for different cores*/
long seedgen(void)
{
	long s, seed, pid, seconds;
	pid = getpid();
	s = time ( &seconds );
	seed = abs(((s*181)*((pid-83)*359))%103929);
	return seed;
}

void take_input()
{
	seedval=seedgen();
	init_genrand(seedval);
}

int find_site(int x, int y)
{
	/* Given x, y it gives site index */
	x=(x+Lx) % Lx;
	y=(y+Ly) % Ly;
	return x+y*Lx;
}

/* Find nearest neighbors */
void neighbour()
{
	int i,j,x,y;
	
	for(y=0;y<Ly;y++)
	{
		for(x=0;x<Lx;x++)
		{
			i=find_site(x,y);
			j=find_site(x-1,y); NN[i][0]=j;
			j=find_site(x+1,y); NN[i][1]=j;
			j=find_site(x,y-1); NN[i][2]=j;
			j=find_site(x,y+1); NN[i][3]=j;
			if(y == 0)
				NN[i][2]=-10;
			if(y == Ly-1)
				NN[i][3]=-10;
		}
	}
}

/* Initializing the lattice*/
void lat_init()
{
	int i,j;

	for(i=0;i<N;i++)
	{
		Lat[i]=0; 
		bind[i]=0;
		diff[i]=0;
	}
	for(i=INC;i<Ly;i=i+2)
	{
		for(j=0;j<Lx;j=j+2)
			bind[j+i*Lx]=1;
	
	}
	for(i=0;i<COM;i++)
		np_type[i]=0;
	
	ntotal=0; nbinds=0; last=0;
	rho=0; rhob=0; rho0=0; rho1=0; rho2=0;
	tot_time=0;
}

/* Initialize diffusion barriers and binding energies */
void initialize()
{
	int i,j;
	for(i=0;i<COM;i++)
	{
		for(j=0;j<COM;j++)
			Eint[i][j]=0;
	}

    /* Diffusion barriers */
	Eint[0][0]=0.0025; Eint[1][1]=0.02; Eint[2][2]=0.03;
 
	Eint[0][1]=(Eint[0][0]+Eint[1][1])/2.0;
	Eint[1][0]=Eint[0][1]; 
	
	Eint[1][2]=(Eint[1][1]+Eint[2][2])/2.0;
	Eint[2][1]=Eint[1][2];

	Eint[0][2]=(Eint[0][0]+Eint[2][2])/2.0;
	Eint[2][0]=Eint[0][2];

    /*Binding energies for three different species*/
	Eb[0]=0.16; Eb[1]=0.75; Eb[2]=2.00;
}

/* Add any of the three species using the Metropolis rule */
void deposit()
{
	int i,m,site,tmps,tmppt,Rtot_new;
	double E_i,E_f,dE,fac;
	site=floor(Lx*genrand_real3());
	if(Lat[site] != 0)
		return;
	i=floor(COM*genrand_real3());

	if(bind[site] == 1)
	{
		Lat[site]=i+1;
		diff[ntotal]=diff[last];
		diff[last]=site;
		np_type[i]++;
		nbinds++;
		ntotal++;
		last++;
		return;
	}

	E_i=0; 
	E_f=0;
	for(m=0;m<dir;m++)
	{
		tmps=NN[site][m];
		if(tmps >= 0)
		{
			if(Lat[tmps] > 0)
			{	
				tmppt=Lat[tmps]-1;
				E_f=E_f+Eint[i][tmppt];
			}
		}
	}
	//E_f=E_f-bind[site]*Eb[i];
	Rtot=1.0*(Lx*Rin+ntotal*Rd+last*Rout);
	Rtot_new=1.0*(Lx*Rin+(ntotal+1)*Rd+(last+1)*Rout);
	
	dE=(E_f-E_i);
	fac=exp(-1.0*beta*dE)*(Rtot/Rtot_new);
	if(fac >= 1.0)
	{
		Lat[site]=i+1;
		diff[ntotal]=diff[last];
		diff[last]=site;
		ntotal++;
		last++;
	}
	else
	{
		if(genrand_real3() < fac)
		{
			Lat[site]=i+1;
			diff[ntotal]=diff[last];
			diff[last]=site;
			ntotal++;
			last++;
		}
	}
}

/* Diffusion move */
void diffuse()
{
	int v,site,next_site,di,pt,tmppt,tmp,swap,tmps,m,val;
	double en_i,en_f,delE,factor,Rtot_new;

	v=floor(ntotal*genrand_real3());
	site=diff[v];
	di=floor(dir*genrand_real3());
	next_site=NN[site][di];
	if(next_site < 0)
		return;
	if(Lat[next_site] != 0)
		return;

	pt=Lat[site]-1;

	en_i=0;
	for(m=0;m<dir;m++)
	{
		tmps=NN[site][m];
		if(tmps >= 0)
		{
			if(Lat[tmps] > 0)
			{	
				tmppt=Lat[tmps]-1;
				en_i=en_i+Eint[pt][tmppt];
			}
		}
	}
	en_i=en_i-bind[site]*Eb[pt];

	en_f=0;
	for(m=0;m<dir;m++)
	{
		tmps=NN[next_site][m];
		if(tmps >= 0)
		{
			if(Lat[tmps] > 0)
			{	
				if(tmps != site)
				{
					tmppt=Lat[tmps]-1;
					en_f=en_f+Eint[pt][tmppt];
				}
			}
		}
	}
	en_f=en_f-bind[next_site]*Eb[pt];
	
	delE=(en_f-en_i);
	val=0;

	Rtot=1.0*(Lx*Rin+ntotal*Rd+last*Rout);
	Rtot_new=Rtot;
	if(site / Lx == 1)
	{
		if(next_site < Lx)
		{
			Rtot_new=1.0*(Lx*Rin+ntotal*Rd+(last+1)*Rout);
			val=1;
		}
	}
	else if(site < Lx)
	{
		if(next_site / Lx == 1)
		{
			Rtot_new=1.0*(Lx*Rin+ntotal*Rd+(last-1)*Rout);
			val=2;
		}
	}
	factor=exp(-1.0*beta*delE)*(Rtot/Rtot_new);
	if(factor >= 1.0)
	{
		Lat[next_site]=Lat[site]; 
		Lat[site]=0; diff[v]=next_site;
		if(bind[next_site] == 1)
		{
			np_type[pt]++;			
			nbinds++;
		}
		else if(bind[site] == 1)
		{
			np_type[pt]--;			
			nbinds--;
		}	
		if(val==1)
		{
			swap=diff[v];
			diff[v]=diff[last];	
			diff[last]=swap;
			last++;
		}
		else if(val==2)
		{
			swap=diff[v];
			diff[v]=diff[last-1];
			diff[last-1]=swap;
			last=last-1;			
		}
		
	}
	else
	{
		if(genrand_real3() < factor)
		{
			Lat[next_site]=Lat[site]; 
			Lat[site]=0; diff[v]=next_site;

			if(bind[next_site] == 1)
			{
				np_type[pt]++;			
				nbinds++;
			}
			else if(bind[site] == 1)
			{
				np_type[pt]--;			
				nbinds--;
			}	
			if(val==1)
			{
				swap=diff[v];
				diff[v]=diff[last];	
				diff[last]=swap;
				last++;
			}
			else if(val==2)
			{
				swap=diff[v];
				diff[v]=diff[last-1];
				diff[last-1]=swap;
				last=last-1;			
			}
		}
	}	
}

/* Remove particles */
void evaporate()
{
	int i,m,site,pt,tmps,tmppt;
	double Ev,dE,fac,Rtot_new;
	i=floor(last*genrand_real3());
	site=diff[i];
	pt=Lat[site]-1;

	Ev=0.0;
	for(m=0;m<dir;m++)
	{
		tmps=NN[site][m];
		if(tmps >= 0)
		{
			if(Lat[tmps] > 0)
			{	
				tmppt=Lat[tmps]-1;
				Ev=Ev+Eint[pt][tmppt];
			}
		}
	}
	Ev=Ev-bind[site]*Eb[pt];
	Rtot=1.0*(Lx*Rin+ntotal*Rd+last*Rout);
	Rtot_new=1.0*(Lx*Rin+(ntotal-1)*Rd+(last-1)*Rout);

	dE=-Ev;
	fac=exp(-1.0*beta*dE)*(Rtot/Rtot_new);
	if(fac >= 1.0)
	{
		Lat[site]=0; 
		diff[i]=diff[last-1];
		diff[last-1]=diff[ntotal-1];
		last=last-1;
		ntotal=ntotal-1;
		if(bind[site] == 1)
		{
			np_type[pt]=np_type[pt]-1;
			nbinds=nbinds-1;
		}
	}
	else
	{
		if(genrand_real3() < fac)
		{
			Lat[site]=0; 
			diff[i]=diff[last-1];
			diff[last-1]=diff[ntotal-1];
			last=last-1;
			ntotal=ntotal-1;
			if(bind[site] == 1)
			{
				np_type[pt]=np_type[pt]-1;
				nbinds=nbinds-1;
			}
		}
	}
}

void evolve()
{
	int t;	
	double ran;
	
	Rtot=(Lx*Rin+ntotal*Rd+last*Rout);
	Pin=1.0*Lx*Rin/Rtot;
	Pd=1.0*ntotal*Rd/Rtot;
	Pout=1-(Pin+Pd);
    /* 1 MC step is defined as the time when on-an-average N diffusion moves are done */
	tot_time=tot_time+1.0*Rd/Rtot;

	ran=genrand_real3();
	if(ran < Pin)
		deposit();
	else
	{
		if(ran < (Pin+Pd))
			diffuse();
		else
			evaporate();
	}
}

void print_data(int inp)
{
    /* prints the occupancies of different species as a function of time: $1- time, $2- total density, $3- density of bound gas molecules, $4- density of species 0, $5- density of species 1, $6- density of species 2 */
	FILE *fp;
	long int t;
	fp=fopen(outfile1,"w");
	fprintf(fp,"# NS = %d\n",inp);
	for(t=1;t<=DA;t++)
		fprintf(fp,"%ld\t%e\t%e\t%e\t%e\t%e\n",t*gap-gap/2,1.0*Rh[t]/inp,1.0*Rhb[t]/inp,1.0*Rh0[t]/inp,1.0*Rh1[t]/inp,1.0*Rh2[t]/inp);
	fclose(fp);
}

int main (int argc, char *argv[])
{
	int ns,tmp;
	long int it;

	MPI_Init(&argc,&argv);
	MPI_Barrier(MPI_COMM_WORLD);
	/*elapsed_time=-MPI_Wtime();*/
	MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
	MPI_Comm_size(MPI_COMM_WORLD,&numtasks);

	initialize();
	take_input();
	neighbour();
	sprintf(outfile1,"mof2DsqT%dLx%dLy%dRi%dRd%dRo%dEH%1.4lfC%1.4lfW%1.4lfSN%d.dat",T,Lx,Ly,Rin,Rd,Rout,Eb[0],Eb[1],Eb[2],taskid);
	for(ns=1;ns<=NS;ns++)
	{	
		lat_init();
		for(it=1;it<=R;it++)
		{
			while(tot_time<it)
				evolve();
			rho=rho+1.0*ntotal/N; rhob=rhob+4.0*nbinds/N; 
			rho0=rho0+4.0*np_type[0]/N; rho1=rho1+4.0*np_type[1]/N; rho2=rho2+4.0*np_type[2]/N;

			if(it % gap == 0)
			{
				tmp=it/gap;
				Rh[tmp]=Rh[tmp]+1.0*rho/gap;
				Rhb[tmp]=Rhb[tmp]+1.0*rhob/gap;
				Rh0[tmp]=Rh0[tmp]+1.0*rho0/gap;
				Rh1[tmp]=Rh1[tmp]+1.0*rho1/gap;
				Rh2[tmp]=Rh2[tmp]+1.0*rho2/gap;

				rho=0; rhob=0; rho0=0; rho1=0; rho2=0;
			}
		}
		print_data(ns);
	}
	MPI_Finalize();
}
