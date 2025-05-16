#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "mpi.h"
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
double ran2(long *idum);
double gasdev(long *idum);
int main(int argc, char **argv)
{
        int rank, size;

        MPI_Status status;

        MPI_Init( &argc, &argv );
        MPI_Comm_rank( MPI_COMM_WORLD, &rank );
        MPI_Comm_size( MPI_COMM_WORLD, &size );

        int spr = 100/size;
        int i,j,m,n,length,p,ptl_num,legs;//variable for loop
        int checkeof;
        int e,N,dump1,dump2,dump3,dump4;
        int ts,*s;
        long long t,T1,T2,T3,T4,plustime;
        double dt, gamma, v, gammar, Dr, sigma;
        double F;
        double ete,sum;
        double di_x,di_y,di_z;
        double xi, eta;//parameters for noise
        double k,R0,f,r,rx,ry,rz,ux,uy,uz;//bond potential
        double a,cos, rx1,ry1,rz1,rx2,ry2,rz2,rsq1,rsq2,r1,r2;//angular
        double fx1,fy1,fz1,fx2,fy2,fz2;//angular
        double *x, *y, *z, *dx, *dy, *dz, vx, vy, vz;//basic 
        double tx,ty,tz;
        double *xt, *yt,*zt,*dx1,*dy1,*dz1,*dx2,*dy2,*dz2;
        long idum;//random seed
        char NAME0[100],NAME1[100],NAME2[100],NAME3[100],NAME4[100];//space to name
        FILE *ofp = NULL, *ofp1,*ofp2,*ofp3,*ofp4, *ofpeta, *ofpxdot, *data;

        char name[MPI_MAX_PROCESSOR_NAME];
        int len=MPI_MAX_PROCESSOR_NAME;
        MPI_Get_processor_name( name, &len );
               
        if(argc<4){
                printf("particle_number a v legs");//50 4
                exit(1);
        }
        ptl_num=atof(argv[1]);
        v=atof(argv[2]);
        legs=atof(argv[3]);
        dt=atof(argv[4]);
        gammar=atof(argv[5]);
        F=atof(argv[6]);

        printf("n=%d,v=%.4f,l=%d,dt=%.4f,gar=%.4f,f=%.4f\n",ptl_num,v,legs,dt,gammar,F);
        N=ptl_num*legs+1;

        k=1.0;
        //dt=0.01;//Dr=1.0/gamma*3.0;//r=0.5,kT=1//gammar=2.0*Dr;//eta=sqrt(4.0*Dr*v*v/dt/3.0);
        //gammar=10.0;
        gamma=1.0;
        eta=sqrt(2.0*gammar*v*v/dt);//amplitude of noise
        dump1=(int)(10000.0/dt);//T=(int)(pow(N,2)*gamma/dt);
        dump2=(int)(100.0/dt);//T=(int)(pow(N,2)*gamma/dt);
        dump3=(int)(1.0/dt);
        dump4=(int)(0.01/dt);//T=(int)(pow(N,2)*gamma/dt);
        T1 = (long long)(10000000/dt);
        T2 = (long long)(100000/dt);
        T3 = (long long)(1000/dt);
        T4 = (long long)(10/dt);
        xi=sqrt(2.0*gamma/dt);
        idum=time(NULL)+rank;
        for(p=0;p<30;p++){
                tx = gasdev(&idum);
        }
               
        for(j=rank*spr;j<(rank+1)*spr;j++){
                printf("%s: %d\n",name,j);
                s=(int*) calloc(N,sizeof(int));
                s[0]=1;
                x=(double*) calloc(N+ptl_num,sizeof(double));
                //y=(double*) calloc(N+ptl_num,sizeof(double));
                //z=(double*) calloc(N+ptl_num,sizeof(double));
                dx=(double*) calloc(N,sizeof(double));
                //dy=(double*) calloc(N,sizeof(double));
                //dz=(double*) calloc(N,sizeof(double));
                xt=(double*) calloc(N,sizeof(double));
                //yt=(double*) calloc(N,sizeof(double));
                //zt=(double*) calloc(N,sizeof(double));
                dx1=(double*) calloc(N,sizeof(double));
                //dy1=(double*) calloc(N,sizeof(double));
                //dz1=(double*) calloc(N,sizeof(double));
                dx2=(double*) calloc(N,sizeof(double));
                //dy2=(double*) calloc(N,sizeof(double));
                //dz2=(double*) calloc(N,sizeof(double));
                                
                vx=gasdev(&idum);
                //vy=gasdev(&idum);
                //vz=gasdev(&idum);//?
                sum=pow(vx,2);//+pow(vy,2)+pow(vz,2);
                sum=sqrt(sum);
                vx=v*vx/sum;
                //vy=v*vy/sum;
                //vz=v*vz/sum;//normalize v

                n=0;
                sprintf(NAME1,"./traj_data/dump100/dump100_N%dtaur%.4fv%.4lfl%ddt%.4fF%.4f_%d.xyz",N,1.0/gammar,v,legs,dt,F,j);
                printf("%s\n",NAME1);
                ofp=fopen(NAME1,"r");
                ofp=NULL;
                if(ofp==NULL){
                        x[0]=0.0;       
                        //y[0]=0.0;       
                        //z[0]=0.0;       
                        for(n=0;n<legs/2;n++){
                                //di_x=gasdev(&idum);
                                //di_y=gasdev(&idum);
                                //di_z=gasdev(&idum);//?
                                di_x=1.0;
                                //sum=pow(di_x,2)+pow(di_y,2)+pow(di_z,2);
                                //sum=sqrt(sum);
                                //di_x=di_x/sum;
                                //di_y=di_y/sum;
                                //di_z=di_z/sum;//normalize direction
                                for(m=1;m<ptl_num+1;m++){
                                        x[n*2*ptl_num + m] = di_x * m; 
                                        //y[n*2*ptl_num + m] = di_y * m; 
                                        //z[n*2*ptl_num + m] = di_z * m; 
                                }
                                for(m=1;m<ptl_num+1;m++){
                                        x[(n*2+1)*ptl_num + m] = -di_x * m; 
                                        //y[(n*2+1)*ptl_num + m] = -di_y * m; 
                                        //z[(n*2+1)*ptl_num + m] = -di_z * m; 
                                }
                        } 
                        plustime = 100000/dt;
                } 
                else{
                        for(n=0;;n++){                                                        
                                checkeof = fscanf(ofp,"%d\t%lf\n",&ts,&tx);
                                if(checkeof==EOF){
                                        break;
                                }
                                vx = tx;
                                for(m=0;m<N;m++){
                                        checkeof = fscanf(ofp,"%d\t%lf\n",&ts,&tx);
                                        if(checkeof==EOF){
                                                break;
                                        }
                                        x[m] = tx;
                                }
                                if(checkeof==EOF){
                                        break;
                                }
                        }//read 
                        length=n;
                        fclose(ofp);//load coordinates
                }

                sprintf(NAME4,"./traj_data/dump100/taetaxdot_N%dtaur%.4fv%.4lfl%ddt%.4fF%.4f_%d.xyz",N,1.0/gammar,v,legs,dt,F,j);
                ofpeta=fopen(NAME4,"w");
                sprintf(NAME2,"./traj_data/dump100/dump100_N%dtaur%.4fv%.4lfl%ddt%.4fF%.4f_%d.xyz",N,1.0/gammar,v,legs,dt,F,j);
                ofp2=fopen(NAME2,"w");
                sprintf(NAME3,"./traj_data/dump1/dump1_N%dtaur%.4fv%.4lfl%ddt%.4fF%.4f_%d.xyz",N,1.0/gammar,v,legs,dt,F,j);
                ofp3=fopen(NAME3,"w");
                sprintf(NAME4,"./traj_data/dump0.01/dump0.01_N%dtaur%.4fv%.4lfl%ddt%.4fF%.4f_%d.xyz",N,1.0/gammar,v,legs,dt,F,j);
                ofp4=fopen(NAME4,"w");
                                
                                
                for(t=1;t<T2+plustime+1;t++){
                        for(i=0;i<N;i++){//dx=(xi/gamma+v)*dt
                                dx1[i]=xi*gasdev(&idum)*dt;
                                //dy1[i]=xi*gasdev(&idum)*dt;
                                //dz1[i]=xi*gasdev(&idum)*dt;
                                if(s[i]==1){
                                        dx1[i]=dx1[i]+vx*dt;
                                        //dy1[i]=dy1[i]+vy*dt;
                                        //dz1[i]=dz1[i]+vz*dt;
                                }
                        }//random and active force
                        vx=vx-gammar*vx*dt+eta*gasdev(&idum)*dt;
                        //vy=vy-gammar*vy*dt+eta*gasdev(&idum)*dt;
                        //vz=vz-gammar*vz*dt+eta*gasdev(&idum)*dt;
                        for(n=0;n<legs;n++){
                                for(m=n*ptl_num;m<(n+1)*ptl_num;m++){
                                        i=m;
                                        if(i%ptl_num==0){
                                                i=0;
                                        }
                                        rx=x[m+1]-x[i];
                                        //ry=y[m+1]-y[i];
                                        //rz=z[m+1]-z[i];

                                        dx1[i]=dx1[i]+rx*k*dt;
                                        //dy1[i]=dy1[i]+ry*k*dt;
                                        //dz1[i]=dz1[i]+rz*k*dt;
                                        dx1[m+1]=dx1[m+1]-rx*k*dt;
                                        //dy1[m+1]=dy1[m+1]-ry*k*dt;
                                        //dz1[m+1]=dz1[m+1]-rz*k*dt;
                                }
                        }//force from bond energy
                        for(m=0;m<N;m++){
                                if(s[m]==1){
                                      dx1[m]=dx1[m] + F*dt; 
                                }
                        }
                        for(i=0;i<N;i++){
                                x[i]=x[i]+dx1[i];
                                //y[i]=y[i]+dy1[i];
                                //z[i]=z[i]+dz1[i];
                        }//position update
                        if(t>plustime && t%dump2==0){
                                fprintf(ofp2,"%d\t%lf\n",2,vx);
                                for(i=0;i<N;i++){
                                       fprintf(ofp2,"%d\t%lf\n",s[i],x[i]);
                                }
                                fprintf(ofpeta,"%lf\t%lf\n",vx,(dx1[0])/dt);
                        }
                        if(t>plustime && t<T3+1+plustime && t%dump3==0){
                                fprintf(ofp3,"%d\t%lf\n",2,vx);
                                for(i=0;i<N;i++){
                                        fprintf(ofp3,"%d\t%lf\n",s[i],x[i]);
                                }
                        }
                        if(t>plustime && t<T4+1+plustime && t%dump4==0){
                                fprintf(ofp4,"%d\t%lf\n",2,vx);
                                for(i=0;i<N;i++){
                                        fprintf(ofp4,"%d\t%lf\n",s[i],x[i]);
                                }
                        }
                }
                //fclose(ofp1);
                fclose(ofp2);
                fclose(ofp4);
                fclose(ofp3);
                fclose(ofpeta);
                free(s);
                free(x);
                //free(y);
                //free(z);
                free(dx);
                //free(dy);
                //free(dz);
                free(xt);
                //free(yt);
                //free(zt);
                free(dx1);
                //free(dy1);
                //free(dz1);
                free(dx2);
                //free(dy2);
                //free(dz2);
        }
        MPI_Finalize();
        return 0;
}

double ran2(long *idum)
{
        int j;
        long k;
        static long idum2=123456789;
        static long iy=0;
        static long iv[NTAB];
        float temp;
        if (*idum <= 0) { //Initialize.
                if (-(*idum) < 1) *idum=1; //Be sure to prevent idum = 0.
                else *idum = -(*idum);
                idum2=(*idum);
                for (j=NTAB+7;j>=0;j--) { //Load the shuffle table (after 8 warm-ups).
                        k=(*idum)/IQ1;
                        *idum=IA1*(*idum-k*IQ1)-k*IR1;
                        if (*idum < 0) *idum += IM1;
                        if (j < NTAB) iv[j] = *idum;
                }
                iy=iv[0];
        }
        k=(*idum)/IQ1; //Start here when not initializing.
        *idum=IA1*(*idum-k*IQ1)-k*IR1;// Compute idum=(IA1*idum) % IM1 without
        if (*idum < 0) *idum += IM1;// overflows by Schrage’s method.
        k=idum2/IQ2;
        idum2=IA2*(idum2-k*IQ2)-k*IR2;// Compute idum2=(IA2*idum) % IM2 likewise.
        if (idum2 < 0) idum2 += IM2;
        j=iy/NDIV; //Will be in the range 0..NTAB-1.
        iy=iv[j]-idum2; //Here idum is shuffled, idum and idum2 are
        iv[j] = *idum; //combined to generate output.
        if (iy < 1) iy += IMM1;
        if ((temp=AM*iy) > RNMX) return RNMX;// Because users don’t expect endpoint values.
        else return temp;
}
double gasdev(long *idum)
{
        double ran2(long *idum);
        static int iset=0;
        static float gset;
        float fac,rsq,v1,v2;
        if (*idum < 0) iset=0; //Reinitialize.
        if (iset == 0) { //We don’t have an extra deviate handy, so
                do {
                        v1=2.0*ran2(idum)-1.0; //pick two uniform numbers in the square ex
                        v2=2.0*ran2(idum)-1.0; //tending from -1 to +1 in each direction,
                        rsq=v1*v1+v2*v2; //see if they are in the unit circle,
                } while (rsq >= 1.0 || rsq == 0.0); //and if they are not, try again.
                fac=sqrt(-2.0*log(rsq)/rsq);
                gset=v1*fac;
                iset=1; //Set flag.
                return v2*fac;
        } else { //We have an extra deviate handy,
                iset=0; //so unset the flag,
                return gset;// and return it.
        }
}

