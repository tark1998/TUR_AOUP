#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "mpi.h"

//int main(){
int main(int argc,char **argv){
        int rank, size;

        MPI_Status status;

        MPI_Init( &argc, &argv );
        MPI_Comm_rank( MPI_COMM_WORLD, &rank );
        MPI_Comm_size( MPI_COMM_WORLD, &size );

        int spr = 1/size;

        int i,j,n,m,p,N,former_size,size_,s,discard_front=0,legs,ptl_num=50,file_cnt=0;
        int eof_check;
        int sete,lag;
        double sum, sumx, sumy, sumz, dx2,k,a,gamma,v,q, scale, gammar;
        double dt,F;
        double *x,*y,*z,*tamsd,*tamsd2,*tamdx2,*tamdy2,*tamdz2;
        double t,dx,dy,dz,tx,ty,tz;
        double time[5000] = {0};
        FILE *data1, *data2, *ifp, *xyz;
        char extra[500];
        char NAME[500], NAME2[500], XYZ[500];
        char dump1[500], dump2[500];
        char OPENF[10000];
        char *file_base[3] = {"10000"};
        char file_base_[500];
        if(argc<1){
                printf("ptln a v legs k scale\n");
                exit(1);
        }
        ptl_num=atof(argv[1]);
        v=atof(argv[2]);
        legs=atoi(argv[3]);
        dt=atof(argv[4]);
        gammar=atof(argv[5]);
        F=atof(argv[6]);
        scale=1;
        //printf("%s\n",extra);
        a/=scale;
        ptl_num/=scale;
        v*=scale;

        N=ptl_num*legs+1;
        
        for (p=0; p<1; p++){
                discard_front = 0;
                size_=100/pow(scale,2);
                if(p==2){
                        discard_front = 0;
                        size_=3000;
                }
                sprintf(file_base_, file_base[p]);
                tamsd2=(double*) calloc(size_-1,sizeof(double));
                tamdx2=(double*) calloc(size_-1,sizeof(double));
                //tamdy2=(double*) calloc(size_-1,sizeof(double));
                //tamdz2=(double*) calloc(size_-1,sizeof(double));
                file_cnt = 0;
                for(i=0;i<100;i++){
                        eof_check = 1;

                        sprintf(OPENF,"./traj_data/dump%s/dump%s_N%dtaur%.4fv%.4lfl%ddt%.4fF%.4f_%d.xyz",file_base_,file_base_,N,1.0/gammar,v,legs,dt,F,i);
                        ifp=fopen(OPENF,"r");
                        if (ifp==NULL){
                               continue; 
                        }

                        x=(double*) calloc(size_,sizeof(double));
                        //y=(double*) calloc(size_,sizeof(double));
                        //z=(double*) calloc(size_,sizeof(double));
                        
                        
                        for(n=0;n<discard_front;n++){
                                for(j=0;j<N+1;j++){
                                        eof_check = fscanf(ifp,"%d\t%lf\n",&s,&tx);
                                        if (eof_check == EOF) break;
                                }
                                if (eof_check == EOF) break;
                        }//read 
                        if (eof_check == EOF) {
                                //free(z);
                                //free(y);
                                free(x);
                                fclose(ifp);
                                continue;
                        }
                        
                        for(n=0;n<size_;n++){
                                for(j=0;j<N+1;j++){
                                        eof_check = fscanf(ifp,"%d\t%lf\n",&s,&tx);
                                        if(s==1){
                                                x[n]=tx;
                                                //y[n]=ty;
                                                //z[n]=tz;
                                        }
                                        if (eof_check == EOF) break;
                                }
                                if (eof_check == EOF) break;
                        }//read 
                        if (eof_check == EOF){
                                //free(z);
                                //free(y);
                                free(x);
                                fclose(ifp);
                                continue;
                        }
                        
                        file_cnt += 1;

                        for(lag=1;lag<size_;lag++){
                                sum=0;
                                sumx=0;
                                //sumy=0;
                                //sumz=0;
                                for(n=0;n<size_-lag;n++){
                                        dx=x[n+lag]-x[n];
                                        dx2 = pow(dx,2);
                                        //dx2 = pow(x[n+lag]-x[n],2)+pow(y[n+lag]-y[n],2)+pow(z[n+lag]-z[n],2);
                                        sum+=dx2;
                                        sumx+=dx;
                                }
                                tamsd2[lag-1]+=sum/(double)(size_-lag);
                                tamdx2[lag-1]+=sumx/(double)(size_-lag);
                                //tamsd[lag-1]+=sum/(double)N;
                        }
                
                        //free(z);
                        //free(y);
                        free(x);
                }
                if (file_cnt == 0){
                        printf("No ensemble\n");
                        printf("%s\n",OPENF);
                        free(tamsd2);                        
                        free(tamdx2);                        
                        //free(tamdy2);                        
                        //free(tamdz2);                        
                        MPI_Finalize();
                        return 0;
                }
                printf("i: %d ensemble, dump: %s sec\n", file_cnt, file_base[p]);

                sprintf(NAME,"./msd_data/msx%s_middle_N%dtaur%.4fv%.4lfl%ddt%.4fF%.4f.dat",file_base[p],(int)(1+legs*ptl_num*scale),1.0/gammar,v/scale,legs,dt,F);
                data1 = fopen(NAME,"w");
                for(lag=1;lag<size_;lag++){
                        fprintf(data1,"%f\t%.12lf\n",lag * atof(file_base[p]) * pow(scale,2),(tamsd2[lag-1])/(double)(file_cnt) * pow(scale,2) - (pow(tamdx2[lag-1],2))/(double)pow(file_cnt,2)); 
                }
                fclose(data1);
                free(tamsd2);
        }
        
        MPI_Finalize();
        return 0;
}
