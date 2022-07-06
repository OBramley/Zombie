// #define PY_SSIZE_T_CLEAN
// #include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// long double *spattospin1(long double *H1ea,int norb){
//     //Converting H1ea from spatial to spin
//     int i, j, ii, jj, nspao;
//     long double *H1ei = calloc(norb*norb,sizeof(long double));

//     if(norb%2 !=0){
//         printf("norb must be even");
//         // Need to put in error catching
//     }
//     nspao = norb/2;
    
//     for(i=0;i<nspao;i++){
//         for(j=0;j<nspao;j++){
//             ii=i*2;
//             jj=j*2;
//             *(H1ei+ii+jj) = *(H1ea+(i*nspao)+j);
//             // *(H1ei+ii+jj+1) = 0.0;
//             *(H1ei+ii+jj+norb+1) = *(H1ea+(i*nspao)+j);
//             // *(H1ei+ii+jj+norb) = 0.0;
//         }
//     }
//     return H1ei;
// }

// long double *spattospin2(long double *H1ea, int norb){
//     int i, j, k, l, ii, jj, kk, ll, nspao;
//     long double *H2ei = calloc(norb*norb*norb*norb,sizeof(long double));

//        if(norb%2 !=0){
//         printf("norb must be even");
//         // Need to put in error catching
//     }
//     nspao = norb/2;
//     for(i=0;i<nspao;i++){
//         for(j=0;j<nspao;j++){
//             ii=i*2;
//             jj=i*2;
//             for(k=0;k<nspao;k++){
//                 for(l=0;l<nspao;l++){
//                     kk=k*2;
//                     ll=l*2;
//                     // Need to work out how the memory has been allocated

//                 }
//             }
//         }
//     }
//     return H2ei;
// }

long double overlap(int norb,long double z1[norb][2],long double z2[norb][2]){
    long double tt, temp = 1.0;
    int i;
    for(i=0;i<norb;i++){
        tt = (z1[i][0]*z2[i][0])+(z1[i][1]*z2[i][1]);
        temp=temp*tt;
    }
    return temp;
}

long double z_an_z3(int norb,long double z1ij[norb][norb][norb][2],long double z2k[norb][norb][2],long double H2ei[norb][norb][norb][norb],int i, int j, int k){
    long double an, vmult[norb][2], gg[norb], hh[norb];
    int n, gmax, hmin;

    for(n=0;n<norb;n++){
        vmult[n][0]=z1ij[i][j][n][0]*z2k[k][n][0];
        vmult[n][1]=z1ij[i][j][n][1]*z2k[k][n][1];
        gg[n]=0;
        hh[n]=0;
    }
    gg[0]=vmult[0][0]-vmult[0][1];
    gmax=norb;
    for(n=1;n<norb;n++){
        gg[n]=gg[n-1]*(vmult[n][0]-vmult[n][1]);
        if(gg[n]==0.0){
            gmax=n;
            break;
        }
    }
    hmin=0;
    hh[norb-1]=vmult[norb-1][0]+vmult[norb-1][1];
    for(n=norb-2;n>0;n--){
        hh[n]=hh[n+1]*(vmult[n][0]+vmult[n][1]);
        if(hh[n]==0.0){
            hmin=n;
            break;
        }
    }
    
    an=0.0;
    if(gmax<hmin){
        return 0.0;
    }
    if(H2ei[i][j][k][0]!=0.0){
        an=an+(z1ij[i][j][0][0]*z2k[k][0][1]*hh[1]*H2ei[i][j][k][0]);
    }
    for(n=1;n<norb-1;n++){
        if(H2ei[i][j][k][1]!=0.0){
            an=an+(gg[n-1]*z1ij[i][j][n][0]*z2k[k][n][1]*hh[n+1]*H2ei[i][j][k][n]);
        }
    }
    if(H2ei[i][j][k][norb-1]!=0.0){
        an = an+(gg[norb-2]*z1ij[i][j][norb-1][0]*z2k[k][norb-1][1]*H2ei[i][j][k][norb-1]);
    }
    return an;
}

int iszero(int norb,long double z1ij[norb][norb][norb][2],int i, int j){
    long double tt;
    int n;

    for(n=0;n<norb;n++){
        tt = (z1ij[i][j][n][0]*z1ij[i][j][n][0])+(z1ij[i][j][n][1]*z1ij[i][j][n][1]);
        if(tt==0.0){
            return 0;
        }
    }

    return 1;


}

long double Ham1z(int norb,long double H1ei[norb][norb],long double zom1[norb][2], long double zom2[norb][2]){
    long double zomt[norb][2], Ht1=0.0;
    int i,j,n;
    
    for(i=0;j<norb;i++){
        for(j=0;j<norb;j++){
            memcpy(zomt,zom2,2*norb*sizeof(long double));
            //anhilation Operator
            zomt[i][0]=zomt[i][1];
            zomt[i][1]=0.0;
            for(n=0;n<i;n++){
                zomt[n][1]=zomt[n][1]*-1;
            }
            //creation Operator
            zomt[j][1]=zomt[j][0];
            zomt[j][0]=0.0;
            for(n=0;n<j;n++){
                zomt[n][1]=zomt[n][1]*-1;
            }
            
            long double ov = overlap(norb,zom1,zom2);
            Ht1 = ov*(H1ei[i][j]);
        }
    }
    return Ht1;
}

long double Ham2z(int norb, long double H2ei[norb][norb][norb][norb], long double zom1[norb][2], long double zom2[norb][2]){
    long double Ht2 = 0.0;
    long double z1ij[norb][norb][norb][2], z2k[norb][norb][2], zomt[norb][2];
    int i, j, k, l, n;

    for(i=0; i<norb; i++){
        for(j=0;j<norb;j++){
            memcpy(zomt,zom1,2*norb*sizeof(long double));
            //anhilation Operator
            zomt[i][0]=zomt[i][1];
            zomt[i][1]=0.0;
            for(n=0;n<i;n++){
                zomt[n][1]=zomt[n][1]*-1;
            }
            //creation Operator
            zomt[j][1]=zomt[j][0];
            zomt[j][0]=0.0;
            for(n=0;n<j;n++){
                zomt[n][1]=zomt[n][1]*-1;
            }
            for(k=0;k<norb;k++){
                for(l=0; k<2;l++){
                    z1ij[i][j][k][l]=zomt[k][l];
                }
            }
        }
    }

    for(k=0;k<norb;k++){
        memcpy(zomt,zom2,2*norb*sizeof(long double));
        //anhilation Operator
        zomt[i][0]=zomt[i][1];
        zomt[i][1]=0.0;
        for(n=0;n<i;n++){
            zomt[n][1]=zomt[n][1]*-1;
        }

        for(i=0;i<norb;i++){
                for(j=0; j<2;j++){
                    z2k[k][i][j]=zomt[i][j];
                }
            }
    }

    for(i=0; i<norb;i++){
        if(zom1[i][1]==0.0){
            continue;
        }
        int ispin = i%2;
        for(j=0;j<norb;j++){
            if(iszero(norb,z1ij,i,j)==0){
                continue;
            }
            // int jspin = j%2;
            for(k=ispin;k<norb;k+=2){
                if(zom2[k][1]==0.0){
                    continue;
                }
                Ht2=Ht2+z_an_z3(norb,z1ij,z2k,H2ei,i,j,k);
            }
        }
    }
    return 0.5*Ht2;
}

long double HTot(long double Hnr, int norb,long double H1ei[norb][norb],long double H2ei[norb][norb][norb][norb],long double zom1[norb][2], long double zom2[norb][2]){
    long double H1et, H2et, HH;
    H1et=Ham1z(norb,H1ei,zom1,zom2);
    H2et=Ham2z(norb,H2ei,zom1,zom2);
    HH= H1et+H2et+(Hnr*overlap(norb,zom1,zom2));
    return HH;
}

int csvwriter(int ndet, long double input[ndet][ndet],char filename[]){
    int i, j;
    FILE *fp1;
    fp1 = fopen(filename, "w");
    if (fp1 == NULL)
    {
        printf("Error while opening the file.\n");
        return 0;
    }
    for(i=0; i<ndet;i++){
        for(j=0; j<ndet; j++){
            fprintf(fp1, "%Lf%s",input[i][j],
                (j<ndet-1?",":""));
        }
        fprintf(fp1,"\n");
    }
    fclose(fp1);
    return 0;
}

int hamiltonian(int ndet, int norb, long double H1ei[norb][norb],long double H2ei[norb][norb][norb][norb], long double Hnr,long double zstore[ndet][norb][2]){
    long double kover[ndet][ndet], bigham[ndet][ndet], zomi[norb][2], zomj[norb][2];
    int i, j;

    for(i=0; i<ndet; i++){
        for(j=i; j<ndet;j++){
            memcpy(zomi, zstore[i], 2*norb*sizeof(long double));
            memcpy(zomi, zstore[j], 2*norb*sizeof(long double));
            kover[i][j]=overlap(norb,zomi,zomj);
            kover[j][i]=kover[i][j];
            bigham[i][j]=HTot(Hnr,norb,H1ei,H2ei,zomi,zomj);
            bigham[j][i]=bigham[i][j];
        }
        printf("\nHamiltonian row %d complete",i);
    }

    csvwriter(ndet, kover, "overlap_2.csv");
    csvwriter(ndet, bigham, "bigham_2.csv");

    return 0;

}

int control(const long double * z_in, const long double * h1in, const long double *h2in, const long double Hnuc, const int ndet, const int norb){
    // int ndet=10;
    // int norb=10;
    long double zstore[ndet][norb][2], H1ei[norb][norb],H2ei[norb][norb][norb][norb];
    // long double temp; //, Hnuc=1.5;
    int i,j,k,l,n,m,p;
    // char *record,*line, *eptr;
    // char buffer[100000] ;

    n=0;
    for(i=0;i<ndet*norb*2;i+=ndet){
        m=0;
        for(j=0;j<norb*2;j+=2){
            zstore[n][m][0]=h1in[i+j];
            zstore[n][m][1]=h1in[i+j+1];
            m++;
        }
    }

    n=0,m=0;
    for(i=0;i<norb*norb;i+=norb){
        for(j=0; j<norb;j++){
            H1ei[n][j]=h1in[i+j];
        }
        n++;
    }

    n=0, m=0, p=0;
    for(i=0;i<norb*norb*norb*norb;i+=(norb*norb*norb)){
        m=0;
        for(j=0;j<norb*norb*norb;j+=(norb*norb)){
            p=0;
            for(k=0;k<norb*norb;k+=norb){
                for(l=0;l<norb;k++){
                    H2ei[n][m][p][l]=h2in[i+j+k+l];
                }
                p++;
            }
            m++;
        }
        n++;
    }

    hamiltonian(ndet,norb,H1ei,H2ei,Hnuc,zstore);


    // Open H1ei
//     FILE *fstream = fopen("H1ei.csv","r");
//     if(fstream == NULL){
//         printf("\n file opening failed ");
//         return -1 ;
//     }
//     i=0,j=0;
//     while((line=fgets(buffer,sizeof(buffer),fstream))!=NULL){
//         record = strtok(line,",");
     
//         while(record != NULL){
//             printf("%s",record);
//             sscanf(record,"lf",&temp);
//         //  temp = strtold(record,&eptr);  
//             printf("%lf",temp); 
//             H1ei[i][j] = temp;
//             record = strtok(NULL,",");
//             j++;
//         }
//         printf("\n");
//         ++i;
//    }

    // for(i=0;i<norb;i++){
    //     for(j=0;j<norb;j++){
    //         fscanf(fstream,"%lf",&temp);
    //         H1ei[i][j]=temp;
    //     }
    // }
    // fclose(fstream);

    // csvwriter(ndet,H1ei,"H1ei_2.csv");

    return 0;
    
}


    