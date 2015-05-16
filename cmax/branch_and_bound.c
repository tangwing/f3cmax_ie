/** @file
 * @author Samia Kouki (Implémentation de l'algorithme Branch @& Bound en C)
 * @author Pierre-Antoine Morin (Adaptation et intégration du code dans ce programme C++)
 * @brief Fichier source définissant la fonction @ref branchAndBound()
 * @note Des commentaires ont été ajoutés au code original (style C++, en début de ligne).
 * @note Ils donnent des indications sur le fonctionnement du programme, et la manière dont le code a été ajusté.
 * @par Remarque à propos des indices
 * La convention choisie pour représenter les indices est différente :
 * @li Dans le code de ce fichier, tous les indices commencent à @b 1.
 * @li Dans le reste du programme, tous les indices commencent à @b 0.
 * @see Code de la fonction @ref branchAndBound() (décalage d'indice effectué pendant la lecture et l'écriture)
 */
//#include "branch_and_bound.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

// Added : optional macro definition (possibly defined in <stdlib.h>)
#ifndef min
#define min(x,y) (((x) < (y)) ? (x) : (y))
#endif /* min */
#ifndef max
#define max(x,y) (((x) < (y)) ? (y) : (x))
#endif /* max */

typedef struct lfs
{
    long j;
    long b_inf;
    long rm1b_b;
    long rm2b_b;
    long qm1b_b;
    long qm2b_b;
    struct lfs * suiv;
    struct nd * Desc;
} lf;
typedef struct Flfs
{
    long Fj;
    long Fb_inf;
    long Findj;
    long FTabCmaxMIN[21];
    long FTabCmaxMOUT[21];
    long FTabJINF[2501];
    long Toul;
    struct Flfs * Fsuiv;
    struct Fnd * FDesc;
} Flf;
typedef struct nd
{
    lf * fils;
    long branch;
    long nb_f;
    long t_sin;
} noeud;
typedef struct Fnd
{
    Flf * Ffils;
    long Fbranch;
    long Fnb_f;
    long Ft_sin;
} Fnoeud;
typedef struct v_job
{
    long r;
    long p;
    long q;
    long job;
} vj;
typedef struct file_int
{
    long j;
    struct file_int * suiv;
} f_int;
typedef struct v_JOB
{
    long r1;
    long r2;
    long q1;
    long q2;
    long alpha;
    long beta;
    long job;
    long p1;
    long p2;
} vJ;
typedef struct JJOB
{
    long j;
    long p;
} jjob;
typedef struct pairMach
{
    long m1;
    long m2;
    long BinfP;
    long BSP;
    long TTlj[51]; /* ???????????????????? */
} PM;
typedef struct OneMach
{
    long m;
    long Binf;
} OM;
#define MAX_MACH 21
#define MAX_JOB 501
#define MAX_PAIR ((MAX_MACH-1)*(MAX_MACH-2)/2)+1
#define IN 1
#define OUT 0
#define FIN 1
#define FOUT 0
lf * P=NULL, *PP=NULL;
long seq[MAX_JOB],Seq[MAX_JOB],SSeq[MAX_JOB],F2SSSeq[MAX_JOB],SSSeq[MAX_JOB],seq_op[MAX_JOB],JR_b[MAX_JOB],n,m,Bborn_sup,born_sup=100000000;
long BSS,BSjak,BSjoh,BestLB=0,BB=0,MINbinf,MINbinfcmax,BP=0,Error,si=0,BOOLF2,UBFinal;
double ERROR;
long nb_iter=0;
long nb_exp=0;
long nbnc=0;
long nbrfcrnd=0;
long Fnb_exp=0;
long Fnbnc=0;
long Fnbnc1=0;
long Fnb_iter=0;
long Nbrjf2;
long FBSS=100000000,Fborn_sup=100000000;
long Reste1;
long FT[MAX_JOB][MAX_MACH];
long FJR_b[MAX_JOB];
long Fseq[MAX_JOB];
long TabCmaxMIN[MAX_MACH],TabCmaxMOUT[MAX_MACH];
long Fseq_op[MAX_JOB];
long Trj[MAX_JOB][MAX_MACH];
long Tqj[MAX_JOB][MAX_MACH];
long Trjqj[MAX_JOB][MAX_MACH];
PM TPM[MAX_PAIR];
OM FTPM[MAX_MACH];
long FbranchF2=0;
long TOUL=0;
long M1,M2,tdisp=0;
long JobR,FTINF,TailIN,TailOUT,Niv1,Niv2,Niv3;
long lj[MAX_JOB];
long Tlj[MAX_JOB];
long llj[MAX_JOB];
vj M1_j_bOne[MAX_JOB];
long SS1[MAX_JOB];
long F2SS1[MAX_JOB];
vJ SEQ[MAX_JOB];
vJ F2SEQ[MAX_JOB];
jjob S1[MAX_JOB],S2[MAX_JOB];
jjob F2S1[MAX_JOB],F2S2[MAX_JOB];
clock_t Sstart,Eend,EEend,TimeBB,FFin,Fdeb,realise;
double TEXE,FTEXE;
FILE *FR;
long double mc=2147483647.0;
long double a=16807.0,b=127773.0,c=2836.0;
long double time_seed;
long born_inf11(long Nb_Job_Restant,long rm1bb,long rm2bb,long qm1bb,long qm2bb,long job);
long born_inf2(long Nb_Job_Restant,long rm1bb,long rm2bb,long qm1bb,long qm2bb,long job,long Tin);
long cmaxfrql(long sequ[MAX_JOB],long R1,long R2,long Q1,long Q2);
long cmaxfrql1(long sequ[MAX_JOB],long R1,long R2,long Q1,long Q2);
long F2rjljqjCmax(long FNb_Job_Restant,long Fbranch,long FTabCmaxMIN[MAX_MACH],long FTabCmaxMOUT[MAX_MACH],long FTin,long FTout,long job);
long F2rjljqjCmax1(long FNb_Job_Restant,long Fbranch,long FTabCmaxMIN[MAX_MACH],long FTabCmaxMOUT[MAX_MACH],long FTin,long FTout,long job);
long F2rjljqjCmax2(long FNb_Job_Restant,long Fbranch,long FTabCmaxMIN[MAX_MACH],long FTabCmaxMOUT[MAX_MACH],long FTin,long FTout,long job);
long Makespan(long SEq[MAX_JOB],long tail,long m);

lf * placer (lf * p , lf * q)
{
    lf * pp;
    pp=p;
    if (p==NULL)
        p=q;
    else if (q->b_inf<p->b_inf)
    {
        q->suiv=p;
        p=q;
    }
    else
    {
        while (pp->suiv!=NULL)
        {
            if (pp->suiv->b_inf<q->b_inf)
                pp=pp->suiv;
            else
            {
                q->suiv=pp->suiv;
                break;
            }
        }
        pp->suiv=q;
    }
    return p;
}
Flf * Fplacer (Flf * p , Flf * q)
{
    Flf * pp;
    pp=p;
    if(p==NULL)
        p=q;
    else if(q->Fb_inf<p->Fb_inf)
    {
        q->Fsuiv=p;
        p=q;
    }
    else
    {
        while(pp->Fsuiv!=NULL)
        {
            if(pp->Fsuiv->Fb_inf<q->Fb_inf)
                pp=pp->Fsuiv;
            else
            {
                q->Fsuiv=pp->Fsuiv;
                break;
            }
        }
        pp->Fsuiv=q;
    }
    return p;
}
noeud * create_noeud (long nbr_f,long new_job,long t_sin,long branch,long rm1b,long rm2b,long qm1b, long qm2b,long binf_pere,long job)
{
    long i,k,rm1bb,rm2bb,qm1bb,qm2bb=0,t_sout,Tin,Tout,kkk,h,hh,maxk,binf;
    long l,Nb_Job_Restant;
    lf * pp;
    noeud * p;
    p =(noeud*)malloc(sizeof(noeud));
    p->fils = NULL;
    p->nb_f=nbr_f;
    p->branch =branch;
    t_sout=t_sin;
    ++nbrfcrnd;
    if (branch==IN)
        seq[++t_sin]=new_job;
    else
        seq[Nbrjf2-t_sin]=new_job;
    p->t_sin=t_sin;
    l=0;
    for (i=1; i<= Nbrjf2; i++)
    {
        for (k=1; k<=t_sout; k++)
            if (FJR_b[i]==seq[k]||FJR_b[i]==seq[Nbrjf2-k]) break;
        if (k > t_sout && FJR_b[i] != new_job && FJR_b[i] !=job)
            JR_b[++l]=FJR_b[i];
    }
    Nb_Job_Restant=Nbrjf2-t_sin-t_sout-1;
    Tin =t_sin;
    Tout =Tin;
    if (branch==OUT)
        Tin =t_sin+1;
    else
        Tout=Tin;
    memcpy(Seq,seq,sizeof(long)*(n+1));
    memcpy(SSeq,seq,sizeof(long)*(n+1));
    for(h=1; h<=Nb_Job_Restant; h++)
    {
        maxk=h;
        for(hh=h+1; hh<=Nb_Job_Restant; hh++)
            if(Trj[JR_b[hh]][M1] < Trj[JR_b[maxk]][M1]) maxk=hh;
        kkk=JR_b[h];
        JR_b[h]=JR_b[maxk];
        JR_b[maxk]=kkk;
    }
    for (i=1; i<=Nb_Job_Restant; i++)
    {
        if(branch==OUT)
        {
            rm1bb=max(rm1b,Trj[JR_b[i]][M1])+FT[JR_b[i]][M1];
            rm2bb=max(rm2b,rm1bb+lj[JR_b[i]])+FT[JR_b[i]][M2];
            qm1bb=qm1b;
            qm2bb=qm2b;
            Seq[Tin]=JR_b[i];
            SSeq[Tin]=JR_b[i];
        }
        else
        {
            qm2bb=max(qm2b,Tqj[JR_b[i]][M2])+FT[JR_b[i]][M2];
            qm1bb=max(qm1b,qm2bb+lj[JR_b[i]])+FT[JR_b[i]][M1];
            rm1bb=rm1b;
            rm2bb=rm2b;
            Seq[Nbrjf2-Tin]=JR_b[i];
            SSeq[Nbrjf2-Tin]=JR_b[i];
        }
        if(p->nb_f!=1)
        {
            binf=born_inf2(Nb_Job_Restant,rm1bb,rm2bb,qm1bb,qm2bb,i,Tin);
            nb_exp++;
        }
        else
        {
            seq[t_sin+1]=JR_b[i];
            binf = cmaxfrql1(seq,rm1bb,rm2bb,qm1bb,qm2bb);
        }
        BSS=BSjoh;
        if (BSS< born_sup )
        {
            born_sup=BSS;
            Bborn_sup=BSS;
            if (binf_pere==born_sup) break;
        }
        if(binf<born_sup)
        {
            pp=(lf*)malloc(sizeof(lf));
            nbnc++;
            pp->j=JR_b[i];
            pp->b_inf=max(binf,binf_pere);
            pp->rm1b_b=rm1bb;
            pp->rm2b_b=rm2bb;
            pp->qm1b_b=qm1bb;
            pp->qm2b_b=qm2bb;
            pp->suiv=NULL;
            pp->Desc=NULL;
            p->fils= placer(p->fils,pp);
        }

    }
    return p;
}
Fnoeud * Fcreate_noeud (long Fnbr_f,long Fnew_job,long Ft_sin,long Fbranch,long FTabCmaxMIN[MAX_MACH],long FTabCmaxMOUT[MAX_MACH],long Fbinf_pere)
{
    long i,k,Ft_sout,FTin,FTout,Fbinf,moujoud=0,TBFJ;
    long l,FNb_Job_Restant, FFTabCmaxMIN[MAX_MACH],FFTabCmaxMOUT[MAX_MACH];
    Flf * Fpp;
    Fnoeud * Fp;
    Fp =(Fnoeud*)malloc(sizeof(Fnoeud));
    Fp->Ffils = NULL;
    Fp->Fnb_f=Fnbr_f;
    Fp->Fbranch =Fbranch;
    Ft_sout=Ft_sin;
    for(i=1; i<=m; i++)
    {
        FFTabCmaxMIN[i]=0;
        FFTabCmaxMOUT[i]=0;
    }
    if (Fbranch==FIN)
        Fseq[++Ft_sin]=Fnew_job;
    else
        Fseq[n-Ft_sin+1]=Fnew_job;
    Fp->Ft_sin=Ft_sin;
    l=0;
    for (i=1; i<= n; i++)
    {
        for (k=1; k<=Ft_sout; k++)
            if (i==Fseq[k]||i==Fseq[n-k+1]) break;
        if (k > Ft_sout && i!= Fnew_job)
            FJR_b[++l]=i;
    }
    TBFJ=l;
    FNb_Job_Restant=n-Ft_sin-Ft_sout;
    FTin =Ft_sin;
    FTout =FTin;
    if(Fbranch==FOUT)
        FTin =Ft_sin+1;
    else
        FTout=FTin;
    TailIN=FTin;
    TailOUT=FTout;
    for (i=1; i<=FNb_Job_Restant; i++)
    {
        if(Fbranch==FOUT)
        {
            FFTabCmaxMIN[1]=FTabCmaxMIN[1]+FT[FJR_b[i]][1];
            for(l=2; l<=m; l++)
                FFTabCmaxMIN[l]=max(FFTabCmaxMIN[l-1],FTabCmaxMIN[l])+FT[FJR_b[i]][l];
            memcpy(FFTabCmaxMOUT,FTabCmaxMOUT,sizeof(long)*(m+1));
        }
        else
        {
            FFTabCmaxMOUT[m]=FTabCmaxMOUT[m]+FT[FJR_b[i]][m];
            for(l=m-1; l>=1; l--)
                FFTabCmaxMOUT[l]=max(FFTabCmaxMOUT[l+1],FTabCmaxMOUT[l])+FT[FJR_b[i]][l];
            memcpy(FFTabCmaxMIN,FTabCmaxMIN,sizeof(long)*(m+1));
        }
        if(Fp->Fnb_f!=1)
        {
            FTINF=FTin+1;

// The following portion of code links the level in the tree with the lower bound to use.
// It is replaced with the systematic use of the same lower bound (constraints relaxed on all machines but one).
// This replacement is made, because the variable 'nbpair' is not initialized when 'm==3' (see below).
            /*
            if((TailIN+TailOUT)<=Niv1)
                Fbinf=F2rjljqjCmax1(FNb_Job_Restant,Fbranch,FFTabCmaxMIN,FFTabCmaxMOUT,FTin,FTout,i);
            if(((TailIN+TailOUT)>Niv1) &&((TailIN+TailOUT)<=Niv2))
                Fbinf=F2rjljqjCmax2(FNb_Job_Restant,Fbranch,FFTabCmaxMIN,FFTabCmaxMOUT,FTin,FTout,i);
            if(((TailIN+TailOUT)>Niv2) &&((TailIN+TailOUT)<=n))
                Fbinf=F2rjljqjCmax(FNb_Job_Restant,Fbranch,FFTabCmaxMIN,FFTabCmaxMOUT,FTin,FTout,i);
            */
            Fbinf=F2rjljqjCmax1(FNb_Job_Restant,Fbranch,FFTabCmaxMIN,FFTabCmaxMOUT,FTin,FTout,i);

            Fnb_exp++;
        }
        else
        {
            Fseq[Ft_sin+1]=FJR_b[i];
            Fbinf = Makespan(Fseq,n,m);
            realise=1;
        }

        if (Fbinf_pere==Fborn_sup) break;
        if(Fbinf<Fborn_sup)
        {
            Fpp=(Flf*)malloc(sizeof(Flf));
            Fnbnc++;
            Fpp->Fj=FJR_b[i];
            if(Fp->Fnb_f!=1)
                Fpp->Fb_inf=max(Fbinf,Fbinf_pere);
            else
                Fpp->Fb_inf=Fbinf;
            memcpy(Fpp->FTabCmaxMIN,FFTabCmaxMIN,sizeof(long)*(m+1));
            memcpy(Fpp->FTabCmaxMOUT,FFTabCmaxMOUT,sizeof(long)*(m+1));
            Fpp->Fsuiv=NULL;
            Fpp->FDesc=NULL;
            Fp->Ffils= Fplacer(Fp->Ffils,Fpp);
        }
    }

    return Fp;
}

void explorer(noeud * N)
{
    lf * p, * q;
    long t_sin;
    double Eroor;
    nb_iter++;
    p= N->fils;
    if (p!=NULL)
    {
        MINbinf=p->b_inf;
        Eroor=(((born_sup-p->b_inf)*100)/p->b_inf);
    }
    t_sin=N->t_sin;
    Eend=clock();
    TEXE=(double)((Eend-Sstart)/CLOCKS_PER_SEC);
    if((TEXE>=2)||(Eroor<=0.05))
    {
        born_sup=-1;
        BOOLF2=1;
    }

    while (p!=NULL)
    {
        if (p->b_inf<born_sup)
        {
            MINbinf=1000000000;
            if (N->nb_f==1)
            {
                seq[t_sin+1]=p->j;
                memcpy(seq_op,seq,sizeof(long)*(n+1));
                born_sup=p->b_inf;
                Bborn_sup=p->b_inf;
            }
            else
            {
                p->Desc=create_noeud(N->nb_f-1,p->j,t_sin,1-N->branch,p->rm1b_b,p->rm2b_b,p->qm1b_b,p->qm2b_b,p->b_inf,JobR);
                explorer(p->Desc);
                free(p->Desc);
                q=p;
                p=p->suiv;
                if(p!=NULL)
                    MINbinf=min(MINbinf,p->b_inf);
                free(q);
            }
        }
        else
        {
            while (p!=NULL)
            {
                q=p;
                p=p->suiv;
                if(p!=NULL)
                    MINbinf=min(MINbinf,p->b_inf);
                free(q);
            }
            break;
        }
    }
}
void Fexplorer(Fnoeud * FN)
{
    Flf * Fp, * Fq;
    long Ft_sin;
    double UB;
    Fnb_iter++;
    Fp= FN->Ffils;
    Ft_sin=FN->Ft_sin;
    if(Fp!=NULL)
    {
        UB=Fborn_sup;
        MINbinfcmax=Fp->Fb_inf;
    }
    FFin=clock();
    FTEXE=(double)((FFin-Fdeb)/CLOCKS_PER_SEC);
    if(FTEXE>=18000)
        Fborn_sup=-1;
    while (Fp!=NULL)
    {
        if (Fp->Fb_inf<Fborn_sup)
        {
            MINbinfcmax=100000000;
            if (FN->Fnb_f==1)
            {
                Fseq[Ft_sin+1]=Fp->Fj;
                memcpy(Fseq_op,Fseq,sizeof(long)*(n+1));
                Fborn_sup=Fp->Fb_inf;
                UBFinal=Fborn_sup;
            }
            else
            {
                Fp->FDesc=Fcreate_noeud(FN->Fnb_f-1,Fp->Fj,Ft_sin,1-FN->Fbranch,Fp->FTabCmaxMIN,Fp->FTabCmaxMOUT,Fp->Fb_inf);
                Fexplorer(Fp->FDesc);
                free(Fp->FDesc);
                Fq=Fp;
                Fp=Fp->Fsuiv;
                if(Fp!=NULL)
                    MINbinfcmax=min(MINbinfcmax,Fp->Fb_inf);
                free(Fq);
            }
        }
        else
        {
            while (Fp!=NULL)
            {
                Fq=Fp;
                Fp=Fp->Fsuiv;
                if(Fp!=NULL)
                    MINbinfcmax=min(MINbinfcmax,Fp->Fb_inf);
                free(Fq);
            }
            break;
        }
    }
}
f_int* placer2 (vj *M_j_b,int i,f_int *Head)
{
    f_int *Ph,*Phead;
    Ph=(f_int*)malloc(sizeof(f_int));
    Ph->j=i;
    if((Head==NULL) ||(M_j_b[i].q>M_j_b[Head->j].q))
    {
        Ph->suiv=Head;
        Head=Ph;
    }
    else
    {
        Phead=Head;
        while ( Phead->suiv!=NULL)
        {
            if(M_j_b[i].q<=M_j_b[Phead->suiv->j].q)
                Phead=Phead->suiv;
            else
                break;
        }
        Ph->suiv=Phead->suiv;
        Phead->suiv=Ph;
    }
    return Head;
}
long jackson(vj *M_j_b, long size)
{
    f_int * Th,*Head=NULL;
    long i,LB=0,J_c=1;
    tdisp=0;
    tdisp = M_j_b[1].r+M_j_b[1].p;
    i=2;
    while(1)
    {
        while( (i<= size) && (M_j_b[i].r< tdisp) )
        {
            if( M_j_b[i].q<=M_j_b[J_c].q)
            {
                Head=placer2(M_j_b,i,Head);
                i++;
            }
            else
            {
                M_j_b[J_c].p = tdisp-M_j_b[i].r;
                Head=placer2(M_j_b,J_c,Head);
                tdisp=tdisp - M_j_b[J_c].p+M_j_b[i].p;
                J_c=i++;
            }
        }
        LB=max(LB,tdisp + M_j_b[J_c].q);
        if (Head!=NULL)
        {
            if (i<= size &&(M_j_b[i].r== tdisp) && (M_j_b[i].q>M_j_b[Head->j].q))
                J_c=i++;
            else
            {
                Th=Head;
                J_c= Head->j;
                Head= Head->suiv;
                free(Th);
            }
            tdisp=tdisp+M_j_b[J_c].p;
        }
        else
        {
            if (i<=size)
            {
                J_c=i++;
                tdisp=M_j_b[J_c].r+M_j_b[J_c].p;
            }
            else return LB;
        }
    }
}
long cmaxfrql(long sequ[MAX_JOB],long R1,long R2,long Q1,long Q2)
{
    long C1=0,C2=0,Cmax=0,i;
    C1=R1;
    C2=R2;
    C1=max(C1,Trj[sequ[1]][M1])+FT[sequ[1]][M1];
    C2= max(C2,C1+Tlj[sequ[1]])+FT[sequ[1]][M2];
    Cmax= C2 + Tqj[sequ[1]][M2];
    for(i=2; i<= Reste1; i++)
    {
        C1 = max(C1,Trj[sequ[i]][M1])+FT[sequ[i]][M1];
        C2 = max(C1+Tlj[sequ[i]],C2)+FT[sequ[i]][M2];
        Cmax = max(Cmax,C2+Tqj[sequ[i]][M2]);
    }
    Cmax=max(Cmax,max((C1+Q1),(C2+Q2)));
    return Cmax;
}
long cmaxfrql1(long sequ[MAX_JOB],long R1,long R2,long Q1,long Q2)
{
    long C1=0,C2=0,Cmax=0,i;
    C1=R1;
    C2=R2;
    C1=max(C1,Trj[sequ[1]][M1])+FT[sequ[1]][M1];
    C2= max(C2,C1+lj[sequ[1]])+FT[sequ[1]][M2];
    Cmax= C2 + Tqj[sequ[1]][M2];
    for(i=2; i<= Nbrjf2-1; i++)
    {
        C1 = max(C1,Trj[sequ[i]][M1])+FT[sequ[i]][M1];
        C2 = max(C1+lj[sequ[i]],C2)+FT[sequ[i]][M2];
        Cmax = max(Cmax,C2+Tqj[sequ[i]][M2]);
    }
    Cmax=max(Cmax,max((C1+Q1),(C2+Q2)));
    return Cmax;
}

long CCMax(long seq[MAX_JOB],long nbjr,long R2,long R1,long Q1,long Q2,long M1, long M2)
{
    long c1=0,c2=(R2 - R1),h;
    for(h=1; h<=nbjr; h++)
    {
        c1= c1+ FT[seq[h]][M1];
        c2= max(c2,c1 + Tlj[seq[h]]) + FT[seq[h]][M2];
    }
    c2=max(c2,c1+(Q1-Q2));
    return c2;
}
long CCMax1(long seq[MAX_JOB],long nbjr,long R2,long R1,long Q1,long Q2,long M1, long M2)
{
    long c1=0,c2=(R2 - R1),h;
    for(h=1; h<=nbjr; h++)
    {
        c1= c1+ FT[seq[h]][M1];
        c2= max(c2,c1 + lj[seq[h]]) + FT[seq[h]][M2];
    }
    c2=max(c2,c1+(Q1-Q2));
    return c2;
}

long born_inf1OneM(long Nb_Job_Restant,long rm1bb,long qm1bb,long job)
{
    long l;
    long LB;
    long MINrj=1000000000,MINqj=1000000000;
    long t,SUM=0;
    t=1;
    for(l=1; l<=Nb_Job_Restant; l++)
    {
        if(FJR_b[l]!=job)
        {
            M1_j_bOne[t].job=FJR_b[l];
            M1_j_bOne[t].r=max(Trj[FJR_b[l]][M1],rm1bb);
            if(MINrj>M1_j_bOne[t].r)
                MINrj=M1_j_bOne[t].r;
            M1_j_bOne[t].p=FT[FJR_b[l]][M1];
            M1_j_bOne[t].q=max(Tqj[FJR_b[l]][M1],qm1bb);
            if(MINqj>M1_j_bOne[t].q)
                MINqj=M1_j_bOne[t].q;
            t++;
        }
    }
    for(l=1; l<=t-1; l++)
        SUM=SUM+M1_j_bOne[l].p;
    LB=MINrj+ SUM + MINqj;
    return (LB);
}
long born_inf2(long Nb_Job_Restant,long rm1bb,long rm2bb,long qm1bb,long qm2bb,long job,long Tin)
{
    long l,t,tail1,tail2,i,j;
    long LB,max,val1,tin,h;
    long minr1j=1000000000;
    long minq2j=1000000000;
    long minr2j=1000000000;
    long minq1j=1000000000;
    jjob k;
    t=1;
    tail1=1;
    tail2=1;
    for(l=1; l<=Nb_Job_Restant; l++)
    {
        if(l!=job)
        {
            SEQ[t].job=JR_b[l];
            SEQ[t].r1=max(Trj[JR_b[l]][M1],rm1bb);
            SEQ[t].p1=FT[JR_b[l]][M1];
            SEQ[t].alpha=FT[JR_b[l]][M1]+lj[JR_b[l]];
            SEQ[t].q2=max(Tqj[JR_b[l]][M2],qm2bb);
            SEQ[t].r2=max(SEQ[t].r1+FT[JR_b[l]][M1]+lj[JR_b[l]],rm2bb);
            SEQ[t].beta=FT[JR_b[l]][M2]+lj[JR_b[l]];
            SEQ[t].q1=max(SEQ[t].q2+FT[JR_b[l]][M2]+lj[JR_b[l]],qm1bb);
            SEQ[t].p2=FT[JR_b[l]][M2];
            if(SEQ[t].r1 < minr1j) minr1j = SEQ[t].r1;
            if(SEQ[t].q2 < minq2j) minq2j = SEQ[t].q2;
            if(SEQ[t].r2 < minr2j) minr2j = SEQ[t].r2;
            if(SEQ[t].q1 < minq1j) minq1j = SEQ[t].q1;
            if(SEQ[t].alpha <= SEQ[t].beta)
            {
                /* S1 ordre croissant sur alpha */
                S1[tail1].j=JR_b[l];
                S1[tail1].p=SEQ[t].alpha;
                tail1++;
            }
            else
            {
                /* S2 ordre decroissant sur beta */
                S2[tail2].j=JR_b[l];
                S2[tail2].p=SEQ[t].beta;
                tail2++;
            }
            t++;
        }
    }
    for(i=1; i<tail1; i++)
    {
        max=i;
        for(j=i+1; j<tail1; j++)
            if(S1[j].p<S1[max].p) max=j;
        k=S1[i];
        S1[i]=S1[max];
        S1[max]=k;
    }
    for(i=1; i<tail2; i++)
    {
        max=i;
        for(j=i+1; j<tail2; j++)
            if(S2[j].p>S2[max].p) max=j;
        k=S2[i];
        S2[i]=S2[max];
        S2[max]=k;
    }
    for(i=1; i<tail1; i++) SS1[i]=S1[i].j;
    for(i=1; i<tail2; i++)
    {
        SS1[tail1]=S2[i].j;
        tail1++;
    }
    tin=Tin+1;
    for(h=1; h <= t-1; h++)
    {
        SSeq[tin]= SS1[h];
        tin++;
    }
    BSjoh=cmaxfrql1(SSeq,rm1bb,rm2bb,qm1bb,qm2bb);
    val1= CCMax1(SS1,t-1,minr2j,minr1j,minq1j,minq2j,M1,M2);
    LB=minr1j + val1 +minq2j;
    return (LB);
}
long born_inf12(long Nb_Job_Restant,long rm1bb,long rm2bb,long qm1bb,long qm2bb,long job)
{
    long l,t,tail1,tail2,i,j;
    long LB,max,val1,h;
    long minr1j=1000000000;
    long minq2j=1000000000;
    long minr2j=1000000000;
    long minq1j=1000000000;
    jjob k;
    t=1;
    tail1=1;
    tail2=1;
    for(l=1; l<=Nb_Job_Restant; l++)
    {
        if(FJR_b[l]!=job)
        {
            SEQ[t].job=FJR_b[l];
            SEQ[t].r1=max(Trj[FJR_b[l]][M1],rm1bb);
            SEQ[t].p1=FT[FJR_b[l]][M1];
            SEQ[t].alpha=FT[FJR_b[l]][M1]+Tlj[FJR_b[l]];
            SEQ[t].q2=max(Tqj[FJR_b[l]][M2],qm2bb);
            SEQ[t].r2=max(SEQ[t].r1+FT[FJR_b[l]][M1]+Tlj[FJR_b[l]],rm2bb);
            SEQ[t].beta=FT[FJR_b[l]][M2]+Tlj[FJR_b[l]];
            SEQ[t].q1=max(SEQ[t].q2+FT[FJR_b[l]][M2]+Tlj[FJR_b[l]],qm1bb);
            SEQ[t].p2=FT[FJR_b[l]][M2];
            if(SEQ[t].r1 < minr1j) minr1j = SEQ[t].r1;
            if(SEQ[t].q2 < minq2j) minq2j = SEQ[t].q2;
            if(SEQ[t].r2 < minr2j) minr2j = SEQ[t].r2;
            if(SEQ[t].q1 < minq1j) minq1j = SEQ[t].q1;
            if(SEQ[t].alpha <= SEQ[t].beta)
            {
                /* S1 ordre croissant sur alpha */
                S1[tail1].j=FJR_b[l];
                S1[tail1].p=SEQ[t].alpha;
                tail1++;
            }
            else
            {
                /* S2 ordre décroissant sur beta */
                S2[tail2].j=FJR_b[l];
                S2[tail2].p=SEQ[t].beta;
                tail2++;
            }
            t++;
        }
    }
    Reste1=t-1;
    for(i=1; i<tail1; i++)
    {
        max=i;
        for(j=i+1; j<tail1; j++)
            if(S1[j].p<S1[max].p) max=j;
        k=S1[i];
        S1[i]=S1[max];
        S1[max]=k;
    }
    for(i=1; i<tail2; i++)
    {
        max=i;
        for(j=i+1; j<tail2; j++)
            if(S2[j].p>S2[max].p) max=j;
        k=S2[i];
        S2[i]=S2[max];
        S2[max]=k;
    }
    for(i=1; i<tail1; i++) SS1[i]=S1[i].j;
    for(i=1; i<tail2; i++)
    {
        SS1[tail1]=S2[i].j;
        tail1++;
    }
    for(h=1; h <= t-1; h++)
        SSSeq[h]= SS1[h];
    BSjoh=cmaxfrql(SSSeq,rm1bb,rm2bb,qm1bb,qm2bb);
    val1= CCMax(SS1,t-1,minr2j,minr1j,minq1j,minq2j,M1,M2);
    LB = minr1j + val1 +minq2j;
    return (LB);
}

long born_inf122(long Nb_Job_Restant,long rm1bb,long rm2bb,long qm1bb,long qm2bb,long job)
{
    long l,t,tail1,tail2,i,j;
    long LB,max,val1;
    long minr1j=1000000000;
    long minq2j=1000000000;
    long minr2j=1000000000;
    long minq1j=1000000000;
    jjob k;
    t=1;
    tail1=1;
    tail2=1;
    for(l=1; l<=Nb_Job_Restant; l++)
    {
        if(FJR_b[l]!=job)
        {
            SEQ[t].job=FJR_b[l];
            SEQ[t].r1=max(Trj[FJR_b[l]][M1],rm1bb);
            SEQ[t].p1=FT[FJR_b[l]][M1];
            SEQ[t].alpha=FT[FJR_b[l]][M1]+Tlj[FJR_b[l]];
            SEQ[t].q2=max(Tqj[FJR_b[l]][M2],qm2bb);
            SEQ[t].r2=max(SEQ[t].r1+FT[FJR_b[l]][M1]+Tlj[FJR_b[l]],rm2bb);
            SEQ[t].beta=FT[FJR_b[l]][M2]+Tlj[FJR_b[l]];
            SEQ[t].q1=max(SEQ[t].q2+FT[FJR_b[l]][M2]+Tlj[FJR_b[l]],qm1bb);
            SEQ[t].p2=FT[FJR_b[l]][M2];
            if(SEQ[t].r1 < minr1j) minr1j = SEQ[t].r1;
            if(SEQ[t].q2 < minq2j) minq2j = SEQ[t].q2;
            if(SEQ[t].r2 < minr2j) minr2j = SEQ[t].r2;
            if(SEQ[t].q1 < minq1j) minq1j = SEQ[t].q1;
            if(SEQ[t].alpha <= SEQ[t].beta)
            {
                /* S1 ordre croissant sur alpha */
                S1[tail1].j=FJR_b[l];
                S1[tail1].p=SEQ[t].alpha;
                tail1++;
            }
            else
            {
                /* S2 ordre decroissant sur beta */
                S2[tail2].j=FJR_b[l];
                S2[tail2].p=SEQ[t].beta;
                tail2++;
            }
            t++;
        }
    }
    Reste1=t-1;
    for(i=1; i<tail1; i++)
    {
        max=i;
        for(j=i+1; j<tail1; j++)
            if(S1[j].p<S1[max].p) max=j;
        k=S1[i];
        S1[i]=S1[max];
        S1[max]=k;
    }
    for(i=1; i<tail2; i++)
    {
        max=i;
        for(j=i+1; j<tail2; j++)
            if(S2[j].p>S2[max].p) max=j;
        k=S2[i];
        S2[i]=S2[max];
        S2[max]=k;
    }
    for(i=1; i<tail1; i++) SS1[i]=S1[i].j;
    for(i=1; i<tail2; i++)
    {
        SS1[tail1]=S2[i].j;
        tail1++;
    }

    val1= CCMax(SS1,t-1,minr2j,minr1j,minq1j,minq2j,M1,M2);
    LB = minr1j + val1 +minq2j;
    return (LB);
}

long F2ljCmax(long Nb_Job_Restant,long rm1bb,long rm2bb,long qm1bb,long qm2bb,long job)
{
    long l,t,tail1,tail2,i,j,UBF2lj;
    long max,h;
    jjob k;
    t=1;
    tail1=1;
    tail2=1;
    for(l=1; l<=Nb_Job_Restant; l++)
    {
        if(FJR_b[l]!=job)
        {
            F2SEQ[t].job=FJR_b[l];
            F2SEQ[t].r1=Trj[FJR_b[l]][M1];
            F2SEQ[t].p1=FT[FJR_b[l]][M1];
            F2SEQ[t].alpha=FT[FJR_b[l]][M1]+Tlj[FJR_b[l]];
            F2SEQ[t].q2=Tqj[FJR_b[l]][M2];
            F2SEQ[t].r2=F2SEQ[t].r1+FT[FJR_b[l]][M1]+Tlj[FJR_b[l]];
            F2SEQ[t].beta=FT[FJR_b[l]][M2]+Tlj[FJR_b[l]];
            F2SEQ[t].q1=F2SEQ[t].q2+FT[FJR_b[l]][M2]+Tlj[FJR_b[l]];
            F2SEQ[t].p2=FT[FJR_b[l]][M2];

            if(F2SEQ[t].alpha <= F2SEQ[t].beta)
            {
                /* S1 ordre croissant sur alpha */
                F2S1[tail1].j=FJR_b[l];
                F2S1[tail1].p=F2SEQ[t].alpha;
                tail1++;
            }
            else
            {
                /* S2 ordre decroissant sur beta */
                F2S2[tail2].j=FJR_b[l];
                F2S2[tail2].p=F2SEQ[t].beta;
                tail2++;
            }
            t++;
        }
    }
    Reste1=t-1;
    for(i=1; i<tail1; i++)
    {
        max=i;
        for(j=i+1; j<tail1; j++)
            if(F2S1[j].p<F2S1[max].p) max=j;
        k=F2S1[i];
        F2S1[i]=F2S1[max];
        F2S1[max]=k;
    }
    for(i=1; i<tail2; i++)
    {
        max=i;
        for(j=i+1; j<tail2; j++)
            if(F2S2[j].p>F2S2[max].p) max=j;
        k=F2S2[i];
        F2S2[i]=F2S2[max];
        F2S2[max]=k;
    }
    for(i=1; i<tail1; i++) F2SS1[i]=F2S1[i].j;
    for(i=1; i<tail2; i++)
    {
        F2SS1[tail1]=F2S2[i].j;
        tail1++;
    }
    for(h=1; h <= t-1; h++)
        F2SSSeq[h]= F2SS1[h];
    UBF2lj=cmaxfrql(F2SSSeq,rm1bb,rm2bb,qm1bb,qm2bb);
    return (UBF2lj);
}
long Makespan(long SEq[MAX_JOB],long tail,long m)
{
    long i,j,h,l,tab1[MAX_MACH],tab2[MAX_MACH],tab3[MAX_MACH],tcmax=0;
    for(j=1; j<=m; j++)
    {
        tcmax = tcmax + FT[SEq[1]][j];
        tab1[j]=tcmax;
    }
    for(i=2; i<=tail; i++)
    {
        for(h=1; h<=m; h++)
            tab2[h]=FT[SEq[i]][h];
        tab3[1]=tab1[1]+tab2[1];
        for(l=2; l<=m; l++)
            tab3[l]=max(tab1[l],tab3[l-1])+tab2[l];
        memcpy(tab1,tab3,sizeof(long)*(m+1));
        tcmax=tab3[m];
    }
    return tcmax;
}

long F2rjljqjCmax(long FNb_Job_Restant,long Fbranch,long FTabCmaxMIN[MAX_MACH],long FTabCmaxMOUT[MAX_MACH],long FTin,long FTout,long job)
{
    noeud * N;
    long born,i,j,k,nbrMC=0,h,M1M2,max,l;
    long BSSS=100000000,BSS=100000000,Bborn=0;
    long MM1,MM2,MaxLBOneM=0,LBOneM,nbpair,m1,m2;
    OM mk;
    Nbrjf2=FNb_Job_Restant;
    Bborn_sup=0;
    for(i=1; i<=n; i++)
    {
        for(j=1; j<=m; j++)
        {
            Trj[i][j]=0;
            Tqj[i][j]=0;
        }
    }
    for(j=1; j<=FNb_Job_Restant; j++)
    {
        if(j!=job)
            Trj[FJR_b[j]][1]=FTabCmaxMIN[1];
    }
    for(i=2; i<=m; i++)
    {
        for(j=1; j<=FNb_Job_Restant; j++)
        {
            if(j!=job)
                Trj[FJR_b[j]][i]=max(FTabCmaxMIN[i],Trj[FJR_b[j]][i-1]+FT[FJR_b[j]][i-1]);
        }
    }
    for(j=1; j<=FNb_Job_Restant; j++)
    {
        if(j!=job)
            Tqj[FJR_b[j]][m]=FTabCmaxMOUT[m];
    }
    for(i=m-1; i>=1; i--)
    {
        for(j=1; j<=FNb_Job_Restant; j++)
        {
            if(j!=job)
                Tqj[FJR_b[j]][i]=max(FTabCmaxMOUT[i],Tqj[FJR_b[j]][i+1]+FT[FJR_b[j]][i+1]);
        }
    }
    JobR=FJR_b[job];
    for(M1=1; M1<=m; M1++)
    {
        LBOneM=born_inf1OneM(FNb_Job_Restant ,FTabCmaxMIN[M1],FTabCmaxMOUT[M1],FJR_b[job]);
        FTPM[M1].m=M1;
        FTPM[M1].Binf=LBOneM;
        if(MaxLBOneM<LBOneM)
        {
            MaxLBOneM=LBOneM;
            if(MaxLBOneM>=FBSS)
                return MaxLBOneM;
        }
    }
    if(MaxLBOneM<FBSS)
    {
        M1=1;
        M2=m;
        born=F2ljCmax(FNb_Job_Restant,FTabCmaxMIN[M1],FTabCmaxMIN[M2],FTabCmaxMOUT[M1],FTabCmaxMOUT[M2],FJR_b[job]);
        if(MaxLBOneM<born)
            MaxLBOneM=born;
    }
    if(MaxLBOneM<FBSS)
    {
        k=1;
        for(h=1; h<=m; h++)
        {
            max=h;
            for(l=h+1; l<=m; l++)
                if(FTPM[l].Binf>FTPM[max].Binf) max=l;
            mk=FTPM[h];
            FTPM[h]=FTPM[max];
            FTPM[max]=mk;
        }

// No value provided for 'm==3' : lower bound not used
        switch(m)
        {
        case 5 : nbpair=3; break;
        case 10: nbpair=4; break;
        case 20: nbpair=8; break;
        default: /* do nothing */ break;
        }

        for(m1=1; m1<=nbpair; m1++)
        {
            for(m2=1; m2<=nbpair; m2++)
            {
                M1=FTPM[m1].m;
                M2=FTPM[m2].m;
                if((M1!=M2)&&(M1<M2) && ((M1!=1)&&(M2!=m)))
                {
                    for(i=1; i<=n; i++)
                    {
                        seq[i]=0;
                        Seq[i]=0;
                        SSeq[i]=0;
                        seq_op[i]=0;
                        JR_b[i]=0;
                        Tlj[i]=0;
                    }
                    for(j=1; j<=FNb_Job_Restant; j++)
                    {
                        if(j!=job)
                        {
                            for(i=M1+1; i<=M2-1; i++)
                                Tlj[FJR_b[j]]=Tlj[FJR_b[j]]+FT[FJR_b[j]][i];
                        }
                    }
                    born= born_inf12(FNb_Job_Restant ,FTabCmaxMIN[M1],FTabCmaxMIN[M2],FTabCmaxMOUT[M1],FTabCmaxMOUT[M2],FJR_b[job]);
                    BSS=BSjoh;
                    TPM[k].m1=M1;
                    TPM[k].m2=M2;
                    TPM[k].BinfP=born;
                    TPM[k].BSP=BSS;
                    memcpy(TPM[k].TTlj,Tlj,sizeof(long)*(n+1));
                    k++;
                    if(born>Bborn)
                    {
                        Bborn=born;
                        memcpy(llj,Tlj,sizeof(long)*(n+1));
                        MM1=M1;
                        MM2=M2;
                        BSSS=BSS;
                    }
                }
            }
        }
        if (Bborn<FBSS)
        {
            BSS=100000000;
            for(h=1; h<=(nbpair*(nbpair-1))/2; h++)
            {
                if((TPM[h].BSP >= Bborn)&&(TPM[h].BinfP >= Bborn))
                {
                    if(TPM[h].BSP<BSS)
                    {
                        BSS=TPM[h].BSP;
                        born=TPM[h].BinfP;
                        M1=TPM[h].m1;
                        M2=TPM[h].m2;
                        memcpy(lj,TPM[h].TTlj,sizeof(long)*(n+1));
                        nbrMC++;
                    }
                }
            }
        }
        if(nbrMC==0)
        {
            born=Bborn;
            M1=MM1;
            M2=MM2;
            BSS=BSSS;
            memcpy(lj,llj,sizeof(long)*(n+1));
        }
        if(BSS==born) Bborn_sup=BSS;
        if((BSS==born+1)||(BSS==born+2)) Bborn_sup=born;
        if((BSS>born+2) && (BSS > MaxLBOneM))
        {
            born_sup=100000000;
            if(BSS<born_sup)
            {
                born_sup=BSS;
                Bborn_sup=BSS;
            }
            if(born<FBSS)
            {
                BB++;
                MINbinf=100000000;
                if(Fbranch==FIN)
                {
                    BOOLF2=0;
                    TOUL=0;
                    nbrfcrnd=0;
                    Sstart=clock();
                    N=create_noeud(FNb_Job_Restant,0,0,0,FTabCmaxMIN[M1],FTabCmaxMIN[M2],FTabCmaxMOUT[M1],FTabCmaxMOUT[M2],born,JobR);
                    explorer(N);
                    EEend=clock();
                    TimeBB=TimeBB+(EEend-Sstart);
                    FbranchF2=1;
                    if(BOOLF2==1)
                        Bborn_sup=MINbinf;
                }
                else
                {
                    BOOLF2=0;
                    memcpy(Trjqj,Trj,sizeof(long)*((m+1)*(n+1)));
                    memcpy(Trj,Tqj,sizeof(long)*((m+1)*(n+1)));
                    memcpy(Tqj,Trjqj,sizeof(long)*((m+1)*(n+1)));
                    M1M2=M1;
                    M1=M2;
                    M2=M1M2;
                    TOUL=0;
                    nbrfcrnd=0;
                    Sstart=clock();
                    N=create_noeud(FNb_Job_Restant,0,0,0,FTabCmaxMOUT[M1],FTabCmaxMOUT[M2],FTabCmaxMIN[M1],FTabCmaxMIN[M2],born,JobR);
                    explorer(N);
                    EEend=clock();
                    TimeBB=TimeBB+(EEend-Sstart);
                    FbranchF2=1;
                    if(BOOLF2==1)
                        Bborn_sup=MINbinf;
                }
            }
            else
            {
                born=max(born,MaxLBOneM);
                return born;
            }
        }
    }
    Bborn_sup=max(Bborn_sup,MaxLBOneM);
    return Bborn_sup;
}
long F2rjljqjCmax1(long FNb_Job_Restant,long Fbranch,long FTabCmaxMIN[MAX_MACH],long FTabCmaxMOUT[MAX_MACH],long FTin,long FTout,long job)
{
    long i,j;
    long BSSS=100000000,BSS=100000000,Bborn=0;
    long MaxLBOneM=0,LBOneM;
    Nbrjf2=FNb_Job_Restant;
    Bborn_sup=0;
    for(i=1; i<=n; i++)
    {
        for(j=1; j<=m; j++)
        {
            Trj[i][j]=0;
            Tqj[i][j]=0;
        }
    }
    for(j=1; j<=FNb_Job_Restant; j++)
    {
        if(j!=job)
            Trj[FJR_b[j]][1]=FTabCmaxMIN[1];
    }
    for(i=2; i<=m; i++)
    {
        for(j=1; j<=FNb_Job_Restant; j++)
        {
            if(j!=job)
                Trj[FJR_b[j]][i]=max(FTabCmaxMIN[i],Trj[FJR_b[j]][i-1]+FT[FJR_b[j]][i-1]);
        }
    }
    for(j=1; j<=FNb_Job_Restant; j++)
    {
        if(j!=job)
            Tqj[FJR_b[j]][m]=FTabCmaxMOUT[m];
    }
    for(i=m-1; i>=1; i--)
    {
        for(j=1; j<=FNb_Job_Restant; j++)
        {
            if(j!=job)
                Tqj[FJR_b[j]][i]=max(FTabCmaxMOUT[i],Tqj[FJR_b[j]][i+1]+FT[FJR_b[j]][i+1]);
        }
    }
    JobR=FJR_b[job];
    for(M1=1; M1<=m; M1++)
    {
        LBOneM=born_inf1OneM(FNb_Job_Restant,FTabCmaxMIN[M1],FTabCmaxMOUT[M1],FJR_b[job]);
        if(MaxLBOneM<LBOneM)
        {
            MaxLBOneM=LBOneM;
            if(MaxLBOneM>=FBSS)
                break;
        }
    }
    return MaxLBOneM;
}
long F2rjljqjCmax2(long FNb_Job_Restant,long Fbranch,long FTabCmaxMIN[MAX_MACH],long FTabCmaxMOUT[MAX_MACH],long FTin,long FTout,long job)
{
    long born,i,j,nbrMC=0,h,l,max,nbpair,m1,m2;
    long BSSS=100000000,BSS=100000000,Bborn=0;
    long MaxLBOneM=0,LBOneM;
    OM mk;
    Nbrjf2=FNb_Job_Restant;
    Bborn_sup=0;
    for(i=1; i<=n; i++)
    {
        for(j=1; j<=m; j++)
        {
            Trj[i][j]=0;
            Tqj[i][j]=0;
        }
    }
    for(j=1; j<=FNb_Job_Restant; j++)
    {
        if(j!=job)
            Trj[FJR_b[j]][1]=FTabCmaxMIN[1];
    }
    for(i=2; i<=m; i++)
    {
        for(j=1; j<=FNb_Job_Restant; j++)
        {
            if(j!=job)
                Trj[FJR_b[j]][i]=max(FTabCmaxMIN[i],Trj[FJR_b[j]][i-1]+FT[FJR_b[j]][i-1]);
        }
    }
    for(j=1; j<=FNb_Job_Restant; j++)
    {
        if(j!=job)
            Tqj[FJR_b[j]][m]=FTabCmaxMOUT[m];
    }
    for(i=m-1; i>=1; i--)
    {
        for(j=1; j<=FNb_Job_Restant; j++)
        {
            if(j!=job)
                Tqj[FJR_b[j]][i]=max(FTabCmaxMOUT[i],Tqj[FJR_b[j]][i+1]+FT[FJR_b[j]][i+1]);
        }
    }
    JobR=FJR_b[job];
    for(M1=1; M1<=m; M1++)
    {
        LBOneM=born_inf1OneM(FNb_Job_Restant ,FTabCmaxMIN[M1],FTabCmaxMOUT[M1],FJR_b[job]);
        FTPM[M1].m=M1;
        FTPM[M1].Binf=LBOneM;
        if(MaxLBOneM<LBOneM)
        {
            MaxLBOneM=LBOneM;
            if(MaxLBOneM>=FBSS)
                return MaxLBOneM;
        }
    }
    if(MaxLBOneM<FBSS)
    {
        M1=1;
        M2=m;
        born=F2ljCmax(FNb_Job_Restant,FTabCmaxMIN[M1],FTabCmaxMIN[M2],FTabCmaxMOUT[M1],FTabCmaxMOUT[M2],FJR_b[job]);
        if(MaxLBOneM<born) MaxLBOneM=born;
    }
    if(MaxLBOneM<FBSS)
    {
        for(h=1; h<=m; h++)
        {
            max=h;
            for(l=h+1; l<=m; l++)
                if(FTPM[l].Binf>FTPM[max].Binf) max=l;
            mk=FTPM[h];
            FTPM[h]=FTPM[max];
            FTPM[max]=mk;
        }

// No value provided for 'm==3' : lower bound not used
        switch(m)
        {
        case 5 : nbpair=3; break;
        case 10: nbpair=4; break;
        case 20: nbpair=8; break;
        default: /* do nothing */ break;
        }

        for(m1=1; m1<=nbpair; m1++)
        {
            for(m2=1; m2<=nbpair; m2++)
            {
                M1=FTPM[m1].m;
                M2=FTPM[m2].m;
                if((M1!=M2)&&(M1<M2) && ((M1!=1)&&(M2!=m)))
                {
                    for(i=1; i<=n; i++)
                    {
                        seq[i]=0;
                        Seq[i]=0;
                        SSeq[i]=0;
                        seq_op[i]=0;
                        JR_b[i]=0;
                        Tlj[i]=0;
                    }
                    for(j=1; j<=FNb_Job_Restant; j++)
                    {
                        if(j!=job)
                        {
                            for(i=M1+1; i<=M2-1; i++)
                                Tlj[FJR_b[j]]=Tlj[FJR_b[j]]+FT[FJR_b[j]][i];
                        }
                    }
                    born=born_inf122(FNb_Job_Restant ,FTabCmaxMIN[M1],FTabCmaxMIN[M2],FTabCmaxMOUT[M1],FTabCmaxMOUT[M2],FJR_b[job]);
                }
            }
        }
    }
    else return MaxLBOneM;
    born=max(born,MaxLBOneM);
    return born;
}

// Replaces the original 'main()' function
unsigned int branchAndBound(const int mym, const int myn, const long* p)
{
    Fnoeud * root; //< A pointer to the root of the research tree
    long lowerBound; //< A parameter of the branch & bound algorithm
    int i, //< A machine index
        j, //< A job ID
        k; //< A job position
    
    m=mym;
    n=myn;
// Initialization : read the processing times
    for(j=0; j<n; ++j)
        for(i=0; i<m; ++i){
            FT[j+1][i+1] = p[j*m + i];
            //printf("%ld,", FT[j+1][i+1]);
        }

// Initialization : set initial bounds value (parameters of the algorithm)
    lowerBound = 0; //< Local variable
    FBSS = 100000000; //< Global variable ; same value as specified during the declaration
    Fborn_sup = 100000000; //< Global variable ; same value as specified during the declaration

// Initialization : set different levels
// When used, these levels indicate which lower bound to use among the 3 implemented.
// Here, they are not used : only the lower bound that relaxes constraints on all machines but one is used.
    switch(n)
    {
    case 20 : Niv1=7  ; Niv2=14 ; break;
    case 50 : Niv1=23 ; Niv2=43 ; break;
    case 100: Niv1=45 ; Niv2=85 ; break;
    case 200: Niv1=85 ; Niv2=170; break;
    case 500: Niv1=230; Niv2=430; break;
    default : /* do nothing */    break;
    }

// Other initializations
    MINbinfcmax=lowerBound;
    ERROR=0.01;
    Fnb_iter=0;
    Fnb_exp=0;
    Fnbnc=0;
    TimeBB=0;
    BB=0;
    Error=0;
    memset(TabCmaxMIN, 0, MAX_MACH * sizeof(long));
    memset(TabCmaxMOUT, 0, MAX_MACH * sizeof(long));

// Branch & Bound algorithm
    root = Fcreate_noeud(n,0,0,0,TabCmaxMIN,TabCmaxMOUT,lowerBound);
    Fexplorer(root);
    free(root);

// Finalization : write the optimal sequence found
    /*for(k=1; k<=n; ++k)
    {
        wrapper_Sequence_addJob(sequenceLocation, dataLocation, Fseq_op[k] - 1);
    }*/

// Return the minimal value of the criterion
    return Fborn_sup;
}
