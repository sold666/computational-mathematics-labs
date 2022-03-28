/* File converted with FORTRAN CONVERTER utility.
   FORTRAN CONVERTER is written by Grigoriev D., 2081/4.*/

/*MARKERS: //! = probably incorrect, //? = possibly incorrect*/

#include "FORSYTHE.H"


void SVD(int M,int N,REAL **A,REAL *W,bool MATU,REAL **U,bool MATV,REAL **V,int &IERR,REAL *RV1)
{

int I,J,K,L,II,I1,KK,K1,LL,L1,ITS,MN;
REAL C,F,G,H,S,X,Y,Z,SCALE,ANORM;


//     �TA �O��PO�PAMMA ECT� TPAHC����� A��O�-�PO�E��P� SVD,
//     O�����KOBAHHO� �O���OM � PA�H�EM B ��PHA�E NUMERISCHE
//     MATHEMATIK, 14, 403-420(1970), A TAK�E B KH��E
//     HANDBOOK FOR AUTOMATIC COMPUTATION  VOL.II-LINEAR
//     ALGEBRA, 134-151 (1971).

//     �TA �O��PO�PAMMA B���C��ET C�H����PHOE PA��O�EH�E
//     A=USV'  �E�CTB�TE��HO� �P�MO��O��HO� MATP��� A C PA�ME-
//     PAM� M � N. �P� �TOM �C�O�����TC� �B�X��A�OHA���A���
//     �OCPE�CTBOM XA�CXO��EPOB�X OTPA�EH��  � BAP�AHT
//     QR-A��OP�TMA.

//     HA BXO�E.

//     NM   C�E��ET �O�O��T� PABH�M CTPO�HO� PA�MEPHOCT�
//          �B�XMEPH�X MACC�BOB ,�A�B�EHHO� B O�EPATOPE PA�-
//          MEPHOCT� B���BA��E� �PO�PAMM�. �AMET�M, �TO NM
//          �O��HO ��T� H  MEH��E,�EM MAKC�M�M �� M � N.

//     M    ��C�O CTPOK A (� U).

//     N    ��C�O CTO���OB A (� U) � �OP��OK V.

//     A    CO�EP��T �P�MO��O��H�� BXO�H�� MATP���, ���
//          KOTOPO� HAXO��TC�  PA��O�EH�E.

//     MATU �O��EH �MET� �HA�EH�E true, EC�� H��HO
//          B���C��T� MATP��� U �� PA��O�EH��, � �HA�EH�E
//          false B �POT�BHOM C���AE.

//     MATV �O��EH �MET� �HA�EH�E true, EC�� H��HO
//          B���C��T� MATP��� V �� PA��O�EH��, � �HA�EH�E
//          false B �POT�BHOM C���AE.

//     HA B�XO�E.

//     A    HE ��MEH�ETC� (EC�� HA EE MECTE HE �A��C�BA�TC�
//          U ���O V).

//     W    CO�EP��T N (HEOTP��ATE��H�X) C�H����PH�X ��CE�
//          A (��A�OHA��H�X ��EMEHTOB S). OH� HE ��OP��O�EH�.
//          EC�� �PO�CXO��T B�XO� �O O���KE, TO ��� �HA�EH��
//          IERR+1, IERR+2,...,N C�H�� ��PH�E ��C�A �O��H�
//          ��T� BEPH�.

//     U    CO�EP��T MATP��� U (C OPTO�OHA��H�M� CTO���AM�)
//          �� PA��O�EH��, EC�� ��� �APAMETPA MATU ���O
//          �A�AHO �HA�EH�E true B �POT�BHOM C���AE HA U
//          �C�O����ETC� KAK BPEMEHH�� MACC�B. U TAK�E MO�ET
//          COB�A�AT� C A. EC�� �PO�CXO��T B�XO�
//          �O O���KE, TO CTO���� U, COOTBETCTB����E �H�EKCAM
//          BEPH�X C�H����PH�X ��CE�, �O��H� ��T� TAK�E BEPH�.

//     V    CO�EP��T MATP��� V (OPTO�OHA��H��) �� PA��O�EH��,
//          EC�� ��� �APAMETPA MATV ���O �A�AHO �HA�EH�E
//          true B �POT�BHOM C���AE HA V HE �PO��BO��TC�
//          CC��OK. V TAK�E MO�ET COB�A�AT� C  A, EC�� U HE
//          B���C��ETC�. EC�� �PO�CXO��T B�XO� �O O���KE,
//          TO CTO���� V, COOTBETCTB����E �H�EKCAM BEPH�X
//          C�H����PH�X ��CE�, �O��H� ��T� TAK�E BEPH�.

//     IERR PABHO
//            0,  EC�� �PO�CXO��T HOPMA��H�� B�XO� �� �O��PO�-
//                PAMM�,
//            K,  EC�� K-E C�H����PHOE ��C�O HE ���O O�PE�E-
//                �EHO �OC�E 30 �TEPA���.

//     RV1  �TO MACC�B �POME��TO�HO�O XPAHEH��.

//     BO�POC� � KOMMEHTAP�� H��HO HA�PAB��T� �O A�PEC�
//     B.S.GARBOW, APPLIED MATEMATICS DIVISION, ARGONNE
//     NATIONAL LABORATORY

//     �O��PO�PAMMA MO������POBAHA C �E��� �CK����T�
//     �EPEMEHHY� MACHEP

IERR=0;
for( I=1; I<=M; I++){   //? target=100
for( J=1; J<=N; J++){   //? target=100
U[I-1][J-1]=A[I-1][J-1];
//_100:;
}
}                         	// CONTINUE

//     XA�CXO��EPOBO �P�BE�EH�E K �B�X��A�OHA��HO� �OPME

G=0.0;
SCALE=0.0;
ANORM=0.0;

for( I=1; I<=N; I++){   //? target=300
L=I+1;
RV1[I-1]=SCALE*G;
G=0.0;
S=0.0;
SCALE=0.0;
if(I>M)goto _210;

for( K=I; K<=M; K++)SCALE=SCALE+ABS(U[K-1][I-1]);
if(SCALE==0.0)goto _210;

for( K=I; K<=M; K++){   //? target=130
U[K-1][I-1]=U[K-1][I-1]/SCALE;
S=S+U[K-1][I-1]*U[K-1][I-1];
//_130:;
}                         	// CONTINUE

F=U[I-1][I-1];
G=-SIGN(SQRT(S),F);
H=F*G-S;
U[I-1][I-1]=F-G;
if(I==N)goto _190;

for( J=L; J<=N; J++){   //? target=150
S=0.0;
for( K=I; K<=M; K++)S=S+U[K-1][I-1]*U[K-1][J-1];
F=S/H;
for( K=I; K<=M; K++){   //? target=150
U[K-1][J-1]=U[K-1][J-1]+F*U[K-1][I-1];
//_150:;
}
}                         	// CONTINUE

_190:;
for( K=I; K<=M; K++)U[K-1][I-1]=SCALE*U[K-1][I-1];
_210:;
W[I-1]=SCALE*G;
G=0.0;
S=0.0;
SCALE=0.0;
if(I>M||I==N)goto _290;

for( K=L; K<=N; K++)SCALE=SCALE+ABS(U[I-1][K-1]);

if(SCALE==0.0)goto _290;

for( K=L; K<=N; K++){   //? target=230
U[I-1][K-1]=U[I-1][K-1]/SCALE;
S=S+U[I-1][K-1]*U[I-1][K-1];
//_230:;
}                         	// CONTINUE

F=U[I-1][L-1];
G=-SIGN(SQRT(S),F);
H=F*G-S;
U[I-1][L-1]=F-G;

for( K=L; K<=N; K++)RV1[K-1]=U[I-1][K-1]/H;

if(I==M)goto _270;

for( J=L; J<=M; J++){   //? target=260
S=0.0;
for( K=L; K<=N; K++)S=S+U[J-1][K-1]*U[I-1][K-1];
for( K=L; K<=N; K++){   //? target=260
U[J-1][K-1]=U[J-1][K-1]+S*RV1[K-1];
//_260:;
}
}                         	// CONTINUE

_270:;
for( K=L; K<=N; K++)U[I-1][K-1]=SCALE*U[I-1][K-1];

_290:;
ANORM=MAX(ANORM,ABS(W[I-1])+ABS(RV1[I-1]));
//_300:;
}                         	// CONTINUE

//     HAKO��EH�E �PABOCTOPOHH�X �PEO�PA�OBAH��

if(!MATV)goto _410;

//     ��� I=N C �A�OM -1 �O 1 B��O�H�T� -

for( II=1; II<=N; II++){   //? target=400
I=N+1-II;
if(I==N)goto _390;
if(G==0.0)goto _360;

for( J=L; J<=N; J++){   //? target=320

//     �BO�HOE �E�EH�E O�XO��T BO�MO�H�� MA��HH�� HY��

//_320:;
V[J-1][I-1]=(U[I-1][J-1]/U[I-1][L-1])/G;
}

for( J=L; J<=N; J++){   //? target=350
S=0.0;
for( K=L; K<=N; K++)S=S+U[I-1][K-1]*V[K-1][J-1];
for( K=L; K<=N; K++){   //? target=350
V[K-1][J-1]=V[K-1][J-1]+S*V[K-1][I-1];
//_350:;
}
}                         	// CONTINUE

_360:;
for( J=L; J<=N; J++){   //? target=380
V[I-1][J-1]=0.0;
V[J-1][I-1]=0.0;
//_380:;
}                         	// CONTINUE

_390:;
V[I-1][I-1]=1.0;
G=RV1[I-1];
L=I;
//_400:;
}                         	// CONTINUE

//     HAKO��EH�E �EBOCTOPOHH�X �PEO�PA�OBAH��

_410:;
if(!MATU)goto _510;

//     ��� I=MIN(N,M) C �A�OM -1 �O 1 B��O�H�T�-

MN=N;
if(M<N)MN=M;

for( II=1; II<=MN; II++){   //? target=500
I=MN+1-II;
L=I+1;
G=W[I-1];
if(I==N)goto _430;

for( J=L; J<=N; J++)U[I-1][J-1]=0.0;

_430:;
if(G==0.0)goto _475;
if(I==MN)goto _460;

for( J=L; J<=N; J++){   //? target=450
S=0.0;
for( K=L; K<=M; K++)S=S+U[K-1][I-1]*U[K-1][J-1];

//     �BO�HOE  �E�EH�E O�XO��T BO�MO�H�� MA��HH�� HY��

F=(S/U[I-1][I-1])/G;

for( K=I; K<=M; K++){   //? target=450
U[K-1][J-1]=U[K-1][J-1]+F*U[K-1][I-1];
//_450:;
}
}                         	// CONTINUE
_460:;
for( J=I; J<=M; J++)U[J-1][I-1]=U[J-1][I-1]/G;

goto _490;

_475:;
for( J=I; J<=M; J++)U[J-1][I-1]=0.0;

_490:;
U[I-1][I-1]=U[I-1][I-1]+1.0;
//_500:;
}                         	// CONTINUE

//     ��A�OHA���A��� �B�X��A�OHA��HO� �OPM� ��� K=N C �A�OM
//     -1 �O 1 B��O�H�T�

_510:;
for( KK=1; KK<=N; KK++){   //? target=700
K1=N-KK;
K=K1+1;
ITS=0;

//     �POBEPKA BO�MO�HOCT� PAC�E��EH��  ��� L=K
//     C �A�OM -1 �O 1 B��O�H�T�

_520:;
for( LL=1; LL<=K; LL++){   //? target=530
L1=K-LL;
L=L1+1;
if(ABS(RV1[L-1])+ANORM==ANORM)goto _565;

//     RV1[0] BCE��A PABHO H���. �O�TOM� B�XO�A
//     �EPE� KOHE� ��K�A HE ���ET

if(ABS(W[L1-1])+ANORM==ANORM)goto _540;
//_530:;
}                         	// CONTINUE

//      EC�� L �O���E, �EM 1, TO RV1[L-1]
//      �P�CBA�BAETC� H��EBOE �HA�EH�E

_540:;
C=0.0;
S=1.0;

for( I=L; I<=K; I++){   //? target=560
F=S*RV1[I-1];
RV1[I-1]=C*RV1[I-1];
if(ABS(F)+ANORM==ANORM)goto _565;
G=W[I-1];
H=SQRT(F*F+G*G);
W[I-1]=H;
C=G/H;
S=-F/H;
if(!MATU)goto _560;

for( J=1; J<=M; J++){   //? target=550
Y=U[J-1][L1-1];
Z=U[J-1][I-1];
U[J-1][L1-1]=Y*C+Z*S;
U[J-1][I-1]=-Y*S+Z*C;
//_550:;
}                         	// CONTINUE
_560:;
}                         	// CONTINUE

//     �POBEPKA CXO��MOCT�

_565:;
Z=W[K-1];
if(L==K)goto _650;

//     C�B�� B���PAETC� �� H��HE�O ���OBO�O
//     M�HOPA �OP��KA 2

if(ITS==30)goto _1000;
ITS=ITS+1;
X=W[L-1];
Y=W[K1-1];
G=RV1[K1-1];
H=RV1[K-1];
F=((Y-Z)*(Y+Z)+(G-H)*(G+H))/(2.0*H*Y);
G=SQRT(F*F+1.0);
F=((X-Z)*(X+Z)+H*(Y/(F+SIGN(G,F))-H))/X;

//     C�E����EE QR-�PEO�PA�OBAH�E

C=1.0;
S=1.0;

for( I1=L; I1<=K1; I1++){   //? target=600
I=I1+1;
G=RV1[I-1];
Y=W[I-1];
H=S*G;
G=C*G;
Z=SQRT(F*F+H*H);
RV1[I1-1]=Z;
C=F/Z;
S=H/Z;
F=X*C+G*S;
G=-X*S+G*C;
H=Y*S;
Y=Y*C;
if(!MATV)goto _575;

for( J=1; J<=N; J++){   //? target=570
X=V[J-1][I1-1];
Z=V[J-1][I-1];
V[J-1][I1-1]=X*C+Z*S;
V[J-1][I-1]=-X*S+Z*C;
//_570:;
}                         	// CONTINUE

_575:;
Z=SQRT(F*F+H*H);
W[I1-1]=Z;

//     BPA�EH�E MO�ET ��T� �PO��BO��H�M, EC�� Z PABHO H���

if(Z==0.0)goto _580;
C=F/Z;
S=H/Z;
_580:;
F=C*G+S*Y;
X=-S*G+C*Y;
if(!MATU)goto _600;

for( J=1; J<=M; J++){   //? target=590
Y=U[J-1][I1-1];
Z=U[J-1][I-1];
U[J-1][I1-1]=Y*C+Z*S;
U[J-1][I-1]=-Y*S+Z*C;
//_590:;
}                         	// CONTINUE
_600:;
}                         	// CONTINUE
RV1[L-1]=0.0;
RV1[K-1]=F;
W[K-1]=X;
goto _520;

//     CXO��MOCT�

_650:;
if(Z>=0.0)goto _700;

//     W[K-1] �E�AETC� HEOTP��ATE��H�M

W[K-1]=-Z;
if(!MATV)goto _700;

for( J=1; J<=N; J++)V[J-1][K-1]=-V[J-1][K-1];
_700:;
}                         	// CONTINUE
goto _1001;

//     �CTAHOB�T� �HA�EH�E �P��HAKA O���K� - �OC�E 30
//     �TEPA��� HET CXOMOCT� K C�H����PHOM� ��C��

_1000:;
IERR=K;
_1001:;
return;
}                         	// END
