/* File converted with FORTRAN CONVERTER utility.
   FORTRAN CONVERTER is written by Grigoriev D., 2081/4.*/

/*MARKERS: //! = probably incorrect, //? = possibly incorrect*/

#include "FORSYTHE.H"

//     METOÄ PÓHÃE-KÓTTA ÔEËÜÁEPÃA ×ETBEPTOÃO-ÏßTOÃO ÏOPßÄKA

//     COCTABÈTEËÈ ÏPOÃPAMMÛ-H.A.WATTS,L.F.SHAMPINE
//                SANDIA LABORATORIES
//              ALBUQUERQUE, NEW MEXICO

//     RKF45 ÏPEÄHAÇHA×EHA ÃËABHÛM OÁPAÇOM ÄËß PEØEHÈß
//     HEÆECTKÈX È CËAÁO ÆECTKÈX ÄÈÔÔEPEHÖÈAËÜHÛX ÓPABHEHÈÉ,
//     KOÃÄA BÛ×ÈCËEHÈE ÏPOÈÇBOÄHÛX HE CËÈØKOM ÄOPOÃOCTOßÙEE.
//     RKF45, BOOÁÙE ÃOBOPß,HE CËEÄÓET ÈCÏOËÜÇOBATÜ
//     ECËÈ ÏOËÜÇOBATEËÞ TPEÁÓETCß BÛCOKAß TO×HOCTÜ

//     PEÇÞME

//     ÏOÄÏPOÃPAMMA RKF45 ÈHTEÃPÈPÓET CÈCTEMÓ ÈÇ NEQN
//     OÁÛKHOBEHHÛX ÄÈÔÔEPEHÖÈAËÜHÛX ÓPABHEHÈÉ ÏEPBOÃO
//     ÏOPßÄKA CËEÄÓÞÙEÃO BÈÄA:

//            DY(I)/DT=F(T,Y(1),Y(2),...,Y(NEQN),

//     ÃÄE Y(I) ÇAÄAHÛ B T.
//     OÁÛ×HO ÏOÄÏPOÃPAMMÓ ÏPÈMEHßÞT ÄËß ÈHTEÃPÈPOBAHÈß
//     OT T ÄO TOUT, OÄHAKO EE MOÆHO ÈCÏOËÜÇOBATÜ È KAK
//     OÄHOØAÃOBÛÉ ÈHTEÃPATOP,×TOÁÛ ÏPOÄOËÆÈTÜ PEØEHÈE HA
//     OÄÈH ØAÃ B HAÏPABËEHÈÈ TOUT.HA BÛXOÄE ÏAPAMETPAM,
//     ÔÈÃÓPÈPÓÞÙÈM B CÏÈCKE BÛÇOBA, ÏPÈCBAÈBAÞTCß ÇHA×EHÈß,
//     HEOÁXOÄÈMÛE ÄËß ÏPOÄOËÆEHÈß ÈHTEÃPÈPOBAHÈß. ÏOËÜÇO-
//     BATEËÞ HÓÆHO ËÈØÜ EÙE PAÇ OÁPATÈTÜCß K RKF45
//     (È,BOÇMOÆHO,OÏPEÄEËÈTÜ HOBOE ÇHA×EHÈE ÄËß TOUT).
//     B ÄEÉCTBÈTEËÜHOCTÈ RKF45-ÝTO ÏPOÃPAMMA ÈHTEPÔEÉCA,
//     KOTOPAß BÛÇÛBAET ÏOÄÏPOÃPAMMÓ RKFS, OCÓÙECTBËßÞÙÓÞ
//     ÏPOÖECC PEØEHÈß.RKFS B CBOÞ O×EPEÄÜ BÛÇÛBAET ÏOÄÏPOÃ-
//     PAMMÓ FEHL, KOTOPAß BÛ×ÈCËßET ÏPÈÁËÈÆEHHOE PEØEHÈE
//     HA OÄÈH ØAÃ.

//     RKF45 ÈCÏOËÜÇÓET METOÄ PÓHÃE-KÓTTA-ÔEËÜÁEPÃA, OÏÈCAHHÛÉ
//     B CËEÄÓÞÙEÉ ÏÓÁËÈKAÖÈÈ:E.FEHLBERG,LOW-ORDER CLASSICAL
//     RUNGE-KUTTA FORMULAS WITH STEPSIZE CONTROL,NASA TR R-315

//     CTÈËÜ PAÁOTÛ ÏPOÃPAMMÛ RKF45 ÈËËÞCTPÈPÓETCß B CËEÄÓÞÙÈX
//     ÏÓÁËÈKAÖÈßX:L.F.SHAMPINE,H.A.WATTS,S.DAVENPORT, SOLVING
//     NON-STIFF ORDINARY DIFFERENTIAL EQUATIONS-THE STAT OF
//     THE ART,SANDIA LABORATORIES REPORT SAND75-0182,SIAM
//     REVIEW,18(1976), N3,376-411.

//     ÏAPAMETPÛ ÏPOÃPAMMÛ:

//     F       -ÏOÄÏPOÃPAMMA F(T,Y,YP) ÄËß BÛ×ÈCËEHÈß
//              ÏPOÈÇBOÄHÛX YP(I)=DY(I)/DT
//     NEQN    -×ÈCËO ÈHTEÃPÈPÓEMÛX ÓPABHEHÈÉ
//     Y(*)    -PEØEHÈE B TO×KE T
//     T       -HEÇABÈCÈMAß ÏEPEMEHHAß
//     TOUT    -TO×KA BÛXOÄA,B KOTOPOÉ HÓÆHO OÏPEÄEËÈTÜ
//              ÇHA×EHÈE PEØEHÈß
//     RELERR  -ÃPAHÈÖA OTHOCÈTEËÜHOÉ ÏOÃPEØHOCTÈ
//              ÄËß TECTA ËOKAËÜHOÉ OØÈÁKÈ.
//     ABSERR  -ÃPAHÈÖA ABCOËÞTHOÉ ÏOÃPEØHOCTÈ
//              ÄËß TECTA ËOKAËÜHOÉ OØÈÁKÈ.
//              HA KAÆÄOM ØAÃE ÏPOÃPAMMA TPEÁÓET BÛÏOËHEHÈß ÓCËOBÈß
//              ABS(LOCAL ERROR)<=RELERR*ABS(Y)+ABSERR
//              ÄËß KAÆÄOÉ KOMÏOHEHTÛ BEKTOPOB ËOKAËÜHOÉ
//              OØÈÁKÈ È PEØEHÈß
//     IFLAG   -ÓKAÇATEËÜ PEÆÈMA ÈHTEÃPÈPOBAHÈß.
//     WORK(*) -MACCÈB,COÄEPÆAÙÈÉ ÈHÔOPMAÖÈÞ,BHÓTPEHHÞÞ ÄËß RKF45,
//              KOTOPAß HEOÁXOÄÈMA ÏPÈ ÏOCËEÄÓÞÙÈX BÛÇOBAX.EÃO
//              PAÇMEPHOCTÜ ÄOËÆHA ÁÛTÜ HE MEHÜØE 3+6*NEQN
//     IWORK(*)-ÖEËÛÉ MACCÈB,COÄEPÆAÙÈÉ ÈHÔOPMAÖÈÞ,BHÓTPEHHÞÞ ÄËß
//              RKF45,KOTOPAß HEOÁXOÄÈMA ÏPÈ ÏOCËEÄÓÞÙÈX BÛÇOBAX.
//              EÃO PAÇMEPHOCTÜ ÄOËÆHA ÁÛTÜ HE MEHÜØE 5.

//     ÏEPBOE OÁPAÙEHÈE K RKF45

//     ÏOËÜÇOBATEËÜ ÄOËÆEH ÏPEÄÓCMOTPETÜ B CBOEÉ BÛÇÛBAÞÙEÉ
//     ÏPOÃPAMME ÏAMßTÜ ÄËß CËEÄÓÞÙÈX MACCÈBOB, ÔÈÃÓPÈPÓÞÙÈX
//     B CÏÈCKE BÛÇOBA- Y(NEQN), WORK(3+6*NEQN), IWORK(5);
//     KPOME TOÃO, OH ÄOËÆEH OÁ'ßBÈTÜ F B OÏEPATOPE EXTERNAL,
//     ÏOÄÃOTOBÈTÜ ÏOÄÏPOÃPAMMÓ F(T,Y,YP) È ÏPÈCBOÈTÜ HA×AËÜ-
//     HÛE ÇHA×EHÈß ÏAPAMETPAM-

//     NEQN  -×ÈCËO ÈHTEÃPÈPÓEMÛX ÓPABHEHÈÉ (NEQN>=1)
//     Y(*)  -BEKTOP HA×AËÜHÛX ÓCËOBÈÉ
//     T     -HA×AËÜHAß TO×KA ÈHTEÃPÈPOBAHÈß,
//            T ÄOËÆHO ÁÛTÜ ÏEPEMEHHOÉ.
//     TOUT  -TO×KA BÛXOÄA,B KOTOPOÉ HÓÆHO HAÉTÈ ÇHA×EHÈE
//            PEØEHÈß. T=TOUT BOÇMOÆHO ËÈØÜ ÏPÈ ÏEPBOM
//            OÁPAÙEHÈÈ.B ÝTOM CËÓ×AE BÛXOÄ ÈÇ RKF45 ÏPOÈ-
//            CXOÄÈT CO ÇHA×EHÈEM ÏAPAMETPA IFLAG=2,ECËÈ
//            MOÆHO ÏPOÄOËÆATÜ ÈHTEÃPÈPOBAHÈE.
//     RELERR-ÃPAHÈÖA ÄËß OTHOCÈTEËÜHOÉ ËOKAËÜHÛÉ ÏOÃPEØHOCTEÈ.
//     ABSERR-ÃPAHÈÖA ÄËß AÁCOËÞTHOÉ    ËOKAËÜHÛÉ ÏOÃPEØHOCTEÈ.
//            ÝTÈ ÃPAHÈÖÛ ÄOËÆHÛ ÁÛTÜ HEOTPÈÖATEËÜHÛ.
//            RELERR ÄOËÆHA ÁÛTÜ ÏEPEMEHHOÉ,A ABSERR MOÆET
//            ÁÛTÜ È KOHCTAHTOÉ.ÏPOÃPAMME, BOOÁÙE ÃOBOPß
//            HE CËEÄÓET ÇAÄABATÜ ÃPAHÈÖÓ ÄËß OTHOCÈTEËÜHOÉ
//            OØÈÁKÈ,MEHÜØÓÞ, ×EM ÏPÈMEPHO 1.E-7. ÄAÁÛ ÈÇÁEÆATÜ
//            TPÓÄHOCTEÉ ,CBßÇAHHÛX C O×EHÜ BÛCOKÈMÈ ÇAÏPOCAMÈ
//            K TO×HOCTÈ, ÏPOÃPAMMA TPEÁÓET,×TOÁÛ RELERR
//            ÁÛËA ÁOËÜØE, ×EM HEKOTOPÛÉ ÏAPAMETP OTHOCÈTEËÜHOÉ
//            OØÈÁKÈ,BÛ×ÈCËßEMÛÉ BHÓTPÈ EE È ÇABÈCßÙÈÉ OT
//            MAØÈHÛ.B ×ACTHOCTÈ,HE PAÇPEØAETCß ÇAÄAHÈE TOËÜKO
//            AÁCOËÞTHOÉ OØÈÁKÈ.ECËÈ ÆE ÇAÄAHO ÇHA×EHÈE RELERR,
//            MEHÜØEE ÄOÏÓCTÈMOÃO, TO RKF45 ÓBEËÈ×ÈBAET RELERR
//            HAÄËEÆAÙÈM OÁPAÇOM È BOÇBPAÙAET ÓÏPABËEHÈE ÏOËÜ-
//            ÇOBATEËÞ, ÏPEÆÄE ×EM ÏPOÄOËÆATÜ ÈHTEÃPÈPOBAHÈE.
//     IFLAG-+1,-1.ÝTO ÓKAÇATEËÜ HACTPOÉKÈ ÏPOÃPAMMÛ ÄËß KAÆÄOÉ
//            HOBOÉ ÇAÄA×È. HOPMAËÜHOE BXOÄHOE ÇHA×EHÈE PABHO+1.
//            ÏOËÜÇOBATEËÜ ÄOËÆEH ÇAÄABATÜ IFLAG=-1 ËÈØÜ B TOM
//            CËÓ×AE,  KOÃÄA HEOÁXOÄÈMO ÓÏPABËEHÈE OÄHOØAÃOBÛM
//            ÈHTEÃPATOPOM.B ÝTOM CËÓ×AE RKF45 ÏÛTAETCß ÏPOÄOËÆÈTÜ
//            PEØEHÈE HA OÄÈH ØAÃ B HAÏPABËEHÈÈ TOUT ÏPÈ KAÆÄOM
//            O×EPEÄHOM BÛÇOBE. ÏOCKOËÜKÓ ÝTOT PEÆÈM PAÁOTÛ
//            BECÜMA HEÝKOHOMÈ×EH, EÃO CËEÄÓET ÏPÈMEHßTÜ
//            ËÈØÜ B CËÓ×AE KPAÉHEÉ HEOÁXOÄÈMOCTÈ.

//     ÈHÔOPMAÖÈß HA BÛXOÄE

//     Y(*)    -PEØEHÈE B TO×KE T
//     T       -ÏOCËEÄHßß TO×KA,ÄOCTÈÃHÓTAß ÏPÈ ÈHTEÃPÈPOBAHÈÈ.
//     IFLAG=2 -ÏPÈÈHTEÃPÈPOBAHÈÈ ÄOCTÈÃHÓTO TOUT.ÝTO ÇHA×EHÈE
//              ÏAPAMETPA ÓKAÇÛBAET HA ÓCÏEØHÛÉ BÛXOÄ È
//              ßBËßETCß HOPMAËÜHÛM PEÆÈMOM ÄËß ÏPOÄOËÆEHÈß
//              ÈHTEÃPÈPOBAHÈß.
//          =3 -ÈHTEÃPÈPOBAHÈE HE ÁÛËO ÇAKOH×EHO ÈÇ-ÇA TOÃO,
//              ×TO ÇAÄAHHOE ÇHA×EHÈE ÃPAHÈÖÛ ÄËß OTHOCÈTEËÜHOÉ
//              OØÈÁKÈ OKAÇAËOCÜ CËÈØKOM MAËO. ÄËß ÏPOÄOËÆEHÈß
//              ÈHTEÃPÈPOBAHÈß RELERR ÁÛËO HAÄËEÆAÙÈM OÁPAÇOM
//              ÓBEËÈ×EHO.
//          =4 -ÈHTEÃPÈPOBAHÈE HE ÁÛËO ÇAKOH×EHO ÈÇ-ÇA TOÃO,
//              ×TO ÏOTPEÁOBAËOCÜ ÁOËEE 3000 BÛ×ÈCËEHÈÉ ÏPO-
//              ÈÇBOÄHOÉ.ÝTO COOTBETCTBYET ÏPÈÁËÈÇÈTEËÜHO
//              500 ØAÃAM.
//          =5 -ÈHTEÃPÈPOBAHÈE HE ÁÛËO ÇAKOH×EHO ÈÇ-ÇA TOÃO,
//              ×TO PEØEHÈE OÁPATÈËOCÜ B HYËÜ,BCËEÄCTBÈE ×EÃO
//              TECT TOËÜKO OTHOCÈTEËÜHOÉ OØÈÁKÈ HE ÏPOXOÄÈT.
//              ÄËß ÏPOÄOËÆEHÈß HEOÁXOÄÈMO HEHYËEBOE ÇHA×EHÈE
//              ÏAPAMETPA ABSERR. ÈCÏOËÜÇOBAHÈE HA OÄÈH ØAÃ
//              PEÆÈMA ÏOØAÃOBOÃO ÈHTEÃPÈPOBAHÈß ßBËßETCß
//              PAÇYMHÛM BÛXOÄOM ÈÇ ÏOËOÆEHÈß.
//          =6 -ÈHTEÃPÈPOBAHÈE  HE ÁÛËO ÇAKOH×EHO ÈÇ-ÇA TOÃO,
//              ×TO TPEÁYEMAß TO×HOCTÜ HE MOÃËA ÁÛTÜ ÄOCTÈÃHÓTA
//              ÄAÆE ÏPÈ HAÈMEHÜØEÉ ÄOÏÓCTÈ MOÉ BEËÈ×ÈHE ØAÃA.
//              ÏOËÜÇOBATEËÜ ÄOËÆEH ÓBEËÈ×ÈTÜ ÃPAHÈÖÓ ÏOÃPEØ-
//              HOCTÈ,ÏPEÆÄE ×EM MOÆHO ÁÓÄET ÏOÏÛTATÜCß
//              ÏPOÄOËÆATÜ ÈHTEÃPÈPOBAHÈE.
//          =7 -ÏO BCEÉ BÈÄÈMOCTÈ, RKF45 HEÝÔÔEKTÈBHA ÏPÈ
//              PEØEHÈÈ ÝTOÉ ÇAÄA×È. CËÈØKOM ÁOËÜØOE ×ÈCËO
//              TPEÁÓEMÛX BÛXOÄHÛX TO×EK ÏPEÏßTCTBÓET BÛÁOPÓ
//              ECTECTBEHHOÉ BEËÈ×ÈHÛ ØAÃA.CËEÄÓET ÈCÏOËÜÇOBATÜ
//              PEÆÈM ÏOØAÃOBOÃO ÈHTEÃPÈPOBAHÈß.
//          =8 -HEÏPABÈËÜHOE ÇAÄAHÈE BXOÄHÛX ÏAPAMETPOB.ÝTO
//              ÇHA×EHÈE ÏOßBËßETCß,ECËÈ ÄOÏÓÙEHA OÄHA ÈÇ
//              CËEÄÓÞÙÈX OØÈÁOK-
//                              NEQN<=0
//                  T=TOUT  È  IFLAG!=+1  ÈËÈ -1
//                  RELERR  ÈËÈ  ABSERR<0
//                  IFLAG==0  ÈËÈ <-2  ÈËÈ >8
//     WORK(*) -ÈHÔOPMAÖÈß, KOTOPAß OÁÛ×HO HE ÏPEÄCTABËßET ÈHTE-
//              PECA ÄËß ÏOËÜÇOBATEËß, HO HEOÁXOÄÈMA ÏPÈ ÏOCËE-
//              ÄÓÞÙÈX BÛÇOBAX. WORK(1),...,WORK(NEQN) COÄEPÆAT
//              ÏEPBÛE ÏPOÈÇBOÄHÛE BEKTOPA PEØEHÈß Y B TO×KE T.
//              WORK(NEQN+1) XPAHÈT BEËÈ×ÈHÓ ØAÃA H,C KOTOPOÉ
//              MOÆHO ÏOÏÛTATÜCß ÏPOBECTÈ CËEÄÓÞÙÈÉ ØAÃ.
//     IWORK(*) -ÈHÔOPMAÖÈß, KOTOPAß OÁÛ×HO HE ÏPEÄCTABËßET ÈHTE-
//               PECA ÄËß ÏOËÜÇOBATEËß, HO HEOÁXOÄÈMA ÏPÈ ÏOCËE-
//               ÄÓÞÙÈX BÛÇOBAX. B IWORK(1) COÄEPÆÈTCß
//               C×ET×ÈK ×ÈCËA BÛ×ÈCËEHÈÉ ÏPOÈÇBOÄHÛX.

//     ÏOCËEÄÓÞÙÈE OÁPAÙEHÈß K RKF45

//          HA BÛXOÄE ÏOÄÏPOÃPAMMÛ RKF45 ÈMEETCß BCß ÈHÔOPMAÖÈß,
//     HEOÁXOÄÈMAß  ÄËß ÏPOÄOËÆEHÈß ÈHTEÃPÈPOBAHÈß.ECËÈ ÏPÈ
//     ÈHTEÃPÈPOBAHÈÈ ÄOCTÈÃHÓTO TOUT,TO ÏOËÜÇOBATEËÞ ÄOCTA-
//     TO×HO OÏPEÄEËÈTÜ HOBOE ÇHA×EHÈE  TOUT È CHOBA OÁPATÈTÜ-
//     Cß K RKF45.
//          B PEÆÈME ÏOØAÃOBOÃO ÈHTEÃPÈPOBAHÈß (IFLAG=-2)
//     ÏOËÜÇOBATEËÜ ÄOËÆEH ÈMETÜ B BÈÄÓ,×TO KAÆÄÛÉ ØAÃ
//     BÛÏOËHßETCß B HAÏPABËEHÈÈ TEKÓÙEÃO ÇHA×EHÈß TOUT
//     (CÈÃHAËÈÇÈPÓEMOM  ÈÇMEHEHÈEM IFLAG HA 2). ÏOËÜÇOBATEËÜ
//     ÄOËÆEH ÇAÄATÜ HOBOE ÇHA×EHÈE TOUT È ÏEPEOÏPEÄEËÈTÜ
//     IFLAG HA -2, ×TOÁÛ ÏPOÄOËÆATÜ B PEÆÈME ÏOØAÃOBOÃO
//     ÈHTEÃPÈPOBAHÈß.
//          ECËÈ ÈHTEÃPÈPOBAHÈE HE ÁÛËO ÇAKOH×EHO,HO
//     ÏOËÜÇOBATEËÜ XO×ET ÏPOÄOËÆATÜ (CËÓ×AÈ IFLAG=3,4), OH
//     ÏOÏPOCTÓ CHOBA OÁPAÙAETCß K RKF45.ÏPÈ IFLAG=3 ÏAPA-
//     METP RELERR ÁÛË  ÈÇMEHEH HAÄËEÆAÙÈM ÄËß ÏPOÄOËÆEHÈß
//     ÈHTEÃPÈPOBAHÈß OÁPAÇOM.B CËÓ×AE IFLAG=4 C×ET×ÈK
//     ×ÈCËA ÇHA×EHÈÉ ÔÓHKÖÈÈ ÁÓÄET ÏEPEOÏPEÄEËEH HA 0, È
//     ÁÓÄÓT PAÇPEØEHÛ EÙE 3000 BÛ×ÈCËEHÈÉ ÔÓHKÖÈÈ.
//          OÄHAKO B CËÓ×AE IFLAG=5, ÏPEÆÄE ×EM MOÆHO ÁÓÄET
//     ÏPOÄOËÆATÜ ÈHTEÃPÈPOBAHÈE,ÏOËÜÇOBATEËÜ ÄOËÆEH CHA×AËA
//     ÈÇMEHÈTÜ KPÈTEPÈÉ OØÈÁKÈ, ÇAÄAB ÏOËOÆÈTEËÜHOE ÇHA×EHÈE
//     ÄËß ABSERR. ECËÈ OH HE CÄEËAET ÝTOÃO, BÛÏOËHEHÈE ÏPO-
//     ÃPAMMÛ ÁÓÄET ÏPEKPAÙEHO.
//        TO×HO TAK ÆE,B CËÓ×AE IFLAG=6,ÏPEÆÄE ×EM ÏPOÄOË-
//     ÆATÜ ÈHTEÃPÈPOBAHÈE,ÏOËÜÇOBATEËÞ HEOÁXOÄÈMO ÏEPEOÏPE-
//     ÄEËÈTÜ IFLAG HA 2 (ÈËÈ -2, ECËÈ ÈCÏOËÜÇÓETCß PEÆÈM
//     ÏOØAÃOBOÃO ÈHTEÃPÈPOBAHÈß) È ÓBEËÈ×ÈTÜ ÇHA×EHÈE ÄËß
//     ABSERR ËÈÁO RELERR,ËÈÁO È ÄËß TOÃO,È ÄËß ÄPÓÃOÃO.
//     ECËÈ ÝTO HE ÁÓÄET CÄEËAHO,BÛÏOËHEHÈE ÏPOÃPAMMÛ
//     ÏPEKPAÙAETCß. ÏOßBËEHÈE IFLAG=6 ÓKAÇÛBAET HA HEPEÃÓ-
//     ËßPHOCTÜ (PEØEHÈE ÁÛCTPO MEHßETCß ÈËÈ, BOÇMOÆHO,
//     ÈMEETCß OCOÁEHHOCTÜ),È ×ACTO B ÏOÄOÁHÛX CËÓ×AßX
//     HE ÈMEET CMÛCËA ÏPOÄOËÆATÜ ÈHTEÃPÈPOBAHÈE.
//          ECËÈ ÁÓÄET ÏOËÓ×EHO ÇHA×EHÈE IFLAG=7,TO ÏOËÜÇO-
//     BATEËÜ ÄOËÆEH ÏEPEÉTÈ K PEÆÈMÓ ÏOØAÃOBOÃO ÈHTEÃPÈPO-
//     BAHÈß C BEËÈ×ÈHOÉ ØAÃA,OÏPEÄEËßEMOÉ ÏPOÃPAMMOÉ, ÈËÈ
//     PACCMOTPETÜ BOÇMOÆHOCTTÜ ÏEPEXOÄA HA ÏPOÃPAMMÛ METOÄOB
//     AÄAMCA.ECËÈ BCE ÆE ÏOËÜÇOBATEËÜ XO×ET ÏPOÄOËÆATÜ
//     ÈHTEÃPÈPOBAHÈE ÏO ÏOÄÏPOÃPAMME RKF45,OH ÄOËÆEH ÄO HOBOÃO
//     OÁPAÙEHÈß K HEÉ ÏEPEOÏPEÄEËÈTÜ IFLAG HA 2.B ÏPOTÈBHOM
//     CËÓ×AE BÛÏOËHEHÈE ÏPOÃPAMMÛ ÁÓÄET ÏPEKPAÙEHO.
//          ECËÈ ÏOËÓ×EHO ÇHA×EHÈE IFLAG=8,TO ÈHTEÃPÈPOBAHÈE
//     HEËÜÇß ÏPOÄOËÆATÜ,ÏOKA HE ÁÓÄÓT ÈCÏPABËEHÛ OØÈÁO×HÛE
//     BXOÄHÛE ÏAPAMETPÛ. HÓÆHO OTMETÈTÜ,×TO MACCÈBÛ WORK È
//     IWORK COÄEPÆAT ÈHÔOPMAÖÈÞ,HEOÁXOÄÈMÓÞ ÄËß ÄAËÜHEÉØEÃO
//     ÈHTEÃPÈPOBAHÈß.ÏOÝTOMÓ B ÝTÈ MACCÈBÛ HEËÜÇß BHOCÈTÜ
//     ÈÇMEHEHÈÉ.

void
FEHL(void(F)(REAL T, REAL *Y, REAL *YP), int NEQN, const REAL *Y, REAL &T, REAL H, const REAL *YP, REAL *F1, REAL *F2, REAL *F3,
     REAL *F4, REAL *F5, REAL *S) {

//     METOÄ PÓHÃE-KÓTTA-ÔEËÜÁEPÃA ×ETBEPTOÃO-ÏßTOÃO ÏOPßÄKA

//     ÏOÄÏPOÃPAMMA FEHL ÈHTEÃPÈPÓET CÈCTEMÓ ÈÇ NEQN
//     OÁÛKHOBEHHÛX ÄÈÔÔEPEHÖÈAËÜHÛX ÓPABHEHÈÉ ÏEPBOÃO
//     ÏOPßÄKA CËEÄÓÞÙEÃO BÈÄA

//          DY(I)/DT=F(T,Y(1),...Y(NEQN)),

//     ÃÄE HA×AËÜHÛE ÇHA×EHÈß Y(I) È HA×AËÜHÛE ÏPOÈÇBOÄHÛE
//     YP(I) ÇAÄAHÛ B HA×AËÜHOÉ TO×KE T. FEHL ÏPOÄOËÆAET
//     PEØEHÈE HA ÔÈKCÈPOBAHHÛÉ ØAÃ H È ÏOMEÙAET B MACCÈB
//     S(I) ÏPÈÁËÈÆEHÈE K PEØEHÈÞ B TO×KE T+H, ÈMEÞÙEE ÏßTÛÉ
//     ÏOPßÄOK TO×HOCTÈ (ËOKAËÜHÛÉ ÏOPßÄOK PABEH ØECTÈ).
//     F1,...F5-MACCÈBÛ PAÇMEPHOCTÈ NEQN,HEOÁXOÄÈMÛE BHÓTPÈ
//     ÏPOÃPAMMÛ.
//          B ÔOPMÓËAX ÏPOÈÇBEÄEHA ÃPÓÏÏÈPOBKA C ÖEËÜÞ
//     ÓMEHÜØÈTÜ ÏOTEPÞ BEPHÛX ÇHAKOB.
//          ×TOÁÛ MOÆHO ÁÛËO PAÇËÈ×ATÜ PAÇHÛE HEÇABÈCÈMÛE
//     APÃÓMEHTÛ, ÏPÈ OÁPAÙEHÈÈ K FEHL HE CËEÄÓET ÇAÄABATÜ
//     ÄËß H ÇHA×EHÈE,MEHÜØEE ÓMHOÆEHHHOÉ HA 13 OØÈÁKÈ
//     OKPÓÃËEHÈß B T.


    REAL CH;
    int K;

    CH = H / 4.0;
    for (K = 1; K <= NEQN; K++)F5[K - 1] = Y[K - 1] + CH * YP[K - 1];

    F(T + CH, F5, F1);
    CH = 3.0 * H / 32.0;
    for (K = 1; K <= NEQN; K++)F5[K - 1] = Y[K - 1] + CH * (YP[K - 1] + 3.0 * F1[K - 1]);
    F(T + 3.0 * H / 8.0, F5, F2);

    CH = H / 2197.0;
    for (K = 1; K <= NEQN; K++)
        F5[K - 1] = Y[K - 1] + CH * (1932.0 * YP[K - 1] + (7296.0 * F2[K - 1] - 7200.0 * F1[K - 1]));
    F(T + 12.0 * H / 13.0, F5, F3);

    CH = H / 4104.0;
    for (K = 1; K <= NEQN; K++)
        F5[K - 1] = Y[K - 1] + CH * ((8341.0 * YP[K - 1] - 845.0 * F3[K - 1]) +
                                     (29440.0 * F2[K - 1] - 32832.0 * F1[K - 1]));
    F(T + H, F5, F4);

    CH = H / 20520.0;
    for (K = 1; K <= NEQN; K++)
        F1[K - 1] = Y[K - 1] + CH * ((-6080.0 * YP[K - 1] + (9295.0 * F3[K - 1] - 5643.0 * F4[K - 1])) +
                                     (41040.0 * F1[K - 1] - 28352.0 * F2[K - 1]));
    F(T + H / 2.0, F1, F5);

//     BÛ×ÈCËÈTÜ ÏPÈÁËÈÆEHHOE PEØEHÈE B TO×KE T+H

    CH = H / 7618050.0;
    for (K = 1; K <= NEQN; K++)
        S[K - 1] = Y[K - 1] + CH * ((902880.0 * YP[K - 1] + (3855735.0 * F3[K - 1] - 1371249.0 * F4[K - 1])) +
                                    (3953664.0 * F2[K - 1] + 277020.0 * F5[K - 1]));
}                            // END


void
RKFS(void(F)(REAL T, REAL *Y, REAL *YP), int NEQN, REAL *Y, REAL &T, REAL TOUT, REAL &RELERR, REAL &ABSERR, int &IFLAG,
     REAL *YP, REAL H, REAL *F1, REAL *F2, REAL *F3, REAL *F4, REAL *F5, REAL SAVRE, REAL SAVAE) {

//     METOÄ PÓHÃE-KÓTTA-ÔEËÜÁEPÃA ×ETBEPTOÃO-ÏßTOÃO ÏOPßÄKA
//     RKFS ÈHTEÃPÈPÓET CÈCTEMÓ OÁÛKHOBEHHÛX ÄÈÔÔE-
//     PEHÖÈAËÜHÛX ÓPABHEHÈÉ ÏEPBOÃO ÏOPßÄKA(CM. KOM-
//     MEHTAPÈÉ K RKF45). MACCÈBÛ YP,F1,F2,F3,F4 È F5
//     (PAÇMEPHOCTÈ ÏO KPAÉHEÉ MEPE NEQN) È ÏEPEMEH-
//     HÛE H,SAVRE,SAVAE,NFE,KOP,INIT,JFLAG È KFLAG
//     ÈCÏOËÜÇÓÞTCß BHÓTPÈ ÏPOÃPAMMÛ È BÛHECEHÛ B
//     CÏÈCOK BÛÇOBA,×TOÁÛ COXPAHÈTÜ ÈX OÏPEÄEËEH-
//     HOCTÜ ÏPÈ ÏOBTOPHOM OÁPAÙEHÈÈ.ÏOÝTOMÓ ÈX ÇHA-
//     ×EHÈß HE ÄOËÆHÛ  ÈÇMEHßTÜCß  ÏOËÜÇOBATEËEM.
//     BOÇMOÆHÛÉ ÈHTEPEC ÏPEÄCTABËßÞT ÏAPAMETPÛ
//     YP  -ÏPOÈÇBOÄHAß BEKTOPA PEØEHÈß B TO×KE T
//     H   -ÏPEÄÏOËAÃAEMÛÉ PAÇMEP ØAÃA ÄËß O×EPEÄHOÃO ÝTAÏA
//     NFE -C×ET×ÈK ×ÈCËA BÛ×ÈCËEHÈÉ ÔÓHKÖÈÈ

    bool HFAILD, OUTPUT;
    static int NFE,KOP, INIT, JFLAG, KFLAG;
    REAL A, AE, DT, EE, EEOET, ESTTOL, ET, HMIN, RER, S, SCALE, TOL, TOLN, U26, EPSP1, EPS, YPK;
    int K, MFLAG;


//     REMIN-ÝTO MÈHÈMAËÜHOE ÄOÏÓCTÈMOE ÇHA×EHÈE ÄËß
//     RELERR.ÏOÏÛTKÈ ÏOËÓ×ÈTÜ ÏO ÝTOÉ ÏOÄÏPOÃPAMME
//     ÁOËEE BÛCOKÓÞ TO×HOCTÜ OÁÛ×HO CTOßT O×EHÜ
//     ÄOPOÃO È ÇA×ACTÓÞ ÁEÇÓCÏEØHÛ.

    REAL REMIN = 1.E-12;

//     CTOÈMOCTÜ C×ETA KOHTPOËÈPÓETCß TPEÁOBAHÈEM,
//     ×TOÁÛ KOËÈ×ECTBO BÛ×ÈCËEHÈÉ ÔÓHKÖÈÈ ÁÛËO OÃ-
//     PAHÈ×EHO BEËÈ×ÈHOÉ,ÏPÈÁËÈÇÈTEËÜHO PABHOÉ ÇHA-
//     ×EHÈÞ ÏAPAMETPA MAXNFE.ÏPÈHßTOE ÇÄECÜ ÇHA×E-
//     HÈE ÏPÈMEPHO COOTBETCTBÓET 500 ØAÃAM.

    int MAXNFE = 3000;

//     ÏPOBEPÈTÜ BXOÄHÛE ÏAPAMETPÛ

    if (NEQN < 1)goto _10;
    if ((RELERR < 0.0) || (ABSERR < 0.0))goto _10;
    MFLAG = ABS(IFLAG);
    if ((MFLAG == 0) || (MFLAG > 8))goto _10;
    if (MFLAG != 1)goto _20;

//     ÏEPBÛÉ BÛÇOB,BÛ×ÈCËÈTÜ MAØÈHHOE ÝÏCÈËOH

    EPS = 1.0;
    _5:;
    EPS = EPS / 2.0;
    EPSP1 = EPS + 1.;
    if (EPSP1 > 1.)goto _5;
    U26 = 26. * EPS;
    goto _50;

//     OØÈÁKA BXOÄHOÉ ÈHÔOPMAÖÈÈ

    _10:;
    IFLAG = 8;
    return;  //?

//     ÏPOBEPÈTÜ BOÇMOÆHOCTÈ ÏPOÄOËÆEHÈß

    _20:;
    if ((T == TOUT) && (KFLAG != 3))goto _10;
    if (MFLAG != 2)goto _25;

//     IFLAG=+2 ÈËÈ -2

    if ((KFLAG == 3) || (INIT == 0))goto _45;
    if (KFLAG == 4)goto _40;
    if ((KFLAG == 5) && (ABSERR == 0.0))goto _30;
    if ((KFLAG == 6) && (RELERR <= SAVRE) && (ABSERR <= SAVAE))goto _30;
    goto _50;

//     IFLAG=3,4,5,6,7 ÈËÈ 8

    _25:;
    if (IFLAG == 3)goto _45;
    if (IFLAG == 4)goto _40;
    if ((IFLAG == 5) && (ABSERR > 0.0))goto _45;

//     ÈHTEÃPÈPOBAHÈE HEËÜÇß ÏPOÄOËÆATÜ,ÏOCKOËÜKÓ ÏOËÜ-
//     ÇOBATEËÜ HE BÛÏOËHÈË ÈHCTPÓKÖÈÉ,COOTBETCTBÓÞÙÈX
//     ÇHA×EHÈßM IFLAG=5,6,7 ÈËÈ 8

    _30:;
//printf("ÈHTEÃPÈPOBAHÈE ÏPEPBAHO, ÏOCKOËÜKÓ ÏOËÜÇOBATEËÜ"
//"HE BÛÏOËHÈË ÈHCTPÓKÖÈÉ RKF45, COOTBETCTBÓÞÙÈX"
//"ÇHA×EHÈßM IFLAG=5,6,7 ÈËÈ 8");
    exit(0);

//     ÏEPEOÏPEÄEËÈTÜ C×ET×ÈK ×ÈCËA BÛ×ÈCËEHÈÉ ÔÓHKÖÈÈ

    _40:;
    NFE = 0;
    if (MFLAG == 2)goto _50;

//     ÏEPEOÏPEÄEËÈTÜ ÇHA×EHÈE FLAG,ÓCTAHOBËEHHOE
//     ÏPÈ ÏPEÄÛÄÓÙEM OÁPAÙEHÈÈ

    _45:;
    IFLAG = JFLAG;
    if (KFLAG == 3)MFLAG = ABS(IFLAG);

//     COXPAHÈTÜ BXOÄHOE ÇHA×EHÈE IFLAG È ÓCTAHOBÈTÜ
//     ÇHA×EHÈE FLAG, COOTBETCTBÓÞÙEE ÏPOÄOËÆEHÈÞ,
//     ÄËß ÁÓÄÓÙEÉ ÏPOBEPKÈ

    _50:;
    JFLAG = IFLAG;
    KFLAG = 0;

//     COXPAHÈTÜ ÇHA×EHÈß RELERR È ABSERR ÄËß BXOÄHOÉ
//     ÏPOBEPKÈ ÏPÈ ÏOCËEÄÓÞÙÈX OÁPAÙEHÈßX

    SAVRE = RELERR;
    SAVAE = ABSERR;

//     ÓCTAHOBÈTÜ ÇHA×EHÈE ÃPAHÈÖÛ ÄËß OTHOCÈTEËÜ-
//     HOÉ ÏOÃPEØHOCTÈ,PABHOE KAK MÈHÈMÓM 2*EPS+
//     REMIN,×TOÁÛ ÈÇÁEÆATÜ TPÓÄHOCTEÉ,CBßÇAHHÛX
//     C TPEÁOBAHÈEM HEÄOCTÈÆÈMOÉ TO×HOCTÈ

    RER = 2. * EPS + REMIN;
    if (RELERR >= RER)goto _55;

//     ÇAÄAHHAß ÃPAHÈÖA OTHOCÈTEËÜHOÉ ÏOÃPEØHOCTÈ
//     CËÈØKOM MAËA

    RELERR = RER;
    IFLAG = 3;
    KFLAG = 3;
    return;  //?

    _55:;
    DT = TOUT - T;

    if (MFLAG == 1)goto _60;
    if (INIT == 0)goto _65;
    goto _80;

//     ÏPÈCBOEHÈE HA×AËÜHÛX ÇHA×EHÈÉ (ÈHÈÖÈÈPOBA-
//             HÈE)-ÓCTAHOBÈTÜ ÇHA×EHÈE ÓKAÇATEËß
//             OKOH×AHÈß ÈHÈÖÈÈPOBAHÈß,INIT
//             ÓCTAHOBÈTÜ ÇHA×EHÈE ÓKAÇATEËß CËÈØ-
//             KOM ÁOËÜØOÃO ÇATPEÁOBAHHOÃO ×ÈCËA BÛ-
//             XOÄHÛX TO×EK,KOP
//             BÛ×ÈCËÈTÜ HA×AËÜHÛE ÏPOÈÇBOÄHÛE
//             ÓCTAHOBÈTÜ ÇHA×EHÈE C×ET×ÈKA ×ÈCËA
//             BÛ×ÈCËEHÈÉ ÔÓHKÖÈÈ,NFE
//             OÖEHÈTÜ HA×EËÜHÓÞ BEËÈ×ÈHÓ ØAÃA

    _60:;
    INIT = 0;
    KOP = 0;

    A = T;
    F(A, Y, YP);
    NFE = 1;
    if (T != TOUT)goto _65;
    IFLAG = 2;
    return;  //?

    _65:;
    INIT = 1;
    H = ABS(DT);
    TOLN = 0.;
    for (K = 1; K <= NEQN; K++) {   //? target=70
        TOL = RELERR * ABS(Y[K - 1]) + ABSERR;
        if (TOL <= 0)goto _70;
        TOLN = TOL;
        YPK = ABS(YP[K - 1]);
        if (YPK * pow(H, 5) > TOL)H = pow((TOL / YPK), 0.2);
        _70:;
    }                            // CONTINUE
    if (TOLN <= 0.0)H = 0.0;
    H = MAX(H, U26 * MAX(ABS(T), ABS(DT)));
    JFLAG = SIGN(2, IFLAG);

//     ÏPÈCBOÈTÜ BEËÈ×ÈHE ØAÃA ÇHAK,COOTBETCTBÓÞÙÈÉ
//     ÈHTEÃPÈPOBAHÈÞ B HAÏPABËEHÈÈ OT T K TOUT

    _80:;
    H = SIGN(H, DT);

//     ÏPOBEPKA, HACKOËÜKO CEPÜEÇHO BËÈßHÈE HA RKF45
//     CËÈØKOM ÁOËÜØOÃO ÇATPEÁOBAHHOÃO ×ÈCËA BÛXOÄ-
//     HÛX TO×EK

    if (ABS(H) >= 2.0 * ABS(DT))KOP = KOP + 1;
    if (KOP != 100)goto _85;

//     ×PEÇMEPHAß ×ACTOTA BÛXOÄOB

    KOP = 0;
    IFLAG = 7;
    return;

    _85:;
    if (ABS(DT) > U26 * ABS(T))goto _95;

//     ECËÈ O×EHÜ ÁËÈÇKO K TO×KE BÛXOÄA,ÏPOÝKCTPAÏO-
//     ËÈPOBATÜ È BEPHÓTÜCß ÏO MECTÓ BÛÇOBA

    for (K = 1; K <= NEQN; K++)Y[K - 1] = Y[K - 1] + DT * YP[K - 1];
    A = TOUT;
    F(A, Y, YP);
    NFE = NFE + 1;
    goto _300;

//     ÏPÈCBOÈTÜ HA×AËÜHOE ÇHA×EHÈE ÈHÄÈKATOPÓ TO×KÈ
//     BÛXOÄA

    _95:;
    OUTPUT = false;

//     ×TOÁÛ ÈÇÁEÆATÜ HEOÏPABÄAHHOÃO MAØÈHHOÃO HÓËß
//     ÏPÈ BÛ×ÈCËEHÈÈ ÔÓHKÖÈÈ OT ÃPAHÈÖ ÏOÃPEØHO-
//     CTEÉ,ÏPOMACØTAÁÈPOBATÜ ÝTÈ ÃPAHÈÖÛ

    SCALE = 2. / RELERR;
    AE = SCALE * ABSERR;

//     ÏOØAÃOBOE ÈHTEÃPÈPOBAHÈE

    _100:;
    HFAILD = false;

//     ÓCTAHOBÈTÜ HAÈMEHÜØÓÞ ÄOÏÓCTÈMÓÞ BEËÈ×ÈHÓ ØAÃA

    HMIN = U26 * ABS(T);

//     ÈCÏPABÈTÜ ÏPÈ HEOÁXOÄÈMOCTÈ BEËÈ×ÈHÓ ØAÃA,
//     ×TOÁÛ ÄOCTÈÃHÓTÜ TO×KÈ BÛXOÄA. PACC×ÈTATÜ HA
//     ÄBA ØAÃA BÏEPEÄ,×TOÁÛ ÈÇÁEÆATÜ CËÈØKOM PEÇKÈX
//     ÈÇMEHEHÈÉ B BEËÈ×ÈHE ØAÃA È TEM CAMÛM ÓMEHÜ-
//     ØÈTÜ BËÈßHÈE BÛXOÄHÛX TO×EK HA ÏPOÃPAMMÓ.

    DT = TOUT - T;
    if (ABS(DT) >= 2. * ABS(H))goto _200;
    if (ABS(DT) > ABS(H))goto _150;

//     CËEÄÓÞÙÈÉ ÓCÏEØHÛÉ ØAÃ ÇABEPØÈT ÈHTEÃPÈPO-
//     BAHÈE ÄO ÓKAÇAHHOÉ TO×KÈ BÛXOÄA

    OUTPUT = true;
    H = DT;
    goto _200;

    _150:;
    H = 0.5 * DT;



//     BHÓTPEHHÈÉ OÄHOØAÃOBÛÉ ÈHTEÃPATOP

//     ÃPAHÈÖÛ ÏOÃPEØHOCTEÉ ÁÛËÈ ÏPOMACØTAÁÈPOBAHÛ,
//     ×TOÁÛ ÈÇÁEÆATÜ HEOÏPABÄAHHOÃO MAØÈHHOÃO HÓËß
//     ÏPÈ BÛ×ÈCËEHÈÈ ÔÓHKÖÈÈ OT HÈX.
//     ×TOÁÛ ÈÇÁEÆATÜ OÁPAÙEHÈß B HÓËÜ ÇHAMEHATEËß
//     B TECTE,OTHOCÈTEËÜHAß OØÈÁKA ÈÇMEPßETCß  ÏO
//     OTHOØEHÈÞ K CPEÄHEMÓ  ÈÇ BEËÈ×ÈH PEØEHÈß
//     B HA×AËE È KOHÖE ØAÃA.
//     B ÔOPMÓËE,OÖEHÈBAÞÙEÉ OØÈÁKÓ,ÏPOÈÇBEÄEHA
//     ÃPÓÏÏÈPOBKA CËAÃAEMÛX,ÓMEHÜØAÞÙAß ÏOTEPÞ
//     BEPHÛX ÇHAKOB.
//     ×TOÁÛ PAÇËÈ×ATÜ MEÆÄÓ COÁOÉ PAÇHÛE APÃÓMEHTÛ,
//     ÄËß H HE ÄOÏÓCKAÞTCß ÇHA×EHÈß,MEHÜØÈE ÓMHO-
//     ÆEHHOÉ HA 26 OØÈÁKÈ OKPÓÃËEHÈß B T.
//     BBEÄEHÛ ÏPAKTÈ×ECKÈE OÃPAHÈ×EHÈß HA CKOPOCTÜ
//     ÈÇMEHEHÈß BEËÈ×ÈHÛ ØAÃA,×TOÁÛ CÃËAÄÈTÜ ÏPO-
//     ÖECC BÛÁOPA ÝTOÉ BEËÈ×ÈHÛ È ÈÇÁEÆATÜ ×PEÇMEP-
//     HOÃO EE PAÇÁPOCA B ÇAÄA×AX C HAPÓØEHÈEM HEÏPE-
//     PÛBHOCTÈ.
//     ÈÇ ÏPEÄOCTOPOÆHOCTÈ ÏPOÃPAMMA ÁEPET 9/10 OT TOÉ
//     BEËÈ×ÈHÛ ØAÃA,KOTOPAß HÓÆHA ÏO EE OÖEHKE.
//     ECËÈHA ÄAHHOM ØAÃE ÁÛËA HEÓÄA×HAß ÏOÏÛTKA
//     TO ÏPÈ ÏËAHÈPOBAHÈÈ CËEÄÓÞÙEÃO ÓBEËÈ×EHÈE
//     ÄËÈHÛ ØAÃA HE ÄOÏÓCKAETCß. ÝTO ÏOBÛØAET ÝÔÔEK-
//     TÈBHOCTÜ ÏPOÃPAMMÛ ÄËß ÇAÄA× C PAÇPÛBAMÈ È
//     B OÁÙEM CËÓ×AE,ÏOCKOËÜKÓ ÈCÏOËÜÇÓETCß ËOKAËÜ-
//     HAß ÝKCTPAÏOËßÖÈß È ÄOÏOËHÈTEËÜHAß ÏPEÄOCTO-
//     POÆHOCTÜ KAÆETCß OÏPABÄAHHOÉ.


//     ÏPOBEPÈTÜ ×ÈCËO BÛ×ÈCËEHÈÉ ÏPOÈÇBOÄHÛX.ECËÈ
//     OHO HE ÏPEBÛØAET ÓCTAHOBËEHHOÃO ÏPEÄEËA,ÏO-
//     ÏPOÁOBATÜ ÏPOÄOËÆÈTÜ ÈHTEÃPÈPOBAHÈE C T ÄO T+H

    _200:;
    if (NFE <= MAXNFE)goto _220;

//     CËÈØKOM ÁOËÜØAß PAÁOTA

    IFLAG = 4;
    KFLAG = 4;
    return;

//     ÏPOÄOËÆÈTÜ ÏPÈÁËÈÆEHHOE PEØEHÈE HA OÄÈH ØAÃ ÄËÈHÛ H

    _220:;
    FEHL(F, NEQN, Y, T, H, YP, F1, F2, F3, F4, F5, F1);
    NFE = NFE + 5;

//     BÛ×ÈCËÈTÜ È CPABHÈTÜ ÄOÏÓCTÈMÛE ÃPAHÈÖÛ È
//     OÖEHKÈ ËOKAËÜHOÉ OØÈÁ,A ÇATEM CHßTÜ MACØTA-
//     ÁÈPOBAHÈE ÃPAHÈÖ.ÇAMETÜTE,×TO OTHOCÈTEËÜHAß
//     OØÈÁKA ÈÇMEPßETCß ÏO OTHOØEHÈÞ K CPEÄHEMÓ ÈÇ
//     BEËÈ×ÈH PEØEHÈß B HA×AËE È KOHÖE ØAÃA.

    EEOET = 0.;
    for (K = 1; K <= NEQN; K++) {   //? target=250
        ET = ABS(Y[K - 1]) + ABS(F1[K - 1]) + AE;
        if (ET > 0.)goto _240;

//     HEÏPABÈËÜHAß ÃPAHÈÖA ÏOÃPEØHOCTÈ

        IFLAG = 5;
        KFLAG = 5;
        return;

        _240:;
        EE = ABS((-2090. * YP[K - 1] + (21970. * F3[K - 1] - 15048. * F4[K - 1])) +
                 (22528. * F2[K - 1] - 27360. * F5[K - 1]));
//_250:;
        EEOET = MAX(EEOET, EE / ET);
    }

    ESTTOL = ABS(H) * EEOET * SCALE / 752400.;

    if (ESTTOL <= 1.0)goto _260;


//     HEÓÄA×HÛÉ ØAÃ
//            ÓMEHÜØÈTÜ BEËÈ×ÈHÓ ØAÃA È CHOBA ÏO-
//            ÏPOÁOBATÜ
//            ÓMEHÜØEHÈE OÃPAHÈ×ÈBAETCß CHÈÇÓ MHO-
//            ÆÈTEËEM 1/10

    HFAILD = true;
    OUTPUT = false;
    S = 0.1;
    if (ESTTOL < 59049.)S = 0.9 / pow(ESTTOL, 0.2);
    H = S * H;
    if (ABS(H) > HMIN)goto _200;

//     ÇAÄAHHAß ÃPAHÈÖA OØÈÁKÈ HEÄOCTÈÆÈMA ÄAÆE ÏPÈ
//     HAÈMEHÜØEÉ ÄOÏÓCTÈMOÉ BEËÈ×ÈHE ØAÃA

    IFLAG = 6;
    KFLAG = 6;
    return;  //?


//     ÓCÏEØHÛÉ ØAÃ
//            ÏOMECTÈTÜ B MACCÈB Y PEØEHÈE B TO×KE
//            T+H È BÛ×ÈCËÈTÜ ÏPOÈÇBOÄHÛE B ÝTOÉ
//            TO×KE

    _260:;
    T = T + H;
    for (K = 1; K <= NEQN; K++)Y[K - 1] = F1[K - 1];
    A = T;
    F(A, Y, YP);
    NFE = NFE + 1;


//     BÛÁPATÜ BEËÈ×ÈHÓ CËEÄÓÞÙEÃO ØAÃA
//     ÓBEËÈ×EHÈE OÃPAHÈ×EHO MHOÆÈTEËEM 5
//     ECËÈ HA ÄAHHOM ØAÃE ÁÛËA HEÓÄA×HAß
//     ÏOÏÛTKA,TO ÄËß CËEÄÓÞÙEÃO HE ÄOÏÓ-
//     CKAETCß BÛÁOP ÁOËÜØEÉ BEËÈ×ÈHÛ ØAÃA

    S = 5.;
    if (ESTTOL > 1.889568E-4)S = 0.9 / pow(ESTTOL, 0.2);
    if (HFAILD)S = MIN(S, 1.0);
    H = SIGN(MAX(S * ABS(H), HMIN), H);

//     KOHEÖ OÄHOØAÃOBOÃO ÈHTEÃPATOPA


//     HÓÆHO ËÈ ÄEËATÜ O×EPEÄHOÉ ØAÃ

    if (OUTPUT)goto _300;
    if (IFLAG > 0)goto _100;


//     ÈHTEÃPÈPOBAHÈE ÓCÏEØHO ÇABEPØEHO

//     PEÆÈM OÄHOØAÃOBOÃO ÈHTEÃPÈPOBAHÈß

    IFLAG = -2;
    return;

//     PEÆÈM ÈHTEÃPÈPOBAHÈß HA ÈHTEPBAËE

    _300:;
    T = TOUT;
    IFLAG = 2;
}                            // END

void
RKF45(void(F)(REAL T, REAL *Y, REAL *YP), int NEQN, REAL *Y, REAL &T, REAL TOUT, REAL &RELERR, REAL &ABSERR, REAL *WORK,
      int &IFLAG) {
//     ECËÈ TPAHCËßTOP ÏPOBEPßET ÈHÄEKCÛ, TO ÇAMEHÈTÜ
//     WORK[0] HA WORK[2+6*NEQN]

    int K1, K2, K3, K4, K5, K6, K1M;

//     BÛ×ÈCËÈTÜ ÈHÄEKCÛ ÄËß PACÙEÏËEHÈß PAÁO×EÃO MACCÈBA

    K1M = NEQN;
    K1 = K1M + 1;
    K2 = K1 + NEQN;
    K3 = K2 + NEQN;
    K4 = K3 + NEQN;
    K5 = K4 + NEQN;
    K6 = K5 + NEQN;

//     ÝTA ÏPOMEÆYTO×HAß ÏPOÃPAMMA ÏPOCTO COKPAÙAET ÄËß
//     ÏOËÜÇOBATEËß ÄËÈHHÛÉ CÏÈCOK BÛÇOBA ÏYTEM PACÙEÏËEHÈß
//     ÄBYX PAÁO×ÈX MACCÈBOB. ECËÈ ÝTO HE COBMECTÈMO C
//     TPAHCËßTOPOM,TO OH ÄOËÆEH OÁPAÙATÜCß HEÏOCPEÄCTBEHHO
//     K ÏOÄÏPOÃPAMME RKFS .

    RKFS(F, NEQN, Y, T, TOUT, RELERR, ABSERR, IFLAG, WORK, WORK[K1M], &WORK[K1], &WORK[K2], &WORK[K3], &WORK[K4],
         &WORK[K5], WORK[K6], WORK[K6 + 1]);
}                            // END


