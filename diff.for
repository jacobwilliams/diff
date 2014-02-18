*DECK DIFF
      SUBROUTINE DIFF(IORD,X0,XMIN,XMAX,F,EPS,ACC,DERIV,ERROR,IFAIL)
C
C             NUMERICAL DIFFERENTIATION OF USER DEFINED FUNCTION
C
C                         DAVID KAHANER, NBS (GAITHERSBURG) 
C
C  THE PROCEDURE DIFFERENTIATE CALCULATES THE FIRST, SECOND OR
C   THIRD ORDER DERIVATIVE OF A FUNCTION BY USING NEVILLE'S PROCESS TO
C   EXTRAPOLATE FROM A SEQUENCE OF SIMPLE POLYNOMIAL APPROXIMATIONS BASED ON
C   INTERPOLATING POINTS DISTRIBUTED SYMMETRICALLY ABOUT X0 (OR LYING ONLY ON
C   ONE SIDE OF X0 SHOULD THIS BE NECESSARY).  IF THE SPECIFIED TOLERANCE IS
C   NON-ZERO THEN THE PROCEDURE ATTEMPTS TO SATISFY THIS ABSOLUTE OR RELATIVE
C   ACCURACY REQUIREMENT, WHILE IF IT IS UNSUCCESSFUL OR IF THE TOLERANCE IS
C   SET TO ZERO THEN THE RESULT HAVING THE MINIMUM ACHIEVABLE ESTIMATED ERROR
C   IS RETURNED INSTEAD.
C
C INPUT PARAMETERS:
C IORD = 1, 2 OR 3 SPECIFIES THAT THE FIRST, SECOND OR THIRD ORDER
C   DERIVATIVE,RESPECTIVELY, IS REQUIRED.
C X0 IS THE POINT AT WHICH THE DERIVATIVE OF THE FUNCTION IS TO BE CALCULATED.
C XMIN, XMAX RESTRICT THE INTERPOLATING POINTS TO LIE IN [XMIN, XMAX], WHICH
C   SHOULD BE THE LARGEST INTERVAL INCLUDING X0 IN WHICH THE FUNCTION IS
C   CALCULABLE AND CONTINUOUS.
C F, A REAL PROCEDURE SUPPLIED BY THE USER, MUST YIELD THE VALUE OF THE
C   FUNCTION AT X FOR ANY X IN [XMIN, XMAX] WHEN CALLED BY F(X).
C EPS DENOTES THE TOLERANCE, EITHER ABSOLUTE OR RELATIVE.  EPS=0 SPECIFIES THAT 
C   THE ERROR IS TO BE MINIMISED, WHILE EPS>0 OR EPS<0 SPECIFIES THAT THE
C   ABSOLUTE OR RELATIVE ERROR, RESPECTIVELY, MUST NOT EXCEED ABS(EPS) IF
C   POSSIBLE.  THE ACCURACY REQUIREMENT SHOULD NOT BE MADE STRICTER THAN
C   NECESSARY, SINCE THE AMOUNT OF COMPUTATION TENDS TO INCREASE AS
C   THE MAGNITUDE OF EPS DECREASES, AND IS PARTICULARLY HIGH WHEN EPS=0.
C ACC DENOTES THAT THE ABSOLUTE (ACC>0) OR RELATIVE (ACC<0) ERRORS IN THE
C   COMPUTED VALUES OF THE FUNCTION ARE MOST UNLIKELY TO EXCEED ABS(ACC), WHICH 
C   SHOULD BE AS SMALL AS POSSIBLE.  IF THE USER CANNOT ESTIMATE ACC WITH
C   COMPLETE CONFIDENCE, THEN IT SHOULD BE SET TO ZERO.
C
C OUTPUT PARAMETERS:
C DERIV IS THE CALCULATED VALUE OF THE DERIVATIVE.
C ERROR IS AN ESTIMATED UPPER BOUND ON THE MAGNITUDE OF THE ABSOLUTE ERROR IN
C   THE CALCULATED RESULT.  IT SHOULD ALWAYS BE EXAMINED, SINCE IN EXTREME CASE 
C   MAY INDICATE THAT THERE ARE NO CORRECT SIGNIFICANT DIGITS IN THE VALUE
C   RETURNED FOR DERIVATIVE.
C IFAIL WILL HAVE ONE OF THE FOLLOWING VALUES ON EXIT:
C   0   THE PROCEDURE WAS SUCCESSFUL.
C   1   THE ESTIMATED ERROR IN THE RESULT EXCEEDS THE (NON-ZERO) REQUESTED
C          ERROR, BUT THE MOST ACCURATE RESULT POSSIBLE HAS BEEN RETURNED.
C   2   INPUT DATA INCORRECT (DERIVATIVE AND ERROR WILL BE UNDEFINED).
C   3   THE INTERVAL [XMIN, XMAX] IS TOO SMALL (DERIVATIVE AND ERROR WILL BE
C          UNDEFINED);
C
      EXTERNAL F
      REAL X0,XMIN,XMAX,ACC,DERIV,ERROR,BETA,BETA4,H,H0,H1,H2,
     +NEWH1,NEWH2,HEVAL,HPREV,BASEH,HACC1,HACC2,NHACC1,
     +NHACC2,MINH,MAXH,MAXH1,MAXH2,TDERIV,F0,TWOF0,F1,F2,F3,F4,FMAX,
     +MAXFUN,PMAXF,DF1,DELTAF,PDELTA,Z,ZPOWER,C0F0,C1,C2,C3,DNEW,DPREV,
     +RE,TE,NEWERR,TEMERR,NEWACC,PACC1,PACC2,FACC1,FACC2,ACC0,
     +ACC1,ACC2,RELACC,TWOINF,TWOSUP,S, 
     +D(10),DENOM(10),E(10),MINERR(10),MAXF(0:10),SAVE(0:13),
     +STOREF(-45:45),FACTOR
C
      INTEGER IORD,IFAIL,ETA,INF,SUP,I,J,K,N,NMAX,METHOD,SIGNH,FCOUNT,
     +INIT
      LOGICAL IGNORE(10),CONTIN,SAVED
C
C
C ETA IS THE MINIMUM NUMBER OF SIGNIFICANT BINARY DIGITS (APART FROM THE
C SIGN DIGIT) USED TO REPRESENT THE MANTISSA OF REAL NUMBERS. IT SHOULD
C BE DEVREASED BY ONE IF THE COMPUTER TRUNCATES RATHER THAN ROUNDS.
C INF, SUP ARE THE LARGEST POSSIBLE POSITIVE INTEGERS SUBJECT TO
C 2**(-INF), -2**(-INF), 2**SUP, AND -2**SUP ALL BEING REPRESENTABLE REAL
C NUMBERS.
      ETA=I1MACH(11) - 1
      INF=-I1MACH(12) - 2
      SUP=I1MACH(13)-1
      IF(IORD.LT.1 .OR. IORD.GT.3 .OR. XMAX.LE.XMIN .OR.
     +  X0.GT.XMAX .OR. X0.LT.XMIN) THEN
          IFAIL = 2 
          RETURN
      ENDIF
C
      TWOINF = 2.**(-INF)
      TWOSUP = 2.**SUP
      FACTOR = 2**(FLOAT((INF+SUP))/30.)
      IF(FACTOR.LT.256.)FACTOR=256.
      MAXH1 = XMAX - X0
      SIGNH = 1
      IF(X0-XMIN .LE. MAXH1)THEN
          MAXH2 = X0 - XMIN
      ELSE
          MAXH2 = MAXH1
          MAXH1 = X0 - XMIN
          SIGNH = -1
      ENDIF
      RELACC = 2.**(1-ETA)
      MAXH1 = (1.-RELACC)*MAXH1
      MAXH2 = (1.-RELACC)*MAXH2
      S=128.*TWOINF 
      IF(ABS(X0).GT.128.*TWOINF*2.**ETA) S = ABS(X0)*2.**(-ETA)
      IF(MAXH1.LT.S)THEN
C         INTERVAL TOO SMALL
          IFAIL =3
          RETURN
      ENDIF
      IF(ACC.LT.0.) THEN
          IF(-ACC.GT.RELACC)RELACC = -ACC
          ACC = 0.
      ENDIF
C
C     DETERMINE THE SMALLEST SPACING AT WHICH THE CALCULATED
C     FUNCTION VALUES ARE UNEQUAL NEAR X0.
C
      F0 = F(X0)
      TWOF0 = F0 + F0
      IF(ABS(X0) .GT. TWOINF*2.**ETA) THEN
          H = ABS(X0)*2.**(-ETA)
          Z = 2.
      ELSE
          H = TWOINF
          Z = 64.
      ENDIF
      DF1 = F(X0+SIGNH*H) - F0
   80 IF(DF1 .NE. 0. .OR. Z*H .GT. MAXH1) GOTO 100
      H = Z*H
      DF1 = F(X0+SIGNH*H) - F0
      IF(Z .NE.2.) THEN
          IF(DF1 .NE. 0.) THEN
              H = H/Z
              Z = 2.
              DF1 = 0.
          ELSE
              IF(Z*H .GT. MAXH1) Z = 2. 
          ENDIF
      ENDIF
      GOTO 80
  100 CONTINUE
C
      IF(DF1 .EQ. 0.) THEN
C         CONSTANT FUNCTION
          DERIV = 0.
          ERROR = 0.
          IFAIL = 0 
          RETURN
      ENDIF
      IF(H .GT. MAXH1/128.) THEN
C         MINIMUM H TOO LARGE 
          IFAIL = 3 
          RETURN
      ENDIF
C
      H = 8.*H
      H1 = SIGNH*H
      H0 = H1
      H2 = -H1
      MINH = 2.**(-MIN(INF,SUP)/IORD)
      IF(MINH.LT.H) MINH = H
      IF(IORD.EQ.1) S = 8.
      IF(IORD.EQ.2) S = 9.*SQRT(3.)
      IF(IORD.EQ.3) S = 27.
      IF(MINH.GT.MAXH1/S) THEN
          IFAIL = 3 
          RETURN
      ENDIF
      IF(MINH.GT.MAXH2/S .OR. MAXH2.LT.128.*TWOINF) THEN
          METHOD = 1
      ELSE
          METHOD = 2
      ENDIF
C
C     METHOD 1 USES 1-SIDED FORMULAE, AND METHOD 2 SYMMETRIC.
C         NOW ESTIMATE ACCURACY OF CALCULATED FUNCTION VALUES.
C
      IF(METHOD.NE.2 .OR. IORD.EQ.2) THEN
          IF(X0.NE.0.) THEN
              CALL FACCUR(0.,-H1,ACC0,X0,F,TWOINF,F0,F1)
          ELSE
              ACC0 = 0.
          ENDIF
      ENDIF
C
      IF(ABS(H1) .GT. TWOSUP/128.) THEN 
          HACC1 = TWOSUP
      ELSE
          HACC1 = 128.*H1
      ENDIF
C
      IF(ABS(HACC1)/4. .LT. MINH) THEN
          HACC1 = 4.*SIGNH*MINH
      ELSEIF(ABS(HACC1) .GT. MAXH1) THEN
          HACC1 = SIGNH*MAXH1 
      ENDIF
      F1 = F(X0+HACC1)
      CALL FACCUR(HACC1,H1,ACC1,X0,F,TWOINF,F0,F1)
      IF(METHOD.EQ.2) THEN
          HACC2 = -HACC1
          IF(ABS(HACC2) .GT. MAXH2) HACC2 = -SIGNH * MAXH2
          F1 = F(X0 + HACC2)
          CALL FACCUR(HACC2,H2,ACC2,X0,F,TWOINF,F0,F1)
      ENDIF
      NMAX = 8
      IF(ETA.GT.36) NMAX = 10 
      N = -1
      FCOUNT = 0
      DERIV = 0.
      ERROR = TWOSUP
      INIT = 3
      CONTIN = .TRUE.
C
  130 N = N+1
      IF(.NOT. CONTIN) GOTO 800
C
      IF(INIT.EQ.3) THEN
C         CALCULATE COEFFICIENTS FOR DIFFERENTIATION FORMULAE
C             AND NEVILLE EXTRAPOLATION ALGORITHM 
          IF(IORD.EQ.1) THEN
              BETA=2.
          ELSEIF(METHOD.EQ.2)THEN
              BETA = SQRT(2.) 
          ELSE
              BETA = SQRT(3.) 
          ENDIF
          BETA4 = BETA**4.
          Z = BETA
          IF(METHOD.EQ.2) Z = Z**2
          ZPOWER = 1.
          DO 150 K = 1,NMAX
              ZPOWER = Z*ZPOWER
              DENOM(K) = ZPOWER-1
  150     CONTINUE
          IF(METHOD.EQ.2 .AND. IORD.EQ.1) THEN
              E(1) = 5.
              E(2) = 6.3
              DO 160 I = 3,NMAX
  160             E(I) = 6.81 
        ELSEIF((METHOD.NE.2.AND.IORD.EQ.1) .OR. (METHOD.EQ.2.AND.
     +            IORD.EQ.2)) THEN
              E(1) = 10.
              E(2) = 16.
              E(3) = 20.36
              E(4) = 23.
              E(5) = 24.46
              DO 165 I = 6,NMAX
  165             E(I) = 26.
              IF(METHOD.EQ.2.AND.IORD.EQ.2) THEN
                  DO 170 I = 1,NMAX
  170                  E(I)=2*E(I)
              ENDIF 
          ELSEIF(METHOD.NE.2.AND.IORD.EQ.2) THEN
              E(1) = 17.78
              E(2) = 30.06
              E(3) = 39.66
              E(4) = 46.16
              E(5) = 50.26
              DO 175 I = 6,NMAX
  175             E(I) = 55.
          ELSEIF(METHOD.EQ.2.AND.IORD.EQ.3) THEN
              E(1) = 25.97
              E(2) = 41.22
              E(3) = 50.95
              E(4) = 56.4
              E(5) = 59.3
              DO 180 I = 6,NMAX
  180             E(I) = 62.
          ELSE
              E(1) = 24.5
              E(2) = 40.4
              E(3) = 52.78
              E(4) = 61.2
              E(5) = 66.55
              DO 185 I = 6,NMAX
  185             E(I) = 73.
              C0F0 = -TWOF0/(3.*BETA)
              C1 = 3./(3.*BETA-1.)
              C2 = -1./(3.*(BETA-1.))
              C3 = 1./(3.*BETA*(5.-2.*BETA))
          ENDIF
      ENDIF
C
C
      IF(INIT.GE.2) THEN
C         INITIALIZATION OF STEPLENGTHS, ACCURACY AND OTHER 
C             PARAMETERS
C
          HEVAL = SIGNH*MINH
          H = HEVAL 
          BASEH = HEVAL
          MAXH = MAXH2
          IF(METHOD.EQ.1)MAXH = MAXH1
          DO 300 K = 1,NMAX
              MINERR(K) = TWOSUP
              IGNORE(K) = .FALSE.
  300     CONTINUE
          IF(METHOD.EQ.1) NEWACC = ACC1 
          IF(METHOD.EQ.-1) NEWACC = ACC2
          IF(METHOD.EQ.2) NEWACC = (ACC1+ACC2)/2. 
          IF(NEWACC.LT.ACC) NEWACC = ACC
          IF((METHOD.NE.2 .OR. IORD.EQ.2) .AND. NEWACC.LT.ACC0)
     +            NEWACC = ACC0
          IF(METHOD.NE.-1) THEN
              FACC1 = ACC1
              NHACC1 = HACC1
              NEWH1 = H1
          ENDIF
          IF(METHOD.NE.1) THEN
              FACC2 = ACC2
              NHACC2 = HACC2
              NEWH2 = H2
          ELSE
              FACC2 = 0.
              NHACC2 = 0.
          ENDIF
          INIT = 1
          J = 0
          SAVED = .FALSE.
      ENDIF
C
C     CALCULATE NEW OR INITIAL FUNCTION VALUES
C
      IF(INIT.EQ.1 .AND. (N.EQ.0 .OR. IORD.EQ.1) .AND.
     +        .NOT.(METHOD.EQ.2 .AND. FCOUNT.GE.45)) THEN
          IF(METHOD.EQ.2) THEN
              FCOUNT = FCOUNT + 1
              F1 = F(X0+HEVAL)
              STOREF(FCOUNT) = F1
              F2 = F(X0-HEVAL)
              STOREF(-FCOUNT) = F2
          ELSE
              J = J+1
              IF(J.LE.FCOUNT) THEN
                  F1 = STOREF(J*METHOD) 
              ELSE
                  F1 = F(X0+HEVAL)
              ENDIF 
          ENDIF
      ELSE
          F1 = F(X0+HEVAL)
          IF(METHOD.EQ.2) F2 = F(X0-HEVAL)
      ENDIF
      IF(N.EQ.0) THEN
          IF(METHOD.EQ.2 .AND. IORD.EQ.3) THEN
              PDELTA = F1-F2
              PMAXF = (ABS(F1)+ABS(F2))/2.
              HEVAL = BETA*HEVAL
              F1 = F(X0+HEVAL)
              F2 = F(X0-HEVAL)
              DELTAF = F1-F2
              MAXFUN = (ABS(F1)+ABS(F2))/2.
              HEVAL = BETA*HEVAL
              F1 = F(X0+HEVAL)
              F2 = F(X0-HEVAL)
          ELSEIF(METHOD.NE.2 .AND. IORD.GE.2) THEN
              IF(IORD.EQ.2) THEN
                  F3 = F1
              ELSE
                  F4 = F1
                  HEVAL = BETA*HEVAL
                  F3 = F(X0+HEVAL)
              ENDIF 
              HEVAL = BETA*HEVAL
              F2 = F(X0+HEVAL)
              HEVAL = BETA*HEVAL
              F1 = F(X0+HEVAL)
          ENDIF
      ENDIF
C
C     EVALUATE A NEW APPROXIMATION DNEW TO THE DERIVATIVE
C
      IF(N.GT.NMAX) THEN
          N = NMAX
          DO 400 I = 1,N
  400         MAXF(I-1) = MAXF(I)
      ENDIF
      IF(METHOD.EQ.2) THEN
          MAXF(N) = (ABS(F1)+ABS(F2))/2.
          IF(IORD.EQ.1) THEN
              DNEW = (F1-F2)/2.
          ELSEIF(IORD.EQ.2) THEN
              DNEW = F1+F2-TWOF0
          ELSE
              DNEW = -PDELTA
              PDELTA = DELTAF 
              DELTAF = F1-F2
              DNEW = DNEW + .5*DELTAF
              IF(MAXF(N).LT.PMAXF) MAXF(N) = PMAXF
              PMAXF = MAXFUN
              MAXFUN = (ABS(F1)+ABS(F2))/2.
          ENDIF
      ELSE
          MAXF(N) = ABS(F1)
          IF(IORD.EQ.1) THEN
              DNEW = F1-F0
          ELSEIF(IORD.EQ.2) THEN
              DNEW = (TWOF0-3*F3+F1)/3. 
              IF(MAXF(N).LT.ABS(F3)) MAXF(N) = ABS(F3)
              F3 = F2
              F2 = F1
          ELSE
              DNEW = C3*F1+C2*F2+C1*F4+C0F0
              IF(MAXF(N).LT.ABS(F2)) MAXF(N) = ABS(F2)
              IF(MAXF(N).LT.ABS(F4)) MAXF(N) = ABS(F4)
              F4 = F3
              F3 = F2
              F2 = F1
          ENDIF
      ENDIF
      IF(ABS(H).GT.1) THEN
          DNEW = DNEW/H**IORD 
      ELSE
          IF(128.*ABS(DNEW).GT.TWOSUP*ABS(H)**IORD) THEN
              DNEW = TWOSUP/128.
          ELSE
              DNEW = DNEW/H**IORD
          ENDIF
      ENDIF
C
      IF(INIT.EQ.0) THEN
C         UPDATE ESTIMATED ACCURACY OF FUNCTION VALUES
          NEWACC = ACC
          IF((METHOD.NE.2 .OR. IORD.EQ.2) .AND. NEWACC.LT.ACC0)
     +        NEWACC = ACC0
          IF(METHOD.NE.-1 .AND. ABS(NHACC1).LE.1.125*ABS(HEVAL)/BETA4)
     +               THEN
              NHACC1 = HEVAL
              PACC1 = FACC1
              CALL FACCUR(NHACC1,NEWH1,FACC1,X0,F,TWOINF,F0,F1)
              IF(FACC1.LT.PACC1) FACC1=(3*FACC1+PACC1)/4.
          ENDIF
          IF(METHOD.NE.1 .AND. ABS(NHACC2).LE.1.125*ABS(HEVAL)/BETA4) 
     +            THEN
              IF(METHOD.EQ.2) THEN
                  F1 = F2
                  NHACC2 = -HEVAL
              ELSE
                  NHACC2 = HEVAL
              ENDIF 
              PACC2 = FACC2
              CALL FACCUR(NHACC2,NEWH2,FACC2,X0,F,TWOINF,F0,F1)
              IF(FACC2.LT.PACC2) FACC2 = (3*FACC2+PACC2)/4. 
          ENDIF
          IF(METHOD.EQ.1 .AND. NEWACC.LT.FACC1) NEWACC = FACC1
          IF(METHOD.EQ.-1 .AND. NEWACC.LT.FACC2) NEWACC = FACC2
          IF(METHOD.EQ.2 .AND. NEWACC.LT.(FACC1+FACC2)/2.)
     +            NEWACC = (FACC1+FACC2)/2.
      ENDIF
C
C     EVALUATE SUCCESSIVE ELEMENTS OF THE CURRENT ROW IN THE NEVILLE
C     ARRAY, ESTIMATING AND EXAMINING THE TRUNCATION AND ROUNDING
C     ERRORS IN EACH
C
      CONTIN = N.LT.NMAX
      HPREV = ABS(H)
      FMAX = MAXF(N)
      IF((METHOD.NE.2 .OR. IORD.EQ.2) .AND. FMAX.LT.ABS(F0))
     +        FMAX = ABS(F0)
C
      DO 500 K = 1,N
          DPREV = D(K)
          D(K) = DNEW
          DNEW = DPREV+(DPREV-DNEW)/DENOM(K)
          TE = ABS(DNEW-D(K)) 
          IF(FMAX.LT.MAXF(N-K)) FMAX = MAXF(N-K)
          HPREV = HPREV/BETA
          IF(NEWACC.GE.RELACC*FMAX) THEN
              RE = NEWACC*E(K)
          ELSE
              RE = RELACC*FMAX*E(K)
          ENDIF
          IF(RE.NE.0.) THEN
              IF(HPREV.GT.1) THEN
                  RE = RE/HPREV**IORD
              ELSEIF(2*RE.GT.TWOSUP*HPREV**IORD) THEN
                  RE = TWOSUP/2.
              ELSE
                  RE = RE/HPREV**IORD
              ENDIF 
          ENDIF
          NEWERR = TE+RE
          IF(TE.GT.RE) NEWERR = 1.25*NEWERR
          IF(.NOT. IGNORE(K)) THEN
              IF((INIT.EQ.0 .OR. (K.EQ.2 .AND. .NOT.IGNORE(1)))
     +                .AND. NEWERR.LT.ERROR) THEN 
                  DERIV = D(K)
                  ERROR = NEWERR
              ENDIF 
              IF(INIT.EQ.1 .AND. N.EQ.1) THEN
              TDERIV = D(1)
                  TEMERR = NEWERR
              ENDIF 
              IF(MINERR(K).LT.TWOSUP/4) THEN
                  S = 4*MINERR(K)
              ELSE
                  S = TWOSUP
              ENDIF 
              IF(TE.GT.RE .OR. NEWERR.GT.S) THEN
                  IGNORE(K) = .TRUE.
              ELSE
                  CONTIN = .TRUE.
              ENDIF 
              IF(NEWERR.LT.MINERR(K)) MINERR(K) = NEWERR
              IF(INIT.EQ.1 .AND. N.EQ.2 .AND. K.EQ.1 .AND.
     +                .NOT.IGNORE(1)) THEN
                  IF(NEWERR.LT.TEMERR) THEN
                      TDERIV = D(1)
                      TEMERR = NEWERR
                  ENDIF
                  IF(TEMERR.LT.ERROR) THEN
                      DERIV = TDERIV
                      ERROR = TEMERR
                  ENDIF
              ENDIF 
          ENDIF
  500 CONTINUE
C
      IF(N.LT.NMAX) D(N+1) = DNEW
                 IF(EPS.LT.0.) THEN
          S = ABS(EPS*DERIV)
      ELSE
          S = EPS
      ENDIF
      IF(ERROR.LE.S) THEN
          CONTIN = .FALSE.
      ELSEIF(INIT.EQ.1 .AND. (N.EQ.2 .OR. IGNORE(1))) THEN
          IF((IGNORE(1) .OR. IGNORE(2)) .AND. SAVED) THEN
              SAVED = .FALSE. 
              N = 2 
              H = BETA * SAVE(0)
              HEVAL = BETA*SAVE(1)
              MAXF(0) = SAVE(2)
              MAXF(1) = SAVE(3)
              MAXF(2) = SAVE(4)
              D(1) = SAVE(5)
              D(2) = SAVE(6)
              D(3) = SAVE(7)
              MINERR(1) = SAVE(8)
              MINERR(2) = SAVE(9)
              IF(METHOD.EQ.2 .AND. IORD.EQ.3) THEN
                  PDELTA = SAVE(10)
                  DELTAF = SAVE(11)
                  PMAXF = SAVE(12)
                  MAXFUN = SAVE(13)
              ELSEIF(METHOD.NE.2 .AND. IORD.GE.2) THEN
                  F2 = SAVE(10)
                  F3 = SAVE(11)
                  IF(IORD.EQ.3) F4 = SAVE(12)
              ENDIF 
              INIT = 0
              IGNORE(1) = .FALSE.
              IGNORE(2) = .FALSE.
          ELSEIF(.NOT. (IGNORE(1) .OR. IGNORE(2)) .AND. N.EQ.2
     +            .AND. BETA4*FACTOR*ABS(HEVAL).LE.MAXH) THEN
C             SAVE ALL CURRENT VALUES IN CASE OF RETURN TO
C                 CURRENT POINT
              SAVED = .TRUE.
              SAVE(0) = H
              SAVE(1) = HEVAL 
              SAVE(2) = MAXF(0)
              SAVE(3) = MAXF(1)
              SAVE(4) = MAXF(2)
              SAVE(5) = D(1)
              SAVE(6) = D(2)
              SAVE(7) = D(3)
              SAVE(8) = MINERR(1)
              SAVE(9) = MINERR (2)
              IF(METHOD.EQ.2 .AND. IORD.EQ.3) THEN
                  SAVE(10) = PDELTA
                  SAVE(11) = DELTAF
                  SAVE(12) = PMAXF
                  SAVE(13) = MAXFUN
              ELSEIF(METHOD.NE.2 .AND. IORD.GE.2) THEN
                  SAVE(10) = F2
                  SAVE(11) = F3
                  IF(IORD.EQ.3) SAVE(12) = F4
              ENDIF 
              H = FACTOR*BASEH
              HEVAL = H
              BASEH = H
              N = -1
          ELSE
              INIT = 0
              H = BETA*H
              HEVAL = BETA*HEVAL
          ENDIF
      ELSEIF(CONTIN .AND. BETA*ABS(HEVAL).LE.MAXH) THEN
          H = BETA*H
          HEVAL = BETA*HEVAL
      ELSEIF(METHOD.NE.1) THEN
          CONTIN = .TRUE.
          IF(METHOD.EQ.2) THEN
              INIT = 3
              METHOD = -1
              IF(IORD.NE.2) THEN
                  IF(X0.NE.0.) THEN
                      CALL FACCUR(0.,-H0,ACC0,X0,F,TWOINF,F0,F1)
                  ELSE
                      ACC0 = 0.
                  ENDIF
              ENDIF 
          ELSE
              INIT = 2
              METHOD = 1
          ENDIF
          N = -1
          SIGNH = -SIGNH
      ELSE
          CONTIN = .FALSE.
      ENDIF
      GOTO 130
  800 IF(EPS.LT.0.) THEN
          S = ABS(EPS*DERIV)
      ELSE
          S = EPS
      ENDIF
      IFAIL = 0
      IF(EPS.NE.0. .AND. ERROR.GT.S) IFAIL = 1
      RETURN
      END 
*DECK FACCUR
      SUBROUTINE FACCUR(H0,H1,FACC,X0,F,TWOINF,F0,F1)
      REAL H0,H1,FACC,A0,A1,F00,F2,DELTAF,T0,T1,X0,F,DF(5),F0,F1
     +        ,TWOINF
      INTEGER J
      EXTERNAL F
      T0 = 0.
      T1 = 0.
      IF(H0.NE.0.) THEN
          IF(X0+H0.NE.0.) THEN
              F00 = F1
          ELSE
              H0 = 0.875*H0
              F00 = F(X0+H0)
          ENDIF
          IF(ABS(H1) .GE. 32.*TWOINF) H1 = H1/8.
          IF(16.*ABS(H1) .GT. ABS(H0)) H1 = SIGN(H1,1.)*ABS(H0)/16.
          IF(F(X0+H0-H1) .EQ. F00) THEN 
              IF(256.*ABS(H1) .LE. ABS(H0)) THEN
                  H1 = 2.*H1
   10             IF(F(X0+H0-H1).NE.F00 .OR. 256.*ABS(H1).GT.ABS(H0)) 
     +                    GOTO 20
                  H1 = 2.*H1
                  GOTO 10
   20             H1 = 8.*H1
  
              ELSE
                  H1 = SIGN(H1,1.)*ABS(H0)/16.
              ENDIF 
          ELSE
              IF(256.*TWOINF .LE. ABS(H0)) THEN
   30             IF(F(X0+H0-H1/2.).EQ.F00 .OR. ABS(H1).LT.4.*TWOINF) 
     +                GOTO 40 
                  H1 = H1/2.
                  GOTO 30
   40             CONTINUE
                  H1 = 8.*H1
                  IF(16.*ABS(H1) .GT. ABS(H0)) H1 = SIGN(H1,1.)
     +                     *ABS(H0)/16. 
              ELSE
                  H1 = SIGN(H1,1.)*ABS(H0)/16.
              ENDIF 
          ENDIF
      ELSE
          F00 = F0
      ENDIF
  
      DO 50 J = 1,5 
          F2 = F(X0+H0-FLOAT(2*J-1)*H1) 
          DF(J) = F2 - F00
          T0 = T0+DF(J)
          T1 = T1+FLOAT(2*J-1)*DF(J)
   50 CONTINUE
      A0 = (33.*T0-5.*T1)/73. 
      A1 = (-5.*T0+1.2*T1)/73.
      FACC = ABS(A0)
      DO 70 J = 1,5 
          DELTAF = ABS(DF(J)-(A0+FLOAT(2*J-1)*A1))
          IF(FACC.LT.DELTAF) FACC = DELTAF
   70 CONTINUE
      FACC = 2.*FACC
      RETURN
      END 
