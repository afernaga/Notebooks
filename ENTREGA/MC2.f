C FENOMENS COL·LECTIUS I TRANSICIONS DE FASE
C SEMESTRE DE TARDOR, CURS 2022-2023
C UNIVERSITAT DE BARCELONA

C ADAM FERNANDEZ GARCIA

C SIMULACIO DEL MODEL D'ISING 2-D
C SIMULACIO RESLITZADA MITJANZANT L'ALGORITME DE METROPOLIS


      program MC2
        implicit none
        integer*4 L,NSEED,SEED0,i,j,MCTOT,a,b,IMC,MCINI,MCD,SEED,T
        parameter (L=48)
        integer*4 PBC(0:L+1),N,DE
        integer*2 S(1:L,1:L),suma
        real*8 genrand_real2,energ,TEMP,DELTA,ENERGIA,MAGNE
        REAL*8 W(-8:8)
        CHARACTER*28 NOM
        REAL*8 SUM,SUME,SUME2,SUMM,SUMM2,SUMAM,VARE,VARM,CV,SUS
        REAL*4 TIME1,TIME2
        CHARACTER*30 DATE
        NAMELIST /DADES/ NOM,TEMP,NSEED,SEED0,MCTOT,MCINI,MCD

C DEFINICIO DE N (TAMANY DEL SISTEMA), I INICIALITZACIO DE LA TEMPERATURA I NOM
        N=L*L
        TEMP=1.4d0

        NOM="SM-L-072-TEMP-1500-MCTOT-10K"

C CONTROL DE TEMPS DE CPU INICIAL
        CALL CPU_TIME(TIME1)

C INICIALITZEM VARIABLES PER ALS BUCLES
        NSEED=1000
        SEED0=117654
        MCTOT=10000
        MCINI=1000
        MCD=10

c LLEGIM LES DADES NOM,CTEMP,NSEED,SEED0,MCTOT,MCINI I MCD DEL FITXER DADES
        OPEN(UNIT=10,FILE="MC2.dat")
        READ(10,DADES)
        CLOSE(10)

C OBRIM EL FITXER ON ESCRIUREM ELS RESULTATS
        OPEN(UNIT=13,FILE=NOM//".dat")
        WRITE(13,*) "# L"," TEMP","    SUM","    SUME","   SUME2",
     &    "    VARE","    SUMM", "     SUMAM","     SUMM2","     VARM",
     &    "    CV","    SUS"
C FEM EL LOOP DE TEMPERATURES
        do T=1,50
          TEMP=TEMP+0.02d0
C ESCRIVIM PER PANTALLA LA TEMPERATURA PER VEURE COM VA EVOLUCIONANT
C L'EXECUCIO DEL PROGRAMA I DETECTAR SI ES PARA
          write(*,*) "Temperatura:",TEMP

C INICIALITZEM LES SUMES PER ALS PROMITJOS
C CONTADOR
        SUM=0.d0
C SUMA ENERGIES
        SUME=0.d0
C SUMA ENERGIES^2
        SUME2=0.d0
C SUMA MAGNETITZACIO
        SUMM=0.d0
C SUMA MAGNETITZACIO^2
        SUMM2=0.d0
C SUMA ABS(MAGNETITZACIO)
        SUMAM=0.d0

C TABULEM ELS VALUR DE L'EXPONENCIAL EXP(DE/CTEMP) PER TAL D'OPTIMITZAR
C EL PROGRAMA
          do de=-8,8
            W(de)=dexp(-dble(de)/TEMP)
          end do

C BUCLE PER A LES DIFERENTS SEEDS

          DO SEED=SEED0,SEED0+NSEED-1,1
C ESCRIVIM PER PANTALLA LES SEEDS PER VEURE COM VA EVOLUCIONANT
C L'EXECUCIO DEL PROGRAMA I DETECTAR SI ES PARA

            WRITE(*,*) "LLAVOR=",SEED

c INICIALITZEM LA MATRIU DE SPINS
            call init_genrand(SEED)
            do i=1,L
              do j=1,L
                if (genrand_real2().lt.0.5d0) then
                  S(i,j)=1
                else
                  S(i,j)=-1
                end if
              end do
            end do
C escrivim el vector PBC (CONDICIONS DE CONTORN PERIODIQUES)
            do i=0,L+1
              if (i.eq.0) then
                PBC(i)=L
              else if (i.eq.(L+1)) then
                PBC(i)=1
              else
                PBC(i)=i
              end if
            end do
            ENERGIA=energ(S,L,PBC)


C SIMULACIO DE MONTECARLO PER AL MODEL D'ISING 2D EN UNA XARXA
C QUADRADA, GENERA UNA SEQUENCIA D'ESTATS {S(I,J)} CORRESPONENTS A UNA
C CADENA DE MARKOV.

            do IMC=1,MCTOT
              do j=1,N
C ELEGIM UN SPIN ARBITRARI
                a=int(genrand_real2()*real(L))+1
                b=int(genrand_real2()*real(L))+1
C EVALUEM EL CANVI D'ENERGIA
                suma=S(a,PBC(b+1))+S(a,PBC(b-1))+S(PBC(a+1),b)
     &                   +S(PBC(a-1),b)
C VARIACIO DE L'ENERGIA
                DE=2*S(a,b)*suma
C DETERMINEM SI ACCEPTEM O NO EL CANVI
                if (DE.le.0) then
                  S(a,b)=-S(a,b)
                  ENERGIA=ENERGIA+de
                else
                  delta=genrand_real2()
                  if (delta.lt.W(de)) then
                    S(a,b)=-S(a,b)
                    ENERGIA=ENERGIA+de
                  end if
                end if
              end do

c EVALUEM SI AGAFEM ELS PROMITJOS I ELS CALCULEM
              IF ((IMC.GT.MCINI).AND.(MCD*(IMC/MCD).EQ.IMC)) THEN

                SUM=SUM+1.d0
                SUME=SUME+ENERGIA
                SUME2=SUME2+ENERGIA*ENERGIA

                SUMM=SUMM+MAGNE(S,L)
                SUMAM=SUMAM+dabs(MAGNE(S,L))
                SUMM2=SUMM2+MAGNE(S,L)*MAGNE(S,L)
              END IF
            end do

          END DO

C NORMALITZACIO DE PROMITJOS
          SUME=SUME/SUM
          SUMM=SUMM/SUM
          SUME2=SUME2/SUM
          SUMAM=SUMAM/SUM
          SUMM2=SUMM2/SUM
          VARE=SUME2-SUME*SUME
          VARM=SUMM2-SUMM*SUMM
C CAPACITAT CALORIFICA/(N*KB) I SUSCEPTIBILITAT/(N*KB)
          CV=VARE/(DBLE(N)*TEMP*TEMP)
          SUS=(SUMM2-SUMAM*SUMAM)/(DBLE(N)*TEMP)
C ESCRIVIM ELS RESULTATS EN UN FITXER


          WRITE(13,*) L,TEMP,SUM,SUME,SUME2,VARE,SUMM,SUMAM,SUMM2,VARM,
     &                CV,SUS
C ACABEM EL LOOP DE TEMPERATURES
        end do


C CONTROL DE TEMPS DE CPU FINAL
        CALL CPU_TIME(TIME2)

C DATA I HORA FINAL:
        CALL FDATE(DATE)



        WRITE(*,*) DATE
        WRITE(13,*) "# CPUTIME=", TIME2-TIME1

C TANQUEM EL FITXER ON HEM ESCRIT ELS RESULTATS
        CLOSE(13)

      end program MC2

c-----------------------------------------------------------------------

C FUNCIO PER A CALCULAR L'ENERGIA

      REAL*8 FUNCTION ENERG(S,L,PBC)
      IMPLICIT NONE
        INTEGER*2 S(1:L,1:L)
        INTEGER*4 I,J,L
        INTEGER*4 PBC(0:L+1)
        REAL*8 ENE
        ENE=0.0D0
        DO I =1,L
          DO J=1,L
            ENE=ENE-S(i,j)*S(PBC(I+1),J)-S(i,j)*S(I, PBC(J+1))
          ENDDO
        ENDDO
        ENERG=ENE
        RETURN
      END

c ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C FUNCTION MAGNETITZACIO
c ----------------------------------------------------------------------
c ----------------------------------------------------------------------

      REAL*8 FUNCTION MAGNE(S,L)
        implicit none
        integer*2 S(1:L,1:L)
        integer*4 i,j,L
        real*8 MAG
        MAG=0.d0
        do i=1,L
          do j=1,L
            MAG=MAG+S(i,j)
          end do
        end do
        MAGNE=MAG
        return
      end function MAGNE

c ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C GENERADOR DE NUMEROS ALEATORIS DEL CAMPUS
c ----------------------------------------------------------------------
c ----------------------------------------------------------------------

c
c  A C-program for MT19937, with initialization improved 2002/1/26.
c  Coded by Takuji Nishimura and Makoto Matsumoto.
c
c  Before using, initialize the state by using init_genrand(seed)
c  or init_by_array(init_key, key_length).
c
c  Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
c  All rights reserved.
c  Copyright (C) 2005, Mutsuo Saito,
c  All rights reserved.
c
c  Redistribution and use in source and binary forms, with or without
c  modification, are permitted provided that the following conditions
c  are met:
c
c    1. Redistributions of source code must retain the above copyright
c       notice, this list of conditions and the following disclaimer.
c
c    2. Redistributions in binary form must reproduce the above copyright
c       notice, this list of conditions and the following disclaimer in the
c       documentation and/or other materials provided with the distribution.
c
c    3. The names of its contributors may not be used to endorse or promote
c       products derived from this software without specific prior written
c       permission.
c
c  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
c  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
c  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
c  A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
c  CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
c  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
c  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
c  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
c  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
c  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
c  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
c
c
c  Any feedback is very welcome.
c  http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
c  email: m-mat @ math.sci.hiroshima-u.ac.jp (remove space)
c
c-----------------------------------------------------------------------
c  FORTRAN77 translation by Tsuyoshi TADA. (2005/12/19)
c
c     ---------- initialize routines ----------
c  subroutine init_genrand(seed): initialize with a seed
c  subroutine init_by_array(init_key,key_length): initialize by an array
c
c     ---------- generate functions ----------
c  integer function genrand_int32(): signed 32-bit integer
c  integer function genrand_int31(): unsigned 31-bit integer
c  double precision function genrand_real1(): [0,1] with 32-bit resolution
c  double precision function genrand_real2(): [0,1) with 32-bit resolution
c  double precision function genrand_real3(): (0,1) with 32-bit resolution
c  double precision function genrand_res53(): (0,1) with 53-bit resolution
c
c  This program uses the following non-standard intrinsics.
c    ishft(i,n): If n>0, shifts bits in i by n positions to left.
c                If n<0, shifts bits in i by n positions to right.
c    iand (i,j): Performs logical AND on corresponding bits of i and j.
c    ior  (i,j): Performs inclusive OR on corresponding bits of i and j.
c    ieor (i,j): Performs exclusive OR on corresponding bits of i and j.
c
c-----------------------------------------------------------------------
c     initialize mt(0:N-1) with a seed
c-----------------------------------------------------------------------
      subroutine init_genrand(s)
      integer s
      integer N
      integer DONE
      integer ALLBIT_MASK
      parameter (N=624)
      parameter (DONE=123456789)
      integer mti,initialized
      integer mt(0:N-1)
      common /mt_state1/ mti,initialized
      common /mt_state2/ mt
      common /mt_mask1/ ALLBIT_MASK
c
      call mt_initln
      mt(0)=iand(s,ALLBIT_MASK)
      do 100 mti=1,N-1
        mt(mti)=1812433253*
     &          ieor(mt(mti-1),ishft(mt(mti-1),-30))+mti
        mt(mti)=iand(mt(mti),ALLBIT_MASK)
  100 continue
      initialized=DONE
c
      return
      end
c-----------------------------------------------------------------------
c     initialize by an array with array-length
c     init_key is the array for initializing keys
c     key_length is its length
c-----------------------------------------------------------------------
      subroutine init_by_array(init_key,key_length)
      integer init_key(0:*)
      integer key_length
      integer N
      integer ALLBIT_MASK
      integer TOPBIT_MASK
      parameter (N=624)
      integer i,j,k
      integer mt(0:N-1)
      common /mt_state2/ mt
      common /mt_mask1/ ALLBIT_MASK
      common /mt_mask2/ TOPBIT_MASK
c
      call init_genrand(19650218)
      i=1
      j=0
      do 100 k=max(N,key_length),1,-1
        mt(i)=ieor(mt(i),ieor(mt(i-1),ishft(mt(i-1),-30))*1664525)
     &           +init_key(j)+j
        mt(i)=iand(mt(i),ALLBIT_MASK)
        i=i+1
        j=j+1
        if(i.ge.N)then
          mt(0)=mt(N-1)
          i=1
        endif
        if(j.ge.key_length)then
          j=0
        endif
  100 continue
      do 200 k=N-1,1,-1
        mt(i)=ieor(mt(i),ieor(mt(i-1),ishft(mt(i-1),-30))*1566083941)-i
        mt(i)=iand(mt(i),ALLBIT_MASK)
        i=i+1
        if(i.ge.N)then
          mt(0)=mt(N-1)
          i=1
        endif
  200 continue
      mt(0)=TOPBIT_MASK
c
      return
      end
c-----------------------------------------------------------------------
c     generates a random number on [0,0xffffffff]-interval
c-----------------------------------------------------------------------
      function genrand_int32()
      integer genrand_int32
      integer N,M
      integer DONE
      integer UPPER_MASK,LOWER_MASK,MATRIX_A
      integer T1_MASK,T2_MASK
      parameter (N=624)
      parameter (M=397)
      parameter (DONE=123456789)
      integer mti,initialized
      integer mt(0:N-1)
      integer y,kk
      integer mag01(0:1)
      common /mt_state1/ mti,initialized
      common /mt_state2/ mt
      common /mt_mask3/ UPPER_MASK,LOWER_MASK,MATRIX_A,T1_MASK,T2_MASK
      common /mt_mag01/ mag01
c
      if(initialized.ne.DONE)then
        call init_genrand(21641)
      endif
c
      if(mti.ge.N)then
        do 100 kk=0,N-M-1
          y=ior(iand(mt(kk),UPPER_MASK),iand(mt(kk+1),LOWER_MASK))
          mt(kk)=ieor(ieor(mt(kk+M),ishft(y,-1)),mag01(iand(y,1)))
  100   continue
        do 200 kk=N-M,N-1-1
          y=ior(iand(mt(kk),UPPER_MASK),iand(mt(kk+1),LOWER_MASK))
          mt(kk)=ieor(ieor(mt(kk+(M-N)),ishft(y,-1)),mag01(iand(y,1)))
  200   continue
        y=ior(iand(mt(N-1),UPPER_MASK),iand(mt(0),LOWER_MASK))
        mt(kk)=ieor(ieor(mt(M-1),ishft(y,-1)),mag01(iand(y,1)))
        mti=0
      endif
c
      y=mt(mti)
      mti=mti+1
c
      y=ieor(y,ishft(y,-11))
      y=ieor(y,iand(ishft(y,7),T1_MASK))
      y=ieor(y,iand(ishft(y,15),T2_MASK))
      y=ieor(y,ishft(y,-18))
c
      genrand_int32=y
      return
      end
c-----------------------------------------------------------------------
c     generates a random number on [0,0x7fffffff]-interval
c-----------------------------------------------------------------------
      function genrand_int31()
      integer genrand_int31
      integer genrand_int32
      genrand_int31=int(ishft(genrand_int32(),-1))
      return
      end
c-----------------------------------------------------------------------
c     generates a random number on [0,1]-real-interval
c-----------------------------------------------------------------------
      function genrand_real1()
      double precision genrand_real1,r
      integer genrand_int32
      r=dble(genrand_int32())
      if(r.lt.0.d0)r=r+2.d0**32
      genrand_real1=r/4294967295.d0
      return
      end
c-----------------------------------------------------------------------
c     generates a random number on [0,1)-real-interval
c-----------------------------------------------------------------------
      function genrand_real2()
      double precision genrand_real2,r
      integer genrand_int32
      r=dble(genrand_int32())
      if(r.lt.0.d0)r=r+2.d0**32
      genrand_real2=r/4294967296.d0
      return
      end
c-----------------------------------------------------------------------
c     generates a random number on (0,1)-real-interval
c-----------------------------------------------------------------------
      function genrand_real3()
      double precision genrand_real3,r
      integer genrand_int32
      r=dble(genrand_int32())
      if(r.lt.0.d0)r=r+2.d0**32
      genrand_real3=(r+0.5d0)/4294967296.d0
      return
      end
c-----------------------------------------------------------------------
c     generates a random number on [0,1) with 53-bit resolution
c-----------------------------------------------------------------------
      function genrand_res53()
      double precision genrand_res53
      integer genrand_int32
      double precision a,b
      a=dble(ishft(genrand_int32(),-5))
      b=dble(ishft(genrand_int32(),-6))
      if(a.lt.0.d0)a=a+2.d0**32
      if(b.lt.0.d0)b=b+2.d0**32
      genrand_res53=(a*67108864.d0+b)/9007199254740992.d0
      return
      end
c-----------------------------------------------------------------------
c     initialize large number (over 32-bit constant number)
c-----------------------------------------------------------------------
      subroutine mt_initln
      integer ALLBIT_MASK
      integer TOPBIT_MASK
      integer UPPER_MASK,LOWER_MASK,MATRIX_A,T1_MASK,T2_MASK
      integer mag01(0:1)
      common /mt_mask1/ ALLBIT_MASK
      common /mt_mask2/ TOPBIT_MASK
      common /mt_mask3/ UPPER_MASK,LOWER_MASK,MATRIX_A,T1_MASK,T2_MASK
      common /mt_mag01/ mag01
CC    TOPBIT_MASK = Z'80000000'
CC    ALLBIT_MASK = Z'ffffffff'
CC    UPPER_MASK  = Z'80000000'
CC    LOWER_MASK  = Z'7fffffff'
CC    MATRIX_A    = Z'9908b0df'
CC    T1_MASK     = Z'9d2c5680'
CC    T2_MASK     = Z'efc60000'
      TOPBIT_MASK=1073741824
      TOPBIT_MASK=ishft(TOPBIT_MASK,1)
      ALLBIT_MASK=2147483647
      ALLBIT_MASK=ior(ALLBIT_MASK,TOPBIT_MASK)
      UPPER_MASK=TOPBIT_MASK
      LOWER_MASK=2147483647
      MATRIX_A=419999967
      MATRIX_A=ior(MATRIX_A,TOPBIT_MASK)
      T1_MASK=489444992
      T1_MASK=ior(T1_MASK,TOPBIT_MASK)
      T2_MASK=1875247104
      T2_MASK=ior(T2_MASK,TOPBIT_MASK)
      mag01(0)=0
      mag01(1)=MATRIX_A
      return
      end
