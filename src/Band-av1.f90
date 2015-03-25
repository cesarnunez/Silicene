! Usar ./a.out < in.nfi
! Compilar en 64bits
!  ifort -parallel -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -O3 -ipo -unroll-aggressive
! EN 32
! ifort -L/opt/intel/Compiler/11.1/046/mkl/lib/32 -O3 -ipo -lmkl_intel -lmkl_intel_thread -lmkl_core -liomp5 -lpthread
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!     VERSION CLUSTER LAPACK
!
	  Program bandas_by
	  implicit none

	complex*16 H0_u(900,900),H0_d(900,900)
	complex*16 VLR_up(900,900),VRL_up(900,900)
	complex*16 VLR_down(900,900),VRL_down(900,900)



	complex*16 Hu(900,900),Hd(900,900)
	real*8 Hru(900,900), Hiu(900,900), Hrd(900,900), Hid(900,900) 


	complex*16 t1p, oo, ii, e0a, e0b, tsvp,tsvn,tsvpx,tsvnx,trsvp,trsvn,tpv
	complex*16 t,t1,t2,Mu,Md,Efield,Ms,aux

	real*8 pi,ek,k0,kf,dk,ef,autou(900),autod(900)
	real*8 er(900),ei(900),fv1(900),fv2(900),fm1(900,900)
  

	INTEGER mx,err,i,j,k,np,nt,l,i1,i2,n,datoN,nn, Ndatos
	character :: str1*6,str2*6,str3*2

	real*4 Efex(50), MM(50)
	INTEGER dato(50),ndimeros(50)



	t  = (-1.d0,0.d0)		!(-Re)Hopping
	t1 = (0.d0,-0.00030d0)		!(-Im)Hopping Rashba
	t2 = (0.d0,+0.03d0)		!(+Im)Hopping Spin-Orbit
	pi = datan(1.d0)*4.d0

	Ndatos = 1



	Open(2,FILE='datosB.dat', ACTION='READ')

	do datoN = 1, Ndatos
	read (2,*) dato(datoN),ndimeros(datoN),Efex(datoN),MM(datoN)
	enddo

	Close(2)


	do datoN = 1, Ndatos

	n=Int(ndimeros(datoN)) ; nn=2*n
	Efield=(-1.d0,0.d0)*real( Efex(datoN) )
	Ms=(1.d0,0.d0)*real( MM(datoN) )


	print* , n , efield, ms

  	write (str1, '(f5.2)') Efex(datoN)
  	write (str2, '(f5.2)') MM(datoN)
  	write (str3, '(I2)') NN


!        open(1, FILE='Band-N'//str3//'EE'//str1//'MM'//str2//'.dat')
	open(1, FILE='Band-zz.dat')

       

	oo = (0.d0,0.d0)				!Cero-Cero
 	e0A = (0d0,0.d0)+Efield				!Energía de sitio A
 	e0B = (0d0,0.d0)-Efield				!Energía de sitio B
	ii = dcmplx(dcos(pi/2.d0),dsin(pi/2.d0))	!Numero Complejo "i"



	!Hopping Spin-Orbita UP y DOWN
	tsvp = t2
	tsvn = conjg(tsvp)

!	!Efecto de Proximidad
	Mu = Ms
	Md = -Ms





! Reemplazmos el HC por el llamado a la matriz H0 y las matrices de transferencia H01 por VLR y H10 por VRL
 
! Hamiltoniano0(HX,eX,eY,M,tpv,tsvp,tsvn,n)
! Hopping(HLR,HRL,tpv,tsvp,tsvn,n)

	call Hamiltoniano0(H0_u,e0A,e0B,Mu,t,tsvp,tsvn,n)		!Hamiltoniano Efecto de proximadad UP
	call Hopping(VLR_up,VRL_up,t,tsvp,tsvn,n)			!Hopping Spin-Orbit UP

	call Hamiltoniano0(H0_d,e0A,e0B,Md,t,tsvn,tsvp,n)		!Hamiltoniano Efecto de proximadad DOWN
	call Hopping(VLR_down,VRL_down,t,tsvn,tsvp,n)		!Hopping Spin-Orbit DOWN

 	do i=1,4*n
 	do j=1,4*n
 	aux = VRL_down(i,j)
 	if(aux.eq.tsvn) then
 	write(*,*) i,j,aux
 	endif
 	enddo
 	enddo
 	STOP


! Barrido en k

      k0 = -1d0; kf = 1d0
      nt = 1001
	  dk  = (kf-k0)/(nt-1)

!     Programa principal
      do np = 1,nt
         ek = k0+(np-1)*dk


	do i = 1, 4*n
	  do j = 1, 4*n
		Hu(i,j) = H0_u(i,j) + VLR_up(i,j)*exp(ii*ek*pi) + VRL_up(i,j)*exp(-ii*ek*pi)
		Hru(i,j) = dreal(Hu(i,j))
		Hiu(i,j) = dimag(Hu(i,j))


		Hd(i,j) = H0_d(i,j) + VLR_down(i,j)*exp(ii*ek*pi) + VRL_down(i,j)*exp(-ii*ek*pi)
		Hrd(i,j) = dreal(Hd(i,j))
		Hid(i,j) = dimag(Hd(i,j))
	  enddo
	enddo


         call ch(900,4*n,Hru,Hiu,autou,0,er,ei,fv1,fv2,fm1,err)
         call ch(900,4*n,Hrd,Hid,autod,0,er,ei,fv1,fv2,fm1,err)

         do i = 1,4*n
         write(1,*) ek,autou(i),autod(i)
         enddo

      	enddo !k


      close(1)




	enddo

	  stop
	  End



!=======================================================================================================
!-------------------------------------------------------------------------------------------------------
! SUBROUTINES  
!-------------------------------------------------------------------------------------------------------
!=======================================================================================================



!-------------------------------------------------------------------------------------------------------
! SUBROUTINES SILICENO 
!-------------------------------------------------------------------------------------------------------


	subroutine Identidad(AA,dimension)

	integer i, dimension
	complex*16 AA(900,900)

	AA=0d0
	do i = 1, dimension
	AA(i,i) = (1d0,0d0)
	enddo

	end subroutine Identidad


!-------------------------------------------------------------------------------------------------------
! HAMILTONIANO H0
!-------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------

! Subroutine Hamiltoniano0

	subroutine Hamiltoniano0(HX,eX,eY,M,tpv,tsvp,tsvn,n)

	integer i,j,k,n
	complex*16 :: eX,eY,tpv,tsvp,tsvn,X1,X2,ii,M
	complex*16 :: id(900,900)
	complex*16 :: EE0X(900,900),EE0Y(900,900)

	complex*16 :: Gab(900,900),Gac(900,900)
	complex*16 :: Gbc(900,900),Gbd(900,900)
	complex*16 :: Gcd(900,900)

	complex*16 :: Gba(900,900),Gca(900,900)
	complex*16 :: Gcb(900,900),Gdb(900,900)
	complex*16 :: Gdc(900,900)

	complex*16 :: HX(900,900)
	
	real*8, parameter :: pi = datan(1.d0)*4.d0

	ii = dcmplx(dcos(pi/2.d0),dsin(pi/2.d0))

	

	!Matrix Hamiltoniano capa 

!Terminos Diagonales

	call Identidad(id,n)
	EE0X = (eX+M)*id
	
	do i=1,n-1
	EE0X(i,i+1)=  tsvn 
	EE0X(i+1,i)=  tsvp
	enddo

	EE0Y = (eY+M)*id

	do i=1,n-1
	EE0Y(i,i+1)=  tsvn !ok
	EE0Y(i+1,i)=  tsvp
	enddo

	!Hopping entre capas

!Primeros Vecinos

	!Capa AB - BA

	Gab = id*tpv
	do i = 1, n-1
	Gab(i,i+1) = tpv
	enddo
	Gba = transpose(conjg(Gab)) !ok

	!Capa BC - CB
	
	Gbc = id*tpv
	Gcb = transpose(conjg(Gbc))

	!Capa CD - DC

	Gcd = id*tpv
	do i = 1, n-1
	Gcd(i+1,i) = tpv
	enddo

	Gdc = transpose(conjg(Gcd)) !ok

!Segundos Vecinos

	!Capa AC - CA (spin-orbit)

	Gac = id*tsvn
	Gca = id*tsvp

	do i = 1, n-1
	Gac(i,i+1) = tsvp
	Gca(i+1,i) = tsvn
	enddo           

	!Capa BD-DB (spin-orbit)

	Gbd = id*tsvp
	Gdb = id*tsvn

	do i = 1, n-1
	Gbd(i+1,i) = tsvn
	Gdb(i,i+1) = tsvp
	enddo          

	!Matriz Hamiltoniana Unitaria

!Contiene las cuatro capas para formar la matriz unitara, es decir, una matriz de HX(4n,4n).

!	 	     | H0x  Gab  Gac  0.0  |  
!		     | Gba  H0y  Gbc  Gbd  |
!	 HX(4n,4n) = | Gca  Gcb  H0x  Gcd  |
!	  	     | 0.0  Gdb  Gdc  H0y  |4nx4n

	HX=(0.d0,0.d0)

	! Terminos Diagonales (4 términos)

	do k = 0,3,2
	do i = 1,n
	do j = 1,n
	HX(i+k*n,j+k*n) = EE0X(i,j) !HX(1,1),HX(3,3)
	enddo
	enddo
	enddo

	do k = 1,4,2
	do i = 1,n
	do j = 1,n
	HX(i+k*n,j+k*n) = EE0Y(i,j) ! HX(2,2),HX(4,4)
	enddo
	enddo
	enddo

	!Hopping Primeros Vecinos (6 términos)	

	do i = 1,n
	do j = 1,n
	HX(i,n+j) = Gab(i,j) !HX(1,2)
	enddo
	enddo
  
	do i = 1,n
	do j = 1,n
	HX(n+i,j) = Gba(i,j) !HX(2,1)
	enddo
	enddo


	do i = 1,n
	do j = 1,n
	HX(i+n,j+2*n) = Gbc(i,j) !HX(2,3)
	enddo
	enddo

	do i = 1,n
	do j = 1,n
	HX(i+2*n,j+n) = Gcb(i,j) !HX(3,2)
	enddo
	enddo

	do i = 1,n
	do j = 1,n
	HX(i+2*n,j+3*n) = Gcd(i,j) !HX(3,4)
	enddo
	enddo

	do i = 1,n
	do j = 1,n
	HX(i+3*n,j+2*n) = Gdc(i,j) !HX(4,3)
	enddo
	enddo

	!Hoppings segundos vecinos ( 4 términos Spin-Orbita)


	do i = 1,n
	do j = 1,n
	HX(i,j+2*n) = Gac(i,j) !HX(1,3)
	enddo
	enddo

	do i = 1,n
	do j = 1,n
	HX(i+2*n,j) = Gca(i,j) !HX(3,1)
	enddo
	enddo

	do i = 1,n
	do j = 1,n
	HX(i+n,j+3*n) = Gbd(i,j) !HX(2,4)
	enddo
	enddo

	do i = 1,n
	do j = 1,n
	HX(i+3*n,j+n) = Gdb(i,j) !HX(4,2)
	enddo
	enddo

	end subroutine Hamiltoniano0



!-------------------------------------------------------------------------------------------------------
! HOPPING ENTRE CAPAS UNITARIAS
!-------------------------------------------------------------------------------------------------------
!Subroutine Hopping entre capas unitarias

	subroutine Hopping(HLR,HRL,tpv,tsvp,tsvn,n)

	integer:: i,j,n
	complex*16 :: tpv,tsvp,tsvn,ii,id(900,900)
	complex*16 :: HLR(900,900),HRL(900,900)
	complex*16 :: V_AD(900,900),V_DA(900,900)
	complex*16 :: V_CA(900,900),V_AC(900,900)
	complex*16 :: V_BD(900,900),V_DB(900,900)

	real*8, parameter :: pi = datan(1.d0)*4.d0

	ii = dcmplx(dcos(pi/2.d0),dsin(pi/2.d0))

	call Identidad(id,n)



	HLR = (0.d0,0.d0); HRL = (0.d0,0.d0)

	V_AC=(0.d0,0.d0); V_AD=(0.d0,0.d0) ; V_BD=(0.d0,0.d0)
	V_CA=(0.d0,0.d0); V_DA=(0.d0,0.d0) ; V_DB=(0.d0,0.d0)



	!Matriz Hopping Primeros Vecinos

!Contiene los elementos que unen las Celdas Unitarias.

!	                  | 0.0(AA)      0.0(AB)      0.0/VRL(AC)  0.0/VRL(AD)  |  
!	        	  | 0.0(BA)      0.0(BB)      0.0(BC)      0.0/VRL(BD)  |
!	 VLR/VRL(4n,4n) = | VLR(CA)/0.0  0.0(CB)      0.0(CC)      0.0(CD)      |
!	                  | VLR(DA)/0.0  VLR(DB)/0.0  0.0(CD)      0.0(DD)      |4nx4n

	!RL
	!Union capa AD VRL(1,4)

	V_DA = id*tpv 	
	V_AD = V_DA

	!Segundos Vecinos RL: AC - BD

	!Union capa AC

	V_CA = id*tsvp
	V_AC=  id*tsvn

	do i=1,n-1
	V_CA(i+1,i) = tsvn	! Union capa CA VLR(3,1)
	V_AC(i,i+1) = tsvp	! Union capa AC VRL(1,3)
	enddo

	!Union capa BD

	V_BD = id*tsvp 	
	V_DB = id*tsvn 	

	do i=1,n-1
	V_BD(i+1,i) = tsvn	! Union capa BD VRL(2,4)
	V_DB(i,i+1) = tsvp	! Union capa DB VLR(4,2)
	enddo

	!Llenado Matriz VLR

!	PRIMEROS VECINOS
	do i = 1, n
	do j = 1, n
	HLR(3*n+i,j) = V_DA(i,j)	!VLR(4,1)
	enddo
	enddo

!	SEGUNDOS VECINOS
	do i = 1, n
	do j = 1, n
	HLR(2*n+i,j) = V_CA(i,j)	!VLR(3,1)
	enddo
	enddo

	do i = 1, n
	do j = 1, n
	HLR(3*n+i,n+j) = V_DB(i,j) 	!VLR(4,2)
	enddo
	enddo

	!Llenado Matriz VRL


	HRL = conjg(transpose(HLR))	




	end subroutine Hopping



!-------------------------------------------------------------------------------------------------------
! HAMILTONIANO RASHBA
!-------------------------------------------------------------------------------------------------------
	subroutine Hamiltoniano0R(HX,tsvp,tsvn,n)

	integer i,j,k,n
	complex*16 :: tsvp,tsvn,X1,X2,ii
	complex*16 :: id(900,900)
	complex*16 :: EE0X(900,900),EE0Y(900,900)

	complex*16 :: Gab(900,900),Gac(900,900)
	complex*16 :: Gbc(900,900),Gbd(900,900)
	complex*16 :: Gcd(900,900)

	complex*16 :: Gba(900,900),Gca(900,900)
	complex*16 :: Gcb(900,900),Gdb(900,900)
	complex*16 :: Gdc(900,900)

	complex*16 :: HX(900,900)
	
	real*8, parameter :: pi = datan(1.d0)*4.d0

	ii = dcmplx(dcos(pi/2.d0),dsin(pi/2.d0))

	

	!Matrix Hamiltoniano capa 


	call Identidad(id,n)
	
	do i=1,n-1
	EE0X(i,i+1)=  trsvp !ok
	EE0X(i+1,i)=  trsvp
	enddo


	do i=1,n-1
	EE0Y(i,i+1)=  trsvn !ok
	EE0Y(i+1,i)=  trsvn
	enddo

	!Hopping entre capas
	!Segundos Vecinos

	!Capa AC - CA (spin-orbit Roshba)

	Gac = id*trsvp
	Gca = id*trsvp

	do i = 1, n-1
	Gac(i,i+1) = trsvp
	Gca(i+1,i) = trsvp
	enddo           

	!Capa BD-DB (spin-orbit Roshba)

	Gbd = id*trsvn
	Gdb = id*trsvn

	do i = 1, n-1
	Gbd(i+1,i) = trsvn
	Gdb(i,i+1) = trsvn
	enddo          

	!Matriz Hamiltoniana Unitaria

!Contiene las cuatro capas para formar la matriz unitara, es decir, una matriz de HX(4n,4n).

!	 	     | H0x  Gab  Gac  0.0  |  
!		     | Gba  H0y  Gbc  Gbd  |
!	 HX(4n,4n) = | Gca  Gcb  H0x  Gcd  |
!	  	     | 0.0  Gdb  Gdc  H0y  |4nx4n

	HX=(0.d0,0.d0)


	!Hoppings segundos vecinos ( 4 t√©rminos Spin-Orbita)

	do i = 1,n
	do j = 1,n
	HX(i,j+2*n) = Gac(i,j) !HX(1,3)
	enddo
	enddo

	do i = 1,n
	do j = 1,n
	HX(i+2*n,j) = Gca(i,j) !HX(3,1)
	enddo
	enddo

	do i = 1,n
	do j = 1,n
	HX(i+n,j+3*n) = Gbd(i,j) !HX(2,4)
	enddo
	enddo

	do i = 1,n
	do j = 1,n
	HX(i+3*n,j+n) = Gdb(i,j) !HX(4,2)
	enddo
	enddo

	end subroutine Hamiltoniano0R


!-------------------------------------------------------------------------------------------------------
! HOPPING RASHBA
!-------------------------------------------------------------------------------------------------------


!Subroutine Hopping Rashba entre capas unitarias

	subroutine HoppingR(HrLR,HrRL,trsvp,trsvn,n)

	integer:: i,j,n
	complex*16 :: trsvp,trsvn,ii
		complex*16  id(900,900)
	complex*16 :: HrLR(900,900),HrRL(900,900)
	complex*16 :: VAD(900,900),VDA(900,900)
	complex*16 :: VCA(900,900),VAC(900,900)
	complex*16 :: VBD(900,900),VDB(900,900)

	real*8, parameter :: pi = datan(1.d0)*4.d0

	ii = dcmplx(dcos(pi/2.d0),dsin(pi/2.d0))

	call Identidad(id,n)



	HrLR = (0.d0,0.d0); HrRL = (0.d0,0.d0)

	VAD=(0.d0,0.d0); VDA=(0.d0,0.d0); VCA=(0.d0,0.d0)
	VAC=(0.d0,0.d0); VDB=(0.d0,0.d0); VBD=(0.d0,0.d0)


	!Matriz Hopping Primeros Vecinos

!Contiene los elementos que unen las Celdas Unitarias.

!	                  | 0.0(AA)      0.0(AB)      0.0/VRL(AC)  0.0/VRL(AD)  |  
!	        	  | 0.0(BA)      0.0(BB)      0.0(BC)      0.0/VRL(BD)  |
!	 VLR/VRL(4n,4n) = | VLR(CA)/0.0  0.0(CB)      0.0(CC)      0.0(CD)      |
!	                  | VLR(DA)/0.0  VLR(DB)/0.0  0.0(CD)      0.0(DD)      |4nx4n

	!RL
	!Union capa AD VRL(1,4)

!	VDA = id*tpv 	
!	VAD = VDA

	!Segundos Vecinos RL: AC - BD

	!Union capa AC

	VAC=  id*trsvp

	VCA = id*trsvp

	do i=1,n-1
	VAC(i,i+1) = trsvp	 ! Union capa AC VRL(1,3)
	VCA(i,i+1) = trsvp	  !Union capa CA VLR(3,1)
	enddo

	!Union capa BD

	VBD = id*trsvn 	
	VDB = id*trsvn 	

	do i=1,n-1
	VBD(i+1,i) = trsvn		!Union capa BD VRL(2,4)
	VDB(i+1,i) = trsvn		!Union capa DB VLR(4,2)
	enddo

	!Llenado Matriz VRL

!	do i = 1, n
!	do j = 1, n
!	HrRL(i,3*n+j) = VDA(i,j)	!VRL(1,4)
!	enddo
!	enddo

	do i = 1, n
	do j = 1, n
	HrRL(i,2*n+j) = VAC(i,j)	!VRL(1,3)
	enddo
	enddo

	do i = 1, n
	do j = 1, n
	HrRL(n+i,3*n+j) = VBD(i,j)	!VRL(2,4)
	enddo
	enddo

!Llenado Matriz VLR

!	do i = 1, n
!	do j = 1, n
!	HrLR(3*n+j,i) = VAD(i,j)	!VLR(4,1)
!	enddo
!	enddo

	do i = 1, n
	do j = 1, n
	HrLR(2*n+j,i) = VCA(i,j)	!VLR(3,1)
	enddo
	enddo

	do i = 1, n
	do j = 1, n
	HrLR(3*n+j,n+i) = VDB(i,j) !VLR(4,2)
	enddo
	enddo

	end subroutine HoppingR


!* =====================================================================================
!* NIST Guide to Available Math Software.
!* Fullsource for module CH from package EISPACK.
!* Retrieved from NETLIB on Thu Apr  6 21:41:50 1995.
!* =====================================================================================
      subroutine ch(nm,n,ar,ai,w,matz,zr,zi,fv1,fv2,fm1,ierr)
!c
      integer i,j,n,nm,ierr,matz
      double precision ar(nm,n),ai(nm,n),w(n),zr(nm,n),zi(nm,n),fv1(n),fv2(n),fm1(2,n)
!c
!c     this subroutine calls the recommended sequence of
!c     subroutines from the eigensystem subroutine package (eispack)
!c     to find the eigenvalues and eigenvectors (if desired)
!c     of a complex hermitian matrix.
!c
!c     on input
!c
!c        nm  must be set to the row dimension of the two-dimensional
!c        array parameters as declared in the calling program
!c        dimension statement.
!c
!c        n  is the order of the matrix  a=(ar,ai).
!c
!c        ar  and  ai  contain the real and imaginary parts,
!c        respectively, of the complex hermitian matrix.
!c
!c        matz  is an integer variable set equal to zero if
!c        only eigenvalues are desired.  otherwise it is set to
!c        any non-zero integer for both eigenvalues and eigenvectors.
!c
!c     on output
!c
!c        w  contains the eigenvalues in ascending order.
!c
!c        zr  and  zi  contain the real and imaginary parts,
!c        respectively, of the eigenvectors if matz is not zero.
!c
!c        ierr  is an integer output variable set equal to an error
!c           completion code described in the documentation for tqlrat
!c           and tql2.  the normal completion code is zero.
!c
!c        fv1, fv2, and  fm1  are temporary storage arrays.
!c
!c     questions and comments should be directed to burton s. garbow,
!c     mathematics and computer science div, argonne national laboratory
!c
!c     this version dated august 1983.
!c
!c     ------------------------------------------------------------------
!c
      if (n .le. nm) go to 10
      ierr = 10 * n
      go to 50
!c
   10 call  htridi(nm,n,ar,ai,w,fv1,fv2,fm1)
      if (matz .ne. 0) go to 20
!c     .......... find eigenvalues only ..........
      call  tqlrat(n,w,fv2,ierr)
      go to 50
!c     .......... find both eigenvalues and eigenvectors ..........
   20 do 40 i = 1, n
!c
         do 30 j = 1, n
            zr(j,i) = 0.0d0
   30    continue
!c
         zr(i,i) = 1.0d0
   40 continue
!c
      call  tql2(nm,n,w,fv1,zr,ierr)
      if (ierr .ne. 0) go to 50
      call  htribk(nm,n,ar,ai,fm1,n,zr,zi)
   50 return
      end
      double precision function epslon (x)
      double precision x
!c
!c     estimate unit roundoff in quantities of size x.
!c
      double precision a,b,c,eps
!c
!c     this program should function properly on all systems
!c     satisfying the following two assumptions,
!c        1.  the base used in representing floating point
!c            numbers is not a power of three.
!c        2.  the quantity  a  in statement 10 is represented to
!c            the accuracy used in floating point variables
!c            that are stored in memory.
!c     the statement number 10 and the go to 10 are intended to
!c     force optimizing compilers to generate code satisfying
!c     assumption 2.
!c     under these assumptions, it should be true that,
!c            a  is not exactly equal to four-thirds,
!c            b  has a zero for its last bit or digit,
!c            c  is not exactly equal to one,
!c            eps  measures the separation of 1.0 from
!c                 the next larger floating point number.
!c     the developers of eispack would appreciate being informed
!c     about any systems where these assumptions do not hold.
!c
!c     this version dated 4/6/83.
!c
      a = 4.0d0/3.0d0
   10 b = a - 1.0d0
      c = b + b + b
      eps = dabs(c-1.0d0)
      if (eps .eq. 0.0d0) go to 10
      epslon = eps*dabs(x)
      return
      end
      subroutine htribk(nm,n,ar,ai,tau,m,zr,zi)
!c
      integer i,j,k,l,m,n,nm
      double precision ar(nm,n),ai(nm,n),tau(2,n),zr(nm,m),zi(nm,m)
      double precision h,s,si
!c
!c     this subroutine is a translation of a complex analogue of
!c     the algol procedure trbak1, num. math. 11, 181-195(1968)
!c     by martin, reinsch, and wilkinson.
!c     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
!c
!c     this subroutine forms the eigenvectors of a complex hermitian
!c     matrix by back transforming those of the corresponding
!c     real symmetric tridiagonal matrix determined by  htridi.
!c
!c     on input
!c
!c        nm must be set to the row dimension of two-dimensional
!c          array parameters as declared in the calling program
!c          dimension statement.
!c
!c        n is the order of the matrix.
!c
!c        ar and ai contain information about the unitary trans-
!c          formations used in the reduction by  htridi  in their
!c          full lower triangles except for the diagonal of ar.
!c
!c        tau contains further information about the transformations.
!c
!c        m is the number of eigenvectors to be back transformed.
!c
!c        zr contains the eigenvectors to be back transformed
!c          in its first m columns.
!c
!c     on output
!c
!c        zr and zi contain the real and imaginary parts,
!c          respectively, of the transformed eigenvectors
!c          in their first m columns.
!c
!c     note that the last component of each returned vector
!c     is real and that vector euclidean norms are preserved.
!c
!c     questions and comments should be directed to burton s. garbow,
!c     mathematics and computer science div, argonne national laboratory
!c
!c     this version dated august 1983.
!c
!c     ------------------------------------------------------------------
!c
      if (m .eq. 0) go to 200
!c     .......... transform the eigenvectors of the real symmetric
!c                tridiagonal matrix to those of the hermitian
!c                tridiagonal matrix. ..........
      do 50 k = 1, n
!c
         do 50 j = 1, m
            zi(k,j) = -zr(k,j) * tau(2,k)
            zr(k,j) = zr(k,j) * tau(1,k)
   50 continue
!c
      if (n .eq. 1) go to 200
!c     .......... recover and apply the householder matrices ..........
      do 140 i = 2, n
         l = i - 1
         h = ai(i,i)
         if (h .eq. 0.0d0) go to 140
!c
         do 130 j = 1, m
            s = 0.0d0
            si = 0.0d0
!c
            do 110 k = 1, l
               s = s + ar(i,k) * zr(k,j) - ai(i,k) * zi(k,j)
               si = si + ar(i,k) * zi(k,j) + ai(i,k) * zr(k,j)
  110       continue
!c     .......... double divisions avoid possible underflow ..........
            s = (s / h) / h
            si = (si / h) / h
!c
            do 120 k = 1, l
               zr(k,j) = zr(k,j) - s * ar(i,k) - si * ai(i,k)
               zi(k,j) = zi(k,j) - si * ar(i,k) + s * ai(i,k)
  120       continue
!c
  130    continue
!c
  140 continue
!c
  200 return
      end
      subroutine htridi(nm,n,ar,ai,d,e,e2,tau)
!c
      integer i,j,k,l,n,ii,nm,jp1
      double precision ar(nm,n),ai(nm,n),d(n),e(n),e2(n),tau(2,n)
      double precision f,g,h,fi,gi,hh,si,scale,pythag
!c
!c     this subroutine is a translation of a complex analogue of
!c     the algol procedure tred1, num. math. 11, 181-195(1968)
!c     by martin, reinsch, and wilkinson.
!c     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
!c
!c     this subroutine reduces a complex hermitian matrix
!c     to a real symmetric tridiagonal matrix using
!c     unitary similarity transformations.
!c
!c     on input
!c
!c        nm must be set to the row dimension of two-dimensional
!c          array parameters as declared in the calling program
!c          dimension statement.
!c
!c        n is the order of the matrix.
!c
!c        ar and ai contain the real and imaginary parts,
!c          respectively, of the complex hermitian input matrix.
!c          only the lower triangle of the matrix need be supplied.
!c
!c     on output
!c
!c        ar and ai contain information about the unitary trans-
!c          formations used in the reduction in their full lower
!c          triangles.  their strict upper triangles and the
!c          diagonal of ar are unaltered.
!c
!c        d contains the diagonal elements of the the tridiagonal matrix.
!c
!c        e contains the subdiagonal elements of the tridiagonal
!c          matrix in its last n-1 positions.  e(1) is set to zero.
!c
!c        e2 contains the squares of the corresponding elements of e.
!c          e2 may coincide with e if the squares are not needed.
!c
!c        tau contains further information about the transformations.
!c
!c     calls pythag for  dsqrt(a*a + b*b) .
!c
!c     questions and comments should be directed to burton s. garbow,
!c     mathematics and computer science div, argonne national laboratory
!c
!c     this version dated august 1983.
!c
!c     ------------------------------------------------------------------
!c
      tau(1,n) = 1.0d0
      tau(2,n) = 0.0d0
!c
      do 100 i = 1, n
  100 d(i) = ar(i,i)
!c     .......... for i=n step -1 until 1 do -- ..........
      do 300 ii = 1, n
         i = n + 1 - ii
         l = i - 1
         h = 0.0d0
         scale = 0.0d0
         if (l .lt. 1) go to 130
!c     .......... scale row (algol tol then not needed) ..........
         do 120 k = 1, l
  120    scale = scale + dabs(ar(i,k)) + dabs(ai(i,k))
!c
         if (scale .ne. 0.0d0) go to 140
         tau(1,l) = 1.0d0
         tau(2,l) = 0.0d0
  130    e(i) = 0.0d0
         e2(i) = 0.0d0
         go to 290
!c
  140    do 150 k = 1, l
            ar(i,k) = ar(i,k) / scale
            ai(i,k) = ai(i,k) / scale
            h = h + ar(i,k) * ar(i,k) + ai(i,k) * ai(i,k)
  150    continue
!c
         e2(i) = scale * scale * h
         g = dsqrt(h)
         e(i) = scale * g
         f = pythag(ar(i,l),ai(i,l))
!c     .......... form next diagonal element of matrix t ..........
         if (f .eq. 0.0d0) go to 160
         tau(1,l) = (ai(i,l) * tau(2,i) - ar(i,l) * tau(1,i)) / f
         si = (ar(i,l) * tau(2,i) + ai(i,l) * tau(1,i)) / f
         h = h + f * g
         g = 1.0d0 + g / f
         ar(i,l) = g * ar(i,l)
         ai(i,l) = g * ai(i,l)
         if (l .eq. 1) go to 270
         go to 170
  160    tau(1,l) = -tau(1,i)
         si = tau(2,i)
         ar(i,l) = g
  170    f = 0.0d0
!c
         do 240 j = 1, l
            g = 0.0d0
            gi = 0.0d0
!c     .......... form element of a*u ..........
            do 180 k = 1, j
               g = g + ar(j,k) * ar(i,k) + ai(j,k) * ai(i,k)
               gi = gi - ar(j,k) * ai(i,k) + ai(j,k) * ar(i,k)
  180       continue
!c
            jp1 = j + 1
            if (l .lt. jp1) go to 220
!c
            do 200 k = jp1, l
               g = g + ar(k,j) * ar(i,k) - ai(k,j) * ai(i,k)
               gi = gi - ar(k,j) * ai(i,k) - ai(k,j) * ar(i,k)
  200       continue
!c     .......... form element of p ..........
  220       e(j) = g / h
            tau(2,j) = gi / h
            f = f + e(j) * ar(i,j) - tau(2,j) * ai(i,j)
  240    continue
!c
         hh = f / (h + h)
!c     .......... form reduced a ..........
         do 260 j = 1, l
            f = ar(i,j)
            g = e(j) - hh * f
            e(j) = g
            fi = -ai(i,j)
            gi = tau(2,j) - hh * fi
            tau(2,j) = -gi
!c
            do 260 k = 1, j
               ar(j,k) = ar(j,k) - f * e(k) - g * ar(i,k) + fi * tau(2,k) + gi * ai(i,k)
               ai(j,k) = ai(j,k) - f * tau(2,k) - g * ai(i,k) - fi * e(k) - gi * ar(i,k)
  260    continue
!c
  270    do 280 k = 1, l
            ar(i,k) = scale * ar(i,k)
            ai(i,k) = scale * ai(i,k)
  280    continue
!c
         tau(2,l) = -si
  290    hh = d(i)
         d(i) = ar(i,i)
         ar(i,i) = hh
         ai(i,i) = scale * dsqrt(h)
  300 continue
!c
      return
      end
      subroutine tqlrat(n,d,e2,ierr)
!c
      integer i,j,l,m,n,ii,l1,mml,ierr
      double precision d(n),e2(n)
      double precision b,c,f,g,h,p,r,s,t,epslon,pythag
!c
!c     this subroutine is a translation of the algol procedure tqlrat,
!c     algorithm 464, comm. acm 16, 689(1973) by reinsch.
!c
!c     this subroutine finds the eigenvalues of a symmetric
!c     tridiagonal matrix by the rational ql method.
!c
!c     on input
!c
!c        n is the order of the matrix.
!c
!c        d contains the diagonal elements of the input matrix.
!c
!c        e2 contains the squares of the subdiagonal elements of the
!c          input matrix in its last n-1 positions.  e2(1) is arbitrary.
!c
!c      on output
!c
!c        d contains the eigenvalues in ascending order.  if an
!c          error exit is made, the eigenvalues are correct and
!c          ordered for indices 1,2,...ierr-1, but may not be
!c          the smallest eigenvalues.
!c
!c        e2 has been destroyed.
!c
!c        ierr is set to
!c          zero       for normal return,
!c          j          if the j-th eigenvalue has not been
!c                     determined after 30 iterations.
!c
!c     calls pythag for  dsqrt(a*a + b*b) .
!c
!c     questions and comments should be directed to burton s. garbow,
!c     mathematics and computer science div, argonne national laboratory
!c
!c     this version dated august 1983.
!c
!c     ------------------------------------------------------------------
!c

      ierr = 0
      if (n .eq. 1) go to 1001

      do 100 i = 2, n
  100 e2(i-1) = e2(i)

      f = 0.0d0
      t = 0.0d0
      e2(n) = 0.0d0

      do 290 l = 1, n
         j = 0
         h = dabs(d(l)) + dsqrt(e2(l))
         if (t .gt. h) go to 105
         t = h
         b = epslon(t)
         c = b * b
!c     .......... look for small squared sub-diagonal element ..........
  105    do 110 m = l, n
            if (e2(m) .le. c) go to 120
!c     .......... e2(n) is always zero, so there is no exit
!c                through the bottom of the loop ..........
  110    continue

  120    if (m .eq. l) go to 210
  130    if (j .eq. 30) go to 1000
         j = j + 1
!c     .......... form shift ..........
         l1 = l + 1
         s = dsqrt(e2(l))
         g = d(l)
         p = (d(l1) - g) / (2.0d0 * s)
         r = pythag(p,1.0d0)
         d(l) = s / (p + dsign(r,p))
         h = g - d(l)

         do 140 i = l1, n
  140    d(i) = d(i) - h

         f = f + h
!c     .......... rational ql transformation ..........
         g = d(m)
         if (g .eq. 0.0d0) g = b
         h = g
         s = 0.0d0
         mml = m - l
!c     .......... for i=m-1 step -1 until l do -- ..........
         do 200 ii = 1, mml
            i = m - ii
            p = g * h
            r = p + e2(i)
            e2(i+1) = s * r
            s = e2(i) / r
            d(i+1) = h + s * (h + d(i))
            g = d(i) - e2(i) / g
            if (g .eq. 0.0d0) g = b
            h = g * p / r
  200    continue

         e2(l) = s * g
         d(l) = h
!c     .......... guard against underflow in convergence test ..........
         if (h .eq. 0.0d0) go to 210
         if (dabs(e2(l)) .le. dabs(c/h)) go to 210
         e2(l) = h * e2(l)
         if (e2(l) .ne. 0.0d0) go to 130
  210    p = d(l) + f
!c     .......... order eigenvalues ..........
         if (l .eq. 1) go to 250
!c     .......... for i=l step -1 until 2 do -- ..........
         do 230 ii = 2, l
            i = l + 2 - ii
            if (p .ge. d(i-1)) go to 270
            d(i) = d(i-1)
  230    continue

  250    i = 1
  270    d(i) = p
  290 continue

      go to 1001
!c     .......... set error -- no convergence to an
!c                eigenvalue after 30 iterations ..........
 1000 ierr = l
 1001 return
      end
      subroutine tql2(nm,n,d,e,z,ierr)

      integer i,j,k,l,m,n,ii,l1,l2,nm,mml,ierr
      double precision d(n),e(n),z(nm,n)
      double precision c,c2,c3,dl1,el1,f,g,h,p,r,s,s2,tst1,tst2,pythag
!c
!c     this subroutine is a translation of the algol procedure tql2,
!c     num. math. 11, 293-306(1968) by bowdler, martin, reinsch, and
!c     wilkinson.
!c     handbook for auto. comp., vol.ii-linear algebra, 227-240(1971).
!c
!c     this subroutine finds the eigenvalues and eigenvectors
!c     of a symmetric tridiagonal matrix by the ql method.
!c     the eigenvectors of a full symmetric matrix can also
!c     be found if  tred2  has been used to reduce this
!c     full matrix to tridiagonal form.
!c
!c     on input
!c
!c        nm must be set to the row dimension of two-dimensional
!c          array parameters as declared in the calling program
!c          dimension statement.
!c
!c        n is the order of the matrix.
!c
!c        d contains the diagonal elements of the input matrix.
!c
!c        e contains the subdiagonal elements of the input matrix
!c          in its last n-1 positions.  e(1) is arbitrary.
!c
!c        z contains the transformation matrix produced in the
!c          reduction by  tred2, if performed.  if the eigenvectors
!c          of the tridiagonal matrix are desired, z must contain
!c          the identity matrix.
!c
!c      on output
!c
!c        d contains the eigenvalues in ascending order.  if an
!c          error exit is made, the eigenvalues are correct but
!c          unordered for indices 1,2,...,ierr-1.
!c
!c        e has been destroyed.
!c
!c        z contains orthonormal eigenvectors of the symmetric
!c          tridiagonal (or full) matrix.  if an error exit is made,
!c          z contains the eigenvectors associated with the stored
!c          eigenvalues.
!c
!c        ierr is set to
!c          zero       for normal return,
!c          j          if the j-th eigenvalue has not been
!c                     determined after 30 iterations.
!c
!c     calls pythag for  dsqrt(a*a + b*b) .
!c
!c     questions and comments should be directed to burton s. garbow,
!c     mathematics and computer science div, argonne national laboratory
!c
!c     this version dated august 1983.
!c
!c     ------------------------------------------------------------------
!c
      ierr = 0
      if (n .eq. 1) go to 1001

      do 100 i = 2, n
  100 e(i-1) = e(i)

      f = 0.0d0
      tst1 = 0.0d0
      e(n) = 0.0d0

      do 240 l = 1, n
         j = 0
         h = dabs(d(l)) + dabs(e(l))
         if (tst1 .lt. h) tst1 = h
!c     .......... look for small sub-diagonal element ..........
         do 110 m = l, n
            tst2 = tst1 + dabs(e(m))
            if (tst2 .eq. tst1) go to 120
!c     .......... e(n) is always zero, so there is no exit
!c                through the bottom of the loop ..........
  110    continue

  120    if (m .eq. l) go to 220
  130    if (j .eq. 30) go to 1000
         j = j + 1
!c     .......... form shift ..........
         l1 = l + 1
         l2 = l1 + 1
         g = d(l)
         p = (d(l1) - g) / (2.0d0 * e(l))
         r = pythag(p,1.0d0)
         d(l) = e(l) / (p + dsign(r,p))
         d(l1) = e(l) * (p + dsign(r,p))
         dl1 = d(l1)
         h = g - d(l)
         if (l2 .gt. n) go to 145

         do 140 i = l2, n
  140    d(i) = d(i) - h

  145    f = f + h
!c     .......... ql transformation ..........
         p = d(m)
         c = 1.0d0
         c2 = c
         el1 = e(l1)
         s = 0.0d0
         mml = m - l
!c     .......... for i=m-1 step -1 until l do -- ..........
         do 200 ii = 1, mml
            c3 = c2
            c2 = c
            s2 = s
            i = m - ii
            g = c * e(i)
            h = c * p
            r = pythag(p,e(i))
            e(i+1) = s * r
            s = e(i) / r
            c = p / r
            p = c * d(i) - s * g
            d(i+1) = h + s * (c * g + s * d(i))
!c     .......... form vector ..........
            do 180 k = 1, n
               h = z(k,i+1)
               z(k,i+1) = s * z(k,i) + c * h
               z(k,i) = c * z(k,i) - s * h
  180       continue
!c
  200    continue
!c
         p = -s * s2 * c3 * el1 * e(l) / dl1
         e(l) = s * p
         d(l) = c * p
         tst2 = tst1 + dabs(e(l))
         if (tst2 .gt. tst1) go to 130
  220    d(l) = d(l) + f
  240 continue
!c     .......... order eigenvalues and eigenvectors ..........
      do 300 ii = 2, n
         i = ii - 1
         k = i
         p = d(i)
!c
         do 260 j = ii, n
            if (d(j) .ge. p) go to 260
            k = j
            p = d(j)
  260    continue
!c
         if (k .eq. i) go to 300
         d(k) = d(i)
         d(i) = p
!c
         do 280 j = 1, n
            p = z(j,i)
            z(j,i) = z(j,k)
            z(j,k) = p
  280    continue
!c
  300 continue
!c
      go to 1001
!c     .......... set error -- no convergence to an
!c                eigenvalue after 30 iterations ..........
 1000 ierr = l
 1001 return
      end



!-------------------------------------------------------------------------------------------------------

!!!! esto lo escribi yo
      double precision function pythag(a,b)
       double precision a,b
       pythag = dsqrt(a*a + b*b)
      end function

!	do i=1,4*n
!	do j=1,4*n

!	AUX = HC(i,j)

!	if(aux.ne.0.d0) then
!	write(*,*) i,j,aux
!	endif

!	enddo
!	enddo
!	STOP
