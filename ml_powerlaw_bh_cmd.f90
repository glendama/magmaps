Program Microlens
use mkl_poisson
    implicit double precision (A-H,O-Z)
	parameter(pi = 3.14159265358979324D0)
	parameter(ratio = 50.D0) ! mass_max/mass_min = 50 in a power law; maybe 100?
	!parameter (Npix = 2000)
	!real f(Npix,Npix)
	double precision, allocatable :: x(:), y(:), p(:,:)
	real, allocatable :: xs(:), ys(:), Re(:), xs2(:), ys2(:), Re2(:) 
	real, allocatable :: f(:,:)
	real start, dtime
	character(len = 8) :: arg
!aportion
    double precision :: du(4) = (/-0.5,0.5,0.5,-0.5/)
    double precision :: dv(4) = (/-0.5,-0.5,0.5,0.5/)
    double precision :: up(4), vp(4)
	double precision :: dy_dx(2)
!potential
	integer state
	integer ipar(128)
	double precision, allocatable :: spar(:), bd_ax(:), bd_bx(:), bd_ay(:), bd_by(:)
	type(DFTI_DESCRIPTOR), pointer :: xhandle
	character(4) BCtype

	start = SECNDS(0.0)
	conv = 0.3; shear = 0.3; smooth = 0.0; L0 = 20; N0pix = 2000; bhfrac = 0.0; bhmass = 0.0; alpha = 0.0 !default parameters
	narg = command_argument_count()
	do i = 1, narg
		call get_command_argument(i, arg)
		if ((arg =='-h').or.(arg=='-help')) then
			stop 'Usage: ml_bh_cmd.out [convergence] [shear] [smooth] [size] [Npix] [bhfrac] [log10(bhmass)] [alpha]'
		endif
		if (i == 1) then
			read(arg,*) conv
		endif
		if (i == 2) then
			read(arg,*) shear
		endif
		if (i == 3) then
			read(arg,*) smooth
		endif
		if (i == 4) then
			read(arg,*) L0
		endif
		if (i == 5) then
			read(arg,*) N0pix
		endif
		if (i == 6) then
			read(arg,*) bhfrac
		endif
		if (i == 7) then
			read(arg,*) bhmass
		endif
		if (i == 8) then
			read(arg,*) alpha
		endif
	enddo
	
	L = L0 + 4
	Npix = int(float(L)/L0*N0pix) ! edges effect
	print*, 'N0pix,Npix', N0pix, Npix
	allocate(f(Npix, Npix))
	bhmass = 10.**bhmass
	conv_s = smooth*conv
	conv1 = conv*(1-bhfrac)
	conv2 = conv*bhfrac
	write(*,'(A,3F5.2,I4,I6,1F5.2,1F8.4)') 'conv, shear, smooth, L(Re), Npix', conv, shear, smooth, L0, N0pix, bhfrac, bhmass
	
	dL = dble(L)/Npix
	xL = 1.5*dble(L)/abs(1-conv-shear)
	yL = 1.5*dble(L)/abs(1-conv+shear)
	dh = dL !(0.5,1.0,2.0,4.0)*dL
	du = du*dh/dL
	dv = dv*dh/dL
	Nx = NINT(xL/dh)
	Ny = NINT(yL/dh)
	print*, 'Nx, Ny', Nx, Ny
	allocate(x(Nx+1), y(Ny+1), p(Nx+1, Ny+1)) 
	do i = 1, Nx+1
		x(i) = (-Nx/2+i-1)*dh
	enddo
	do j = 1, Ny+1
		y(j) = (-Ny/2+j-1)*dh
	enddo
! Random stars coordinates
	R = max(xL,yL)
	print*, 'R', R
	Ns = conv1*(1-smooth)/pi*R**2
	allocate(xs(Ns),ys(Ns),Re(Ns))
	call RANDOM_SEED()
	call RANDOM_NUMBER(xs)
	call RANDOM_NUMBER(ys)
	xs = R*(xs-0.5)
	ys = R*(ys-0.5)
	!Re = 1.
	print*, 'Nstar', Ns
! Power law mass distribution
	if (alpha.le.0.01) then ! identical
		Re = 1.
	else ! power-law
		xmean = 1.0 ! mean mass
		xmin = xmean*(2-alpha)/(1-alpha)*(ratio**(1-alpha)-1)/(ratio**(2-alpha)-1)
		xmax = ratio*xmin
		call RANDOM_NUMBER(Re) ! Re ~ mass, really Re2
		Re = ((xmax**(1-alpha)-xmin**(1-alpha)) * Re + xmin**(1-alpha))**(1./(1-alpha))
	endif
	
! Random BH coordinates
	Ns2 = conv2*(1-smooth)/pi*R**2/bhmass
	allocate(xs2(Ns2),ys2(Ns2),Re2(Ns2))
	call RANDOM_NUMBER(xs2)
	call RANDOM_NUMBER(ys2)
	xs2 = R*(xs2-0.5)
	ys2 = R*(ys2-0.5)
	Re2 = bhmass
	print*, 'Nbh', Ns2

! mass distribution
	do k = 1, Ns
		x1 = xs(k)/dh+Nx/2+1
		y1 = ys(k)/dh+Ny/2+1
		i = floor(x1)
		j = floor(y1)
		if (i.ge.1.and.i.le.Nx.and.j.ge.1.and.j.le.Ny) then
			dx = x1 - i
			dy = y1 - j
			p(i,j) = p(i,j) + (1-dx)*(1-dy)
			p(i+1,j) = p(i+1,j) + dx*(1-dy)
			p(i,j+1) = p(i,j+1) + (1-dx)*dy
			p(i+1,j+1) = p(i+1,j+1) + dx*dy
			p(i:i+1,j:j+1) = Re(k)*p(i:i+1,j:j+1)
		endif
	enddo
	do k = 1, Ns2
		x1 = xs2(k)/dh+Nx/2+1
		y1 = ys2(k)/dh+Ny/2+1
		i = floor(x1)
		j = floor(y1)
		if (i.ge.1.and.i.le.Nx.and.j.ge.1.and.j.le.Ny) then
			dx = x1 - i
			dy = y1 - j
			p(i,j) = p(i,j) + (1-dx)*(1-dy)
			p(i+1,j) = p(i+1,j) + dx*(1-dy)
			p(i,j+1) = p(i,j+1) + (1-dx)*dy
			p(i+1,j+1) = p(i+1,j+1) + dx*dy
			p(i:i+1,j:j+1) = Re2(k)*p(i:i+1,j:j+1)
		endif
	enddo
	p = p*(-2.D0*pi/dh/dh)
	allocate(spar(13*Nx/2+7), bd_ax(Ny+1), bd_bx(Ny+1), bd_ay(Nx+1), bd_by(Nx+1))

! boundary conditions
	ax=-Nx*dh/2
	bx= Nx*dh/2
	ay=-Ny*dh/2
	by= Ny*dh/2
	q = 0.0D0

	BCtype = 'DDDD' !Dirichlet conditions
	do k = 1, Ns
		do j = 1,Ny+1
		   r2 = (x(1)-xs(k))**2 + (y(j)-ys(k))**2
		   bd_ax(j) = Re(k)*dlog(r2)/2.D0 + bd_ax(j)
		   r2 = (x(Nx+1)-xs(k))**2 + (y(j)-ys(k))**2
		   bd_bx(j) = Re(k)*dlog(r2)/2.D0 + bd_bx(j)
		enddo
		do i = 1,Nx+1
		   r2 = (x(i)-xs(k))**2 + (y(1)-ys(k))**2
		   bd_ay(i) = Re(k)*dlog(r2)/2.D0 + bd_ay(i)
		   r2 = (x(i)-xs(k))**2 + (y(Ny+1)-ys(k))**2
		   bd_by(i) = Re(k)*dlog(r2)/2.D0 + bd_by(i) 
		enddo
	enddo
	do k = 1, Ns2
		do j = 1,Ny+1
		   r2 = (x(1)-xs2(k))**2 + (y(j)-ys2(k))**2
		   bd_ax(j) = Re2(k)*dlog(r2)/2.D0 + bd_ax(j)
		   r2 = (x(Nx+1)-xs2(k))**2 + (y(j)-ys2(k))**2
		   bd_bx(j) = Re2(k)*dlog(r2)/2.D0 + bd_bx(j)
		enddo
		do i = 1,Nx+1
		   r2 = (x(i)-xs2(k))**2 + (y(1)-ys2(k))**2
		   bd_ay(i) = Re2(k)*dlog(r2)/2.D0 + bd_ay(i)
		   r2 = (x(i)-xs2(k))**2 + (y(Ny+1)-ys2(k))**2
		   bd_by(i) = Re2(k)*dlog(r2)/2.D0 + bd_by(i) 
		enddo
	enddo

	ipar = 0
	call d_init_Helmholtz_2D(ax, bx, ay, by, Nx, Ny, BCtype, q, ipar, spar, state)
	call d_commit_Helmholtz_2D(p, bd_ax, bd_bx, bd_ay, bd_by, xhandle, ipar, spar, state)
	call d_Helmholtz_2D(p, bd_ax, bd_bx, bd_ay, bd_by, xhandle, ipar, spar, state)
	call free_Helmholtz_2D(xhandle, ipar, state)

	dtime = SECNDS(start)
	write(*,*) 'Potential calculation takes', dtime, 'seconds'
	print*, 'Potential: min, max', minval(p), maxval(p)

	f = 0.0
!$OMP PARALLEL DO COLLAPSE(2) PRIVATE(i,j,alphaX,alphaY,p11,p12,p22,a11,a12,a22,detA,uc,vc,up,vp,imin,imax,jmin,jmax,jm,im,dy_dx)
	do i = 2, Nx-1
		do j = 2, Ny-1
			alphaX = (p(i+1,j)-p(i-1,j))/(2.D0*dh)
			alphaY = (p(i,j+1)-p(i,j-1))/(2.D0*dh)
			p11 = (p(i+1,j) - 2.D0*p(i,j) + p(i-1,j))/(dh*dh)
			p22 = (p(i,j+1) - 2.D0*p(i,j) + p(i,j-1))/(dh*dh)
			p12 = (p(i+1,j+1)-p(i+1,j-1)-p(i-1,j+1)+p(i-1,j-1))/(4.D0*dh*dh)
			a11 = 1.D0-conv_s-shear-p11
			a22 = 1.D0-conv_s+shear-p22
			a12 = -p12
			detA = a11*a22 - a12*a12  

			uc = x(i)*(1.D0-conv_s-shear) - alphaX
			vc = y(j)*(1.D0-conv_s+shear) - alphaY
	        uc = uc/dL + Npix/2 + 1
       		vc = vc/dL + Npix/2 + 1
			if (detA >= 0) then
				up = a11*du + a12*dv + uc
           		vp = a12*du + a22*dv + vc 
				dy_dx(1) = a12/a11
				dy_dx(2) = a22/a12
			else
				up = a12*du + a11*dv + uc
 		        vp = a22*du + a12*dv + vc 
				dy_dx(1) = a22/a12
				dy_dx(2) = a12/a11
				detA = -detA 
			endif
!			uwidth = (abs(a11)+abs(a12))/2
!			vwidth = (abs(a22)+abs(a12))/2

!			width_limit = 5.D0
!			uwidth = min(uwidth,width_limit)
!			vwidth = min(vwidth,width_limit)

!			imin = uc - uwidth + 0.5
!			imax = uc + uwidth + 0.5
!			jmin = vc - vwidth + 0.5
!			jmax = vc + vwidth + 0.5

			imin = min(up(1),up(2),up(3),up(4)) + 0.5
			imax = max(up(1),up(2),up(3),up(4)) + 0.5
			jmin = min(vp(1),vp(2),vp(3),vp(4)) + 0.5
			jmax = max(vp(1),vp(2),vp(3),vp(4)) + 0.5

			if ((imin.ge.1).and.(imax.le.Npix).and.(jmin.ge.1).and.(jmax.le.Npix)) then
				do jm = jmin, jmax
					do im = imin, imax
						f(im,jm) = f(im,jm) + area(im, jm, up, vp, dy_dx) / detA
					enddo
				enddo
			endif
		enddo
	enddo
!$OMP END PARALLEL DO
	dtime = SECNDS(start)
	write(*,*) 'Program has used', dtime, 'seconds'

	N1 = (Npix - N0pix)/2
	N2 = N1 + N0pix - 1
	print*, 'N1, N2', N1, N2
	print*, 'f min, max', minval(f(N1:N2, N1:N2)), maxval(f(N1:N2, N1:N2))
!    open(unit = 10, status='replace',file='map1.bin',access='stream') 
    open(unit = 20, status='replace',file='map1.bin', action='write',access='stream')
    do i = N1, N2
    	write(20) f(i,N1:N2)
    enddo
	close(20)
	
	ftheor =  1./((1-conv)**2-shear**2)
	fmean = sum(f)/Npix**2
	print*, 'f theor, mean:', ftheor, fmean
end Program	

function area(im, jm, x, y, dy_dx)
    implicit double precision (A-H,O-Z)
	double precision area
	integer im, jm
	double precision x(4), y(4)
	double precision dy_dx(2)
	double precision   xv(8), yv(8) !vertices
	byte code, code1, code2
	byte codes(5)
	byte, parameter :: INSIDE = 0, LEFT = 2#0001, RIGHT = 2#0010, BOTTOM = 2#0100, TOP = 2#1000
	byte, parameter :: LEFT_BOTTOM = 2#0101, BOTTOM_RIGHT = 2#0110, RIGHT_TOP = 2#1010, TOP_LEFT = 2#1001
	
	xl = im-0.5D0
	xr = im+0.5D0
	yb = jm-0.5D0
	yt = jm+0.5D0	
	do k = 1, 4
		x0 = x(k); y0 = y(k)
		code = 0
		if (x0 < xl) then
			code = LEFT
		elseif (x0 > xr) then
			code = RIGHT
		endif
		if (y0 < yb) then
			code = ior(code, BOTTOM)
		elseif (y0 > yt) then
			code = ior(code, TOP)
		endif
		codes(k) = code
	enddo
	codes(5) = codes(1)

	area = 0
	kv = 0
	do k = 1, 4
		x1 = x(k); y1 = y(k) 
		s = dy_dx(MOD(k-1,2)+1) ! even, odd
		code1 = codes(k); code2 = codes(k+1)
!		print*, k, x1, y1, code1

		select case(code1)
			case(INSIDE)
				kv = kv + 1; xv(kv) = x1; yv(kv) = y1
				select case(code2)
					case(INSIDE)
					case(LEFT)
						yc = y1 + s*(xl-x1)
						kv = kv + 1; xv(kv) = xl; yv(kv) = yc
					case(BOTTOM)
						xc = x1 + (yb-y1)/s
						kv = kv + 1; xv(kv) = xc; yv(kv) = yb
					case(RIGHT)
						yc = y1 + s*(xr-x1)
						kv = kv + 1; xv(kv) = xr; yv(kv) = yc
					case(TOP)
						xc = x1 + (yt-y1)/s
						kv = kv + 1; xv(kv) = xc; yv(kv) = yt
					case(LEFT_BOTTOM)
						yc = y1 + s*(xl-x1)
						if (yc >= yb) then !left cross
							kv = kv + 1; xv(kv) = xl; yv(kv) = yc
						else !bottom cross
							xc = x1 + (yb-y1)/s
							kv = kv + 1; xv(kv) = xc; yv(kv) = yb
						endif
					case(BOTTOM_RIGHT)
						xc = x1 + (yb-y1)/s
						if (xc <= xr) then !bottom cross
							kv = kv + 1; xv(kv) = xc; yv(kv) = yb
						else !right cross
							yc = y1 + s*(xr-x1)
							kv = kv + 1; xv(kv) = xr; yv(kv) = yc
						endif
					case(RIGHT_TOP)
						yc = y1 + s*(xr-x1)
						if (yc <= yt) then !right cross
							kv = kv + 1; xv(kv) = xr; yv(kv) = yc
						else !top cross
							xc = x1 + (yt-y1)/s
							kv = kv + 1; xv(kv) = xc; yv(kv) = yt
						endif
					case(TOP_LEFT)
						xc = x1 + (yt-y1)/s
						if (xc >= xl) then !top cross
							kv = kv + 1; xv(kv) = xc; yv(kv) = yt
						else !left cross
							yc = y1 + s*(xl-x1)
							kv = kv + 1; xv(kv) = xl; yv(kv) = yc
						endif
				end select
			case(LEFT)
				select case(code2)
					case(LEFT, LEFT_BOTTOM, TOP_LEFT)
					case(INSIDE)
						yc = y1 + s*(xl-x1)
						kv = kv + 1; xv(kv) = xl; yv(kv) = yc
					case(BOTTOM)
						yc = y1 + s*(xl-x1)
						if (yc <= yb) then
							kv = kv + 1; xv(kv) = xl; yv(kv) = yb
						else ! add the two vertices
							kv = kv + 1; xv(kv) = xl; yv(kv) = yc
							xc = x1 + (yb-y1)/s
							kv = kv + 1; xv(kv) = xc; yv(kv) = yb
						endif
					case(BOTTOM_RIGHT)
						yc = y1 + s*(xl-x1)
						if (yc <= yb) then
							kv = kv + 1; xv(kv) = xl; yv(kv) = yb
						else !add two vertices
							kv = kv + 1; xv(kv) = xl; yv(kv) = yc
							xc = x1 + (yb-y1)/s
							if (xc <= xr) then
								kv = kv + 1; xv(kv) = xc; yv(kv) = yb
							else
								yc = y1 + s*(xr-x1)
								kv = kv + 1; xv(kv) = xr; yv(kv) = yc
							endif
						endif
					case(RIGHT) ! add two
						yc = y1 + s*(xl-x1)
						kv = kv + 1; xv(kv) = xl; yv(kv) = yc
						yc = y1 + s*(xr-x1)
						kv = kv + 1; xv(kv) = xr; yv(kv) = yc
					case(TOP)
						yc = y1 + s*(xl-x1)
						if (yc >= yt) then
							return
						else ! add the two vertices
							kv = kv + 1; xv(kv) = xl; yv(kv) = yc
							xc = x1 + (yt-y1)/s
							kv = kv + 1; xv(kv) = xc; yv(kv) = yt
						endif
					case(RIGHT_TOP)
						yc = y1 + s*(xl-x1)
						if (yc >= yt) then
							return
						else ! add the two vertices
							kv = kv + 1; xv(kv) = xl; yv(kv) = yc
							xc = x1 + (yt-y1)/s
							if (xc <= xr) then
								kv = kv + 1; xv(kv) = xc; yv(kv) = yt
							else
								yc = y1 + s*(xr-x1)
								kv = kv + 1; xv(kv) = xr; yv(kv) = yc
							endif
						endif
				end select
			case(BOTTOM)
				select case(code2)
					case(BOTTOM, BOTTOM_RIGHT, LEFT_BOTTOM)
					case(INSIDE)
						xc = x1 + (yb-y1)/s
						kv = kv + 1; xv(kv) = xc; yv(kv) = yb
					case(LEFT)
						xc = x1 + (yb-y1)/s
						if (xc <= xl) then
							return
						else ! add the two vertices
							kv = kv + 1; xv(kv) = xc; yv(kv) = yb
							yc = y1 + s*(xl-x1)
							kv = kv + 1; xv(kv) = xl; yv(kv) = yc
						endif
					case(TOP_LEFT)
						xc = x1 + (yb-y1)/s
						if (xc <= xl) then
							return
						else ! add the two vertices
							kv = kv + 1; xv(kv) = xc; yv(kv) = yb
							yc = y1 + s*(xl-x1)
							if (yc <= yt) then
								kv = kv + 1; xv(kv) = xl; yv(kv) = yc
							else
								xc = x1 + (yt-y1)/s
								kv = kv + 1; xv(kv) = xc; yv(kv) = yt
							endif
						endif
					case(RIGHT)
						xc = x1 + (yb-y1)/s
						if (xc >= xr) then
							kv = kv + 1; xv(kv) = xr; yv(kv) = yb
						else ! add the two vertices
							kv = kv + 1; xv(kv) = xc; yv(kv) = yb
							yc = y1 + s*(xr-x1)
							kv = kv + 1; xv(kv) = xr; yv(kv) = yc
						endif
					case(RIGHT_TOP)
						xc = x1 + (yb-y1)/s
						if (xc >= xr) then
							kv = kv + 1; xv(kv) = xr; yv(kv) = yb
						else ! add the two vertices
							kv = kv + 1; xv(kv) = xc; yv(kv) = yb
							yc = y1 + s*(xr-x1)
							if (yc <= yt) then
								kv = kv + 1; xv(kv) = xr; yv(kv) = yc
							else
								xc = x1 + (yt-y1)/s
								kv = kv + 1; xv(kv) = xc; yv(kv) = yt
							endif
						endif
					case(TOP)! add two
						xc = x1 + (yb-y1)/s
						kv = kv + 1; xv(kv) = xc; yv(kv) = yb
						xc = x1 + (yt-y1)/s
						kv = kv + 1; xv(kv) = xc; yv(kv) = yt
				end select
			case(RIGHT)
				select case(code2)
					case(RIGHT, RIGHT_TOP, BOTTOM_RIGHT)
					case(INSIDE)
						yc = y1 + s*(xr-x1)
						kv = kv + 1; xv(kv) = xr; yv(kv) = yc
					case(LEFT) ! add two
						yc = y1 + s*(xr-x1)
						kv = kv + 1; xv(kv) = xr; yv(kv) = yc
						yc = y1 + s*(xl-x1)
						kv = kv + 1; xv(kv) = xl; yv(kv) = yc
					case(BOTTOM)
						yc = y1 + s*(xr-x1)
						if (yc <= yb) then
							return
						else
							kv = kv + 1; xv(kv) = xr; yv(kv) = yc
							xc = x1 + (yb-y1)/s
							kv = kv + 1; xv(kv) = xc; yv(kv) = yb
						endif	
					case(LEFT_BOTTOM)
						yc = y1 + s*(xr-x1)
						if (yc <= yb) then
							return
						else
							kv = kv + 1; xv(kv) = xr; yv(kv) = yc
							xc = x1 + (yb-y1)/s
							if (xc >= xl) then
								kv = kv + 1; xv(kv) = xc; yv(kv) = yb
							else
								yc = y1 + s*(xl-x1)
								kv = kv + 1; xv(kv) = xl; yv(kv) = yc
							endif
						endif
					case(TOP)
						yc = y1 + s*(xr-x1)
						if (yc >= yt) then
							kv = kv + 1; xv(kv) = xr; yv(kv) = yt
						else ! add two
							kv = kv + 1; xv(kv) = xr; yv(kv) = yc
							xc = x1 + (yt-y1)/s
							kv = kv + 1; xv(kv) = xc; yv(kv) = yt
						endif
					case(TOP_LEFT)
						yc = y1 + s*(xr-x1)
						if (yc >= yt) then
							kv = kv + 1; xv(kv) = xr; yv(kv) = yt
						else ! add two
							kv = kv + 1; xv(kv) = xr; yv(kv) = yc
							xc = x1 + (yt-y1)/s
							if (xc >= xl) then
								kv = kv + 1; xv(kv) = xc; yv(kv) = yt
							else
								yc = y1 + s*(xl-x1)
								kv = kv + 1; xv(kv) = xl; yv(kv) = yc
							endif
						endif
				end select
			case(TOP)
				select case(code2)
					case(TOP, TOP_LEFT, RIGHT_TOP)
					case(INSIDE)
						xc = x1 + (yt-y1)/s
						kv = kv + 1; xv(kv) = xc; yv(kv) = yt
					case(LEFT)
						xc = x1 + (yt-y1)/s
						if (xc <= xl) then
							kv = kv + 1; xv(kv) = xl; yv(kv) = yt
						else
							kv = kv + 1; xv(kv) = xc; yv(kv) = yt
							yc = y1 + s*(xl-x1)
							kv = kv + 1; xv(kv) = xl; yv(kv) = yc
						endif  
					case(LEFT_BOTTOM)
						xc = x1 + (yt-y1)/s
						if (xc <= xl) then
							kv = kv + 1; xv(kv) = xl; yv(kv) = yt
						else
							kv = kv + 1; xv(kv) = xc; yv(kv) = yt
							yc = y1 + s*(xl-x1)
							if (yc >= yb) then
								kv = kv + 1; xv(kv) = xl; yv(kv) = yc
							else
								xc = x1 + (yb-y1)/s
								kv = kv + 1; xv(kv) = xc; yv(kv) = yb
							endif
						endif 
					case(BOTTOM) ! add two
						xc = x1 + (yt-y1)/s
						kv = kv + 1; xv(kv) = xc; yv(kv) = yt
						xc = x1 + (yb-y1)/s
						kv = kv + 1; xv(kv) = xc; yv(kv) = yb
					case(RIGHT)
						xc = x1 + (yt-y1)/s
						if (xc >= xr) then
							return
						else
							kv = kv + 1; xv(kv) = xc; yv(kv) = yt
							yc = y1 + s*(xr-x1)
							kv = kv + 1; xv(kv) = xr; yv(kv) = yc
						endif  
					case(BOTTOM_RIGHT)
						xc = x1 + (yt-y1)/s
						if (xc >= xr) then
							return
						else
							kv = kv + 1; xv(kv) = xc; yv(kv) = yt
							yc = y1 + s*(xr-x1)
							if (yc >= yb) then
								kv = kv + 1; xv(kv) = xr; yv(kv) = yc
							else
								xc = x1 + (yb-y1)/s
								kv = kv + 1; xv(kv) = xc; yv(kv) = yb
							endif
						endif  
				end select
			! Two bits corners
			case(LEFT_BOTTOM)
				select case(code2)
					case(LEFT_BOTTOM, LEFT)
					case(INSIDE)
						yc = y1 + s*(xl-x1)
						if (yc >= yb) then
							kv = kv + 1; xv(kv) = xl; yv(kv) = yc
						else
							kv = kv + 1; xv(kv) = xl; yv(kv) = yb
							xc = x1 + (yb-y1)/s
							kv = kv + 1; xv(kv) = xc; yv(kv) = yb
						endif
					case(BOTTOM,BOTTOM_RIGHT)
						kv = kv + 1; xv(kv) = xl; yv(kv) = yb
					case(RIGHT)
						yc = y1 + s*(xl-x1)
						if (yc >= yb) then
							kv = kv + 1; xv(kv) = xl; yv(kv) = yc
							yc = y1 + s*(xr-x1)
							kv = kv + 1; xv(kv) = xr; yv(kv) = yc
						else
							kv = kv + 1; xv(kv) = xl; yv(kv) = yb
							xc = x1 + (yb-y1)/s
							if (xc <= xr) then
								kv = kv + 1; xv(kv) = xc; yv(kv) = yb
								yc = y1 + s*(xr-x1)
								kv = kv + 1; xv(kv) = xr; yv(kv) = yc
							else
								kv = kv + 1; xv(kv) = xr; yv(kv) = yb
							endif
						endif
					case(RIGHT_TOP)
						yc = y1 + s*(xl-x1)
						if (yc >= yt) then
							return
						elseif (yc < yb) then
							kv = kv + 1; xv(kv) = xl; yv(kv) = yb
							xc = x1 + (yb-y1)/s
							if (xc >= xr) then
								kv = kv + 1; xv(kv) = xr; yv(kv) = yb
							else
								kv = kv + 1; xv(kv) = xc; yv(kv) = yb
								yc = y1 + s*(xr-x1)						
								if (yc <= yt) then
									kv = kv + 1; xv(kv) = xr; yv(kv) = yc
								else
									xc = x1 + (yt-y1)/s
									kv = kv + 1; xv(kv) = xc; yv(kv) = yt
								endif
							endif
						else !(yb <= yc < yt)
							kv = kv + 1; xv(kv) = xl; yv(kv) = yc
							xc = x1 + (yt-y1)/s
							if (xc <= xr) then
								kv = kv + 1; xv(kv) = xc; yv(kv) = yt 
							else
								yc = y1 + s*(xr-x1)
								kv = kv + 1; xv(kv) = xr; yv(kv) = yc
							endif
						endif
					case(TOP)
						yc = y1 + s*(xl-x1)
						if (yc >= yt) then
							return
						elseif (yc < yb) then
							kv = kv + 1; xv(kv) = xl; yv(kv) = yb
							xc = x1 + (yb-y1)/s
							kv = kv + 1; xv(kv) = xc; yv(kv) = yb
							xc = x1 + (yt-y1)/s
							kv = kv + 1; xv(kv) = xc; yv(kv) = yt
						else !(yb <= yc < yt)
							kv = kv + 1; xv(kv) = xl; yv(kv) = yc
							xc = x1 + (yt-y1)/s
							kv = kv + 1; xv(kv) = xc; yv(kv) = yt
						endif
					case(TOP_LEFT)
						return
				end select
			case(BOTTOM_RIGHT)
				select case(code2)
					case(BOTTOM_RIGHT, BOTTOM)
					case(INSIDE)
						xc = x1 + (yb-y1)/s
						if (xc <= xr) then
							kv = kv + 1; xv(kv) = xc; yv(kv) = yb
						else
							kv = kv + 1; xv(kv) = xr; yv(kv) = yb
							yc = y1 + s*(xr-x1)
							kv = kv + 1; xv(kv) = xr; yv(kv) = yc
						endif
					case(LEFT)
						xc = x1 + (yb-y1)/s
						if (xc >= xr) then
							kv = kv + 1; xv(kv) = xr; yv(kv) = yb
							yc = y1 + s*(xr-x1)
							kv = kv + 1; xv(kv) = xr; yv(kv) = yc
							yc = y1 + s*(xl-x1)
							kv = kv + 1; xv(kv) = xl; yv(kv) = yc
						elseif (xc <= xl) then
							return
						else
							kv = kv + 1; xv(kv) = xc; yv(kv) = yb
							yc = y1 + s*(xl-x1)
							kv = kv + 1; xv(kv) = xl; yv(kv) = yc
						endif
					case(RIGHT,RIGHT_TOP)
						kv = kv + 1; xv(kv) = xr; yv(kv) = yb
					case(TOP)
						xc = x1 + (yb-y1)/s
						if (xc >= xr) then
							kv = kv + 1; xv(kv) = xr; yv(kv) = yb
							yc = y1 + s*(xr-x1)
							if (yc >= yt) then
								kv = kv + 1; xv(kv) = xr; yv(kv) = yt
							else
								kv = kv + 1; xv(kv) = xr; yv(kv) = yc
								xc = x1 + (yt-y1)/s
								kv = kv + 1; xv(kv) = xc; yv(kv) = yt
							endif
						else
							kv = kv + 1; xv(kv) = xc; yv(kv) = yb
							xc = x1 + (yt-y1)/s
							kv = kv + 1; xv(kv) = xc; yv(kv) = yt
						endif
					case(LEFT_BOTTOM)
						return
					case(TOP_LEFT)
						xc = x1 + (yb-y1)/s
						if (xc <= xl) then
							return
						elseif (xc > xr) then
							kv = kv + 1; xv(kv) = xr; yv(kv) = yb
							yc = y1 + s*(xr-x1)
							if (yc >= yt) then
								kv = kv + 1; xv(kv) = xr; yv(kv) = yt
							else
								kv = kv + 1; xv(kv) = xr; yv(kv) = yc
								xc = x1 + (yt-y1)/s
								if (xc >= xl) then
									kv = kv + 1; xv(kv) = xc; yv(kv) = yt
								else
									yc = y1 + s*(xl-x1)
									kv = kv + 1; xv(kv) = xl; yv(kv) = yc
								endif
							endif
						else !(xl < xc <=xr)
							kv = kv + 1; xv(kv) = xc; yv(kv) = yb
							yc = y1 + s*(xl-x1)
							if (yc <= yt) then
								kv = kv + 1; xv(kv) = xl; yv(kv) = yc
							else
								xc = x1 + (yt-y1)/s
								kv = kv + 1; xv(kv) = xc; yv(kv) = yt
							endif
						endif
 				end select
			case(RIGHT_TOP)
				select case(code2)
					case(RIGHT_TOP, RIGHT)
					case(INSIDE)
						yc = y1 + s*(xr-x1)
						if (yc <= yt) then
							kv = kv + 1; xv(kv) = xr; yv(kv) = yc
						else
							kv = kv + 1; xv(kv) = xr; yv(kv) = yt
							xc = x1 + (yt-y1)/s
							kv = kv + 1; xv(kv) = xc; yv(kv) = yt
						endif
					case(LEFT)
						xc = x1 + (yt-y1)/s
						if (xc <= xl) then
							kv = kv + 1; xv(kv) = xr; yv(kv) = yt
							kv = kv + 1; xv(kv) = xl; yv(kv) = yt
						elseif (xc >= xr) then
							yc = y1 + s*(xr-x1)
							kv = kv + 1; xv(kv) = xr; yv(kv) = yc
							yc = y1 + s*(xl-x1)
							kv = kv + 1; xv(kv) = xl; yv(kv) = yc
						else
							kv = kv + 1; xv(kv) = xr; yv(kv) = yt
							kv = kv + 1; xv(kv) = xc; yv(kv) = yt
							yc = y1 + s*(xl-x1)
							kv = kv + 1; xv(kv) = xl; yv(kv) = yc
						endif
					case(BOTTOM)
						yc = y1 + s*(xr-x1)
						if (yc <= yb) then
							return
						elseif (yc > yt) then
							kv = kv + 1; xv(kv) = xr; yv(kv) = yt
							xc = x1 + (yt-y1)/s
							kv = kv + 1; xv(kv) = xc; yv(kv) = yt
							xc = x1 + (yb-y1)/s
							kv = kv + 1; xv(kv) = xc; yv(kv) = yb
						else !(yb < yc <= yt)
							kv = kv + 1; xv(kv) = xr; yv(kv) = yc
							xc = x1 + (yb-y1)/s
							kv = kv + 1; xv(kv) = xc; yv(kv) = yb
						endif
					case(TOP,TOP_LEFT)
						kv = kv + 1; xv(kv) = xr; yv(kv) = yt
					case(LEFT_BOTTOM)
						yc = y1 + s*(xr-x1)
						if (yc <= yb) then
							return
						elseif (yc > yt) then
							kv = kv + 1; xv(kv) = xr; yv(kv) = yt
							xc = x1 + (yt-y1)/s
							if (xc <= xl) then
								kv = kv + 1; xv(kv) = xl; yv(kv) = yt
							else
								kv = kv + 1; xv(kv) = xc; yv(kv) = yt
								yc = y1 + s*(xl-x1)
								if (yc >= yb) then
									kv = kv + 1; xv(kv) = xl; yv(kv) = yc
								else
									xc = x1 + (yb-y1)/s
									kv = kv + 1; xv(kv) = xc; yv(kv) = yb
								endif
							endif
						else !(yb < yc < yt)
							kv = kv + 1; xv(kv) = xr; yv(kv) = yc
							xc = x1 + (yb-y1)/s
							if (xc >= xl) then
								kv = kv + 1; xv(kv) = xc; yv(kv) = yb
							else
								yc = y1 + s*(xl-x1)
								kv = kv + 1; xv(kv) = xl; yv(kv) = yc
							endif
						endif
					case(BOTTOM_RIGHT)
						return
				end select
			case(TOP_LEFT)
				select case(code2)
					case(TOP_LEFT, TOP)
					case(INSIDE)
						xc = x1 + (yt-y1)/s
						if (xc >= xl) then
							kv = kv + 1; xv(kv) = xc; yv(kv) = yt
						else
							kv = kv + 1; xv(kv) = xl; yv(kv) = yt
							yc = y1 + s*(xl-x1)
							kv = kv + 1; xv(kv) = xl; yv(kv) = yc
						endif
					case(LEFT,LEFT_BOTTOM)
						kv = kv + 1; xv(kv) = xl; yv(kv) = yt
					case(BOTTOM)
						xc = x1 + (yt-y1)/s
						if (xc >= xl) then
							kv = kv + 1; xv(kv) = xc; yv(kv) = yt
							xc = x1 + (yb-y1)/s
							kv = kv + 1; xv(kv) = xc; yv(kv) = yb
						else
							kv = kv + 1; xv(kv) = xl; yv(kv) = yt
							yc = y1 + s*(xl-x1)
							if (yc <= yb) then
								kv = kv + 1; xv(kv) = xl; yv(kv) = yb
							else
								kv = kv + 1; xv(kv) = xl; yv(kv) = yc
								xc = x1 + (yb-y1)/s
								kv = kv + 1; xv(kv) = xc; yv(kv) = yb
							endif
						endif
					case(RIGHT)
						xc = x1 + (yt-y1)/s
						if (xc >= xr) then
							return
						elseif (xc < xl) then
							kv = kv + 1; xv(kv) = xl; yv(kv) = yt
							yc = y1 + s*(xl-x1)
							kv = kv + 1; xv(kv) = xl; yv(kv) = yc
							yc = y1 + s*(xr-x1)
							kv = kv + 1; xv(kv) = xr; yv(kv) = yc
						else
							kv = kv + 1; xv(kv) = xc; yv(kv) = yt
							yc = y1 + s*(xr-x1)
							kv = kv + 1; xv(kv) = xr; yv(kv) = yc
						endif
					case(BOTTOM_RIGHT)
						xc = x1 + (yt-y1)/s
						if (xc >= xr) then
							return
						elseif (xc < xl) then
							kv = kv + 1; xv(kv) = xl; yv(kv) = yt
							yc = y1 + s*(xl-x1)
							if (yc <= yb) then
								kv = kv + 1; xv(kv) = xl; yv(kv) = yb
							else
								kv = kv + 1; xv(kv) = xl; yv(kv) = yc
								xc = x1 + (yb-y1)/s
								if (xc <= xr) then
									kv = kv + 1; xv(kv) = xc; yv(kv) = yb
								else
									yc = y1 + s*(xr-x1)
									kv = kv + 1; xv(kv) = xr; yv(kv) = yc
								endif
							endif
						else !(xl <= xc < xr)
							kv = kv + 1; xv(kv) = xc; yv(kv) = yt
							yc = y1 + s*(xr-x1)
							if (yc >= yb) then
								kv = kv + 1; xv(kv) = xr; yv(kv) = yc
							else
								xc = x1 + (yb-y1)/s
								kv = kv + 1; xv(kv) = xc; yv(kv) = yb
							endif
						endif
					case(RIGHT_TOP)
						return
				end select
		end select
	enddo

!	print*, 'kv', kv
!	print*, 'xv', xv(1:kv)
!	print*, 'yv', yv(1:kv)

	do k = 1, kv-1	
		area = area + xv(k)*yv(k+1)-xv(k+1)*yv(k)
	enddo
	area = area + xv(kv)*yv(1) - xv(1)*yv(kv)
	area = 0.5D0*area
end function area
