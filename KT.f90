module eskt

	implicit none
	!!! PARAMETERS !!!
	public :: initial, metropolis, DeltaE, prob, config, accum
	integer, public, parameter :: L=32
	integer, public, parameter :: N_E_val=4000
	integer, public, parameter :: double = selected_real_kind(13)
	
	
	!!! PUBLIC VARIABLES !!!
	real , public :: tmax !max angle for trial
	real, public, dimension(7) :: cum
	integer, public :: nmcs=5000, nequil=1000,accept,cw,ccw
	real, public :: T,E,dE
	real, public, dimension(2) :: M
	real, public, dimension(L,L) :: s
	real, public, dimension(N_E_val) :: w
	
contains
!--FUN1----FUN1----FUN1----FUN1----FUN1----FUN1----FUN1----FUN1--
	function DeltaE(x,y,ang,N) result (DeltaE_result)
		integer, intent (in) :: x,y, N
		real , intent (in) :: ang
		real :: E0, Etr, DeltaE_result
		real :: left, right, up, down
		if (x == 1) then
		   left = s(N,y)
		   right = s(2,y)
		else if (x == N) then
		   left = s(N-1,y)
		   right = s(1,y)
		else
		   left = s(x-1,y)
		   right = s(x+1,y)
		end if
		if (y == 1) then
		   up = s(x,2)
		   down = s(x,N)
		else if (y == N) then
		   up = s(x,1)
		   down = s(x,N-1)
		else
		   up = s(x,y+1)
		   down = s(x,y-1)
		end if
		Etr = -cos(ang-up)-cos(ang-down)-cos(ang-right)-cos(ang-left)
		E0= -cos(s(x,y)-up)-cos(s(x,y)-down)-cos(s(x,y)-right)-cos(s(x,y)-left)
		DeltaE_result= (Etr-E0) 
	end function DeltaE
	

!--SUB1----SUB1----SUB1----SUB1----SUB1----SUB1----SUB1----SUB1--
	subroutine initial(N)
	integer, intent(in) :: N
	character(len=20) :: in_mode
	integer :: x,y,up,right,i
	real :: rnd
		
	!L lattice 
	!T temperature
	!nequil equlibration
	!nmcs mc steps
		
	in_mode='random'
	
	!!! INITIALISATION !!!
		s=0.0
	if(in_mode == 'random') then
		do y=1,N
			do x=1,N
				call random_number(rnd)
				s(x,y)= 2*acos(-1.)*rnd  
			end do
		end do
	end if
	
	if(in_mode == 'zero') then
		s=0.0 
		
	end if

	
	!!! ENERGY COMPUTATION !!!
	E=0.0
	M=0.0
	do y=1,N
		if(y==N) then
			up=1
		else
			up=y+1
		end if
		do x=1,N
			if(x==N) then
				right=1
			else
				right=x+1
			end if
			E=E-cos(s(x,y)-s(x,up))-cos(s(x,y)-s(right,y))
			M(1)=M(1)+cos(s(x,y))
			M(2)=M(2)+sin(s(x,y))
		end do
	end do
	
	accept=0 
	cum=0.0
	
	end subroutine initial
	
!--SUB2----SUB2----SUB2----SUB2----SUB2----SUB2----SUB2----SUB2--
	subroutine prob()
	integer :: j
	
	!!! TRANSITION PROBABILITY !!!
	do j=1,N_E_val
		w(j)=exp(-j*8./(N_E_val*T))
	end do
	
	end subroutine prob
	
!--SUB3----SUB3----SUB3----SUB3----SUB3----SUB3----SUB3----SUB3--
	subroutine metropolis(N)
		integer, intent(in) :: N
		integer :: x,y,i
		real :: rnd, str
		!real (kind=double) :: dE
		
		do i=1,N**2
			call random_number(rnd)
			x = int(N*rnd) + 1
			call random_number(rnd)
			y = int(N*rnd) + 1
			call random_number(rnd)
			str=tmax*(rnd-.5)*2. 
			dE = DeltaE(x,y,s(x,y)+str,N)
			if(dE<=0.0 .and. dE>=-8.5) then
				accept = accept + 1
				E = E + dE 
				M(1)=M(1)+cos(s(x,y)+str)-cos(s(x,y))
				M(2)=M(2)+sin(s(x,y)+str)-sin(s(x,y))
				s(x,y) = s(x,y) + str
				if(s(x,y)>2*acos(-1.)) then
					s(x,y)=s(x,y)-2*acos(-1.)
				else if(s(x,y)<0.0) then
					s(x,y)=s(x,y)+2*acos(-1.)
				end if
			else
				call random_number(rnd)
				if (rnd <= w(int(dE*N_E_val/8.)) .and. dE<=8.5) then
					accept = accept + 1
					E = E + dE
					M(1)=M(1)+cos(s(x,y)+str)-cos(s(x,y))
					M(2)=M(2)+sin(s(x,y)+str)-sin(s(x,y))
					s(x,y) = s(x,y) + str
					if(s(x,y)>2*acos(-1.)) then
						s(x,y)=s(x,y)-2*acos(-1.)
					else if(s(x,y)<0.0) then
						s(x,y)=s(x,y)+2*acos(-1.)
					end if
				end if
			end if
			
		end do
		
	end subroutine metropolis
	
!--SUB4----SUB4----SUB4----SUB4----SUB4----SUB4----SUB4----SUB4--
	subroutine config(file,N)
		integer ,intent(in) :: file, N
		integer :: x,y
		character(len=20) :: name1
		
		write(name1, '(A, I0, A)') "spin",file,".dat"
		
		open(unit=1,file=trim(name1),status='unknown')
		do x=1,N
			do y=1,N
				write(1,*) x,y,.4*cos(s(x,y)),.4*sin(s(x,y)), .4*M(1)/N**2, .4*M(2)/N**2
			end do
		end do
		close(1)
		
	end subroutine config
	
!--SUB5----SUB5----SUB5----SUB5----SUB5----SUB5----SUB5----SUB5--
	
	subroutine vortex(file)
		integer, intent(in) :: file
		integer :: x,y, right,up 
		real :: a,b,c,d,v,a1,b1,c1,d1
		character(len=20) :: name1,name2,name3
		
		write(name1, '(A, I0, A)') "cwv",file,".dat"
		write(name2, '(A, I0, A)') "ccwv",file,".dat"
		write(name3, '(A, I0, A)') "dthetav",file,".dat"
		
		if(file /= 0) then
			open(unit=2,file=trim(name1),status='unknown')
			open(unit=3,file=trim(name2),status='unknown')
		end if
		
		!open(unit=4,file=trim(name3),status='unknown')
		
			cw=0
			ccw=0
			do y=1,L
				if(y==L) then
						up=1
				else
						up=y+1
				end if
				do x=1,L
					if(x==L) then
						right=1
					else
						right=x+1
					end if
					a=(s(x,y)-s(right,y))
					if(a>acos(-1.)) then
						a1=a-2*acos(-1.)
					else if(a<-acos(-1.)) then
						a1=a+2*acos(-1.)
					else 
						a1=a
					end if
					b=(s(right,y)-s(right,up))
					if(b>acos(-1.)) then
						b1=b-2*acos(-1.)
					else if(b<-acos(-1.)) then
						b1=b+2*acos(-1.)
					else
						b1=b
					end if
					c=(s(right,up)-s(x,up))
					if(c>acos(-1.)) then
						c1=c-2*acos(-1.)
					else if(c<-acos(-1.)) then
						c1=c+2*acos(-1.)
					else 
						c1=c
					end if
					d=(s(x,up)-s(x,y))
					if(d>acos(-1.)) then
						d1=d-2*acos(-1.)
					else if(d<-acos(-1.)) then
						d1=d+2*acos(-1.)
					else
						d1=d
					end if
					v=a1+b1+c1+d1
					
						
						if(abs(v-2*acos(-1.))<0.001) then
							if(file /= 0) write(2,*) x+.5, y+.5
							cw=cw+1
						end if
						if(abs(v+2*acos(-1.))<0.001) then
							if(file /= 0) write(3,*) x+.5, y+.5
							ccw=ccw+1
						end if
						!write(4,*) L*x+y,v
						
					
				end do	
			end do
			
			print*, cw, 'clock wise vortices'
			print*, ccw, 'counter clock wise vortices'
			
			close(2)
			close(3)
			!close(4)
	
	end subroutine vortex
	
!--SUB6----SUB6----SUB6----SUB6----SUB6----SUB6----SUB6----SUB6--
		
	subroutine accum(N)
		integer, intent(in) :: N
		integer :: x,y
		real :: dc
		
		cum(1) = cum(1) + M(1)
		cum(2) = cum(2) + M(2)
		cum(4) = cum(4) + (M(1)**2+M(2)**2)
		cum(5) = cum(5) + sqrt(M(1)**2+M(2)**2)
		cum(6) = cum(6) + E**2
		cum(7) = cum(7) + E
		
		do x=1,N 
			do y=1,N
				dc=(sin(s(x,y))*M(2)/N**2+cos(s(x,y))*M(1)/N**2)/(sqrt(M(1)**2+M(2)**2)/N**2)
				if(dc<=1. .and. dc>=-1.) then
					cum(3)=cum(3)+(acos(dc))**2
				!else 
					!print*, 'error'
				end if
			end do
		end do

	end subroutine accum
	
end module eskt


!!!!! PROGRAM !!!!!
		
program esa
use eskt

	integer :: a, b
	integer :: nquench=200
	integer :: nt_q=30, nt_t=10, nt_f=25 !n. of temp steps
	real :: c, chi
	real ::  Tmin_q=.5, dt_q=20., Tmin_t=.5, Tmax_t=1.1, Tmin_f=1., Tmax_f=1.2, scale=.5
	character(len=20) :: mode='quench'
	integer, dimension(8) :: seed
	
	seed=(/1,2,3,4,5,6,7,8/)
    call random_seed(put=seed)

	call initial(L)
	print*, 'energy per site:',E/L**2,', Magnetization per site:', M/L**2
	
	!call config(1,L)
	   
	accept=0
	cum=0.
		
	if(mode == 'quench') then
		
		call config(10,L)
		
		!T=Tmin_q/scale**nt_q
		T=Tmin_q+nt_q*dt_q   
		
		tmax=1.2
		
		do b=1,nt_q+1
			accept=0
			call prob()
			do a=1,nquench
				call metropolis(L)
			end do
			print*, 'acc. prob.:',real(accept)/nquench/L**2,', Temperature:',T,',L:',L
			!T=T*scale
			T=T-dt_q
		end do
		
		!cum=0.0
		!	do a=1,nquench
		!		call metropolis(L)
		!		call accum(L)
		!		!write(98,*) a, E/L**2
		!	end do
		
		call vortex(11)
		call config(11,L)
		
	else if(mode=='temp') then
		
		T=Tmin_t
		dt=(Tmax_t-Tmin_t)/nt_t
		tmax=1.4
		
		
		do b=1,nt_t+1
			accept=0
			call prob()
			do a=1,nequil
				call metropolis(L)
			end do
			cum=0.0
			do a=1,nmcs
				call metropolis(L)
				call accum(L)
			end do
			call vortex(2000+b)
			call config(2000+b,L)
			
			c=(cum(6)/nmcs-(cum(7)/nmcs)**2)/T**2/L**2
			
			write(200,*) T,(cum(4)/nmcs-(cum(5)/nmcs)**2)/T/L**2,c     !,sqrt(cum(1)**2+cum(2)**2)/L**2/nmcs
			write(201,*) T, E/L**2, cum(7)/nmcs/L**2, cw/real(L**2),1./T, real(accept)/(nequil+nmcs)/L**2
			print*, 'acc. prob.:',real(accept)/(nequil+nmcs)/L**2,', Temperature:',T,',L:',L
			T=T+dt
			!tmax=tmax+.05
		end do
		
	else if(mode=='fit') then
		
		T=Tmin_f
		dt=(Tmax_f-Tmin_f)/nt_f
		tmax=1.8
		
		do b=1,nt_f+1
			accept=0
			call prob()
			do a=1,nequil
				call metropolis(L)
			end do
			cum=0.0
			do a=1,nmcs
				call metropolis(L)
				call accum(L)
				if(mod(a,100)==0) write(5000+b,*) a, E/L**2, cum(7)/a/L**2
			end do
			
			c=(cum(6)/nmcs-(cum(7)/nmcs)**2)/T**2/L**2
			chi=(cum(4)/nmcs)/T/L**2
			
			
			write(500,*) T,chi,log(chi), log(chi-(cum(1)**2+cum(2)**2)/nmcs**2/T/L**2)  
			write(501,*) T,c,sqrt(cum(1)**2+cum(2)**2)/L**2/nmcs      ! E/L**2, cum(7)/nmcs/L**2, cw/real(L**2),1./T, real(accept)/(nequil+nmcs)/L**2
			print*, T,real(accept)/(nequil+nmcs)/L**2
			T=T+dt
			!tmax=tmax+.05
		end do
		
		
		
	else if(mode == 'theta') then
		
		T=.1
		call prob()
		
		do i=4,L,4
			tmax=.6
			accept=0
			call initial(i)
		
			do a=1,nequil
				call metropolis(i)
				!write(300+i,*) a, E/real(i**2), M(1)/i**2,M(2)/i**2, dE !!! CHECK EQ
			end do
			cum=0.0
			do a=1,nmcs
				call metropolis(i)
				call accum(i)
			end do
			chi=(cum(4)/nmcs)/T/i**2
			print*, 'acc. prob.:',real(accept)/(nequil+nmcs)/i**2,', N:',i,', tmax:',tmax,',T:',T
			write(30,*) i, cum(3)/nmcs/i**2, chi, real(accept)/(nequil+nmcs)/i**2
			!tmax=tmax*1
		end do
		
	else
		
		T=.5
		call prob()
		tmax=.9
		
		do a=1,nequil
			call metropolis(L)
			write(40,*) a, E/real(L**2), M(1)/L**2,M(2)/L**2, dE
		end do
		print*, 'acc. prob.:',real(accept)/(nequil)/L**2,', Temperature:',T,',L:',L
		
		do a=1,nmcs
			call metropolis(L)
			!call vortex(0)
			call accum(L) 
		end do
	
		
		!print*, cum/L**2/nmcs, sqrt(cum(1)**2+cum(2)**2)/L**2/nmcs
		
	end if
	
	   
	!call config(2,L)
	do a=1,L
		do b=1,L
			if(s(a,b)<0.0 .or. s(a,b)>2*acos(-1.)) print*, a,b,s(a,b)  !check
		end do
	end do
	!call vortex

end program esa