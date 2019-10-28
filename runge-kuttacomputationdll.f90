!
!  Dll.f90
!
!  Fortran DLL-Exported Subroutine
!  Generated by PGI Visual Fortran(R)
!  4/29/2018 1:24:17 PM
!
  subroutine integrate(cp,x,y,z,vx,vy,vz,m,md,iter,datastorage,time,startenergy,merge,mergeenergy)
	  !DEC$ ATTRIBUTES DLLEXPORT:: INTEGRATE
        integer i,j,l
	    integer cp
		integer iter
		integer indexi
		integer indexj
		integer r 
		real*8 dvdt(cp,15)
		real*8 xyzstart(cp,3)
		real*8 dxyzdt(cp,9)
		real*8 dtmiddle
		real*8 mergeenergy
		real*8 deltaenergy 
		real*8 datastorage(3*cp*(iter+1))
		real*8 g
		real*8 startenergy
		real*8 kineticenergy
		real*8 potentialenergy
		real*8 energycoeff
		real*8 coeff
		real*8 dt
		real*8 sumtime
		real*8 sigma
		real*8 maxv
		real*8 vnew
		real*8 mindist
		real*8 m(cp)
		real*8 time(iter+1)
		real*8 temp(cp)
		real*8 tempx(cp)
		real*8 tempy(cp)
		real*8 tempz(cp)
		real*8 xacc(cp)
		real*8 yacc(cp)
		real*8 zacc(cp)
		real*8 md (cp,cp)
		real*8 x(cp)
		real*8 y(cp)
		real*8 z(cp)
		real*8 vx(cp)
		real*8 vy(cp)
		real*8 vz(cp)
		real*8 velocitycentermass(3)
		logical merge

		g=6.67408d-11	
		time(1)=0d0
		sumtime=0d0
		indexi=0
		indexj=0	
		correction=0d0

		do i=1, cp
			datastorage(3*i-2) = x(i)
			datastorage(3*i-1) = y(i)
			datastorage(3*i) = z(i)
		end do

		grandloop:do l = 1, iter
				   
		   mindist=1d295
		   do i=1,cp
				do j=i+1,cp
					if(mindist.GT.md(i,j)) then
						mindist=md(i,j)
						indexi=i
						indexj=j
					end if
				end do
		   end do

		   maxv=0
		   do i=1,cp
				
				if(maxv<dsqrt(vx(i)*vx(i)+vy(i)*vy(i)+vz(i)*vz(i))) then
						maxv=dsqrt(vx(i)*vx(i)+vy(i)*vy(i)+vz(i)*vz(i))
				end if

		   end do
	   		   
		   if(mindist<1d9) then
				
				merge = .TRUE.

				mergeenergy = mergeenergy -g*m(indexi)*m(indexj)/mindist

				vx(indexi)=(m(indexi) * vx(indexi) + m(indexj) * vx(indexj))/(m(indexi) + m(indexj))
			    vy(indexi)=(m(indexi) * vy(indexi) + m(indexj) * vy(indexj))/(m(indexi) + m(indexj))
				vz(indexi)=(m(indexi) * vz(indexi) + m(indexj) * vz(indexj))/(m(indexi) + m(indexj))
				
				
				x(indexi) = (m(indexi)*x(indexi)+m(indexj)*x(indexj))/(m(indexi)+m(indexj))
				y(indexi) = (m(indexi)*y(indexi)+m(indexj)*y(indexj))/(m(indexi)+m(indexj))
				z(indexi) = (m(indexi)*z(indexi)+z(indexj)*z(indexj))/(m(indexi)+m(indexj))

				m(indexi)=m(indexi) + m(indexj)

			    do i=indexj,cp-1

				x(i)=x(i+1)
				y(i)=y(i+1)
				z(i)=z(i+1)
				vx(i)=vx(i+1)
				vy(i)=vy(i+1)
				vz(i)=vz(i+1)

				end do

				cp=cp-1 	

				do i=l, iter
					
					time(i+1)=time(l)

				end do

				do i = l, iter
					
					do j = 1, cp+1

						datastorage(1+3*(i*cp+j-1))=x(j)
						datastorage(2+3*(i*cp+j-1))=y(j)
						datastorage(3+3*(i*cp+j-1))=z(j)

					end do

				end do	

				return 

		   end if

		  dt = mindist/(maxv*1d2)
		  sumtime=dt+sumtime
		  time(l+1)=sumtime
		  dtmiddle = dt/4d0

		  do i = 1, cp

			xyzstart(i,1) = x(i)
			xyzstart(i,2) = y(i)
			xyzstart(i,3) = z(i)
			dxyzdt(i,1) = vx(i)
			dxyzdt(i,2) = vy(i)
			dxyzdt(i,3) = vz(i)
		

		  end do

		  do r = 1, 5
		  	 
		    do i = 1, cp
		      do j = 1, cp
			    if(j.NE.i) then
				  temp(j) = g*m(j)/md(i,j)**2			  		
			    end if 
		      end do 

		      do j = 1, cp
		        if(i.NE.j) then
		          coeff = temp(j)/md(i,j)			 		  
			    end if 
			    tempx(j) = (x(j)-x(i))*coeff
			    tempy(j) = (y(j)-y(i))*coeff
			    tempz(j) = (z(j)-z(i))*coeff
			  end do 

			  xacc(i) = 0d0
			  yacc(i) = 0d0
			  zacc(i) = 0d0

			  do j = 1, cp
			    xacc(i) = xacc(i)+tempx(j)
			    yacc(i) = yacc(i)+tempy(j)
			    zacc(i) = zacc(i)+tempz(j)
			  end do 
			  						
		    end do 	

			do i = 1, cp
			  dvdt(i,(r-1)*3+1) = xacc(i)
			  dvdt(i,(r-1)*3+2) = yacc(i)
			  dvdt(i,(r-1)*3+3) = zacc(i)
			end do
				  		  	
		    do i = 1, cp
		      x(i) = x(i)+vx(i)*dtmiddle+xacc(i)*dtmiddle**2/2
		      y(i) = y(i)+vy(i)*dtmiddle+yacc(i)*dtmiddle**2/2
		      z(i) = z(i)+vz(i)*dtmiddle+zacc(i)*dtmiddle**2/2
		      vx(i) = vx(i)+xacc(i)*dtmiddle
		      vy(i) = vy(i)+yacc(i)*dtmiddle
		      vz(i) = vz(i)+zacc(i)*dtmiddle
		    end do

		  end do

		  do i = 1, cp

		    dxyzdt(i,4) = dxyzdt(i,1) + 0.5d0*dt*(dvdt(i,1) + 4*dvdt(i,4) + dvdt(i,7))/6d0
		    dxyzdt(i,5) = dxyzdt(i,2) + 0.5d0*dt*(dvdt(i,2) + 4*dvdt(i,5) + dvdt(i,8))/6d0
		    dxyzdt(i,6) = dxyzdt(i,3) + 0.5d0*dt*(dvdt(i,3) + 4*dvdt(i,6) + dvdt(i,9))/6d0
		    dxyzdt(i,7) = dxyzdt(i,4) + 0.5d0*dt*(dvdt(i,7) + 4*dvdt(i,10) + dvdt(i,13))/6d0
		    dxyzdt(i,8) = dxyzdt(i,5) + 0.5d0*dt*(dvdt(i,8) + 4*dvdt(i,11) + dvdt(i,14))/6d0
		    dxyzdt(i,9) = dxyzdt(i,6) + 0.5d0*dt*(dvdt(i,9) + 4*dvdt(i,12) + dvdt(i,15))/6d0

		  end do
		  		
		  do i =1, cp

		    x(i) = xyzstart(i,1) + dt*(dxyzdt(i,1) + 4*dxyzdt(i,4) + dxyzdt(i,7))/6d0
			y(i) = xyzstart(i,2) + dt*(dxyzdt(i,2) + 4*dxyzdt(i,5) + dxyzdt(i,8))/6d0
			z(I) = xyzstart(i,3) + dt*(dxyzdt(i,3) + 4*dxyzdt(i,6) + dxyzdt(i,9))/6d0
			vx(i) = dxyzdt(i,7)
			vy(i) = dxyzdt(i,8)
			vz(i) = dxyzdt(i,9)

		  end do

		    	 open(2,file = 'Out.txt')
		  
		  do i = 1, cp
			
			do j = 1, 3

			write(2,*) xyzstart(i,j)

		  end do

		  write(2,*) x(i)
		  write(2,*) y(i)
		  write(2,*) z(i)

		  end do

		  close(2)

		  do i = 1, cp
		    do j = i+1, cp
		      md(i,j) = dsqrt((x(i)-x(j))**2+(y(i)-y(j))**2+(z(i)-z(j))**2)
			  md(j,i) = md(i,j)
		    end do 
		  end do 	

		   do i=1, cp
		    datastorage(1+3*(l*cp+i-1))=x(i)
			datastorage(2+3*(l*cp+i-1))=y(i)
			datastorage(3+3*(l*cp+i-1))=z(i)
		  end do

		  kineticenergy =0d0
		  potentialenergy =0d0
		  do i=1,cp

		  kineticenergy = kineticenergy +m(i)*(vx(i)*vx(i)+vy(i)*vy(i)+vz(i)*vz(i))/2

		  end do

		  !energyconservation
		  do i=1,cp
			
			do j=i+1, cp

			potentialenergy = potentialenergy - g*m(i)*m(j)/md(i,j)

			end do
		  
		  end do

		  deltaenergy = startenergy - (potentialenergy+kineticenergy+mergeenergy)

		  energycoeff=dsqrt(1+deltaenergy/kineticenergy)

		  do i=1, cp

          vx(i)=vx(i)*energycoeff
		  vy(i)=vy(i)*energycoeff
		  vz(i)=vz(i)*energycoeff

		  end do
		  !energyconservation

		  !impulseconservation
		  sigma =0d0

		  velocitycentermass(1)=0d0
		  velocitycentermass(2)=0d0
		  velocitycentermass(3)=0d0

		  do i=1,cp

		  sigma=sigma+m(i)
		  velocitycentermass(1)=velocitycentermass(1)+m(i)*vx(i)
		  velocitycentermass(2)=velocitycentermass(2)+m(i)*vy(i)
		  velocitycentermass(3)=velocitycentermass(3)+m(i)*vz(i)

		  end do

		  velocitycentermass(1)=velocitycentermass(1)/sigma
		  velocitycentermass(2)=velocitycentermass(2)/sigma
		  velocitycentermass(3)=velocitycentermass(3)/sigma

		  do i=1, cp

		  vx(i)=vx(i)-velocitycentermass(1)
          vy(i)=vy(i)-velocitycentermass(2)
		  vz(i)=vz(i)-velocitycentermass(3)


		  end do

		  !impulseconservation


	   end do grandloop
		

	end subroutine integrate


		
	
!qsrtit		

 
