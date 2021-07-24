


module symmQRBM
  use debugging
  use matrices
 
  

  implicit none

  private :: logfile, right_translate, compute_theta, wf, compute_VarOp
  public :: s_init, s_delete,s_SR

  
  character*20 :: logfile="symm_outlog.txt"


  type :: sRBM
     integer :: Nv, alpha, dmn
     double precision ::  spin
     double complex, allocatable :: a(:),b(:)
     double complex, allocatable :: W(:,:)
  end type sRBM


contains

  subroutine s_init(self, Nv, alpha, spin)
    implicit none
    type(sRBM):: self
    integer, intent(in) :: Nv, alpha
    integer ::  status, ii, jj
    double precision, intent(in) :: spin
    double precision ::  re, im
    logical :: ca = .false., cb = .false., cW = .false.

    
    self%Nv = Nv
    self%alpha = alpha
    self%dmn = 2*alpha +Nv*alpha
    self%spin = spin

    allocate(self%a(alpha), stat=status)
    if (status .eq. 0) ca = .true.
       
    allocate(self%b(alpha), stat = status)
    if (status .eq. 0) cb = .true.

    allocate(self%W(Nv,alpha), stat = status)
    if (status .eq. 0) cW = .true.

    if (ca .and. cb .and. cW .eqv. .true.) then
       print *,"RBM succesfully initialized"
    else
       print *,"RBM allocation: FAILED"
       stop
    endif

    do ii=1,alpha
       call random_number(re)
       call random_number(im)
       self%a(ii) = complex(re*0.2-0.1,im*0.2-0.1)
    enddo

    do jj=1,alpha
       call random_number(re)
       call random_number(im)
       self%b(jj) = complex(re*0.2-0.1,im*0.2-0.1)
    enddo

    do jj=1,alpha
       do ii=1,Nv
          call random_number(re)
          call random_number(im)
          self%W(ii,jj) = complex(re*0.2-0.1,im*0.2-0.1)
       enddo
    enddo

    
    
  end subroutine s_init

  subroutine s_delete(self)
    implicit none
    type(sRBM) :: self
    deallocate(self%a, self%b, self%W)

    print *, "RBM succesfully deleted"

  end subroutine s_delete
  

   subroutine right_translate(self,S)
    implicit none
    type(sRBM), intent(in) :: self
    integer :: ii
    double complex :: S(self%Nv), S_temp(self%Nv)

    do  ii=1,self%Nv-1
       S_temp(ii+1) = S(ii)
    enddo

    S_temp(1) = S(self%Nv)
    
    S = S_temp
  
  end subroutine right_translate
  


  
  
  

    

  subroutine compute_theta(self, S, theta)
    implicit none
    type(sRBM), intent(in) :: self
    double complex, intent(in) :: S(self%Nv)
    double complex, intent(out) :: theta(self%alpha)

    theta = self%b+matmul(S, self%W)

  end subroutine compute_theta

  double complex function wf(self,S)
    implicit none
    type(sRBM), intent(in) :: self
    integer :: ss,jj
    double complex, dimension(self%Nv) :: S
    double complex, dimension(self%alpha) :: theta, act
    

    

    ! Initialize wave function
    wf = exp(sum(self%a)*sum(S)) 
    
    

       ! Add product of cosh
       do ss=1,self%Nv
       
          ! Compute theta
          call compute_theta(self, S, theta)

          !print *, "theta", double precision(sum((conjg(theta)*theta)**2))

          ! Compute activation
          act = cosh(theta)

          do jj=1,self%alpha
             !print *, "act", double precision(sum((conjg(act)*act)**2))
             wf = wf*act(jj) 
          enddo

          !Translate spin configuration
          call right_translate(self,S)
          
       enddo

       ! We translate S Nv times, so in the end it should be equal to original one


  end function wf

 
  double complex function Ising_loc(self, S, h, pbc)
    implicit none
    type(sRBM), intent(in) :: self
    integer :: ii, jj,ss, iph
    double precision, intent(in) :: h
    double complex:: frac
    double complex, dimension(self%Nv) :: S
    double complex, dimension(self%alpha) :: theta, theta_new, act, act_new
    logical, intent(in) :: pbc

    ! Initialize local energy
    Ising_loc = complex(0d0,0d0)

   

    ! Linear term
    do  ii=1,self%Nv

      
       ! Take the difference of exponential part
       frac = exp(-2.*sum(self%a)*S(ii)) 

       ! Add product of cosh
     
       do ss=0,self%Nv-1

          
          ! Compute theta
          call compute_theta(self, S, theta)

          !print *, "theta", double precision(sum((conjg(theta)*theta)**2))

          ! Compute activation
          act = cosh(theta)


          ! Copmute permutation index
          iph = mod(ii+ss-1,self%Nv)+1

          !print *, "ii", ii, "iph", iph
          
          ! Compute theta of new configuration
          theta_new = theta - 2.*self%W(iph,:)*S(iph)

         

          !print *, S(iph)-temp_spin

          !Compute  act_new
          act_new = cosh(theta_new)

          do jj=1, self%alpha
             !print *, "act", double precision(sum((conjg(act)*act)**2))
             frac = frac*act_new(jj)/act(jj) 
          enddo

          !Translate spin configuration
          call right_translate(self,S)
         
       enddo
       
       
       ! Update local energy
       Ising_loc = Ising_loc - frac

    enddo

    ! Add field
    Ising_loc = Ising_loc*h

    ! Interaction term
    do ii=1,self%Nv-1
       Ising_loc = Ising_loc -S(ii)*S(ii+1)
    end do
   

    ! Periodic boundary conditions

    if (pbc .eqv. .true.) then
       Ising_loc = Ising_loc - S(1)*S(self%Nv)
    endif

  end function Ising_loc
  
 

    

  double precision function sample_wf(self, S, N_sweeps)
    implicit none
    type(sRBM), intent(in) :: self
    integer, intent(in) :: N_sweeps
    integer :: sweep,jj, ss, site, iph
    double precision :: uu, prob
    double complex:: frac
    double complex, dimension(self%Nv):: S 
    double complex, dimension(self%alpha) :: theta, theta_new, act, act_new


    
    do sweep=1,self%Nv*N_sweeps

       ! Generate a site in (1,.., Nv)
       call random_number(uu)
       site = int(floor(self%Nv*uu)+1.)
 

    
       ! Initialize the fraction
       frac = exp(-2.*sum(self%a)*S(site)) 

       ! Add product of cosh
     
       do ss=0,self%Nv-1

  
          ! Compute theta
          call compute_theta(self, S, theta)

          !print *, "theta", double precision(sum((conjg(theta)*theta)**2))

          ! Compute activation
          act = cosh(theta)

          ! Compute index of permutation site
          iph = mod(site+ss-1,self%Nv)+1
        
          ! Compute theta of new configuration
          theta_new = theta - 2.*self%W(iph,:)*S(iph)

          !Compute act_new
          act_new = cosh(theta_new)

          do jj=1, self%alpha
             !print *, "act", double precision(sum((conjg(act)*act)**2))
             frac = frac*act_new(jj)/act(jj)
          enddo

          !Translate spin configuration
          call right_translate(self,S)
         
       enddo

       
       ! Compute the probability
       prob = min(1., dble(conjg(frac)*frac))

       !print *, "prob", prob
      
    

       ! Metropolis algorithm
       call random_number(uu)
       if (uu.le.prob) then
          S(site) = -1.*S(site)
       endif

    enddo
    
       
    ! Compute the square wave function of sampled configuration
    sample_wf = dble(conjg(wf(self,S))*wf(self,S))

  end function sample_wf

  
  subroutine compute_VarOp(self, S, O)
    implicit none
    type(sRBM), intent(in) :: self
    integer :: ii,jj,ss
    double complex, intent(in) :: S(self%Nv)
    double complex :: Oh(self%alpha), theta(self%alpha), act(self%alpha)
    double complex, intent(out) :: O(self%dmn)
    

   
    ! Copmute variational operator
    O(1:self%alpha) = sum(S)

  
    
    Oh = complex(0d0,0d0)
    do  ii=1,self%Nv
       call  compute_theta(self, S, theta)
       Oh = Oh+tanh(theta)
       call right_translate(self,S)
    enddo

 
    O(self%alpha+1:2*self%alpha) = Oh

    
    do jj=1,self%alpha
       do ii=1,self%Nv
          O(self%alpha+self%alpha+(jj-1)*self%Nv+ii)=complex(0d0,0d0)
          do ss=0,self%Nv-1
             call compute_theta(self,S,theta)
             act= tanh(theta)
             O(self%alpha+self%alpha+(jj-1)*self%Nv+ii) =  O(self%alpha+self%alpha+(jj-1)*self%Nv+ii) + S(ii)*act(jj)
             call right_translate(self,S)
          enddo
       enddo
    enddo

   end subroutine compute_VarOp

  subroutine s_SR(self, model, h, pbc, steps, N_sweeps, N_samples, E, debug)
    implicit none
    type(sRBM) :: self
    integer :: step, sample, ii,jj
    integer, intent(in) :: model, steps, N_sweeps, N_samples
    double precision :: norm, square_wf, lr, uu, threshold, start, finish
    double precision, intent(in) :: h
    double precision, intent(out) :: E(steps)
    double complex :: Eloc_av, Eloc, S(self%Nv)
    double complex, dimension(self%dmn) :: O, Odag, O_av, Odag_av, ElocOdag, ElocOdag_av, Force
    double complex, dimension(self%dmn, self%dmn) :: OdagO,OdagO_av, t_dot, Fisher, Inverse
    logical, intent(in) :: pbc , debug

    threshold = 0.00001

    open(9, file=logfile, status="unknown")

    if(debug) then
       write(9,*) "Logfile for symmetric RBM ansatz"
       write(9,*) "    "
       if(model.eq.0) then
          write(9,*) "Model: Heisenberg"
          write(9,*) "Parameters:"
          write(9,*) "Number of spins:   ", self%Nv
          write(9,*) "Alpha:   ", self%alpha
          write(9,*) "Number of variational parameters:   ", self%dmn
          write(9,*) "Periodic boundary conditions:   ", pbc
          write(9,*) "Number  of SR iterations:   ", steps
          write(9,*) "Number of MC sweep through the lattice:   ", N_sweeps
          write(9,*) "Number of samples:   ", N_samples
       elseif(model.eq.1) then
          write(9,*) "Model: Ising"
          write(9,*) "Parameters:"
          write(9,*) "Number of spins:   ", self%Nv
          write(9,*) "Alpha:   ", self%alpha
          write(9,*) "Number of variational parameters:   ", self%dmn
          write(9,*) "H_field:   ", h
          write(9,*) "Periodic boundary conditions:   ", pbc
          write(9,*) "Number  of SR iterations:   ", steps
          write(9,*) "Number of MC sweep through the lattice:   ", N_sweeps
          write(9,*) "Number of samples:   ", N_samples
       else
          write(9,*) "Unrecognize model"
          write(9,*) "Computation stopped."
          stop
       endif

       write(9,*) "*************************************************************************************************"

       write(9,*) "Starting Sochastic Reconfiguration..."
    endif
    
 

    do step=1,steps

       if(debug) then
          write(9,*) "Epoch: ", step
       end if
   

       ! Initialize tensors
       OdagO_av = complex(0d0,0d0)
       Odag_av = complex(0d0,0d0)
       O_av = complex(0d0,0d0)
       ElocOdag_av = complex(0d0,0d0)
       Eloc_av = complex(0d0,0d0)
       norm = 0.0

       !Generate the standard sequence
       do ii=1,self%Nv
          call random_number(uu)
          uu = (2.*floor(2*uu)-1.)*self%spin
          S(ii) = uu
       enddo

       if(debug) then
          write(9,*) "Square a:   ", dble(sum(conjg(self%a)*self%a))
          write(9,*) "Square b:   ", dble(sum(conjg(self%b)*self%b))
          write(9,*) "Square W:   ", dble(sum(conjg(self%W)*self%W))
       endif
       
       ! Sample sequence and probability according to thermalize MCMC
       square_wf = sample_wf(self, S, 100)

       if(debug) then
          write(9,*)   "-------------------------------------------------------------------------------------------------"
          write(9,*) "Starting sampling..."
       endif
       

       call cpu_time(start)
       do sample=1,N_samples

         
          
          ! Sample sequence and probability according to MCMC
          square_wf = sample_wf(self, S, N_sweeps)

         
          

          ! Compute O and Odag
          call compute_VarOp(self, S, O)
          Odag = conjg(O)

          !write(9,*) "Variational operator:   ", O

          ! Compute OdagO
          do jj=1,self%dmn
             do ii=1,self%dmn
                OdagO(ii,jj) = Odag(ii)*O(jj)
             enddo
          enddo

          ! Compute Eloc
          if (model .eq. 0) then
             !Eloc = Heisenberg_loc(self,S, pbc)
          elseif (model .eq. 1) then
             Eloc = Ising_loc(self, S, h, pbc)
          else
             print *, "Undefined model"
             stop
          endif
          
             

          ! Copmute ElocOdag
          ElocOdag = Eloc*Odag

          ! Update averages
          norm = norm+square_wf
          ElocOdag_av = ElocOdag_av+ElocOdag*square_wf
          Eloc_av =  Eloc_av+Eloc*square_wf

          OdagO_av = OdagO_av + OdagO*square_wf
          O_av = O_av +O*square_wf
          Odag_av = Odag_av +Odag*square_wf

       enddo

       call cpu_time(finish)

       if(debug) then
          write(9,*) "Sampling recap:"
          write(9,*) "Sampling time (s):   ", finish-start
          write(9,*) "Last sampled square_wf:   ", square_wf
          write(9,*) "Normalization constant:   ", norm
          write(9,*) "End sampling."
          write(9,*)  "-------------------------------------------------------------------------------------------------"
       endif
       
     



       !Compute t_dot
       do jj=1,self%dmn
          do ii=1,self%dmn
             t_dot(ii,jj) = Odag_av(ii)*O_av(jj)
          enddo
       enddo

       ! Compute Fisher matrix
       Fisher = OdagO_av/norm-t_dot/norm**2

       if(debug) then
          write(9,*) "Sum squared Fisher matrix:   ", sum(conjg(Fisher)*Fisher)
          !Check hermiticity
          call check_hermiticity(Fisher, threshold, logfile)
       endif

       
       
       

       ! Compute learning rate
       lr =max(0.01*(0.95**step), 0.001)       

       ! Regularize Fisher matrix
       do ii=1,self%dmn
          Fisher(ii,ii) = Fisher(ii,ii)+ 0.0001
       enddo

       !call print_matrix(Fisher)

       if(debug) then
          write(9,*) "Sum squared regularized Fisher matrix:   ", sum(conJG(Fisher)*Fisher)
       endif
       

       !print *, "Fisher", Fisher
       

       Inverse = Fisher
       call cpu_time(start)
       ! Invert Fisher matrix
       call inv_herm(Inverse)
       call cpu_time(finish)
   

       
       if(debug) then
          write(9,*) "Inversion time(s):   ", finish-start
          ! Check inverse
          call check_inverse(Fisher, Inverse, threshold, logfile)
       endif

       ! Compute the force
       Force = ElocOdag_av/norm-Eloc_av*Odag_av/norm**2

       if(debug) then
          write(9,*) "Sum squared Force:   ", sum(conjg(Force)*Force)
          write(9,*) "Sum squared Inverse:   ", sum(conjg(Inverse)*Inverse)
       end if
      

       ! Compute new weights (save in Force)
       Force = lr*matmul(Inverse, Force)


       if(debug) then
          write(9,*) "Sum squared learning force:   ", sum(conjg(Force)*Force)
          !write(9,*) "Learning force:   ", Force
       endif
       
      
      
       ! Update weights and biases
       self%a = self%a -Force(1:self%alpha)
       self%b = self%b - Force(self%alpha+1: self%alpha+self%alpha)

     

       do jj=1,self%alpha
          do ii=1,self%Nv
             self%W(ii,jj) = self%W(ii,jj) - Force(self%alpha+self%alpha+(jj-1)*self%Nv+ii)
          enddo
       enddo


       if(debug) then
          write(9,*) "Expected energy:   ", Eloc_av/norm
       endif
       

       E(step) = dble(Eloc_av/norm)

       if(debug) then
          write(9,*) "End of epoch:   ", step, "."
          write(9,*)  "*************************************************************************************************"
       endif
       

      

    enddo

  end subroutine s_SR
  
    

end module symmQRBM


  
     
  

   

