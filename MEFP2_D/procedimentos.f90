    module procedimentos
    use VAR
    use t10
    use algebra
    implicit none
    contains
 
    !-----------------------------------------------------A0_el----------------------------------------------------!
        subroutine A0_el()
       
           do i=1,10
                X1(i) = X((inc(el,i))*2-1)
                X2(i) = X((inc(el,i))*2)
           end do
              
           do i=1,10
           
                call dfi1()
                call dfi2()
            
                A0(1,1) = A0(1,1) + dfis1(i)*X1(i) 
                A0(1,2) = A0(1,2) + dfis2(i)*X1(i) 
                A0(2,1) = A0(2,1) + dfis1(i)*X2(i)
                A0(2,2) = A0(2,2) + dfis2(i)*X2(i)
   
           end do 
    
        end subroutine A0_el
    !-----------------------------------------------------A1_el----------------------------------------------------!
       
        subroutine A1_el()
       
           do i=1,10
                Y1(i) = Y((inc(el,i))*2-1)
                Y2(i) = Y((inc(el,i))*2)
           end do
              
           do i=1,10
            
                call dfi1()
                call dfi2()
            
                A1(1,1) = A1(1,1) + dfis1(i)*Y1(i) 
                A1(1,2) = A1(1,2) + dfis2(i)*Y1(i)
                A1(2,1) = A1(2,1) + dfis1(i)*Y2(i)
                A1(2,2) = A1(2,2) + dfis2(i)*Y2(i)
   
           end do 
    
        end subroutine A1_el
    
    !-----------------------------------------------------DA1i----------------------------------------------------! 
        subroutine DA1_ab_xsi1xsi2 !1 pra cada p.i.
            
                call dfi1()
                call dfi2()        
        
            DA1(1,1,gl) = (2-alfa)*dfis1(beta)
            DA1(1,2,gl) = (2-alfa)*dfis2(beta)
            
            DA1(2,1,gl) = (alfa-1)*dfis1(beta)
            DA1(2,2,gl) = (alfa-1)*dfis2(beta)
        
        end subroutine DA1_ab_xsi1xsi2
    
    !-----------------------------------------------------S_pi----------------------------------------------------!    
        subroutine S_el()
            !Já foi chamada       
            !call A0_el
            !call A1_el
            
            ID=0.d0
            ID(1,1) = 1.d0
            ID(2,2) = 1.d0
            
            maux1 = matmul(transpose(inv2(A0(:,:))), transpose(A1(:,:)))
            
            E = (matmul(maux1,matmul(A1(:,:),inv2(A0(:,:))))-Id)*0.5d0
            
            v = prop(el,2)
            
            G = prop(el,1)/(2*(1+v))
            
            if (ep=="EPD") then
                S(1,1) = (G/(1-2*v))*(2*E(1,1)*(1-v)+2*E(2,2)*v)
                S(1,2) = 2*E(1,2)*G
                S(2,1) = 2*E(2,1)*G
                S(2,2) = (G/(1-2*v))*(2*E(2,2)*(1-v)+2*E(1,1)*v)
            else !EPT
                S(1,1) = (G/(1-v))*(2*E(1,1)+2*E(2,2)*v)
                S(1,2) = 2*E(1,2)*G
                S(2,1) = 2*E(2,1)*G
                S(2,2) = (G/(1-v))*(2*E(2,2)+2*E(1,1)*v)
            end if
                    
        end subroutine S_el
        
     !-----------------------------------------------------DS_pi----------------------------------------------------! 
        
        subroutine DS_pi
            ! Já foi chamada
            !call A0_el
            !call A1_el
          
            v = prop(el,2)
            
            G = prop(el,1)/(2*(1+v))
            
            if (ep=="EPD") then
                DS(1,1) = (G/(1-2*v))*(2*DE(1,1,gl)*(1-v)+2*DE(2,2,gl)*v)
                DS(1,2) = 2*DE(1,2,gl)*G
                DS(2,1) = 2*DE(2,1,gl)*G
                DS(2,2) = (G/(1-2*v))*(2*DE(2,2,gl)*(1-v)+2*DE(1,1,gl)*v)
            else !EPT
                DS(1,1) = (G/(1-v))*(2*DE(1,1,gl)+2*DE(2,2,gl)*v)
                DS(1,2) = 2*DE(1,2,gl)*G
                DS(2,1) = 2*DE(2,1,gl)*G
                DS(2,2) = (G/(1-v))*(2*DE(2,2,gl)+2*DE(1,1,gl)*v)
            end if
                    
        end subroutine DS_pi
        
    !-----------------------------------------------------tensoes_el----------------------------------------------------!
        
        subroutine tensoes_el()
            MFI = MFIs_pi
            maux3 = matmul(transpose(MFI),MFI)
            S_nodal(:,1) = matmul(transpose(MFI),SX) 
            S_nodal(:,2) = matmul(transpose(MFI),SY)
            S_nodal(:,3) = matmul(transpose(MFI),SXY)
            
            call dgesv(10,3,maux3,10,ip,S_nodal,10,info)

        end subroutine tensoes_el
        
    !-----------------------------------------------------CC_el----------------------------------------------------! 
        subroutine CC_el()
            do i=1,10
                X1(i) = X((inc(el,i))*2-1)
                X2(i) = X((inc(el,i))*2)
            end do
            
            do i=1,10
                ccar(i*2-1) = X1(i)*2 + X2(i)*1.5 !força de corpo na dir x1
                ccar(i*2) = X1(i)*2 + X2(i)*1.5   !força de corpo na dir x2
            end do        
        end subroutine CC_el
    
    !-----------------------------------------------------TM_pi----------------------------------------------------!
        subroutine TM_pi()
            do i=1,12
                xsi1 = pintegra(1,i)
                xsi2 = pintegra(2,i)
                
                call fi()
                
                do j=1,10 
                    do k=1,10
                        TM(i,2*j-1,2*k-1) = fis(j)*fis(k) 
                        TM(i,2*j,2*k-1) = 0.d0
                        
                        TM(i,2*j-1,2*k) = 0.d0
                        TM(i,2*j,2*k) = fis(j)*fis(k) 
                    end do
                end do
            end do
        
        end subroutine TM_pi
        
    !-----------------------------------------------------calc_Clocal----------------------------------------------------!
        subroutine calc_C()
            C_G = b1*HESS + b2*MASS_G      
        end subroutine calc_C
    
    !-----------------------------------------------------MASS_el----------------------------------------------------!  
        subroutine MASS_el()
        MASS=0.D0
            do i=1,12
                MASS = MASS + prop(el,3)*TM(i,:,:)*pintegra(3,i)*det2(A0)*prop(el,4)
            end do
        
        end subroutine MASS_el
        
    !-----------------------------------------------------FUNÇÕES----------------------------------------------------!
    !-----------------------------------------------------var_t(t)----------------------------------------------------!
        function var_t(t)
            real, intent(in) :: t
            real :: var_t    
            var_t = 1.d0
        end function    
        
    end module procedimentos
    
 