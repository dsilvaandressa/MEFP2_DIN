    module algebra
    implicit none
    contains
        
     !-----------------------------------------------------det2x2----------------------------------------------------!
        function det2(M) 
        
            real(8),dimension(2,2),intent(in)::M
            real(8)::det2

                det2=m(1,1)*m(2,2)-m(1,2)*m(2,1)

        end function
    
    !-----------------------------------------------------inv2x2------------------------------------------------------!
        function inv2(M)
            real(8),dimension(2,2),intent(in)::M
            real(8),dimension(2,2) :: inv2
            real(8) :: det
            
            det = (-M(1,2)*M(2,1)+M(1,1)*M(2,2))
            inv2(1,1) = M(2,2)/det
            inv2(1,2) = -M(1,2)/det
            inv2(2,1) = -M(2,1)/det
            inv2(2,2) = M(1,1)/det
        
        end function
                
        
    end module algebra