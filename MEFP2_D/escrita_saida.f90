    module saida
    use VAR
    implicit none
    contains

    !-----------------------------------------------------acadview_init----------------------------------------------------!
    ! Início do arquivo de saída - fora dos loops
        subroutine acadview_init()
        
            open (20,FILE="saida_acadview.ogl",STATUS="REPLACE")
            write(20,"('#')")

            write(20,"(3(i0,3x))")nnos,nel,pct*5

            write(20,"('#')")
            
            do i=1,nnos
                write(20,"(2(f0.15,3x),4(i0,3x))")x(2*i-1),x(2*i),0,0,0,0
            end do
            
            write(20,"('#')")

            do i=1,nel
                write(20,"(12(i0,3x))")2,3,inc(i,1),inc(i,2),inc(i,3),inc(i,4),inc(i,5), &
                inc(i,6),inc(i,7),inc(i,8),inc(i,9),inc(i,10)
            end do
        
        end subroutine acadview_init
        
    !-----------------------------------------------------pc_acadview----------------------------------------------------!
    ! Saída a cada passo de carga - logo após fim do N.R.
        subroutine pc_acadview
            write(20,"('#')")
            write(20,"('Passo ',i0,'- Desloc. x')")pc

            do i=1,nnos
                write(20,"(4(f0.15,3x))")u(2*(i-1)+1),u(2*(i-1)+2),0.0d0,u(2*(i-1)+1)
            end do

            write(20,"('#')")
            write(20,"('Passo ',i0,'- Desloc. y')")pc

            do i=1,nnos
                write(20,"(4(f0.15,3x))")u(2*(i-1)+1),u(2*(i-1)+2),0.0d0,u(2*(i-1)+2)
            end do
            
            write(20,"('#')")
            write(20,"('Passo ',i0,'- Sigma X')")pc

            do i=1,nnos
                write(20,"(4(f0.15,3x))")u(2*(i-1)+1),u(2*(i-1)+2),0.0d0,tensoes(i,1)
            end do

            write(20,"('#')")
            write(20,"('Passo ',i0,'- Sigma Y')")pc

            do i=1,nnos
                write(20,"(4(f0.15,3x))")u(2*(i-1)+1),u(2*(i-1)+2),0.0d0,tensoes(i,2)
            end do

            write(20,"('#')")
            write(20,"('Passo ',i0,'- Tau XY')")pc

            do i=1,nnos
                write(20,"(4(f0.15,3x))")u(2*(i-1)+1),u(2*(i-1)+2),0.0d0,tensoes(i,3)
            end do

        end subroutine pc_acadview
    
    
    end module saida
    