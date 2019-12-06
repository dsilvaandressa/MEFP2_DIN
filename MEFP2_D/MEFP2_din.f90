!  MEFP2.f90 

!****************************************************************************
!
!  PROGRAM: MEFP2
!
!  PURPOSE:  Análise de chapa 2D pelo MEFP.
!****************************************************************************

    program MEFP2

    use VAR
    use leitura
    use procedimentos
    use saida
    use algebra
    
    implicit none

    call acadm() !leitura da saida do acadmesh
    
    !Pontos de integração
    pi=12
    allocate(pintegra(3,12))
    !xsi1    
    pintegra(1,:) = (/0.501426509658179d0, 0.249286745170910d0, 0.249286745170910d0, 0.873821971016996d0, 0.063089014491502d0, 0.063089014491502d0, 0.053145049844816d0, 0.310352451033785d0, 0.636502499121399d0, 0.310352451033785d0, 0.636502499121399d0, 0.053145049844816d0/)
    !xsi2
    pintegra(2,:) = (/0.249286745170910d0, 0.249286745170910d0, 0.501426509658179d0, 0.063089014491502d0, 0.063089014491502d0, 0.873821971016996d0, 0.310352451033785d0, 0.636502499121399d0, 0.053145049844816d0, 0.053145049844816d0, 0.310352451033785d0, 0.636502499121399d0/)
    !peso
    pintegra(3,:) = (/0.116786275726379d0/2.0d0, 0.116786275726379d0/2.0d0, 0.116786275726379d0/2.0d0, 0.050844906370207d0/2.0d0, 0.050844906370207/2.0d0, 0.050844906370207d0/2.0d0, 0.082851075618374d0/2.0d0, 0.082851075618374d0/2.0d0, 0.082851075618374d0/2.0d0, 0.082851075618374d0/2.0d0, 0.082851075618374d0/2.0d0, 0.082851075618374d0/2.0d0/)
    
    !definição de MFI mfi(i,j) = fi_j(pi_i)
    do p=1,12
        xsi1=pintegra(1,p)
        xsi2=pintegra(2,p)
        call fi
        MFIs_pi(p,:) = fis
    end do
       
    !Definição do tensor m para o cálculo da matriz de massa MASS
    call TM_pi() ! (pi,20,20)
    
    !tentativa inicial - inicializações
    Y=X
    fint=0.d0
    pct=200 !nº total de passos de tempo
    tol=1d-6
    err = 1d0
    it = 0
    ep = "EPD" ! "EPD" ou "EPT"
    
    !Inicializações Dinâmica
    ace=0.d0
    vel=0.d0
    dt = 0.0001d0
    t = 0.d0
    gama_newmark = 0.5d0
    beta_newmark = 0.25d0
    b1 = 1d0 !parâmetros para o cálculo de C - multiplica HESS
    b2 = 1d0 !parâmetros para o cálculo de C - multiplica matriz de massa
    
    !Início do processamento
    !Passos de tempo - pc
    
    !inicia o arquivo de saida
    call acadview_init()

    do pc=1,pct !passos de tempo
        
        t = t + dt
        fnodal = fnodal*var_t(t)
        dnodal = dnodal*var_t(t)
        
        Q = Y/(beta_newmark*dt*dt)+vel/(beta_newmark*dt) + ace*(1d0/(2d0*beta_newmark)-1d0)
        R = vel + dt*(1-gama_newmark)*ace

        print*,"Passo = ",pc,"de",pct
        it=0

nr:        do !while (err>tol) !newton raphson
            it = it + 1
            HESS=0.d0
            fint=0.d0
            MASS_G = 0.D0
            
            do el=1,nel  
                do i=1,20
                end do
                do p=1,pi
                    alfa = 1.d0
                    beta = 1.d0
                    
                    xsi1 = pintegra(1,p)
                    xsi2 = pintegra(2,p)
                    A1=0.d0
                    A0=0.d0
                    call A0_el !retorna A0(2x2) p/cada p
                    call A1_el !retorna A1(2x2) p/cada p
                    call S_el ! retorna S(2x2) p/cada p
                    call MASS_el !retorna a matriz de massa do elemento
                    
                    do gl=1,20
                        
                        alfa = 2.d0-mod(gl,2)
                        beta = int((gl+1)/2)
                        glg=inc(el,beta)*2-2+alfa

                        call DA1_ab_xsi1xsi2 !retorna DA1(2x2) p/cada alfa beta e gl
                        
                        maux1 = matmul(transpose(inv2(A0(:,:))),transpose(DA1(:,:,gl)))
                        maux2 = matmul(transpose(inv2(A0(:,:))),transpose(A1(:,:)))
                        
                        DE(:,:,gl) = 0.5d0*(matmul(matmul(maux1,A1(:,:)),inv2(A0(:,:))) + matmul(matmul(maux2,DA1(:,:,gl)),inv2(A0(:,:))))
                        
                        auxr = DE(1,1,gl)*S(1,1) + DE(1,2,gl)*S(1,2) + DE(2,1,gl)*S(2,1) + DE(2,2,gl)*S(2,2)
                          
                        fint(glg) = fint(glg) + auxr*pintegra(3,p)*det2(A0)*prop(el,4) !força interna
                        
                        ind(gl) = glg

                    end do !gl
                    
                    !forças de corpo
                    call CC_el()
                    fint(ind) = fint(ind) + matmul(TM(p,:,:),ccar(:)*var_t(t))*pintegra(3,p)*det2(A0)*prop(el,4) 
             
                    do gl=1,20
                        do gl2=1,20
                            glg=inc(el,int((gl+1)/2))*2-mod(gl,2)
                            glg2=inc(el,int((gl2+1)/2))*2-mod(gl2,2)

                            call DS_pi()
                            maux1 = matmul(transpose(inv2(A0(:,:))),transpose(DA1(:,:,gl)))
                            maux2 = matmul(transpose(inv2(A0(:,:))),transpose(DA1(:,:,gl2)))
                            
                            D2E = 0.5d0*(matmul(matmul(maux1,DA1(:,:,gl2)),inv2(A0(:,:))) + matmul(matmul(maux2,DA1(:,:,gl)),inv2(A0(:,:)))) 
                        
                            hessianinha =  DE(1,1,gl2)*DS(1,1) + DE(1,2,gl2)*DS(1,2) + DE(2,1,gl2)*DS(2,1) + DE(2,2,gl2)*DS(2,2) + &
                                D2E(1,1)*S(1,1) + D2E(1,2)*S(1,2) + D2E(2,1)*S(2,1) + D2E(2,2)*S(2,2)

                            HESS(glg,glg2) = HESS(glg,glg2) + hessianinha*det2(A0)*pintegra(3,p)*prop(el,4) 

                        end do !gl2
                    end do !gl
                end do !pi
                
                !adicionar parte dinâmica à força interna 
                !fint(ind) = fint(ind) + matmul(MASS,Y(ind))/(beta_newmark*dt*dt) - matmul(MASS,Q) + gama_newmark*matmul(C_loc,Y(ind))/(beta_newmark*dt) + matmul(C_loc,R) - matmul(C_loc,Q)*gama_newmark*dt
                MASS_G(ind,ind) = MASS
            end do !el
            
            call calc_C()
            HESS = HESS + MASS_G/(beta_newmark*dt*dt) + C_G*gama_newmark/(beta_newmark*dt)
            fint = fint + matmul(MASS_G,Y)/(beta_newmark*dt*dt) - matmul(MASS_G,Q) + gama_newmark*matmul(C_G,Y)/(beta_newmark*dt) + matmul(C_G,R) - matmul(C_G,Q)*gama_newmark*dt
            
            dbm = fnodal - fint ! (-1)*desbalanciamento mecânico           
            
            !condições de contorno
            
            do i=1,2*nnos
                if(HESS(i,i)<=0) then
                    print*,"erro"
                    pause
                else
                end if
                
                if(rest(i)==1) then
                    dbm(i)=0.d0
                    HESS(:,i) = 0.d0
                    HESS(i,:) = 0.d0
                    HESS(i,i) = 1.d0   
                else
                end if                
            end do
            
            !solução do sistema
            ip=0.d0
            call dgesv(2*nnos,1,HESS,2*nnos,ip,dbm,2*nnos,info)
            
            Y=Y+dbm
            
            err = norm2(dbm)/norm2(x)
            
            !atualização da velocidade 
            ace = Y/(beta_newmark*dt*dt) - Q
            vel = gama_newmark*Y/(beta_newmark*dt) + R - gama_newmark*dt*Q
            
            print*, "It = ",it,"err = ", err
            
            if(it>100) then
                print*,"Problema de convergencia"
                pause
            else if(err<tol) then
                exit nr
            else
            end if                       
        end do nr !fim do newton raphson
        !deslocamento
        u=y-x
        
        !calcular tensões
        tensoes=0.d0
        do el=1,nel
            do p=1,12
                xsi1 = pintegra(1,p)
                xsi2 = pintegra(2,p)
                A1=0.d0
                A0=0.d0
                call A0_el !retorna A0(2x2) p/cada p
                call A1_el !retorna A1(2x2) p/cada p
                call S_el ! retorna S(2x2) p/cada p
                
                A = matmul(A1,inv2(A0))
                
                S = matmul(matmul(A,S),transpose(A))*(1.0d0/det2(A))
                    SX(p) = S(1,1)
                    SY(p) = S(2,2)
                    SXY(p) = S(1,2)
            end do
            call tensoes_el
            tensoes(inc(el,:),1:3) = tensoes(inc(el,:),1:3)  + S_nodal
            tensoes(inc(el,:),4) = tensoes(inc(el,:),4) + 1.d0
        end do
            tensoes(:,1) = tensoes(:,1)/tensoes(:,4)
            tensoes(:,2) = tensoes(:,2)/tensoes(:,4)
            tensoes(:,3) = tensoes(:,3)/tensoes(:,4)

        call pc_acadview()
    end do !passos de carga

    
    close(20)
    end program MEFP2

