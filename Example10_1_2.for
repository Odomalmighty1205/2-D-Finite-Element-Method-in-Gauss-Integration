      PROGRAM Example10_1
      
      INTEGER*4   Num_node            ! Number of nodes
      INTEGER*4   Num_element         ! Number of elements
      INTEGER*4   Num_Force           ! Number of forces
      INTEGER*4   tmp                 ! Used for calculation
      INTEGER*4   Ele_nodes(4,50)     ! Nodes in each element

      
      REAL*8  Coor(2,50)              ! Node coordinates
      REAL*8  F_node(2,50)            ! Forces on the nodes
      REAL*8  thickness               ! Thickness of the element
      REAL*8  E                       ! Modulus of elasticity
      REAL*8  Poisson                 ! Poisson ration of element
      REAL*8  D_mtx(3,3)              ! Constitutive matrix
      REAL*8  K_glob(50,50)           ! Global stiffness matrix
      REAL*8  Ke(8,8)                 ! Element stiffness matrix
      REAL*8  F_glob(50)              ! Array of global forces
      REAL*8  F_ele(8)                ! Array of element forces
      REAL*8  Dis(50)                 ! Array of global displacements
      REAL*8  Dis_ele(8)              ! Array of element displacements

      
c     Read the data
      !Open the data file and solution sheet
      OPEN(1, FILE='Data.txt', STATUS='UNKNOWN')
      OPEN(2, FILE='Solution10_1.out', STATUS='UNKNOWN')
      
      ! Read the title without storing value into any variable
      READ(1,*)

      ! Read number of element
      READ(1,*)
      READ(1,*) Num_element     
      
      ! Read number of nodes
      READ(1,*)
      READ(1,*) Num_node  
       
      ! Read the node coordinates
      READ(1,*)
      
      DO i=1,Num_node
          READ(1,*) m, (Coor(j,i),j=1,2)
      ENDDO
      
      ! Read the nodes of each element
      READ(1,*)
      DO i=1,Num_element
          READ(1,*) m,(Ele_nodes(j,i),j=1,4)
      ENDDO
      
      ! Read number of force
      READ(1,*)
      READ(1,*) Num_Force
      
      ! Read forces on the nodes
      READ(1,*)
      
      tmp=Num_Force
      DO WHILE(tmp.GT.0) 
          READ(1,*) n,(F_node(j,n),j=1,2)
          DO k=1,2
              IF(F_node(k,n).NE.0) then
                  tmp=tmp-1
              ENDIF
          ENDDO
      ENDDO
      
      ! Read thickness
      READ(1,*)
      READ(1,*) thickness
      
      ! Read the modulus of elasticity
      READ(1,*)
      READ(1,*) E
          
      ! Read the Poisson's ratio
      READ(1,*)
      READ(1,*) Poisson

c     Yielding [D] matrix
      CALL sol_D_mtx(D_mtx,E,Poisson)

c     Yielding global stiffness matrix
      CALL kglobal(Num_element,D_mtx,Coor,thickness,Ele_nodes,K_glob)
      
c     Yielding force matrix
      CALL yield_force_mtx(Num_node,F_node,F_glob)

c     Solving displacements of each nodes
      CALL sol_dis(Num_node,K_glob,F_glob,Dis)

c     Outputting the solutions
      WRITE(2,*),'The followings are the solutions of Example 10.1'
      WRITE(2,*)
      
      ! Solving (a)
      WRITE(2,*),'(a) The displacements of each node'
      DO i=1,Num_node
          WRITE(2,10),'     u',i,' =',Dis(2*i-1),' in.',
     +     '     v',i,'=',Dis(2*i),' in.'
10        FORMAT(A,I1,A,F9.5,A,A,I1,A,F9.5,A)    
      ENDDO
      WRITE(2,*)
      
      ! Solving (b)
      WRITE(2,*),'(b) The element forces'
      WRITE(2,*)
      
      DO i=1,Num_element
          CALL sol_ele_stiffness(D_mtx,Coor,thickness,i,Ele_nodes,Ke)
          CALL sol_ele_dis(i,Dis,Dis_ele,Ele_nodes)
          CALL sol_ele_force(Ke,Dis_ele,F_ele)
          WRITE(2,20),' (',i,')',' Element ',i
20        FORMAT(A,I1,A,A,I1)

          WRITE(2,30),'      fx_',Ele_nodes(1,i),
     +                '       fy_',Ele_nodes(1,i),
     +                '        fx_',Ele_nodes(2,i),
     +                '         fy_',Ele_nodes(2,i),
     +                '         fx_',Ele_nodes(3,i),
     +                '        fy_',Ele_nodes(3,i), 
     +                '        fx_',Ele_nodes(4,i),
     +                '        fy_',Ele_nodes(4,i)
30        FORMAT(4(A,I1,A,I1))      
          WRITE(2,40),(F_ele(k),k=1,8)
40        FORMAT(8(F12.2))
          WRITE(2,*)
      ENDDO
      
      STOP
      END
      
      
c
c
c     The followings are the subprograms used in the main program      
c
c    
      
*     This subroutine is aimed at solving the [D] matrix
      SUBROUTINE sol_D_mtx(D_mtx,E,Poisson)
      REAL*8  D_mtx(3,*),E,Poisson,tmp
      tmp=E/1-Poisson**2
      D_mtx(1,1)=1 ; D_mtx(2,2)=1
      D_mtx(1,2)=Poisson ; D_mtx(2,1)=Poisson
      D_mtx(3,3)=1-Poisson/2
      DO i=1,3
          DO j=1,3
              D_mtx(i,j)=tmp*D_mtx(i,j)
          ENDDO
      ENDDO
      
      RETURN
      END
      

*     This subroutine is aimed at solving [B] matrix
      SUBROUTINE q4_B_mtx(Ele_num,Ele_nodes,s,t,B_4,Coor,det_J)
      INTEGER*4   Ele_num,Ele_nodes(4,*)
      REAL*8  s,t,B_4(3,8),Coor(2,*),Xc(4),Yc(4),a,b,c,d
      REAL*8  part_x_s,part_x_t,part_y_s,part_y_t,part_N(2,50),det_J

      Xc(1)=Coor(1,Ele_nodes(1,Ele_num))
      Xc(2)=Coor(1,Ele_nodes(2,Ele_num))
      Xc(3)=Coor(1,Ele_nodes(3,Ele_num))
      Xc(4)=Coor(1,Ele_nodes(4,Ele_num))
      
      Yc(1)=Coor(2,Ele_nodes(1,Ele_num))
      Yc(2)=Coor(2,Ele_nodes(2,Ele_num))
      Yc(3)=Coor(2,Ele_nodes(3,Ele_num))
      Yc(4)=Coor(2,Ele_nodes(4,Ele_num))
      
      a=(1/4)*(Yc(1)*(s-1)+Yc(2)*(-1-s)+Yc(3)*(1+s)+Yc(4)*(1-s))
      b=(1/4)*(Yc(1)*(t-1)+Yc(2)*(1-t)+Yc(3)*(1+t)+Yc(4)*(-1-t))
      c=(1/4)*(Xc(1)*(t-1)+Xc(2)*(1-t)+Xc(3)*(1+t)+Xc(4)*(-1-t))
      d=(1/4)*(Xc(1)*(s-1)+Xc(2)*(-1-s)+Xc(3)*(1+s)+Xc(4)*(1-s))
      
      part_N(1,1)=-(1-t)/4;part_N(2,1)=-(1-s)/4
      part_N(1,2)=(1-t)/4;part_N(2,2)=-(1+s)/4
      part_N(1,3)=(1+t)/4;part_N(2,3)=(1+s)/4
      part_N(1,4)=-(1+t)/4;part_N(2,4)=(1-s)/4
      
      part_x_s=part_N(1,1)*Xc(1)+part_N(1,2)*Xc(2)+part_N(1,3)*Xc(3)
     +         +part_N(1,4)*Xc(4)
      part_x_t=part_N(2,1)*Xc(1)+part_N(2,2)*Xc(2)+part_N(2,3)*Xc(3)
     +         +part_N(2,4)*Xc(4)
      
      part_y_s=part_N(1,1)*Yc(1)+part_N(1,2)*Yc(2)+part_N(1,3)*Yc(3)
     +         +part_N(1,4)*Yc(4)
      part_y_t=part_N(2,1)*Yc(1)+part_N(2,2)*Yc(2)+part_N(2,3)*Yc(3)
     +         +part_N(2,4)*Yc(4)
      
      det_J=part_x_s*part_y_t-part_y_s*part_x_t
      
      j=1
      DO i=1,4
          B_4(1,j)=part_y_t*part_N(1,i)-part_y_s*part_N(2,i)
          B_4(2,j+1)=part_x_s*part_N(2,i)-part_x_t*part_N(1,i)
          B_4(3,j)=part_x_s*part_N(2,i)-part_x_t*part_N(1,i)
          B_4(3,j+1)=part_y_t*part_N(1,i)-part_y_s*part_N(2,i)      
          j=j+2
      ENDDO      
      
      DO i=1,3
          DO j=1,8
              B_4(i,j)=B_4(i,j)/det_J
          ENDDO
      ENDDO
      
      RETURN
      END
      

*     This subroutine is aimed at yielding stiffness matrix
      SUBROUTINE sol_ele_stiffness(D_mtx,Coor,thickness,Ele_num,
     + Ele_nodes,Ke)
      INTEGER*4   Ele_num,Ele_nodes(4,*)
      REAL*8  D_mtx(3,3),Coor(2,*),Ke(8,8)
      REAL*8  B_4(3,8),B_T(8,3),det_J,thickness
      REAL*8  tmp(8,3),sum,s,t
      
      ! Dot 1 (-0.5773,-0.5773)
      s=-0.5773
      t=-0.5773
      CALL q4_B_mtx(Ele_num,Ele_nodes,s,t,B_4,Coor,det_J)
      DO i=1,3
          DO j=1,8
              B_T(j,i)=B_4(i,j)
          ENDDO
      ENDDO
      
      sum=0
      DO i=1,8
          DO j=1,3
              DO k=1,3
                  sum=sum+B_T(i,k)*D_mtx(k,j)
              ENDDO
          tmp(i,j)=sum
          sum=0
          ENDDO
      ENDDO
      
      sum=0
      DO i=1,8
          DO j=1,8
              DO k=1,3
                  sum=sum+tmp(i,k)*B_4(k,j)
              ENDDO
          sum=sum*thickness*det_J*1*1
          Ke(i,j)=Ke(i,j)+sum
          sum=0
          ENDDO
      ENDDO
      
      
      ! Dot 2 (-0.5773,0.5773)
      s=-0.5773
      t=0.5773
      CALL q4_B_mtx(Ele_num,Ele_nodes,s,t,B_4,Coor,det_J)
      DO i=1,3
          DO j=1,8
              B_T(j,i)=B_4(i,j)
          ENDDO
      ENDDO
      
      sum=0
      DO i=1,8
          DO j=1,3
              DO k=1,3
                  sum=sum+B_T(i,k)*D_mtx(k,j)
              ENDDO
          tmp(i,j)=sum
          sum=0
          ENDDO
      ENDDO
      
      sum=0
      DO i=1,8
          DO j=1,8
              DO k=1,3
                  sum=sum+tmp(i,k)*B_4(k,j)
              ENDDO
          sum=sum*thickness*det_J*1*1
          Ke(i,j)=Ke(i,j)+sum
          sum=0
          ENDDO
      ENDDO
      
      ! Dot 3 (0.5773,-0.5773)
      s=0.5773
      t=-0.5773
      CALL q4_B_mtx(Ele_num,Ele_nodes,s,t,B_4,Coor,det_J)
      DO i=1,3
          DO j=1,8
              B_T(j,i)=B_4(i,j)
          ENDDO
      ENDDO
      
      sum=0
      DO i=1,8
          DO j=1,3
              DO k=1,3
                  sum=sum+B_T(i,k)*D_mtx(k,j)
              ENDDO
          tmp(i,j)=sum
          sum=0
          ENDDO
      ENDDO
      
      sum=0
      DO i=1,8
          DO j=1,8
              DO k=1,3
                  sum=sum+tmp(i,k)*B_4(k,j)
              ENDDO
          sum=sum*thickness*det_J*1*1
          Ke(i,j)=Ke(i,j)+sum
          sum=0
          ENDDO
      ENDDO
      

      ! Dot 4 (0.5773,0.5773)
      s=0.5773
      t=0.5773
      CALL q4_B_mtx(Ele_num,Ele_nodes,s,t,B_4,Coor,det_J)
      DO i=1,3
          DO j=1,8
              B_T(j,i)=B_4(i,j)
          ENDDO
      ENDDO
      
      sum=0
      DO i=1,8
          DO j=1,3
              DO k=1,3
                  sum=sum+B_T(i,k)*D_mtx(k,j)
              ENDDO
          tmp(i,j)=sum
          sum=0
          ENDDO
      ENDDO
      
      sum=0
      DO i=1,8
          DO j=1,8
              DO k=1,3
                  sum=sum+tmp(i,k)*B_4(k,j)
              ENDDO
          sum=sum*thickness*det_J*1*1
          Ke(i,j)=Ke(i,j)+sum
          sum=0
          ENDDO
      ENDDO
      
      RETURN
      END
     
      
*     This subroutine is aimed at yielding global stiffness matrix
      SUBROUTINE kglobal(Num_element,D_mtx,Coor,thickness,Ele_nodes,
     + K_glob)
      INTEGER*4   Num_element,Ele_nodes(4,*),tmp1,tmp2
      REAL*8  D_mtx(3,3),Ke(8,8),Coor(2,*),K_glob(50,*),thickness
      
      DO i=1,Num_element
       CALL sol_ele_stiffness(D_mtx,Coor,thickness,i,Ele_nodes,Ke)
          
          tmp1=2*Ele_nodes(1,i)-1
          tmp2=tmp1
          DO j=1,2
              DO k=1,2
                  K_glob(tmp1,tmp2)=K_glob(tmp1,tmp2)+Ke(j,k)
                  tmp2=tmp2+1
              ENDDO
              tmp1=tmp1+1
          ENDDO
          
          tmp1=2*Ele_nodes(2,i)-1
          tmp2=tmp1
          DO j=3,4
              DO k=3,4
                  K_glob(tmp1,tmp2)=K_glob(tmp1,tmp2)+Ke(j,k)
                  tmp2=tmp2+1
              ENDDO
              tmp1=tmp1+1
          ENDDO
          
          tmp1=2*Ele_nodes(3,i)-1
          tmp2=tmp1
          DO j=5,6
              DO k=5,6
                  K_glob(tmp1,tmp2)=K_glob(tmp1,tmp2)+Ke(j,k)
                  tmp2=tmp2+1
              ENDDO
              tmp1=tmp1+1
          ENDDO      
                  
          tmp1=2*Ele_nodes(4,i)-1
          tmp2=tmp1
          DO j=7,8
              DO k=7,8
                  K_glob(tmp1,tmp2)=K_glob(tmp1,tmp2)+Ke(j,k)
                  tmp2=tmp2+1
              ENDDO
              tmp1=tmp1+1
          ENDDO   
          
                    
          tmp1=2*Ele_nodes(1,i)-1
          tmp2=2*Ele_nodes(2,i)-1
          DO j=1,2
              DO k=3,4
                  K_glob(tmp1,tmp2)=K_glob(tmp1,tmp2)+Ke(j,k)
                  tmp2=tmp2+1
              ENDDO
              tmp1=tmp1+1
          ENDDO
          
          tmp1=2*Ele_nodes(2,i)-1
          tmp2=2*Ele_nodes(1,i)-1
          DO j=3,4
              DO k=1,2
                  K_glob(tmp1,tmp2)=K_glob(tmp1,tmp2)+Ke(j,k)
                  tmp2=tmp2+1
              ENDDO
              tmp1=tmp1+1
          ENDDO
          
          
          tmp1=2*Ele_nodes(1,i)-1
          tmp2=2*Ele_nodes(3,i)-1
          DO j=1,2
              DO k=5,6
                  K_glob(tmp1,tmp2)=K_glob(tmp1,tmp2)+Ke(j,k)
                  tmp2=tmp2+1
              ENDDO
              tmp1=tmp1+1
          ENDDO
          
          tmp1=2*Ele_nodes(3,i)-1
          tmp2=2*Ele_nodes(1,i)-1
          DO j=5,6
              DO k=1,2
                  K_glob(tmp1,tmp2)=K_glob(tmp1,tmp2)+Ke(j,k)
                  tmp2=tmp2+1
              ENDDO
              tmp1=tmp1+1
          ENDDO
          
          
          tmp1=2*Ele_nodes(1,i)-1
          tmp2=2*Ele_nodes(4,i)-1
          DO j=1,2
              DO k=7,8
                  K_glob(tmp1,tmp2)=K_glob(tmp1,tmp2)+Ke(j,k)
                  tmp2=tmp2+1
              ENDDO
              tmp1=tmp1+1
          ENDDO
          
          tmp1=2*Ele_nodes(4,i)-1
          tmp2=2*Ele_nodes(1,i)-1
          DO j=7,8
              DO k=1,2
                  K_glob(tmp1,tmp2)=K_glob(tmp1,tmp2)+Ke(j,k)
                  tmp2=tmp2+1
              ENDDO
              tmp1=tmp1+1
          ENDDO
          
          
          tmp1=2*Ele_nodes(2,i)-1
          tmp2=2*Ele_nodes(3,i)-1
          DO j=2,3
              DO k=5,6
                  K_glob(tmp1,tmp2)=K_glob(tmp1,tmp2)+Ke(j,k)
                  tmp2=tmp2+1
              ENDDO
              tmp1=tmp1+1
          ENDDO
          
          tmp1=2*Ele_nodes(3,i)-1
          tmp2=2*Ele_nodes(2,i)-1
          DO j=5,6
              DO k=2,3
                  K_glob(tmp1,tmp2)=K_glob(tmp1,tmp2)+Ke(j,k)
                  tmp2=tmp2+1
              ENDDO
              tmp1=tmp1+1
          ENDDO
          
          
          tmp1=2*Ele_nodes(2,i)-1
          tmp2=2*Ele_nodes(4,i)-1
          DO j=2,3
              DO k=7,8
                  K_glob(tmp1,tmp2)=K_glob(tmp1,tmp2)+Ke(j,k)
                  tmp2=tmp2+1
              ENDDO
              tmp1=tmp1+1
          ENDDO
          
          tmp1=2*Ele_nodes(4,i)-1
          tmp2=2*Ele_nodes(2,i)-1
          DO j=7,8
              DO k=2,3
                  K_glob(tmp1,tmp2)=K_glob(tmp1,tmp2)+Ke(j,k)
                  tmp2=tmp2+1
              ENDDO
              tmp1=tmp1+1
          ENDDO      

        
          tmp1=2*Ele_nodes(3,i)-1
          tmp2=2*Ele_nodes(4,i)-1
          DO j=5,6
              DO k=7,8
                  K_glob(tmp1,tmp2)=K_glob(tmp1,tmp2)+Ke(j,k)
                  tmp2=tmp2+1
              ENDDO
              tmp1=tmp1+1
          ENDDO
          
          tmp1=2*Ele_nodes(4,i)-1
          tmp2=2*Ele_nodes(3,i)-1
          DO j=7,8
              DO k=5,6
                  K_glob(tmp1,tmp2)=K_glob(tmp1,tmp2)+Ke(j,k)
                  tmp2=tmp2+1
              ENDDO
              tmp1=tmp1+1
          ENDDO
          
      ENDDO
      
      RETURN 
      END
      
*     This subroutine is aimed at yielding force array 
      SUBROUTINE yield_force_mtx(Num_node,F_node,F_glob)
      INTEGER*4   Num_node
      REAL*8  F_node(2,*),F_glob(*)
      DO i=1,Num_node
          DO j=1,2
              IF(F_node(j,i).NE.0) THEN
                  F_glob(2*i+j-2)=F_node(j,i)
              ELSE
                  F_glob(2*i+j-2)=0
              ENDIF
          ENDDO
      ENDDO
      
      RETURN
      END

*     This subroutine is aimed at solving displacement of each node
      SUBROUTINE sol_dis(Num_node,K_glob,F_glob,Dis)
      INTEGER*4   Num_node
      REAL*8  K_glob(50,*),F_glob(2*Num_node),Dis(2*Num_node)
      REAL*8  tmp(2*Num_node,4*Num_node),Inv(2*Num_node,2*Num_node)
      REAL*8  Ratio,sum
      
      DO i=1,2*Num_node
          DO j=1,2*Num_node
              tmp(i,j)=K_glob(i,j)
          ENDDO
      ENDDO

      ! Augmenting identity matrix of order n
      DO i=1,2*Num_node
          DO j=1,2*Num_node
              IF (i.EQ.j) THEN
                  tmp(i,j+2*Num_node)=1
              ELSE
                  tmp(i,j+2*Num_node)=0
              ENDIF
          ENDDO
      ENDDO
      
      ! Applying Gauss-Jordan Elimination
      DO i=1,2*Num_node
          IF (tmp(i,i).NE.0) THEN
              DO j=1,2*Num_node
                  IF(i.NE.j) THEN
                      Ratio=tmp(j,i)/tmp(i,i)
                      DO p=1,2*2*Num_node
                          tmp(j,p)=tmp(j,p)-Ratio*tmp(i,p)
                      ENDDO
                  ENDIF
              ENDDO
          ENDIF
      ENDDO
      
      ! Row operation to make principal diagonal to 1
      DO i=1,2*Num_node
          DO j=2*Num_node+1,2*2*Num_node
              tmp(i,j)=tmp(i,j)/tmp(i,i)
          ENDDO
      ENDDO
      
      DO i=1,2*Num_node
          m=2*Num_node
          DO j=1,2*Num_node
              Inv(i,j)=tmp(i,m+1)
              m=m+1
          ENDDO
      ENDDO
      
      sum=0
      DO i=1,2*Num_node
          DO j=1,2*Num_node
              sum=sum+Inv(i,j)*F_glob(j)
          ENDDO
      Dis(i)=sum
      sum=0
      ENDDO
      
      RETURN
      END
      
*     This subroutine is aimed at yielding element displacement arrary
      SUBROUTINE sol_ele_dis(Ele_num,Dis,Dis_ele,Ele_nodes)
      INTEGER*4   Ele_nodes(4,*),Ele_num,Node(4)
      REAL*8  Dis(*),Dis_ele(*)
      
      DO i=1,4
          Node(i)=Ele_nodes(i,Ele_num)
      ENDDO
      
      DO i=1,4
          Dis_ele(2*i-1)=Dis(2*Node(i)-1)
          Dis_ele(2*i)=Dis(2*Node(i))
      ENDDO
      
      RETURN
      END
      
*     This subroutine is aimed at solving element forces
      SUBROUTINE sol_ele_force(Ke,Dis_ele,F_ele)
      REAL*8  Ke(8,8),Dis_ele(8),F_ele(8)
      REAL*8  sum
      
      sum=0
      DO i=1,8
          DO j=1,8
              sum=sum+Ke(i,j)*Dis_ele(j)
          ENDDO
          F_ele(i)=sum
          sum=0
      ENDDO
      
      RETURN
      END