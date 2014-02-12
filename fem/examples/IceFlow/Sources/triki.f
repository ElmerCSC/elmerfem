c
c   sens = 1 : tri des ki pour avoir k1 < k2 < k3
c   sens =-1 : tri des visc apres calcul avec k1 < k2 < k3 
c
c
       Subroutine triki(ki0,ki,visc,ordre,sens)

       use defgrid
       
       implicit none
       
       Real(kind=dp), Dimension(3) ::  ki0,ki
       real(kind=dp), Dimension(6) :: visc,b
       real(kind=dp) :: a
       Integer :: sens,i,j
       Integer, Dimension(3) :: ordre 
       
       

c
c    Passage pour trier les ki
c
       If (sens.EQ.1) Then
         Do i=1,3
           ki(i)=ki0(i)
           ordre(i)=i
         End Do
         Do j=2,3
          a=ki(j)
          Do i=j-1,1,-1
          If (ki(i).LE.a) Goto 20
          ki(i+1)=ki(i)
          ordre(i+1)=ordre(i)
          End Do
  20      Continue
         ki(i+1)=a
         ordre(i+1)=j
         End Do
c
c   Passage pour remettre les viscosite dans le bon ordre
c
       ElseIf (sens.EQ.-1) Then

         Do i=1,6
         b(i)=visc(i)
         End Do
         Do i=1,3
          visc(ordre(i))=b(i)
          visc(ordre(i)+3)=b(i+3)
        End Do

       Else
         Write(*,*)'triki.f : sens <> 1 ou -1' 
       Stop
       End If
       End
