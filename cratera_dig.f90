!*************************************************************
module cratera_dig
	use parametros
	implicit none
	integer :: Ixcell_hole, Iycell_hole
	integer :: left_hole, right_hole, bottom_hole, upper_hole !define as fronteiras da região a ser cavada
	double precision :: Rcell

contains
	
subroutine cavar_cratera(altura_CI, e_bed, e_parede, flag_dummy, aspect_ratio_dummy, contdig, N_resto, flag_digtype_dummy)
	use parametros
	implicit none
	integer :: pt_A, pt_B, pt_C, pt_D			 !pontos que delimitam as dimensões horizontais da cratera
	integer :: cutoff_up 					 !delimita a dimensão vertical superior da cavidade
	integer :: cutoff_2   					 !auxiliar para calcular o número restante de partículas no sistema
	integer :: aux_1, aux_2					 !auxiliares para comparativo - condicionais
	integer :: H_comp, L_comp				 !H e L previstos para o caso H/L = 1.0
	integer :: H, L						 !comprimento da cavidade, largura da cavidade
	double precision :: H_real
	integer, intent(in) :: altura_CI 			 !altura da CI vista da coluna central de células de Verlet
	integer, intent(in) :: e_bed 				 !espessura do fundo da cavidade até o fundo do recipiente - delimita a dimensão vertical inferior da cavidade
	integer, intent(in) :: e_parede 			 !espessura horizontal para minimizar influência da rugosidade da parede
	integer, intent(out) :: contdig		 		 !quantidade de partículas cavadas
	integer, intent(inout) :: N_resto 			 !número restante de partículas móveis
	integer, intent(inout) :: flag_digtype_dummy		 !especifica o tipo de cavamento a ser realizado
	integer :: contdig_2 					 !quantidade de partículas do topo retiradas para cutoff
	integer, dimension(N), intent(inout) :: flag_dummy
	double precision, intent(in) :: aspect_ratio_dummy 	 !razão de aspecto

	contdig = 0 		!contador de partículas que ainda permanecem no sistema
	contdig_2 = 0

	Rcell = 2.0d0*raiomax   !o raio da célula é o diâmetro das partículas das paredes (fronteiras) do domínio 

	!puxamos flag_digtype_dummy conforme o cutoff scheme (uniformização da CI para novo cavamento) já foi realizada ou não
	!se já foi realizada, N_resto é automaticamente puxado da call da subroutine (inout)
	if (flag_digtype_dummy .eq. 0) then !não houve uniformização
		N_resto = N
	end if

	write(*,*) "aspect ratio", aspect_ratio_dummy

	!Razão de aspecto MENOR do que 1.0
	if (aspect_ratio_dummy .lt. 1.0d0) then	
		!define L primeiro, pois será maior - esse scheme busca o melhor valor de H para o valor de L especificado
		L = altura_CI - e_bed 				!a largura da cavidade é igual à espessura da bed subtraída da altura da CI 
		H = floor(aspect_ratio_dummy*L) 		!o comprimento da cavidade é igual à razão de aspecto multiplicado pela largura da cavidade
		aux_1 = 3*H + 2 + e_parede			!valor auxiliar que é igual a 3 vezes o comprimento da cavidade + colunas das paredes esq. e dir. +  2 vezes a espessura horizontal

		!se o comprimento do recipiente é igual ou maior do que aux_1
22		if (bw .ge. aux_1) then
			!define o tamanho do bloco a ser separado em degrau/cavidade/degrau - todos com o mesmo comprimento de bloco
			if (mod((bw - 3*H),2) .eq. 0) then	!caso par
				aux_2 = (bw - 3*H)/2		
				!vestigial = 0
			else					!caso ímpar
				aux_2 = floor((bw - 3*H)/2.0)
				!vestigial = 1			
			end if

			!calculando a posição da coluna do ponto A
			pt_A = aux_2
			!calculando a posição da coluna do ponto A
			pt_B = pt_A + H
			!calculando a posição da coluna do ponto A
			pt_C = pt_B + H !+ vestigial
			!calculando a posição da coluna do ponto A
			pt_D = pt_C + H

		!se o comprimento do recipiente é menor do que aux_1
		else if (bw .lt. aux_1) then
			!decremente o valor de L atual até que caibam 3 vezes o comprimento da cavidade + colunas das paredes esq.e dir. + 2 vezes a espessura horizontal
			do while (bw .lt. aux_1)
			
				L = L - 1
				H = floor(aspect_ratio_dummy*L) 	!recalculando comprimento da cavidade

				aux_1 = 3*H + 2 + e_parede		!recalculando aux_1 para reiteração
			end do
			go to 22	!volte para o condicional que checa se o comprimento do recipiente é igual ou maior do que aux_1 - ela vai ser satisfeita
					!P.S.: se não for satisfeita, há uma inconsistência! Checar essa inconsistência!
		end if
			
	!Razão de aspecto IGUAL a 1.0
	else if (aspect_ratio_dummy .eq. 1.0d0) then
		write(*,*) "aspect ratio = 1.0"
		!comparando valores previstos para H e L - cavaremos pelo menor
		H_comp = floor((bw - 2 - 2*e_parede)/3.0)	!o comprimento da cavidade previsto é igual ao comprimento do recipiente - colunas das paredes - 2 vezes a espessura horizontal
		L_comp = altura_CI - e_bed			!a largura da cavidade prevista é igual à espessura da bed subtraída da altura da CI
		
		write(*,*) "H_comp, L_comp", H_comp, L_comp, flag_digtype_dummy

		!se H_comp > L_comp, cave em função do menor, L_comp
		!se H_comp = L_comp, é o mesmo procedimento pois não interessa se cavamos em função de H_comp ou L_comp
		if (H_comp .ge. L_comp) then
			L = L_comp 				!a largura da cavidade é igual à largura da cavidade prevista
			H = L					!razão de aspecto é igual a 1.0
			aux_1 = 3*H + 2 + 2*e_parede		!valor auxiliar que é igual a 3 vezes o comprimento da cavidade + colunas das paredes esq. e dir. + 2 vezes a espessura horizontal
			cutoff_up = 0
		!se H_comp < L_comp, cave em função do menor, H_comp	
		else if (H_comp .lt. L_comp) then
			H = H_comp				!o comprimento da cavidade é igual ao comprimento da cavidade previsto
			L = H					!razão de aspecto é igual a 1.0
			cutoff_up = L_comp - H			!valor de corte superior - acima desta linha de células, todas as partículas vão ser cavadas - nivelando o corte para melhor adequar a cavidade!
		end if
	
		write(*,*) "H_comp < L_comp", H, L, cutoff_up

		!se não há corte a ser feito OU se o corte já foi efetuado, apenas continue			
		if (cutoff_up .eq. 0) then 
			!define o tamanho do bloco a ser separado em degrau/cavidade/degrau - todos com o mesmo comprimento de bloco
			if (mod((bw - 3*H),2) .eq. 0) then	!caso par
				aux_2 = (bw - 3*H)/2		
				!vestigial = 0
			else					!caso ímpar
				aux_2 = floor((bw - 3*H)/2.0)
				!vestigial = 1			
			end if	
			
			!calculando a posição da coluna do ponto A
			pt_A = aux_2
			!calculando a posição da coluna do ponto A
			pt_B = pt_A + H
			!calculando a posição da coluna do ponto A
			pt_C = pt_B + H !+ vestigial
			!calculando a posição da coluna do ponto A
			pt_D = pt_C + H

			write(*,*) "pontos", pt_A, pt_B, pt_C, pt_D, H
			
		!se o corte ainda não foi realizado, pule para o momento em que definimos flag_dig - em que cavamos partículas propriamente
		else if (cutoff_up .gt. 0) then
			flag_digtype_dummy = 1
			write(*,*) "indicating cutoff type", flag_digtype_dummy
			go to 23		
		end if

	!Razão de aspecto MAIOR que 1.0
	else if (aspect_ratio_dummy .gt. 1.0d0) then
		write(*,*) "aspect ratio > 1.0"

		!define H primeiro, pois será maior - esse scheme busca o melhor valor de L para o valor de H especificado
		H_real = (bw - 2 - 2*e_parede)/3.0

		!define o tamanho do bloco a ser separado em degrau/cavidade/degrau - todos com o mesmo comprimento de bloco
		if (mod(H, 2) .eq. 0) then !caso par
			H = int(H_real)
		else 
			H = floor(H_real)
		end if

		L = floor(H/aspect_ratio_dummy)
		aux_2 = altura_CI - e_bed

		write(*,*) "L, aux_2", L, aux_2, altura_CI, e_bed
		read(*,*)

		!se a largura do recipiente é
25		if (aux_2 .ge. L) then
			write(*,*) "aux_2 >= L"
			cutoff_up = aux_2 - L
 
			!se não há corte a ser feito OU se o corte já foi efetuado, apenas continue			
			if (cutoff_up .eq. 0) then 
				!define o tamanho do bloco a ser separado em degrau/cavidade/degrau - todos com o mesmo comprimento de bloco
				if (mod((bw - 3*H),2) .eq. 0) then	!caso par
					aux_2 = (bw - 3*H)/2		
					!vestigial = 0
				else					!caso ímpar
					aux_2 = floor((bw - 3*H)/2.0)
					!vestigial = 1			
				end if	
			
				!calculando a posição da coluna do ponto A
				pt_A = aux_2
				!calculando a posição da coluna do ponto A
				pt_B = pt_A + H
				!calculando a posição da coluna do ponto A
				pt_C = pt_B + H !+ vestigial
				!calculando a posição da coluna do ponto A
				pt_D = pt_C + H


				flag_digtype_dummy = 0
				go to 23
			
			!se o corte ainda não foi realizado, pule para o momento em que definimos flag_dig - em que cavamos partículas propriamente
			else if (cutoff_up .gt. 0) then
				flag_digtype_dummy = 1
				go to 23		
			end if			

		!se a largura acima de e_bed é menor do que aux_1		
		else if (aux_2 .lt. L) then
			write(*,*) "aux_2 < L", aux_2

			!decremente o valor de H atual até que caiba a largura altura_CI - e_bed
			do while (aux_2 .lt. L)

				write(*,*) "aux_2 < L loop", aux_2, L
			
				H = H - 1
				L = floor(H/aspect_ratio_dummy) 	!recalculando comprimento da cavidade

				write(*,*) "H, L recalculados", H, L

			end do
			read(*,*)
			go to 25	!volte para o condicional que checa se o comprimento do recipiente é igual ou maior do que aux_1 - ela vai ser satisfeita
					!P.S.: se não for satisfeita, há uma inconsistência! Checar essa inconsistência!			
		end if
	end if

!	!definindo dimensões da cratera em função das células de Verlet (H da razão de aspecto - comprimento da cratera)
!	K_1 = int(((bw - 2) - 2*e_parede)/3.0)

!	pt_B = (1 + e_parede) + K_1		!parede esquerda DA CAVIDADE. P.S.: refere-se à COLUNA de células de Verlet
!	pt_C = pt_B + K_1 			!parede direita DA CAVIDADE. P.S.: refere-se à COLUNA de células de Verlet

!	cutoff_up = K_1 + e_bed  		!delimitação superior DA CAVIDADE - cutoff: todas as partículas acima desta linha de células deve ser cavada. P.S.: refere-se à LINHA de células de Verlet

	!identificando quais partículas serão cavadas

23	if (flag_digtype_dummy .eq. 1) then

		write(*,*) "scheme for cutoff flag = 1", flag_digtype_dummy, cutoff_up
		do j = 1, N_resto
			!determinando a célula de Verlet para a partícula j
			!Ixcell_hole = int((xold(j) - (xinicial-raiomax))/(Rcell) + 1) 
			Iycell_hole = int((yold(j) - (yinicial-raiomax))/(Rcell)) + 1 
				!write(*,*) "j, Iycell_hole", j, Iycell_hole

			!se a linha de células correspondente à partícula j é maior do que o cutoff_up, cave!
			if (Iycell_hole .gt. (altura_CI - cutoff_up)) then
				!write(*,*) "b, j, Iycell_hole", j, Iycell_hole
				flag_dummy(j) = -1
				contdig = contdig + 1	
			end if

			if ((Iycell_hole .lt. (altura_CI - cutoff_up + 5)) .and. (j .le. paredes)) then
				flag_dummy(j) = 1
				contdig = contdig - 1
			end if

		end do
		
		!pule para o final da subroutine, já cavamos e vamos estabilizar a NOVA CI (CI uniformizada para novo corte)
		go to 26
		
	
	else if (flag_digtype_dummy .eq. 0) then
		write(*,*) "scheme for cutoff flag = 0", flag_digtype_dummy, cutoff_up
		do j = 1, N_resto	
			Ixcell_hole = int((xold(j) - (xinicial-raiomax))/(Rcell)) + 1 
			Iycell_hole = int((yold(j) - (yinicial-raiomax))/(Rcell)) + 1 
				
			!se a partícula está dentro da região a ser cavada, mude seu flag_dig
			!  1. se o número da coluna de células é maior que pt_B
			!E 2. se o número da coluna de células é menor que pt_C
			!E 3. se o número da linha de células é maior que e_bed
			if ((Ixcell_hole .gt. pt_B) .AND. (Ixcell_hole .lt. pt_C) .AND. (Iycell_hole .gt. e_bed)) then
				write(*,*) "j", j, Ixcell_hole, Iycell_hole
				flag_dummy(j) = -1
				contdig = contdig + 1
			end if

			!OU SE, ALÉM dessas 3 condições,
			!	4. o número da linha de células é maior que cutoff_up
			if (Iycell_hole .ge. (altura_CI - cutoff_up)) then
				flag_dummy(j) = -1
			      	contdig = contdig + 1
			end if
		end do
	end if
				

	!P.S.: separamos os condicionais pois é interessante saber:
		!1. o volume da cavidade obtida - para o qual temos uma estimativa pois é proporcional ao número de células cavadas
		!2. o volume de fato ocupado pelas partículas - afinal, sabemos quais partículas tem flag_dig = -1


26	if (flag_digtype_dummy .eq. 0) then
  		N_resto = N_resto - contdig 
	else
		N_resto = N - contdig
	end if
	write(*,*) "26", N_resto, contdig, contdig_2

	write(*,*) "terminou dig"

	return
end subroutine cavar_cratera

end module
