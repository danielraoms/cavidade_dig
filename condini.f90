include 'parametros.f90'		!parametros necessários para todo o main
include 'forces.f90'			!define as subroutines de forças a serem usadas
include 'verlet_list.f90'		!contém subroutine para início da lista de verlet
include 'integracao_eps.f90'		!contém subroutines para evolução das partículas e salvar eps
include 'cratera_dig.f90'		!contém subroutine para cavar uma região a partir da condição inicial

!*************************************************************
!*************************************************************
PROGRAM Colisaoverlet2D
	use parametros 		
	use integracao_eps	
	use cratera_dig	

	implicit none
	integer :: l
	integer :: a, b !contadores para célula teste (a,b)
	integer :: p, q !contadores para partículas teste (p) e comparada (q)
	integer :: cont !contador do tempo 
	integer :: cont_chains !contador para force chains
	double precision :: start, finish !auxiliares para contar tempo de simulação

	!arquivo para dados de entrada
	write(entrada,"(a,i0,a)") "dados_entrada.dat"

	!DEFININDO RAZÃO DE ASPECTO E DIMENSÕES DA CAVIDADE
	!arquivo para dados da cavidade
	write(dados_cavidade,"(a,i0,a)") "dados_cavidade.dat"

	!abre arquivo para dados da cavidade
	open(unit = 40, file=trim(dados_cavidade), status = "old")

	!puxa valor da razão de aspecto desejada, espessura_bed e espessura_parede
	read(40,*) aspect_ratio, espessura_bed, espessura_parede

	!abre arquivo para dados de entrada
	open(unit = 1, file = trim(entrada), status = "unknown")

	!arquivo para quantidades da cratera
	write(cratera, "(a,i0,a)") "cratera.dat"

	!arquivo para energias
	write(kinenergy_uniform, "(a,i0,a)") "energiacinetica_u.dat"
	write(potgenergy_uniform, "(a,i0,a)") "energiapotg_u.dat"
	write(rotenergy_uniform, "(a,i0,a)") "energiarot_u.dat"
	write(totenergy_uniform, "(a,i0,a)") "energiatotal_u.dat"
	write(kinenergy_cav, "(a,i0,a)") "energiacinetica_cav.dat"
	write(potgenergy_cav, "(a,i0,a)") "energiapotg_cv.dat"
	write(rotenergy_cav, "(a,i0,a)") "energiarot_cav.dat"
	write(totenergy_cav, "(a,i0,a)") "energiatotal_cav.dat"

	!abre arquivo para energias
	open (unit = 50, file = trim(kinenergy_uniform), status = "unknown")
	open (unit = 60, file = trim(potgenergy_uniform), status = "unknown")
	open (unit = 70, file = trim(rotenergy_uniform), status = "unknown")
	open (unit = 80, file = trim(totenergy_uniform), status = "unknown")

	!abre arquivo para energias
	open (unit = 51, file = trim(kinenergy_cav), status = "unknown")
	open (unit = 61, file = trim(potgenergy_cav), status = "unknown")
	open (unit = 71, file = trim(rotenergy_cav), status = "unknown")
	open (unit = 81, file = trim(totenergy_cav), status = "unknown")
	
call cpu_time(start)

!*******************************************************************************************************************************************************************************************
!*****PRÉ-EVOLUÇÃO DO SISTEMA (GERANDO A CONFIGURAÇÃO INICIAL)******************************************************************************************************************************

	!puxa dados do arquivo de ENTRADA para definir parâmetros de entrada da SIMULAÇÃO
	read(1,*) bw			!número de partículas na BOTTOM WALL
	read(1,*) lw			!número de partículas na LEFT WALL, sem contar a bottom wall
	read(1,*) rw			!número de partículas na RIGHT WALL, sem contar a bottom wall e a left wall
	read(1,*) contmod		!salva resultados a cada contmod passos
	read(1,*) contmodeps		!salva arquivo .eps a cada contmodeps passos
	read(1,*) contmodsafe		!salva energia a cada contmodsafe passos
	read(1,*) dt			!passo de tempo
	read(1,*) g 			!constante gravitacional
	read(1,*) xinicial		!posição x do centro da partícula mais perto da origem (0,0)
	read(1,*) yinicial		!posição y do centro da partícula mais perto da origem (0,0)
	read(1,*) raiomed		!raio médio das partículas do SISTEMA	


	!LENDO A CONDIÇÃO INICIAL A PARTIR DE UM .dat 
	!arquivo para dados da CI
	write(CI_dados,"(a,i0,a)") "dados_condini.dat"
	!abre arquivo para dados da CI
	open(unit=105, file = trim(CI_dados), status = "old")

	!LENDO DADOS DA CI
	!read(105,*)  t,N,bw,lw,rw,xinicial,yinicial,raiomed,raiomax,dt,gama_n,mi_t,gama_s,mi_roll_1
	read(105,*) t,N,bw,lw,rw,xinicial,yinicial,raiomed,raiomax,dt,g,paredes,maxIxcell,maxIycell,&
	  	    gama_n,mi_t,gama_s,mi_roll_1,mi_roll_2, E_young_p, v_poisson_p, G_shear_p

	!alocando posições e outros parâmetros físicos importantes do SISTEMA corretamente:
	!arrays que devem estar alocados a todo o momento
	allocate(r(N), m(N), inertia(N))
	allocate(xold(N), xnew(N))
	allocate(yold(N), ynew(N))
	allocate(xnewer(N), ynewer(N))
	allocate(vxold(N), vyold(N))
	allocate(vxnew(N), vynew(N))
	allocate(theta_old(N), theta_new(N))
	allocate(omega_old(N), omega_new(N))
	allocate(forcax(N), forcay(N), torque(N))
	allocate(tdl(N,N), links(N)) 

	!arrays necessários para gerar arquivos .eps
	allocate(F_elastica(N))
	allocate(velocidade_total(2,N))
	
	!definindo as dimensões dos arrays para a história do atrito
	allocate(detector_old(N,N), detector_new(N,N))
	allocate(dx_history_x(N,N), dx_history_y(N,N))
	allocate(E_young(N), v_poisson(N), G_shear(N)) 

	!arquivo para CI
	write(CI_file,"(a,i0,a)") "condini.dat"
	!abre arquivo para CI
	open(unit = 106, file = trim(CI_file), status = "old")

	!escreve o número total de partículas do SISTEMA obtidas da CI, número máximo de colunas de células de Verlet, número máximo de linhas de células de Verlet
	write(*,*) "N", N, maxIycell, maxIxcell

	!LENDO CI (PARTÍCULAS)
	do i = 1, N
		!read(106,*)	r(i),m(i),inertia(i),xold(i),xnew(i),vxold(i),vxnew(i),forcax(i),yold(i),ynew(i),&
				!vyold(i),vynew(i),forcay(i),theta_old(i),theta_new(i),omega_old(i),omega_new(i),torque(i),&
				!E_young(i),v_poisson(i),G_shear(i)
		read(106,*) r(i),m(i),inertia(i),xold(i),xnew(i),vxold(i),vxnew(i),forcax(i),yold(i),ynew(i),&
			    vyold(i),vynew(i),forcay(i),theta_old(i),theta_new(i),omega_old(i),omega_new(i),torque(i)
	end do	

	!definindo parâmetros físicos das partículas
	do i = 1, N
		E_young(i) = E_young_p
		v_poisson(i) = v_poisson_p
		G_shear(i) = G_shear_p
	end do

	!define a lista de Verlet para calcular corretamente número máximo de células de Verlet
	call lista_verlet(tdl)

	!define o modelo de rolling friction
	mi_roll = mi_roll_1

	!a tolerância para o equilíbrio de energia cinética 
	!é igual à energia cinética de uma partícula da parede quase no repouso multiplicado pelo número de partículas móveis
	toleranciakin = (N - paredes)*m(1)*0.4d0*(1.0e-1)/2.0d0

	!a tolerância para o equilíbrio de energia rotacional
	!é igual à energia rotacional de uma partícula da parede quase no repouso multiplicado pelo número de partículas móveis
	toleranciarot = (N - paredes)*m(1)*0.4d0*(1.0e-1)/2.0d0

	!escreve a tolerância para as energias cinética e rotacional
	write(*,*) "tolerancia energia cinetica", toleranciakin
	write(*,*) "tolerancia energia rotacional", toleranciarot
	
	!iniciar a Friction History 
	detector_old(:,:) = 0
	detector_new(:,:) = 0

	!determinando o tempo máximo de simulação para obter a cavidade 
	maxtime_cav = t + 100000

	tmax_u = 350000

	write(*,*) "Dados da CI puxados."
	
!**********************************************************************************************************************************************************************************************
!*****OBTENDO A CAVIDADE (CAVANDO) *****************************************************************************************************************************************
	write(*,*) "Obtendo a cavidade..."
	allocate (flag_dig(N))

	flag_dig(:) = 1

	!ainda não houve nenhum cavamento ==> ainda não houve possível uniformização da CI ==> flag_digtype = 0
	flag_digtype = 0

	!alocando array para calcular partícula mais alta da coluna correspondente - highest_of_cell(x e y, linha, coluna)
	allocate(highest_of_cell(2,maxIycell,maxIxcell))

	!determina qual linha de células possui, na coluna central de células, a partícula com altura máxima
	!isto deve ser feito de forma a cavarmos automaticamente a cavidade

		!determina a coluna central de células
30		central_column = ceiling(bw/2.0)

		!inicializando array para calcular a altura máxima
		highest_central = 0.0d0

		!itere todas as linhas de células do sistema
		do a = 2, maxIycell-1
			!tome a partícula head da célula teste (a,central_column) como a partícula teste
			p = tdl(a,central_column)
	
			do while (p .gt. 0)

				!considere somente as partículas móveis
				if (p .gt. paredes) then

					if (flag_dig(p) .gt. 0) then

						!calculando a altura máxima entre as partículas da célula (a,central_column)
						if (yold(p) .gt. highest_central) then
							highest_central = yold(p) !partícula com maior altura até agora
							highest_central_line = a  !linha de células correspondente à partícula com maior altura até agora
						end if
					end if

					p = links(p)
				else 
					p = links(p)
					continue
				end if
			end do	
		end do

	write(*,*) "highest central line", highest_central, highest_central_line, espessura_bed, espessura_parede

	!salvando imagem .eps com configuração do sistema no passo de tempo atual
	!call salva_eps_cratera(int(cont/contmodeps),ywall,N,paredes,r,xnew,ynew,0,theta_new,F_elastica,velocidade_total,flag_dig)
		
	!1. possuindo o valor da altura da CI em vista da coluna central das células de Verlet (highest_central_line), podemos então automatizar a obtenção da cavidade, CALCULADO acima
	!2. possuindo o valor da espessura do fundo da cavidade até o fundo do recipiente (espessura_bed), DADO 
	!3. possuindo o valor da espessura horizontal para minimizar influência da rugosidade da parede (espessura_parede), DADO
	!cava uma região a partir da condição inicial
	write(*,*) "before dig", highest_central_line, flag_digtype
	call cavar_cratera(highest_central_line, espessura_bed, espessura_parede, flag_dig, aspect_ratio,&
			   qtde_cels_cav, N_restante, flag_digtype) 


	!se o scheme de dig não envolve uniformização da cavidade, pule para renomeação das partículas
	if (flag_digtype .eq. 0) then
		write(*,*) "CI uniformizada pronta para ser cavada."
		go to 32
	else if (flag_digtype .eq. 1) then
		write(*,*) "CI pronta para uniformização."
		go to 31
	end if


	!SCHEME PARA UNIFORMIZAÇÃO DA CRATERA	
	!loop do tempo 
31	DO cont = floor(t/dt), tmax_u	
		!write(*,*) "Dentro loop tempo - uniformização", cont

		!calcula as energias do sistema no tempo atual 
		call energias(flag_dig)

		if (mod(cont,100) .eq. 0) then
			write(*,*) "t", cont*dt, kinetic, rotational_en 
		end if


		!escreve o tempo atual na tela
		!escreve as energias em seus respectivos arquivos .dat
		if (mod(cont,100) .eq. 0) then
			write(50,*) cont*dt, kinetic
			write(60,*) cont*dt, potentialg
			write(70,*) cont*dt, rotational_en
			write(80,*) cont*dt, totalenergy
		end if

		if (mod(cont,contmodeps) .eq. 0) then 
			!write(*,*) "cont", cont
			!itere enquanto (kinetic .LT. toleranciastop)
			if ((kinetic .LT. toleranciakin*percent_digged) .AND. (rotational_en .LT. toleranciarot*percent_digged)&
										.AND. (cont .GT. maxtime_cav)) then
					
				write(*,*) "Equilíbrio (uniformização da CI):", t, kinetic, rotational_en
					
				!A CI foi uniformizada, pule para a obtenção da cavidade
				flag_digtype = 0 !mudando flag_digtype para scheme a partir da uniformização da cratera
				go to 30
			end if
		end if

		!redefine as forças e torques do SISTEMA sofrendo apenas forças de campo para novo cálculo da iteração atual
		forcax(:) = 0.0d0
		forcay(:) = -g*m(:)
		torque(:) = 0.0d0

		do j = 1, N
			if (flag_dig(j) .eq. -1) then
				forcax(j) = 0.0d0
				forcay(j) = 0.0d0		
				torque(j) = 0.0d0
			end if
		end do

		!reinicia a lista de Verlet e reposiciona as partículas nas células
		call lista_verlet(tdl) 

		!write(*,*) "depois de refazer lista de Verlet"
		
		!itere todas as células interiores do sistema
		do b = 1, maxIxcell
		do a = 1, maxIycell

			!write(*,*) "dentro 1"
			!read(*,*)
		
			!tome a partícula head da célula teste (a,b) (Head Of Cell, Hoc; Tête de Liste, TDL) como a partícula teste
			p = tdl(a,b)

			!itere todas as partículas da célula teste (a,b) 
			do while (p .GT. 0)

				!write(*,*) "dentro 2"
				!read(*,*)

				!compare a célula atual (a,b) com ela mesma e com as células adjacentes
				do j = b-1, b+1
					!ignore as ghost cells (colunas) que seriam iteradas além das fronteiras do sistema
					if ((j .le. 0) .OR. (j .gt. maxIxcell)) then 
						cycle
					end if	
				do i = a-1, a+1
					!ignore as ghost cells (linhas) que seriam iteradas além das fronteiras do sistema
					if ((i .le. 0) .OR. (i .gt. maxIycell)) then
						cycle
					end if

				!tome a partícula head da célula (i,j) como a partícula comparada
				q = tdl(i,j)		
					
					!write(*,*) "dentro 3"
					!read(*,*)

					!compare todas as partículas da célula comparada (i,j) 
					do while (q .NE. 0)

						!write(*,*) "dentro 4"
						!read(*,*)

						!não calcule as quantidades das partículas da fronteira
						if (p .gt. paredes) then

						if (((flag_dig(p) .eq. 1) .AND. (flag_dig(q) .eq. 1)) .OR.& !duas partículas móveis não-cavadas
						   ((flag_dig(p) .eq. 1) .AND. (q .le. paredes)))     then  !uma partícula móvel e uma da parede
						   !((p .le. paredes) .AND. (flag_dig(q) .eq. 1))) then	    !uma partícula da parede e uma partícula móvel

							
						!write(*,*) "dentro 5"
						!read(*,*)

						!não compare a partícula teste com ela mesma
						!critério de colisão entre a partícula teste (p) e a partícula comparada (q)
						if ((p.ne.q).AND.((((xnew(p)-xnew(q))**2.0d0)+(ynew(p)-ynew(q))**2.0d0).LE.((r(p)+r(q))**2.0))) then

							!write(*,*) "dentro 6"
							!read(*,*)

							detector_new(p,q) = 1
			
							!Friction History Detection 
							!caso 1: início de uma colisão
							if ((detector_old(p,q) .eq. 0) .and. (detector_new(p,q) .eq. 1)) then						
								dx_history_x(p,q) = 0.0d0
								dx_history_y(p,q) = 0.0d0
							!caso 2: no meio de uma colisão
							else if ((detector_old(p,q) .eq. 1) .and. (detector_new(p,q) .eq. 1)) then
								dx_history_x(p,q) = dx_history_x(p,q)
								dx_history_y(p,q) = dx_history_y(p,q)
								!mantenha dx_history_x(p,q) e dx_history_y(p,q) unchanged
							end if
								 
							!if ((p .le. paredes) .and. (q .le. paredes)) then
								!não calcule a força resultante para qualquer partícula da parede
							!	continue
							!else
							  !calcule as forças 
							call all_forces(gama_n,mi_t,gama_s,mi_roll,E_young(p),E_young(q),v_poisson(p),v_poisson(q),&
										G_shear(p),G_shear(q),m(p),m(q),r(p),r(q),xold(p),xold(q),yold(p),yold(q),&
										vxold(p),vxold(p),vyold(p),vyold(q),omega_old(p),omega_old(q),Fx_elastica,Fy_elastica,&
										Fx_viscosa,Fy_viscosa,Fat_x,Fat_y,Fs_tangencial,T_rolling,ndx,ndy,&
										dx_history_x(p,q),dx_history_y(p,q),sinal_vrel)

							!end if

							!soma as contribuições de força que a partícula comparada (q) realiza na partícula teste (p) à força total 
							!que a partícula teste (p) sofre de todas suas partículas vizinhas
							forcax(p) = forcax(p) + Fx_elastica + Fx_viscosa + Fat_x
							forcay(p) = forcay(p) + Fy_elastica + Fy_viscosa + Fat_y
							torque(p) = torque(p) - sinal_vrel*(Fs_tangencial)*r(p) + T_rolling

							!write(*,*) "p", p, sinal_vrel, Fs_tangencial, r(p), T_rolling, torque(p)

							!write(*,*) "p", p, Fx_elastica, Fx_viscosa, Fat_x, Fy_elastica,&
								       ! Fy_viscosa, Fat_y, Fs_tangencial

							!força elástica resultante na partícula p
							F_elastica(p) = F_elastica(p) + dsqrt(Fx_elastica**2.0d0 + Fy_elastica**2.0d0)

							!reseta as contribuições que a partícula comparada (q) realiza na partícula teste (p) - já foram contabilizadas!
							Fx_elastica = 0.0d0
							Fy_elastica = 0.0d0
							Fx_viscosa = 0.0d0
							Fy_viscosa = 0.0d0
							Fs_tangencial = 0.0d0
							ndx = 0.0d0
							ndy = 0.0d0

						else
							detector_new(p,q) = 0
						end if !fecha condicional de critério de colisão

						end if

						end if

						!atualiza a partícula comparada (q) como sendo a próxima da lista de vizinhos da célula comparada (i,j)
						q = links(q)

					end do !fecha loop de partículas da célula comparada (i,j)

				end do !fecha double loop de células comparadas (i,j), vizinhas da célula teste (a,b)
				end do

				!atualiza a partícula atual (p) como sendo a próxima da lista de vizinhos da célula teste (a,b)
				p = links(p)

			end do !fecha loop de partículas da célula teste (a,b)

		end do !fecha double loop de células teste (a,b)
		end do

		!write(*,*) "depois loop partículas"

	!evolui no tempo as posições, velocidades e ângulos das partículas
	call integracao_verlet_cratera()

!!!!!!!!!!!!!!!!!!!!SETUP para POST-PROCESSING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	if (mod(cont,contmodeps) .eq. 0) then
 
		!calculando as energias do tempo atual do sistema
		call energias(flag_dig)

		!salvando imagem .eps com configuração do sistema no passo de tempo atual
		call salva_eps_cratera(int(cont/contmodeps),ywall,N,paredes,r,xnew,ynew,0,theta_new,F_elastica,velocidade_total,flag_dig)
	
	end if

		!write(*,*) "depois eps"
		!read(*,*) 

END DO

	write(*,*) "Equilíbrio (uniformização da CI):", t, kinetic, rotational_en
					
	!A CI foi uniformizada, pule para a obtenção da cavidade
	flag_digtype = 0 !mudando flag_digtype para scheme a partir da uniformização da cratera
	go to 30

	!determinando percentual de partículas que foram cavadas com relação ao total de partículas na CI (para recalcular as tolerâncias de energia)
	percent_digged = N_restante/N


	write(*,*) "Evoluindo um pouco a cavidade..."
	!evoluindo a cavidade em um passo de tempo 
	!loop do tempo 
32	DO cont = tmax_u, tmax_u + 1000	
		!write(*,*) "Dentro loop tempo - uniformização", cont

		!calcula as energias do sistema no tempo atual 
		call energias(flag_dig)

		if (mod(cont,1) .eq. 0) then
			write(*,*) "t", cont*dt, kinetic, rotational_en 
		end if


		!escreve o tempo atual na tela
		!escreve as energias em seus respectivos arquivos .dat
		if (mod(cont,1) .eq. 0) then
			write(51,*) cont*dt, kinetic
			write(61,*) cont*dt, potentialg
			write(71,*) cont*dt, rotational_en
			write(81,*) cont*dt, totalenergy
		end if

		if (mod(cont,1) .eq. 0) then 
			!write(*,*) "cont", cont
			!itere enquanto (kinetic .LT. toleranciastop)
			if ((kinetic .LT. toleranciakin*percent_digged) .AND. (rotational_en .LT. toleranciarot*percent_digged)&
										.AND. (cont .GT. maxtime_cav)) then
					
				write(*,*) "Equilíbrio (obtenção da cavidade):", t, kinetic, rotational_en
					
				!A CI foi uniformizada, pule para a obtenção da cavidade
				flag_digtype = 0 !mudando flag_digtype para scheme a partir da uniformização da cratera
				go to 35
			end if
		end if

		!redefine as forças e torques do SISTEMA sofrendo apenas forças de campo para novo cálculo da iteração atual
		forcax(:) = 0.0d0
		forcay(:) = -g*m(:)
		torque(:) = 0.0d0

		do j = 1, N
			if (flag_dig(j) .eq. -1) then
				forcax(j) = 0.0d0
				forcay(j) = 0.0d0		
				torque(j) = 0.0d0
			end if
		end do

		!reinicia a lista de Verlet e reposiciona as partículas nas células
		call lista_verlet(tdl) 

		!write(*,*) "depois de refazer lista de Verlet"
		
		!itere todas as células interiores do sistema
		do b = 1, maxIxcell
		do a = 1, maxIycell

			!write(*,*) "dentro 1"
			!read(*,*)
		
			!tome a partícula head da célula teste (a,b) (Head Of Cell, Hoc; Tête de Liste, TDL) como a partícula teste
			p = tdl(a,b)

			!itere todas as partículas da célula teste (a,b) 
			do while (p .GT. 0)

				!write(*,*) "dentro 2"
				!read(*,*)

				!compare a célula atual (a,b) com ela mesma e com as células adjacentes
				do j = b-1, b+1
					!ignore as ghost cells (colunas) que seriam iteradas além das fronteiras do sistema
					if ((j .le. 0) .OR. (j .gt. maxIxcell)) then 
						cycle
					end if	
				do i = a-1, a+1
					!ignore as ghost cells (linhas) que seriam iteradas além das fronteiras do sistema
					if ((i .le. 0) .OR. (i .gt. maxIycell)) then
						cycle
					end if

				!tome a partícula head da célula (i,j) como a partícula comparada
				q = tdl(i,j)		
					
					!write(*,*) "dentro 3"
					!read(*,*)

					!compare todas as partículas da célula comparada (i,j) 
					do while (q .NE. 0)

						!write(*,*) "dentro 4"
						!read(*,*)

						!não calcule as quantidades das partículas da fronteira
						if (p .gt. paredes) then

						if (((flag_dig(p) .eq. 1) .AND. (flag_dig(q) .eq. 1)) .OR.& !duas partículas móveis não-cavadas
						   ((flag_dig(p) .eq. 1) .AND. (q .le. paredes)))     then  !uma partícula móvel e uma da parede
						   !((p .le. paredes) .AND. (flag_dig(q) .eq. 1))) then	    !uma partícula da parede e uma partícula móvel

							
						!write(*,*) "dentro 5"
						!read(*,*)

						!não compare a partícula teste com ela mesma
						!critério de colisão entre a partícula teste (p) e a partícula comparada (q)
						if ((p.ne.q).AND.((((xnew(p)-xnew(q))**2.0d0)+(ynew(p)-ynew(q))**2.0d0).LE.((r(p)+r(q))**2.0))) then

							!write(*,*) "dentro 6"
							!read(*,*)

							detector_new(p,q) = 1
			
							!Friction History Detection 
							!caso 1: início de uma colisão
							if ((detector_old(p,q) .eq. 0) .and. (detector_new(p,q) .eq. 1)) then						
								dx_history_x(p,q) = 0.0d0
								dx_history_y(p,q) = 0.0d0
							!caso 2: no meio de uma colisão
							else if ((detector_old(p,q) .eq. 1) .and. (detector_new(p,q) .eq. 1)) then
								dx_history_x(p,q) = dx_history_x(p,q)
								dx_history_y(p,q) = dx_history_y(p,q)
								!mantenha dx_history_x(p,q) e dx_history_y(p,q) unchanged
							end if
								 
							!if ((p .le. paredes) .and. (q .le. paredes)) then
								!não calcule a força resultante para qualquer partícula da parede
							!	continue
							!else
							  !calcule as forças 
							call all_forces(gama_n,mi_t,gama_s,mi_roll,E_young(p),E_young(q),v_poisson(p),v_poisson(q),&
										G_shear(p),G_shear(q),m(p),m(q),r(p),r(q),xold(p),xold(q),yold(p),yold(q),&
										vxold(p),vxold(p),vyold(p),vyold(q),omega_old(p),omega_old(q),Fx_elastica,Fy_elastica,&
										Fx_viscosa,Fy_viscosa,Fat_x,Fat_y,Fs_tangencial,T_rolling,ndx,ndy,&
										dx_history_x(p,q),dx_history_y(p,q),sinal_vrel)

							!end if

							!soma as contribuições de força que a partícula comparada (q) realiza na partícula teste (p) à força total 
							!que a partícula teste (p) sofre de todas suas partículas vizinhas
							forcax(p) = forcax(p) + Fx_elastica + Fx_viscosa + Fat_x
							forcay(p) = forcay(p) + Fy_elastica + Fy_viscosa + Fat_y
							torque(p) = torque(p) - sinal_vrel*(Fs_tangencial)*r(p) + T_rolling

							!write(*,*) "p", p, sinal_vrel, Fs_tangencial, r(p), T_rolling, torque(p)

							!write(*,*) "p", p, Fx_elastica, Fx_viscosa, Fat_x, Fy_elastica,&
								       ! Fy_viscosa, Fat_y, Fs_tangencial

							!força elástica resultante na partícula p
							F_elastica(p) = F_elastica(p) + dsqrt(Fx_elastica**2.0d0 + Fy_elastica**2.0d0)

							!reseta as contribuições que a partícula comparada (q) realiza na partícula teste (p) - já foram contabilizadas!
							Fx_elastica = 0.0d0
							Fy_elastica = 0.0d0
							Fx_viscosa = 0.0d0
							Fy_viscosa = 0.0d0
							Fs_tangencial = 0.0d0
							ndx = 0.0d0
							ndy = 0.0d0

						else
							detector_new(p,q) = 0
						end if !fecha condicional de critério de colisão

						end if

						end if

						!atualiza a partícula comparada (q) como sendo a próxima da lista de vizinhos da célula comparada (i,j)
						q = links(q)

					end do !fecha loop de partículas da célula comparada (i,j)

				end do !fecha double loop de células comparadas (i,j), vizinhas da célula teste (a,b)
				end do

				!atualiza a partícula atual (p) como sendo a próxima da lista de vizinhos da célula teste (a,b)
				p = links(p)

			end do !fecha loop de partículas da célula teste (a,b)

		end do !fecha double loop de células teste (a,b)
		end do

		!write(*,*) "depois loop partículas"

	!evolui no tempo as posições, velocidades e ângulos das partículas
	call integracao_verlet_cratera()

!!!!!!!!!!!!!!!!!!!!SETUP para POST-PROCESSING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
		!calculando as energias do tempo atual do sistema
		call energias(flag_dig)

		!salvando imagem .eps com configuração do sistema no passo de tempo atual
		call salva_eps_cratera(int(cont/contmodeps),ywall,N,paredes,r,xnew,ynew,0,theta_new,F_elastica,velocidade_total,flag_dig)

		!write(*,*) "depois eps"
		!read(*,*) 

END DO


35	write(*,*) "Cavidade obtida."
	
		massa_media = 0.0d0
	do j = 1, N
		if (flag_dig(j) .eq. 1) then
			massa_media = massa_media + m(j)
		else 
			CONTINUE 
		end if
	end do

	massa_media = massa_media/N_restante

	!arquivos .dat com dados das partículas para evolução do colapso.
	!arquivo para dados da cavidade obtida
	write(dados_cavidade, "(a,i0,a)") "dados_out_cavidade.dat"
	
	!arquivo para quantidades da condição inicial
	write(cavidade, "(a,i0,a)") "cavidade.dat"

	open (unit = 101, file=trim(cavidade), status = "unknown")

	!escreve num arquivo .dat as quantidades das partículas na condição inicial
	open (unit = 101, file=trim(cavidade), status = "unknown")


	cont_dig_end = 0
	do i = 1, N
		!salve somente as partículas cavadas
		if (flag_dig(i) .gt. 0) then
			cont_dig_end = cont_dig_end + 1
			write(101,*) r(i),m(i),inertia(i),xold(i),xnew(i),vxold(i),vxnew(i),forcax(i),yold(i),ynew(i),&
		    		    vyold(i),vynew(i),forcay(i),theta_old(i),theta_new(i),omega_old(i),omega_new(i),torque(i)

			write(101,*) ""
		end if
	end do

	close(unit=101)	


	!escreve num arquivo .dat os parâmetros da cavidade
	open (unit = 100, file=trim(dados_cavidade), status = "unknown")

	write(100, *) cont*dt,N_restante,N_cavidade,cont_dig_end,massa_media,bw,lw,rw,xinicial,yinicial,raiomed,&
		      raiomax,dt,g,paredes,maxIxcell,maxIycell,gama_n,mi_t,gama_s,mi_roll_1,mi_roll_2, E_young,v_poisson,G_shear
	write(100,*) " "

	close(unit=100)


	
	!abre arquivo para altura máxima e packing fraction do sistema ao longo da obtenção da cavidade
	write(altura_maxima, "(a,i0,a)") "altura_max.dat"
	open (unit = 212, file=trim(altura_maxima), status = "unknown")

	!alocando arrays - POST-PROCESSING
	allocate(sum_packcell(maxIycell,maxIxcell))

	allocate(packing_fraction(maxIycell,maxIxcell))
	allocate(highest_height(maxIxcell))


	write(*,*) "Scheme de cavamento completo. Pronto para proceder para evolução do colapso"

	call cpu_time(finish)

	!escreve o tempo de simulação na tela
	write(101,*) '("Time = ", f6.3," minutes.")', (finish-start)/60, (finish-start)/3600
	write(*,*) '("Time = ", f6.3," minutes.")', (finish-start)/60, (finish-start)/3600

	write(*,*) "The simulation reached its end time."
end program Colisaoverlet2D
