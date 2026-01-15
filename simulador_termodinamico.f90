module thermo_engine
    implicit none
    ! Coeficientes: a, b, c (Entalpia kJ/kg) | R (kJ/kg.K) | Antoine (A, B, C)
    
    ! R134a (Perfil 1 e 3)
    real, parameter :: r134_f(3) = (/200.0, 1.20, 0.0008/), r134_g(3) = (/398.5, 0.61, -0.0003/), R_134 = 0.0815
    real, parameter :: ant_134(3) = (/14.4421, 2074.63, 230.31/)
    
    ! R404A (Perfil 2)
    real, parameter :: r404_f(3) = (/200.0, 1.58, 0.002/), r404_g(3) = (/370.2, 0.41, -0.001/), R_404 = 0.0839
    real, parameter :: ant_404(3) = (/14.3443, 1969.55, 232.24/)

    

contains
    function get_hf(T, f_type)
        real, intent(in) :: T
        integer, intent(in) :: f_type
        real :: get_hf
        select case(f_type)
            case(1); get_hf = r134_f(1) + r134_f(2)*T + r134_f(3)*(T**2)
            case(2); get_hf = r404_f(1) + r404_f(2)*T + r404_f(3)*(T**2)
            case(3); get_hf = r134_f(1) + r134_f(2)*T + r134_f(3)*(T**2)
        end select
    end function get_hf

    function get_hg(T, f_type)
        real, intent(in) :: T
        integer, intent(in) :: f_type
        real :: get_hg
        select case(f_type)
            case(1); get_hg = r134_g(1) + r134_g(2)*T + r134_g(3)*(T**2)
            case(2); get_hg = r404_g(1) + r404_g(2)*T + r404_g(3)*(T**2)
            case(3); get_hg = r134_g(1) + r134_g(2)*T + r134_g(3)*(T**2)  
        end select
    end function get_hg

    function get_pressure(T, f_type)
        real, intent(in) :: T
        integer, intent(in) :: f_type
        real :: get_pressure, a, b, c
        select case(f_type)
            case(1); a=ant_134(1); b=ant_134(2); c=ant_134(3)
            case(2); a=ant_404(1); b=ant_404(2); c=ant_404(3)
            case(3); a=ant_134(1); b=ant_134(2); c=ant_134(3)
            end select
        get_pressure = exp(a - b / (T + c))
    end function get_pressure
end module thermo_engine

program refrig_sim_final
    use thermo_engine
    implicit none

    ! --- 1. DECLARAÇÕES DE VARIÁVEIS ---
    real, dimension(4) :: x, f, dx, fp
    real, dimension(4,4) :: J
    real :: Te, Tc, n_is, Q_atual, h_step, erro, p_evap, p_cond, cop1st, cop2nd
    real :: potencia_kW, energia_total_kJ, rand_val, valor_energia_total, cop_carnot
    real :: Tk_cond, Tk_evap, soma_cop2nd, contador_pontos, media_final_cop2nd, tempo_total
    real :: v_especifico, vazao_m3s, deslocamento_cm3, soma_deslocamento
    real :: Q_base, tarifa, R_gas
    integer :: it, c, t_passo, i, i_semente, tamanho_semente, perfil
    character(len=10) :: nome_fluido
    logical :: ok

    ! --- 2. CONFIGURAÇÃO INICIAL DOS PARÂMETROS ---
    energia_total_kJ = 0.0
    h_step = 1.0e-4         ! Passo para derivada numérica (Jacobiana)
    tempo_total = 100       ! Tempo total de iterações do Newton-Raphson
    soma_cop2nd = 0         ! Contador COP 2nd Law para fazer a média
    contador_pontos = 0     ! Contador de passos para fazer as médias
    soma_deslocamento = 0   ! Contador do deslocamento cm^3

    print *, "|  Selecione o perfil para a aplicacao desejada:  | "
    print *, "  -----------------------------------------------"
    print *, "|  1 - Domestico (R134a) | 2 - Comercial (R404A)  | 3- Comercial (R134a)"
    read *, perfil

    select case (perfil)
    case (1) ! DOMÉSTICO
        Te = -10.0         ! Temperatura do evaporador fixa °C
        n_is = 0.80        ! Eficiência do compressor 
        Q_base = 0.15      ! Carga térmica base de operação kW
        tarifa = 0.85      ! Tarifa de energia média residencial R$
        R_gas = 0.0815     ! Constante R do fluido R134a
        nome_fluido = "R134a"
        ! Chutes h1, h2s, h3 e m_ponto
        x = (/ 390.0, 430.0, 240.0, 0.001 /) 

    case (2) ! COMERCIAL (R404A)
        Te = -7.0         ! Temperatura do evaporador fixa °C
        n_is = 0.75        ! Eficiência do compressor 
        Q_base = 1.2       ! Carga térmica base de operação kW
        tarifa = 0.60      ! Tarifa de energia média comercial R$
        R_gas = 0.0839     ! Constante R do fluido R404a 
        nome_fluido = "R404A"
        ! Chutes h1, h2s, h3 e m_ponto
        x = (/ 380.0, 440.0, 230.0, 0.020 /) 
    case (3) ! COMERCIAL (R134a)
        Te = -10.0         ! Temperatura do evaporador fixa °C
        n_is = 0.8        ! Eficiência do compressor 
        Q_base = 1.3       ! Carga térmica base de operação kW
        tarifa = 0.60      ! Tarifa de energia média comercial R$
        R_gas = 0.0815     ! Constante R do fluido R404a 
        nome_fluido = "R134a"
        ! Chutes h1, h2s, h3 e m_ponto
        x = (/ 380.0, 440.0, 230.0, 0.020 /)     

    end select
    

    ! --- 3. INICIALIZAÇÃO DA SEMENTE ALEATÓRIA ---
    ! Usa o relógio do sistema para garantir que cada simulação seja única
    call random_seed(size=tamanho_semente)
    call system_clock(count=i_semente)
    call random_seed(put=[(i_semente + i, i=1, tamanho_semente)])

    ! --- 4. PREPARAÇÃO DE SAÍDAS (TERMINAL E ARQUIVO) ---
    open(unit=10, file='resultado_simulacao.csv', status='replace')
    write(10, '(A)') "Passo, Carga_kW, Tc_C, m_dot, COP, Potencia_kW, p_cond, p_evap, cop2nd"

    print *, "========================================================================================================"
    print *, "            SIMULADOR DE SISTEMA DE REFRIGERACAO: MONITORAMENTO DE ESTABILIDADE E CUSTO"
    print *, "========================================================================================================"
    print *, "Tempo | Carga(kW) | Tc(C) | Vazao(kg/s) | COP  | Potencia(kW) | Pressao Condensador | Pressao Evaporador"
    print *, "--------------------------------------------------------------------------------------------------------"

    ! --- 5. LOOP PRINCIPAL DE TEMPO (100 MINUTOS) ---
    do t_passo = 1, 100
        
        ! Geração de perturbações estocásticas (Vida Real)
        call random_number(rand_val)
        if (rand_val > 0.8) then 
            ! CENÁRIO DE ESTRESSE (Proporcional ao porte do sistema)
            call random_number(rand_val)
            ! A carga sobe para 150% a 200% da base (ex: porta aberta ou nova carga)
            Q_atual = Q_base * (1.5 + rand_val * 0.5)  
            
            call random_number(rand_val)
            Tc = 45.0 + (rand_val * 12.0)  ! Condensador sobrecarregado
        else
            ! CENÁRIO NOMINAL (Estável)
            call random_number(rand_val)
            ! Oscila levemente em torno da carga base (90% a 110%)
            Q_atual = Q_base * (0.9 + rand_val * 0.2)  
            
            call random_number(rand_val)
            Tc = 38.0 + (rand_val * 3.0)   ! Temperatura ambiente controlada
        end if


        ! --- 6. SOLVER NUMÉRICO (NEWTON-RAPHSON) ---
        ok = .false.
        do it = 1, 100
            call calc_F(x, f, Te, Tc, n_is, Q_atual, perfil)
            erro = sqrt(sum(f**2))

            if (erro < 1.0e-4) then
                ok = .true.
                exit
            end if

            ! PROTEÇÃO: Se o valor de entalpia ficar negativo ou absurdo, interrompe
            if (any(x /= x) .or. erro > 1.0e10) exit

            ! Cálculo da Jacobiana (Matriz de Sensibilidade)
            do c = 1, 4
                x(c) = x(c) + h_step
                call calc_F(x, fp, Te, Tc, n_is, Q_atual, perfil)
                J(:, c) = (fp - f) / h_step
                x(c) = x(c) - h_step 
            end do

            ! Sistema Linear: dx = -J^-1 * f
            call eliminacao_gauss(J, -f, dx)
            ! APLICAÇÃO DE RELAXAÇÃO PARA ESTABILIDADE
            x = x + dx 
                       


        end do



        ! --- 7. PROCESSAMENTO ECONÔMICO E SAÍDA ---
        if (ok) then

            ! Equação de Antoine para o R134a (P em kPa, T em Celsius)
            p_evap = get_pressure(Te, perfil)
            p_cond = get_pressure(Tc, perfil)

            ! Cálculo da potência instantânea do compressor (W = m_dot * delta_H)
            potencia_kW = x(4) * (x(2) - x(1))
            
            ! Integração da energia (Potência * Tempo) acumulada em kJ
            energia_total_kJ = energia_total_kJ + (potencia_kW * 60.0)
            
            ! Conversão para Custo (kWh = kJ/3600 | Tarifa média R$ 0,85)
            valor_energia_total = (energia_total_kJ / 3600.0) * tarifa

            ! COP 1st LAW
            cop1st = ((x(1)-x(3))/(x(2)-x(1)))

            ! Temperaturas em Kelvin para a 2nda Lei
            Tk_evap = Te + 273.15
            Tk_cond = Tc + 273.15

            ! COP máximo teórico (Carnot)
            cop_carnot = Tk_evap / (Tk_cond - Tk_evap)

            ! COP 2nd LAW
            cop2nd = cop1st / cop_carnot 

            ! --- CÁLCULO VOLUMÉTRICO ---
            ! R para R134a aprox. 0.0815 kJ/kg.K. Usamos Te para sucção.
            v_especifico = (R_gas * (Te + 273.15)) / p_evap
            vazao_m3s = x(4) * v_especifico
            ! Deslocamento em cm3/rev considerando 3500 RPM (padrão Brasil 60Hz)
            deslocamento_cm3 = (vazao_m3s * 1.0e6 * 60.0) / 3500.0
            
            ! Acumula para média final
            soma_deslocamento = soma_deslocamento + deslocamento_cm3

            ! ACUMULADOR PARA A MÉDIA FINAL
            soma_cop2nd = soma_cop2nd + cop2nd
            contador_pontos = contador_pontos + 1

            ! Gravação no CSV para análise em gráfico (Excel)
            write(10, '(I3, ",", F6.2, ",", F6.1, ",", F8.5, ",", F5.2, ",", F6.2, ",", F8.2, ",", F8.2, ",", F6.2)') &
            t_passo, Q_atual, Tc, x(4), cop1st, potencia_kW, p_cond, p_evap, cop2nd*100

            ! Exibição monitorada no terminal
            print '(I3, " | Q: ", F4.2, " | Tc: ", F4.1, " | m: ", F7.5, " | COP: ", F4.2, &
            & " | ", F4.2, " kW | ", F7.2, " kPa | ", F7.2 " kPa", " | COP2nd ", F5.2, "%")', &
            t_passo, Q_atual, Tc, x(4), cop1st, potencia_kW, &
            p_cond, p_evap, cop2nd*100
        
        else
            print '(I3, " | ERRO: O sistema divergiu neste ponto!")', t_passo
        end if
        
    end do

    close(10)

    ! --- 8. RELATÓRIO FINAL DE CONSUMO ---

    if (contador_pontos > 0) then
        media_final_cop2nd = (soma_cop2nd / contador_pontos) * 100.0
    else
        media_final_cop2nd = 0.0
    end if

    print *, ""
    print *, "=========================================================="
    print *, "          RELATORIO DE CERTIFICACAOO REFRIGERADOR"
    print *, "=========================================================="
    print '(A, F10.2, A)', " EFICIENCIA DE 2nda LEI MEDIA: ", media_final_cop2nd, "%"
    print '(A, F5.1, A)', " CONDICAO DE TESTE (Temperatura Evaporador):  ", Te, " C"

    ! Lógica de Classificação
    if (media_final_cop2nd >= 42.0) then
        print *, " CLASSIFICAO: [CLASSE A] - ALTA PERFORMANCE"
    else if (media_final_cop2nd >= 38.0) then
        print *, " CLASSIFICACAO: [CLASSE B]  - EFICIENTE"
    else if (media_final_cop2nd >= 33.0) then
        print *, " CLASSIFICACAO: [CLASSE C]   - BAIXA PERFORMANCE"
    else
        print *, " CLASSIFICACAO: [CLASSE D]    - MUITO BAIXA PERFORMANCE"
    end if

    print *, "----------------------------------------------------------------------------------"
    print *, "Arquivo 'resultado_simulacao.csv' gerado com sucesso!"
    print '(A, F8.4, A)', " CONSUMO TOTAL DO CICLO (100 min): ", energia_total_kJ/3600.0, " kWh"
    ! Projeção Mensal: 2 eventos de estresse por dia durante 30 dias
    print '(A, F8.2)',    " ESTIMATIVA DE CUSTO MENSAL: R$ ", valor_energia_total * 30 * 2
    print *, "----------------------------------------------------------------------------------"

    ! --- 9. DIMENSIONAMENTO DO HARDWARE (APÓS O END DO) ---
    print *, ""
    print *, "=========================================================="
    print *, "          DIMENSIONAMENTO DO COMPRESSOR"
    print *, "=========================================================="
    
    if (contador_pontos > 0) then
        deslocamento_cm3 = soma_deslocamento / contador_pontos
        print '(A, F8.2, A)', " DESLOCAMENTO MEDIO CALCULADO: ", deslocamento_cm3, " cm3/rev"
        
        if (deslocamento_cm3 <= 12.0) then
            print *, " CATEGORIA: [COMPRESSOR RESIDENCIAL]"
            print *, " SUGESTAO: Modelos de 1/4 a 1/3 HP (Ex: Embraco FFI)"
        else if (deslocamento_cm3 <= 40.0) then
            print *, " CATEGORIA: [COMPRESSOR COMERCIAL MEDIO PORTE]"
            print *, " SUGESTAO: Unidades Condensadoras (Ex: Embraco NEK)"
        else if(deslocamento_cm3 <= 90.0) then
            print *, " CATEGORIA: [COMPRESSOR COMERCIAL PESADO]"
            print *, " SUGESTAO: Unidades Condensadoras (Ex: Bitzer 4PES-15)"    
        else
            print *, " CATEGORIA: [COMPRESSOR COMERCIAL/INDUSTRIAL]"
            print *, " SUGESTAO: Sistemas de Grande Porte"
        end if
    end if
    print *, "=========================================================="
    print*, "Simulacao concluida. Aperte ENTER para sair..."
read(*,*)
contains

    subroutine calc_F(v, res, Te_val, Tc_val, n_val, Q_val, perfil)
        ! Define as equações de balanço de massa e energia do ciclo
        real, dimension(4), intent(in) :: v
        real, dimension(4), intent(out) :: res
        real, intent(in) :: Te_val, Tc_val, n_val, Q_val
        real :: h1, h2, h3, m, h2s, cp_v
        integer :: perfil
        
        ! Cálculo simplificado da entalpia isentrópica (ideal)
        select case(perfil)
        case(1); cp_v = 1.02  ! R134a
        case(2); cp_v = 1.2   ! R404A
        case(3); cp_v = 1.02  ! R134a
        end select

        !Entalpia de compressão
        h1 = v(1); h2 = v(2); h3 = v(3); m = v(4)
               
        !Entalpia Isentrópica (ideal)
        h2s = h1 + cp_v * (Tc_val - Te_val)

        ! Equações: 1.Carga Térmica | 2.Eficiência Compressor | 3.Sucção | 4.Expansão
        res(1) = m * (h1 - h3) - Q_val               
        res(2) = (h2 - h1) - ((h2s - h1) / n_val)  
        res(3) = h1 - get_hg(Te_val, perfil)
        res(4) = h3 - get_hf(Tc_val, perfil)          
    end subroutine calc_F

    subroutine eliminacao_gauss(A, b, s)
        ! Resolve sistemas lineares por Eliminação de Gauss com Retrosubstituição
        real, dimension(4,4), intent(in) :: A
        real, dimension(4), intent(in) :: b
        real, dimension(4), intent(out) :: s
        real, dimension(4,5) :: m_aug
        real :: fac
        integer :: k, i

        m_aug(1:4, 1:4) = A
        m_aug(1:4, 5)   = b

        ! Eliminação progressiva
        do k = 1, 3
            do i = k+1, 4
                fac = m_aug(i,k) / m_aug(k,k)
                m_aug(i, k:5) = m_aug(i, k:5) - fac * m_aug(k, k:5)
            end do
        end do

        ! Retrosubstituição
        s(4) = m_aug(4,5) / m_aug(4,4)
        do i = 3, 1, -1
            s(i) = (m_aug(i,5) - sum(m_aug(i, i+1:4) * s(i+1:4))) / m_aug(i,i)
        end do
    end subroutine eliminacao_gauss

end program refrig_sim_final
