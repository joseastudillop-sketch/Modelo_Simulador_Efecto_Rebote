##############################################################################
# Modelo Dinamico de Simulacion de Efecto Rebote Energetico
##############################################################################
# Horizonte: 60 meses
# Escenarios: No Tax | With Tax | With Tax + Compensation
##############################################################################

using DifferentialEquations
using PlotlyJS
using Statistics

##############################################################################
# 0. CARPETA DE RESULTADOS Y UTILIDADES
##############################################################################

RESULTS_DIR_LIB = joinpath(@__DIR__, "resultados")
isdir(RESULTS_DIR_LIB) || mkpath(RESULTS_DIR_LIB)
DISPLAY_FIGURES = true # Mostrar figuras en VS Code

function rgba(hex::String, alpha::Real)
    r = parse(Int, hex[2:3], base = 16)
    g = parse(Int, hex[4:5], base = 16)
    b = parse(Int, hex[6:7], base = 16)
    return "rgba($r,$g,$b,$alpha)"
end

function save_plot(fig, filename::String; width = 980, height = 760)
    savefig(fig, joinpath(RESULTS_DIR_LIB, filename * ".png"); width = width, height = height)
    savefig(fig, joinpath(RESULTS_DIR_LIB, filename * ".pdf"); width = width, height = height)
    return nothing
end

function finite_values(vs...)
    vals = Float64[]
    for v in vs
        append!(vals, filter(isfinite, Float64.(collect(v))))
    end
    return vals
end

function auto_range(vs...; pad = 0.08, minwidth = 0.02, include_zero = false)
    vals = finite_values(vs...)
    if isempty(vals)
        return [0.0, 1.0]
    end
    mn = minimum(vals)
    mx = maximum(vals)
    if include_zero
        mn = min(mn, 0.0)
        mx = max(mx, 0.0)
    end
    width = max(mx - mn, minwidth)
    return [mn - pad * width, mx + pad * width]
end

function line_trace(x, y, name, color; dash = "solid", width = 2.4,
                    xaxis = "", yaxis = "", fillopacity = 0.0, showlegend = true)
    tr = PlotlyJS.scatter(
        x = x,
        y = y,
        mode = "lines",
        name = name,
        line = attr(color = color, width = width, dash = dash),
        showlegend = showlegend
    )
    if xaxis != ""
        tr[:xaxis] = xaxis
    end
    if yaxis != ""
        tr[:yaxis] = yaxis
    end
    if fillopacity > 0.0
        tr[:fill] = "tozeroy"
        tr[:fillcolor] = rgba(color, fillopacity)
    end
    return tr
end

function horizontal_line(x, yval, name, color; dash = "dash", width = 2,
                         xaxis = "", yaxis = "", showlegend = true)
    tr = PlotlyJS.scatter(
        x = [x[1], x[end]],
        y = [yval, yval],
        mode = "lines",
        name = name,
        line = attr(color = color, width = width, dash = dash),
        showlegend = showlegend
    )
    if xaxis != ""
        tr[:xaxis] = xaxis
    end
    if yaxis != ""
        tr[:yaxis] = yaxis
    end
    return tr
end

function axis_id(base::String, panel::Int)
    panel == 1 ? base : base * string(panel)
end

clamp01(x) = clamp(x, 0.0, 1.0)

############################################################################################################
# 1. PARAMETROS                                                                                             
############################################################################################################

############################################################################################################
#                        PARAMETROS MODIFICABLES PARA SIMULAR OTROS ESCENARIOS                             #
############################################################################################################
# El modelo trabaja en p.u. y usa bases fisicas para interpretar resultados.                               #
# Para simular otro caso, se pueden modificar los parametros señalados en esta seccion de comentarios.     #
# Mantener intacta la estructura causal del modelo.                                                        #    
#                                                                                                          #
# Variables que se pueden modificar para simular otros escenarios:                                         #
#       Horizonte temporal: t_final_months                                                                 # 
#       Sector socioeconomico: Income_base_CLP y energy_budget_share_target                                #
#       Escala del sistema: GP_base_physical_MW y FEC_base_physical_GWh_month                              # 
#       Politica de eficiencia: M_IE_policy, t_IE y tau_IE_policy                                          #
#       Rebote: eta_epsilon_S, DRE_max, indirect_rebound_max y k_indirect                                  #
#       Impuesto: P_T, t_T y tau_tax_implementation                                                        #
#       Compensacion: compensation_enabled y compensation_recycling_share                                  #
#       Inversion: pipelines iniciales, marginal_capacity_credit y tiempos de ajuste                       #
#                                                                                                          #     
# No modificar sin revisar ecuaciones:                                                                     #
#       Sustentabilidad = crecimiento GP / crecimiento S                                                   #
#       Seguridad = margen de reserva obtenido                                                             #   
#       Equidad = UCUW inicial / UCUW actual                                                               #
#       La compensacion no reduce el precio marginal, solo el gasto neto                                   #
#       GP crece por inversion endogena y retiros                                                          #
############################################################################################################

p_lib = (
    # --- control de simulacion ---
    t_final_months    = 60.0,        # Horizonte total
    dt_save_months    = 0.1,         # Paso guardado

    # --- bases internas p.u. ---
    S_base            = 100.0,       # Demanda util base
    epsilon_0         = 1.0,         # Eficiencia base
    P_E_base          = 1.0,         # Precio base

    # --- referencias fisicas ---
    FEC_base_physical_GWh_month = 7000.0,  # FEC base fisico (Dato generación y demanda a escala Pais - Coordinador)
    P_E_base_physical_USD_MWh   = 100.0,   # Precio fisico base de energía (referencia conservadora de precio energía a nivel generación, y
                                           # CNE publica precios de nudo promedio como referencia regulada del componente energía y potencia)
    GP_base_physical_MW         = 39200.0, # Capacidad base SEN (capacidad instalada de referencia del Sistema Eléctrico Nacional)

    # --- sector socioeconomico ---
    Income_base_CLP_vulnerable = 644000.0,  # Ingreso sector vulnerable
    Income_base_CLP_middle     = 1250000.0, # Ingreso sector medio
    Income_base_CLP_high       = 3000000.0, # Ingreso sector alto
    Income_base_CLP            = 644000.0,  # Ingreso simulado para defensa de modelo
    Income_base                = 1000.0,    # Ingreso interno
    energy_budget_share_target = 0.10,      # Umbral gasto energia
    EnergyBudget_base          = 100.0,     # Presupuesto energia base

    # --- trayectoria de eficiencia ---
    M_IE_policy       = 0.10,        # Mejora de eficiencia objetivo
    t_IE              = 2.0,         # Mes de inicio eficiencia
    tau_IE_policy     = 8.0,         # Velocidad politica EE  
    tau_IE            = 6.0,         # Ajuste eficiencia        
    IE_max            = 0.22,        # Techo eficiencia
    k_GP_to_IE        = 0.018,       # GP a eficiencia
    k_price_to_IE     = 0.018,       # Precio a eficiencia
    k_reserve_to_IE   = 0.012,       # Reserva a eficiencia

    # --- rebote directo ---
    eta_epsilon_S      = 1.20,       # Elasticidad eficiencia         
    DRE_max            = 0.130,      # Techo rebote directo             
    t_DRE              = 6.0,        # Mes de inicio rebote
    tau_DRE_activation = 3.0,        # Activacion rebote
    tau_D              = 7.0,        # Retardo rebote
    N_delay            = 3,          # Orden retardo

    # --- impuesto correctivo ---
    P_T                     = 0.0,   # Tasa impuesto
    t_T                     = 12.0,  # Mes de inicio de impuesto
    tau_tax_implementation  = 5.0,   # Rampa impuesto     
    tax_demand_pass_through = 0.45,   # Traspaso demanda

    # --- compensacion focalizada ---
    compensation_enabled = 0.0,            # Activa compensacion
    t_compensation = 18.0,                 # Mes de inicio de compensacion
    tau_compensation_implementation = 5.0, # Rampa compensacion     
    compensation_recycling_share = 0.60,   # Fraccion reciclaje
    compensation_cap_share = 0.28,         # Tope compensacion
    compensation_full_duration = 32.0,     # Duracion apoyo fuerte
    compensation_phaseout_tau = 12.0,      # Velocidad retiro
    compensation_phaseout_fraction = 0.40, # Reduccion final
    compensation_trigger_width = 0.35,     # Suavizado activacion
    compensation_trigger_floor = 0.55,     # Piso activacion

    # --- demanda macro ---
    eta_P              = -0.3,       # Elasticidad precio
    tau_Q              = 4.0,        # Ajuste demanda
    demand_growth_rate = 0.00045,    # Crecimiento demanda (permite que el bloque de Olsina no quede plano y tenga dinamica)

    # --- socioeconomico y rebote indirecto ---
    savings_ref_months    = 12.0,    # Normalizacion CMS
    indirect_rebound_max  = 0.045,   # Techo rebote indirecto
    k_indirect            = 18.0,    # Sensibilidad IR
    tau_IR                = 6.0,     # Ajuste IR
    beta_budget           = 1.0,     # Peso restriccion presupuestaria
    beta_stress           = 0.8,     # Peso estres 
    budget_pressure_scale = 1.2,     # Escala presion presupuestaria
    budget_constraint_max = 0.32,    # Techo restriccion presupuestaria
    tau_BC                = 5.0,     # Ajuste restriccion presupuestaria
    monthly_slack_weight  = 0.45,    # Peso holgura mensual 
    monthly_stress_weight = 0.20,    # Peso estres mensual

    # --- reserva ---
    reserve_ratio        = 0.15,     # Reserva requerida
    reserve_demand_floor = 0.965,    # Piso demanda reserva
    reserve_floor_blend  = 0.35,     # Mezcla piso reserva

    # --- costos y precios ---
    RRR                       = 0.08,  # Retorno requerido
    variable_cost_base        = 0.98,  # Costo variable base
    expected_fuel_price_index = 1.04,  # Precio combustible
    fixed_cost_adder          = 0.0,   # Costo fijo adicional
    equil_reserve_add         = 0.0,   # Ajuste reserva costo
    k_pexp_to_prodcost        = 0.04,  # Precio esperado a costo
    tau_energy_price          = 7.0,   # Ajuste precio energia
    tau_expect_price          = 7.0,   # Ajuste precio esperado
    w_spot_to_energy_price    = 0.40,  # Peso spot
    w_expected_to_energy_price = 0.60, # Peso esperado
    w_energy_price_to_expect   = 0.52, # Energia a expectativa
    w_expected_cost_to_expect  = 0.32, # Costo a expectativa
    w_fuel_price_to_expect     = 0.16, # Combustible a expectativa
    w_anchor_to_P_E_base       = 0.35, # Ancla precio base
    expectation_scarcity_gain  = 0.42, # Prima escasez
    surplus_discount_scale     = 0.080, # Descuento excedente
    max_surplus_discount       = 0.075, # Techo descuento
    tax_investment_wedge       = 0.22,  # Brecha impuesto inversion
    tax_market_adaptation_gain = 0.055, # Adaptacion mercado
    tax_expansion_brake_k      = 16.0,  # Freno inversion impuesto

    # --- inversion y Generation Park ---
    tau_ID_up                    = 7.0,             # Ajuste ID subida
    tau_ID_down                  = 5.0,             # Ajuste ID bajada
    tau_construction             = 12.0,            # Tiempo construccion
    tau_gap                      = 10.0,            # Cierre brecha
    asset_life_months            = 300.0,           # Vida util activos
    pipeline_credit              = 0.55,            # Credito pipeline
    investment_forecast_horizon  = 9.0,             # Horizonte inversion
    over_reserve_brake_k         = 38.0,            # Freno sobre reserva
    over_reserve_brake_threshold = 1.18,            # Umbral freno reserva
    marginal_capacity_credit     = 0.55,            # Capacidad efectiva
    initial_approval_pipeline_MW = 3500.0,          # Pipeline aprobacion
    initial_construction_pipeline_MW = 6500.0,      # Pipeline construccion
    capacity_gap_response_gain   = 0.95,            # Respuesta brecha              
    demand_growth_capacity_gain  = 1.40,            # Demanda a capacidad
    scarcity_capacity_gain       = 0.15,            # Escasez a capacidad           
    investment_pi_threshold      = 0.98,            # Umbral PI
    investment_pi_width          = 0.12,            # Ancho puerta PI

    # --- precio spot y PI ---
    VOLL              = 8.0,        # Precio escasez max
    price_spike_scale = 2.0,        # Escala spike
    k_spike           = 3.0,        # Sensibilidad spike
    invest_sat        = 1.75,       # Saturacion inversion
    pi_slope          = 2.4,        # Pendiente PI alta
    pi_decay_slope    = 3.0,        # Pendiente PI baja
    pi_decay_floor    = 0.02,       # Piso decaimiento PI

    # --- planificacion ---
    planning_reserve_ratio     = 0.150, # Reserva planificada
    demand_planning_weight     = 0.15,  # Peso demanda base
    peak_demand_elasticity     = 1.35,  # Elasticidad punta
    expectation_demand_gain    = 0.17,  # Prima demanda esperada
    demand_investment_gate_gain = 0.90  # Puerta demanda inversion
)
p_lib_tax = merge(p_lib, (; P_T = 0.10, tax_expansion_brake_k = 18.0)) # Escenario con impuesto

p_lib_comp = merge(p_lib, (; 
    P_T = 0.10, # Tasa impuesto
    compensation_enabled = 1.0, # Activa la compensacion
    compensation_recycling_share = 0.60, # Fraccion reciclaje
    tax_expansion_brake_k = 14.0
)) # Escenario impuesto y compensacion

##############################################################################
# 2. FUNCIONES AUXILIARES DEL MODELO
##############################################################################

function policy_ramp(t, t0, tau)                    # Calcula una activación progresiva entre 0 y 1.
    if t < t0                                       # Devuelve el grado de activación de la política
        return 0.0                                  # Sirve para que eficiencia, impuesto, rebote y compensación no aparezcan de golpe.
    elseif tau <= 0.0
        return 1.0
    else
        return 1.0 - exp(-(t - t0) / tau)
    end
end

function scarcity_component(reserve_margin, p)                                              # Calcula una señal de escasez cuando la reserva cae bajo el nivel requerido.
    if reserve_margin >= p.reserve_ratio                                                    # Devuelve un componente que presiona al alza el precio spot.
        return 0.0                                                                          # Sirve para conectar seguridad de suministro con precios eléctricos.
    elseif reserve_margin >= 0.0
        x = (p.reserve_ratio - reserve_margin) / max(p.reserve_ratio, 1e-6)
        return p.price_spike_scale * x^2
    else
        x = -reserve_margin / max(p.reserve_ratio, 1e-6)
        return p.price_spike_scale * exp(p.k_spike * x)
    end
end

function surplus_discount_component(reserve_margin, p)                                         # Calcula un descuento cuando existe exceso de reserva.
    if reserve_margin <= p.reserve_ratio                                                       # Devuelve una reducción aplicada al precio spot.                 
        return 0.0                                                                             # Sirve para representar que, cuando hay mucha capacidad disponible, el precio tiende a moderarse.    
    else
        x = (reserve_margin - p.reserve_ratio) / max(p.reserve_ratio, 1e-6)
        return min(p.max_surplus_discount, p.surplus_discount_scale * x)
    end
end

function spot_price_liberalized(production_cost, reserve_margin, p)                         # Calcula el precio spot considerando costo, escasez y descuento por exceso de capacidad.
    scarcity = scarcity_component(reserve_margin, p)                                        # Devuelve precio spot, señal de escasez y descuento por excedente.
    discount = surplus_discount_component(reserve_margin, p)                                # Sirve como parte del bloque de mercado eléctrico inspirado en la lógica de Olsina.
    spot = production_cost * (1.0 - discount) + scarcity
    return clamp(spot, 0.05, p.VOLL), scarcity, discount
end

function investment_multiplier(PI, p)                                                               # Calcula cuánto se amplifica o reduce la inversión según la rentabilidad esperada.
    if PI >= 1.0                                                                                    # Devuelve un multiplicador de inversión.
        return 1.0 + (p.invest_sat - 1.0) * (1.0 - exp(-p.pi_slope * (PI - 1.0)))                   # Sirve para que la inversión no sea fija, sino endógena.
    else
        return max(0.0, exp(p.pi_decay_slope * (PI - 1.0)) - p.pi_decay_floor)
    end
end

function over_reserve_brake(reserve_margin, p)                                                      # Calcula un freno a la inversión cuando existe sobrecapacidad.
    excess = max(0.0, reserve_margin - p.over_reserve_brake_threshold * p.reserve_ratio)            # Devuelve un factor entre 0 y 1.    
    return exp(-p.over_reserve_brake_k * excess)                                                    # Sirve para evitar que el parque generador crezca sin control. 
end

##### Funciones auxiliares de peso (weight)#####
function normalize_price_weights(p)                                                             # Esta función normaliza los pesos usados para calcular el precio objetivo de energía.
    w_total = p.w_spot_to_energy_price + p.w_expected_to_energy_price                           # Devuelven pesos ajustados en proporciones que suman 1.   
    return w_total <= 0.0 ? (0.5, 0.5) :                                                        # Sirven para construir el precio de energía y el precio esperado.   
           (p.w_spot_to_energy_price / w_total, p.w_expected_to_energy_price / w_total)
end

function normalize_expectation_weights(p)                                                                 # Esta función normaliza los pesos usados para calcular el precio esperado de largo plazo       
    w_total = p.w_energy_price_to_expect + p.w_expected_cost_to_expect + p.w_fuel_price_to_expect       
    return w_total <= 0.0 ? (0.70, 0.20, 0.10) :                                                        
           (p.w_energy_price_to_expect / w_total,           
            p.w_expected_cost_to_expect / w_total,
            p.w_fuel_price_to_expect / w_total)
end


##############################################################################
# 3. ESTADO INICIAL - Desde qué valores parte el modelo en el mes 0
##############################################################################

function initial_state_liberalized(p)               
    # Recibe los parámetros de un escenario.
    # Calcula el vector inicial de variables de estado.
    # Devuelve u0, que es el estado inicial que usará el solver para resolver las respectivas ecuaciones diferenciales.

    N  = Int(p.N_delay)             # Orden del retardo para rebote directo
    u0 = zeros(N + 10)              # Estado inicial

    # Como N = 3, 3 posiciones para la cadena de retardo del rebote directo
    # 10 posiciones para los otros estados dinámicos del modelo.
    # entonces u0 = 13 estados que van cambiando en el tiempo

    F_EC_base = p.S_base / p.epsilon_0                  # Consumo energetico base
            # = 100      /   1.0
            # = 100 p.u.
    GP0       = F_EC_base * (1.0 + p.reserve_ratio)     # Generation Park inicial 
            # = 100 p.u. * (1 + 0.15)
            # = 115 p.u.

    retirements0 = GP0 / p.asset_life_months                                 # cuánta capacidad se retiraría inicialmente por envejecimiento de activos.
                # = 115 / 300 = 0.3833 p.u./mes 

    mw_to_model_unit = GP0 / max(p.GP_base_physical_MW, 1e-6)                # Conversion MW a escala interna del modelo.
                   # = 115 p.u. / 39200 MW = 0.00293 p.u./MW 

    committed_uc0 = p.initial_construction_pipeline_MW * mw_to_model_unit    # Capacidad inicial en construccion 
                # = 8500 MW * 0.00293 p.u./MW = 24.9 p.u.

    committed_id0 = p.initial_approval_pipeline_MW * mw_to_model_unit        # Capacidad inicial aprobada  
                # = 5000 MW * 0.00293 p.u./MW = 14.65 p.u.

    UC0 = retirements0 * p.tau_construction + committed_uc0                  # Under Construction inicial 
                # = 0.3833 p.u./mes * 12 meses + 24.9 p.u. = 29.5 p.u.

    ID0 = retirements0 * p.tau_ID_up + committed_id0                         # Investment Decision stock inicial en unidad modelo
      # = 0.3833 p.u./mes * 7 meses + 14.65 p.u. = 17.35 p.u.
      
# Asignación de valores iniciales al vector de estado u0 #

    u0[N + 1]  = 0.0          # CMS
    u0[N + 2]  = GP0          # Generation Park
    u0[N + 3]  = p.S_base     # Q suavizado de demanda util para reserva
    u0[N + 4]  = p.P_E_base   # EnergyPrice
    u0[N + 5]  = p.P_E_base   # Pexp
    u0[N + 6]  = ID0          # Investment Decision stock
    u0[N + 7]  = UC0          # Under Construction stock
    u0[N + 8]  = 0.0          # IE state
    u0[N + 9]  = 0.0          # Budget Constraint
    u0[N + 10] = 0.0          # Indirect Rebound

    return u0
end

#######################################################################################################################################################################################################
# 4. ALGEBRAICAS
#######################################################################################################################################################################################################

#############################################################################################
# 4.1 PRIMERA SECCION DE ALGEBRAICAS
# Esta seccion lee los estados actuales del modelo desde el vector u.
# Luego calcula las referencias base del sistema, como demanda util,
# consumo energetico, parque generador, gasto energetico y costo unitario inicial.
# Tambien define ingreso, demanda base, impuesto efectivo y margen de reserva.
# Finalmente calcula costos, precio spot, precio de energia, precio al consumidor
# y precio esperado de largo plazo.
#############################################################################################


 function liberalized_algebraics(u, p, t)    # toma el estado actual del modelo en el mes t y prepara las variables base para calcular demanda, precios, reserva, consumo, gasto e índices. 

    # "u" Es el vector de estados actuales del modelo. Aquí vienen las variables que cambian dinámicamente    
    # "p" Son los parámetros del escenario que se están simulando. Sirven para calcular las ecuaciones del modelo según las características de cada escenario.
    # "t" Es el mes actual de la simulación, que se usa para activar políticas y calcular variables dependientes del tiempo. 

    N = Int(p.N_delay)          # número de etapas del retardo del rebote directo.
                                # Como N = 3, las primeras tres posiciones del vector u se usan para representar el retardo del rebote directo.


    L                = u[1:N]                                                # Variables de retardo dinámico para el rebote directo. 
    CMS              = u[N + 1]                                              # Ahorro monetario acumulado respecto al presupuesto energético base.     
    GP               = max(0.0, u[N + 2])                                    # Capacidad instalada, se asegura que no sea negativa. Se actualiza por inversión y retiros.   
    Q                = max(1e-6, u[N + 3])                                   # Demanda util suavizada para planificación y cálculo de reserva.
    EnergyPrice      = max(0.0, u[N + 4])                                    # Precio energético actual del modelo. Es el precio base de la energía consumida expresado en p.u.  
    Pexp             = max(0.0, u[N + 5])                                    # Precio esperado de largo plazo
    ID               = max(0.0, u[N + 6])                                    # ID representa proyectos de generación aprobados, pero que aún no pasan a construcción.
    UC               = max(0.0, u[N + 7])                                    # Capacidad que está en construcción (No aumenta el GP hasta que se completa la construcción)
    IE               = clamp(u[N + 8], 0.0, p.IE_max)                        # Representa la mejora de eficiencia energética.
    BudgetConstraint = clamp(u[N + 9], 0.0, p.budget_constraint_max)         # Representa cuánto la presión financiera limita la demanda útil.
    IndirectRebound  = clamp(u[N + 10], 0.0, p.indirect_rebound_max)         # Rebote indirecto


    S0     = p.S_base                           # Referencia inicial de demanda de trabajo útil - Sirve para medir crecimiento de demanda y como base para calcular elasticidades y efectos de precio.
        #  = 100 p.u.

    FEC0   = p.S_base / p.epsilon_0             # Consumo energético final base - Sirve para medir el crecimiento del consumo energético y para calcular el gasto energético bruto base.
        #  = 100 p.u. / 1.0 = 100 p.u.

    GP0    = FEC0 * (1.0 + p.reserve_ratio)     # GP inicial - Sirve como referencia para medir crecimiento del parque generador y seguridad de suministro.
        #  = 100 p.u. * (1 + 0.15) = 115 p.u.

    EEXP0  = p.P_E_base * FEC0                  # Gasto energético bruto base  - Sirve para comparar si el gasto sube o baja durante la simulación.
        #  = 1.0 * 100 p.u. = 100 p.u.

    UCUW0  = p.P_E_base / p.epsilon_0           # Costo unitario inicial del trabajo útil - Sirve para para medir equidad y rebote directo.
        #  = 1.0 / 1.0 = 1.0

    epsilon_t = p.epsilon_0 + IE                # Eficiencia total actual
            # = 1 + IE


    # Parte socioeconómica - Esta parte fija el ingreso del agente vulnerable representativo
       # para las simulaciones se definió el ingreso fijo y que no aumentase con el tiempo

    income_t        = p.Income_base                     # Ingreso normalizado (se usa en el modelo y plots)
    income_CLP_t    = p.Income_base_CLP                 # Ingreso en CLP (se usa para cálculos en el diagnostico)
    energy_budget_available_t = p.EnergyBudget_base     # Presupuesto energético disponible de referencia.



    # Crecimiento exogeno suave de demanda util base.
        # Permite que la demanda útil crezca levemente por evolución estructural del sistema, a pesar de que no haya rebote.
        
    demand_growth_factor = exp(p.demand_growth_rate * t)            # factor de crecimiento de la demanda en el tiempo.
    S_base_t             = p.S_base * demand_growth_factor          # demanda útil base ajustada por ese crecimiento.


    # Activación del impuesto con rampa temporal.
    tax_activation_t = policy_ramp(t, p.t_T, p.tau_tax_implementation)      # Toma valores entre 0 (no activado) y crece gradualmente hasta 1 (activación completa)
    tax_val_t        = p.P_T * tax_activation_t                             # Impuesto efectivo aplicado al precio de la energía.    

    # Margen de reserva.
    reserve_floor_t = p.reserve_demand_floor * p.S_base                                                              # Piso mínimo de demanda para el cálculo de reserva, sirve para que el margen de reserva no caiga tanto
    smoothed_capacity_requirement_t = p.S_base * (max(Q, 1e-6) / max(p.S_base, 1e-6))^p.peak_demand_elasticity       # Se calcula la demanda relevante para capacidad.    
                    # Usa Q (demanda suavizada) y la transforma en una necesidad de capacidad.
                    # peak_demand_elasticity hace que la necesidad de capacidad responda de forma no lineal.
                    # Si Q aumenta, la capacidad requerida aumenta más que proporcionalmente si la elasticidad es mayor que 1.

    # Esta condición asegura que la demanda usada para calcular reserva no sea excesivamente baja. Si la demanda suavizada cae bajo el piso, el modelo usa una corrección parcial               
    if smoothed_capacity_requirement_t >= reserve_floor_t                                                                                     
        reserve_demand_for_margin_t = smoothed_capacity_requirement_t                                                                       
    else
        reserve_demand_for_margin_t = reserve_floor_t + p.reserve_floor_blend * (smoothed_capacity_requirement_t - reserve_floor_t)
    end

    # Esta línea evita que la demanda para calcular reserva sea cero o negativa.
    reserve_demand_for_margin_t = max(1e-6, reserve_demand_for_margin_t)

    # Aquí se calcula la capacidad efectiva usada para el margen de reserva.
    if GP >= GP0                                                                                # Esta primera restriccion significa que la capacidad nueva sobre el nivel inicial no se reconoce completamente como capacidad firme.
        GP_effective_for_margin_t = GP0 + p.marginal_capacity_credit * (GP - GP0)               # Sólo se reconoce una fracción definida por marginal_capacity_credit    
    else                                                                                        # Si GP es menor que GP0, se usa directamente:
        GP_effective_for_margin_t = GP
    end

    # Aquí se calcula el Margen de Reserva - mide cuánta capacidad efectiva sobra respecto a la demanda relevante del sistema.
    reserve_margin_t = (GP_effective_for_margin_t - reserve_demand_for_margin_t) / reserve_demand_for_margin_t

    # Costos y Precios.

    # Primero se calcula el costo de producción de energía
    production_cost_t = p.variable_cost_base * p.expected_fuel_price_index *                # Depende de Costo variable base de generación, el precio esperado del combustible
                        (1.0 + p.k_pexp_to_prodcost * (Pexp / p.P_E_base - 1.0))            # y de la sensibilidad del costo de producción a las expectativas de precio.

    production_cost_t = max(0.85, production_cost_t)                                        # impone un piso mínimo al costo de producción

    # Luego se calcula el precio esperado de combustibles.
    expected_fuel_price_t = p.P_E_base * p.expected_fuel_price_index
                        # = 1.0 * 1.04 = 1.04 p.u. = 4% sobre el precio base de energía

    # Aquí se calcula el costo esperado total.                    
    expected_cost_t = max(0.05, production_cost_t + p.fixed_cost_adder + p.equil_reserve_add) 
          # Se suma =  costo de produccion + costo fijo adicional + ajuste para que el costo refleje el efecto de la reserva.
    
    # Luego se calcula precio spot, señal de escasez cuando falta reserva y descuento por exceso de reserva     
    spot_price_t, scarcity_t, surplus_discount_t = spot_price_liberalized(production_cost_t, reserve_margin_t, p)               # Si hay escasez, el precio spot sube

    # Se normalizan los pesos para combinar precio spot y precio esperado.
    w_spot, w_expected = normalize_price_weights(p)    
    market_price_target_raw_t = max(0.05, w_spot * spot_price_t + w_expected * Pexp)   # Se calcula el precio objetivo bruto de mercado.

    # Aquí se calcula el precio objetivo final de energía.
    energy_price_target_t = (1.0 - p.w_anchor_to_P_E_base) * market_price_target_raw_t +            # w_anchor_to_P_E_base evita que el precio de energía se mueva demasiado rápido o demasiado lejos de la base
                             p.w_anchor_to_P_E_base * p.P_E_base
    energy_price_target_t = max(0.05, energy_price_target_t)

    # Se calcula el precio que enfrenta el consumidor luego de agregar el impuesto.
    consumer_price_t        = EnergyPrice * (1.0 + tax_val_t)

    # Se calcula el precio que afecta la respuesta de demanda.
        # No necesariamente es igual al precio consumidor. No usa todo el impuesto, porque se asume que la respuesta del consumo al impuesto es parcial
    demand_response_price_t = EnergyPrice * (1.0 + p.tax_demand_pass_through * tax_val_t)

    # Aquí se normalizan los pesos para calcular el precio esperado de largo plazo.
    # Los pesos corresponden a precio actual de la energia, costo esperado y combustible esperado
    w_energy, w_cost, w_fuel = normalize_expectation_weights(p)

    # Se calcula un recargo adicional de escasez para el precio esperado. 
        # Si hay escasez actual, los agentes esperan precios futuros más altos.
    scarcity_expectation_premium_t = p.expectation_scarcity_gain * scarcity_t

    # se calcula un recargo adicional por crecimiento de demanda.
    demand_expectation_premium_t = p.expectation_demand_gain * max(0.0, Q / max(S0, 1e-6) - 1.0)    
        # Si Q está por encima de S0, significa que la demanda suavizada creció respecto a la base.
        # Eso aumenta el precio esperado.

    # Se calcula el objetivo del precio esperado de largo plazo (Pexp)    
    # precio esperado objetivo = mezcla ponderada de precio actual, costo esperado y combustible esperado, más primas por escasez y demanda.
    pexp_target_t = max(0.05,
        w_energy * EnergyPrice +
        w_cost   * expected_cost_t +
        w_fuel   * expected_fuel_price_t +
        scarcity_expectation_premium_t +
        demand_expectation_premium_t
    )


#########################################################################################
# 4.2 SEGUNDA SECCION DE ALGEBRAICAS
# Esta seccion calcula el rebote directo a partir de la reduccion del costo util.
# Luego aplica el efecto precio y la restriccion presupuestaria sobre la demanda util.
# Con la demanda util y la eficiencia calcula el consumo energetico.
# Despues calcula gasto bruto, compensacion focalizada y gasto neto.
# A partir del gasto neto calcula presion presupuestaria, CMS normalizado,
# holgura financiera, estres financiero, restriccion presupuestaria objetivo
# y rebote indirecto objetivo.
# Finalmente inicia el bloque de inversion y parque generador, calculando retiros,
# completaciones, rentabilidad esperada y multiplicador de inversion.
#############################################################################################

    # Rebote directo: la señal usa precio efectivo para demanda, no compensacion.
    U_CUW_signal_t = demand_response_price_t / max(epsilon_t, 1e-6)                                      # Aquí se calcula el costo unitario percibido del trabajo útil.
    useful_cost_reduction_t = max(0.0, 1.0 - U_CUW_signal_t / max(UCUW0, 1e-6))                          # se calcula cuánto bajó el costo útil respecto al costo inicial. Esta variable mide cuánto se abarató el uso de energía útil respecto al inicio
    rebound_activation_t = policy_ramp(t, p.t_DRE, p.tau_DRE_activation)                                 # Activación temporal del rebote
    X_DRE_t = clamp(p.eta_epsilon_S * useful_cost_reduction_t * rebound_activation_t, 0.0, p.DRE_max)    # Se calcula el objetivo del rebote directo, limitado a un DRE_max
            # = ( Sensibilidad de la demanda útil frente a la mejora de eficiencia + 
            #      Cuánto bajó el costo útil + Qué tan activado está el rebote en ese mes)
    DRE_t   = clamp(L[N], 0.0, p.DRE_max)                                                                # Se obtiene el rebote directo efectivo despues del retardo

    # Efecto precio sobre demanda.
    price_dev_t    = (demand_response_price_t / p.P_E_base) - 1.0           # Aquí se calcula cuánto cambió el precio que afecta la demanda respecto al precio base.                         
    price_effect_t = max(0.0, 1.0 + p.eta_P * price_dev_t)                  # Se calcula el efecto del precio sobre la demanda. Como eta_P es negativo, si el precio sube, la demanda baja.

    # Restriccion presupuestaria.
        # Esta línea convierte la restricción presupuestaria en un multiplicador que reduce la demanda.
        # Si BudgetConstraint = 0, entonces budget_constraint_effect_t = 1, y no se reduce la demanda.
        # Si BudgetConstraint aumenta, el denominador aumenta y el efecto baja de 1.
    budget_constraint_effect_t = 1.0 / (1.0 + p.beta_budget * max(0.0, BudgetConstraint))

    # Demanda util efectiva.
        # Aquí se calcula la demanda útil antes de aplicar restricción presupuestaria.
        # La eficiencia puede aumentar demanda por rebote, pero el precio puede reducirla.
         S_before_budget_t = S_base_t * (1.0 + DRE_t) * (1.0 + IndirectRebound) * price_effect_t         
                         # = demanda base , rebote directo, rebote indirecto, efecto precio 

         S_actual_t        = max(0.0, S_before_budget_t * budget_constraint_effect_t)      # Aquí se calcula la demanda útil final despues de considerar la restricción presupuestaria.              

    # Consumo energetico FEC.
    F_EC_t  = max(0.0, S_actual_t / max(epsilon_t, 1e-6))

    # Gasto energetico bruto.
    E_Exp_GROSS_t = consumer_price_t * F_EC_t       # precio que enfrenta el consumidor por consumo energético, sin considerar la compensación


    # COMPENSACION FOCALIZADA POR RECICLAJE PARCIAL DEL IMPUESTO - No reduce el precio marginal ni la señal de demanda; solo reduce gasto neto #

     # Aquí se estima una recaudación mensual aproximada del impuesto. Es un proxy interno para calcular la compensación.
        tax_revenue_proxy_t = max(0.0, tax_val_t * EnergyPrice * F_EC_t)                
                          # = (impuesto efectivo, precio de la energia, consumo energetico)      

     # Aquí se mide qué parte del ingreso se va a gasto energético antes de compensación
        gross_budget_share_t = E_Exp_GROSS_t / max(income_t, 1e-6)      

     # Se calcula si el gasto energético supera el umbral de asequibilidad                                   
        affordability_excess_t = max(0.0, gross_budget_share_t / max(p.energy_budget_share_target, 1e-6) - 1.0)
            # El umbral del modelo está en 10%
            # Si el gasto está bajo el 10%, esta variable vale 0
            # Si está sobre el 10%, aparece exceso de asequibilidad.

     # Aquí se suaviza el gatillo de compensación.       
        # Esta línea transforma el problema de asequibilidad en una señal suave de activación de compensación.
        affordability_trigger_raw_t = tanh(affordability_excess_t / max(p.compensation_trigger_width, 1e-6))
     
     # Se calcula el gatillo final de compensación.   
        # define qué tan fuerte se activa la compensación según el problema de asequibilidad
        # compensation_trigger_floor da un piso mínimo de activación.
        affordability_trigger_t = p.compensation_trigger_floor + (1.0 - p.compensation_trigger_floor) * affordability_trigger_raw_t
            
     # Aquí se activa temporalmente la compensación desde el mes 18.
        compensation_activation_t = policy_ramp(t, p.t_compensation, p.tau_compensation_implementation)

     # Aquí se define cuándo empieza el retiro gradual de la compensación.
        # La compensación tiene un periodo completo de aplicación y luego empieza a reducirse   
        phaseout_start_t = p.t_compensation + p.compensation_full_duration
        # Aquí se calcula la rampa de retiro de la compensación.
        phaseout_t = policy_ramp(t, phaseout_start_t, p.compensation_phaseout_tau)
     
     # Aquí se calcula el perfil de compensación en el tiempo. 
        # Esta línea reduce gradualmente la compensación después de su periodo principal.    
        compensation_profile_t = max(0.0, 1.0 - p.compensation_phaseout_fraction * phaseout_t)

     # Aquí se calcula la compensación bruta, es decir,
     # la compensación antes del tope, reciclando una parte del impuesto hacia el agente vulnerable representativo
        compensation_raw_t = p.compensation_enabled *               # define la compensación está habilitada
                             p.compensation_recycling_share *       # porcentaje de reciclaje del impuesto
                             tax_revenue_proxy_t *                  # recaudación aproximada
                             affordability_trigger_t *              # gatillo de asequibilidad
                             compensation_activation_t *            # activación temporal
                             compensation_profile_t                 # perfil de retiro

     # Aquí se calcula el tope máximo de compensación - No puede superar cierto porcentaje del gasto bruto           
        # Esta línea evita que la compensación sea excesiva respecto al gasto energético        
        compensation_cap_t = p.compensation_cap_share * E_Exp_GROSS_t
     
     # Aquí se aplica el tope - define la compensación efectiva, respetando el límite máximo   
        compensation_t = min(compensation_raw_t, compensation_cap_t)


    # CALCULO GASTO ENERGETICO NETO - Gasto final que enfrenta el agente vulnerable después de la compensación
        E_Exp_NET_t = max(0.0, E_Exp_GROSS_t - compensation_t)      

     # Se define que el gasto energético principal será el gasto neto.
        E_Exp_t = E_Exp_NET_t                


    # Se calcula el costo unitario BRUTO del trabajo útil
        U_CUW_gross_t = S_actual_t > 1e-9 ? E_Exp_GROSS_t / S_actual_t : NaN
                    # = gasto bruto dividido por demanda útil
                    # Mide cuánto cuesta cada unidad de trabajo útil antes de compensación
                    # Si la demanda útil es prácticamente cero, devuelve NaN para evitar error.

    # Se calcula el costo unitario neto del trabajo útil - Importante para el índice de equidad.
        U_CUW_net_t   = S_actual_t > 1e-9 ? E_Exp_NET_t   / S_actual_t : NaN
                    # = gasto neto dividido por demanda útil+
                    # Mide cuánto cuesta cada unidad de trabajo útil después de compensación.

    # Variables socioeconomicas internas.

        # Aquí se calcula el porcentaje del ingreso destinado a energía -  Usa gasto neto, no gasto bruto.
             energy_budget_share_t = E_Exp_NET_t / max(income_t, 1e-6)   

        # Aquí se calcula la presión presupuestaria.     
             budget_pressure_raw_t = max(0.0, energy_budget_share_t / max(p.energy_budget_share_target, 1e-6) - 1.0)
                # Si el gasto energético neto está bajo el umbral de 10%, vale 0.
                # Si lo supera, aparece presión presupuestaria.
                # Mide cuánto se supera el umbral de asequibilidad energética

    # Aquí se normaliza el CMS y se convierte en una señal proporcional comparable            
        CMS_norm_t = CMS / max(energy_budget_available_t * p.savings_ref_months, 1e-6)
             # = CMS / presupuesto energético de referencia * meses de referencia(12)
    
    # Aquí se calcula el ahorro o pérdida mensual respecto al presupuesto energético disponible.         
        monthly_savings_signal_t = (energy_budget_available_t - E_Exp_NET_t) / max(energy_budget_available_t, 1e-6)
                    # Si el gasto neto es menor al presupuesto, la señal es positiva.
                    # Si el gasto neto es mayor, la señal es negativa.

    # Aquí se calcula la holgura financiera.
        financial_slack_t  = max(0.0, CMS_norm_t) + p.monthly_slack_weight * max(0.0, monthly_savings_signal_t)
                     # = CMS positivo acumulado, ahorro mensual positivo
                     # Si hay ahorro acumulado y ahorro mensual, aumenta la holgura.

    # Aquí se calcula el estrés financiero.
        financial_stress_t = max(0.0, -CMS_norm_t) + p.monthly_stress_weight * max(0.0, -monthly_savings_signal_t)
                     # = CMS negativo acumulado, pérdida mensual.
                     # Si el gasto supera el presupuesto, aumenta el estrés.

    # Aquí se calcula la presión total que alimenta la restricción presupuestaria.
        budget_constraint_raw_t    = budget_pressure_raw_t + p.beta_stress * financial_stress_t
                             # = presión presupuestaria directa, estrés financiero ponderado.
                             
    # Aquí se calcula el objetivo de la restricción presupuestaria.                         
        budget_constraint_target_t = p.budget_constraint_max * tanh(budget_constraint_raw_t / max(p.budget_pressure_scale, 1e-6))
                    # tanh suaviza la respuesta y evita que crezca infinitamente.
                    # budget_constraint_max fija el máximo posible.
                    # Se calcula hacia dónde debería moverse la restricción presupuestaria según la presión financiera
                    
    # Aquí se calcula el objetivo del rebote indirecto - Se convierte la holgura financiera en rebote indirecto potencial.
        indirect_rebound_target_t = p.indirect_rebound_max * tanh(p.k_indirect * financial_slack_t)
                    # Si hay holgura financiera, aumenta el rebote indirecto.
                    # tanh hace que el rebote se acerque a un máximo, pero no lo supere.

    # Inversion y parque generador.
    # GP no contiene crecimiento exogeno. La capacidad cambia solo por completacion menos retiros.
    # Inversion por logica de Olsina et al. : I = m(PI) * I_ref, con decision y construccion retrasadas.

        # Aquí se calculan los retiros de capacidad.
            retirements_t = GP / p.asset_life_months
                # Mientras mayor sea GP, mayor será el retiro mensual promedio.
        
        # Aquí se calcula cuánta capacidad en construcción se completa y entra a GP en cada mes.        
            completion_t  = UC / p.tau_construction
                # Si hay más capacidad en construcción, mayor será la completación.

        # Aquí se calcula el precio esperado neto para inversión.
            net_expected_price_for_investment_t = max(0.05, Pexp / (1.0 + p.tax_investment_wedge * tax_val_t))
                # Parte desde Pexp, pero descuenta una brecha asociada al impuesto.
                # La idea es que el impuesto puede afectar la rentabilidad percibida para invertir.

        # Aquí se calcula el costo esperado que se compara contra el precio esperado para evaluar rentabilidad      
            expected_cost_for_investment_t = max(0.05, expected_cost_t - 0.090)
                # Se toma expected_cost_t y se ajusta levemente hacia abajo.
                # El max(0.05, ...) evita costos negativos o demasiado bajos.

        # Aquí se calcula el margen de rentabilidad - mide qué tan rentable sería invertir en nueva capacidad. 
            profitability_margin_t = (net_expected_price_for_investment_t - expected_cost_for_investment_t) /max(expected_cost_for_investment_t, 1e-6)
                                 # = (precio esperado neto - costo esperado) / dividido por costo esperado
                                 # Si el precio esperado supera el costo, el margen es positivo.
                                 # Si el costo supera el precio, el margen cae.

        # Aquí se mide la intensidad relativa del impuesto en ese momento - calcula qué tan activado está el impuesto respecto a su valor máximo.
            tax_policy_strength_t = p.P_T > 0.0 ? tax_val_t / max(p.P_T, 1e-6) : 0.0
                # Si no hay impuesto, queda en cero.
                # Si el impuesto está completamente activado, se acerca a 1.
        
        # Aquí se calcula una adaptación gradual del mercado despues del impuesto.        
            tax_market_adaptation_t = p.tax_market_adaptation_gain * tax_policy_strength_t * (1.0 - exp(-max(0.0, t - p.t_T) / 24.0))
                # el mercado puede ajustar parcialmente sus expectativas o decisiones frente al nuevo contexto tributario.

        # Aquí se calcula el índice de rentabilidad esperado.
            PI_raw_t = 1.0 + (profitability_margin_t - p.RRR) / 0.30 + tax_market_adaptation_t
                   # RRR es el retorno requerido (0,07, Olsina usaba 0,125 para los 24 años).
                   # Si el margen de rentabilidad supera el retorno requerido, el PI sube.

        # Aquí se limita el PI entre 0.05 y 2.50 para evitar valores extremos     
            PI_raw_t = clamp(PI_raw_t, 0.05, 2.50)

        # Aquí se convierte el PI en un multiplicador de inversión.    
            m_t      = investment_multiplier(PI_raw_t, p)
                            # Si PI es alto, el multiplicador aumenta.
                            # Si PI es bajo, el multiplicador disminuye.


##########################################################################################################################
# 4.3 TERCERA SECCION DE ALGEBRAICAS
# Esta seccion completa el bloque de inversion y parque generador.
# Calcula la demanda usada para planificacion, la capacidad objetivo, la brecha de capacidad y la expansion de referencia.
# Luego aplica frenos y condiciones de activacion para obtener la expansion endogena, aprobaciones, inicios de construccion y cancelaciones.
# Despues calcula las señales que ajustan la eficiencia dinamica.
# Finalmente calcula los indices del trilema y normaliza variables para diagnosticos y graficos.
##########################################################################################################################

  # Aquí se calcula la demanda que el sistema usa para planificar capacidad futura. 
    planning_demand_t = p.demand_planning_weight * S_base_t +
                        (1.0 - p.demand_planning_weight) * reserve_demand_for_margin_t
                    # = demand_planning_weight define cuánto pesa cada una , Demanda útil base del sistema, Demanda usada para calcular reserva   
                    # calcula una demanda de planificación combinando la demanda base con la demanda relevante para reserva

  # Aquí se proyecta la demanda de planificación hacia adelante.                  
    forecast_Q_t = planning_demand_t * exp(p.demand_growth_rate * p.investment_forecast_horizon)
                # = Crecimiento de demanda, Horizonte de previsión para inversión.
                # La idea es que los inversionistas no miran solo la demanda actual, sino la demanda esperada en el futuro.

  # Aquí se calcula la capacidad objetivo por planificación.
    planning_capacity_target_t = forecast_Q_t * (1.0 + p.planning_reserve_ratio)
        # demanda futura esperada más un margen de reserva de planificación
        # calcula cuánta capacidad se necesitaría para cubrir la demanda futura con reserva

  # Aquí se calcula la capacidad necesaria para cumplir la reserva actual requerida.      
    reserve_target_capacity_t = reserve_demand_for_margin_t * (1.0 + p.reserve_ratio)

  # Aquí se define la capacidad objetivo final.  
    capacity_target_t = max(planning_capacity_target_t, reserve_target_capacity_t)
            # Toma el valor mayor entre:
            # capacidad requerida por demanda futura, capacidad requerida por reserva actual.
            # El modelo elige la exigencia más alta entre planificación futura y reserva actual.

  # Aquí se asume que la tasa de inversión de reemplazo debe igualar los retiros. 
    replacement_start_rate_t = retirements_t
            # Si se retira capacidad, se necesita iniciar capacidad nueva para reemplazarla.

  # Aquí se calcula cuánto stock de decisiones de inversión debería existir solo para reemplazar retiros.          
    replacement_id_ref_t     = replacement_start_rate_t * p.tau_ID_up
            # tau_ID_up representa el tiempo asociado a la etapa de decisión o aprobación.
            # calcula el nivel de decisiones de inversión necesario para sostener el reemplazo
    
  # Aquí se calcula cuánto stock en construcción debería existir solo para reemplazar retiros.                                                       
    replacement_uc_ref_t     = replacement_start_rate_t * p.tau_construction


  # EXPANSION COMPROMETIDA
     # Aquí se separa qué parte de ID corresponde a expansión real y no solo reemplazo.
        expansion_id_stock_t  = max(0.0, ID - replacement_id_ref_t)
            # Si ID es mayor que el stock necesario para reemplazo, el excedente se considera expansión.
            # identifica cuánta inversión aprobada corresponde realmente a expansión adicional

     # Aqui se separa construcción de reemplazo y construcción de expansión.    
        expansion_uc_stock_t  = max(0.0, UC - replacement_uc_ref_t)
            # identifica cuánta capacidad en construcción corresponde a expansión adicional.

     # Aquí se calcula la expansión comprometida nominal - estima cuánta expansión futura ya está comprometida en el pipeline.      
        committed_expansion_nominal_t = p.pipeline_credit * (expansion_id_stock_t + expansion_uc_stock_t)
            # Suma la expansión aprobada y la expansión en construcción, pero aplica pipeline_credit.
            # pipeline_credit representa que no todo el pipeline se considera completamente seguro o comprometido.

     # Aquí se convierte la expansión nominal en expansión efectiva.      
        committed_expansion_effective_t = p.marginal_capacity_credit * committed_expansion_nominal_t
            # Aplica marginal_capacity_credit, porque no toda nueva capacidad aporta como capacidad firme completa.
            # Reconoce solo una parte de la expansión comprometida como capacidad efectiva para seguridad.


   # BRECHA DE CAPACIDAD
     # Aquí se calcula la brecha futura de capacidad.         
       future_capacity_gap_effective_t = max(0.0, capacity_target_t - GP_effective_for_margin_t - committed_expansion_effective_t)
                                     # = capacidad objetivo, GP efectiva actual, expansion ya comprometida
                                     # Si todavía falta capacidad, aparece una brecha positiva.
                                     # mide cuánta capacidad efectiva faltaría en el futuro, descontando lo que ya está comprometido

     # Aquí se calcula una brecha inmediata de reserva.                                
       immediate_reserve_gap_effective_t = max(0.0, reserve_target_capacity_t - GP_effective_for_margin_t)
            # No mira la demanda futura. 
            # Mide si actualmente falta capacidad para cumplir el margen de reserva

     # Aquí se suma la brecha futura y la brecha inmediata - junta la falta de capacidad futura y la falta de reserva actual.
       capacity_gap_effective_t = future_capacity_gap_effective_t + immediate_reserve_gap_effective_t

     # Aquí se convierte la brecha efectiva en brecha nominal.   
       capacity_gap_nominal_t = capacity_gap_effective_t / max(p.marginal_capacity_credit, 1e-6)
            # no toda nueva capacidad cuenta como capacidad efectiva
            # se necesita más capacidad nominal para cubrir una brecha efectiva.
            # Entonces, esta linea calcula cuánta capacidad nominal habría que construir para cubrir la brecha efectiva.

   # REFERENCIAS DE EXPANSIÓN 
     # Aquí se calcula la necesidad adicional de capacidad causada por crecimiento futuro de demanda.         
        forecast_growth_effective_need_t = max(0.0, planning_capacity_target_t - reserve_target_capacity_t)
                # Si la capacidad de planificación futura supera la capacidad requerida por reserva actual, aparece una necesidad por crecimiento.

     # Aquí se calcula una tasa de inversión de referencia asociada al crecimiento esperado de demanda.   
        # transforma la necesidad futura de capacidad por crecimiento de demanda en una tasa mensual de inversión.”        
        demand_growth_reference_t = p.demand_growth_capacity_gain * forecast_growth_effective_need_t /
                                    max(p.investment_forecast_horizon, 1e-6) /                                  
                                    max(p.marginal_capacity_credit, 1e-6)
                            # Divide por el horizonte de planificación para convertir una necesidad futura en una tasa mensual
                            # También corrige por marginal_capacity_credit.

     # Aquí se calcula una tasa de inversión necesaria para cerrar la brecha de capacidad.
        gap_closure_reference_t = p.capacity_gap_response_gain * capacity_gap_nominal_t / max(p.tau_gap, 1e-6)
                # Si falta mucha capacidad, esta referencia aumenta.
                # tau_gap indica en cuánto tiempo se intenta cerrar la brecha.

     # Aquí se calcula una referencia de inversión inducida por escasez - agrega inversión inducida por escasez de capacidad
        scarcity_reference_t = p.scarcity_capacity_gain * scarcity_t * S0 /
                           max(p.tau_gap * p.marginal_capacity_credit, 1e-6)
                # Si la escasez aumenta, el sistema genera una señal adicional para invertir.
                        
     # Aquí se suma toda la expansión de referencia.           
        expansion_ref_t = demand_growth_reference_t + gap_closure_reference_t + scarcity_reference_t
                # Incluye tres motivos para invertir:
                    # crecimiento esperado de demanda, cierre de brecha de capacidad, escasez.
                # junta las tres razones principales para expandir capacidad.
                
   # FRENOS Y PUERTAS DE INVERSION             
     # Aquí se calcula un freno por sobre reserva - Si el sistema tiene demasiada capacidad sobrante, la inversión se frena.
        brake_t = over_reserve_brake(reserve_margin_t, p)
            # Evita que el modelo siga invirtiendo cuando ya existe mucha reserva

     # Aquí se calcula un freno asociado al impuesto.   
        tax_expansion_brake_t = 1.0 / (1.0 + p.tax_expansion_brake_k * tax_val_t)
            # Si el impuesto está activo, puede reducir la expansión de capacidad, al reducir la señal económica de inversión.

     # Aquí se calcula una puerta de inversión por rentabilidad - activa la inversión cuando la rentabilidad esperada supera cierto umbral 
        profitability_gate_t = clamp((PI_raw_t - p.investment_pi_threshold) / max(p.investment_pi_width, 1e-6), 0.0, 1.0)
            # Si PI_raw_t está bajo el umbral, esta puerta se acerca a 0.
            # Si PI_raw_t supera el umbral, se acerca a 1.
      
     # Aquí se calcula una puerta por necesidad de reserva - activa inversión cuando la reserva está por debajo del nivel requerido
        reserve_need_gate_t  = clamp((p.reserve_ratio - reserve_margin_t) / 0.06, 0.0, 1.0)
            # Si el margen de reserva está bajo el requerido, esta puerta se activa.
            # Si hay suficiente reserva, queda cerca de 0.

     # Aquí se calcula una puerta por presión de demanda - activa inversión cuando la demanda planificada crece sobre la base    
        demand_pressure_gate_t = clamp(p.demand_investment_gate_gain * max(0.0, planning_demand_t / max(S0, 1e-6) - 1.0), 0.0, 1.0)
            # Si la demanda de planificación supera la demanda base, esta puerta se activa.

     # Aquí se calcula la puerta total de expansión - se decide si existen condiciones para activar expansión de capacidad.
        expansion_gate_t = clamp(max(profitability_gate_t, reserve_need_gate_t, demand_pressure_gate_t, capacity_gap_effective_t > 1e-6 ? 1.0 : 0.0), 0.0, 1.0)
            # El modelo permite expansión si se cumple cualquiera de estas condiciones:
                # hay rentabilidad suficiente, falta reserva, hay presión de demanda, existe brecha efectiva de capacidad.

  # EXPANSION ENDÓGENA, APROBACIÓN Y CONSTRUCCIÓN 
    # Aquí se calcula la tasa de expansión endógena.
        endogenous_expansion_rate_t = tax_expansion_brake_t * brake_t * expansion_gate_t * m_t * expansion_ref_t
                                  # = freno por impuesto, freno por sobre reserva, puerta de expansión, multiplicador de inversión, expansión de referencia.
            # calcula cuánta nueva expansión se activa realmente según señales de mercado, rentabilidad y capacidad
    
    # Aquí se calcula la tasa deseada de inicio de proyectos.       
        desired_start_rate_t = replacement_start_rate_t + endogenous_expansion_rate_t
                           # = reemplazo por retiros, expansion nueva
            # suma la inversión necesaria para reemplazar retiros y la inversión para expandir capacidad.               

    # Aquí se calcula el stock deseado de decisiones de inversión.         
        desired_ID_t = desired_start_rate_t * p.tau_ID_up
            # Si se quiere iniciar más proyectos, se necesita un mayor stock de decisiones aprobadas.
            # Esto calcula cuántas decisiones de inversión deberían existir para sostener la tasa deseada de proyectos.

    # Aquí se define la tasa de aprobación de nuevos proyectos.
        approval_rate_t = desired_start_rate_t     
            #  # convierte la tasa deseada de inicio en aprobación de proyectos.

    # Aquí se calcula cuánta inversión aprobada pasa a construcción en cada mes.    
        I_start_t = ID / p.tau_ID_up
            # Depende del stock ID y del tiempo de decisión tau_ID_up.+

    # Aquí se calculan cancelaciones de proyectos.        
        cancellation_rate_t = max(0.0, (ID - desired_ID_t) / p.tau_ID_down)
            # Si hay más decisiones aprobadas de las que se desean, se cancela parte del exceso.}
            # cancela proyectos cuando hay más inversión aprobada de la necesaria.


  # SEÑALES PARA EFICIENCIA DINÁMICA
    # Aquí se calcula presión por baja reserva.
        pressure_reserve_t   = max(0.0, 1.0 - reserve_margin_t / max(p.reserve_ratio, 1e-6))
            # Si la reserva está bajo el nivel requerido, aparece presión para mejorar eficiencia.

    # Aquí se calcula presión por precio alto.        
        pressure_price_t     = max(0.0, demand_response_price_t / p.P_E_base - 1.0)
            # Si el precio que afecta demanda supera el precio base, aparece presión para eficiencia.
            # se mide presión por aumento del precio energético

    # Aquí se calcula déficit de capacidad.        
        capacity_shortfall_t = max(0.0, 1.0 - GP_effective_for_margin_t / max(capacity_target_t, 1e-6))
            # Si la capacidad efectiva es menor que la capacidad objetivo, aparece presión para eficiencia.
            # mide cuánto falta de capacidad respecto al objetivo del sistema

    # Aquí se calcula el objetivo de eficiencia por política pública.        
        IE_policy_target = p.M_IE_policy * policy_ramp(t, p.t_IE, p.tau_IE_policy)
            # Se activa gradualmente desde el mes de política de eficiencia.
            # define la mejora de eficiencia impulsada directamente por política.

    # Aquí se calcula el objetivo total de eficiencia - hace que la eficiencia responda tanto a la política como a presiones internas del sistema
        IE_target_t = IE_policy_target +                            # meta de política,
                    p.k_GP_to_IE      * capacity_shortfall_t +      # presión por falta de capacidad,
                    p.k_price_to_IE   * pressure_price_t +          # presión por precio alto,
                    p.k_reserve_to_IE * pressure_reserve_t          # presión por baja reserva.

    # Aquí se limita el objetivo de eficiencia entre 0 y el máximo permitido.      
        IE_target_t = clamp(IE_target_t, 0.0, p.IE_max)
            # evita que la mejora de eficiencia supere un límite técnico o supuesto máximo

  #  RESERVA ABSOLUTA  
    # Aquí se calcula la reserva absoluta en unidades internas.
    # No es el margen porcentual, sino la diferencia entre capacidad y demanda relevante.
    GR_t = GP - reserve_demand_for_margin_t
    # GR mide cuánta capacidad sobra en términos absolutos, antes de dividir por la demanda.
        

 # INDICES DEL TRILEMA ENERGÉTICO
    # Se define el INDICE DE SEGURIDAD
        Security_Index_t = reserve_margin_t         # Es directamente el margen de reserva.

    # Aquí se calcula el crecimiento relativo del parque generador.
        growth_GP_t = GP / max(GP0, 1e-6)
      # Aquí se calcula el crecimiento relativo de la demanda útil - mide cuánto creció la demanda útil respecto al inicio.  
        growth_S_t  = S_actual_t / max(S0, 1e-6)
      # Aquí se calcula el INDICE DE SOSTENIBILIDAD 
        Sustainability_Index_t = growth_S_t > 1e-9 ? growth_GP_t / growth_S_t : NaN
                # Si la demanda útil crece más rápido que el GP, el índice baja.
                # Si GP crece más que la demanda útil, el índice sube.

    # Aquí se calcula el INDICE DE EQUIDAD            
        Equity_Index_t = isfinite(U_CUW_net_t) && U_CUW_net_t > 0.0 ? UCUW0 / U_CUW_net_t : NaN 
                     # = costo unitario inicial del trabajo útil, costo unitario neto actual 
                     # Si el costo unitario neto baja, la equidad mejora.
                     # Si sube, la equidad empeora.

    # Aquí se calcula una eficiencia ajustada por rebote.
        rebound_adjusted_efficiency_t = F_EC_t > 1e-9 ? p.S_base / F_EC_t : NaN
            # Compara la demanda base contra el consumo energético efectivo.
            # muestra qué tan efectiva fue la eficiencia después de considerar el consumo resultante

    # Aquí se normaliza el CMS respecto al presupuesto energético base mensual.        
        CMS_index_t        = energy_budget_available_t > 0.0 ? CMS / energy_budget_available_t : NaN

    # Aquí se calcula CMS como proporción del ingreso.    
        CMS_income_ratio_t = income_t > 0.0 ? CMS / income_t : NaN

    # Demanda útil en p.u. - normaliza la demanda útil respecto a su valor inicial    
        S_pu_t          = S_actual_t / max(S0, 1e-6)

    # Consumo energético en p.u.    
        FEC_pu_t        = F_EC_t / max(FEC0, 1e-6)
        
    # Parque generador en p.u.    
        GP_pu_t         = GP / max(GP0, 1e-6)

    # Gasto bruto en p.u.    
        EEXP_gross_pu_t = E_Exp_GROSS_t / max(EEXP0, 1e-6)

    # Gasto neto en p.u.
        EEXP_net_pu_t   = E_Exp_NET_t / max(EEXP0, 1e-6)

    # Costo unitario bruto en p.u.   
        UCUW_gross_pu_t = U_CUW_gross_t / max(UCUW0, 1e-6)

    # Costo unitario neto en p.u.    
        UCUW_net_pu_t   = U_CUW_net_t / max(UCUW0, 1e-6)



## BLOQUE DE SALIDA DE LA FUNCION ##
#### Todo lo que se calculó dentro de liberalized_algebraics se empaqueta y se devuelve con nombres.

    return (
        L = L,
        CMS = CMS, CMS_INDEX = CMS_index_t, CMS_INCOME = CMS_income_ratio_t,
        CMS_NORM = CMS_norm_t, INCOME = income_t, INCOME_CLP = income_CLP_t,
        ENERGY_BUDGET_AVAILABLE = energy_budget_available_t,
        GP = GP, GP_PU = GP_pu_t, GP_EFFECTIVE = GP_effective_for_margin_t, GP_EFFECTIVE_PU = GP_effective_for_margin_t / max(GP0, 1e-6), Q = Q,
        ENERGYPRICE = EnergyPrice, ENERGYPRICE_TARGET = energy_price_target_t,
        CONSUMERPRICE = consumer_price_t, DEMANDPRICE = demand_response_price_t,
        PEXP = Pexp, PEXP_TARGET = pexp_target_t,
        SCARCITY_EXPECTATION_PREMIUM = scarcity_expectation_premium_t,
        DEMAND_EXPECTATION_PREMIUM = demand_expectation_premium_t,
        NET_EXPECTED_PRICE_INVESTMENT = net_expected_price_for_investment_t,
        ID = ID, UC = UC, IE = IE, IE_TARGET = IE_target_t, IE_POLICY_TARGET = IE_policy_target,
        EPS = epsilon_t, X_DRE = X_DRE_t, DRE = DRE_t,
        REBOUND_ACTIVATION = rebound_activation_t,
        TAX = tax_val_t, TAX_ACTIVATION = tax_activation_t,
        RESERVE = reserve_margin_t, RESERVE_DEMAND = reserve_demand_for_margin_t,
        PRODCOST = production_cost_t, EXPECTEDFUEL = expected_fuel_price_t,
        EXPECTEDCOST = expected_cost_t, SCARCITY = scarcity_t,
        SURPLUS_DISCOUNT = surplus_discount_t, SPOTPRICE = spot_price_t,
        S = S_actual_t, S_PU = S_pu_t, S_BEFORE_BUDGET = S_before_budget_t, S_BASE_T = S_base_t,
        FEC = F_EC_t, FEC_PU = FEC_pu_t,
        EEXP = E_Exp_t, EEXP_NET = E_Exp_NET_t, EEXP_GROSS = E_Exp_GROSS_t,
        EEXP_NET_PU = EEXP_net_pu_t, EEXP_GROSS_PU = EEXP_gross_pu_t,
        COMPENSATION = compensation_t,
        COMPENSATION_ACTIVATION = compensation_activation_t,
        COMPENSATION_TRIGGER = affordability_trigger_t,
        COMPENSATION_PROFILE = compensation_profile_t,
        TAX_REVENUE = tax_revenue_proxy_t,
        UCUW = U_CUW_net_t, UCUW_GROSS = U_CUW_gross_t, UCUW_NET = U_CUW_net_t,
        UCUW_GROSS_PU = UCUW_gross_pu_t, UCUW_NET_PU = UCUW_net_pu_t,
        UCUW_SIGNAL = U_CUW_signal_t,
        COST_REDUCTION = useful_cost_reduction_t,
        EBS = energy_budget_share_t, BPRESS = budget_pressure_raw_t,
        BCONS = BudgetConstraint, BCONS_TARGET = budget_constraint_target_t,
        BCEFF = budget_constraint_effect_t,
        FSLACK = financial_slack_t, FSTRESS = financial_stress_t,
        IR = IndirectRebound, IR_TARGET = indirect_rebound_target_t,
        RETIREMENTS = retirements_t, COMPLETION = completion_t,
        PI = PI_raw_t, M_PI = m_t, BRAKE = brake_t,
        CAP_TARGET = capacity_target_t, CAP_GAP = capacity_gap_effective_t,
        CAP_GAP_NOMINAL = capacity_gap_nominal_t,
        CAP_TARGET_PLANNING = planning_capacity_target_t,
        CAP_TARGET_RESERVE = reserve_target_capacity_t,
        CAP_GAP_FUTURE = future_capacity_gap_effective_t,
        CAP_GAP_IMMEDIATE = immediate_reserve_gap_effective_t,
        COMMITTED_EXPANSION = committed_expansion_effective_t,
        COMMITTED_EXPANSION_NOMINAL = committed_expansion_nominal_t,
        I_REF = replacement_start_rate_t + expansion_ref_t,
        REPLACEMENT_START = replacement_start_rate_t,
        EXPANSION_REF = expansion_ref_t,
        ENDOGENOUS_EXPANSION = endogenous_expansion_rate_t,
        DEMAND_GROWTH_INVESTMENT = demand_growth_reference_t,
        GAP_CLOSURE_INVESTMENT = gap_closure_reference_t,
        SCARCITY_INVESTMENT = scarcity_reference_t,
        EXPANSION_GATE = expansion_gate_t,
        TAX_EXPANSION_BRAKE = tax_expansion_brake_t,
        APPROVAL = approval_rate_t, ISTART = I_start_t,
        CANCELLATION = cancellation_rate_t,
        DESIRED_START = desired_start_rate_t, DESIRED_ID = desired_ID_t,
        GR = GR_t,
        SEC = Security_Index_t, SECURITY_NORM = Security_Index_t,
        SUST = Sustainability_Index_t,
        EQUITY = Equity_Index_t,
        REB_EFF = rebound_adjusted_efficiency_t,
        PRESSURE_RESERVE = pressure_reserve_t,
        PRESSURE_PRICE = pressure_price_t,
        CAPACITY_SHORTFALL = capacity_shortfall_t,
        DEMAND_GROWTH_FACTOR = demand_growth_factor,
        GROWTH_GP = growth_GP_t,
        GROWTH_S = growth_S_t
    )
end

##############################################################################
# 5. DINAMICA
##############################################################################

function energy_rebound_model_liberalized!(du, u, p, t)
    N = Int(p.N_delay)
    a = liberalized_algebraics(u, p, t)

    flow_rate = p.tau_D > 0.0 ? N / p.tau_D : 1.0e9

    du[1] = flow_rate * (a.X_DRE - a.L[1])
    for i in 2:N
        du[i] = flow_rate * (a.L[i - 1] - a.L[i])
    end

    du[N + 1] = p.EnergyBudget_base - a.EEXP_NET

    du[N + 2]  = a.COMPLETION - a.RETIREMENTS
    du[N + 3]  = (a.S - a.Q) / p.tau_Q
    du[N + 4]  = (a.ENERGYPRICE_TARGET - a.ENERGYPRICE) / p.tau_energy_price
    du[N + 5]  = (a.PEXP_TARGET - a.PEXP) / p.tau_expect_price
    du[N + 6]  = a.APPROVAL - a.ISTART - a.CANCELLATION
    du[N + 7]  = a.ISTART - a.COMPLETION
    du[N + 8]  = (a.IE_TARGET - a.IE) / p.tau_IE
    du[N + 9]  = (a.BCONS_TARGET - a.BCONS) / p.tau_BC
    du[N + 10] = (a.IR_TARGET - a.IR) / p.tau_IR

    return nothing
end

##############################################################################
# 6. RECONSTRUCCION DE SERIES
##############################################################################

function SolLiberalized(sol_ode, p)
    nT = length(sol_ode.t)
    names = [
        :S, :S_PU, :S_BEFORE_BUDGET, :S_BASE_T,
        :FEC, :FEC_PU,
        :EEXP, :EEXP_NET, :EEXP_GROSS, :EEXP_NET_PU, :EEXP_GROSS_PU,
        :COMPENSATION, :COMPENSATION_ACTIVATION, :COMPENSATION_TRIGGER,
        :COMPENSATION_PROFILE, :TAX_REVENUE,
        :UCUW, :UCUW_GROSS, :UCUW_NET, :UCUW_GROSS_PU, :UCUW_NET_PU,
        :UCUW_SIGNAL, :COST_REDUCTION,
        :SPOTPRICE, :ENERGYPRICE, :ENERGYPRICE_TARGET, :CONSUMERPRICE, :DEMANDPRICE,
        :PEXP, :PEXP_TARGET, :SCARCITY_EXPECTATION_PREMIUM, :DEMAND_EXPECTATION_PREMIUM, :NET_EXPECTED_PRICE_INVESTMENT,
        :EPS, :IE, :IE_TARGET, :IE_POLICY_TARGET, :DRE, :REBOUND_ACTIVATION,
        :CMS, :CMS_INDEX, :CMS_INCOME, :CMS_NORM, :INCOME, :INCOME_CLP, :ENERGY_BUDGET_AVAILABLE,
        :GP, :GP_PU, :GP_EFFECTIVE, :GP_EFFECTIVE_PU, :Q, :ID, :UC, :ISTART,
        :APPROVAL, :CANCELLATION, :COMPLETION, :RETIREMENTS,
        :RESERVE, :RESERVE_DEMAND, :REQRES,
        :PRODCOST, :EXPECTEDFUEL, :EXPECTEDCOST,
        :SCARCITY, :SURPLUS_DISCOUNT,
        :PI, :M_PI, :BRAKE, :GR, :SEC, :SECURITY_NORM, :SUST, :EQUITY, :REB_EFF,
        :DESIRED_START, :DESIRED_ID, :CAP_TARGET, :CAP_GAP, :CAP_GAP_NOMINAL,
        :CAP_TARGET_PLANNING, :CAP_TARGET_RESERVE, :CAP_GAP_FUTURE, :CAP_GAP_IMMEDIATE, :I_REF,
        :REPLACEMENT_START, :EXPANSION_REF, :ENDOGENOUS_EXPANSION,
        :DEMAND_GROWTH_INVESTMENT, :GAP_CLOSURE_INVESTMENT, :SCARCITY_INVESTMENT,
        :EXPANSION_GATE, :TAX_EXPANSION_BRAKE,
        :COMMITTED_EXPANSION, :COMMITTED_EXPANSION_NOMINAL,
        :TAX, :TAX_ACTIVATION, :PRESSURE_RESERVE, :PRESSURE_PRICE, :CAPACITY_SHORTFALL,
        :DEMAND_GROWTH_FACTOR, :GROWTH_GP, :GROWTH_S,
        :EBS, :BPRESS, :BCONS, :BCONS_TARGET, :BCEFF,
        :FSLACK, :FSTRESS, :IR, :IR_TARGET
    ]

    series = Dict{Symbol, Vector{Float64}}()
    for name in names
        series[name] = zeros(nT)
    end

    for i in 1:nT
        u = sol_ode.u[i]
        t = sol_ode.t[i]
        a = liberalized_algebraics(u, p, t)
        for name in names
            if name == :REQRES
                series[name][i] = p.reserve_ratio
            else
                series[name][i] = getfield(a, name)
            end
        end
    end

    return NamedTuple{Tuple(names)}(Tuple(series[name] for name in names))
end

##############################################################################
# 7. SIMULACION
##############################################################################

tspan_lib = (0.0, p_lib.t_final_months)
times_lib = range(tspan_lib[1], tspan_lib[2], step = p_lib.dt_save_months)

u0_lib  = initial_state_liberalized(p_lib)
prob_lib = ODEProblem(energy_rebound_model_liberalized!, u0_lib, tspan_lib, p_lib)
sol_lib  = solve(prob_lib, Tsit5(), reltol = 1e-8, abstol = 1e-8, saveat = times_lib)
ser_lib  = SolLiberalized(sol_lib, p_lib)

u0_lib_tax  = initial_state_liberalized(p_lib_tax)
prob_lib_tax = ODEProblem(energy_rebound_model_liberalized!, u0_lib_tax, tspan_lib, p_lib_tax)
sol_lib_tax  = solve(prob_lib_tax, Tsit5(), reltol = 1e-8, abstol = 1e-8, saveat = times_lib)
ser_lib_tax  = SolLiberalized(sol_lib_tax, p_lib_tax)

u0_lib_comp  = initial_state_liberalized(p_lib_comp)
prob_lib_comp = ODEProblem(energy_rebound_model_liberalized!, u0_lib_comp, tspan_lib, p_lib_comp)
sol_lib_comp  = solve(prob_lib_comp, Tsit5(), reltol = 1e-8, abstol = 1e-8, saveat = times_lib)
ser_lib_comp  = SolLiberalized(sol_lib_comp, p_lib_comp)

##############################################################################
# 8. DIAGNOSTICO
##############################################################################

function print_scenario_final(label, ser)
    gp_mw = ser.GP_PU[end] * p_lib.GP_base_physical_MW
    delta_gp_mw = (ser.GP_PU[end] - 1.0) * p_lib.GP_base_physical_MW
    horizon_years = max(p_lib.t_final_months / 12.0, 1e-6)
    delta_gp_mw_year = delta_gp_mw / horizon_years
    delta_gp_pct_total = (ser.GP_PU[end] - 1.0) * 100.0
    delta_gp_pct_year = delta_gp_pct_total / horizon_years

    println(label)
    println("  S [p.u.]              = ", round(ser.S_PU[end], digits = 4))
    println("  FEC [p.u.]            = ", round(ser.FEC_PU[end], digits = 4))
    println("  GP [p.u.]             = ", round(ser.GP_PU[end], digits = 4))
    println("  GP [MW]               = ", round(gp_mw, digits = 0))
    println("  Delta GP [MW]         = ", round(delta_gp_mw, digits = 0))
    println("  Delta GP [MW/year]    = ", round(delta_gp_mw_year, digits = 0))
    println("  Delta GP [% 5y]       = ", round(delta_gp_pct_total, digits = 2))
    println("  Delta GP [%/year]     = ", round(delta_gp_pct_year, digits = 2))
    println("  GP effective [p.u.]   = ", round(ser.GP_EFFECTIVE_PU[end], digits = 4))
    println("  IE state              = ", round(ser.IE[end], digits = 4))
    println("  DRE [%]               = ", round(ser.DRE[end] * 100, digits = 2))
    println("  IR  [%]               = ", round(ser.IR[end] * 100, digits = 2))
    println("  UCUW net [p.u.]       = ", round(ser.UCUW_NET_PU[end], digits = 4))
    println("  EEXP net [p.u.]       = ", round(ser.EEXP_NET_PU[end], digits = 4))
    println("  Compensation [p.u.]   = ", round(ser.COMPENSATION[end], digits = 4))
    println("  CMS [p.u.]            = ", round(ser.CMS_INDEX[end], digits = 4))
    println("  CMS norm [p.u.]       = ", round(ser.CMS_NORM[end], digits = 4))
    println("  Financial Slack       = ", round(ser.FSLACK[end], digits = 4))
    println("  Financial Stress      = ", round(ser.FSTRESS[end], digits = 4))
    println("  EBS                   = ", round(ser.EBS[end], digits = 4))
    println("  BCons                 = ", round(ser.BCONS[end], digits = 4))
    println("  Reserve Margin [p.u.] = ", round(ser.RESERVE[end], digits = 4))
    println("  Reserve Margin [%]    = ", round(ser.RESERVE[end] * 100, digits = 2))
    println("  Security Index        = ", round(ser.SEC[end], digits = 4))
    println("  Sustainability Index  = ", round(ser.SUST[end], digits = 4))
    println("  Equity Index          = ", round(ser.EQUITY[end], digits = 4))
    println("  PI                    = ", round(ser.PI[end], digits = 4))
    println("")
end

function model_unit_to_mw_factor(p)
    GP0 = (p.S_base / p.epsilon_0) * (1.0 + p.reserve_ratio)
    return p.GP_base_physical_MW / max(GP0, 1e-6)
end

function print_rebound_diagnostic(label, ser)
    total_rebound = (1.0 + ser.DRE[end]) * (1.0 + ser.IR[end]) - 1.0
    backfire_flag = ser.FEC_PU[end] > 1.0 ? "YES" : "NO"

    println(label)
    println("  DRE final [%]              = ", round(ser.DRE[end] * 100, digits = 2))
    println("  IR final [%]               = ", round(ser.IR[end] * 100, digits = 2))
    println("  Rebound total aprox. [%]   = ", round(total_rebound * 100, digits = 2))
    println("  Cost reduction signal      = ", round(ser.COST_REDUCTION[end], digits = 4))
    println("  Financial Slack            = ", round(ser.FSLACK[end], digits = 4))
    println("  Financial Stress           = ", round(ser.FSTRESS[end], digits = 4))
    println("  S final [p.u.]             = ", round(ser.S_PU[end], digits = 4))
    println("  FEC final [p.u.]           = ", round(ser.FEC_PU[end], digits = 4))
    println("  Backfire detected          = ", backfire_flag)
    println("")
end

function time_of_extreme(ser, v; mode = :max)
    idx = mode == :min ? argmin(v) : argmax(v)
    return round(times_lib[idx], digits = 1), v[idx]
end

function print_metric_line(name, values; digits = 4, unit = "")
    idx_max = argmax(values)
    idx_min = argmin(values)
    final_v = values[end]
    max_v = values[idx_max]
    min_v = values[idx_min]
    println("  ", rpad(name, 30),
        " final = ", round(final_v, digits = digits), unit,
        " | max = ", round(max_v, digits = digits), unit, " at m ", round(times_lib[idx_max], digits = 1),
        " | min = ", round(min_v, digits = digits), unit, " at m ", round(times_lib[idx_min], digits = 1))
end

function print_main_variable_diagnostic(label, ser, p)
    base_budget_clp = p.Income_base_CLP * p.energy_budget_share_target
    efficiency_index = 1.0 .+ ser.IE
    gp_mw = ser.GP_PU .* p.GP_base_physical_MW
    fec_gwh_month = ser.FEC_PU .* p.FEC_base_physical_GWh_month
    cms_clp = ser.CMS_INDEX .* base_budget_clp

    println(label)
    print_metric_line("Efficiency index", efficiency_index)
    print_metric_line("S [p.u.]", ser.S_PU)
    print_metric_line("FEC [p.u.]", ser.FEC_PU)
    print_metric_line("FEC [GWh/month]", fec_gwh_month; digits = 0, unit = " GWh/month")
    print_metric_line("UCUW net [p.u.]", ser.UCUW_NET_PU)
    print_metric_line("Energy expenditure [p.u.]", ser.EEXP_NET_PU)
    print_metric_line("CMS [p.u.]", ser.CMS_INDEX)
    print_metric_line("CMS [CLP acumulado]", cms_clp; digits = 0, unit = " CLP")
    print_metric_line("Equity index", ser.EQUITY)
    print_metric_line("GP [p.u.]", ser.GP_PU)
    print_metric_line("GP [MW]", gp_mw; digits = 0, unit = " MW")
    print_metric_line("Reserve margin", ser.RESERVE)
    print_metric_line("Security index", ser.SEC)
    print_metric_line("Sustainability index", ser.SUST)
    print_metric_line("Profitability index", ser.PI)
    println("")
end

function print_olsina_diagnostic(label, ser, p)
    mw_factor = model_unit_to_mw_factor(p)
    cap_base = max(p.S_base, 1e-6)

    id_mw = ser.ID[end] * mw_factor
    uc_mw = ser.UC[end] * mw_factor
    desired_start_mw_year = ser.DESIRED_START[end] * mw_factor * 12.0
    approval_mw_year = ser.APPROVAL[end] * mw_factor * 12.0
    construction_start_mw_year = ser.ISTART[end] * mw_factor * 12.0
    completion_mw_year = ser.COMPLETION[end] * mw_factor * 12.0
    retirement_mw_year = ser.RETIREMENTS[end] * mw_factor * 12.0
    expansion_ref_mw_year = ser.EXPANSION_REF[end] * mw_factor * 12.0
    endogenous_expansion_mw_year = ser.ENDOGENOUS_EXPANSION[end] * mw_factor * 12.0
    demand_growth_inv_mw_year = ser.DEMAND_GROWTH_INVESTMENT[end] * mw_factor * 12.0
    gap_closure_inv_mw_year = ser.GAP_CLOSURE_INVESTMENT[end] * mw_factor * 12.0
    scarcity_inv_mw_year = ser.SCARCITY_INVESTMENT[end] * mw_factor * 12.0
    committed_expansion_mw = ser.COMMITTED_EXPANSION[end] * mw_factor
    committed_expansion_nominal_mw = ser.COMMITTED_EXPANSION_NOMINAL[end] * mw_factor

    t_min_reserve, min_reserve = time_of_extreme(ser, ser.RESERVE; mode = :min)
    t_max_scarcity, max_scarcity = time_of_extreme(ser, ser.SCARCITY; mode = :max)
    t_max_gap, max_gap = time_of_extreme(ser, ser.CAP_GAP; mode = :max)
    t_max_expansion, max_expansion = time_of_extreme(ser, ser.ENDOGENOUS_EXPANSION; mode = :max)

    println(label)
    println("  PI                                  = ", round(ser.PI[end], digits = 4))
    println("  Investment multiplier               = ", round(ser.M_PI[end], digits = 4))
    println("  Energy price [p.u.]                 = ", round(ser.ENERGYPRICE[end], digits = 4))
    println("  Expected price [p.u.]               = ", round(ser.PEXP[end], digits = 4))
    println("  Net expected price investment [p.u.] = ", round(ser.NET_EXPECTED_PRICE_INVESTMENT[end], digits = 4))
    println("  Expected cost [p.u.]                = ", round(ser.EXPECTEDCOST[end], digits = 4))
    println("  Spot price [p.u.]                   = ", round(ser.SPOTPRICE[end], digits = 4))
    println("  Scarcity signal final               = ", round(ser.SCARCITY[end], digits = 4))
    println("  Scarcity signal max                 = ", round(max_scarcity, digits = 4), " at month ", t_max_scarcity)
    println("  Surplus discount                    = ", round(ser.SURPLUS_DISCOUNT[end], digits = 4))
    println("  Reserve margin final [p.u.]         = ", round(ser.RESERVE[end], digits = 4))
    println("  Reserve margin min [p.u.]           = ", round(min_reserve, digits = 4), " at month ", t_min_reserve)
    println("  Required capacity [p.u.]            = ", round(ser.RESERVE_DEMAND[end] / cap_base, digits = 4))
    println("  Capacity target [p.u.]              = ", round(ser.CAP_TARGET[end] / cap_base, digits = 4))
    println("  Capacity gap final [p.u.]           = ", round(ser.CAP_GAP[end] / cap_base, digits = 4))
    println("  Capacity gap max [p.u.]             = ", round(max_gap / cap_base, digits = 4), " at month ", t_max_gap)
    println("  Immediate gap [p.u.]                = ", round(ser.CAP_GAP_IMMEDIATE[end] / cap_base, digits = 4))
    println("  Future gap [p.u.]                   = ", round(ser.CAP_GAP_FUTURE[end] / cap_base, digits = 4))
    println("  GP nominal [p.u.]                   = ", round(ser.GP_PU[end], digits = 4))
    println("  GP effective [p.u.]                 = ", round(ser.GP_EFFECTIVE_PU[end], digits = 4))
    println("  Marginal capacity credit            = ", round(p.marginal_capacity_credit, digits = 4))
    println("  Approval pipeline ID [MW]           = ", round(id_mw, digits = 0))
    println("  Under construction UC [MW]          = ", round(uc_mw, digits = 0))
    println("  Desired start [MW/year]             = ", round(desired_start_mw_year, digits = 0))
    println("  Approval rate [MW/year]             = ", round(approval_mw_year, digits = 0))
    println("  Construction start [MW/year]        = ", round(construction_start_mw_year, digits = 0))
    println("  Capacity completion [MW/year]       = ", round(completion_mw_year, digits = 0))
    println("  Retirement [MW/year]                = ", round(retirement_mw_year, digits = 0))
    println("  Exogenous GP growth [MW/year]       = 0.0")
    println("  Reference investment [MW/year]      = ", round(expansion_ref_mw_year, digits = 0))
    println("  Endogenous expansion [MW/year]      = ", round(endogenous_expansion_mw_year, digits = 0))
    println("  Endogenous expansion max [MW/year]  = ", round(max_expansion * mw_factor * 12.0, digits = 0), " at month ", t_max_expansion)
    println("  Demand-growth investment [MW/year]  = ", round(demand_growth_inv_mw_year, digits = 0))
    println("  Capacity-gap investment [MW/year]   = ", round(gap_closure_inv_mw_year, digits = 0))
    println("  Scarcity investment [MW/year]       = ", round(scarcity_inv_mw_year, digits = 0))
    println("  Committed pipeline effective [MW]   = ", round(committed_expansion_mw, digits = 0))
    println("  Committed pipeline nominal [MW]     = ", round(committed_expansion_nominal_mw, digits = 0))
    println("  Expansion gate                      = ", round(ser.EXPANSION_GATE[end], digits = 4))
    println("  Tax expansion brake                 = ", round(ser.TAX_EXPANSION_BRAKE[end], digits = 4))
    println("  Over-reserve brake                  = ", round(ser.BRAKE[end], digits = 4))
    println("")
end

println("══════════════════════════════════════════════════════════")
println("  DIAGNOSTICO MODELO ENDOGENO")
println("══════════════════════════════════════════════════════════")
println("Income_base CLP             = ", round(p_lib.Income_base_CLP, digits = 0))
println("Income growth               = 0.0  (ingreso fijo)")
println("P_T escenario impuesto       = ", p_lib_tax.P_T)
println("Compensacion                 = reciclaje parcial focalizado del impuesto")
println("Recycling share              = ", p_lib_comp.compensation_recycling_share)
println("GP growth driver            = endogenous investment + committed pipeline + retirements")
println("Horizonte                    = ", p_lib.t_final_months, " meses")
println("")

print_scenario_final("Sin impuesto final:", ser_lib)
print_scenario_final("Con impuesto final:", ser_lib_tax)
print_scenario_final("Con impuesto + compensacion final:", ser_lib_comp)

println("DIAGNOSTICO FINAL, PEAK Y MINIMO DE VARIABLES GRAFICADAS:")
print_main_variable_diagnostic("Sin impuesto:", ser_lib, p_lib)
print_main_variable_diagnostic("Con impuesto:", ser_lib_tax, p_lib_tax)
print_main_variable_diagnostic("Con impuesto + compensacion:", ser_lib_comp, p_lib_comp)

println("DIAGNOSTICO REBOTE Y BACKFIRE:")
print_rebound_diagnostic("Sin impuesto:", ser_lib)
print_rebound_diagnostic("Con impuesto:", ser_lib_tax)
print_rebound_diagnostic("Con impuesto + compensacion:", ser_lib_comp)

println("DIAGNOSTICO BLOQUE OLSINA / INVERSION-CAPACIDAD:")
print_olsina_diagnostic("Sin impuesto:", ser_lib, p_lib)
print_olsina_diagnostic("Con impuesto:", ser_lib_tax, p_lib_tax)
print_olsina_diagnostic("Con impuesto + compensacion:", ser_lib_comp, p_lib_comp)

println("Deltas With Tax menos No Tax:")
println("  Delta S [p.u.]             = ", round(ser_lib_tax.S_PU[end] - ser_lib.S_PU[end], digits = 4))
println("  Delta FEC [p.u.]           = ", round(ser_lib_tax.FEC_PU[end] - ser_lib.FEC_PU[end], digits = 4))
println("  Delta UCUW net [p.u.]      = ", round(ser_lib_tax.UCUW_NET_PU[end] - ser_lib.UCUW_NET_PU[end], digits = 4))
println("  Delta Security             = ", round(ser_lib_tax.SEC[end] - ser_lib.SEC[end], digits = 4))
println("  Delta Sustainability       = ", round(ser_lib_tax.SUST[end] - ser_lib.SUST[end], digits = 4))
println("  Delta Equity               = ", round(ser_lib_tax.EQUITY[end] - ser_lib.EQUITY[end], digits = 4))
println("  Delta CMS [p.u.]           = ", round(ser_lib_tax.CMS_INDEX[end] - ser_lib.CMS_INDEX[end], digits = 4))
println("")
println("Deltas With Tax + Comp menos With Tax:")
println("  Delta EEXP net [p.u.]      = ", round(ser_lib_comp.EEXP_NET_PU[end] - ser_lib_tax.EEXP_NET_PU[end], digits = 4))
println("  Delta UCUW net [p.u.]      = ", round(ser_lib_comp.UCUW_NET_PU[end] - ser_lib_tax.UCUW_NET_PU[end], digits = 4))
println("  Delta Equity               = ", round(ser_lib_comp.EQUITY[end] - ser_lib_tax.EQUITY[end], digits = 4))
println("  Delta FEC [p.u.]           = ", round(ser_lib_comp.FEC_PU[end] - ser_lib_tax.FEC_PU[end], digits = 4))
println("  Delta CMS [p.u.]           = ", round(ser_lib_comp.CMS_INDEX[end] - ser_lib_tax.CMS_INDEX[end], digits = 4))
println("══════════════════════════════════════════════════════════")

##############################################################################
# 9. PALETA
##############################################################################

lib_c_main   = "#264653"
lib_c_fill   = "#1D3557"
lib_c_index  = "#3B0910"
lib_c_green  = "#2A9D8F"
lib_c_amber  = "#E9C46A"
lib_c_gray   = "#6C757D"
lib_c_orange = "#E76F51"
lib_c_purple = "#7B2CBF"

##############################################################################
# 10. FUNCIONES DE GRAFICACION FINAL
##############################################################################


function phase_event_shapes_three_panel()
    events = [p_lib.t_IE, p_lib.t_DRE, p_lib.t_T, p_lib.t_compensation]
    xrefs = ["x", "x2", "x3"]
    domains = [(0.72, 1.00), (0.38, 0.66), (0.04, 0.32)]
    shapes = Any[]
    for (xref, domain) in zip(xrefs, domains)
        for ev in events
            push!(shapes, attr(
                type = "line",
                xref = xref,
                yref = "paper",
                x0 = ev,
                x1 = ev,
                y0 = domain[1],
                y1 = domain[2],
                line = attr(color = "rgba(80,80,80,0.35)", width = 1.1, dash = "dot")
            ))
        end
    end
    return shapes
end

function phase_event_shapes_single()
    events = [p_lib.t_IE, p_lib.t_DRE, p_lib.t_T, p_lib.t_compensation]
    shapes = Any[]
    for ev in events
        push!(shapes, attr(
            type = "line",
            xref = "x",
            yref = "paper",
            x0 = ev,
            x1 = ev,
            y0 = 0.0,
            y1 = 1.0,
            line = attr(color = "rgba(80,80,80,0.35)", width = 1.1, dash = "dot")
        ))
    end
    return shapes
end

function make_three_panel_figure(t, y_no_tax, y_tax, y_comp;
                                 title::String,
                                 ytitle::String,
                                 trace_name::String,
                                 filename::String,
                                 color = "#264653",
                                 baseline = nothing,
                                 baseline_name = "Base = 1",
                                 baseline_color = "#6C757D",
                                 baseline_dash = "dot",
                                 width = 980,
                                 height = 760,
                                 pad = 0.08)

    baseline_vec = baseline === nothing ? Float64[] : [Float64(baseline)]
    y_range = auto_range(y_no_tax, y_tax, y_comp, baseline_vec; pad = pad, minwidth = 0.02)

    ys = [y_no_tax, y_tax, y_comp]
    panels = ["No Tax", "With Tax", "With Tax + Compensation"]

    first_trace = line_trace(t, ys[1], trace_name, color;
                             xaxis = "",
                             yaxis = "",
                             fillopacity = 0.18,
                             showlegend = true)
    traces = [first_trace]

    if baseline !== nothing
        push!(traces, horizontal_line(t, baseline, baseline_name, baseline_color;
                                      xaxis = "",
                                      yaxis = "",
                                      dash = baseline_dash,
                                      width = 1.8,
                                      showlegend = true))
    end

    for i in 2:3
        xaxis_name = "x" * string(i)
        yaxis_name = "y" * string(i)
        push!(traces, line_trace(t, ys[i], trace_name, color;
                                 xaxis = xaxis_name,
                                 yaxis = yaxis_name,
                                 fillopacity = 0.18,
                                 showlegend = false))
        if baseline !== nothing
            push!(traces, horizontal_line(t, baseline, baseline_name, baseline_color;
                                          xaxis = xaxis_name,
                                          yaxis = yaxis_name,
                                          dash = baseline_dash,
                                          width = 1.8,
                                          showlegend = false))
        end
    end

    lyt = Layout(
        title = title,
        height = height,
        width = width,
        paper_bgcolor = "white",
        plot_bgcolor = "white",

        xaxis  = attr(domain = [0.10, 0.92], anchor = "y",  showticklabels = false, showgrid = true, gridcolor = "lightgray", showline = true, linecolor = "black"),
        yaxis  = attr(domain = [0.72, 1.00], title = "", showline = true, linecolor = "black", range = y_range),

        xaxis2 = attr(domain = [0.10, 0.92], anchor = "y2", showticklabels = false, showgrid = true, gridcolor = "lightgray", showline = true, linecolor = "black"),
        yaxis2 = attr(domain = [0.38, 0.66], title = ytitle, showline = true, linecolor = "black", range = y_range),

        xaxis3 = attr(domain = [0.10, 0.92], anchor = "y3", title = "Time (months)", showgrid = true, gridcolor = "lightgray", showline = true, linecolor = "black"),
        yaxis3 = attr(domain = [0.04, 0.32], title = "", showline = true, linecolor = "black", range = y_range),

        annotations = [
            attr(text = "(a) " * panels[1], x = 0.51, y = 1.04, xref = "paper", yref = "paper", showarrow = false, font = attr(size = 13, color = lib_c_index)),
            attr(text = "(b) " * panels[2], x = 0.51, y = 0.695, xref = "paper", yref = "paper", showarrow = false, font = attr(size = 13, color = lib_c_index)),
            attr(text = "(c) " * panels[3], x = 0.51, y = 0.355, xref = "paper", yref = "paper", showarrow = false, font = attr(size = 13, color = lib_c_index))
        ],
        shapes = phase_event_shapes_three_panel(),
        legend = attr(orientation = "h", x = 0.5, y = -0.14, xanchor = "center", yanchor = "center")
    )

    fig = Plot(traces, lyt)
    if DISPLAY_FIGURES
        display(fig)
    end
    save_plot(fig, filename; width = width, height = height)
    return nothing
end

function make_efficiency_figure(t, ie_policy_target, ie_state)
    y_policy = 1.0 .+ ie_policy_target
    y_state  = 1.0 .+ ie_state
    y_range = auto_range(y_policy, y_state, [1.0]; pad = 0.10, minwidth = 0.03)

    tr_policy = line_trace(t, y_policy, "IE policy target", lib_c_index; width = 2.8, fillopacity = 0.10)
    tr_state  = line_trace(t, y_state,  "IE state",         lib_c_main;  width = 2.8, fillopacity = 0.10)
    tr_base   = horizontal_line(t, 1.0, "Base = 1", lib_c_gray; dash = "dot", width = 1.8)

    lyt = Layout(
        title = "Dynamic Efficiency Improvement",
        height = 520,
        width = 980,
        paper_bgcolor = "white",
        plot_bgcolor = "white",
        xaxis = attr(title = "Time (months)", showgrid = true, gridcolor = "lightgray", showline = true, linecolor = "black"),
        yaxis = attr(title = "Efficiency index [p.u.]", showline = true, linecolor = "black", range = y_range),
        shapes = phase_event_shapes_single(),
        legend = attr(orientation = "h", x = 0.5, y = -0.20, xanchor = "center", yanchor = "center")
    )

    fig = Plot([tr_policy, tr_state, tr_base], lyt)
    if DISPLAY_FIGURES
        display(fig)
    end
    save_plot(fig, "Fig_Final_Dynamic_Efficiency"; width = 980, height = 520)
    return nothing
end

##############################################################################
# 11. FIGURAS FINALES
##############################################################################

fig_eff = make_efficiency_figure(sol_lib.t, ser_lib.IE_POLICY_TARGET, ser_lib.IE)

fig_S = make_three_panel_figure(sol_lib.t, ser_lib.S_PU, ser_lib_tax.S_PU, ser_lib_comp.S_PU;
    title = "Useful Work Demand S",
    ytitle = "S [p.u.]",
    trace_name = "S",
    filename = "Fig_Final_S_3Panels",
    color = lib_c_green,
    baseline = 1.0,
    baseline_name = "Base = 1"
)

fig_FEC = make_three_panel_figure(sol_lib.t, ser_lib.FEC_PU, ser_lib_tax.FEC_PU, ser_lib_comp.FEC_PU;
    title = "Flow of Energy Consumption FEC",
    ytitle = "FEC [p.u.]",
    trace_name = "FEC",
    filename = "Fig_Final_FEC_3Panels",
    color = lib_c_main,
    baseline = 1.0,
    baseline_name = "Base = 1"
)

fig_UCUW = make_three_panel_figure(sol_lib.t, ser_lib.UCUW_NET_PU, ser_lib_tax.UCUW_NET_PU, ser_lib_comp.UCUW_NET_PU;
    title = "Unitary Cost of Useful Work",
    ytitle = "UCUW net [p.u.]",
    trace_name = "UCUW net",
    filename = "Fig_Final_UCUW_3Panels",
    color = lib_c_index,
    baseline = 1.0,
    baseline_name = "Base = 1"
)

fig_EEXP = make_three_panel_figure(sol_lib.t, ser_lib.EEXP_NET_PU, ser_lib_tax.EEXP_NET_PU, ser_lib_comp.EEXP_NET_PU;
    title = "Energy Expenditure",
    ytitle = "Net expenditure [p.u.]",
    trace_name = "Net energy expenditure",
    filename = "Fig_Final_Energy_Expenditure_3Panels",
    color = lib_c_orange,
    baseline = 1.0,
    baseline_name = "Base = 1"
)

fig_SEC = make_three_panel_figure(sol_lib.t, ser_lib.SEC, ser_lib_tax.SEC, ser_lib_comp.SEC;
    title = "Energy Security Index",
    ytitle = "Reserve margin [p.u.]",
    trace_name = "Security index",
    filename = "Fig_Final_Security_Index_3Panels",
    color = lib_c_fill,
    baseline = p_lib.reserve_ratio,
    baseline_name = "Required reserve = 0.15"
)

fig_SUST = make_three_panel_figure(sol_lib.t, ser_lib.SUST, ser_lib_tax.SUST, ser_lib_comp.SUST;
    title = "Energy Sustainability Index",
    ytitle = "Sustainability index [p.u.]",
    trace_name = "Sustainability index",
    filename = "Fig_Final_Sustainability_Index_3Panels",
    color = lib_c_green,
    baseline = 1.0,
    baseline_name = "Base = 1"
)

fig_EQUITY = make_three_panel_figure(sol_lib.t, ser_lib.EQUITY, ser_lib_tax.EQUITY, ser_lib_comp.EQUITY;
    title = "Energy Equity Index",
    ytitle = "Equity index [p.u.]",
    trace_name = "Equity index",
    filename = "Fig_Final_Equity_Index_3Panels",
    color = lib_c_purple,
    baseline = 1.0,
    baseline_name = "Base = 1"
)

fig_CMS = make_three_panel_figure(sol_lib.t, ser_lib.CMS_INDEX, ser_lib_tax.CMS_INDEX, ser_lib_comp.CMS_INDEX;
    title = "Cumulative Monetary Savings",
    ytitle = "CMS [p.u.]",
    trace_name = "CMS",
    filename = "Fig_Final_CMS_3Panels",
    color = lib_c_purple,
    baseline = 0.0,
    baseline_name = "CMS = 0",
    pad = 0.10
)

fig_GP = make_three_panel_figure(sol_lib.t, ser_lib.GP_PU, ser_lib_tax.GP_PU, ser_lib_comp.GP_PU;
    title = "Generation Park",
    ytitle = "GP [p.u.]",
    trace_name = "Generation Park",
    filename = "Fig_Final_Generation_Park_3Panels",
    color = lib_c_fill,
    baseline = 1.0,
    baseline_name = "Base = 1"
)

fig_RESERVE = make_three_panel_figure(sol_lib.t, ser_lib.RESERVE, ser_lib_tax.RESERVE, ser_lib_comp.RESERVE;
    title = "Generation Reserve Margin",
    ytitle = "Reserve margin [p.u.]",
    trace_name = "Reserve margin",
    filename = "Fig_Final_Generation_Reserve_3Panels",
    color = lib_c_index,
    baseline = p_lib.reserve_ratio,
    baseline_name = "Required reserve = 0.15"
)

fig_PI = make_three_panel_figure(sol_lib.t, ser_lib.PI, ser_lib_tax.PI, ser_lib_comp.PI;
    title = "Profitability Index",
    ytitle = "PI [p.u.]",
    trace_name = "PI",
    filename = "Fig_Final_Profitability_Index_3Panels",
    color = lib_c_amber,
    baseline = 1.0,
    baseline_name = "PI = 1"
)

println("\nGraficas guardadas en:")
println(RESULTS_DIR_LIB)
println("Archivos principales generados:")
println("  Fig_Final_Dynamic_Efficiency")
println("  Fig_Final_S_3Panels")
println("  Fig_Final_FEC_3Panels")
println("  Fig_Final_UCUW_3Panels")
println("  Fig_Final_Energy_Expenditure_3Panels")
println("  Fig_Final_Security_Index_3Panels")
println("  Fig_Final_Sustainability_Index_3Panels")
println("  Fig_Final_Equity_Index_3Panels")
println("  Fig_Final_CMS_3Panels")
println("  Fig_Final_Generation_Park_3Panels")
println("  Fig_Final_Generation_Reserve_3Panels")
println("  Fig_Final_Profitability_Index_3Panels")

##############################################################################
# FIN DEL CODIGO
##############################################################################
