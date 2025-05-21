import math

def calculate_conduit_trial_flow_momentum_eq(
    q_old_barrel: float,                # 이전 시간 스텝에서의 배럴당 유량 (cfs)
    q_last_iter_barrel: float,          # 이전 반복에서의 배럴당 유량 (cfs) - 현재 속도 계산에 사용
    h1_node: float,                     # 상류 노드의 수리 수두 (ft)
    h2_node: float,                     # 하류 노드의 수리 수두 (ft)
    z1_conduit: float,                  # 도관 상류단 관저고 (ft) 
    z2_conduit: float,                  # 도관 하류단 관저고 (ft) 
    a_old_barrel: float,                # 이전 시간 스텝에서의 배럴당 단면적 (ft²)
    length_modified: float,             # Courant 수정된 도관 길이 (ft)
    length_actual: float,               # 도관의 실제 물리적 길이 (ft)
    conduit_roughness_factor: float,    # 도관 조도 계수 (SWMM의 Conduit[k].roughFactor 값)
    dt: float,                          # 시간 스텝 (s)
    sigma_inertial_damping: float,      # 관성 감쇠 계수 (0 ~ 1)
    rho_factor_upstream_weighting: float, # 마찰/압력 항에서 면적 및 동수반경에 대한 상류 가중 계수
    has_local_losses: bool,             # 국부 손실 적용 여부 (True/False)
    coeff_loss_inlet: float,            # 입구 국부 손실 계수 (단위 없음, K)
    coeff_loss_outlet: float,           # 출구 국부 손실 계수 (단위 없음, K)
    coeff_loss_avg: float,              # 평균 국부 손실 계수 (단위 없음, K)
    GRAVITY: float,                     # 중력 가속도 (예: 32.2 ft/s²)
    FUDGE: float,                       # 0으로 나누기 방지를 위한 작은 수 (예: 1e-6)
    y1_current_iter: float,             # 현재 반복에서의 도관 상류단 수심 (ft) (h1_node - z1_conduit)
    y2_current_iter: float,             # 현재 반복에서의 도관 하류단 수심 (ft) (h2_node - z2_conduit)
    a1_current_iter: float,             # 현재 반복에서의 도관 상류단 단면적 (ft²)
    a2_current_iter: float,             # 현재 반복에서의 도관 하류단 단면적 (ft²)
    a_mid_current_iter: float,          # 현재 반복에서의 도관 중간 단면적 (ft²)
    r1_current_iter: float,             # 현재 반복에서의 도관 상류단 동수반경 (ft)
    r_mid_current_iter: float           # 현재 반복에서의 도관 중간 동수반경 (ft)
) -> dict:
    """
    SWMM의 `dwflow_findConduitFlow` 함수 내 핵심 운동량 방정식 로직을 구현하여
    동적파 도관 유량에 대한 시험 유량을 계산하고, 중간 계산값들을 반환합니다. (배럴 1개 기준)
    이 함수는 Preissmann Slot 로직이 적용된 후의 유효 단면적/동수반경을 입력받는다고 가정합니다.

    반환값:
        딕셔너리:
            'q_trial_momentum_barrel': 시험 유량 (cfs)
            'v_current_iter_barrel': 현재 반복 유속 (ft/s)
            'aWtd': 가중 단면적 (ft^2)
            'rWtd': 가중 동수반경 (ft)
            'dq1' - 'dq6': 운동량 방정식의 각 항
            'denominator': 최종 유량 계산식의 분모
    """

    # --- 1. 현재 반복에서의 배럴당 유속 계산 (v_current_iter_barrel) ---
    if a_mid_current_iter > FUDGE:
        v_current_iter_barrel = q_last_iter_barrel / a_mid_current_iter
    else:
        v_current_iter_barrel = 0.0

    # --- 2. 가중 단면적 (aWtd) 계산 ---
    aWtd = a1_current_iter + (a_mid_current_iter - a1_current_iter) * rho_factor_upstream_weighting
    aWtd = max(aWtd, FUDGE) 

    # --- 3. 가중 동수반경 (rWtd) 계산 ---
    rWtd = r1_current_iter + (r_mid_current_iter - r1_current_iter) * rho_factor_upstream_weighting
    rWtd = max(rWtd, FUDGE) 

    # --- 4. dq1 (마찰 경사 항) 계산 ---
    if rWtd > FUDGE:
        try:
            rWtd_pow = rWtd**1.3333333333333333 # 4/3
            if abs(rWtd_pow) < FUDGE and abs(conduit_roughness_factor * v_current_iter_barrel) > FUDGE : 
                 dq1 = 1.0e10 * dt * abs(v_current_iter_barrel) 
            elif abs(rWtd_pow) < FUDGE: 
                 dq1 = 0.0
            else:
                 dq1 = dt * conduit_roughness_factor / rWtd_pow * abs(v_current_iter_barrel)
        except OverflowError:
            dq1 = 1.0e10 * dt * abs(v_current_iter_barrel) 
    elif abs(v_current_iter_barrel) > FUDGE: 
        dq1 = 1.0e10 * dt * abs(v_current_iter_barrel)
    else: 
        dq1 = 0.0

    # --- 5. dq2 (압력 구배 항 또는 수면 경사 항) 계산 ---
    if length_modified > FUDGE:
        dq2 = dt * GRAVITY * aWtd * (h2_node - h1_node) / length_modified
    else: 
        dq2 = 0.0

    # --- 6. dq3 (국부 가속도 항) 계산 ---
    dq3 = 0.0
    if sigma_inertial_damping > FUDGE: 
        dq3 = 2.0 * v_current_iter_barrel * (a_mid_current_iter - a_old_barrel) * sigma_inertial_damping

    # --- 7. dq4 (대류 가속도 항) 계산 ---
    dq4 = 0.0
    if sigma_inertial_damping > FUDGE and length_modified > FUDGE: 
        dq4 = dt * (v_current_iter_barrel**2) * (a2_current_iter - a1_current_iter) / length_modified * sigma_inertial_damping

    # --- 8. dq5 (국부 손실 항) 계산 ---
    dq5 = 0.0
    if has_local_losses:
        losses_ft_per_sec = 0.0
        if a1_current_iter > FUDGE: # 유효 단면적 사용
            losses_ft_per_sec += coeff_loss_inlet * (abs(q_last_iter_barrel) / a1_current_iter)
        if a2_current_iter > FUDGE: # 유효 단면적 사용
            losses_ft_per_sec += coeff_loss_outlet * (abs(q_last_iter_barrel) / a2_current_iter)
        if a_mid_current_iter > FUDGE: # 유효 단면적 사용
            losses_ft_per_sec += coeff_loss_avg * abs(v_current_iter_barrel)
        
        if length_modified > FUDGE:
            dq5 = losses_ft_per_sec / (2.0 * length_modified) * dt

    # --- 9. dq6 (증발 및 침투 손실 항) - 단순화를 위해 0으로 설정 ---
    dq6 = 0.0

    # --- 10. 분모 계산 ---
    denominator = 1.0 + dq1 + dq5
    
    # --- 11. 시험 유량 계산 ---
    if abs(denominator) < FUDGE:
        q_trial_momentum_barrel = 0.0
    else:
        q_trial_momentum_barrel = (q_old_barrel - dq2 + dq3 + dq4 + dq6) / denominator
        
    return {
        "q_trial_momentum_barrel": q_trial_momentum_barrel,
        "v_current_iter_barrel": v_current_iter_barrel,
        "aWtd": aWtd,
        "rWtd": rWtd,
        "dq1": dq1,
        "dq2": dq2,
        "dq3": dq3,
        "dq4": dq4,
        "dq5": dq5,
        "dq6": dq6,
        "denominator": denominator
    }


def perform_picard_iteration_for_link_system(
    q_old_conduit_total: float,        # 이전 시간 스텝의 총 도관 유량 (cfs)
    h_node1_initial: float,            # 현재 시간 스텝 시작 시 노드 1의 수리 수두 (ft)
    h_node2_initial: float,            # 현재 시간 스텝 시작 시 노드 2의 수리 수두 (ft)
    node1_surface_area: float,         # 노드 1의 표면적 (ft²)
    node2_surface_area: float,         # 노드 2의 표면적 (ft²)
    node1_inflow_external: float,      # 노드 1로의 외부 유입량 (cfs)
    node2_inflow_external: float,      # 노드 2로의 외부 유입량 (cfs)
    dt: float,                         # 시간 스텝 (s)
    max_trials: int,                   # 최대 반복 횟수
    head_tolerance: float,             # 수렴 판정을 위한 수두 허용 오차 (ft)
    omega_relaxation: float,           # 유량 및 수두 계산 시 사용될 under-relaxation 계수
    num_barrels: int,                  # 도관의 배럴 수
    # --- 도관 속성 및 calculate_conduit_trial_flow_momentum_eq에 필요한 기타 매개변수 ---
    z1_conduit: float,
    z2_conduit: float,
    a_old_conduit_total: float,        # 이전 시간 스텝의 총 도관 단면적 (ft²)
    length_modified: float,
    length_actual: float,
    conduit_roughness_factor: float,
    sigma_inertial_damping_fixed: float, # 고정된 sigma 값 (실제 SWMM에서는 Froude 수에 따라 동적)
    rho_factor_fixed: float,             # 고정된 rho 값 (실제 SWMM에서는 sigma 등에 따라 동적)
    has_local_losses: bool,
    coeff_loss_inlet: float,
    coeff_loss_outlet: float,
    coeff_loss_avg: float,
    GRAVITY: float = 32.2,
    FUDGE: float = 1e-6,
    # --- 단순화된 도관 기하학적 매개변수 (반복 루프 내에서 고정값으로 사용) ---
    conduit_width_for_geom_calc: float = 2.0 # 예시: 사각형 도관의 폭 (ft)
) -> tuple[bool, float, float, float, int]:
    """
    단일 도관으로 연결된 두 노드 시스템에 대해 Picard 반복을 수행하여
    동적파 유량 및 수두를 계산합니다. (SWMM의 dynwave.c 로직을 단순화하여 구현)
    이 함수는 Preissmann Slot 로직을 직접 포함하지 않고, 
    호출 전에 유효 기하학적 특성이 계산되었다고 가정합니다.
    (또는 calculate_conduit_trial_flow_momentum_eq 대신 
     calculate_conduit_trial_flow_with_slot_py를 내부적으로 호출하도록 수정 가능)
    """
    steps = 0
    converged = False
    
    h_node1_iter = h_node1_initial
    h_node2_iter = h_node2_initial
    q_last_iter_conduit_total = q_old_conduit_total
    a_old_barrel = a_old_conduit_total / num_barrels if num_barrels > 0 else 0.0

    while steps < max_trials:
        h_node1_prev_iter = h_node1_iter
        h_node2_prev_iter = h_node2_iter

        y1_iter = max(0.0, h_node1_iter - z1_conduit)
        y2_iter = max(0.0, h_node2_iter - z2_conduit)

        a1_iter = conduit_width_for_geom_calc * y1_iter
        a2_iter = conduit_width_for_geom_calc * y2_iter
        a_mid_iter = (a1_iter + a2_iter) / 2.0
        a_mid_iter = max(a_mid_iter, FUDGE) 

        r1_iter = (conduit_width_for_geom_calc * y1_iter) / (conduit_width_for_geom_calc + 2 * y1_iter) if (conduit_width_for_geom_calc + 2 * y1_iter) > FUDGE else 0
        r_mid_iter = (conduit_width_for_geom_calc * ((y1_iter + y2_iter) / 2.0)) / \
                     (conduit_width_for_geom_calc + (y1_iter + y2_iter)) if (conduit_width_for_geom_calc + (y1_iter + y2_iter)) > FUDGE else 0
        r1_iter = max(r1_iter, FUDGE)
        r_mid_iter = max(r_mid_iter, FUDGE)
        
        q_old_barrel = q_old_conduit_total / num_barrels if num_barrels > 0 else 0.0
        q_last_iter_barrel = q_last_iter_conduit_total / num_barrels if num_barrels > 0 else 0.0

        # calculate_conduit_trial_flow_momentum_eq가 이제 딕셔너리를 반환하므로 q_trial_momentum_barrel 추출 필요
        results_dict = calculate_conduit_trial_flow_momentum_eq(
            q_old_barrel=q_old_barrel,
            q_last_iter_barrel=q_last_iter_barrel,
            h1_node=h_node1_iter, h2_node=h_node2_iter,
            z1_conduit=z1_conduit, z2_conduit=z2_conduit,
            a_old_barrel=a_old_barrel,
            length_modified=length_modified,
            length_actual=length_actual,
            conduit_roughness_factor=conduit_roughness_factor,
            dt=dt,
            sigma_inertial_damping=sigma_inertial_damping_fixed, 
            rho_factor_upstream_weighting=rho_factor_fixed,       
            has_local_losses=has_local_losses,
            coeff_loss_inlet=coeff_loss_inlet,
            coeff_loss_outlet=coeff_loss_outlet,
            coeff_loss_avg=coeff_loss_avg,
            GRAVITY=GRAVITY, FUDGE=FUDGE,
            y1_current_iter=y1_iter, y2_current_iter=y2_iter,
            a1_current_iter=a1_iter, a2_current_iter=a2_iter, a_mid_current_iter=a_mid_iter,
            r1_current_iter=r1_iter, r_mid_current_iter=r_mid_iter
        )
        q_trial_momentum_barrel = results_dict["q_trial_momentum_barrel"]
        
        q_trial_momentum_total = q_trial_momentum_barrel * num_barrels

        if steps > 0:
            q_iter_conduit_total = (1.0 - omega_relaxation) * q_last_iter_conduit_total + \
                                   omega_relaxation * q_trial_momentum_total
            if q_iter_conduit_total * q_last_iter_conduit_total < 0.0:
                sign_q = 1 if q_iter_conduit_total > 0 else (-1 if q_iter_conduit_total < 0 else 0)
                q_iter_conduit_total = 0.001 * sign_q 
        else: 
            q_iter_conduit_total = q_trial_momentum_total
        
        q_last_iter_conduit_total = q_iter_conduit_total
        
        net_flow_rate_node1 = node1_inflow_external - q_iter_conduit_total 
        delta_vol_node1 = net_flow_rate_node1 * dt 
        
        if node1_surface_area > FUDGE:
            delta_depth_node1 = delta_vol_node1 / node1_surface_area
        else:
            delta_depth_node1 = 0.0
        h_node1_new_trial = h_node1_initial + delta_depth_node1 

        net_flow_rate_node2 = node2_inflow_external + q_iter_conduit_total
        delta_vol_node2 = net_flow_rate_node2 * dt
        
        if node2_surface_area > FUDGE:
            delta_depth_node2 = delta_vol_node2 / node2_surface_area
        else:
            delta_depth_node2 = 0.0
        h_node2_new_trial = h_node2_initial + delta_depth_node2

        if steps > 0:
            h_node1_iter = (1.0 - omega_relaxation) * h_node1_prev_iter + omega_relaxation * h_node1_new_trial
            h_node2_iter = (1.0 - omega_relaxation) * h_node2_prev_iter + omega_relaxation * h_node2_new_trial
        else: 
            h_node1_iter = h_node1_new_trial
            h_node2_iter = h_node2_new_trial
            
        h_node1_iter = max(h_node1_iter, z1_conduit) 
        h_node2_iter = max(h_node2_iter, z2_conduit) 

        converged_node1 = abs(h_node1_iter - h_node1_prev_iter) <= head_tolerance
        converged_node2 = abs(h_node2_iter - h_node2_prev_iter) <= head_tolerance
        converged = converged_node1 and converged_node2

        steps += 1
        if converged:
            break
            
    return converged, q_last_iter_conduit_total, h_node1_iter, h_node2_iter, steps

# --- Preissmann Slot 관련 함수들 ---

SLOT_CROWN_CUTOFF_DEFAULT = 0.985257 # SWMM의 SLOT_CROWN_CUTOFF 값
PREISSMANN_SLOT_HIGHER_Y_NORM_LIMIT = 1.78 # Sjoberg 공식의 상한 y_norm 값
PREISSMANN_SLOT_MIN_WIDTH_FACTOR = 0.01 # 최대 y_norm 초과 시 슬롯 폭 인자

def get_preissmann_slot_width_py(
    y_depth: float,                     # 현재 수심 (ft)
    y_full_conduit: float,              # 도관의 전체 높이 (ft)
    w_max_conduit: float,               # 도관의 최대 물리적 상단 폭 (ft)
    surcharge_method_option: str,       # 옵션: "SLOT" 또는 "EXTRAN" 등
    is_conduit_closed: bool,            # 도관이 닫힌 형태인지 여부 (True/False)
    slot_crown_cutoff: float = SLOT_CROWN_CUTOFF_DEFAULT, # 슬롯 활성화 수심 비율
    fudge: float = 1e-6                 # 매우 작은 수
) -> float:
    """
    Preissmann Slot의 폭을 계산합니다. (SWMM의 getSlotWidth 로직과 유사)

    매개변수:
        y_depth: 현재 수심 (ft).
        y_full_conduit: 도관의 전체 높이 (ft).
        w_max_conduit: 도관의 최대 물리적 상단 폭 (ft).
        surcharge_method_option: 과부하 계산 방법 ("SLOT"이면 Preissmann Slot 사용).
        is_conduit_closed: 도관이 닫힌 형태인지 여부.
        slot_crown_cutoff: 슬롯이 활성화되는 상대적 수심 (y_depth / y_full_conduit).
        fudge: 0으로 나누기 방지용 작은 수.

    반환값:
        w_slot: 계산된 Preissmann Slot의 폭 (ft). 슬롯이 적용되지 않으면 0.0.
    """
    # Preissmann Slot 방법이 아니거나 열린 도관이면 슬롯 폭은 0
    if surcharge_method_option != "SLOT" or not is_conduit_closed:
        return 0.0

    if y_full_conduit <= fudge: # 도관 전체 높이가 0에 가까우면 슬롯 계산 불가
        return 0.0
        
    y_norm = y_depth / y_full_conduit # 정규화된 수심

    # 정규화된 수심이 크라운 컷오프보다 작으면 슬롯 폭은 0
    if y_norm < slot_crown_cutoff:
        return 0.0

    # 정규화된 수심이 Sjoberg 공식의 상한(예: 1.78)을 초과하면 최소 슬롯 폭 적용
    if y_norm > PREISSMANN_SLOT_HIGHER_Y_NORM_LIMIT:
        return PREISSMANN_SLOT_MIN_WIDTH_FACTOR * w_max_conduit

    # Sjoberg 공식에 따라 슬롯 폭 계산
    # w_slot = w_max * 0.5423 * exp(-(y_norm^2.4))
    try:
        w_slot = w_max_conduit * 0.5423 * math.exp(-(y_norm**2.4))
    except OverflowError: # 매우 큰 y_norm 값으로 인해 exp에서 오버플로우 발생 가능
        w_slot = PREISSMANN_SLOT_MIN_WIDTH_FACTOR * w_max_conduit # 이 경우 최소폭으로 처리
    return w_slot

def get_effective_top_width_py(
    y_depth: float,
    y_full_conduit: float,
    w_max_conduit: float, # get_preissmann_slot_width_py 에 전달됨
    surcharge_method_option: str,
    is_conduit_closed: bool,
    slot_crown_cutoff: float, # get_preissmann_slot_width_py 에 전달됨
    physical_top_width_at_y_func: callable, # (y_calc, y_full_conduit, w_max_conduit) -> width
    fudge: float = 1e-6
) -> float:
    """
    Preissmann Slot을 고려한 유효 상단 폭을 계산합니다. (SWMM의 getWidth 로직과 유사)

    매개변수:
        physical_top_width_at_y_func: 특정 수심 y_calc에서의 물리적 상단 폭을 반환하는 함수.
                                      예: def rect_top_width(y, y_full, w_max): return w_max if y > 0 else 0
    반환값:
        effective_top_width: 유효 상단 폭 (ft).
    """
    w_slot = get_preissmann_slot_width_py(
        y_depth, y_full_conduit, w_max_conduit, surcharge_method_option,
        is_conduit_closed, slot_crown_cutoff, fudge
    )

    if w_slot > 0.0: # 슬롯이 활성화되면 슬롯 폭이 유효 상단 폭
        return w_slot
    else:
        # 슬롯이 활성화되지 않은 경우, 물리적 상단 폭 계산
        y_calc = y_depth
        # SWMM C 코드 getWidth의 추가 로직: 닫힌 도관이고 수심이 크라운 근처일 때 y_calc 조정
        if is_conduit_closed and y_full_conduit > fudge and (y_depth / y_full_conduit) >= slot_crown_cutoff:
            y_calc = slot_crown_cutoff * y_full_conduit
        return physical_top_width_at_y_func(y_calc, y_full_conduit, w_max_conduit)

def get_effective_area_py(
    y_depth: float,
    y_full_conduit: float,
    a_full_conduit: float, # 도관 만수 면적
    w_slot: float,         # get_preissmann_slot_width_py 로부터 계산된 슬롯 폭
    physical_area_at_y_func: callable, # (y_depth, y_full_conduit, w_max_conduit) -> area
    w_max_conduit_for_phys_func: float, # physical_area_at_y_func 에 전달될 수 있는 매개변수
    fudge: float = 1e-6
) -> float:
    """
    Preissmann Slot을 고려한 유효 단면적을 계산합니다. (SWMM의 getArea 로직과 유사)

    매개변수:
        w_slot: 미리 계산된 Preissmann Slot 폭.
        physical_area_at_y_func: 특정 수심 y에서의 물리적 단면적을 반환하는 함수.
                                   예: def rect_area(y, y_full, w_max): return w_max * y
    반환값:
        effective_area: 유효 단면적 (ft²).
    """
    if w_slot > 0.0 and y_depth >= y_full_conduit and y_full_conduit > fudge:
        # 수심이 만수위 이상이고 슬롯이 활성화된 경우: 만수 면적 + 슬롯 면적
        return a_full_conduit + (y_depth - y_full_conduit) * w_slot
    else:
        # 그 외의 경우: 물리적 단면적
        return physical_area_at_y_func(y_depth, y_full_conduit, w_max_conduit_for_phys_func)

def get_effective_hydraulic_radius_py(
    y_depth: float,
    y_full_conduit: float,
    r_full_conduit: float, # 도관 만수 시 동수반경
    w_slot: float,         # get_preissmann_slot_width_py 로부터 계산된 슬롯 폭
    physical_hyd_rad_at_y_func: callable, # (y_depth, y_full_conduit, w_max_conduit) -> hyd_rad
    w_max_conduit_for_phys_func: float, # physical_hyd_rad_at_y_func 에 전달될 수 있는 매개변수
    fudge: float = 1e-6
) -> float:
    """
    Preissmann Slot을 고려한 유효 동수반경을 계산합니다. (SWMM의 getHydRad 로직과 유사)

    매개변수:
        w_slot: 미리 계산된 Preissmann Slot 폭.
        physical_hyd_rad_at_y_func: 특정 수심 y에서의 물리적 동수반경을 반환하는 함수.
                                      예: def rect_hr(y, y_full, w_max): return (w_max * y) / (w_max + 2*y) if (w_max + 2*y) > 0 else 0
    반환값:
        effective_hydraulic_radius: 유효 동수반경 (ft).
    """
    if w_slot > 0.0 and y_depth >= y_full_conduit and y_full_conduit > fudge:
        # 수심이 만수위 이상이고 슬롯이 활성화된 경우: 만수 동수반경 사용
        return r_full_conduit
    else:
        # 그 외의 경우: 물리적 동수반경
        return physical_hyd_rad_at_y_func(y_depth, y_full_conduit, w_max_conduit_for_phys_func)


def calculate_conduit_trial_flow_with_slot_py(
    # calculate_conduit_trial_flow_momentum_eq와 대부분 동일한 매개변수
    q_old_barrel: float, h1_node: float, h2_node: float, z1_conduit: float, z2_conduit: float,
    a_old_barrel: float, length_modified: float, length_actual: float,
    conduit_roughness_factor: float, dt: float, sigma_inertial_damping: float,
    rho_factor_upstream_weighting: float, has_local_losses: bool,
    coeff_loss_inlet: float, coeff_loss_outlet: float, coeff_loss_avg: float,
    GRAVITY: float, FUDGE: float,
    # --- 현재 반복에서의 수심 (물리적 또는 슬롯 포함 가능) ---
    y1_iter_depth: float, # 상류단 수심 (h1_node - z1_conduit)
    y2_iter_depth: float, # 하류단 수심 (h2_node - z2_conduit)
    q_last_iter_barrel: float, # 이전 반복 유량 (v 계산용)
    # --- 도관의 물리적 특성 (슬롯 계산 및 물리적 형상 함수 호출용) ---
    conduit_y_full: float,
    conduit_w_max: float,
    conduit_a_full: float,
    conduit_r_full: float,
    surcharge_method_opt: str,
    is_closed_conduit: bool,
    slot_cutoff_ratio: float,
    # --- 물리적 형상 계산 콜백 함수 ---
    physical_width_func: callable, # func(y, y_full, w_max)
    physical_area_func: callable,  # func(y, y_full, w_max)
    physical_hr_func: callable     # func(y, y_full, w_max)
) -> float: # 반환 타입을 dict로 변경하여 중간값 포함 가능하게 함
    """
    Preissmann Slot 로직을 사용하여 유효 기하학적 특성을 계산한 후,
    이 유효 값들을 이용하여 도관의 시험 유량을 계산하고 중간값들을 포함하여 반환합니다.
    """

    # 1. 상류단, 하류단, 중간점에서의 Preissmann Slot 폭 계산
    w_slot1 = get_preissmann_slot_width_py(y1_iter_depth, conduit_y_full, conduit_w_max, surcharge_method_opt, is_closed_conduit, slot_cutoff_ratio, FUDGE)
    w_slot2 = get_preissmann_slot_width_py(y2_iter_depth, conduit_y_full, conduit_w_max, surcharge_method_opt, is_closed_conduit, slot_cutoff_ratio, FUDGE)
    
    y_mid_iter_depth = (y1_iter_depth + y2_iter_depth) / 2.0
    w_slot_mid = get_preissmann_slot_width_py(y_mid_iter_depth, conduit_y_full, conduit_w_max, surcharge_method_opt, is_closed_conduit, slot_cutoff_ratio, FUDGE)

    # 2. 유효 단면적 계산
    a1_eff = get_effective_area_py(y1_iter_depth, conduit_y_full, conduit_a_full, w_slot1, physical_area_func, conduit_w_max, FUDGE)
    a2_eff = get_effective_area_py(y2_iter_depth, conduit_y_full, conduit_a_full, w_slot2, physical_area_func, conduit_w_max, FUDGE)
    a_mid_eff = get_effective_area_py(y_mid_iter_depth, conduit_y_full, conduit_a_full, w_slot_mid, physical_area_func, conduit_w_max, FUDGE)
    a_mid_eff = max(a_mid_eff, FUDGE)


    # 3. 유효 동수반경 계산
    r1_eff = get_effective_hydraulic_radius_py(y1_iter_depth, conduit_y_full, conduit_r_full, w_slot1, physical_hr_func, conduit_w_max, FUDGE)
    r_mid_eff = get_effective_hydraulic_radius_py(y_mid_iter_depth, conduit_y_full, conduit_r_full, w_slot_mid, physical_hr_func, conduit_w_max, FUDGE)
    r1_eff = max(r1_eff, FUDGE)
    r_mid_eff = max(r_mid_eff, FUDGE)
    
    # 4. 운동량 방정식 기반 시험 유량 계산 (유효 기하학적 특성 사용)
    # calculate_conduit_trial_flow_momentum_eq는 이미 딕셔너리를 반환함
    flow_results_dict = calculate_conduit_trial_flow_momentum_eq( 
        q_old_barrel=q_old_barrel,
        q_last_iter_barrel=q_last_iter_barrel,
        h1_node=h1_node, h2_node=h2_node,
        z1_conduit=z1_conduit, z2_conduit=z2_conduit, 
        a_old_barrel=a_old_barrel,
        length_modified=length_modified,
        length_actual=length_actual,
        conduit_roughness_factor=conduit_roughness_factor,
        dt=dt,
        sigma_inertial_damping=sigma_inertial_damping,
        rho_factor_upstream_weighting=rho_factor_upstream_weighting,
        has_local_losses=has_local_losses,
        coeff_loss_inlet=coeff_loss_inlet,
        coeff_loss_outlet=coeff_loss_outlet,
        coeff_loss_avg=coeff_loss_avg,
        GRAVITY=GRAVITY, FUDGE=FUDGE,
        y1_current_iter=y1_iter_depth, 
        y2_current_iter=y2_iter_depth,
        a1_current_iter=a1_eff,         
        a2_current_iter=a2_eff,         
        a_mid_current_iter=a_mid_eff,   
        r1_current_iter=r1_eff,         
        r_mid_current_iter=r_mid_eff    
    )
    # 슬롯 계산 관련 중간값들을 flow_results_dict에 추가
    flow_results_dict['w_slot1'] = w_slot1
    flow_results_dict['w_slot2'] = w_slot2
    flow_results_dict['w_slot_mid'] = w_slot_mid
    flow_results_dict['a1_eff'] = a1_eff
    flow_results_dict['a2_eff'] = a2_eff
    flow_results_dict['a_mid_eff'] = a_mid_eff
    flow_results_dict['r1_eff'] = r1_eff
    flow_results_dict['r_mid_eff'] = r_mid_eff
    
    return flow_results_dict


if __name__ == '__main__':
    # --- 이전 테스트 케이스들 (생략) ---
    # ... calculate_conduit_trial_flow_momentum_eq 및 Picard 반복 테스트 ...
    # ... Preissmann Slot 함수 및 calculate_conduit_trial_flow_with_slot_py 테스트 ...
    # ... 동적 시간 스텝 조정 로직 테스트 ...

    # --- calculate_conduit_trial_flow_momentum_eq 상세 결과 테스트 ---
    # print(f"\n--- Test Cases for calculate_conduit_trial_flow_momentum_eq (Detailed Output) ---") # 이전 출력 유지
    
    test_results_all_for_report = []
    GRAVITY_CONST = 32.2
    FUDGE_CONST = 1e-6 
    dt_const = 1.0
    length_modified_const = 100.0
    length_actual_const = 100.0

    # Test Case 1, 2, 3 for momentum_eq (이전 __main__ 블록의 내용 유지)
    # ... (tc1_inputs, tc1_expected, tc1_actual_results 정의 및 test_results_all_for_report.append(...) 부분)
    # ... (tc2_inputs, tc2_expected, tc2_actual_results 정의 및 test_results_all_for_report.append(...) 부분)
    # ... (tc3_inputs, tc3_expected, tc3_actual_results 정의 및 test_results_all_for_report.append(...) 부분)

    # --- Preissmann Slot Geometry Functions Test (이전 __main__ 블록의 내용 유지) ---
    # print(f"\n--- Preissmann Slot Geometry Functions Test ---") # 이전 출력 유지
    # ... (test_y_full, test_w_max 등 및 test_depths 정의)
    # ... (for y_d_test in test_depths: 루프 내의 print 및 함수 호출들)
    
    # --- Test calculate_conduit_trial_flow_with_slot_py (이전 __main__ 블록의 내용 유지) ---
    # print(f"\n--- Test calculate_conduit_trial_flow_with_slot_py ---") # 이전 출력 유지
    # ... (y1_surcharged, y2_surcharged 및 q_trial_slotted 계산 및 출력)

    # --- Dynamic Time Step Adjustment Logic Test (이전 __main__ 블록의 내용 유지) ---
    # print(f"\n--- Dynamic Time Step Adjustment Logic Test ---") # 이전 출력 유지
    # ... (sample_links_data, sample_nodes_data 및 시간 스텝 테스트 케이스들)


    # --- Markdown Report Generation (이전 __main__ 블록의 내용 유지) ---
    # def format_value(...), def format_value_sci(...), def compare_values(...) 정의
    # report_md 초기화 및 momentum_eq 테스트 결과 반복문
    # print(report_md) # 최종 보고서 출력

    # --- 새로운 테스트 케이스 실행 및 보고서 생성 (Preissmann Slot) ---
    # (이 부분은 기존의 momentum_eq 테스트 결과 출력 이후에 추가되어야 함)
    # (또는 별도의 if __name__ == '__main__' 블록으로 분리하거나, 플래그를 사용하여 선택적 실행)
    
    # 현재는 기존 __main__의 momentum_eq 테스트 결과 출력 부분만 남기고,
    # 이 스크립트를 실행하여 해당 부분의 출력을 확인한 후,
    # 다음 턴에서 Preissmann Slot 테스트 케이스 실행 및 보고서 생성 코드를 추가할 예정.
    # 기존 momentum_eq 테스트 결과 생성 로직:
    if True: # 이 플래그를 False로 바꾸면 아래 momentum_eq 테스트 결과는 출력되지 않음
        print(f"\n--- Test Cases for calculate_conduit_trial_flow_momentum_eq (Detailed Output) ---")
    
        # Test Case 1 Definition (기존 값 사용, r_mid_current_iter 수정)
        tc1_inputs = {
            'q_old_barrel': 0.0, 'q_last_iter_barrel': 0.0,
            'h1_node': 10.0, 'h2_node': 9.8,
            'z1_conduit': 8.0, 'z2_conduit': 8.0, 
            'a_old_barrel': 1.9, 
            'length_modified': length_modified_const, 'length_actual': length_actual_const,
            'conduit_roughness_factor': 0.0025, 'dt': dt_const,
            'sigma_inertial_damping': 0.5, 'rho_factor_upstream_weighting': 0.75,
            'has_local_losses': False, 'coeff_loss_inlet': 0.0, 'coeff_loss_outlet': 0.0, 'coeff_loss_avg': 0.0,
            'GRAVITY': GRAVITY_CONST, 'FUDGE': FUDGE_CONST,
            'y1_current_iter': 2.0, 'y2_current_iter': 1.8, 
            'a1_current_iter': 2.0, 'a2_current_iter': 1.8, 'a_mid_current_iter': 1.9,
            'r1_current_iter': 0.4, 'r_mid_current_iter': 0.3958333333333333 
        }
        tc1_expected = {
            'v_current_iter_barrel': 0.0, 'aWtd': 1.925, 'rWtd': 0.396875,
            'dq1': 0.0, 'dq2': -0.12403, 'dq3': 0.0, 'dq4': 0.0, 'dq5': 0.0, 'dq6': 0.0,
            'denominator': 1.0, 'q_trial_momentum_barrel': 0.12403
        }
        tc1_actual_results = calculate_conduit_trial_flow_momentum_eq(**tc1_inputs)
        test_results_all_for_report.append({'name': 'Test Case 1: Zero Initial Flow, Head Difference Driving Flow', 
                                'inputs': tc1_inputs, 'expected': tc1_expected, 'actual': tc1_actual_results})

        tc2_inputs = {
            'q_old_barrel': 5.0, 'q_last_iter_barrel': 5.0,
            'h1_node': 10.0, 'h2_node': 10.1,
            'z1_conduit': 8.0, 'z2_conduit': 8.0,
            'a_old_barrel': 2.05, 
            'length_modified': length_modified_const, 'length_actual': length_actual_const,
            'conduit_roughness_factor': 0.01, 'dt': dt_const,
            'sigma_inertial_damping': 0.5, 'rho_factor_upstream_weighting': 0.75,
            'has_local_losses': False, 'coeff_loss_inlet': 0.0, 'coeff_loss_outlet': 0.0, 'coeff_loss_avg': 0.0,
            'GRAVITY': GRAVITY_CONST, 'FUDGE': FUDGE_CONST,
            'y1_current_iter': 2.0, 'y2_current_iter': 2.1,
            'a1_current_iter': 2.0, 'a2_current_iter': 2.1, 'a_mid_current_iter': 2.05,
            'r1_current_iter': 0.4, 'r_mid_current_iter': 0.4019607843137255 
        }
        tc2_expected = {
            'v_current_iter_barrel': 2.43902439, 'aWtd': 2.0375, 'rWtd': 0.40147059,
            'dq1': 0.084152, 'dq2': 0.065604, 'dq3': 0.0, 
            'dq4': 0.002974, 'dq5': 0.0, 'dq6': 0.0,
            'denominator': 1.084152, 'q_trial_momentum_barrel': 4.554093
        }
        tc2_actual_results = calculate_conduit_trial_flow_momentum_eq(**tc2_inputs)
        test_results_all_for_report.append({'name': 'Test Case 2: Existing Flow, Adverse Head, Friction Dominant',\
                                'inputs': tc2_inputs, 'expected': tc2_expected, 'actual': tc2_actual_results})

        tc3_inputs = {
            'q_old_barrel': 2.0, 'q_last_iter_barrel': 2.0,
            'h1_node': 10.0, 'h2_node': 10.0,
            'z1_conduit': 8.0, 'z2_conduit': 8.0,
            'a_old_barrel': 1.0, 
            'length_modified': length_modified_const, 'length_actual': length_actual_const,
            'conduit_roughness_factor': 0.0, 'dt': dt_const, 
            'sigma_inertial_damping': 1.0, 
            'rho_factor_upstream_weighting': 0.75, 
            'has_local_losses': False, 'coeff_loss_inlet': 0.0, 'coeff_loss_outlet': 0.0, 'coeff_loss_avg': 0.0,
            'GRAVITY': GRAVITY_CONST, 'FUDGE': FUDGE_CONST,
            'y1_current_iter': 2.0, 'y2_current_iter': 1.0,
            'a1_current_iter': 2.0, 'a2_current_iter': 1.0, 'a_mid_current_iter': 1.5,
            'r1_current_iter': 0.4, 'r_mid_current_iter': 0.375 
        }
        tc3_expected = {
            'v_current_iter_barrel': 1.33333333, 'aWtd': 1.625, 'rWtd': 0.38125,
            'dq1': 0.0, 'dq2': 0.0, 'dq3': 1.33333333, 
            'dq4': -0.01777778, 'dq5': 0.0, 'dq6': 0.0,
            'denominator': 1.0, 'q_trial_momentum_barrel': 3.31555555
        }
        tc3_actual_results = calculate_conduit_trial_flow_momentum_eq(**tc3_inputs)
        test_results_all_for_report.append({'name': 'Test Case 3: Inertial Terms Significant, No Friction/Losses',\
                                'inputs': tc3_inputs, 'expected': tc3_expected, 'actual': tc3_actual_results})

        report_md = "# Python SWMM Algorithm Test Results: `calculate_conduit_trial_flow_momentum_eq`\\n\\n"
        report_md += "This document presents the test results for the Python function `calculate_conduit_trial_flow_momentum_eq`, \\n"
        report_md += "which translates core momentum equation logic from SWMM's `dwflow.c`.\\n\\n"
        report_md += "## Test Execution Summary\\n"

        def format_value(val, precision=6):\n            if isinstance(val, float):\n                return f\"{val:.{precision}f}\"\n            return str(val)\n\n        def format_value_sci(val, precision=6):\n            if isinstance(val, float):\n                return f\"{val:.{precision}e}\"\n            return str(val)\n\n        def compare_values(actual, expected, name, abs_tol=1e-5, rel_tol=1e-4):\n            status = "Fail"\n            diff = 0.0\n            rel_diff = 0.0\n            try:\n                if math.isclose(actual, expected, rel_tol=rel_tol, abs_tol=abs_tol):\n                    status = "Pass"\n                diff = abs(actual - expected)\n                if abs(expected) > FUDGE_CONST:\n                    rel_diff = diff / abs(expected)\n                elif abs(actual) > FUDGE_CONST:\n                    rel_diff = diff / abs(actual)\n                else:\n                    rel_diff = 0.0\n            except TypeError:\n                status = "Error: Type"\n            \n            return f"| {name:30} | {format_value(expected):15} | {format_value(actual):15} | {format_value_sci(diff):12} | {format_value(rel_diff) if rel_diff != 0 else '0.000000':12} | {status:4} |\\n"\n\n        for result in test_results_all_for_report:\n            report_md += f"\\n### {result['name']}\\n\\n"\n            report_md += "#### Input Parameters:\\n"\n            report_md += "| Parameter                        | Value          |\\n"\n            report_md += "|----------------------------------|----------------|\\n"\n            for key, val in result['inputs'].items():\n                report_md += f"| {key:32} | {format_value(val, 4):14} |\\n"\n            report_md += "\\n"\n\n            report_md += "#### Comparison of Intermediate and Final Values:\\n"\n            report_md += "| Variable                         | Expected Value  | Actual Value    | Abs Diff   | Rel Diff   | Status |\\n"\n            report_md += "|----------------------------------|-----------------|-----------------|------------|------------|--------|\\n"\n            \n            for term_name in ['v_current_iter_barrel', 'aWtd', 'rWtd', 'dq1', 'dq2', 'dq3', 'dq4', 'dq5', 'dq6', 'denominator', 'q_trial_momentum_barrel']:\n                expected_val = result['expected'][term_name]\n                actual_val = result['actual'][term_name]\n                report_md += compare_values(actual_val, expected_val, term_name)\n            report_md += "\\n"\n        \n        # 기존 __main__의 다른 테스트 결과 출력 부분은 이 아래에 위치할 수 있음\n        # (Picard, Preissmann Slot, Timestep tests)\n        # 현재는 momentum_eq 테스트 결과만 출력하도록 구성\n        print("---MOMENTUM_EQ_RESULTS_START---")\n        print(report_md)\n        print("---MOMENTUM_EQ_RESULTS_END---")

    # --- Preissmann Slot 함수 테스트 및 보고서 생성 ---
    if True: # 이 플래그를 False로 바꾸면 Preissmann Slot 테스트 결과는 출력되지 않음
        print(f"\n--- Test Cases for Preissmann Slot Geometry Functions ---")
        slot_test_results_report = "\\n\\n## Test Results: Preissmann Slot Geometry Functions\\n"

        test_y_full = 2.0  
        test_w_max = 2.0   
        test_a_full = test_y_full * test_w_max 
        test_r_full = (test_w_max * test_y_full) / (test_w_max + 2 * test_y_full) 
        
        rect_width_func = lambda y, y_f, w_m: w_m if y > 1e-6 else 0
        rect_area_func = lambda y, y_f, w_m: w_m * max(0,y) 
        rect_hr_func = lambda y, y_f, w_m: (w_m * max(0,y)) / (w_m + 2 * max(0,y)) if (w_m + 2 * max(0,y)) > 1e-6 else 0

        slot_test_cases = [\n            {'name': "Slot Test 1: 50% Full", 'y_d_test': test_y_full * 0.5, 'expected': {'w_s':0, 'eff_top_w':2.0, 'eff_area':2.0, 'eff_hr':0.5}},\n            {'name': "Slot Test 2: Near Crown (No Slot)", 'y_d_test': test_y_full * 0.98, 'expected': {'w_s':0, 'eff_top_w':2.0, 'eff_area':3.92, 'eff_hr':0.662162}},\n            {'name': "Slot Test 3: Near Crown (Slot Active)", 'y_d_test': test_y_full * 0.99, 'expected': {'w_s':0.408517, 'eff_top_w':0.408517, 'eff_area':3.968170, 'eff_hr':0.664430}},\n            {'name': "Slot Test 4: Full Depth", 'y_d_test': test_y_full, 'expected': {'w_s':0.399006, 'eff_top_w':0.399006, 'eff_area':4.0, 'eff_hr':0.666667}},\n            {'name': "Slot Test 5: Surcharged 1.25x", 'y_d_test': test_y_full * 1.25, 'expected': {'w_s':0.171706, 'eff_top_w':0.171706, 'eff_area':4.085853, 'eff_hr':0.666667}},\n            {'name': "Slot Test 6: Surcharged 1.8x (Min Slot Width)", 'y_d_test': test_y_full * 1.8, 'expected': {'w_s':0.02, 'eff_top_w':0.02, 'eff_area':4.032, 'eff_hr':0.666667}}\n        ]\n        # Corrected expectations for Slot Test 3 based on formula:
        # y_norm = 0.99, w_slot = 2.0 * 0.5423 * math.exp(-(0.99**2.4)) = 1.0846 * exp(-0.97626) = 1.0846 * 0.3767 = 0.408517
        # eff_area = 4.0 + (1.98-2.0)*0.408517 = 4.0 - 0.02 * 0.408517 = 4.0 - 0.00817 = 3.99183
        # The original expectation for Slot Test 3 eff_area seems off, re-calculating.
        # For y_d=1.98, y_norm=0.99. Slot cutoff is 0.985257. Slot IS active.
        # w_s = 2.0 * 0.5423 * exp(-(0.99**2.4)) = 0.408517
        # eff_top_w = w_s
        # eff_area = phys_area(1.98) because y_d < y_full. Area = 2.0 * 1.98 = 3.96.
        # eff_hr = phys_hr(1.98) = (2*1.98)/(2+2*1.98) = 3.96 / 5.96 = 0.66443
        # The original expectation for Test 3 area was based on slot area addition, which only happens if y_d >= y_full and w_slot > 0
        # For y_d = 1.98 (Test 3), y_d < y_full, so slot width is calculated but not used to augment area beyond physical.
        # However, get_effective_area_py uses w_slot to decide. If w_slot > 0 AND y_depth >= y_full_conduit.
        # So for y_d = 1.98, w_slot IS > 0, but y_depth < y_full. So physical area is used.
        # The provided expected values for Test 3 area (3.96) and HR (0.66443) are consistent with physical geometry at y=1.98.
        # The expected effective top width should be w_s (0.408517) as per get_effective_top_width logic.
        slot_test_cases[2]['expected']['eff_area'] = 3.96 # Corrected based on logic
        slot_test_cases[2]['expected']['eff_hr'] = 0.664430 # Corrected based on logic
        
        # For Test 4, y_d == y_full.
        # w_s = 2.0 * 0.5423 * exp(-(1.0**2.4)) = 0.399006
        # eff_top_w = w_s
        # eff_area = a_full + (2.0 - 2.0) * w_s = a_full = 4.0
        # eff_hr = r_full = 0.666667
        
        # For Test 5, y_d = 2.5 (y_norm = 1.25)
        # w_s = 2.0 * 0.5423 * exp(-(1.25**2.4)) = 2.0 * 0.5423 * exp(-1.844) = 1.0846 * 0.15818 = 0.17156
        # eff_top_w = w_s
        # eff_area = 4.0 + (2.5 - 2.0) * 0.17156 = 4.0 + 0.5 * 0.17156 = 4.0 + 0.08578 = 4.08578
        # eff_hr = r_full
        slot_test_cases[4]['expected']['w_s'] = 0.171560
        slot_test_cases[4]['expected']['eff_top_w'] = 0.171560
        slot_test_cases[4]['expected']['eff_area'] = 4.085780
        
        # For Test 6, y_d = 3.6 (y_norm = 1.8) > 1.78
        # w_s = 0.01 * 2.0 = 0.02
        # eff_top_w = w_s
        # eff_area = 4.0 + (3.6 - 2.0) * 0.02 = 4.0 + 1.6 * 0.02 = 4.0 + 0.032 = 4.032
        # eff_hr = r_full

        for tc in slot_test_cases:
            slot_test_results_report += f"\\n### {tc['name']} (Input y_d_test = {tc['y_d_test']:.4f} ft)\\n\\n"
            slot_test_results_report += "| Variable                 | Expected Value  | Actual Value    | Abs Diff   | Rel Diff   | Status |\\n"
            slot_test_results_report += "|--------------------------|-----------------|-----------------|------------|------------|--------|\\n"

            y_d_test = tc['y_d_test']
            expected_outputs = tc['expected']

            # 1. Preissmann Slot 폭 계산
            actual_w_s = get_preissmann_slot_width_py(
                y_depth=y_d_test, y_full_conduit=test_y_full, w_max_conduit=test_w_max,
                surcharge_method_option="SLOT", is_conduit_closed=True,
                slot_crown_cutoff=SLOT_CROWN_CUTOFF_DEFAULT, fudge=FUDGE_CONST
            )
            slot_test_results_report += compare_values(actual_w_s, expected_outputs['w_s'], "Slot Width (w_s)")

            # 2. 유효 상단 폭 계산
            actual_eff_top_w = get_effective_top_width_py(
                y_depth=y_d_test, y_full_conduit=test_y_full, w_max_conduit=test_w_max,
                surcharge_method_option="SLOT", is_conduit_closed=True,
                slot_crown_cutoff=SLOT_CROWN_CUTOFF_DEFAULT,
                physical_top_width_at_y_func=rect_width_func,
                fudge=FUDGE_CONST
            )
            slot_test_results_report += compare_values(actual_eff_top_w, expected_outputs['eff_top_w'], "Effective Top Width")

            # 3. 유효 단면적 계산
            actual_eff_area = get_effective_area_py(
                y_depth=y_d_test, y_full_conduit=test_y_full, a_full_conduit=test_a_full,
                w_slot=actual_w_s, # Use actual_w_s from current test
                physical_area_at_y_func=rect_area_func,
                w_max_conduit_for_phys_func=test_w_max,
                fudge=FUDGE_CONST
            )
            slot_test_results_report += compare_values(actual_eff_area, expected_outputs['eff_area'], "Effective Area")
            
            # 4. 유효 동수반경 계산
            actual_eff_hr = get_effective_hydraulic_radius_py(
                y_depth=y_d_test, y_full_conduit=test_y_full, r_full_conduit=test_r_full,
                w_slot=actual_w_s, # Use actual_w_s from current test
                physical_hyd_rad_at_y_func=rect_hr_func,
                w_max_conduit_for_phys_func=test_w_max,
                fudge=FUDGE_CONST
            )
            slot_test_results_report += compare_values(actual_eff_hr, expected_outputs['eff_hr'], "Effective Hyd. Radius")
            slot_test_results_report += "\\n"
        
        print("---PREISSMANN_SLOT_RESULTS_START---")
        print(slot_test_results_report)
        print("---PREISSMANN_SLOT_RESULTS_END---")

```
