import numpy as np

# Define global variable DEBUG
DEBUG = False

def compute_p(m_exp, k):
    log2_m = - m_exp / (np.log2(k) + 1)
    return float(np.exp2(log2_m))

def chernoff_bounds(p, n, confidence):
    """
    Calculate Chernoff bounds for given probability p, number of trials n, and confidence level.
    
    :param p: Expected probability (empirical success probability)
    :param n: Number of trials
    :param confidence: Desired confidence level (e.g., 0.99 for 99% confidence)
    :return: Lower and upper bounds
    """
    epsilon = np.sqrt(-2 * np.log(1 - confidence) / (n * p))
    lower_bound = max(0, p - epsilon * p)
    upper_bound = min(1, p + epsilon * p)
    return lower_bound, upper_bound

def SetRange(s):
    return 2.0*(s//2)+1

def ProbUnif(s):
    '''
    Returns the probability mass on an element in the support of the uniform distribution over the range defined by range_size
    '''
    return (1.0/SetRange(s))

def ProbSumToZ(s, z):
    denom = 2 * (s // 2) + 1
    result = 1.0 / denom - float(abs(z)) / (denom ** 2)
    return float(result)

def ProbSumInRange(s, p):
    if not (0 <= p <= 1):
        raise ValueError("p must be in [0, 1]")
    thres = (s * p) // 2
    denom = 2 * (s // 2) + 1
    
    result_0 = (thres / denom) - ((thres + 1) * thres / 2.0) / (denom ** 2)
    
    result = 2.0 * result_0 + 1.0 / denom
    return result

def ProbSumModSInRange(s, p):
    if not (0 <= p <= 1):
        raise ValueError("p must be in [0, 1]")
    t = (s * p) // 2
    n = 2 * t + 1
    d = 2 * (s // 2) + 1
    return float(n)/d

def MRDistFromUnif(s, p):
    range = s * p
    numerator = ProbSumToZ(s, 0)/ProbSumInRange(s, p)
    denominator = ProbSumToZ(s, (range // 2))/ProbSumInRange(s, p)
    prob_unit = ProbUnif(range)
    if denominator == 0.0 or prob_unit == 0.0:
        raise ZeroDivisionError("Denominator is zero in MRDistFromUnif")
    return max(numerator/prob_unit, prob_unit/denominator)

def SoS(n):
    return n * (n + 1) * (2 * n + 1) / 6

def ProbSumWithTwoRVInRange(s, p):
    if not (0 <= p <= 1):
        raise ValueError("p must be in [0, 1]")
    
    s2 = s // 2
    t = (s*p // 2)
    range_s = SetRange(s)
    range_sp = SetRange(s*p)
    sum1 = (2 * (s2 - t) / range_s) * (range_sp/range_s)**2
    sum2_1 = SoS(range_sp - 1) - SoS(range_sp - t - 1)
    sum2 = (2.0 / range_s) * ((sum2_1/range_s)/range_s)
    
    return sum1 + sum2

def MRDistFromPairUnif(s, p):
    if not (0 <= p <= 0.5):
        raise ValueError("p must be in [0, 0.5]")
    sp_half = s * p // 2
    s_half = s // 2
    range_s = SetRange(s)
    range_sp = SetRange(s*p)
    pr2 = ProbSumWithTwoRVInRange(s, p)
    max_term = (range_sp/range_s)**2 / pr2
    min_term = (range_s/range_sp)**2/(2.0*(s_half - sp_half)) * pr2 * range_s
    return max(max_term, min_term)

def FirstMomentInductFactors(m_exp, k, p, d, zm):
    # check whether k is a power of 2
    if not (k & (k - 1) == 0):
        raise ValueError("k must be a power of 2")
    log2_k = int(np.log2(k))
    if not (0 <= d <= log2_k - 1):
        raise ValueError("d must be in [0, log2(k) - 1]")
    new_m_exp = m_exp * (1 - d / (log2_k + 1))
    s = 2.0 ** new_m_exp
    temp_exp = k / (2.0 ** (d + 1))
    if zm == 1 and d == 0:  # First iteration for Z_m version
        term1_noexp = ProbSumModSInRange(s, p)
        term2_ub_noexp = 1.0
        term2_lb_noexp = 1.0
    else:
        term1_noexp = ProbSumInRange(s, p)
        term2_ub_noexp = MRDistFromUnif(s, p)
        term2_lb_noexp = 1.0 / term2_ub_noexp
    log2_ub = temp_exp * float((np.log2(term1_noexp) + np.log2(term2_ub_noexp)))
    log2_lb = temp_exp * float((np.log2(term1_noexp) + np.log2(term2_lb_noexp)))
    return log2_ub, log2_lb

def FirstMomentBounds(m_exp, k, p, zm):
    log2_k = int(np.log2(k))
    log2_ub = float(np.log2(ProbUnif(2.0**(m_exp / (log2_k + 1)))))
    log2_lb = log2_ub
    for d in range(0, log2_k):
        log2_ub_temp, log2_lb_temp = FirstMomentInductFactors(m_exp, k, p, d, zm)
        log2_ub += log2_ub_temp
        log2_lb += log2_lb_temp
    return log2_ub, log2_lb

def FirstMomentSumBounds(m_exp, k, n, zm):
    log2_k = float(np.log2(k))
    p = 2.0 ** (- m_exp / (log2_k + 1))
    log2_ub, log2_lb = FirstMomentBounds(m_exp, k, p, zm)
    ub = np.exp2(log2_ub + k * np.log2(n))
    lb = np.exp2(log2_lb + k * np.log2(n))
    return float(ub), float(lb)

def ProbUnifXMRDist(m_exp, log_p, zm):
    s = 2.0 ** m_exp
    p = 2.0 ** log_p
    return ProbSumInRange(s, p) * MRDistFromUnif(s, p) if zm == 0 else ProbSumModSInRange(s, p)

def ProbUnifXMRDistPair(m_exp, log_p, zm):
    s = 2.0 ** m_exp
    p = 2.0 ** log_p
    return ProbSumWithTwoRVInRange(s, p) * MRDistFromPairUnif(s, p) if zm == 0 else (ProbSumModSInRange(s, p))**2

def SecondMomentUBRecurse(m_exp, k, log2_n, log_p, zm):
    """_summary_

    Args:
        m_exp (int): the exponent of m
        k (int): the total number of lists
        n (float): the number of elements in the input list
        log_p (float): log_p = - m_exp / (log2(k) + 1)
    """
    # check k is a power of 2
    k = int(k)
    if not (k & (k - 1) == 0):
        raise ValueError("k must be a power of 2")
    if k > 1:
        term_1 = ProbUnifXMRDist(m_exp, log_p, zm)
        term_2 = ProbUnifXMRDistPair(m_exp, log_p, zm)
        n = 2.0 ** (log2_n + 2*float(np.log2(term_1)))
        log2_n = 3/2 * log2_n + np.log2(n + 2*term_2)/2
        return SecondMomentUBRecurse(m_exp + log_p, k // 2, log2_n, log_p, 0)
    else:
        m = 2.0 ** m_exp
        prob_unif = ProbUnif(m)
        log2_np = log2_n + float(np.log2(prob_unif))
        return log2_np

def SecondMomentUB(m_exp, k, n, zm):
    if not (k & (k - 1) == 0):
        raise ValueError("k must be a power of 2")
    log2_k = float(np.log2(k))
    log2_p = - m_exp / (log2_k + 1)
    
    # make log2_n
    log2_n = float(np.log2(n))
    log2_np = SecondMomentUBRecurse(m_exp, k, log2_n, log2_p, zm)
    
    term_np = np.exp2(log2_np)
    
    return float(term_np * (term_np + 1))

def SizeBoundInductFactors(m_exp, k, p, d, t, zm):
    log2_k = int(np.log2(k))
    if not (0 <= d <= log2_k):
        raise ValueError("d must be in [0, log2(k) - 1]")
    if not (0 <= t <= d - 2):
        raise ValueError("t must be in [0, d - 2]")
    new_m_exp = m_exp * (1 - float(t) / (log2_k + 1))
    s = 2.0 ** new_m_exp
    temp_exp = 2.0 ** (d - t - 1)
    if zm == 1 and t == 0:  # First iteration for Z_m version
        term1_noexp = ProbSumModSInRange(s, p)
        term2_ub_noexp = 1.0
        term2_lb_noexp = 1.0
    else:
        term1_noexp = ProbSumInRange(s, p)
        term2_ub_noexp = MRDistFromUnif(s, p)
        term2_lb_noexp = 1.0 / term2_ub_noexp
    log2_ub = temp_exp * float((np.log2(term1_noexp) + np.log2(term2_ub_noexp)))
    log2_lb = temp_exp * float((np.log2(term1_noexp) + np.log2(term2_lb_noexp)))
    return log2_ub, log2_lb

def SizeBounds(m_exp, k, n, zm):
    log2_k = int(np.log2(k))
    log2_n = float(np.log2(n))
    p = 2.0 ** ( - m_exp / (log2_k + 1) )
    ub = 0.0
    lb = 0.0
    for d in range(0, log2_k + 1):
        if d == 0:
            ub += k * n
            lb += k * n
        elif d == 1:
            if zm == 1:
                additive_term = k/2 * ProbSumModSInRange(2.0**m_exp, p) * n**2
            else:
                additive_term = k/2 * ProbSumInRange(2.0**m_exp, p) * n**2
            ub += additive_term
            lb += additive_term
        else:
            new_m_exp_1 = m_exp * (1 - (d - 1) / (log2_k + 1))
            s_1 = 2.0 ** new_m_exp_1
            term_0 = ProbSumInRange(s_1, p)
            log2_term_0_ub = np.log2(term_0)
            log2_term_0_lb = np.log2(term_0)
            for t in range(0, d - 1):
                log2_ub_temp, log2_lb_temp = SizeBoundInductFactors(m_exp, k, p, d, t, zm)
                log2_term_0_ub += log2_ub_temp
                log2_term_0_lb += log2_lb_temp
            ub += k/(2.0**d) * 2.0 ** ((2.0**d)*log2_n + log2_term_0_ub)
            lb += k/(2.0**d) * 2.0 ** ((2.0**d)*log2_n + log2_term_0_lb)
    return ub, lb

def main_theorem(k, n, m_exp, mode="success_prob", zm=0, remain_float=False):
    # Ensure k is a power of 2 and at least 4
    k = int(k)
    n = float(n)
    m_exp = float(m_exp)

    log2_k = float(np.log2(k))

    # Compute p
    p = 2.0 ** ( - m_exp / (log2_k + 1) )

    # Compute c
    c = p * n

    if mode == "success_prob":
        n = float(n) if remain_float else np.round(max(n, 1.0))
        # Compute upper and lower bounds from FirstMomentSumBounds
        ub, lb = FirstMomentSumBounds(m_exp, k, n, zm)

        # Compute second moment upper bound from prop_2_4
        second_moment_ub = SecondMomentUB(m_exp, k, n, zm)

        if DEBUG:
            print(f"m_exp={m_exp}, k={k}, n={n}, p={p}, c={c}, ub={ub}, lb={lb}, second_moment_ub={second_moment_ub}")

        # Compute success probability bounds
        success_prob_lb = max((float(lb) ** 2) / second_moment_ub, 0.0) if second_moment_ub > 0.0 else 0.0
        success_prob_ub = min(ub, 1.0)
        
        return success_prob_ub, success_prob_lb
    
    elif mode == "size":
        n = float(n) if remain_float else np.round(max(n, 1.0))
        # Compute size bounds from SizeBounds
        size_ub, size_lb = SizeBounds(m_exp, k, n, zm)
        
        if DEBUG:
            print(f"success_prob_ub={success_prob_ub}, success_prob_lb={success_prob_lb}, size_ub={size_ub}, size_lb={size_lb}")

        return size_ub, size_lb
    
    else:
        raise ValueError("Invalid mode. Use either 'success_prob' or 'size'")

def ProbBoundsAnaly(m, k, n):
    n = max(np.round(n), 1)
    # Check if k is a power of 2 and >= 4
    if not (k & (k-1) == 0) or k < 4:
        raise ValueError("k must be a power of 2 and >= 4")
    
    log_k = np.log2(k)
    
    # Check if m > 30^(log(k)+1)
    if m <= 30**(log_k + 1):
        raise ValueError("m must be > 30^(log(k)+1)")
    
    # Compute p and c
    p = m**(-1 / (log_k + 1))
    c = p * n
    
    # Compute lower bound
    lb = 1 / (c**(-k) + (1 + k/n)**k)
    
    # Compute upper bound
    ub = c**k * (1 + 37*p)**k
    
    return min(ub, 1.0), max(lb, 0.0)

def SizeBoundsAnaly(m, k, n):
    # Check if k is a power of 2 and >= 4
    if not (k & (k-1) == 0) or k < 4:
        raise ValueError("k must be a power of 2 and >= 4")
    
    log_k = np.log2(k)
    
    # Check if m > 30^(log(k)+1)
    if m <= 30**(log_k + 1):
        raise ValueError("m must be > 30^(log(k)+1)")
    
    # Compute p and c
    p = m**(-1 / (log_k + 1))
    c = p * n
    
    # Compute the sum term
    sum_term = sum(c**(2**d - 1) / 2**d for d in range(int(log_k) + 1))
    
    # Compute common factor
    common_factor = k * n * (1 + sum_term)
    
    # Compute lower bound
    lb = common_factor * (1 - 37*p)**(k-1)
    
    # Compute upper bound
    ub = common_factor * (1 + 37*p)**(k-1)
    
    return ub, lb

def ShallueProbBounds(m, k, n):
    return 1, 0

def ShallueSizeBounds(m, k, n):
    m = float(m)
    k = float(k)
    #alpha = max( (1024.0*np.log2(m))/(np.log2(np.e)*np.log2(k)) , 1024.0, k)
    alpha = max( (np.log2(m))/(np.log2(np.e)*np.log2(k)) , 1024.0, k)
    #alpha = max( 1024.0, k)
    threshold = alpha * (m**(1.0/np.log2(k)))
    complexity = k * threshold
    return complexity, complexity

def JouxProbBounds(m, k, n):
    return 1, 0

def JouxSizeBounds(m, k, n):
    alpha = 0.79591
    c = 1/6.0
    l = np.log2(float(k))
    threshold = ( (1.0/alpha)**( (l-1)/(l+1) ) ) * ( (1.0/c)**( (l-1)*(l+2)/(2*(l+1)) ) ) * ( m**( 1.0 / (np.log2(k)+1) ) )
    complexity = k * threshold
    return complexity, complexity
