def concentration(m_eff, T):
    return 2*(2*pi*m_eff*k*T/(2*pi*h_bar)**2)**(3/2)

# w depletion area width
def depletion_width(eps, n, d_phi):
    """Ширина зоны между металлом и полупроводником"""
    return sqrt(eps*d_phi/(2*pi*e*n))

# L_D screening length
def debye_length(eps, n, T):
    """Длина Дебая"""
    return sqrt(eps*k*T/(4*pi*e**2*n))
